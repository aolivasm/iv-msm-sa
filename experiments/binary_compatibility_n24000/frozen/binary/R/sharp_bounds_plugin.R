# Complete standalone cross-fitted plug-in sharp-bound estimator
#
# This file contains everything needed to:
#   1. Generate the simulation data.
#   2. Compute oracle e(X), pi_z(X), gamma_0(Y,X), and gamma_1(Y,X).
#   3. Estimate e(X), pi_z(X) with stable cross-fitted GLM/spline models.
#   4. Estimate gamma_0 and gamma_1 with the direct RKHS density-ratio method.
#   5. Estimate the sharp bounds for mu^1, mu^0, and the ATE.
#   6. Compare:
#        - oracle_all:    oracle e, pi_z, and gamma_t;
#        - oracle_gamma:  estimated e, pi_z and oracle gamma_t;
#        - rkhs:          estimated e, pi_z and normalized RKHS gamma_t.
#
# The production default contains no classification estimator. An optional
# fit-provider hook permits explicitly sourced audit learners; the separate
# classifier-ratio implementation does not change the direct-RKHS default.
#
# The file returns cross-fitted plug-in estimates and layered cross-fitted
# one-step estimates. The production one-step entry point currently requires
# cross-validated natural-spline thresholds; the RKHS label in its fully
# estimated diagnostic mode refers only to the gamma_t learner. The
# threshold functions Q_t^- and Q_t^+ are fitted by minimizing the appropriate
# observation-specific check-loss objectives. The default is the original
# fixed natural-spline class; cross-validated spline and regularized RBF-RKHS
# learners are optional. The primary bounds evaluate the original functional
# R(lambda) directly on held-out observations, and the check-loss pseudo-bound
# plus the exact optimizer-residual discrepancy are returned as diagnostics.
#
# No external R package is required.
#
# Recommended debugging run:
#
#   sim_data <- generate_iv_data(n = 2000, seed = 2026)
#
#   bound_results <- estimate_sharp_bounds_plugin(
#     data = sim_data,
#     Lambda0_T1 = 5,
#     Lambda1_T1 = 5,
#     Lambda0_T0 = 5,
#     Lambda1_T0 = 5,
#     x_cols = "X",
#     K = 5,
#     rkhs_lambda_by_t = NULL,
#     seed = 2026
#   )
#
#   bound_results$summary
#   bound_results$pseudo_summary
#   bound_results$truth
#   bound_results$unconfounded_checks
#   bound_results$normalization_by_x_bin
#   bound_results$fold_diagnostics
#
# Cross-validation is the default. Supplying a named vector to
# rkhs_lambda_by_t uses fixed penalties for faster debugging runs.
#
# =============================================================================
# 1. Basic utilities and simulation DGP
# =============================================================================

clip_probability <- function(p, epsilon = 0.05) {
  pmin(1 - epsilon, pmax(epsilon, p))
}

clip_positive <- function(x, floor = 1e-6, cap = Inf) {
  x[!is.finite(x)] <- floor
  x <- pmax(floor, x)
  if (is.finite(cap)) {
    x <- pmin(cap, x)
  }
  x
}

check_binary <- function(x, name) {
  values <- sort(unique(stats::na.omit(x)))
  if (!all(values %in% c(0, 1))) {
    stop(sprintf("%s must be coded 0/1.", name))
  }
}

rnormtc <- function(n, mu, s = 1) {
  lower <- mu - 2 * s
  upper <- mu + 2 * s

  F_lower <- stats::pnorm(lower, mean = mu, sd = s)
  F_upper <- stats::pnorm(upper, mean = mu, sd = s)

  u <- stats::runif(n, min = F_lower, max = F_upper)
  stats::qnorm(u, mean = mu, sd = s)
}

dtnormtc <- function(x, mu = 0, s = 1) {
  standardized <- (x - mu) / s
  normalizing_constant <- stats::pnorm(2) - stats::pnorm(-2)

  ifelse(
    standardized >= -2 & standardized <= 2,
    stats::dnorm(standardized) / (s * normalizing_constant),
    0
  )
}

get_dgp_outcome_support <- function(treatment_effect = 2) {
  if (length(treatment_effect) != 1L || !is.finite(treatment_effect)) {
    stop("treatment_effect must be one finite number.")
  }
  # X, U, and each outcome error lie in [-2, 2]. Therefore
  # Y^0 = 1.5 X + 1.5 U + error_0 lies in [-8, 8], while
  # Y^1 = treatment_effect + 1.5 X + 1.5 U + error_1 has the same
  # width and is translated by the additive treatment effect.
  list(`0` = c(-8, 8), `1` = treatment_effect + c(-8, 8))
}

generate_iv_data <- function(
    n = 5000,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    seed = NULL
) {
  if (length(treatment_effect) != 1L || !is.finite(treatment_effect)) {
    stop("treatment_effect must be one finite number.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  X <- rnormtc(n, mu = 0, s = 1)
  U <- rnormtc(n, mu = 0, s = 1)

  prob_Z <- stats::plogis(0.5 * X)
  Z <- stats::rbinom(n, 1, prob_Z)

  prob_T <- stats::plogis(
    treatment_intercept + X + instrument_strength * Z +
      confounder_strength * U
  )
  T_obs <- stats::rbinom(n, 1, prob_T)

  # Symmetric truncated-normal outcome errors keep the causal means unchanged
  # while giving known bounded supports for IV-only partial identification.
  mean_Y0 <- 1.5 * X + 1.5 * U
  mean_Y1 <- treatment_effect + 1.5 * X + 1.5 * U
  Y0 <- rnormtc(n, mu = mean_Y0, s = 1)
  Y1 <- rnormtc(n, mu = mean_Y1, s = 1)
  Y <- ifelse(T_obs == 1, Y1, Y0)

  data.frame(
    Y = Y,
    T_obs = T_obs,
    Z = Z,
    X = X,
    U = U,
    Y0 = Y0,
    Y1 = Y1,
    prob_Z = prob_Z,
    prob_T = prob_T
  )
}

d_u <- function(u) {
  dtnormtc(u, mu = 0, s = 1)
}

# =============================================================================
# 2. Oracle quantities by Gauss-Legendre quadrature
# =============================================================================

gauss_legendre_rule <- function(
    n_nodes = 80,
    lower = -2,
    upper = 2
) {
  if (n_nodes < 2) {
    stop("n_nodes must be at least 2.")
  }

  index <- seq_len(n_nodes - 1)
  beta <- index / sqrt(4 * index^2 - 1)

  jacobi <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  jacobi[cbind(index, index + 1)] <- beta
  jacobi[cbind(index + 1, index)] <- beta

  eig <- eigen(jacobi, symmetric = TRUE)
  ordering <- order(eig$values)

  nodes_standard <- eig$values[ordering]
  weights_standard <- 2 * eig$vectors[1, ordering]^2

  list(
    nodes = (upper - lower) * nodes_standard / 2 +
      (upper + lower) / 2,
    weights = (upper - lower) * weights_standard / 2
  )
}

compute_oracle_quantities <- function(
    data,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    n_nodes = 80,
    batch_size = 1000
) {
  required <- c("Y", "X")
  missing_columns <- setdiff(required, names(data))

  if (length(missing_columns) > 0) {
    stop(
      "Missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  rule <- gauss_legendre_rule(
    n_nodes = n_nodes,
    lower = -2,
    upper = 2
  )
  density_rule <- gauss_legendre_rule(
    n_nodes = n_nodes,
    lower = -1,
    upper = 1
  )

  u <- rule$nodes
  integration_weights <- rule$weights * d_u(u)

  n <- nrow(data)

  pi0_true <- rep(NA_real_, n)
  pi1_true <- rep(NA_real_, n)
  e_true <- rep(NA_real_, n)
  gamma0_oracle <- rep(NA_real_, n)
  gamma1_oracle <- rep(NA_real_, n)

  starts <- seq.int(1, n, by = batch_size)

  for (start in starts) {
    finish <- min(start + batch_size - 1, n)
    index <- start:finish

    x <- data$X[index]
    y <- data$Y[index]
    n_batch <- length(index)

    x_matrix <- matrix(
      x,
      nrow = n_batch,
      ncol = n_nodes
    )
    u_matrix <- matrix(
      u,
      nrow = n_batch,
      ncol = n_nodes,
      byrow = TRUE
    )

    eta_z0 <- treatment_intercept +
      x_matrix +
      confounder_strength * u_matrix

    eta_z1 <- eta_z0 + instrument_strength

    p_t1_z0 <- stats::plogis(eta_z0)
    p_t1_z1 <- stats::plogis(eta_z1)

    weighted_integral <- function(matrix_value) {
      as.vector(matrix_value %*% integration_weights)
    }

    pi0_batch <- weighted_integral(p_t1_z0)
    pi1_batch <- weighted_integral(p_t1_z1)

    pi0_true[index] <- pi0_batch
    pi1_true[index] <- pi1_batch

    # e(X) = P(Z=0 | X).
    e_true[index] <- stats::plogis(-0.5 * x)

    for (t_value in c(0, 1)) {
      outcome_location <- treatment_effect * t_value + 1.5 * x

      # The truncated outcome density is nonzero only when
      # |y - outcome_location - 1.5 u| <= 2. Integrating on that
      # observation-specific interval avoids applying a global quadrature rule
      # to an integrand with moving discontinuities.
      density_u_lower <- pmax(
        -2,
        (y - outcome_location - 2) / 1.5
      )
      density_u_upper <- pmin(
        2,
        (y - outcome_location + 2) / 1.5
      )
      density_half_width <- pmax(
        0,
        (density_u_upper - density_u_lower) / 2
      )
      density_midpoint <- (density_u_upper + density_u_lower) / 2
      density_u_matrix <- matrix(
        density_midpoint,
        nrow = n_batch,
        ncol = n_nodes
      ) + outer(
        density_half_width,
        density_rule$nodes
      )
      density_quadrature_weights <- outer(
        density_half_width,
        density_rule$weights
      ) * d_u(density_u_matrix)

      density_x_matrix <- matrix(
        x,
        nrow = n_batch,
        ncol = n_nodes
      )
      density_y_matrix <- matrix(
        y,
        nrow = n_batch,
        ncol = n_nodes
      )
      mean_y <- treatment_effect * t_value +
        1.5 * density_x_matrix +
        1.5 * density_u_matrix

      density_y_given_xu <- dtnormtc(
        density_y_matrix,
        mu = mean_y,
        s = 1
      )

      density_eta_z0 <- treatment_intercept +
        density_x_matrix +
        confounder_strength * density_u_matrix
      density_eta_z1 <- density_eta_z0 + instrument_strength
      density_p_t1_z0 <- stats::plogis(density_eta_z0)
      density_p_t1_z1 <- stats::plogis(density_eta_z1)

      treatment_weight_z0 <- if (t_value == 1) {
        density_p_t1_z0
      } else {
        1 - density_p_t1_z0
      }

      treatment_weight_z1 <- if (t_value == 1) {
        density_p_t1_z1
      } else {
        1 - density_p_t1_z1
      }

      joint_y_t_z0 <- rowSums(
        density_y_given_xu *
          treatment_weight_z0 *
          density_quadrature_weights
      )

      joint_y_t_z1 <- rowSums(
        density_y_given_xu *
          treatment_weight_z1 *
          density_quadrature_weights
      )

      probability_t_z0 <- if (t_value == 1) {
        pi0_batch
      } else {
        1 - pi0_batch
      }

      probability_t_z1 <- if (t_value == 1) {
        pi1_batch
      } else {
        1 - pi1_batch
      }

      density_y_t_z0 <- joint_y_t_z0 /
        pmax(probability_t_z0, 1e-12)

      density_y_t_z1 <- joint_y_t_z1 /
        pmax(probability_t_z1, 1e-12)

      gamma_value <- ifelse(
        density_y_t_z0 > 1e-300,
        density_y_t_z1 / density_y_t_z0,
        1
      )

      if (t_value == 1) {
        gamma1_oracle[index] <- gamma_value
      } else {
        gamma0_oracle[index] <- gamma_value
      }
    }
  }

  data.frame(
    pi0_true = pi0_true,
    pi1_true = pi1_true,
    e_true = e_true,
    gamma0_oracle = gamma0_oracle,
    gamma1_oracle = gamma1_oracle
  )
}

add_oracle_quantities <- function(
    data,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    n_nodes = 80,
    batch_size = 1000
) {
  oracle <- compute_oracle_quantities(
    data = data,
    confounder_strength = confounder_strength,
    treatment_intercept = treatment_intercept,
    instrument_strength = instrument_strength,
    treatment_effect = treatment_effect,
    n_nodes = n_nodes,
    batch_size = batch_size
  )

  cbind(data, oracle)
}

# Slow scalar oracle, useful for spot checks.
get_oracle_gamma_analytical <- function(
    x,
    y,
    t_value,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2
) {
  probability_t <- function(z, u) {
    p_t1 <- stats::plogis(
      treatment_intercept + x + instrument_strength * z +
        confounder_strength * u
    )

    if (t_value == 1) {
      p_t1
    } else {
      1 - p_t1
    }
  }

  density_y_given_xu <- function(u) {
    dtnormtc(
      y,
      mu = treatment_effect * t_value + 1.5 * x + 1.5 * u,
      s = 1
    )
  }

  conditional_density <- function(z) {
    denominator <- stats::integrate(
      function(u) {
        probability_t(z, u) * d_u(u)
      },
      lower = -2,
      upper = 2
    )$value

    outcome_location <- treatment_effect * t_value + 1.5 * x
    density_u_lower <- max(-2, (y - outcome_location - 2) / 1.5)
    density_u_upper <- min(2, (y - outcome_location + 2) / 1.5)

    if (density_u_lower >= density_u_upper) {
      return(0)
    }

    numerator <- stats::integrate(
      function(u) {
        density_y_given_xu(u) *
          probability_t(z, u) *
          d_u(u)
      },
      lower = density_u_lower,
      upper = density_u_upper
    )$value

    numerator / denominator
  }

  density_z0 <- conditional_density(0)
  density_z1 <- conditional_density(1)
  if (density_z0 <= 0) 1 else density_z1 / density_z0
}

# =============================================================================
# 3. Cross-fitting folds
# =============================================================================

make_tz_stratified_folds <- function(
    data,
    K = 5,
    seed = 1
) {
  check_binary(data$T_obs, "T_obs")
  check_binary(data$Z, "Z")

  if (K < 2) {
    stop("K must be at least 2.")
  }

  set.seed(seed)

  strata <- interaction(
    data$T_obs,
    data$Z,
    drop = TRUE
  )

  fold_id <- integer(nrow(data))

  for (stratum in levels(strata)) {
    index <- which(strata == stratum)

    if (length(index) < K) {
      stop(sprintf(
        "Stratum %s has only %d observations, fewer than K=%d.",
        stratum,
        length(index),
        K
      ))
    }

    fold_id[index] <- sample(
      rep(seq_len(K), length.out = length(index))
    )
  }

  fold_id
}

make_z_stratified_folds <- function(
    data,
    K = 3,
    seed = 1
) {
  if (K < 2) {
    stop("K must be at least 2.")
  }

  set.seed(seed)
  fold_id <- integer(nrow(data))

  for (z_value in c(0, 1)) {
    index <- which(data$Z == z_value)

    if (length(index) < K) {
      stop(sprintf(
        "Only %d observations with Z=%d, fewer than K=%d.",
        length(index),
        z_value,
        K
      ))
    }

    fold_id[index] <- sample(
      rep(seq_len(K), length.out = length(index))
    )
  }

  fold_id
}

# =============================================================================
# 4. Stable nuisance estimation
# =============================================================================

make_additive_spline_formula <- function(
    response,
    predictors,
    spline_df = 5,
    data = NULL
) {
  terms <- vapply(predictors, function(name) {
    use_linear <- !is.null(data) &&
      length(unique(data[[name]][is.finite(data[[name]])])) <= spline_df + 1L
    if (use_linear) name else sprintf("splines::ns(%s, df=%d)", name, spline_df)
  }, character(1))

  stats::as.formula(
    paste(response, "~", paste(terms, collapse = " + "))
  )
}

fit_nuisance_models_glm <- function(
    train,
    x_cols = "X",
    positivity_epsilon = 0.05,
    treatment_spline_df = 5
) {
  # Exact instrument model for this DGP:
  # P(Z=1|X) = expit(0.5 X).
  instrument_formula <- stats::reformulate(
    x_cols,
    response = "Z"
  )

  instrument_fit <- stats::glm(
    formula = instrument_formula,
    data = train,
    family = stats::binomial(),
    control = stats::glm.control(maxit = 100)
  )

  # Separate flexible marginal treatment models by instrument level.
  # This allows the marginal pi_z(X) curves to have different shapes.
  treatment_formula <- make_additive_spline_formula(
    response = "T_obs",
    predictors = x_cols,
    spline_df = treatment_spline_df,
    data = train
  )

  treatment_fit_z0 <- stats::glm(
    formula = treatment_formula,
    data = train[train$Z == 0, , drop = FALSE],
    family = stats::binomial(),
    control = stats::glm.control(maxit = 100)
  )

  treatment_fit_z1 <- stats::glm(
    formula = treatment_formula,
    data = train[train$Z == 1, , drop = FALSE],
    family = stats::binomial(),
    control = stats::glm.control(maxit = 100)
  )

  predict_nuisance <- function(newdata) {
    p_z1 <- stats::predict(
      instrument_fit,
      newdata = newdata,
      type = "response"
    )

    pi0 <- stats::predict(
      treatment_fit_z0,
      newdata = newdata,
      type = "response"
    )

    pi1 <- stats::predict(
      treatment_fit_z1,
      newdata = newdata,
      type = "response"
    )

    list(
      e = clip_probability(
        1 - p_z1,
        positivity_epsilon
      ),
      pi0 = clip_probability(
        pi0,
        positivity_epsilon
      ),
      pi1 = clip_probability(
        pi1,
        positivity_epsilon
      )
    )
  }

  list(
    predict = predict_nuisance,
    instrument_fit = instrument_fit,
    treatment_fit_z0 = treatment_fit_z0,
    treatment_fit_z1 = treatment_fit_z1,
    learner = "glm_spline",
    diagnostics = data.frame(
      component = c("e", "pi0", "pi1"),
      learner = "glm_spline",
      converged = c(
        isTRUE(instrument_fit$converged),
        isTRUE(treatment_fit_z0$converged),
        isTRUE(treatment_fit_z1$converged)
      ),
      model_rank = c(
        instrument_fit$rank,
        treatment_fit_z0$rank,
        treatment_fit_z1$rank
      ),
      effective_df = c(
        instrument_fit$rank,
        treatment_fit_z0$rank,
        treatment_fit_z1$rank
      ),
      stringsAsFactors = FALSE
    )
  )
}

resolve_nuisance_gam_basis_dimension <- function(
    data,
    x_cols,
    basis_dimension = NULL,
    minimum = 5L,
    maximum = 20L,
    growth_exponent = 1 / 3
) {
  if (!is.null(basis_dimension)) {
    basis_dimension <- as.integer(basis_dimension)
    if (length(basis_dimension) != 1L || is.na(basis_dimension) ||
        basis_dimension < 3L) {
      stop("nuisance_gam_k must be NULL or one integer of at least 3.")
    }
    return(basis_dimension)
  }
  minimum <- as.integer(minimum)
  maximum <- as.integer(maximum)
  if (length(minimum) != 1L || is.na(minimum) || minimum < 3L ||
      length(maximum) != 1L || is.na(maximum) || maximum < minimum ||
      length(growth_exponent) != 1L || !is.finite(growth_exponent) ||
      growth_exponent <= 0 || growth_exponent > 1) {
    stop("Invalid adaptive nuisance-GAM basis controls.")
  }
  observations_per_covariate <- nrow(data) / max(1L, length(x_cols))
  proposed <- ceiling(observations_per_covariate^growth_exponent)
  as.integer(min(maximum, max(minimum, proposed)))
}

fit_nuisance_models_gam <- function(
    train,
    x_cols = "X",
    positivity_epsilon = 0.05,
    basis_dimension = NULL,
    minimum_basis_dimension = 5L,
    maximum_basis_dimension = 20L,
    basis_growth_exponent = 1 / 3,
    method = "REML"
) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop(
      "The optional nuisance_learner='gam' requires the recommended mgcv package."
    )
  }
  non_numeric <- x_cols[
    !vapply(train[, x_cols, drop = FALSE], is.numeric, logical(1))
  ]
  if (length(non_numeric) > 0L) {
    stop(
      "GAM nuisance predictors must be numeric: ",
      paste(non_numeric, collapse = ", ")
    )
  }
  if (length(method) != 1L || !nzchar(method)) {
    stop("nuisance_gam_method must be one nonempty string.")
  }

  basis_for <- function(data) {
    resolve_nuisance_gam_basis_dimension(
      data = data,
      x_cols = x_cols,
      basis_dimension = basis_dimension,
      minimum = minimum_basis_dimension,
      maximum = maximum_basis_dimension,
      growth_exponent = basis_growth_exponent
    )
  }
  make_formula <- function(response, data, selected_basis_dimension) {
    smooth_terms <- vapply(x_cols, function(variable) {
      unique_count <- length(unique(data[[variable]][
        is.finite(data[[variable]])
      ]))
      if (unique_count <= 2L) {
        return(sprintf("`%s`", variable))
      }
      local_k <- min(selected_basis_dimension, unique_count - 1L)
      local_k <- max(3L, local_k)
      sprintf("s(`%s`, bs='cr', k=%d)", variable, local_k)
    }, character(1))
    formula <- stats::as.formula(
      paste0("`", response, "` ~ ", paste(smooth_terms, collapse = " + "))
    )
    environment(formula) <- asNamespace("mgcv")
    formula
  }

  e_basis_dimension <- basis_for(train)
  z0_training <- train[train$Z == 0, , drop = FALSE]
  z1_training <- train[train$Z == 1, , drop = FALSE]
  pi0_basis_dimension <- basis_for(z0_training)
  pi1_basis_dimension <- basis_for(z1_training)
  instrument_fit <- mgcv::gam(
    formula = make_formula("Z", train, e_basis_dimension),
    data = train,
    family = stats::binomial(),
    method = method,
    select = TRUE
  )
  treatment_fit_z0 <- mgcv::gam(
    formula = make_formula("T_obs", z0_training, pi0_basis_dimension),
    data = z0_training,
    family = stats::binomial(),
    method = method,
    select = TRUE
  )
  treatment_fit_z1 <- mgcv::gam(
    formula = make_formula("T_obs", z1_training, pi1_basis_dimension),
    data = z1_training,
    family = stats::binomial(),
    method = method,
    select = TRUE
  )

  predict_nuisance <- function(newdata) {
    p_z1 <- stats::predict(
      instrument_fit,
      newdata = newdata,
      type = "response"
    )
    pi0 <- stats::predict(
      treatment_fit_z0,
      newdata = newdata,
      type = "response"
    )
    pi1 <- stats::predict(
      treatment_fit_z1,
      newdata = newdata,
      type = "response"
    )
    list(
      e = clip_probability(1 - p_z1, positivity_epsilon),
      pi0 = clip_probability(pi0, positivity_epsilon),
      pi1 = clip_probability(pi1, positivity_epsilon)
    )
  }

  gam_diagnostics <- function(fit, component, selected_basis_dimension) {
    data.frame(
      component = component,
      learner = "gam",
      converged = isTRUE(fit$converged),
      model_rank = fit$rank,
      effective_df = sum(fit$edf),
      basis_dimension = selected_basis_dimension,
      stringsAsFactors = FALSE
    )
  }

  list(
    predict = predict_nuisance,
    instrument_fit = instrument_fit,
    treatment_fit_z0 = treatment_fit_z0,
    treatment_fit_z1 = treatment_fit_z1,
    learner = "gam",
    diagnostics = rbind(
      gam_diagnostics(instrument_fit, "e", e_basis_dimension),
      gam_diagnostics(treatment_fit_z0, "pi0", pi0_basis_dimension),
      gam_diagnostics(treatment_fit_z1, "pi1", pi1_basis_dimension)
    ),
    configuration = list(
      basis_selection = if (is.null(basis_dimension)) {
        "adaptive_training_stratum_size"
      } else {
        "fixed"
      },
      minimum_basis_dimension = as.integer(minimum_basis_dimension),
      maximum_basis_dimension = as.integer(maximum_basis_dimension),
      basis_growth_exponent = basis_growth_exponent,
      selected_basis_dimensions = c(
        e = e_basis_dimension,
        pi0 = pi0_basis_dimension,
        pi1 = pi1_basis_dimension
      )
    )
  )
}

fit_nuisance_models_ranger <- function(
    train,
    x_cols = "X",
    positivity_epsilon = 0.05,
    num_trees = 500L,
    min_node_size = NULL,
    seed = 1L
) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop(
      "The optional nuisance_learner='ranger' requires the ranger package."
    )
  }
  non_numeric <- x_cols[
    !vapply(train[, x_cols, drop = FALSE], is.numeric, logical(1))
  ]
  if (length(non_numeric) > 0L) {
    stop(
      "Ranger nuisance predictors must be numeric: ",
      paste(non_numeric, collapse = ", ")
    )
  }
  num_trees <- as.integer(num_trees)
  if (length(num_trees) != 1L || is.na(num_trees) || num_trees < 1L) {
    stop("num_trees must be one positive integer.")
  }

  fit_probability_forest <- function(data, response, local_seed) {
    outcome <- as.numeric(data[[response]])
    check_binary(outcome, response)
    if (length(unique(outcome)) == 1L) {
      probability <- (sum(outcome) + 0.5) / (length(outcome) + 1)
      return(list(
        predict = function(newdata) rep(probability, nrow(newdata)),
        fit = NULL,
        min_node_size = NA_integer_
      ))
    }
    local_min_node_size <- if (is.null(min_node_size)) {
      max(5L, as.integer(floor(nrow(data)^(1 / 4))))
    } else {
      as.integer(min_node_size)
    }
    if (length(local_min_node_size) != 1L ||
        is.na(local_min_node_size) || local_min_node_size < 1L) {
      stop("min_node_size must be NULL or one positive integer.")
    }
    ranger_data <- data[, c(response, x_cols), drop = FALSE]
    ranger_data[[response]] <- factor(outcome, levels = c(0, 1))
    fit <- ranger::ranger(
      formula = stats::reformulate(x_cols, response = response),
      data = ranger_data,
      probability = TRUE,
      num.trees = num_trees,
      mtry = max(1L, floor(sqrt(length(x_cols)))),
      min.node.size = local_min_node_size,
      seed = as.integer(local_seed),
      respect.unordered.factors = "order",
      num.threads = 1L
    )
    list(
      predict = function(newdata) {
        prediction <- predict(fit, data = newdata)$predictions
        if (is.null(dim(prediction))) return(as.numeric(prediction))
        if ("1" %in% colnames(prediction)) {
          return(as.numeric(prediction[, "1"]))
        }
        as.numeric(prediction[, ncol(prediction)])
      },
      fit = fit,
      min_node_size = local_min_node_size
    )
  }

  z0_training <- train[train$Z == 0, , drop = FALSE]
  z1_training <- train[train$Z == 1, , drop = FALSE]
  instrument_fit <- fit_probability_forest(train, "Z", seed)
  treatment_fit_z0 <- fit_probability_forest(
    z0_training, "T_obs", seed + 1L
  )
  treatment_fit_z1 <- fit_probability_forest(
    z1_training, "T_obs", seed + 2L
  )

  predict_nuisance <- function(newdata) {
    list(
      e = clip_probability(
        1 - instrument_fit$predict(newdata), positivity_epsilon
      ),
      pi0 = clip_probability(
        treatment_fit_z0$predict(newdata), positivity_epsilon
      ),
      pi1 = clip_probability(
        treatment_fit_z1$predict(newdata), positivity_epsilon
      )
    )
  }

  list(
    predict = predict_nuisance,
    instrument_fit = instrument_fit$fit,
    treatment_fit_z0 = treatment_fit_z0$fit,
    treatment_fit_z1 = treatment_fit_z1$fit,
    learner = "ranger",
    diagnostics = data.frame(
      component = c("e", "pi0", "pi1"),
      learner = "ranger",
      converged = TRUE,
      model_rank = rep(length(x_cols), 3L),
      effective_df = rep(NA_real_, 3L),
      n_training = c(nrow(train), nrow(z0_training), nrow(z1_training)),
      min_node_size = c(
        instrument_fit$min_node_size,
        treatment_fit_z0$min_node_size,
        treatment_fit_z1$min_node_size
      ),
      num_trees = rep(num_trees, 3L),
      stringsAsFactors = FALSE
    )
  )
}

binary_cross_entropy <- function(outcome, probability) {
  outcome <- as.numeric(outcome)
  probability <- pmin(1 - 1e-8, pmax(1e-8, as.numeric(probability)))
  if (length(outcome) != length(probability) ||
      any(!is.finite(outcome)) || any(!is.finite(probability))) {
    stop("Invalid inputs to binary_cross_entropy().")
  }
  mean(-outcome * log(probability) - (1 - outcome) * log1p(-probability))
}

fit_nuisance_models_library_cv <- function(
    train,
    x_cols = "X",
    positivity_epsilon = 0.05,
    candidates = c("glm_spline", "gam", "ranger"),
    cv_folds = 3L,
    selection_rule = c("one_se", "minimum"),
    treatment_spline_df = 5L,
    gam_basis_dimension = NULL,
    gam_minimum_basis_dimension = 5L,
    gam_maximum_basis_dimension = 20L,
    gam_basis_growth_exponent = 1 / 3,
    gam_method = "REML",
    ranger_num_trees = 500L,
    seed = 1L
) {
  selection_rule <- match.arg(selection_rule)
  candidates <- unique(as.character(candidates))
  allowed <- c("glm_spline", "gam", "ranger")
  if (length(candidates) < 1L || any(!candidates %in% allowed)) {
    stop(
      "Primitive nuisance library candidates must be drawn from: ",
      paste(allowed, collapse = ", "), "."
    )
  }
  cv_folds <- as.integer(cv_folds)
  if (length(cv_folds) != 1L || is.na(cv_folds) || cv_folds < 2L) {
    stop("Primitive nuisance library CV requires at least two folds.")
  }

  fit_candidate <- function(data, candidate, local_seed) {
    switch(
      candidate,
      glm_spline = fit_nuisance_models_glm(
        train = data,
        x_cols = x_cols,
        positivity_epsilon = positivity_epsilon,
        treatment_spline_df = treatment_spline_df
      ),
      gam = fit_nuisance_models_gam(
        train = data,
        x_cols = x_cols,
        positivity_epsilon = positivity_epsilon,
        basis_dimension = gam_basis_dimension,
        minimum_basis_dimension = gam_minimum_basis_dimension,
        maximum_basis_dimension = gam_maximum_basis_dimension,
        basis_growth_exponent = gam_basis_growth_exponent,
        method = gam_method
      ),
      ranger = fit_nuisance_models_ranger(
        train = data,
        x_cols = x_cols,
        positivity_epsilon = positivity_epsilon,
        num_trees = ranger_num_trees,
        seed = local_seed
      )
    )
  }

  fold_id <- make_tz_stratified_folds(
    train, K = cv_folds, seed = as.integer(seed)
  )
  risk_rows <- list()
  counter <- 1L
  for (candidate_index in seq_along(candidates)) {
    candidate <- candidates[candidate_index]
    for (fold in seq_len(cv_folds)) {
      validation <- fold_id == fold
      fit <- fit_candidate(
        train[!validation, , drop = FALSE],
        candidate,
        as.integer(seed + 10000L * candidate_index + fold)
      )
      prediction <- fit$predict(train[validation, , drop = FALSE])
      local_data <- train[validation, , drop = FALSE]
      component_specifications <- list(
        e = list(
          index = rep(TRUE, nrow(local_data)),
          outcome = 1 - local_data$Z,
          prediction = prediction$e
        ),
        pi0 = list(
          index = local_data$Z == 0,
          outcome = local_data$T_obs,
          prediction = prediction$pi0
        ),
        pi1 = list(
          index = local_data$Z == 1,
          outcome = local_data$T_obs,
          prediction = prediction$pi1
        )
      )
      for (component in names(component_specifications)) {
        specification <- component_specifications[[component]]
        index <- specification$index
        if (!any(index)) {
          stop("A primitive nuisance CV fold contains an empty stratum.")
        }
        risk_rows[[counter]] <- data.frame(
          component = component,
          candidate = candidate,
          fold = fold,
          risk = binary_cross_entropy(
            specification$outcome[index],
            specification$prediction[index]
          ),
          observations = sum(index),
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }
  }
  fold_risk <- do.call(rbind, risk_rows)
  groups <- split(
    fold_risk,
    interaction(
      fold_risk$component, fold_risk$candidate,
      drop = TRUE, lex.order = TRUE
    )
  )
  cv_table <- do.call(rbind, lapply(groups, function(group) {
    data.frame(
      component = group$component[1L],
      candidate = group$candidate[1L],
      mean_risk = mean(group$risk),
      standard_error = if (nrow(group) > 1L) {
        stats::sd(group$risk) / sqrt(nrow(group))
      } else {
        0
      },
      stringsAsFactors = FALSE
    )
  }))
  rownames(cv_table) <- NULL
  complexity <- c(glm_spline = 1L, gam = 2L, ranger = 3L)
  selected <- stats::setNames(character(3L), c("e", "pi0", "pi1"))
  cv_table$selected <- FALSE
  for (component in names(selected)) {
    local <- which(cv_table$component == component)
    minimum <- local[which.min(cv_table$mean_risk[local])]
    choice <- minimum
    if (selection_rule == "one_se") {
      cutoff <- cv_table$mean_risk[minimum] +
        cv_table$standard_error[minimum]
      eligible <- local[cv_table$mean_risk[local] <= cutoff]
      choice <- eligible[which.min(complexity[cv_table$candidate[eligible]])]
    }
    cv_table$selected[choice] <- TRUE
    selected[component] <- cv_table$candidate[choice]
  }

  selected_candidates <- unique(unname(selected))
  final_fits <- stats::setNames(lapply(seq_along(selected_candidates), function(index) {
    fit_candidate(
      train,
      selected_candidates[index],
      as.integer(seed + 900000L + index)
    )
  }), selected_candidates)

  predict_nuisance <- function(newdata) {
    candidate_predictions <- lapply(final_fits, function(fit) {
      fit$predict(newdata)
    })
    list(
      e = candidate_predictions[[selected[["e"]]]]$e,
      pi0 = candidate_predictions[[selected[["pi0"]]]]$pi0,
      pi1 = candidate_predictions[[selected[["pi1"]]]]$pi1
    )
  }

  selected_diagnostics <- do.call(rbind, lapply(names(selected), function(component) {
    candidate <- selected[[component]]
    diagnostic <- final_fits[[candidate]]$diagnostics
    diagnostic <- diagnostic[diagnostic$component == component, , drop = FALSE]
    local_cv <- cv_table[
      cv_table$component == component & cv_table$candidate == candidate,
      , drop = FALSE
    ]
    value_or <- function(name, default) {
      if (name %in% names(diagnostic)) diagnostic[[name]][1L] else default
    }
    data.frame(
      component = component,
      learner = "library_cv",
      selected_learner = candidate,
      converged = value_or("converged", TRUE),
      model_rank = value_or("model_rank", NA_real_),
      effective_df = value_or("effective_df", NA_real_),
      n_training = value_or(
        "n_training",
        if (component == "e") nrow(train) else sum(
          train$Z == if (component == "pi0") 0 else 1
        )
      ),
      cv_mean_risk = local_cv$mean_risk[1L],
      cv_standard_error = local_cv$standard_error[1L],
      stringsAsFactors = FALSE
    )
  }))
  rownames(selected_diagnostics) <- NULL

  list(
    predict = predict_nuisance,
    fits = final_fits,
    learner = "library_cv",
    selected_learners = selected,
    diagnostics = selected_diagnostics,
    cv_table = cv_table,
    configuration = list(
      candidates = candidates,
      cv_folds = cv_folds,
      selection_rule = selection_rule,
      selected_learners = selected
    )
  )
}

# Convex Super Learner for the three primitive Bernoulli nuisances used by
# both continuous- and discrete-outcome estimators.  The implementation of
# the individual library fits and the Bernoulli-log-loss metalearner lives in
# R/binary_outcome_bounds.R because it is also used for the binary outcome
# regressions.  Nothing here estimates gamma_t through classification odds:
# gamma_t remains the direct RKHS/Riesz density-ratio fit below.
fit_nuisance_models_super_learner <- function(
    train,
    x_cols = "X",
    positivity_epsilon = 0.05,
    candidates = NULL,
    cv_folds = 3L,
    ranger_num_trees = 500L,
    hal_max_degree = 2L,
    seed = 1L
) {
  required <- c(
    "binary_fit_regression_super_learner",
    "binary_default_super_learner_candidates"
  )
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing)) {
    stop(
      "Source R/binary_outcome_bounds.R before requesting the primitive ",
      "Super Learner; missing: ", paste(missing, collapse = ", "), "."
    )
  }
  if (is.null(candidates)) {
    candidates <- binary_default_super_learner_candidates()
  }
  candidates <- unique(as.character(candidates))
  fit_component <- function(data, response, local_seed) {
    binary_fit_regression_super_learner(
      data = data,
      response = response,
      x_cols = x_cols,
      candidates = candidates,
      cv_folds = cv_folds,
      ranger_num_trees = ranger_num_trees,
      hal_max_degree = hal_max_degree,
      seed = local_seed
    )
  }

  z0_training <- train[train$Z == 0, , drop = FALSE]
  z1_training <- train[train$Z == 1, , drop = FALSE]
  fits <- list(
    e = fit_component(train, "Z", seed + 1L),
    pi0 = fit_component(z0_training, "T_obs", seed + 2L),
    pi1 = fit_component(z1_training, "T_obs", seed + 3L)
  )
  training_sizes <- c(
    e = nrow(train), pi0 = nrow(z0_training), pi1 = nrow(z1_training)
  )

  diagnostics <- do.call(rbind, lapply(names(fits), function(component) {
    fit <- fits[[component]]
    row <- data.frame(
      component = component,
      learner = "super_learner",
      selected_learner = fit$selected_learner,
      converged = isTRUE(fit$converged),
      model_rank = as.integer(fit$rank),
      effective_df = NA_real_,
      n_training = as.integer(training_sizes[[component]]),
      cv_mean_risk = as.numeric(fit$cv_mean_risk),
      stringsAsFactors = FALSE
    )
    weights <- fit$ensemble_weights
    risk_table <- fit$cv_table
    for (candidate in candidates) {
      row[[paste0("weight_", candidate)]] <- if (
        !is.null(weights) && candidate %in% names(weights)
      ) as.numeric(weights[[candidate]]) else 0
      local_risk <- if (
        !is.null(risk_table) && nrow(risk_table) &&
          candidate %in% risk_table$candidate
      ) {
        risk_table$mean_risk[match(candidate, risk_table$candidate)]
      } else {
        NA_real_
      }
      row[[paste0("risk_", candidate)]] <- as.numeric(local_risk)
    }
    row
  }))
  rownames(diagnostics) <- NULL

  cv_tables <- do.call(rbind, lapply(names(fits), function(component) {
    table <- fits[[component]]$cv_table
    if (is.null(table) || !nrow(table)) return(NULL)
    table$component <- component
    table
  }))

  predict_nuisance <- function(newdata) {
    list(
      e = clip_probability(
        1 - fits$e$predict(newdata), positivity_epsilon
      ),
      pi0 = clip_probability(
        fits$pi0$predict(newdata), positivity_epsilon
      ),
      pi1 = clip_probability(
        fits$pi1$predict(newdata), positivity_epsilon
      )
    )
  }

  list(
    predict = predict_nuisance,
    fits = fits,
    learner = "super_learner",
    diagnostics = diagnostics,
    cv_table = cv_tables,
    configuration = list(
      candidates = candidates,
      cv_folds = as.integer(cv_folds),
      ranger_num_trees = as.integer(ranger_num_trees),
      hal_max_degree = as.integer(hal_max_degree),
      metalearner = "simplex_probability_mixture_bernoulli_log_loss"
    )
  )
}

fit_nuisance_models_integrated_parametric <- function(
    train,
    x_cols = "X",
    positivity_epsilon = 0.05,
    confounder_strength = 1,
    quadrature_nodes = 40L,
    maxit = 500L
) {
  # Simulation benchmark for the correctly specified primitive nuisances.
  # The instrument model is logistic and linear in X. The marginal treatment
  # propensity integrates the known truncated-normal U distribution out of
  #   logit P(T=1 | X,Z,U) = alpha + beta'X + eta Z + c U.
  # The confounding strength c and the distribution of U are treated as known;
  # alpha, beta, and eta are estimated by marginal Bernoulli likelihood.
  non_numeric <- x_cols[
    !vapply(train[, x_cols, drop = FALSE], is.numeric, logical(1))
  ]
  if (length(non_numeric) > 0L) {
    stop(
      "Integrated-parametric nuisance predictors must be numeric: ",
      paste(non_numeric, collapse = ", ")
    )
  }
  if (
    length(confounder_strength) != 1L ||
      !is.finite(confounder_strength)
  ) {
    stop("confounder_strength must be one finite number.")
  }
  if (
    length(quadrature_nodes) != 1L ||
      !is.finite(quadrature_nodes) ||
      quadrature_nodes < 2L ||
      abs(quadrature_nodes - round(quadrature_nodes)) >
        sqrt(.Machine$double.eps)
  ) {
    stop("nuisance_integrated_nodes must be one integer of at least 2.")
  }
  if (
    length(maxit) != 1L ||
      !is.finite(maxit) ||
      maxit < 1L ||
      abs(maxit - round(maxit)) > sqrt(.Machine$double.eps)
  ) {
    stop("nuisance_integrated_maxit must be one positive integer.")
  }

  instrument_formula <- stats::reformulate(x_cols, response = "Z")
  instrument_fit <- stats::glm(
    formula = instrument_formula,
    data = train,
    family = stats::binomial(),
    control = stats::glm.control(maxit = 100)
  )

  rule <- gauss_legendre_rule(
    n_nodes = as.integer(quadrature_nodes),
    lower = -2,
    upper = 2
  )
  integration_weights <- rule$weights * d_u(rule$nodes)
  integration_weights <- integration_weights / sum(integration_weights)
  treatment <- as.numeric(train$T_obs)
  check_binary(treatment, "T_obs")
  treatment_design <- cbind(
    `(Intercept)` = 1,
    as.matrix(train[, x_cols, drop = FALSE]),
    Z = as.numeric(train$Z)
  )
  storage.mode(treatment_design) <- "double"

  marginal_probability <- function(coefficients, design) {
    linear_predictor <- as.vector(design %*% coefficients)
    probability_matrix <- stats::plogis(
      outer(
        linear_predictor,
        confounder_strength * rule$nodes,
        "+"
      )
    )
    as.vector(probability_matrix %*% integration_weights)
  }

  objective <- function(coefficients) {
    probability <- clip_probability(
      marginal_probability(coefficients, treatment_design),
      1e-10
    )
    -sum(
      treatment * log(probability) +
        (1 - treatment) * log1p(-probability)
    )
  }

  gradient <- function(coefficients) {
    linear_predictor <- as.vector(treatment_design %*% coefficients)
    conditional_probability <- stats::plogis(
      outer(
        linear_predictor,
        confounder_strength * rule$nodes,
        "+"
      )
    )
    probability <- as.vector(
      conditional_probability %*% integration_weights
    )
    probability_derivative <- as.vector(
      (conditional_probability * (1 - conditional_probability)) %*%
        integration_weights
    )
    probability <- clip_probability(probability, 1e-10)
    score_by_row <-
      (probability - treatment) * probability_derivative /
      (probability * (1 - probability))
    as.vector(crossprod(treatment_design, score_by_row))
  }

  initial_formula <- stats::reformulate(
    c(x_cols, "Z"),
    response = "T_obs"
  )
  initial_fit <- suppressWarnings(stats::glm(
    formula = initial_formula,
    data = train,
    family = stats::binomial(),
    control = stats::glm.control(maxit = 100)
  ))
  initial_coefficients <- stats::coef(initial_fit)
  initial_coefficients[!is.finite(initial_coefficients)] <- 0
  if (length(initial_coefficients) != ncol(treatment_design)) {
    initial_coefficients <- rep(0, ncol(treatment_design))
  }
  names(initial_coefficients) <- colnames(treatment_design)

  treatment_fit <- stats::optim(
    par = initial_coefficients,
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(maxit = as.integer(maxit), reltol = 1e-10)
  )
  names(treatment_fit$par) <- colnames(treatment_design)

  make_treatment_design <- function(newdata, z_value) {
    design <- cbind(
      `(Intercept)` = 1,
      as.matrix(newdata[, x_cols, drop = FALSE]),
      Z = rep(z_value, nrow(newdata))
    )
    storage.mode(design) <- "double"
    design
  }

  predict_nuisance <- function(newdata) {
    p_z1 <- stats::predict(
      instrument_fit,
      newdata = newdata,
      type = "response"
    )
    pi0 <- marginal_probability(
      treatment_fit$par,
      make_treatment_design(newdata, 0)
    )
    pi1 <- marginal_probability(
      treatment_fit$par,
      make_treatment_design(newdata, 1)
    )
    list(
      e = clip_probability(1 - p_z1, positivity_epsilon),
      pi0 = clip_probability(pi0, positivity_epsilon),
      pi1 = clip_probability(pi1, positivity_epsilon)
    )
  }

  treatment_converged <- identical(treatment_fit$convergence, 0L)
  treatment_rank <- ncol(treatment_design)
  list(
    predict = predict_nuisance,
    instrument_fit = instrument_fit,
    treatment_fit = treatment_fit,
    treatment_fit_initial_glm = initial_fit,
    learner = "integrated_parametric",
    diagnostics = data.frame(
      component = c("e", "pi0", "pi1"),
      learner = "integrated_parametric",
      converged = c(
        isTRUE(instrument_fit$converged),
        treatment_converged,
        treatment_converged
      ),
      model_rank = c(
        instrument_fit$rank,
        treatment_rank,
        treatment_rank
      ),
      effective_df = c(
        instrument_fit$rank,
        treatment_rank,
        treatment_rank
      ),
      n_training = rep(nrow(train), 3L),
      stringsAsFactors = FALSE
    ),
    configuration = list(
      confounder_strength = confounder_strength,
      quadrature_nodes = as.integer(quadrature_nodes),
      optimization_convergence = treatment_fit$convergence,
      optimization_message = treatment_fit$message,
      estimated_coefficients = treatment_fit$par
    )
  )
}

fit_nuisance_models <- function(
    train,
    x_cols = "X",
    positivity_epsilon = 0.05,
    learner = c(
      "gam", "ranger", "glm_spline", "library_cv", "super_learner",
      "integrated_parametric"
    ),
    treatment_spline_df = 5L,
    gam_basis_dimension = NULL,
    gam_minimum_basis_dimension = 5L,
    gam_maximum_basis_dimension = 20L,
    gam_basis_growth_exponent = 1 / 3,
    gam_method = "REML",
    integrated_confounder_strength = 1,
    integrated_quadrature_nodes = 40L,
    integrated_maxit = 500L,
    library_candidates = c("glm_spline", "gam", "ranger"),
    library_cv_folds = 3L,
    library_selection_rule = c("one_se", "minimum"),
    library_ranger_num_trees = 500L,
    library_seed = 1L,
    super_learner_candidates = NULL,
    super_learner_cv_folds = 3L,
    super_learner_ranger_num_trees = 500L,
    super_learner_hal_max_degree = 2L
) {
  learner <- match.arg(learner)
  library_selection_rule <- match.arg(library_selection_rule)
  if (learner == "glm_spline") {
    return(fit_nuisance_models_glm(
      train = train,
      x_cols = x_cols,
      positivity_epsilon = positivity_epsilon,
      treatment_spline_df = treatment_spline_df
    ))
  }
  if (learner == "gam") {
    return(fit_nuisance_models_gam(
      train = train,
      x_cols = x_cols,
      positivity_epsilon = positivity_epsilon,
      basis_dimension = gam_basis_dimension,
      minimum_basis_dimension = gam_minimum_basis_dimension,
      maximum_basis_dimension = gam_maximum_basis_dimension,
      basis_growth_exponent = gam_basis_growth_exponent,
      method = gam_method
    ))
  }
  if (learner == "ranger") {
    return(fit_nuisance_models_ranger(
      train = train,
      x_cols = x_cols,
      positivity_epsilon = positivity_epsilon
    ))
  }
  if (learner == "library_cv") {
    return(fit_nuisance_models_library_cv(
      train = train,
      x_cols = x_cols,
      positivity_epsilon = positivity_epsilon,
      candidates = library_candidates,
      cv_folds = library_cv_folds,
      selection_rule = library_selection_rule,
      treatment_spline_df = treatment_spline_df,
      gam_basis_dimension = gam_basis_dimension,
      gam_minimum_basis_dimension = gam_minimum_basis_dimension,
      gam_maximum_basis_dimension = gam_maximum_basis_dimension,
      gam_basis_growth_exponent = gam_basis_growth_exponent,
      gam_method = gam_method,
      ranger_num_trees = library_ranger_num_trees,
      seed = library_seed
    ))
  }
  if (learner == "super_learner") {
    return(fit_nuisance_models_super_learner(
      train = train,
      x_cols = x_cols,
      positivity_epsilon = positivity_epsilon,
      candidates = super_learner_candidates,
      cv_folds = super_learner_cv_folds,
      ranger_num_trees = super_learner_ranger_num_trees,
      hal_max_degree = super_learner_hal_max_degree,
      seed = library_seed
    ))
  }
  fit_nuisance_models_integrated_parametric(
    train = train,
    x_cols = x_cols,
    positivity_epsilon = positivity_epsilon,
    confounder_strength = integrated_confounder_strength,
    quadrature_nodes = integrated_quadrature_nodes,
    maxit = integrated_maxit
  )
}

subset_nuisance <- function(nuisance, index) {
  list(
    e = nuisance$e[index],
    pi0 = nuisance$pi0[index],
    pi1 = nuisance$pi1[index]
  )
}

compute_riesz_weight <- function(
    nuisance,
    t_value
) {
  if (t_value == 1) {
    nuisance$e * nuisance$pi0 /
      ((1 - nuisance$e) * nuisance$pi1)
  } else {
    nuisance$e * (1 - nuisance$pi0) /
      ((1 - nuisance$e) * (1 - nuisance$pi1))
  }
}

# =============================================================================
# 6. Direct RKHS estimator of gamma_t
# =============================================================================

rkhs_pairwise_distance <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)

  A2 <- rowSums(A^2)
  B2 <- rowSums(B^2)

  distance_squared <- outer(A2, B2, "+") -
    2 * tcrossprod(A, B)

  distance_squared[distance_squared < 0] <- 0

  sqrt(distance_squared)
}

rkhs_kernel <- function(
    A, B, sigma,
    kernel_type = c("gaussian", "sobolev_matern32", "sobolev_matern52")
) {
  kernel_type <- match.arg(kernel_type)
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("sigma must be one positive finite number.")
  }
  scaled_distance <- rkhs_pairwise_distance(A, B) / sigma

  if (kernel_type == "gaussian") {
    return(exp(-scaled_distance^2 / 2))
  }
  if (kernel_type == "sobolev_matern32") {
    argument <- sqrt(3) * scaled_distance
    return((1 + argument) * exp(-argument))
  }
  argument <- sqrt(5) * scaled_distance
  (1 + argument + argument^2 / 3) * exp(-argument)
}

rbf_kernel <- function(A, B, sigma) {
  rkhs_kernel(A, B, sigma, kernel_type = "gaussian")
}

first_primes <- function(count) {
  count <- as.integer(count)
  if (length(count) != 1L || is.na(count) || count < 1L) {
    stop("count must be one positive integer.")
  }
  primes <- integer()
  candidate <- 2L
  while (length(primes) < count) {
    divisor <- primes[primes <= floor(sqrt(candidate))]
    is_prime <- !length(divisor) || all(candidate %% divisor != 0L)
    if (is_prime) primes <- c(primes, candidate)
    candidate <- candidate + 1L
  }
  primes
}

radical_inverse <- function(index, base) {
  index <- as.integer(index)
  base <- as.integer(base)
  output <- numeric(length(index))
  fraction <- 1 / base
  remaining <- index
  while (any(remaining > 0L)) {
    output <- output + fraction * (remaining %% base)
    remaining <- remaining %/% base
    fraction <- fraction / base
  }
  output
}

make_fixed_rbf_centers <- function(
    rank,
    dimension,
    lower = -2,
    upper = 2,
    skip = 20L,
    marginal = c("truncated_normal", "uniform")
) {
  rank <- as.integer(rank)
  dimension <- as.integer(dimension)
  skip <- as.integer(skip)
  marginal <- match.arg(marginal)
  if (
    length(rank) != 1L || is.na(rank) || rank < 1L ||
      length(dimension) != 1L || is.na(dimension) || dimension < 1L ||
      length(skip) != 1L || is.na(skip) || skip < 0L
  ) {
    stop("rank and dimension must be positive; skip must be nonnegative.")
  }
  lower <- rep_len(as.numeric(lower), dimension)
  upper <- rep_len(as.numeric(upper), dimension)
  if (
    any(!is.finite(lower)) || any(!is.finite(upper)) ||
      any(lower >= upper)
  ) {
    stop("Every fixed-RBF lower limit must be smaller than its upper limit.")
  }
  index <- skip + seq_len(rank)
  bases <- first_primes(dimension)
  unit <- vapply(
    bases,
    function(base) radical_inverse(index, base),
    numeric(rank)
  )
  if (dimension == 1L) unit <- matrix(unit, ncol = 1L)
  centers <- if (marginal == "uniform") {
    sweep(sweep(unit, 2L, upper - lower, "*"), 2L, lower, "+")
  } else {
    lower_probability <- stats::pnorm(lower)
    upper_probability <- stats::pnorm(upper)
    probability <- sweep(
      sweep(unit, 2L, upper_probability - lower_probability, "*"),
      2L,
      lower_probability,
      "+"
    )
    stats::qnorm(probability)
  }
  colnames(centers) <- paste0("X", seq_len(dimension))
  centers
}

prepare_rkhs_coordinates <- function(
    train_t,
    x_cols = "X",
    bandwidth_sample = 1000,
    seed = 1
) {
  non_numeric <- x_cols[
    !vapply(
      train_t[, x_cols, drop = FALSE],
      is.numeric,
      logical(1)
    )
  ]

  if (length(non_numeric) > 0) {
    stop(
      "RKHS coordinates must be numeric: ",
      paste(non_numeric, collapse = ", ")
    )
  }

  raw_coordinates <- cbind(
    Y = train_t$Y,
    as.matrix(
      train_t[, x_cols, drop = FALSE]
    )
  )

  center <- colMeans(raw_coordinates)

  scale_value <- apply(
    raw_coordinates,
    2,
    stats::sd
  )

  scale_value[
    !is.finite(scale_value) |
      scale_value < 1e-8
  ] <- 1

  transform_coordinates <- function(newdata) {
    coordinates <- cbind(
      Y = newdata$Y,
      as.matrix(
        newdata[, x_cols, drop = FALSE]
      )
    )

    sweep(
      sweep(
        coordinates,
        2,
        center,
        "-"
      ),
      2,
      scale_value,
      "/"
    )
  }

  transformed_train <- transform_coordinates(
    train_t
  )

  set.seed(seed)

  sample_index <- sample(
    seq_len(nrow(transformed_train)),
    min(
      bandwidth_sample,
      nrow(transformed_train)
    )
  )

  distances <- as.numeric(
    stats::dist(
      transformed_train[
        sample_index,
        ,
        drop = FALSE
      ]
    )
  )

  distances <- distances[
    is.finite(distances) &
      distances > 0
  ]

  sigma <- if (length(distances) == 0) {
    1
  } else {
    stats::median(distances)
  }

  sigma <- max(sigma, 1e-4)

  list(
    transformed_train = transformed_train,
    transform = transform_coordinates,
    center = center,
    scale = scale_value,
    sigma = sigma
  )
}

resolve_gamma_nystrom_rank <- function(rank, n_training) {
  if (
    length(n_training) != 1L ||
      !is.finite(n_training) ||
      n_training < 1L
  ) {
    stop("n_training must be a positive integer.")
  }
  n_training <- as.integer(n_training)
  if (is.function(rank)) {
    rank <- rank(n_training)
  }
  if (is.null(rank)) {
    # A fixed finite cap can leave non-vanishing approximation error.  Keep
    # the automatic rank increasing with the treatment-arm training size.
    rank <- ceiling(10 * sqrt(n_training))
  }
  if (
    length(rank) != 1L ||
      !is.finite(rank) ||
      rank < 1L ||
      abs(rank - round(rank)) > sqrt(.Machine$double.eps)
  ) {
    stop(
      paste0(
        "gamma_nystrom_rank must be NULL, one positive integer, or a ",
        "function returning one positive integer from n_training."
      )
    )
  }
  min(as.integer(round(rank)), n_training)
}

solve_gamma_exact_reduced_system <- function(
    kernel,
    baseline_indicator,
    right_hand_side,
    ridge,
    eigen_tolerance = 1e-12
) {
  kernel <- as.matrix(kernel)
  n_training <- nrow(kernel)
  if (
    n_training < 1L ||
      ncol(kernel) != n_training ||
      length(baseline_indicator) != n_training ||
      length(right_hand_side) != n_training
  ) {
    stop("The exact gamma system has incompatible dimensions.")
  }
  if (length(ridge) != 1L || !is.finite(ridge) || ridge <= 0) {
    stop("ridge must be one positive finite number.")
  }
  if (
    length(eigen_tolerance) != 1L ||
      !is.finite(eigen_tolerance) ||
      eigen_tolerance <= 0
  ) {
    stop("eigen_tolerance must be one positive finite number.")
  }

  baseline <- which(as.numeric(baseline_indicator) > 0.5)
  comparison <- setdiff(seq_len(n_training), baseline)
  alpha <- numeric(n_training)

  # In comparison-stratum rows diag(W0) is zero, so the reduced equation
  # gives ridge * alpha_1 = w_1 exactly.
  if (length(comparison) > 0L) {
    alpha[comparison] <- right_hand_side[comparison] / ridge
  }
  if (length(baseline) == 0L) {
    return(list(
      solution = alpha,
      method = "comparison_closed_form",
      minimum_eigenvalue = NA_real_,
      eigenvalue_cutoff = NA_real_
    ))
  }

  baseline_kernel <- kernel[baseline, baseline, drop = FALSE]
  baseline_kernel <- (baseline_kernel + t(baseline_kernel)) / 2
  diag(baseline_kernel) <- diag(baseline_kernel) + ridge
  baseline_rhs <- right_hand_side[baseline]
  if (length(comparison) > 0L) {
    baseline_rhs <- baseline_rhs - as.vector(
      kernel[baseline, comparison, drop = FALSE] %*%
        alpha[comparison]
    )
  }

  factor <- tryCatch(chol(baseline_kernel), error = function(error) NULL)
  if (!is.null(factor)) {
    alpha[baseline] <- backsolve(
      factor,
      forwardsolve(t(factor), baseline_rhs)
    )
    return(list(
      solution = alpha,
      method = "block_cholesky",
      minimum_eigenvalue = NA_real_,
      eigenvalue_cutoff = NA_real_
    ))
  }

  # Ridge regularization should make this block positive definite.  Retain a
  # reported symmetric fallback for floating-point failures rather than
  # silently applying a large arbitrary jitter.
  decomposition <- eigen(baseline_kernel, symmetric = TRUE)
  maximum_eigenvalue <- max(decomposition$values)
  cutoff <- eigen_tolerance * max(1, maximum_eigenvalue)
  stabilized_values <- pmax(decomposition$values, cutoff)
  alpha[baseline] <- as.vector(
    decomposition$vectors %*%
      (crossprod(decomposition$vectors, baseline_rhs) /
        stabilized_values)
  )
  list(
    solution = alpha,
    method = "block_eigenvalue_floor",
    minimum_eigenvalue = min(decomposition$values),
    eigenvalue_cutoff = cutoff
  )
}

prepare_gamma_nystrom_landmarks <- function(
    coordinate_spec,
    rank = NULL,
    seed = 1
) {
  n_training <- nrow(coordinate_spec$transformed_train)
  resolved_rank <- resolve_gamma_nystrom_rank(rank, n_training)
  had_random_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  if (had_random_seed) {
    saved_random_seed <- get(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  }
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", saved_random_seed, envir = .GlobalEnv)
    } else if (exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  landmark_indices <- if (resolved_rank == n_training) {
    seq_len(n_training)
  } else {
    sort(sample(seq_len(n_training), resolved_rank, replace = FALSE))
  }
  list(
    rank_requested = resolved_rank,
    landmark_indices = landmark_indices,
    transformed_landmarks = coordinate_spec$transformed_train[
      landmark_indices,
      ,
      drop = FALSE
    ]
  )
}

make_gamma_nystrom_features <- function(
    coordinate_spec,
    landmark_spec,
    bandwidth_multiplier = 1,
    kernel_type = c("gaussian", "sobolev_matern32", "sobolev_matern52"),
    eigen_tolerance = 1e-10,
    numerical_jitter = 1e-10
) {
  kernel_type <- match.arg(kernel_type)
  if (
    length(eigen_tolerance) != 1L ||
      !is.finite(eigen_tolerance) ||
      eigen_tolerance < 0
  ) {
    stop("gamma_nystrom_eigen_tolerance must be nonnegative and finite.")
  }
  if (
    length(numerical_jitter) != 1L ||
      !is.finite(numerical_jitter) ||
      numerical_jitter < 0
  ) {
    stop("gamma_nystrom_jitter must be nonnegative and finite.")
  }
  sigma <- coordinate_spec$sigma * bandwidth_multiplier
  landmarks <- landmark_spec$transformed_landmarks
  landmark_kernel <- rkhs_kernel(
    landmarks, landmarks, sigma, kernel_type = kernel_type
  )
  landmark_kernel <- (landmark_kernel + t(landmark_kernel)) / 2
  if (numerical_jitter > 0) {
    diag(landmark_kernel) <- diag(landmark_kernel) + numerical_jitter
  }
  decomposition <- eigen(landmark_kernel, symmetric = TRUE)
  maximum_eigenvalue <- max(decomposition$values)
  cutoff <- max(
    eigen_tolerance * max(maximum_eigenvalue, 0),
    0
  )
  retained <- which(decomposition$values > cutoff)
  if (length(retained) == 0L) {
    retained <- which.max(decomposition$values)
  }
  vectors <- decomposition$vectors[, retained, drop = FALSE]
  values <- decomposition$values[retained]
  if (any(!is.finite(values)) || any(values <= 0)) {
    stop("Nyström landmark kernel has no stable positive eigenspace.")
  }
  feature_map <- function(newdata) {
    kernel_to_landmarks <- rkhs_kernel(
      coordinate_spec$transform(newdata),
      landmarks,
      sigma,
      kernel_type = kernel_type
    )
    sweep(kernel_to_landmarks %*% vectors, 2L, sqrt(values), "/")
  }
  list(
    features_train = sweep(
      rkhs_kernel(
        coordinate_spec$transformed_train,
        landmarks,
        sigma,
        kernel_type = kernel_type
      ) %*% vectors,
      2L,
      sqrt(values),
      "/"
    ),
    predict_features = feature_map,
    rank_requested = landmark_spec$rank_requested,
    landmark_count = nrow(landmarks),
    effective_rank = length(values),
    min_retained_eigenvalue = min(values),
    max_eigenvalue = maximum_eigenvalue,
    landmark_indices = landmark_spec$landmark_indices,
    sigma = sigma,
    kernel_type = kernel_type
  )
}

fit_gamma_rkhs_reduced <- function(
    train_t,
    nuisance_t,
    t_value,
    x_cols = "X",
    lambda,
    bandwidth_multiplier = 1 / sqrt(2),
    coordinate_spec = NULL,
    bandwidth_sample = 1000,
    numerical_jitter = 1e-10,
    backend = c("exact", "nystrom"),
    nystrom_rank = NULL,
    nystrom_landmark_spec = NULL,
    nystrom_feature_spec = NULL,
    nystrom_eigen_tolerance = 1e-10,
    nystrom_jitter = 1e-10,
    kernel_type = c("gaussian", "sobolev_matern32", "sobolev_matern52"),
    seed = 1
) {
  backend <- match.arg(backend)
  kernel_type <- match.arg(kernel_type)

  if (
    length(bandwidth_multiplier) != 1 ||
      !is.finite(bandwidth_multiplier) ||
      bandwidth_multiplier <= 0
  ) {
    stop("bandwidth_multiplier must be one positive finite value.")
  }

  if (is.null(coordinate_spec)) {
    coordinate_spec <- prepare_rkhs_coordinates(
      train_t = train_t,
      x_cols = x_cols,
      bandwidth_sample = bandwidth_sample,
      seed = seed
    )
  }

  selected_sigma <- coordinate_spec$sigma *
    bandwidth_multiplier

  W0 <- as.numeric(train_t$Z == 0)

  w1 <- as.numeric(train_t$Z == 1) *
    compute_riesz_weight(
      nuisance = nuisance_t,
      t_value = t_value
    )

  n_t <- nrow(train_t)

  if (backend == "nystrom") {
    if (is.null(nystrom_feature_spec)) {
      if (is.null(nystrom_landmark_spec)) {
        nystrom_landmark_spec <- prepare_gamma_nystrom_landmarks(
          coordinate_spec = coordinate_spec,
          rank = nystrom_rank,
          seed = seed
        )
      }
      nystrom_feature_spec <- make_gamma_nystrom_features(
        coordinate_spec = coordinate_spec,
        landmark_spec = nystrom_landmark_spec,
        bandwidth_multiplier = bandwidth_multiplier,
        kernel_type = kernel_type,
        eigen_tolerance = nystrom_eigen_tolerance,
        numerical_jitter = nystrom_jitter
      )
    }

    features <- nystrom_feature_spec$features_train
    weighted_features <- sweep(features, 1L, W0, "*")
    system_matrix <- crossprod(features, weighted_features)
    diag(system_matrix) <- diag(system_matrix) +
      n_t * lambda + numerical_jitter
    right_hand_side <- as.vector(crossprod(features, w1))
    theta <- tryCatch(
      solve(system_matrix, right_hand_side),
      error = function(e) {
        qr.solve(system_matrix, right_hand_side, tol = 1e-10)
      }
    )

    return(list(
      predict_raw = function(newdata) {
        as.vector(
          nystrom_feature_spec$predict_features(newdata) %*% theta
        )
      },
      theta = theta,
      sigma = selected_sigma,
      base_sigma = coordinate_spec$sigma,
      bandwidth_multiplier = bandwidth_multiplier,
      coordinate_spec = coordinate_spec,
      backend = backend,
      kernel_type = kernel_type,
      nystrom_rank_requested =
        nystrom_feature_spec$rank_requested,
      nystrom_landmark_count =
        nystrom_feature_spec$landmark_count,
      nystrom_effective_rank =
        nystrom_feature_spec$effective_rank,
      nystrom_min_retained_eigenvalue =
        nystrom_feature_spec$min_retained_eigenvalue,
      nystrom_max_eigenvalue =
        nystrom_feature_spec$max_eigenvalue,
      nystrom_landmark_indices =
        nystrom_feature_spec$landmark_indices,
      nystrom_feature_spec = nystrom_feature_spec
    ))
  }

  K <- rkhs_kernel(
    coordinate_spec$transformed_train,
    coordinate_spec$transformed_train,
    selected_sigma,
    kernel_type = kernel_type
  )

  exact_solution <- solve_gamma_exact_reduced_system(
    kernel = K,
    baseline_indicator = W0,
    right_hand_side = w1,
    ridge = n_t * lambda + numerical_jitter
  )
  alpha <- exact_solution$solution

  list(
    predict_raw = function(newdata) {
      K_new <- rkhs_kernel(
        coordinate_spec$transform(newdata),
        coordinate_spec$transformed_train,
        selected_sigma,
        kernel_type = kernel_type
      )

      as.vector(K_new %*% alpha)
    },
    alpha = alpha,
    sigma = selected_sigma,
    base_sigma = coordinate_spec$sigma,
    bandwidth_multiplier = bandwidth_multiplier,
    coordinate_spec = coordinate_spec,
    backend = backend,
    kernel_type = kernel_type,
    nystrom_rank_requested = NA_integer_,
    nystrom_landmark_count = NA_integer_,
    nystrom_effective_rank = NA_integer_,
    nystrom_min_retained_eigenvalue = NA_real_,
    nystrom_max_eigenvalue = NA_real_,
    nystrom_landmark_indices = integer(),
    exact_solver_method = exact_solution$method,
    exact_solver_minimum_eigenvalue =
      exact_solution$minimum_eigenvalue,
    exact_solver_eigenvalue_cutoff =
      exact_solution$eigenvalue_cutoff
  )
}

rkhs_validation_risk <- function(
    validation_t,
    nuisance_validation_t,
    gamma_raw,
    t_value
) {
  w_tilde <- compute_riesz_weight(
    nuisance = nuisance_validation_t,
    t_value = t_value
  )

  mean(
    0.5 *
      as.numeric(validation_t$Z == 0) *
      gamma_raw^2 -
      as.numeric(validation_t$Z == 1) *
      w_tilde *
      gamma_raw
  )
}

make_rkhs_cv_plan <- function(
    train_t,
    x_cols = "X",
    cv_folds = 3,
    cv_repeats = 2L,
    bandwidth_sample = 1000,
    backend = c("exact", "nystrom"),
    nystrom_rank = NULL,
    seed = 1
) {
  backend <- match.arg(backend)
  if (
    length(cv_repeats) != 1L || !is.finite(cv_repeats) ||
      cv_repeats < 1L || cv_repeats != as.integer(cv_repeats)
  ) {
    stop("cv_repeats must be one positive integer.")
  }
  cv_repeats <- as.integer(cv_repeats)

  repeat_fold_id <- lapply(seq_len(cv_repeats), function(repeat_index) {
    make_z_stratified_folds(
      train_t,
      K = cv_folds,
      seed = seed + 100000L * (repeat_index - 1L)
    )
  })
  splits <- vector("list", cv_repeats * cv_folds)
  split_index <- 0L
  for (repeat_index in seq_len(cv_repeats)) {
    for (fold in seq_len(cv_folds)) {
      split_index <- split_index + 1L
      validation_index <- which(repeat_fold_id[[repeat_index]] == fold)
      splits[[split_index]] <- list(
        repeat_id = repeat_index,
        fold = fold,
        training_index = which(repeat_fold_id[[repeat_index]] != fold),
        validation_index = validation_index
      )
    }
  }

  coordinate_specs <- vector(
    "list",
    length(splits)
  )
  nystrom_landmark_specs <- vector(
    "list",
    length(splits)
  )

  for (split_index in seq_along(splits)) {
    split <- splits[[split_index]]
    repeat_offset <- 100000L * (split$repeat_id - 1L)
    training_index <- split$training_index

    coordinate_specs[[split_index]] <-
      prepare_rkhs_coordinates(
        train_t = train_t[
          training_index,
          ,
          drop = FALSE
        ],
        x_cols = x_cols,
        bandwidth_sample = bandwidth_sample,
        seed = seed + repeat_offset + 1000L * split$fold
      )

    if (backend == "nystrom") {
      nystrom_landmark_specs[[split_index]] <-
        prepare_gamma_nystrom_landmarks(
          coordinate_spec = coordinate_specs[[split_index]],
          rank = nystrom_rank,
          seed = seed + repeat_offset + 2000L * split$fold
        )
    }
  }

  final_coordinate_spec <-
    prepare_rkhs_coordinates(
      train_t = train_t,
      x_cols = x_cols,
      bandwidth_sample = bandwidth_sample,
      seed = seed + 50000
    )

  final_nystrom_landmark_spec <- if (backend == "nystrom") {
    prepare_gamma_nystrom_landmarks(
      coordinate_spec = final_coordinate_spec,
      rank = nystrom_rank,
      seed = seed + 60000
    )
  } else {
    NULL
  }

  list(
    # Keep the first partition under the historical name for backward
    # compatibility. New code should use splits and repeat_fold_id.
    fold_id = repeat_fold_id[[1L]],
    repeat_fold_id = repeat_fold_id,
    splits = splits,
    cv_folds = as.integer(cv_folds),
    cv_repeats = cv_repeats,
    coordinate_specs = coordinate_specs,
    final_coordinate_spec = final_coordinate_spec,
    backend = backend,
    nystrom_rank = nystrom_rank,
    nystrom_landmark_specs = nystrom_landmark_specs,
    final_nystrom_landmark_spec = final_nystrom_landmark_spec
  )
}

select_rkhs_lambda <- function(
    train_t,
    nuisance_t,
    t_value,
    x_cols = "X",
    lambda_grid = 10^seq(-5, -1, length.out = 9),
    bandwidth_multipliers = 1 / sqrt(2),
    cv_plan,
    selection_rule = c("minimum", "one_se"),
    bandwidth_sample = 1000,
    backend = c("exact", "nystrom"),
    nystrom_rank = NULL,
    nystrom_eigen_tolerance = 1e-10,
    nystrom_jitter = 1e-10,
    kernel_type = c("gaussian", "sobolev_matern32", "sobolev_matern52")
) {
  backend <- match.arg(backend)
  kernel_type <- match.arg(kernel_type)
  selection_rule <- match.arg(
    selection_rule
  )

  lambda_grid <- sort(
    unique(
      as.numeric(lambda_grid)
    )
  )

  if (
    length(lambda_grid) == 0 ||
      any(!is.finite(lambda_grid)) ||
      any(lambda_grid <= 0)
  ) {
    stop(
      "lambda_grid must contain positive finite values."
    )
  }

  bandwidth_multipliers <- sort(
    unique(
      as.numeric(bandwidth_multipliers)
    )
  )

  if (
    length(bandwidth_multipliers) == 0 ||
      any(!is.finite(bandwidth_multipliers)) ||
      any(bandwidth_multipliers <= 0)
  ) {
    stop(
      "bandwidth_multipliers must contain positive finite values."
    )
  }

  candidate_grid <- expand.grid(
    lambda = lambda_grid,
    bandwidth_multiplier = bandwidth_multipliers,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  splits <- if (!is.null(cv_plan$splits)) {
    cv_plan$splits
  } else {
    lapply(seq_len(max(cv_plan$fold_id)), function(fold) {
      list(
        repeat_id = 1L,
        fold = fold,
        training_index = which(cv_plan$fold_id != fold),
        validation_index = which(cv_plan$fold_id == fold)
      )
    })
  }
  cv_resamples <- length(splits)

  fold_risk <- matrix(
    NA_real_,
    nrow = nrow(candidate_grid),
    ncol = cv_resamples
  )
  fold_effective_rank <- matrix(
    NA_integer_,
    nrow = length(bandwidth_multipliers),
    ncol = cv_resamples
  )

  for (resample in seq_along(splits)) {
    validation_index <- splits[[resample]]$validation_index
    training_index <- splits[[resample]]$training_index

    cv_train <- train_t[
      training_index,
      ,
      drop = FALSE
    ]

    cv_validation <- train_t[
      validation_index,
      ,
      drop = FALSE
    ]

    nuisance_cv_train <- subset_nuisance(
      nuisance_t,
      training_index
    )

    nuisance_cv_validation <- subset_nuisance(
      nuisance_t,
      validation_index
    )

    nystrom_feature_specs <- if (backend == "nystrom") {
      lapply(
        bandwidth_multipliers,
        function(multiplier) {
          make_gamma_nystrom_features(
            coordinate_spec = cv_plan$coordinate_specs[[resample]],
            landmark_spec = cv_plan$nystrom_landmark_specs[[resample]],
            bandwidth_multiplier = multiplier,
            kernel_type = kernel_type,
            eigen_tolerance = nystrom_eigen_tolerance,
            numerical_jitter = nystrom_jitter
          )
        }
      )
    } else {
      NULL
    }

    if (backend == "nystrom") {
      fold_effective_rank[, resample] <- vapply(
        nystrom_feature_specs,
        function(specification) specification$effective_rank,
        integer(1)
      )
    }

    for (candidate in seq_len(nrow(candidate_grid))) {
      bandwidth_index <- match(
        candidate_grid$bandwidth_multiplier[candidate],
        bandwidth_multipliers
      )
      fit <- fit_gamma_rkhs_reduced(
        train_t = cv_train,
        nuisance_t = nuisance_cv_train,
        t_value = t_value,
        x_cols = x_cols,
        lambda = candidate_grid$lambda[candidate],
        bandwidth_multiplier =
          candidate_grid$bandwidth_multiplier[candidate],
        coordinate_spec =
          cv_plan$coordinate_specs[[resample]],
        bandwidth_sample = bandwidth_sample,
        backend = backend,
        nystrom_rank = nystrom_rank,
        nystrom_landmark_spec =
          cv_plan$nystrom_landmark_specs[[resample]],
        nystrom_feature_spec = if (backend == "nystrom") {
          nystrom_feature_specs[[bandwidth_index]]
        } else {
          NULL
        },
        nystrom_eigen_tolerance = nystrom_eigen_tolerance,
        nystrom_jitter = nystrom_jitter,
        kernel_type = kernel_type
      )

      raw_validation <- fit$predict_raw(
        cv_validation
      )

      fold_risk[
        candidate,
        resample
      ] <- rkhs_validation_risk(
        validation_t = cv_validation,
        nuisance_validation_t =
          nuisance_cv_validation,
        gamma_raw = raw_validation,
        t_value = t_value
      )
    }
  }

  mean_risk <- rowMeans(
    fold_risk
  )

  standard_error <- apply(
    fold_risk,
    1,
    stats::sd
  ) / sqrt(cv_resamples)

  minimum_index <- which.min(
    mean_risk
  )

  selected_index <- if (
    selection_rule == "minimum"
  ) {
    minimum_index
  } else {
    threshold <- mean_risk[
      minimum_index
    ] + standard_error[
      minimum_index
    ]

    acceptable <- which(
      mean_risk <= threshold
    )

    # Among candidates within one standard error, prefer stronger ridge
    # regularization and then the smoother (larger) kernel bandwidth.
    ordering <- order(
      candidate_grid$lambda[acceptable],
      candidate_grid$bandwidth_multiplier[acceptable],
      decreasing = TRUE
    )
    acceptable[ordering[1]]
  }

  list(
    selected_lambda =
      candidate_grid$lambda[selected_index],
    selected_bandwidth_multiplier =
      candidate_grid$bandwidth_multiplier[selected_index],
    cv_table = data.frame(
      lambda = candidate_grid$lambda,
      bandwidth_multiplier =
        candidate_grid$bandwidth_multiplier,
      rkhs_backend = backend,
      rkhs_kernel = kernel_type,
      nystrom_effective_rank_min = if (backend == "nystrom") {
        apply(fold_effective_rank, 1L, min)[
          match(
            candidate_grid$bandwidth_multiplier,
            bandwidth_multipliers
          )
        ]
      } else {
        NA_integer_
      },
      nystrom_effective_rank_max = if (backend == "nystrom") {
        apply(fold_effective_rank, 1L, max)[
          match(
            candidate_grid$bandwidth_multiplier,
            bandwidth_multipliers
          )
        ]
      } else {
        NA_integer_
      },
      cv_folds = if (is.null(cv_plan$cv_folds)) {
        max(cv_plan$fold_id)
      } else {
        cv_plan$cv_folds
      },
      cv_repeats = if (is.null(cv_plan$cv_repeats)) {
        1L
      } else {
        cv_plan$cv_repeats
      },
      cv_resamples = cv_resamples,
      mean_risk = mean_risk,
      standard_error = standard_error,
      selected =
        seq_len(nrow(candidate_grid)) ==
        selected_index
    ),
    cv_resample_table = do.call(
      rbind,
      lapply(seq_len(nrow(candidate_grid)), function(candidate) {
        data.frame(
          lambda = candidate_grid$lambda[candidate],
          bandwidth_multiplier =
            candidate_grid$bandwidth_multiplier[candidate],
          cv_repeat = vapply(splits, `[[`, integer(1L), "repeat_id"),
          fold = vapply(splits, `[[`, integer(1L), "fold"),
          risk = as.numeric(fold_risk[candidate, ]),
          selected = candidate == selected_index,
          stringsAsFactors = FALSE
        )
      })
    ),
    fold_risk = fold_risk
  )
}

default_gamma_lambda_grid <- function(
    n_training,
    treatment = NULL,
    constants = 10^seq(-1.5, 1, by = 0.5)
) {
  if (
    length(n_training) != 1L ||
      !is.finite(n_training) ||
      n_training < 1L
  ) {
    stop("n_training must be a positive integer.")
  }
  constants <- sort(unique(as.numeric(constants)))
  if (
    length(constants) == 0L ||
      any(!is.finite(constants)) ||
      any(constants <= 0)
  ) {
    stop("constants must contain positive finite values.")
  }
  constants / sqrt(n_training)
}

resolve_gamma_lambda_grid <- function(
    lambda_grid,
    n_training,
    treatment
) {
  if (
    length(n_training) != 1L ||
      !is.finite(n_training) ||
      n_training < 1L
  ) {
    stop("n_training must be a positive integer.")
  }

  if (is.function(lambda_grid)) {
    lambda_grid <- lambda_grid(
      n_training = as.integer(n_training),
      treatment = treatment
    )
  }

  lambda_grid <- sort(unique(as.numeric(lambda_grid)))
  if (
    length(lambda_grid) == 0L ||
      any(!is.finite(lambda_grid)) ||
      any(lambda_grid <= 0)
  ) {
    stop(
      paste0(
        "The resolved gamma lambda grid must contain positive finite ",
        "values."
      )
    )
  }
  lambda_grid
}

fit_gamma_rkhs <- function(
    train,
    nuisance_train,
    t_value,
    x_cols = "X",
    lambda_grid = default_gamma_lambda_grid,
    bandwidth_multipliers = 1 / sqrt(2),
    cv_folds = 3,
    cv_repeats = 2L,
    selection_rule = c("minimum", "one_se"),
    gamma_floor = 1e-6,
    gamma_cap = 100,
    bandwidth_sample = 1000,
    cv_plan = NULL,
    backend = c("exact", "nystrom"),
    nystrom_rank = NULL,
    nystrom_eigen_tolerance = 1e-10,
    nystrom_jitter = 1e-10,
    kernel_type = c("gaussian", "sobolev_matern32", "sobolev_matern52"),
    seed = 1
) {
  backend <- match.arg(backend)
  kernel_type <- match.arg(kernel_type)
  selection_rule <- match.arg(
    selection_rule
  )

  treatment_index <- which(
    train$T_obs == t_value
  )

  train_t <- train[
    treatment_index,
    ,
    drop = FALSE
  ]

  nuisance_t <- subset_nuisance(
    nuisance_train,
    treatment_index
  )

  lambda_grid <- resolve_gamma_lambda_grid(
    lambda_grid = lambda_grid,
    n_training = nrow(train_t),
    treatment = t_value
  )

  if (is.null(cv_plan)) {
    cv_plan <- make_rkhs_cv_plan(
      train_t = train_t,
      x_cols = x_cols,
      cv_folds = cv_folds,
      cv_repeats = cv_repeats,
      bandwidth_sample = bandwidth_sample,
      backend = backend,
      nystrom_rank = nystrom_rank,
      seed = seed
    )
  } else if (!identical(cv_plan$backend, backend)) {
    stop("cv_plan backend does not match the requested gamma backend.")
  }

  tuning <- select_rkhs_lambda(
    train_t = train_t,
    nuisance_t = nuisance_t,
    t_value = t_value,
    x_cols = x_cols,
    lambda_grid = lambda_grid,
    bandwidth_multipliers = bandwidth_multipliers,
    cv_plan = cv_plan,
    selection_rule = selection_rule,
    bandwidth_sample = bandwidth_sample,
    backend = backend,
    nystrom_rank = nystrom_rank,
    nystrom_eigen_tolerance = nystrom_eigen_tolerance,
    nystrom_jitter = nystrom_jitter,
    kernel_type = kernel_type
  )

  final_nystrom_feature_spec <- if (backend == "nystrom") {
    make_gamma_nystrom_features(
      coordinate_spec = cv_plan$final_coordinate_spec,
      landmark_spec = cv_plan$final_nystrom_landmark_spec,
      bandwidth_multiplier = tuning$selected_bandwidth_multiplier,
      kernel_type = kernel_type,
      eigen_tolerance = nystrom_eigen_tolerance,
      numerical_jitter = nystrom_jitter
    )
  } else {
    NULL
  }

  final_fit <- fit_gamma_rkhs_reduced(
    train_t = train_t,
    nuisance_t = nuisance_t,
    t_value = t_value,
    x_cols = x_cols,
    lambda = tuning$selected_lambda,
    bandwidth_multiplier =
      tuning$selected_bandwidth_multiplier,
    coordinate_spec =
      cv_plan$final_coordinate_spec,
    bandwidth_sample = bandwidth_sample,
    backend = backend,
    nystrom_rank = nystrom_rank,
    nystrom_landmark_spec =
      cv_plan$final_nystrom_landmark_spec,
    nystrom_feature_spec = final_nystrom_feature_spec,
    nystrom_eigen_tolerance = nystrom_eigen_tolerance,
    nystrom_jitter = nystrom_jitter,
    kernel_type = kernel_type
  )

  raw_train <- final_fit$predict_raw(
    train_t
  )

  positive_train <- clip_positive(
    raw_train,
    floor = gamma_floor,
    cap = gamma_cap
  )

  normalization_constant <- mean(
    positive_train[
      train_t$Z == 0
    ]
  )

  normalization_constant <- max(
    normalization_constant,
    gamma_floor
  )

  list(
    predict = function(newdata) {
      raw <- final_fit$predict_raw(
        newdata
      )

      positive <- clip_positive(
        raw,
        floor = gamma_floor,
        cap = gamma_cap
      )

      normalized <- clip_positive(
        positive /
          normalization_constant,
        floor = gamma_floor,
        cap = gamma_cap
      )

      list(
        raw = raw,
        positive = positive,
        normalized = normalized
      )
    },
    selected_lambda =
      tuning$selected_lambda,
    selected_bandwidth_multiplier =
      tuning$selected_bandwidth_multiplier,
    base_bandwidth = final_fit$base_sigma,
    selected_bandwidth = final_fit$sigma,
    cv_table = tuning$cv_table,
    cv_resample_table = tuning$cv_resample_table,
    normalization_constant =
      normalization_constant,
    training_negative_rate =
      mean(raw_train <= 0),
    rkhs_backend = final_fit$backend,
    rkhs_kernel = final_fit$kernel_type,
    nystrom_rank_requested =
      final_fit$nystrom_rank_requested,
    nystrom_landmark_count =
      final_fit$nystrom_landmark_count,
    nystrom_effective_rank =
      final_fit$nystrom_effective_rank,
    nystrom_min_retained_eigenvalue =
      final_fit$nystrom_min_retained_eigenvalue,
    nystrom_max_eigenvalue =
      final_fit$nystrom_max_eigenvalue,
    nystrom_landmark_indices =
      final_fit$nystrom_landmark_indices,
    exact_solver_method = final_fit$exact_solver_method,
    exact_solver_minimum_eigenvalue =
      final_fit$exact_solver_minimum_eigenvalue,
    exact_solver_eigenvalue_cutoff =
      final_fit$exact_solver_eigenvalue_cutoff,
    cv_plan = cv_plan
  )
}

make_cached_gamma_rkhs_provider <- function(
    lambda_grid = default_gamma_lambda_grid,
    bandwidth_multipliers = 1 / sqrt(2),
    cv_folds = 3L,
    cv_repeats = 2L,
    selection_rule = c("minimum", "one_se"),
    gamma_floor = 1e-6,
    gamma_cap = 100,
    bandwidth_sample = 1000L,
    backend = c("exact", "nystrom"),
    nystrom_rank = NULL,
    nystrom_eigen_tolerance = 1e-10,
    nystrom_jitter = 1e-10,
    kernel_type = c("gaussian", "sobolev_matern32", "sobolev_matern52")
) {
  selection_rule <- match.arg(selection_rule)
  backend <- match.arg(backend)
  kernel_type <- match.arg(kernel_type)
  cache <- new.env(parent = emptyenv())

  training_signature <- function(train, x_cols) {
    list(
      row_id = if (".row_id" %in% names(train)) train$.row_id else {
        seq_len(nrow(train))
      },
      Y = train$Y,
      T_obs = train$T_obs,
      Z = train$Z,
      X = as.matrix(train[, x_cols, drop = FALSE])
    )
  }

  provider <- function(
      fold,
      treatment,
      train,
      nuisance_train,
      x_cols,
      seed,
      ...
  ) {
    key <- paste0("fold", fold, "_t", treatment)
    signature <- training_signature(train, x_cols)
    if (exists(key, envir = cache, inherits = FALSE)) {
      entry <- get(key, envir = cache, inherits = FALSE)
      if (!identical(entry$training_signature, signature)) {
        stop(
          "A cached gamma fit was requested for a different training sample."
        )
      }
      entry$cache_hits <- entry$cache_hits + 1L
      assign(key, entry, envir = cache)
      return(entry$fit)
    }

    fit <- fit_gamma_rkhs(
      train = train,
      nuisance_train = nuisance_train,
      t_value = treatment,
      x_cols = x_cols,
      lambda_grid = lambda_grid,
      bandwidth_multipliers = bandwidth_multipliers,
      cv_folds = cv_folds,
      cv_repeats = cv_repeats,
      selection_rule = selection_rule,
      gamma_floor = gamma_floor,
      gamma_cap = gamma_cap,
      bandwidth_sample = bandwidth_sample,
      backend = backend,
      nystrom_rank = nystrom_rank,
      nystrom_eigen_tolerance = nystrom_eigen_tolerance,
      nystrom_jitter = nystrom_jitter,
      kernel_type = kernel_type,
      seed = seed
    )
    assign(
      key,
      list(
        fit = fit,
        training_signature = signature,
        cache_hits = 0L,
        fold = fold,
        treatment = treatment
      ),
      envir = cache
    )
    fit
  }
  attr(provider, "gamma_cache") <- cache
  class(provider) <- c("cached_gamma_rkhs_provider", class(provider))
  provider
}

compact_gamma_nystrom_fit_for_checkpoint <- function(fit) {
  if (!is.list(fit) || !identical(fit$rkhs_backend, "nystrom") ||
      !is.function(fit$predict)) {
    stop("Only a fitted Nyström gamma object can be compacted.")
  }
  fit_environment <- environment(fit$predict)
  required_fit_bindings <- c(
    "final_fit", "gamma_floor", "gamma_cap", "normalization_constant",
    "x_cols"
  )
  if (!all(vapply(
    required_fit_bindings,
    exists,
    logical(1),
    envir = fit_environment,
    inherits = TRUE
  ))) {
    stop("The Nyström gamma prediction closure has an unknown structure.")
  }
  final_fit <- get("final_fit", envir = fit_environment, inherits = TRUE)
  feature_spec <- final_fit$nystrom_feature_spec
  if (!is.list(feature_spec) || !is.function(feature_spec$predict_features)) {
    stop("The Nyström gamma fit is missing its feature specification.")
  }
  feature_environment <- environment(feature_spec$predict_features)
  required_feature_bindings <- c("landmarks", "vectors", "values")
  if (!all(vapply(
    required_feature_bindings,
    exists,
    logical(1),
    envir = feature_environment,
    inherits = TRUE
  ))) {
    stop("The Nyström feature closure has an unknown structure.")
  }
  coordinate_spec <- final_fit$coordinate_spec
  required_coordinate_fields <- c("center", "scale")
  if (!is.list(coordinate_spec) ||
      !all(required_coordinate_fields %in% names(coordinate_spec))) {
    stop("The Nyström gamma fit is missing its coordinate transformation.")
  }

  prediction_environment <- new.env(parent = environment(rkhs_kernel))
  bindings <- list(
    x_cols = get("x_cols", envir = fit_environment, inherits = TRUE),
    center = coordinate_spec$center,
    scale_value = coordinate_spec$scale,
    landmarks = get("landmarks", envir = feature_environment, inherits = TRUE),
    vectors = get("vectors", envir = feature_environment, inherits = TRUE),
    values = get("values", envir = feature_environment, inherits = TRUE),
    sigma = feature_spec$sigma,
    kernel_type = feature_spec$kernel_type,
    theta = final_fit$theta,
    gamma_floor = get("gamma_floor", envir = fit_environment, inherits = TRUE),
    gamma_cap = get("gamma_cap", envir = fit_environment, inherits = TRUE),
    normalization_constant = get(
      "normalization_constant", envir = fit_environment, inherits = TRUE
    )
  )
  list2env(bindings, envir = prediction_environment)
  compact_predict <- function(newdata) {
    coordinates <- cbind(
      Y = newdata$Y,
      as.matrix(newdata[, x_cols, drop = FALSE])
    )
    transformed <- sweep(
      sweep(coordinates, 2L, center, "-"),
      2L, scale_value, "/"
    )
    kernel_to_landmarks <- rkhs_kernel(
      transformed, landmarks, sigma, kernel_type = kernel_type
    )
    features <- sweep(
      kernel_to_landmarks %*% vectors,
      2L, sqrt(values), "/"
    )
    raw <- as.vector(features %*% theta)
    positive <- clip_positive(raw, floor = gamma_floor, cap = gamma_cap)
    normalized <- clip_positive(
      positive / normalization_constant,
      floor = gamma_floor,
      cap = gamma_cap
    )
    list(raw = raw, positive = positive, normalized = normalized)
  }
  environment(compact_predict) <- prediction_environment
  compact <- fit
  compact$predict <- compact_predict
  compact$cv_plan <- NULL
  compact$checkpoint_format <- "minimal_nystrom_prediction_v1"
  compact
}

collect_cached_gamma_rkhs_fits <- function(provider) {
  if (!inherits(provider, "cached_gamma_rkhs_provider")) {
    stop("provider must be returned by make_cached_gamma_rkhs_provider().")
  }
  cache <- attr(provider, "gamma_cache")
  keys <- sort(ls(envir = cache, all.names = TRUE))
  lapply(keys, function(key) get(key, envir = cache, inherits = FALSE))
}


# =============================================================================
# 7. Fixed-lambda RKHS fit used for fast bound debugging
# =============================================================================

fit_gamma_rkhs_fixed <- function(
    train,
    nuisance_train,
    t_value,
    x_cols = "X",
    lambda,
    bandwidth_multiplier = 1 / sqrt(2),
    gamma_floor = 1e-6,
    gamma_cap = 100,
    bandwidth_sample = 1000,
    backend = c("exact", "nystrom"),
    nystrom_rank = NULL,
    nystrom_eigen_tolerance = 1e-10,
    nystrom_jitter = 1e-10,
    kernel_type = c("gaussian", "sobolev_matern32", "sobolev_matern52"),
    seed = 1
) {
  backend <- match.arg(backend)
  kernel_type <- match.arg(kernel_type)
  treatment_index <- which(train$T_obs == t_value)

  train_t <- train[
    treatment_index,
    ,
    drop = FALSE
  ]

  nuisance_t <- subset_nuisance(
    nuisance_train,
    treatment_index
  )

  coordinate_spec <- prepare_rkhs_coordinates(
    train_t = train_t,
    x_cols = x_cols,
    bandwidth_sample = bandwidth_sample,
    seed = seed
  )

  nystrom_landmark_spec <- if (backend == "nystrom") {
    prepare_gamma_nystrom_landmarks(
      coordinate_spec = coordinate_spec,
      rank = nystrom_rank,
      seed = seed + 10000
    )
  } else {
    NULL
  }

  fit <- fit_gamma_rkhs_reduced(
    train_t = train_t,
    nuisance_t = nuisance_t,
    t_value = t_value,
    x_cols = x_cols,
    lambda = lambda,
    bandwidth_multiplier = bandwidth_multiplier,
    coordinate_spec = coordinate_spec,
    bandwidth_sample = bandwidth_sample,
    backend = backend,
    nystrom_rank = nystrom_rank,
    nystrom_landmark_spec = nystrom_landmark_spec,
    nystrom_eigen_tolerance = nystrom_eigen_tolerance,
    nystrom_jitter = nystrom_jitter,
    kernel_type = kernel_type,
    seed = seed
  )

  raw_train <- fit$predict_raw(train_t)

  positive_train <- clip_positive(
    raw_train,
    floor = gamma_floor,
    cap = gamma_cap
  )

  normalization_constant <- mean(
    positive_train[train_t$Z == 0]
  )

  normalization_constant <- max(
    normalization_constant,
    gamma_floor
  )

  list(
    predict = function(newdata) {
      raw <- fit$predict_raw(newdata)

      positive <- clip_positive(
        raw,
        floor = gamma_floor,
        cap = gamma_cap
      )

      normalized <- clip_positive(
        positive / normalization_constant,
        floor = gamma_floor,
        cap = gamma_cap
      )

      list(
        raw = raw,
        positive = positive,
        normalized = normalized
      )
    },
    selected_lambda = lambda,
    selected_bandwidth_multiplier = bandwidth_multiplier,
    base_bandwidth = fit$base_sigma,
    selected_bandwidth = fit$sigma,
    cv_table = NULL,
    normalization_constant = normalization_constant,
    training_negative_rate = mean(raw_train <= 0),
    sigma = fit$sigma,
    rkhs_backend = fit$backend,
    rkhs_kernel = fit$kernel_type,
    nystrom_rank_requested = fit$nystrom_rank_requested,
    nystrom_landmark_count = fit$nystrom_landmark_count,
    nystrom_effective_rank = fit$nystrom_effective_rank,
    nystrom_min_retained_eigenvalue =
      fit$nystrom_min_retained_eigenvalue,
    nystrom_max_eigenvalue = fit$nystrom_max_eigenvalue,
    nystrom_landmark_indices = fit$nystrom_landmark_indices,
    exact_solver_method = fit$exact_solver_method,
    exact_solver_minimum_eigenvalue =
      fit$exact_solver_minimum_eigenvalue,
    exact_solver_eigenvalue_cutoff =
      fit$exact_solver_eigenvalue_cutoff
  )
}

# =============================================================================
# 8. Spline regressions and check-loss optimization
# =============================================================================

make_spline_design <- function(
    data,
    x_cols = "X",
    spline_df = 5
) {
  non_numeric <- x_cols[
    !vapply(
      data[, x_cols, drop = FALSE],
      is.numeric,
      logical(1)
    )
  ]

  if (length(non_numeric) > 0) {
    stop(
      "All x_cols must be numeric for the spline learners: ",
      paste(non_numeric, collapse = ", ")
    )
  }

  basis_specs <- vector(
    "list",
    length(x_cols)
  )

  names(basis_specs) <- x_cols

  training_parts <- list(
    intercept = rep(
      1,
      nrow(data)
    )
  )

  for (name in x_cols) {
    unique_values <- sort(unique(data[[name]][is.finite(data[[name]])]))
    if (length(unique_values) < 2L) {
      warning(
        "Dropping constant spline predictor in this training subset: ",
        name,
        call. = FALSE
      )
      basis_specs[[name]] <- list(
        type = "constant",
        training_value = if (length(unique_values)) unique_values[1L] else NA_real_
      )
      next
    }
    if (length(unique_values) <= spline_df + 1L) {
      # Binary indicators and other low-cardinality numeric predictors should
      # enter as ordinary linear columns. Forcing natural splines through them
      # creates repeated boundary knots and is not a meaningful expansion.
      basis <- matrix(data[[name]], ncol = 1L)
      basis_specs[[name]] <- list(type = "linear")
    } else {
      basis <- splines::ns(
        data[[name]],
        df = spline_df,
        intercept = FALSE
      )

      basis_specs[[name]] <- list(
        type = "spline",
        knots = attr(
          basis,
          "knots"
        ),
        boundary_knots = attr(
          basis,
          "Boundary.knots"
        )
      )
    }

    colnames(basis) <- paste0(
      name,
      "_s",
      seq_len(ncol(basis))
    )

    training_parts[[name]] <- basis
  }

  training_matrix <- do.call(
    cbind,
    training_parts
  )

  storage.mode(training_matrix) <-
    "double"

  predict_matrix <- function(newdata) {
    new_parts <- list(
      intercept = rep(
        1,
        nrow(newdata)
      )
    )

    for (name in x_cols) {
      spec <- basis_specs[[name]]
      if (identical(spec$type, "constant")) next
      basis <- if (identical(spec$type, "linear")) {
        matrix(newdata[[name]], ncol = 1L)
      } else {
        splines::ns(
          newdata[[name]],
          knots = spec$knots,
          Boundary.knots =
            spec$boundary_knots,
          intercept = FALSE
        )
      }

      colnames(basis) <- paste0(
        name,
        "_s",
        seq_len(ncol(basis))
      )

      new_parts[[name]] <- basis
    }

    output <- do.call(
      cbind,
      new_parts
    )

    storage.mode(output) <-
      "double"

    output
  }

  list(
    matrix = training_matrix,
    predict_matrix = predict_matrix,
    basis_specs = basis_specs
  )
}

fit_spline_mean <- function(
    data,
    outcome,
    x_cols = "X",
    spline_df = 5
) {
  design <- make_spline_design(
    data = data,
    x_cols = x_cols,
    spline_df = spline_df
  )

  outcome_value <- data[[outcome]]

  fit <- stats::lm.fit(
    x = design$matrix,
    y = outcome_value
  )

  coefficients <- fit$coefficients
  coefficients[!is.finite(coefficients)] <- 0

  list(
    predict = function(newdata) {
      as.vector(
        design$predict_matrix(newdata) %*%
          coefficients
      )
    },
    coefficients = coefficients,
    rank = fit$rank
  )
}

fit_spline_means <- function(
    data,
    outcomes,
    x_cols = "X",
    spline_df = 5
) {
  if (length(outcomes) == 0L) {
    stop("outcomes must contain at least one column name.")
  }

  missing_outcomes <- setdiff(outcomes, names(data))

  if (length(missing_outcomes) > 0L) {
    stop(
      "Missing spline outcomes: ",
      paste(missing_outcomes, collapse = ", ")
    )
  }

  design <- make_spline_design(
    data = data,
    x_cols = x_cols,
    spline_df = spline_df
  )

  outcome_matrix <- as.matrix(
    data[, outcomes, drop = FALSE]
  )

  storage.mode(outcome_matrix) <- "double"

  fit <- stats::lm.fit(
    x = design$matrix,
    y = outcome_matrix
  )

  coefficients <- fit$coefficients

  if (is.null(dim(coefficients))) {
    coefficients <- matrix(
      coefficients,
      ncol = 1L,
      dimnames = list(NULL, outcomes)
    )
  }

  coefficients[!is.finite(coefficients)] <- 0

  list(
    predict = function(newdata) {
      output <- design$predict_matrix(newdata) %*%
        coefficients

      output <- as.matrix(output)
      colnames(output) <- outcomes
      output
    },
    coefficients = coefficients,
    rank = fit$rank,
    outcomes = outcomes
  )
}

default_eif_conditional_control <- function() {
  list(
    spline_df_grid = c(3L, 5L, 7L, 9L),
    cv_folds = 3L,
    selection_rule = "one_se",
    gam_k = 6L,
    gam_method = "REML",
    ranger_num_trees = 500L,
    ranger_min_node_size = NULL,
    library_candidates = c("linear", "gam", "ranger"),
    library_selection_rule = "one_se"
  )
}

resolve_eif_conditional_control <- function(control = NULL) {
  output <- default_eif_conditional_control()
  if (!is.null(control)) {
    unknown <- setdiff(names(control), names(output))
    if (length(unknown) > 0L) {
      stop(
        "Unknown eif_conditional_control entries: ",
        paste(unknown, collapse = ", ")
      )
    }
    output[names(control)] <- control
  }
  output$spline_df_grid <- sort(unique(as.integer(output$spline_df_grid)))
  output$cv_folds <- as.integer(output$cv_folds)
  output$selection_rule <- match.arg(
    output$selection_rule,
    c("one_se", "minimum")
  )
  output$gam_k <- as.integer(output$gam_k)
  output$ranger_num_trees <- as.integer(output$ranger_num_trees)
  output$library_candidates <- unique(as.character(
    output$library_candidates
  ))
  output$library_selection_rule <- match.arg(
    output$library_selection_rule,
    c("one_se", "minimum")
  )
  if (!is.null(output$ranger_min_node_size)) {
    output$ranger_min_node_size <- as.integer(
      output$ranger_min_node_size
    )
  }
  if (
    length(output$spline_df_grid) == 0L ||
      anyNA(output$spline_df_grid) ||
      any(output$spline_df_grid < 1L) ||
      length(output$cv_folds) != 1L ||
      is.na(output$cv_folds) ||
      output$cv_folds < 2L ||
      length(output$gam_k) != 1L ||
      is.na(output$gam_k) ||
      output$gam_k < 3L ||
      length(output$gam_method) != 1L ||
      !nzchar(output$gam_method) ||
      length(output$ranger_num_trees) != 1L ||
      is.na(output$ranger_num_trees) ||
      output$ranger_num_trees < 1L ||
      length(output$library_candidates) < 1L ||
      any(!output$library_candidates %in% c("linear", "gam", "ranger")) ||
      (!is.null(output$ranger_min_node_size) && (
        length(output$ranger_min_node_size) != 1L ||
          is.na(output$ranger_min_node_size) ||
          output$ranger_min_node_size < 1L
      ))
  ) {
    stop("Invalid eif_conditional_control values.")
  }
  output
}

fit_spline_means_cv <- function(
    data,
    outcomes,
    x_cols = "X",
    spline_df_grid = c(3L, 5L, 7L, 9L),
    cv_folds = 3L,
    selection_rule = c("one_se", "minimum"),
    seed = 1L
) {
  selection_rule <- match.arg(selection_rule)
  spline_df_grid <- sort(unique(as.integer(spline_df_grid)))
  if (length(spline_df_grid) == 0L || any(spline_df_grid < 1L)) {
    stop("spline_df_grid must contain positive integers.")
  }
  cv_folds <- as.integer(cv_folds)
  if (cv_folds < 2L || nrow(data) < cv_folds) {
    stop("Conditional-mean CV requires at least two nonempty folds.")
  }
  set.seed(seed)
  fold_id <- sample(rep(seq_len(cv_folds), length.out = nrow(data)))
  outcome_matrix <- as.matrix(data[, outcomes, drop = FALSE])
  outcome_scale <- apply(outcome_matrix, 2L, stats::sd)
  outcome_scale[!is.finite(outcome_scale) | outcome_scale < 1e-8] <- 1
  risk_rows <- list()
  counter <- 1L
  for (spline_df in spline_df_grid) {
    for (fold in seq_len(cv_folds)) {
      validation <- fold_id == fold
      fit <- fit_spline_means(
        data = data[!validation, , drop = FALSE],
        outcomes = outcomes,
        x_cols = x_cols,
        spline_df = spline_df
      )
      prediction <- fit$predict(data[validation, , drop = FALSE])
      residual <- sweep(
        outcome_matrix[validation, , drop = FALSE] - prediction,
        2L,
        outcome_scale,
        "/"
      )
      risk_rows[[counter]] <- data.frame(
        spline_df = spline_df,
        fold = fold,
        risk = mean(residual^2),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }
  fold_risk <- do.call(rbind, risk_rows)
  cv_table <- do.call(rbind, lapply(
    split(fold_risk, fold_risk$spline_df),
    function(group) data.frame(
      spline_df = group$spline_df[1L],
      mean_risk = mean(group$risk),
      standard_error = stats::sd(group$risk) / sqrt(nrow(group)),
      stringsAsFactors = FALSE
    )
  ))
  rownames(cv_table) <- NULL
  minimum <- which.min(cv_table$mean_risk)
  selected <- if (selection_rule == "minimum") {
    minimum
  } else {
    cutoff <- cv_table$mean_risk[minimum] + cv_table$standard_error[minimum]
    eligible <- which(cv_table$mean_risk <= cutoff)
    eligible[which.min(cv_table$spline_df[eligible])]
  }
  cv_table$selected <- seq_len(nrow(cv_table)) == selected
  final <- fit_spline_means(
    data = data,
    outcomes = outcomes,
    x_cols = x_cols,
    spline_df = cv_table$spline_df[selected]
  )
  final$learner <- "spline_cv"
  final$selected_spline_df <- cv_table$spline_df[selected]
  final$effective_df <- final$rank
  final$cv_table <- cv_table
  final
}

fit_gam_means <- function(
    data,
    outcomes,
    x_cols = "X",
    basis_dimension = 6L,
    method = "REML"
) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("The optional conditional EIF learner 'gam' requires mgcv.")
  }
  non_numeric <- x_cols[
    !vapply(data[, x_cols, drop = FALSE], is.numeric, logical(1))
  ]
  if (length(non_numeric) > 0L) {
    stop(
      "Conditional EIF GAM predictors must be numeric: ",
      paste(non_numeric, collapse = ", ")
    )
  }
  unique_counts <- vapply(
    data[, x_cols, drop = FALSE],
    function(value) length(unique(value[is.finite(value)])),
    integer(1)
  )
  smooth_cols <- x_cols[unique_counts >= as.integer(basis_dimension)]
  linear_cols <- setdiff(x_cols, smooth_cols)
  smooth_terms <- sprintf(
    "s(`%s`, bs='cr', k=%d)",
    gsub("`", "", smooth_cols, fixed = TRUE),
    as.integer(basis_dimension)
  )
  linear_terms <- sprintf(
    "`%s`", gsub("`", "", linear_cols, fixed = TRUE)
  )
  model_terms <- c(linear_terms, smooth_terms)
  if (!length(model_terms)) model_terms <- "1"
  fits <- stats::setNames(lapply(outcomes, function(outcome) {
    response <- data[[outcome]]
    center <- mean(response)
    scale <- stats::sd(response)
    if (!is.finite(scale) || scale < 1e-10) {
      return(list(
        fit = NULL,
        center = center,
        scale = 0,
        rank = 1,
        effective_df = 1
      ))
    }
    local_data <- data
    local_data$.eif_gam_response <- (response - center) / scale
    formula <- stats::as.formula(paste0(
      ".eif_gam_response ~ ", paste(model_terms, collapse = " + ")
    ))
    environment(formula) <- asNamespace("mgcv")
    fit <- tryCatch(
      mgcv::gam(
        formula = formula,
        data = local_data,
        family = stats::gaussian(),
        method = method,
        select = TRUE
      ),
      error = function(error) {
        stop(
          "Conditional EIF GAM failed for outcome ", outcome, ": ",
          conditionMessage(error),
          call. = FALSE
        )
      }
    )
    list(
      fit = fit,
      center = center,
      scale = scale,
      rank = fit$rank,
      effective_df = sum(fit$edf)
    )
  }), outcomes)
  list(
    predict = function(newdata) {
      output <- vapply(fits, function(specification) {
        if (is.null(specification$fit)) {
          return(rep(specification$center, nrow(newdata)))
        }
        specification$center + specification$scale * as.numeric(
          stats::predict(
            specification$fit,
            newdata = newdata,
            type = "response"
          )
        )
      }, numeric(nrow(newdata)))
      output <- as.matrix(output)
      colnames(output) <- outcomes
      output
    },
    fits = fits,
    rank = sum(vapply(fits, function(specification) {
      specification$rank
    }, numeric(1))),
    outcomes = outcomes,
    learner = "gam",
    selected_spline_df = NA_integer_,
    smooth_predictors = smooth_cols,
    linear_predictors = linear_cols,
    effective_df = sum(vapply(
      fits,
      function(specification) specification$effective_df,
      numeric(1)
    )),
    cv_table = data.frame()
  )
}

fit_ranger_means <- function(
    data,
    outcomes,
    x_cols = "X",
    num_trees = 500L,
    min_node_size = NULL,
    seed = 1L
) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("The optional conditional EIF learner 'ranger' requires ranger.")
  }
  non_numeric <- x_cols[
    !vapply(data[, x_cols, drop = FALSE], is.numeric, logical(1))
  ]
  if (length(non_numeric) > 0L) {
    stop(
      "Conditional EIF ranger predictors must be numeric: ",
      paste(non_numeric, collapse = ", ")
    )
  }
  num_trees <- as.integer(num_trees)
  local_min_node_size <- if (is.null(min_node_size)) {
    max(5L, as.integer(floor(nrow(data)^(1 / 4))))
  } else {
    as.integer(min_node_size)
  }
  if (length(num_trees) != 1L || is.na(num_trees) || num_trees < 1L ||
      length(local_min_node_size) != 1L ||
      is.na(local_min_node_size) || local_min_node_size < 1L) {
    stop("Invalid conditional EIF ranger controls.")
  }

  fits <- stats::setNames(lapply(seq_along(outcomes), function(index) {
    outcome <- outcomes[index]
    response <- as.numeric(data[[outcome]])
    if (all(abs(response - response[1L]) < 1e-12)) {
      return(list(fit = NULL, constant = response[1L]))
    }
    ranger_data <- data[, x_cols, drop = FALSE]
    ranger_data$.eif_ranger_response <- response
    fit <- ranger::ranger(
      formula = stats::reformulate(
        x_cols, response = ".eif_ranger_response"
      ),
      data = ranger_data,
      num.trees = num_trees,
      mtry = max(1L, floor(sqrt(length(x_cols)))),
      min.node.size = local_min_node_size,
      seed = as.integer(seed + index - 1L),
      num.threads = 1L
    )
    list(fit = fit, constant = NA_real_)
  }), outcomes)

  list(
    predict = function(newdata) {
      output <- vapply(fits, function(specification) {
        if (is.null(specification$fit)) {
          return(rep(specification$constant, nrow(newdata)))
        }
        as.numeric(predict(
          specification$fit, data = newdata
        )$predictions)
      }, numeric(nrow(newdata)))
      output <- as.matrix(output)
      colnames(output) <- outcomes
      output
    },
    fits = fits,
    rank = length(x_cols) * length(outcomes),
    outcomes = outcomes,
    learner = "ranger",
    selected_spline_df = NA_integer_,
    effective_df = NA_real_,
    cv_table = data.frame(),
    num_trees = num_trees,
    min_node_size = local_min_node_size
  )
}

fit_eif_conditional_library_cv <- function(
    data,
    outcomes,
    x_cols = "X",
    candidates = c("linear", "gam", "ranger"),
    cv_folds = 3L,
    selection_rule = c("one_se", "minimum"),
    gam_k = 6L,
    gam_method = "REML",
    ranger_num_trees = 500L,
    ranger_min_node_size = NULL,
    seed = 1L
) {
  selection_rule <- match.arg(selection_rule)
  candidates <- unique(as.character(candidates))
  allowed <- c("linear", "gam", "ranger")
  if (length(candidates) < 1L || any(!candidates %in% allowed)) {
    stop(
      "Conditional EIF library candidates must be drawn from: ",
      paste(allowed, collapse = ", "), "."
    )
  }
  cv_folds <- as.integer(cv_folds)
  if (length(cv_folds) != 1L || is.na(cv_folds) ||
      cv_folds < 2L || nrow(data) < cv_folds) {
    stop("Conditional EIF library CV requires at least two nonempty folds.")
  }
  missing_outcomes <- setdiff(outcomes, names(data))
  if (length(missing_outcomes) > 0L) {
    stop(
      "Missing conditional EIF outcomes: ",
      paste(missing_outcomes, collapse = ", ")
    )
  }

  fit_candidate <- function(training, candidate, local_seed) {
    switch(
      candidate,
      linear = {
        fit <- fit_spline_means(
          data = training,
          outcomes = outcomes,
          x_cols = x_cols,
          spline_df = 1L
        )
        fit$learner <- "linear"
        fit$effective_df <- fit$rank
        fit
      },
      gam = fit_gam_means(
        data = training,
        outcomes = outcomes,
        x_cols = x_cols,
        basis_dimension = gam_k,
        method = gam_method
      ),
      ranger = fit_ranger_means(
        data = training,
        outcomes = outcomes,
        x_cols = x_cols,
        num_trees = ranger_num_trees,
        min_node_size = ranger_min_node_size,
        seed = local_seed
      )
    )
  }

  set.seed(as.integer(seed))
  fold_id <- sample(rep(seq_len(cv_folds), length.out = nrow(data)))
  outcome_matrix <- as.matrix(data[, outcomes, drop = FALSE])
  storage.mode(outcome_matrix) <- "double"
  outcome_scale <- apply(outcome_matrix, 2L, stats::sd)
  outcome_scale[!is.finite(outcome_scale) | outcome_scale < 1e-8] <- 1
  risk_rows <- list()
  counter <- 1L
  for (candidate_index in seq_along(candidates)) {
    candidate <- candidates[candidate_index]
    for (fold in seq_len(cv_folds)) {
      validation <- fold_id == fold
      fit <- fit_candidate(
        data[!validation, , drop = FALSE],
        candidate,
        as.integer(seed + 10000L * candidate_index + fold)
      )
      prediction <- fit$predict(data[validation, , drop = FALSE])
      residual <- sweep(
        outcome_matrix[validation, , drop = FALSE] - prediction,
        2L,
        outcome_scale,
        "/"
      )
      for (outcome in outcomes) {
        risk_rows[[counter]] <- data.frame(
          outcome = outcome,
          candidate = candidate,
          fold = fold,
          risk = mean(residual[, outcome]^2),
          observations = sum(validation),
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }
  }
  fold_risk <- do.call(rbind, risk_rows)
  groups <- split(
    fold_risk,
    interaction(
      fold_risk$outcome, fold_risk$candidate,
      drop = TRUE, lex.order = TRUE
    )
  )
  cv_table <- do.call(rbind, lapply(groups, function(group) {
    data.frame(
      outcome = group$outcome[1L],
      candidate = group$candidate[1L],
      mean_risk = mean(group$risk),
      standard_error = if (nrow(group) > 1L) {
        stats::sd(group$risk) / sqrt(nrow(group))
      } else {
        0
      },
      stringsAsFactors = FALSE
    )
  }))
  rownames(cv_table) <- NULL
  complexity <- c(linear = 1L, gam = 2L, ranger = 3L)
  selected <- stats::setNames(character(length(outcomes)), outcomes)
  cv_table$selected <- FALSE
  for (outcome in outcomes) {
    local <- which(cv_table$outcome == outcome)
    minimum <- local[which.min(cv_table$mean_risk[local])]
    choice <- minimum
    if (selection_rule == "one_se") {
      cutoff <- cv_table$mean_risk[minimum] +
        cv_table$standard_error[minimum]
      eligible <- local[cv_table$mean_risk[local] <= cutoff]
      choice <- eligible[which.min(complexity[cv_table$candidate[eligible]])]
    }
    cv_table$selected[choice] <- TRUE
    selected[outcome] <- cv_table$candidate[choice]
  }

  selected_candidates <- unique(unname(selected))
  final_fits <- stats::setNames(lapply(seq_along(selected_candidates), function(index) {
    fit_candidate(
      data,
      selected_candidates[index],
      as.integer(seed + 900000L + index)
    )
  }), selected_candidates)

  list(
    predict = function(newdata) {
      candidate_predictions <- lapply(final_fits, function(fit) {
        fit$predict(newdata)
      })
      output <- matrix(
        NA_real_, nrow = nrow(newdata), ncol = length(outcomes),
        dimnames = list(NULL, outcomes)
      )
      for (outcome in outcomes) {
        output[, outcome] <- candidate_predictions[[
          selected[[outcome]]
        ]][, outcome]
      }
      output
    },
    fits = final_fits,
    rank = NA_real_,
    outcomes = outcomes,
    learner = "library_cv",
    selected_learners = selected,
    selected_spline_df = NA_integer_,
    effective_df = NA_real_,
    cv_table = cv_table,
    configuration = list(
      candidates = candidates,
      cv_folds = cv_folds,
      selection_rule = selection_rule,
      selected_learners = selected
    )
  )
}

fit_eif_conditional_means <- function(
    data,
    outcomes,
    x_cols = "X",
    learner = c(
      "spline_fixed", "spline_cv", "gam", "ranger", "library_cv"
    ),
    spline_df = 5L,
    control = NULL,
    seed = 1L
) {
  learner <- match.arg(learner)
  control <- resolve_eif_conditional_control(control)
  if (learner == "spline_fixed") {
    fit <- fit_spline_means(
      data = data,
      outcomes = outcomes,
      x_cols = x_cols,
      spline_df = spline_df
    )
    fit$learner <- learner
    fit$selected_spline_df <- as.integer(spline_df)
    fit$effective_df <- fit$rank
    fit$cv_table <- data.frame()
    return(fit)
  }
  if (learner == "spline_cv") {
    return(fit_spline_means_cv(
      data = data,
      outcomes = outcomes,
      x_cols = x_cols,
      spline_df_grid = control$spline_df_grid,
      cv_folds = control$cv_folds,
      selection_rule = control$selection_rule,
      seed = seed
    ))
  }
  if (learner == "gam") {
    return(fit_gam_means(
      data = data,
      outcomes = outcomes,
      x_cols = x_cols,
      basis_dimension = control$gam_k,
      method = control$gam_method
    ))
  }
  if (learner == "library_cv") {
    return(fit_eif_conditional_library_cv(
      data = data,
      outcomes = outcomes,
      x_cols = x_cols,
      candidates = control$library_candidates,
      cv_folds = control$cv_folds,
      selection_rule = control$library_selection_rule,
      gam_k = control$gam_k,
      gam_method = control$gam_method,
      ranger_num_trees = control$ranger_num_trees,
      ranger_min_node_size = control$ranger_min_node_size,
      seed = seed
    ))
  }
  fit_ranger_means(
    data = data,
    outcomes = outcomes,
    x_cols = x_cols,
    num_trees = control$ranger_num_trees,
    min_node_size = control$ranger_min_node_size,
    seed = seed
  )
}

bind_conditional_cv_tables <- function(baseline, comparison) {
  rows <- list()
  if (nrow(baseline) > 0L) {
    baseline$stratum <- "baseline_z0"
    rows[[length(rows) + 1L]] <- baseline
  }
  if (nrow(comparison) > 0L) {
    comparison$stratum <- "comparison_z1"
    rows[[length(rows) + 1L]] <- comparison
  }
  if (length(rows) == 0L) return(data.frame())
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

check_loss <- function(
    y,
    q,
    level
) {
  residual <- y - q

  residual * (
    level -
      as.numeric(residual < 0)
  )
}

fit_threshold_by_check_loss <- function(
    data,
    y,
    level,
    weight,
    x_cols = "X",
    spline_df = 5,
    tie_break_ridge = 0,
    maxit = 4000,
    reltol = 1e-10,
    retry_on_iteration_limit = TRUE,
    retry_maxit_multiplier = 4,
    retry_reltol = 1e-8
) {
  if (
    length(y) != nrow(data) ||
      length(level) != nrow(data) ||
      length(weight) != nrow(data)
  ) {
    stop(
      "y, level, and weight must have one value per row."
    )
  }

  if (
    any(!is.finite(y)) ||
      any(!is.finite(level)) ||
      any(!is.finite(weight))
  ) {
    stop(
      "Non-finite values supplied to the threshold optimizer."
    )
  }

  if (
    any(level <= 0) ||
      any(level >= 1)
  ) {
    stop(
      "All check-loss levels must lie strictly between zero and one."
    )
  }

  if (
    any(weight < 0) ||
      sum(weight) <= 0
  ) {
    stop(
      "The check-loss weights must be nonnegative with positive sum."
    )
  }
  if (
    length(tie_break_ridge) != 1L ||
      !is.finite(tie_break_ridge) ||
      tie_break_ridge < 0
  ) {
    stop("tie_break_ridge must be one nonnegative finite number.")
  }
  if (
    length(maxit) != 1L || !is.finite(maxit) || maxit < 1 ||
      length(reltol) != 1L || !is.finite(reltol) || reltol <= 0 ||
      length(retry_on_iteration_limit) != 1L ||
      is.na(retry_on_iteration_limit) ||
      length(retry_maxit_multiplier) != 1L ||
      !is.finite(retry_maxit_multiplier) ||
      retry_maxit_multiplier < 1 ||
      length(retry_reltol) != 1L ||
      !is.finite(retry_reltol) || retry_reltol <= 0
  ) {
    stop("Invalid threshold-optimizer controls.")
  }
  maxit <- as.integer(maxit)
  retry_maxit <- as.integer(ceiling(maxit * retry_maxit_multiplier))

  design <- make_spline_design(
    data = data,
    x_cols = x_cols,
    spline_df = spline_df
  )

  # Optimize in a standardized spline coordinate system.  The intercept is
  # left on its original scale, while every nonconstant spline column is
  # centered and divided by its training standard deviation.  This makes the
  # coefficient ridge comparable across spline degrees of freedom without
  # changing the represented threshold class.  The ridge below applies to the
  # complete standardized coefficient vector, including the intercept, so it
  # selects a unique vector even when constant thresholds are tied.
  design_center <- rep(0, ncol(design$matrix))
  design_scale <- rep(1, ncol(design$matrix))
  if (ncol(design$matrix) > 1L) {
    nonintercept_columns <- 2:ncol(design$matrix)
    design_center[nonintercept_columns] <- colMeans(
      design$matrix[, nonintercept_columns, drop = FALSE]
    )
    design_scale[nonintercept_columns] <- apply(
      design$matrix[, nonintercept_columns, drop = FALSE],
      2L,
      stats::sd
    )
    design_scale[
      !is.finite(design_scale) | design_scale < sqrt(.Machine$double.eps)
    ] <- 1
  }
  standardize_design <- function(matrix) {
    sweep(sweep(matrix, 2L, design_center, "-"), 2L, design_scale, "/")
  }
  optimization_matrix <- standardize_design(design$matrix)

  initial_fit <- stats::lm.wfit(
    x = optimization_matrix,
    y = y,
    w = pmax(weight, 1e-8)
  )

  initial <- initial_fit$coefficients
  initial[!is.finite(initial)] <- 0

  objective <- function(beta) {
    q <- as.vector(
      optimization_matrix %*% beta
    )

    empirical_risk <- sum(
      weight *
        check_loss(
          y = y,
          q = q,
          level = level
        )
    ) / sum(weight)

    empirical_risk +
      0.5 * tie_break_ridge * sum(beta^2)
  }

  initial_optimization <- stats::optim(
    par = initial,
    fn = objective,
    method = "Nelder-Mead",
    control = list(
      maxit = maxit,
      reltol = reltol
    )
  )

  # Nelder-Mead occasionally reaches its iteration limit on the nonsmooth,
  # piecewise-linear check-loss objective even after arriving at a stable
  # solution.  Continue deterministically from that solution with a larger
  # iteration budget and a still-stringent numerical tolerance.  This retry
  # does not change the threshold class, ridge, or empirical objective.
  retry_used <- isTRUE(retry_on_iteration_limit) &&
    identical(as.integer(initial_optimization$convergence), 1L)
  retry_optimization <- NULL
  if (retry_used) {
    retry_optimization <- stats::optim(
      par = initial_optimization$par,
      fn = objective,
      method = "Nelder-Mead",
      control = list(
        maxit = retry_maxit,
        reltol = max(reltol, retry_reltol)
      )
    )
  }

  retry_objective_tolerance <- if (retry_used) {
    sqrt(.Machine$double.eps) * max(
      1,
      abs(initial_optimization$value),
      abs(retry_optimization$value)
    )
  } else {
    NA_real_
  }
  retry_accepted <- retry_used &&
    is.finite(retry_optimization$value) &&
    retry_optimization$value <=
      initial_optimization$value + retry_objective_tolerance
  optimization <- if (retry_accepted) {
    retry_optimization
  } else {
    initial_optimization
  }

  combine_counts <- function(first, second = NULL) {
    if (is.null(second)) return(first)
    count_names <- union(names(first), names(second))
    output <- vapply(count_names, function(name) {
      values <- c(first[[name]], second[[name]])
      if (all(is.na(values))) NA_integer_ else sum(values, na.rm = TRUE)
    }, integer(1L))
    output
  }
  optimizer_counts <- combine_counts(
    initial_optimization$counts,
    if (retry_used) retry_optimization$counts else NULL
  )

  coefficients <- optimization$par
  fitted_q <- as.vector(optimization_matrix %*% coefficients)
  empirical_objective <- sum(
    weight * check_loss(y, fitted_q, level)
  ) / sum(weight)

  list(
    predict = function(newdata) {
      as.vector(
        standardize_design(design$predict_matrix(newdata)) %*%
          coefficients
      )
    },
    coefficients = coefficients,
    design_center = design_center,
    design_scale = design_scale,
    objective = empirical_objective,
    penalized_objective = optimization$value,
    tie_break_ridge = tie_break_ridge,
    convergence = optimization$convergence,
    counts = optimizer_counts,
    initial_convergence = initial_optimization$convergence,
    initial_penalized_objective = initial_optimization$value,
    retry_used = retry_used,
    retry_accepted = retry_accepted,
    retry_convergence = if (retry_used) {
      retry_optimization$convergence
    } else {
      NA_integer_
    },
    retry_penalized_objective = if (retry_used) {
      retry_optimization$value
    } else {
      NA_real_
    },
    retry_maxit = if (retry_used) retry_maxit else NA_integer_,
    retry_reltol = if (retry_used) max(reltol, retry_reltol) else NA_real_
  )
}

default_threshold_control <- function() {
  list(
    spline_df_grid = 1:5,
    cv_folds = 3L,
    selection_rule = "one_se",
    # spline_tie_break_ridge is retained for explicitly requested fixed-ridge
    # sensitivity analyses.  The official vanishing tie-break is specified by
    # spline_tie_break_ridge_numerator = c and applies c/m in a fit with m
    # training observations.  At most one of the two entries may be positive.
    # Held-out bounds are always evaluated through the original, unpenalized
    # fixed-q functional.
    spline_tie_break_ridge = 0,
    spline_tie_break_ridge_numerator = 0,
    # The exact full-kernel backend is the scientific reference.  Nyström is
    # retained as an explicitly named computational approximation.
    rkhs_backend = "exact",
    rkhs_lambda_grid = 10^seq(-3, 0, length.out = 4),
    rkhs_bandwidth_multipliers = c(0.5, 1, 2),
    rkhs_rank = 40L,
    rkhs_bandwidth_sample = 1000L,
    # The check-loss weights are normalized to sum to one.  A modest ADMM
    # penalty therefore balances the loss and consensus residual much better
    # than rho = 1, especially once the baseline stratum contains hundreds of
    # observations.
    rkhs_admm_rho = 0.01,
    rkhs_admm_maxit = 3000L,
    rkhs_admm_abs_tol = 1e-5,
    rkhs_admm_rel_tol = 1e-4,
    rkhs_admm_check_interval = 5L,
    rkhs_condition_tolerance = 1e-12,
    rkhs_numerical_jitter = 1e-10,
    # A fixed finite-rank Gaussian-RBF class provides a deterministic
    # restricted threshold target.  The centers and bandwidth define the
    # class and must therefore be supplied rather than learned from the same
    # data.  The c/m ridge is a vanishing numerical tie-break and is omitted
    # when evaluating the held-out bound functional.
    fixed_rbf_centers = NULL,
    fixed_rbf_bandwidth = 1,
    fixed_rbf_ridge_numerator = 1,
    fixed_rbf_admm_rho = 0.01,
    fixed_rbf_admm_maxit = 6000L,
    fixed_rbf_admm_abs_tol = 1e-5,
    fixed_rbf_admm_rel_tol = 1e-4,
    # The growing RBF sieve is an automatic, increasingly rich threshold
    # learner rather than a user-specified restricted class.  Its rank is
    # ceil(multiplier * m^exponent), where m is the realized threshold
    # training size.  Deterministic Halton points are mapped through the
    # training empirical marginals, bandwidth is selected by honest CV, and
    # c/m is used only as a vanishing numerical tie-break.
    growing_rbf_rank_multiplier = 6,
    growing_rbf_rank_exponent = 1 / 3,
    growing_rbf_min_rank = 12L,
    growing_rbf_center_skip = 20L,
    growing_rbf_bandwidth_rule = "absolute",
    growing_rbf_bandwidth_grid = c(0.5, 1, 2),
    growing_rbf_bandwidth_sample = 1000L,
    growing_rbf_ridge_numerator = 1,
    growing_rbf_admm_rho = 0.01,
    growing_rbf_admm_maxit = 6000L,
    growing_rbf_admm_abs_tol = 1e-5,
    growing_rbf_admm_rel_tol = 1e-4
  )
}

resolve_threshold_control <- function(control = NULL) {
  output <- default_threshold_control()
  if (!is.null(control)) {
    unknown <- setdiff(names(control), names(output))
    if (length(unknown) > 0L) {
      stop(
        "Unknown threshold_control entries: ",
        paste(unknown, collapse = ", ")
      )
    }
    output[names(control)] <- control
  }
  output$selection_rule <- match.arg(
    output$selection_rule,
    c("one_se", "minimum")
  )
  output$rkhs_backend <- match.arg(
    output$rkhs_backend,
    c("exact", "nystrom")
  )
  output$growing_rbf_bandwidth_rule <- match.arg(
    output$growing_rbf_bandwidth_rule,
    c("absolute", "median_multiplier")
  )
  output$cv_folds <- as.integer(output$cv_folds)
  output$spline_df_grid <- sort(unique(as.integer(
    output$spline_df_grid
  )))
  if (!is.null(output$rkhs_rank)) {
    output$rkhs_rank <- as.integer(output$rkhs_rank)
  }
  output$rkhs_bandwidth_sample <- as.integer(
    output$rkhs_bandwidth_sample
  )
  output$rkhs_admm_maxit <- as.integer(output$rkhs_admm_maxit)
  output$rkhs_admm_check_interval <- as.integer(
    output$rkhs_admm_check_interval
  )
  if (!is.null(output$fixed_rbf_centers)) {
    output$fixed_rbf_centers <- as.matrix(output$fixed_rbf_centers)
    storage.mode(output$fixed_rbf_centers) <- "double"
  }
  output$fixed_rbf_admm_maxit <- as.integer(
    output$fixed_rbf_admm_maxit
  )
  output$growing_rbf_min_rank <- as.integer(
    output$growing_rbf_min_rank
  )
  output$growing_rbf_center_skip <- as.integer(
    output$growing_rbf_center_skip
  )
  output$growing_rbf_bandwidth_sample <- as.integer(
    output$growing_rbf_bandwidth_sample
  )
  output$growing_rbf_admm_maxit <- as.integer(
    output$growing_rbf_admm_maxit
  )
  if (
    output$spline_tie_break_ridge > 0 &&
      output$spline_tie_break_ridge_numerator > 0
  ) {
    stop(
      paste0(
        "Specify at most one of spline_tie_break_ridge and ",
        "spline_tie_break_ridge_numerator."
      )
    )
  }
  if (
    length(output$cv_folds) != 1L ||
      is.na(output$cv_folds) ||
      output$cv_folds < 2L ||
      length(output$spline_df_grid) == 0L ||
      any(!is.finite(output$spline_df_grid)) ||
      any(output$spline_df_grid < 1L) ||
      length(output$spline_tie_break_ridge) != 1L ||
      !is.finite(output$spline_tie_break_ridge) ||
      output$spline_tie_break_ridge < 0 ||
      length(output$spline_tie_break_ridge_numerator) != 1L ||
      !is.finite(output$spline_tie_break_ridge_numerator) ||
      output$spline_tie_break_ridge_numerator < 0 ||
      length(output$rkhs_lambda_grid) == 0L ||
      any(!is.finite(output$rkhs_lambda_grid)) ||
      any(output$rkhs_lambda_grid <= 0) ||
      length(output$rkhs_bandwidth_multipliers) == 0L ||
      any(!is.finite(output$rkhs_bandwidth_multipliers)) ||
      any(output$rkhs_bandwidth_multipliers <= 0) ||
      (!is.null(output$rkhs_rank) && (
        length(output$rkhs_rank) != 1L ||
          is.na(output$rkhs_rank) ||
          output$rkhs_rank < 1L
      )) ||
      length(output$rkhs_bandwidth_sample) != 1L ||
      is.na(output$rkhs_bandwidth_sample) ||
      output$rkhs_bandwidth_sample < 2L ||
      length(output$rkhs_admm_rho) != 1L ||
      !is.finite(output$rkhs_admm_rho) ||
      output$rkhs_admm_rho <= 0 ||
      length(output$rkhs_admm_maxit) != 1L ||
      is.na(output$rkhs_admm_maxit) ||
      output$rkhs_admm_maxit < 1L ||
      length(output$rkhs_admm_abs_tol) != 1L ||
      !is.finite(output$rkhs_admm_abs_tol) ||
      output$rkhs_admm_abs_tol <= 0 ||
      length(output$rkhs_admm_rel_tol) != 1L ||
      !is.finite(output$rkhs_admm_rel_tol) ||
      output$rkhs_admm_rel_tol <= 0 ||
      length(output$rkhs_admm_check_interval) != 1L ||
      is.na(output$rkhs_admm_check_interval) ||
      output$rkhs_admm_check_interval < 1L ||
      length(output$rkhs_condition_tolerance) != 1L ||
      !is.finite(output$rkhs_condition_tolerance) ||
      output$rkhs_condition_tolerance <= 0 ||
      length(output$rkhs_numerical_jitter) != 1L ||
      !is.finite(output$rkhs_numerical_jitter) ||
      output$rkhs_numerical_jitter < 0 ||
      (!is.null(output$fixed_rbf_centers) && (
        nrow(output$fixed_rbf_centers) < 1L ||
          ncol(output$fixed_rbf_centers) < 1L ||
          any(!is.finite(output$fixed_rbf_centers))
      )) ||
      length(output$fixed_rbf_bandwidth) != 1L ||
      !is.finite(output$fixed_rbf_bandwidth) ||
      output$fixed_rbf_bandwidth <= 0 ||
      length(output$fixed_rbf_ridge_numerator) != 1L ||
      !is.finite(output$fixed_rbf_ridge_numerator) ||
      output$fixed_rbf_ridge_numerator <= 0 ||
      length(output$fixed_rbf_admm_rho) != 1L ||
      !is.finite(output$fixed_rbf_admm_rho) ||
      output$fixed_rbf_admm_rho <= 0 ||
      length(output$fixed_rbf_admm_maxit) != 1L ||
      is.na(output$fixed_rbf_admm_maxit) ||
      output$fixed_rbf_admm_maxit < 1L ||
      length(output$fixed_rbf_admm_abs_tol) != 1L ||
      !is.finite(output$fixed_rbf_admm_abs_tol) ||
      output$fixed_rbf_admm_abs_tol <= 0 ||
      length(output$fixed_rbf_admm_rel_tol) != 1L ||
      !is.finite(output$fixed_rbf_admm_rel_tol) ||
      output$fixed_rbf_admm_rel_tol <= 0 ||
      length(output$growing_rbf_rank_multiplier) != 1L ||
      !is.finite(output$growing_rbf_rank_multiplier) ||
      output$growing_rbf_rank_multiplier <= 0 ||
      length(output$growing_rbf_rank_exponent) != 1L ||
      !is.finite(output$growing_rbf_rank_exponent) ||
      output$growing_rbf_rank_exponent <= 0 ||
      output$growing_rbf_rank_exponent >= 1 ||
      length(output$growing_rbf_min_rank) != 1L ||
      is.na(output$growing_rbf_min_rank) ||
      output$growing_rbf_min_rank < 1L ||
      length(output$growing_rbf_center_skip) != 1L ||
      is.na(output$growing_rbf_center_skip) ||
      output$growing_rbf_center_skip < 0L ||
      length(output$growing_rbf_bandwidth_grid) < 1L ||
      any(!is.finite(output$growing_rbf_bandwidth_grid)) ||
      any(output$growing_rbf_bandwidth_grid <= 0) ||
      length(output$growing_rbf_bandwidth_sample) != 1L ||
      is.na(output$growing_rbf_bandwidth_sample) ||
      output$growing_rbf_bandwidth_sample < 2L ||
      length(output$growing_rbf_ridge_numerator) != 1L ||
      !is.finite(output$growing_rbf_ridge_numerator) ||
      output$growing_rbf_ridge_numerator <= 0 ||
      length(output$growing_rbf_admm_rho) != 1L ||
      !is.finite(output$growing_rbf_admm_rho) ||
      output$growing_rbf_admm_rho <= 0 ||
      length(output$growing_rbf_admm_maxit) != 1L ||
      is.na(output$growing_rbf_admm_maxit) ||
      output$growing_rbf_admm_maxit < 1L ||
      length(output$growing_rbf_admm_abs_tol) != 1L ||
      !is.finite(output$growing_rbf_admm_abs_tol) ||
      output$growing_rbf_admm_abs_tol <= 0 ||
      length(output$growing_rbf_admm_rel_tol) != 1L ||
      !is.finite(output$growing_rbf_admm_rel_tol) ||
      output$growing_rbf_admm_rel_tol <= 0
  ) {
    stop("Invalid threshold_control tuning or optimizer values.")
  }
  output
}

resolve_threshold_tie_break_ridge <- function(
    training_size,
    fixed_ridge = 0,
    ridge_numerator = 0
) {
  training_size <- as.integer(training_size)
  if (
    length(training_size) != 1L || is.na(training_size) ||
      training_size < 1L
  ) {
    stop("training_size must be one positive integer.")
  }
  if (
    length(fixed_ridge) != 1L || !is.finite(fixed_ridge) ||
      fixed_ridge < 0 ||
      length(ridge_numerator) != 1L || !is.finite(ridge_numerator) ||
      ridge_numerator < 0
  ) {
    stop("Threshold ridge inputs must be nonnegative finite scalars.")
  }
  if (fixed_ridge > 0 && ridge_numerator > 0) {
    stop("A fixed ridge and a vanishing ridge cannot both be positive.")
  }
  if (fixed_ridge > 0) fixed_ridge else ridge_numerator / training_size
}

make_threshold_cv_folds <- function(n, K = 3L, seed = 1L) {
  K <- as.integer(K)
  if (K < 2L || n < K) {
    stop("Threshold CV requires at least two folds and one row per fold.")
  }
  set.seed(seed)
  sample(rep(seq_len(K), length.out = n))
}

select_threshold_cv_row <- function(
    cv_table,
    selection_rule = c("one_se", "minimum"),
    complexity_columns
) {
  selection_rule <- match.arg(selection_rule)
  minimum_index <- which.min(cv_table$mean_risk)
  if (selection_rule == "minimum") return(minimum_index)
  cutoff <- cv_table$mean_risk[minimum_index] +
    cv_table$standard_error[minimum_index]
  eligible <- which(cv_table$mean_risk <= cutoff)
  ordering_arguments <- lapply(
    complexity_columns,
    function(specification) {
      value <- cv_table[[specification$name]][eligible]
      if (isTRUE(specification$decreasing)) -value else value
    }
  )
  ordering <- do.call(order, ordering_arguments)
  eligible[ordering[1L]]
}

fit_threshold_spline_cv <- function(
    data,
    y,
    level,
    weight,
    x_cols = "X",
    spline_df_grid = 1:5,
    tie_break_ridge = 0,
    tie_break_ridge_numerator = 0,
    cv_folds = 3L,
    selection_rule = c("one_se", "minimum"),
    maxit = 4000L,
    seed = 1L
) {
  selection_rule <- match.arg(selection_rule)
  spline_df_grid <- sort(unique(as.integer(spline_df_grid)))
  if (any(spline_df_grid < 1L)) {
    stop("All spline_df_grid values must be positive.")
  }
  fold_id <- make_threshold_cv_folds(nrow(data), cv_folds, seed)
  risk_rows <- list()
  counter <- 1L

  for (spline_df in spline_df_grid) {
    for (fold in seq_len(cv_folds)) {
      validation <- fold_id == fold
      training <- !validation
      local_tie_break_ridge <- resolve_threshold_tie_break_ridge(
        training_size = sum(training),
        fixed_ridge = tie_break_ridge,
        ridge_numerator = tie_break_ridge_numerator
      )
      fit <- fit_threshold_by_check_loss(
        data = data[training, , drop = FALSE],
        y = y[training],
        level = level[training],
        weight = weight[training],
        x_cols = x_cols,
        spline_df = spline_df,
        tie_break_ridge = local_tie_break_ridge,
        maxit = maxit
      )
      prediction <- fit$predict(data[validation, , drop = FALSE])
      local_loss <- check_loss(
        y[validation],
        prediction,
        level[validation]
      )
      risk_rows[[counter]] <- data.frame(
        spline_df = spline_df,
        fold = fold,
        risk = sum(weight[validation] * local_loss) /
          sum(weight[validation]),
        training_size = sum(training),
        tie_break_ridge = local_tie_break_ridge,
        convergence = fit$convergence,
        initial_convergence = fit$initial_convergence,
        retry_used = fit$retry_used,
        retry_accepted = fit$retry_accepted,
        retry_convergence = fit$retry_convergence,
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }
  fold_risk <- do.call(rbind, risk_rows)
  groups <- split(fold_risk, fold_risk$spline_df)
  cv_table <- do.call(rbind, lapply(groups, function(group) {
    data.frame(
      spline_df = group$spline_df[1L],
      mean_risk = mean(group$risk),
      standard_error = stats::sd(group$risk) / sqrt(nrow(group)),
      minimum_training_size = min(group$training_size),
      maximum_training_size = max(group$training_size),
      minimum_tie_break_ridge = min(group$tie_break_ridge),
      maximum_tie_break_ridge = max(group$tie_break_ridge),
      maximum_convergence = max(group$convergence),
      maximum_initial_convergence = max(group$initial_convergence),
      any_retry_used = any(group$retry_used),
      all_retries_accepted = if (any(group$retry_used)) {
        all(group$retry_accepted[group$retry_used])
      } else {
        NA
      },
      maximum_retry_convergence = if (any(group$retry_used)) {
        max(group$retry_convergence[group$retry_used])
      } else {
        NA_integer_
      },
      stringsAsFactors = FALSE
    )
  }))
  rownames(cv_table) <- NULL
  selected_index <- select_threshold_cv_row(
    cv_table,
    selection_rule,
    complexity_columns = list(list(
      name = "spline_df",
      decreasing = FALSE
    ))
  )
  cv_table$selected <- seq_len(nrow(cv_table)) == selected_index
  selected_df <- cv_table$spline_df[selected_index]
  final_tie_break_ridge <- resolve_threshold_tie_break_ridge(
    training_size = nrow(data),
    fixed_ridge = tie_break_ridge,
    ridge_numerator = tie_break_ridge_numerator
  )
  final <- fit_threshold_by_check_loss(
    data = data,
    y = y,
    level = level,
    weight = weight,
    x_cols = x_cols,
    spline_df = selected_df,
    tie_break_ridge = final_tie_break_ridge,
    maxit = maxit
  )
  final$learner <- "spline_cv"
  final$tie_break_ridge_numerator <- tie_break_ridge_numerator
  final$tie_break_ridge_training_size <- nrow(data)
  final$selected_spline_df <- selected_df
  final$selected_lambda <- NA_real_
  final$selected_bandwidth_multiplier <- NA_real_
  final$selected_bandwidth <- NA_real_
  final$base_bandwidth <- NA_real_
  final$selected_rank <- NA_integer_
  final$cv_table <- cv_table
  final$primal_residual <- NA_real_
  final$dual_residual <- NA_real_
  final$final_rho <- NA_real_
  final
}

prepare_threshold_rkhs_coordinates <- function(
    data,
    x_cols = "X",
    bandwidth_sample = 1000L,
    seed = 1L
) {
  coordinates <- as.matrix(data[, x_cols, drop = FALSE])
  center <- colMeans(coordinates)
  scale_value <- apply(coordinates, 2, stats::sd)
  scale_value[!is.finite(scale_value) | scale_value < 1e-8] <- 1
  transform <- function(newdata) {
    value <- as.matrix(newdata[, x_cols, drop = FALSE])
    sweep(sweep(value, 2, center, "-"), 2, scale_value, "/")
  }
  transformed <- transform(data)
  set.seed(seed)
  sample_index <- sample(
    seq_len(nrow(transformed)),
    min(as.integer(bandwidth_sample), nrow(transformed))
  )
  distances <- as.numeric(stats::dist(
    transformed[sample_index, , drop = FALSE]
  ))
  distances <- distances[is.finite(distances) & distances > 0]
  base_bandwidth <- if (length(distances) == 0L) {
    1
  } else {
    stats::median(distances)
  }
  list(
    transformed = transformed,
    transform = transform,
    center = center,
    scale = scale_value,
    base_bandwidth = max(base_bandwidth, 1e-4)
  )
}

make_threshold_nystrom_features <- function(
    data,
    x_cols = "X",
    bandwidth_multiplier = 1,
    rank = 40L,
    bandwidth_sample = 1000L,
    seed = 1L
) {
  coordinate_spec <- prepare_threshold_rkhs_coordinates(
    data,
    x_cols,
    bandwidth_sample,
    seed
  )
  sigma <- coordinate_spec$base_bandwidth * bandwidth_multiplier
  rank <- min(as.integer(rank), nrow(data))
  set.seed(seed + 1L)
  center_index <- if (rank == nrow(data)) {
    seq_len(nrow(data))
  } else {
    sample(seq_len(nrow(data)), rank)
  }
  centers <- coordinate_spec$transformed[
    center_index,
    ,
    drop = FALSE
  ]
  K_mm <- rbf_kernel(centers, centers, sigma)
  eigen_decomposition <- eigen(K_mm, symmetric = TRUE)
  keep <- eigen_decomposition$values >
    max(eigen_decomposition$values) * 1e-10
  vectors <- eigen_decomposition$vectors[, keep, drop = FALSE]
  values <- eigen_decomposition$values[keep]
  feature_map <- function(newdata) {
    transformed <- coordinate_spec$transform(newdata)
    kernel <- rbf_kernel(transformed, centers, sigma)
    sweep(kernel %*% vectors, 2, sqrt(values), "/")
  }
  list(
    features = feature_map(data),
    predict_features = feature_map,
    base_bandwidth = coordinate_spec$base_bandwidth,
    bandwidth_multiplier = bandwidth_multiplier,
    bandwidth = sigma,
    rank = length(values)
  )
}

make_threshold_exact_kernel <- function(
    data,
    x_cols = "X",
    bandwidth_multiplier = 1,
    bandwidth_sample = 1000L,
    seed = 1L
) {
  coordinate_spec <- prepare_threshold_rkhs_coordinates(
    data,
    x_cols,
    bandwidth_sample,
    seed
  )
  sigma <- coordinate_spec$base_bandwidth * bandwidth_multiplier
  kernel <- rbf_kernel(
    coordinate_spec$transformed,
    coordinate_spec$transformed,
    sigma
  )
  kernel <- (kernel + t(kernel)) / 2
  decomposition <- eigen(kernel, symmetric = TRUE)
  original_values <- decomposition$values
  numerical_scale <- max(abs(original_values), 1)
  values <- pmax(original_values, 0)
  effective_rank <- sum(values > max(values) * 1e-10)

  list(
    kernel = kernel,
    eigenvectors = decomposition$vectors,
    eigenvalues = values,
    minimum_raw_eigenvalue = min(original_values),
    negative_eigenvalue_mass = sum(abs(
      original_values[original_values < -1e-12 * numerical_scale]
    )),
    effective_rank = effective_rank,
    base_bandwidth = coordinate_spec$base_bandwidth,
    bandwidth_multiplier = bandwidth_multiplier,
    bandwidth = sigma,
    transform = coordinate_spec$transform,
    training_coordinates = coordinate_spec$transformed,
    predict_kernel = function(newdata) {
      rbf_kernel(
        coordinate_spec$transform(newdata),
        coordinate_spec$transformed,
        sigma
      )
    }
  )
}

fit_weighted_quantile_kernel_admm <- function(
    kernel_spec,
    y,
    level,
    weight,
    lambda,
    rho = 0.01,
    maxit = 3000L,
    abs_tol = 1e-5,
    rel_tol = 1e-4,
    check_interval = 5L,
    condition_tolerance = 1e-12,
    numerical_jitter = 1e-10
) {
  kernel <- kernel_spec$kernel
  vectors <- kernel_spec$eigenvectors
  values <- kernel_spec$eigenvalues
  n <- nrow(kernel)
  if (
    ncol(kernel) != n ||
      length(y) != n ||
      length(level) != n ||
      length(weight) != n
  ) {
    stop("Exact-kernel ADMM inputs must have a common row count.")
  }
  if (
    any(!is.finite(y)) ||
      any(!is.finite(level)) ||
      any(level <= 0 | level >= 1) ||
      any(!is.finite(weight)) ||
      any(weight < 0) ||
      sum(weight) <= 0 ||
      length(lambda) != 1L ||
      !is.finite(lambda) ||
      lambda <= 0
  ) {
    stop("Invalid exact-kernel weighted quantile inputs.")
  }
  normalized_weight <- weight / sum(weight)
  one <- rep(1, n)
  projected_one <- as.vector(crossprod(vectors, one))
  current_rho <- rho
  cached_rho <- NA_real_
  inverse_one <- NULL
  effective_c <- NA_real_
  current_reciprocal_condition <- NA_real_
  minimum_reciprocal_condition <- Inf
  maximum_jitter <- 0

  update_regularized_inverse <- function() {
    c_value <- lambda / current_rho
    denominator <- values + c_value
    denominator_scale <- max(denominator)
    reciprocal_condition <- min(denominator) / denominator_scale
    jitter <- 0
    if (
      !is.finite(reciprocal_condition) ||
        reciprocal_condition < condition_tolerance
    ) {
      jitter <- max(
        numerical_jitter * max(denominator_scale, 1),
        condition_tolerance * denominator_scale - min(denominator)
      )
    }
    denominator <- denominator + jitter
    local_inverse_one <- as.vector(
      vectors %*% (projected_one / denominator)
    )
    denominator_one <- sum(local_inverse_one)
    if (!is.finite(denominator_one) || abs(denominator_one) < 1e-12) {
      stop("The exact-kernel intercept Schur complement is singular.")
    }
    list(
      denominator = denominator,
      inverse_one = local_inverse_one,
      denominator_one = denominator_one,
      effective_c = c_value + jitter,
      reciprocal_condition = min(denominator) / max(denominator),
      jitter = jitter
    )
  }

  inverse_state <- update_regularized_inverse()
  cached_rho <- current_rho
  inverse_one <- inverse_state$inverse_one
  effective_c <- inverse_state$effective_c
  minimum_reciprocal_condition <- min(
    minimum_reciprocal_condition,
    inverse_state$reciprocal_condition
  )
  maximum_jitter <- max(maximum_jitter, inverse_state$jitter)

  intercept <- sum(weight * y) / sum(weight)
  alpha <- numeric(n)
  fitted <- rep(intercept, n)
  residual <- y - fitted
  dual <- numeric(n)
  convergence <- 1L
  primal_norm <- dual_norm <- Inf
  primal_tolerance <- dual_tolerance <- NA_real_

  for (iteration in seq_len(as.integer(maxit))) {
    if (current_rho != cached_rho) {
      inverse_state <- update_regularized_inverse()
      cached_rho <- current_rho
      inverse_one <- inverse_state$inverse_one
      effective_c <- inverse_state$effective_c
      minimum_reciprocal_condition <- min(
        minimum_reciprocal_condition,
        inverse_state$reciprocal_condition
      )
      maximum_jitter <- max(maximum_jitter, inverse_state$jitter)
    }

    residual_previous <- residual
    target <- y - residual - dual
    inverse_target <- as.vector(
      vectors %*% (
        as.vector(crossprod(vectors, target)) /
          inverse_state$denominator
      )
    )
    intercept <- sum(inverse_target) /
      inverse_state$denominator_one
    alpha <- inverse_target - intercept * inverse_one
    # From (K + c I) alpha = target - intercept, this is algebraically
    # identical to intercept + K alpha and avoids another dense multiplication.
    fitted <- target - effective_c * alpha

    proximal_center <- y - fitted - dual
    threshold_positive <- normalized_weight * level / current_rho
    threshold_negative <-
      normalized_weight * (1 - level) / current_rho
    residual <- ifelse(
      proximal_center > threshold_positive,
      proximal_center - threshold_positive,
      ifelse(
        proximal_center < -threshold_negative,
        proximal_center + threshold_negative,
        0
      )
    )
    primal <- fitted + residual - y
    dual <- dual + primal

    check_now <- iteration %% as.integer(check_interval) == 0L ||
      iteration == 1L || iteration == as.integer(maxit)
    if (check_now) {
      residual_change <- residual - residual_previous
      kernel_residual_change <- as.vector(kernel %*% residual_change)
      dual_norm <- current_rho * sqrt(
        sum(residual_change)^2 +
          sum(kernel_residual_change^2)
      )
      primal_norm <- sqrt(sum(primal^2))
      primal_tolerance <- sqrt(n) * abs_tol + rel_tol * max(
        sqrt(sum(fitted^2)),
        sqrt(sum(residual^2)),
        sqrt(sum(y^2))
      )
      kernel_dual <- as.vector(kernel %*% dual)
      dual_tolerance <- sqrt(n + 1L) * abs_tol +
        rel_tol * current_rho * sqrt(
          sum(dual)^2 + sum(kernel_dual^2)
        )
      if (
        primal_norm <= primal_tolerance &&
          dual_norm <= dual_tolerance
      ) {
        convergence <- 0L
        break
      }
    }

    if (iteration %% 25L == 0L) {
      new_rho <- current_rho
      if (primal_norm > 10 * dual_norm) {
        new_rho <- min(current_rho * 2, 100)
      } else if (dual_norm > 10 * primal_norm) {
        new_rho <- max(current_rho / 2, 1e-5)
      }
      if (new_rho != current_rho) {
        dual <- dual * current_rho / new_rho
        current_rho <- new_rho
      }
    }
  }

  # Recompute from the actual kernel for the returned objective and fitted
  # values. Conditional jitter is therefore visible rather than silently
  # changing the target functional.
  fitted <- intercept + as.vector(kernel %*% alpha)
  rkhs_norm_squared <- sum(
    values * as.vector(crossprod(vectors, alpha))^2
  )
  unpenalized_risk <- sum(
    weight * check_loss(y, fitted, level)
  ) / sum(weight)

  list(
    intercept = intercept,
    alpha = alpha,
    fitted = fitted,
    convergence = convergence,
    iterations = iteration,
    primal_residual = primal_norm,
    dual_residual = dual_norm,
    primal_tolerance = primal_tolerance,
    dual_tolerance = dual_tolerance,
    final_rho = current_rho,
    objective = unpenalized_risk,
    penalized_objective = unpenalized_risk +
      lambda * rkhs_norm_squared / 2,
    rkhs_norm_squared = rkhs_norm_squared,
    minimum_reciprocal_condition = minimum_reciprocal_condition,
    maximum_jitter = maximum_jitter,
    effective_rank = kernel_spec$effective_rank,
    predict_kernel = function(new_kernel) {
      intercept + as.vector(as.matrix(new_kernel) %*% alpha)
    }
  )
}

fit_weighted_quantile_ridge_admm <- function(
    features,
    y,
    level,
    weight,
    lambda,
    rho = 1,
    maxit = 1000L,
    abs_tol = 1e-5,
    rel_tol = 1e-4
) {
  features <- as.matrix(features)
  design <- cbind(intercept = 1, features)
  n <- nrow(design)
  p <- ncol(design)
  normalized_weight <- weight / sum(weight)
  penalty <- diag(c(0, rep(lambda, p - 1L)), p)
  current_rho <- rho
  update_factor <- function() {
    system <- current_rho * crossprod(design) +
      penalty + diag(1e-10, p)
    chol(system)
  }
  factor <- update_factor()
  solve_system <- function(rhs) {
    backsolve(factor, forwardsolve(t(factor), rhs))
  }
  beta <- stats::lm.wfit(
    design,
    y,
    pmax(weight, 1e-8)
  )$coefficients
  beta[!is.finite(beta)] <- 0
  residual <- y - as.vector(design %*% beta)
  dual <- numeric(n)
  convergence <- 1L
  primal_norm <- dual_norm <- Inf

  for (iteration in seq_len(as.integer(maxit))) {
    residual_previous <- residual
    rhs <- current_rho * crossprod(
      design,
      y - residual - dual
    )
    beta <- as.vector(solve_system(rhs))
    proximal_center <- y - as.vector(design %*% beta) - dual
    threshold_positive <- normalized_weight * level / current_rho
    threshold_negative <-
      normalized_weight * (1 - level) / current_rho
    residual <- ifelse(
      proximal_center > threshold_positive,
      proximal_center - threshold_positive,
      ifelse(
        proximal_center < -threshold_negative,
        proximal_center + threshold_negative,
        0
      )
    )
    primal <- as.vector(design %*% beta) + residual - y
    dual <- dual + primal
    primal_norm <- sqrt(sum(primal^2))
    dual_vector <- current_rho * crossprod(
      design,
      residual - residual_previous
    )
    dual_norm <- sqrt(sum(dual_vector^2))
    primal_tolerance <- sqrt(n) * abs_tol + rel_tol * max(
      sqrt(sum((design %*% beta)^2)),
      sqrt(sum(residual^2)),
      sqrt(sum(y^2))
    )
    dual_tolerance <- sqrt(p) * abs_tol +
      rel_tol * sqrt(sum(
        (current_rho * crossprod(design, dual))^2
      ))
    if (
      primal_norm <= primal_tolerance &&
        dual_norm <= dual_tolerance
    ) {
      convergence <- 0L
      break
    }
    # Residual balancing substantially reduces the number of iterations over
    # the range of lambda values considered by CV.  Because `dual` is the
    # scaled multiplier, rescale it whenever rho changes.
    if (iteration %% 25L == 0L) {
      new_rho <- current_rho
      if (primal_norm > 10 * dual_norm) {
        new_rho <- min(current_rho * 2, 100)
      } else if (dual_norm > 10 * primal_norm) {
        new_rho <- max(current_rho / 2, 1e-5)
      }
      if (new_rho != current_rho) {
        dual <- dual * current_rho / new_rho
        current_rho <- new_rho
        factor <- update_factor()
      }
    }
  }
  fitted <- as.vector(design %*% beta)
  unpenalized_risk <- sum(
    weight * check_loss(y, fitted, level)
  ) / sum(weight)
  penalized_objective <- unpenalized_risk +
    lambda * sum(beta[-1L]^2) / 2
  list(
    coefficients = beta,
    fitted = fitted,
    convergence = convergence,
    iterations = iteration,
    primal_residual = primal_norm,
    dual_residual = dual_norm,
    final_rho = current_rho,
    objective = unpenalized_risk,
    penalized_objective = penalized_objective,
    predict_features = function(new_features) {
      as.vector(cbind(1, as.matrix(new_features)) %*% beta)
    }
  )
}

polish_weighted_quantile_ridge_fit <- function(
    features,
    y,
    level,
    weight,
    lambda,
    initial,
    maxit = 6000L,
    reltol = 1e-9
) {
  features <- as.matrix(features)
  design <- cbind(intercept = 1, features)
  if (
    nrow(design) != length(y) || length(y) != length(level) ||
      length(y) != length(weight) || length(initial) != ncol(design)
  ) {
    stop("Invalid weighted-quantile polish dimensions.")
  }
  objective <- function(beta) {
    fitted <- as.vector(design %*% beta)
    sum(weight * check_loss(y, fitted, level)) / sum(weight) +
      lambda * sum(beta[-1L]^2) / 2
  }
  optimization <- stats::optim(
    par = as.numeric(initial),
    fn = objective,
    method = "Nelder-Mead",
    control = list(maxit = as.integer(maxit), reltol = reltol)
  )
  coefficients <- optimization$par
  fitted <- as.vector(design %*% coefficients)
  unpenalized_risk <- sum(
    weight * check_loss(y, fitted, level)
  ) / sum(weight)
  list(
    coefficients = coefficients,
    fitted = fitted,
    convergence = optimization$convergence,
    counts = optimization$counts,
    objective = unpenalized_risk,
    penalized_objective = optimization$value,
    predict_features = function(new_features) {
      as.vector(cbind(1, as.matrix(new_features)) %*% coefficients)
    }
  )
}

make_threshold_fixed_rbf_features <- function(
    data,
    x_cols,
    centers,
    bandwidth
) {
  non_numeric <- x_cols[
    !vapply(data[, x_cols, drop = FALSE], is.numeric, logical(1L))
  ]
  if (length(non_numeric)) {
    stop(
      "Fixed-RBF threshold coordinates must be numeric: ",
      paste(non_numeric, collapse = ", "), "."
    )
  }
  centers <- as.matrix(centers)
  storage.mode(centers) <- "double"
  if (!is.null(colnames(centers)) && all(x_cols %in% colnames(centers))) {
    centers <- centers[, x_cols, drop = FALSE]
  }
  if (
    nrow(centers) < 1L || ncol(centers) != length(x_cols) ||
      any(!is.finite(centers)) ||
      length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0
  ) {
    stop(
      "Fixed-RBF centers must be a finite rank-by-length(x_cols) matrix ",
      "and bandwidth must be positive."
    )
  }
  feature_map <- function(newdata) {
    coordinates <- as.matrix(newdata[, x_cols, drop = FALSE])
    if (any(!is.finite(coordinates))) {
      stop("Fixed-RBF prediction coordinates must be finite.")
    }
    rbf_kernel(coordinates, centers, bandwidth) / sqrt(nrow(centers))
  }
  list(
    features = feature_map(data),
    predict_features = feature_map,
    centers = centers,
    bandwidth = bandwidth,
    rank = nrow(centers)
  )
}

fit_threshold_fixed_rbf <- function(
    data,
    y,
    level,
    weight,
    x_cols,
    centers,
    bandwidth,
    ridge_numerator = 1,
    rho = 0.01,
    maxit = 6000L,
    abs_tol = 1e-5,
    rel_tol = 1e-4
) {
  feature_specification <- make_threshold_fixed_rbf_features(
    data = data,
    x_cols = x_cols,
    centers = centers,
    bandwidth = bandwidth
  )
  ridge <- as.numeric(ridge_numerator) / nrow(data)
  fit <- fit_weighted_quantile_ridge_admm(
    features = feature_specification$features,
    y = y,
    level = level,
    weight = weight,
    lambda = ridge,
    rho = rho,
    maxit = maxit,
    abs_tol = abs_tol,
    rel_tol = rel_tol
  )
  list(
    predict = function(newdata) {
      fit$predict_features(feature_specification$predict_features(newdata))
    },
    coefficients = fit$coefficients,
    objective = fit$objective,
    penalized_objective = fit$penalized_objective,
    tie_break_ridge = ridge,
    tie_break_ridge_numerator = ridge_numerator,
    tie_break_ridge_training_size = nrow(data),
    convergence = fit$convergence,
    counts = c(iterations = fit$iterations),
    learner = "fixed_rbf",
    rkhs_backend = "fixed_rbf",
    selected_spline_df = NA_integer_,
    selected_lambda = ridge,
    selected_bandwidth_multiplier = NA_real_,
    selected_bandwidth = bandwidth,
    base_bandwidth = bandwidth,
    selected_rank = feature_specification$rank,
    cv_table = data.frame(),
    primal_residual = fit$primal_residual,
    dual_residual = fit$dual_residual,
    final_rho = fit$final_rho,
    minimum_reciprocal_condition = NA_real_,
    maximum_numerical_jitter = NA_real_,
    kernel_minimum_raw_eigenvalue = NA_real_,
    centers = feature_specification$centers
  )
}

resolve_growing_rbf_rank <- function(
    n_training,
    multiplier = 6,
    exponent = 1 / 3,
    minimum_rank = 12L
) {
  n_training <- as.integer(n_training)
  minimum_rank <- as.integer(minimum_rank)
  if (
    length(n_training) != 1L || is.na(n_training) || n_training < 1L ||
      length(multiplier) != 1L || !is.finite(multiplier) ||
      multiplier <= 0 ||
      length(exponent) != 1L || !is.finite(exponent) ||
      exponent <= 0 || exponent >= 1 ||
      length(minimum_rank) != 1L || is.na(minimum_rank) ||
      minimum_rank < 1L
  ) {
    stop("Invalid growing-RBF rank rule.")
  }
  min(
    n_training,
    max(minimum_rank, as.integer(ceiling(multiplier * n_training^exponent)))
  )
}

resolve_growing_rbf_base_bandwidth <- function(
    data,
    x_cols,
    rule = c("absolute", "median_multiplier"),
    bandwidth_sample = 1000L
) {
  rule <- match.arg(rule)
  if (rule == "absolute") return(1)
  coordinates <- as.matrix(data[, x_cols, drop = FALSE])
  center <- colMeans(coordinates)
  scale_value <- apply(coordinates, 2L, stats::sd)
  scale_value[!is.finite(scale_value) | scale_value < 1e-8] <- 1
  transformed <- sweep(
    sweep(coordinates, 2L, center, "-"), 2L, scale_value, "/"
  )
  sample_size <- min(as.integer(bandwidth_sample), nrow(transformed))
  sample_index <- unique(as.integer(round(seq(
    1, nrow(transformed), length.out = sample_size
  ))))
  distance <- as.numeric(stats::dist(
    transformed[sample_index, , drop = FALSE]
  ))
  distance <- distance[is.finite(distance) & distance > 0]
  if (length(distance)) max(stats::median(distance), 1e-4) else 1
}

make_threshold_growing_rbf_features <- function(
    data,
    x_cols,
    bandwidth,
    base_bandwidth = 1,
    rank_multiplier = 6,
    rank_exponent = 1 / 3,
    minimum_rank = 12L,
    center_skip = 20L
) {
  non_numeric <- x_cols[
    !vapply(data[, x_cols, drop = FALSE], is.numeric, logical(1L))
  ]
  if (length(non_numeric)) {
    stop(
      "Growing-RBF threshold coordinates must be numeric: ",
      paste(non_numeric, collapse = ", "), "."
    )
  }
  coordinates <- as.matrix(data[, x_cols, drop = FALSE])
  if (!nrow(coordinates) || any(!is.finite(coordinates))) {
    stop("Growing-RBF training coordinates must be finite and nonempty.")
  }
  coordinate_center <- colMeans(coordinates)
  coordinate_scale <- apply(coordinates, 2L, stats::sd)
  coordinate_scale[
    !is.finite(coordinate_scale) | coordinate_scale < 1e-8
  ] <- 1
  transform <- function(newdata) {
    value <- as.matrix(newdata[, x_cols, drop = FALSE])
    if (any(!is.finite(value))) {
      stop("Growing-RBF prediction coordinates must be finite.")
    }
    sweep(
      sweep(value, 2L, coordinate_center, "-"),
      2L, coordinate_scale, "/"
    )
  }
  transformed <- transform(data)
  rank <- resolve_growing_rbf_rank(
    n_training = nrow(data),
    multiplier = rank_multiplier,
    exponent = rank_exponent,
    minimum_rank = minimum_rank
  )
  center_skip <- as.integer(center_skip)
  if (length(center_skip) != 1L || is.na(center_skip) || center_skip < 0L) {
    stop("center_skip must be one nonnegative integer.")
  }
  index <- center_skip + seq_len(rank)
  bases <- first_primes(length(x_cols))
  unit <- vapply(
    bases,
    function(base) radical_inverse(index, base),
    numeric(rank)
  )
  if (length(x_cols) == 1L) unit <- matrix(unit, ncol = 1L)
  centers <- vapply(
    seq_along(x_cols),
    function(column) {
      as.numeric(stats::quantile(
        transformed[, column],
        probs = unit[, column],
        names = FALSE,
        type = 8
      ))
    },
    numeric(rank)
  )
  if (length(x_cols) == 1L) centers <- matrix(centers, ncol = 1L)
  colnames(centers) <- x_cols
  feature_map <- function(newdata) {
    rbf_kernel(transform(newdata), centers, bandwidth) / sqrt(rank)
  }
  list(
    features = feature_map(data),
    predict_features = feature_map,
    centers = centers,
    bandwidth = as.numeric(bandwidth),
    base_bandwidth = as.numeric(base_bandwidth),
    rank = rank,
    coordinate_center = coordinate_center,
    coordinate_scale = coordinate_scale
  )
}

make_threshold_growing_rbf_plan <- function(
    data,
    x_cols = "X",
    bandwidth_grid = c(0.5, 1, 2),
    rank_multiplier = 6,
    rank_exponent = 1 / 3,
    minimum_rank = 12L,
    center_skip = 20L,
    bandwidth_rule = c("absolute", "median_multiplier"),
    bandwidth_sample = 1000L,
    cv_folds = 3L,
    seed = 1L
) {
  bandwidth_rule <- match.arg(bandwidth_rule)
  fold_id <- make_threshold_cv_folds(nrow(data), cv_folds, seed)
  make_specifications <- function(local_data) {
    base_bandwidth <- resolve_growing_rbf_base_bandwidth(
      data = local_data,
      x_cols = x_cols,
      rule = bandwidth_rule,
      bandwidth_sample = bandwidth_sample
    )
    output <- lapply(bandwidth_grid, function(bandwidth_candidate) {
      bandwidth <- if (bandwidth_rule == "median_multiplier") {
        base_bandwidth * bandwidth_candidate
      } else {
        bandwidth_candidate
      }
      make_threshold_growing_rbf_features(
        data = local_data,
        x_cols = x_cols,
        bandwidth = bandwidth,
        base_bandwidth = base_bandwidth,
        rank_multiplier = rank_multiplier,
        rank_exponent = rank_exponent,
        minimum_rank = minimum_rank,
        center_skip = center_skip
      )
    })
    names(output) <- vapply(
      bandwidth_grid, threshold_bandwidth_key, character(1L)
    )
    output
  }
  folds <- lapply(seq_len(cv_folds), function(fold) {
    validation <- fold_id == fold
    training <- !validation
    list(
      training = training,
      validation = validation,
      specifications = make_specifications(data[training, , drop = FALSE])
    )
  })
  list(
    backend = "growing_rbf_sieve",
    n = nrow(data),
    fold_id = fold_id,
    folds = folds,
    final_specifications = make_specifications(data),
    bandwidth_grid = as.numeric(bandwidth_grid),
    bandwidth_rule = bandwidth_rule,
    bandwidth_sample = as.integer(bandwidth_sample),
    rank_multiplier = as.numeric(rank_multiplier),
    rank_exponent = as.numeric(rank_exponent),
    minimum_rank = as.integer(minimum_rank),
    center_skip = as.integer(center_skip),
    seed = as.integer(seed)
  )
}

fit_threshold_growing_rbf_cv <- function(
    data,
    y,
    level,
    weight,
    x_cols = "X",
    bandwidth_grid = c(0.5, 1, 2),
    rank_multiplier = 6,
    rank_exponent = 1 / 3,
    minimum_rank = 12L,
    center_skip = 20L,
    bandwidth_rule = c("absolute", "median_multiplier"),
    bandwidth_sample = 1000L,
    ridge_numerator = 1,
    cv_folds = 3L,
    selection_rule = c("minimum", "one_se"),
    rho = 0.01,
    maxit = 6000L,
    abs_tol = 1e-5,
    rel_tol = 1e-4,
    plan = NULL,
    seed = 1L
) {
  selection_rule <- match.arg(selection_rule)
  bandwidth_rule <- match.arg(bandwidth_rule)
  if (is.null(plan)) {
    plan <- make_threshold_growing_rbf_plan(
      data = data,
      x_cols = x_cols,
      bandwidth_grid = bandwidth_grid,
      rank_multiplier = rank_multiplier,
      rank_exponent = rank_exponent,
      minimum_rank = minimum_rank,
      center_skip = center_skip,
      bandwidth_rule = bandwidth_rule,
      bandwidth_sample = bandwidth_sample,
      cv_folds = cv_folds,
      seed = seed
    )
  }
  if (
    !identical(plan$backend, "growing_rbf_sieve") ||
      plan$n != nrow(data) || length(plan$folds) != cv_folds ||
      !identical(as.numeric(plan$bandwidth_grid), as.numeric(bandwidth_grid)) ||
      !identical(plan$bandwidth_rule, bandwidth_rule)
  ) {
    stop("The supplied growing-RBF plan does not match this fit.")
  }
  risk_rows <- list()
  counter <- 1L
  for (fold in seq_len(cv_folds)) {
    training <- plan$folds[[fold]]$training
    validation <- plan$folds[[fold]]$validation
    local_ridge <- as.numeric(ridge_numerator) / sum(training)
    for (bandwidth in bandwidth_grid) {
      key <- threshold_bandwidth_key(bandwidth)
      specification <- plan$folds[[fold]]$specifications[[key]]
      fit <- fit_weighted_quantile_ridge_admm(
        features = specification$features,
        y = y[training],
        level = level[training],
        weight = weight[training],
        lambda = local_ridge,
        rho = rho,
        maxit = maxit,
        abs_tol = abs_tol,
        rel_tol = rel_tol
      )
      prediction <- fit$predict_features(
        specification$predict_features(data[validation, , drop = FALSE])
      )
      local_loss <- check_loss(
        y[validation], prediction, level[validation]
      )
      risk_rows[[counter]] <- data.frame(
        bandwidth = bandwidth,
        fold = fold,
        risk = sum(weight[validation] * local_loss) /
          sum(weight[validation]),
        training_size = sum(training),
        rank = specification$rank,
        base_bandwidth = specification$base_bandwidth,
        realized_bandwidth = specification$bandwidth,
        tie_break_ridge = local_ridge,
        convergence = fit$convergence,
        iterations = fit$iterations,
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }
  fold_risk <- do.call(rbind, risk_rows)
  groups <- split(fold_risk, fold_risk$bandwidth)
  cv_table <- do.call(rbind, lapply(groups, function(local) {
    data.frame(
      bandwidth = local$bandwidth[1L],
      mean_risk = mean(local$risk),
      standard_error = stats::sd(local$risk) / sqrt(nrow(local)),
      minimum_training_size = min(local$training_size),
      maximum_training_size = max(local$training_size),
      minimum_rank = min(local$rank),
      maximum_rank = max(local$rank),
      mean_base_bandwidth = mean(local$base_bandwidth),
      mean_realized_bandwidth = mean(local$realized_bandwidth),
      minimum_tie_break_ridge = min(local$tie_break_ridge),
      maximum_tie_break_ridge = max(local$tie_break_ridge),
      maximum_convergence = max(local$convergence),
      mean_iterations = mean(local$iterations),
      rkhs_backend = "growing_rbf_sieve",
      stringsAsFactors = FALSE
    )
  }))
  rownames(cv_table) <- NULL
  selected_index <- select_threshold_cv_row(
    cv_table,
    selection_rule,
    complexity_columns = list(list(
      name = "bandwidth", decreasing = TRUE
    ))
  )
  cv_table$selected <- seq_len(nrow(cv_table)) == selected_index
  selected_bandwidth <- cv_table$bandwidth[selected_index]
  selected_key <- threshold_bandwidth_key(selected_bandwidth)
  final_specification <- plan$final_specifications[[selected_key]]
  final_ridge <- as.numeric(ridge_numerator) / nrow(data)
  final <- fit_weighted_quantile_ridge_admm(
    features = final_specification$features,
    y = y,
    level = level,
    weight = weight,
    lambda = final_ridge,
    rho = rho,
    maxit = maxit,
    abs_tol = abs_tol,
    rel_tol = rel_tol
  )
  polish <- if (final$convergence == 0L) {
    NULL
  } else {
    polish_weighted_quantile_ridge_fit(
      features = final_specification$features,
      y = y,
      level = level,
      weight = weight,
      lambda = final_ridge,
      initial = final$coefficients,
      maxit = maxit
    )
  }
  polish_tolerance <- sqrt(.Machine$double.eps) * max(
    1,
    abs(final$penalized_objective),
    if (is.null(polish)) 0 else abs(polish$penalized_objective)
  )
  polish_accepted <- !is.null(polish) &&
    is.finite(polish$penalized_objective) &&
    polish$penalized_objective <= final$penalized_objective + polish_tolerance
  selected_fit <- if (polish_accepted) polish else final
  list(
    predict = function(newdata) {
      selected_fit$predict_features(
        final_specification$predict_features(newdata)
      )
    },
    coefficients = selected_fit$coefficients,
    objective = selected_fit$objective,
    penalized_objective = selected_fit$penalized_objective,
    tie_break_ridge = final_ridge,
    tie_break_ridge_numerator = ridge_numerator,
    tie_break_ridge_training_size = nrow(data),
    convergence = selected_fit$convergence,
    counts = c(iterations = final$iterations),
    initial_convergence = final$convergence,
    initial_penalized_objective = final$penalized_objective,
    retry_used = !is.null(polish),
    retry_accepted = polish_accepted,
    retry_convergence = if (is.null(polish)) {
      NA_integer_
    } else {
      polish$convergence
    },
    retry_penalized_objective = if (is.null(polish)) {
      NA_real_
    } else {
      polish$penalized_objective
    },
    learner = "growing_rbf_sieve",
    rkhs_backend = "growing_rbf_sieve",
    selected_spline_df = NA_integer_,
    selected_lambda = final_ridge,
    selected_bandwidth_multiplier = if (
      bandwidth_rule == "median_multiplier"
    ) selected_bandwidth else NA_real_,
    selected_bandwidth = final_specification$bandwidth,
    base_bandwidth = final_specification$base_bandwidth,
    selected_rank = final_specification$rank,
    cv_table = cv_table,
    primal_residual = final$primal_residual,
    dual_residual = final$dual_residual,
    final_rho = final$final_rho,
    minimum_reciprocal_condition = NA_real_,
    maximum_numerical_jitter = NA_real_,
    kernel_minimum_raw_eigenvalue = NA_real_,
    centers = final_specification$centers,
    coordinate_center = final_specification$coordinate_center,
    coordinate_scale = final_specification$coordinate_scale,
    plan = plan
  )
}

threshold_bandwidth_key <- function(value) {
  format(as.numeric(value), digits = 17, scientific = FALSE, trim = TRUE)
}

resolve_threshold_nystrom_rank <- function(rank, n_training) {
  resolve_gamma_nystrom_rank(rank, n_training)
}

make_threshold_rkhs_plan <- function(
    data,
    x_cols = "X",
    backend = c("exact", "nystrom"),
    bandwidth_multipliers = c(0.5, 1, 2),
    rank = 40L,
    cv_folds = 3L,
    bandwidth_sample = 1000L,
    seed = 1L
) {
  backend <- match.arg(backend)
  fold_id <- make_threshold_cv_folds(nrow(data), cv_folds, seed)
  make_specification <- function(local_data, multiplier, local_seed) {
    if (backend == "exact") {
      make_threshold_exact_kernel(
        data = local_data,
        x_cols = x_cols,
        bandwidth_multiplier = multiplier,
        bandwidth_sample = bandwidth_sample,
        seed = local_seed
      )
    } else {
      local_rank <- resolve_threshold_nystrom_rank(rank, nrow(local_data))
      make_threshold_nystrom_features(
        data = local_data,
        x_cols = x_cols,
        bandwidth_multiplier = multiplier,
        rank = local_rank,
        bandwidth_sample = bandwidth_sample,
        seed = local_seed
      )
    }
  }
  folds <- vector("list", cv_folds)
  for (fold in seq_len(cv_folds)) {
    validation <- fold_id == fold
    training <- !validation
    specifications <- list()
    for (multiplier in bandwidth_multipliers) {
      key <- threshold_bandwidth_key(multiplier)
      specifications[[key]] <- make_specification(
        data[training, , drop = FALSE],
        multiplier,
        seed + 1000L * fold
      )
    }
    folds[[fold]] <- list(
      training = training,
      validation = validation,
      specifications = specifications
    )
  }
  final_specifications <- list()
  for (multiplier in bandwidth_multipliers) {
    key <- threshold_bandwidth_key(multiplier)
    final_specifications[[key]] <- make_specification(
      data,
      multiplier,
      seed + 900000L
    )
  }
  list(
    backend = backend,
    n = nrow(data),
    fold_id = fold_id,
    folds = folds,
    final_specifications = final_specifications,
    bandwidth_multipliers = bandwidth_multipliers,
    rank = rank,
    seed = seed
  )
}

fit_threshold_rkhs_cv <- function(
    data,
    y,
    level,
    weight,
    x_cols = "X",
    backend = c("exact", "nystrom"),
    lambda_grid = 10^seq(-3, 0, length.out = 4),
    bandwidth_multipliers = c(0.5, 1, 2),
    rank = 40L,
    cv_folds = 3L,
    selection_rule = c("one_se", "minimum"),
    bandwidth_sample = 1000L,
    admm_rho = 0.01,
    admm_maxit = 3000L,
    admm_abs_tol = 1e-5,
    admm_rel_tol = 1e-4,
    admm_check_interval = 5L,
    condition_tolerance = 1e-12,
    numerical_jitter = 1e-10,
    plan = NULL,
    seed = 1L
) {
  backend <- match.arg(backend)
  selection_rule <- match.arg(selection_rule)
  if (is.null(plan)) {
    plan <- make_threshold_rkhs_plan(
      data = data,
      x_cols = x_cols,
      backend = backend,
      bandwidth_multipliers = bandwidth_multipliers,
      rank = rank,
      cv_folds = cv_folds,
      bandwidth_sample = bandwidth_sample,
      seed = seed
    )
  }
  if (
    !identical(plan$backend, backend) ||
      plan$n != nrow(data) ||
      length(plan$folds) != cv_folds
  ) {
    stop("The supplied threshold RKHS plan does not match this fit.")
  }
  risk_rows <- list()
  counter <- 1L
  for (fold in seq_len(cv_folds)) {
    validation <- plan$folds[[fold]]$validation
    training <- plan$folds[[fold]]$training
    for (bandwidth_multiplier in bandwidth_multipliers) {
      key <- threshold_bandwidth_key(bandwidth_multiplier)
      specification <- plan$folds[[fold]]$specifications[[key]]
      if (is.null(specification)) {
        stop("The threshold RKHS plan omits a bandwidth candidate.")
      }
      validation_representation <- if (backend == "exact") {
        specification$predict_kernel(data[validation, , drop = FALSE])
      } else {
        specification$predict_features(data[validation, , drop = FALSE])
      }
      for (lambda in lambda_grid) {
        fit <- if (backend == "exact") {
          fit_weighted_quantile_kernel_admm(
            kernel_spec = specification,
            y = y[training],
            level = level[training],
            weight = weight[training],
            lambda = lambda,
            rho = admm_rho,
            maxit = admm_maxit,
            abs_tol = admm_abs_tol,
            rel_tol = admm_rel_tol,
            check_interval = admm_check_interval,
            condition_tolerance = condition_tolerance,
            numerical_jitter = numerical_jitter
          )
        } else {
          fit_weighted_quantile_ridge_admm(
            specification$features,
            y[training],
            level[training],
            weight[training],
            lambda,
            admm_rho,
            admm_maxit,
            admm_abs_tol,
            admm_rel_tol
          )
        }
        prediction <- if (backend == "exact") {
          fit$predict_kernel(validation_representation)
        } else {
          fit$predict_features(validation_representation)
        }
        local_loss <- check_loss(
          y[validation],
          prediction,
          level[validation]
        )
        risk_rows[[counter]] <- data.frame(
          lambda = lambda,
          bandwidth_multiplier = bandwidth_multiplier,
          fold = fold,
          risk = sum(weight[validation] * local_loss) /
            sum(weight[validation]),
          convergence = fit$convergence,
          iterations = fit$iterations,
          reciprocal_condition = if (backend == "exact") {
            fit$minimum_reciprocal_condition
          } else {
            NA_real_
          },
          numerical_jitter = if (backend == "exact") {
            fit$maximum_jitter
          } else {
            NA_real_
          },
          effective_rank = if (backend == "exact") {
            fit$effective_rank
          } else {
            specification$rank
          },
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }
  }
  fold_risk <- do.call(rbind, risk_rows)
  group <- interaction(
    fold_risk$lambda,
    fold_risk$bandwidth_multiplier,
    drop = TRUE
  )
  groups <- split(fold_risk, group)
  cv_table <- do.call(rbind, lapply(groups, function(local) {
    data.frame(
      lambda = local$lambda[1L],
      bandwidth_multiplier = local$bandwidth_multiplier[1L],
      mean_risk = mean(local$risk),
      standard_error = stats::sd(local$risk) / sqrt(nrow(local)),
      maximum_convergence = max(local$convergence),
      mean_iterations = mean(local$iterations),
      minimum_reciprocal_condition = if (
        all(is.na(local$reciprocal_condition))
      ) {
        NA_real_
      } else {
        min(local$reciprocal_condition, na.rm = TRUE)
      },
      maximum_numerical_jitter = if (
        all(is.na(local$numerical_jitter))
      ) {
        NA_real_
      } else {
        max(local$numerical_jitter, na.rm = TRUE)
      },
      mean_effective_rank = mean(local$effective_rank),
      rkhs_backend = backend,
      stringsAsFactors = FALSE
    )
  }))
  rownames(cv_table) <- NULL
  selected_index <- select_threshold_cv_row(
    cv_table,
    selection_rule,
    complexity_columns = list(
      list(name = "lambda", decreasing = TRUE),
      list(name = "bandwidth_multiplier", decreasing = TRUE)
    )
  )
  cv_table$selected <- seq_len(nrow(cv_table)) == selected_index
  selected_lambda <- cv_table$lambda[selected_index]
  selected_multiplier <-
    cv_table$bandwidth_multiplier[selected_index]
  selected_key <- threshold_bandwidth_key(selected_multiplier)
  final_specification <- plan$final_specifications[[selected_key]]
  final <- if (backend == "exact") {
    fit_weighted_quantile_kernel_admm(
      kernel_spec = final_specification,
      y = y,
      level = level,
      weight = weight,
      lambda = selected_lambda,
      rho = admm_rho,
      maxit = admm_maxit,
      abs_tol = admm_abs_tol,
      rel_tol = admm_rel_tol,
      check_interval = admm_check_interval,
      condition_tolerance = condition_tolerance,
      numerical_jitter = numerical_jitter
    )
  } else {
    fit_weighted_quantile_ridge_admm(
      final_specification$features,
      y,
      level,
      weight,
      selected_lambda,
      admm_rho,
      admm_maxit,
      admm_abs_tol,
      admm_rel_tol
    )
  }
  list(
    predict = function(newdata) {
      if (backend == "exact") {
        final$predict_kernel(final_specification$predict_kernel(newdata))
      } else {
        final$predict_features(final_specification$predict_features(newdata))
      }
    },
    coefficients = if (backend == "exact") {
      c(intercept = final$intercept, final$alpha)
    } else {
      final$coefficients
    },
    objective = final$objective,
    penalized_objective = final$penalized_objective,
    tie_break_ridge = NA_real_,
    tie_break_ridge_numerator = NA_real_,
    tie_break_ridge_training_size = nrow(data),
    convergence = final$convergence,
    counts = c(iterations = final$iterations),
    learner = "rkhs",
    rkhs_backend = backend,
    selected_spline_df = NA_integer_,
    selected_lambda = selected_lambda,
    selected_bandwidth_multiplier = selected_multiplier,
    selected_bandwidth = final_specification$bandwidth,
    base_bandwidth = final_specification$base_bandwidth,
    selected_rank = if (backend == "exact") {
      final_specification$effective_rank
    } else {
      final_specification$rank
    },
    cv_table = cv_table,
    primal_residual = final$primal_residual,
    dual_residual = final$dual_residual,
    final_rho = final$final_rho,
    minimum_reciprocal_condition = if (backend == "exact") {
      final$minimum_reciprocal_condition
    } else {
      NA_real_
    },
    maximum_numerical_jitter = if (backend == "exact") {
      final$maximum_jitter
    } else {
      NA_real_
    },
    kernel_minimum_raw_eigenvalue = if (backend == "exact") {
      final_specification$minimum_raw_eigenvalue
    } else {
      NA_real_
    },
    plan = plan
  )
}

fit_threshold_learner <- function(
    data,
    y,
    level,
    weight,
    x_cols = "X",
    learner = c(
      "spline_fixed", "spline_cv", "rkhs", "fixed_rbf",
      "growing_rbf_sieve"
    ),
    spline_df = 5L,
    maxit = 4000L,
    control = NULL,
    rkhs_plan = NULL,
    seed = 1L
) {
  learner <- match.arg(learner)
  control <- resolve_threshold_control(control)
  if (learner == "spline_fixed") {
    realized_tie_break_ridge <- resolve_threshold_tie_break_ridge(
      training_size = nrow(data),
      fixed_ridge = control$spline_tie_break_ridge,
      ridge_numerator = control$spline_tie_break_ridge_numerator
    )
    fit <- fit_threshold_by_check_loss(
      data,
      y,
      level,
      weight,
      x_cols,
      spline_df,
      realized_tie_break_ridge,
      maxit
    )
    fit$learner <- learner
    fit$tie_break_ridge_numerator <-
      control$spline_tie_break_ridge_numerator
    fit$tie_break_ridge_training_size <- nrow(data)
    fit$selected_spline_df <- spline_df
    fit$selected_lambda <- NA_real_
    fit$selected_bandwidth_multiplier <- NA_real_
    fit$selected_bandwidth <- NA_real_
    fit$base_bandwidth <- NA_real_
    fit$selected_rank <- NA_integer_
    fit$cv_table <- data.frame()
    fit$primal_residual <- NA_real_
    fit$dual_residual <- NA_real_
    fit$final_rho <- NA_real_
    fit$rkhs_backend <- NA_character_
    fit$minimum_reciprocal_condition <- NA_real_
    fit$maximum_numerical_jitter <- NA_real_
    fit$kernel_minimum_raw_eigenvalue <- NA_real_
    return(fit)
  }
  if (learner == "spline_cv") {
    fit <- fit_threshold_spline_cv(
      data = data,
      y = y,
      level = level,
      weight = weight,
      x_cols = x_cols,
      spline_df_grid = control$spline_df_grid,
      tie_break_ridge = control$spline_tie_break_ridge,
      tie_break_ridge_numerator =
        control$spline_tie_break_ridge_numerator,
      cv_folds = control$cv_folds,
      selection_rule = control$selection_rule,
      maxit = maxit,
      seed = seed
    )
    fit$rkhs_backend <- NA_character_
    fit$minimum_reciprocal_condition <- NA_real_
    fit$maximum_numerical_jitter <- NA_real_
    fit$kernel_minimum_raw_eigenvalue <- NA_real_
    return(fit)
  }
  if (learner == "fixed_rbf") {
    if (is.null(control$fixed_rbf_centers)) {
      stop(
        "threshold_control$fixed_rbf_centers is required for ",
        "threshold_learner='fixed_rbf'."
      )
    }
    return(fit_threshold_fixed_rbf(
      data = data,
      y = y,
      level = level,
      weight = weight,
      x_cols = x_cols,
      centers = control$fixed_rbf_centers,
      bandwidth = control$fixed_rbf_bandwidth,
      ridge_numerator = control$fixed_rbf_ridge_numerator,
      rho = control$fixed_rbf_admm_rho,
      maxit = control$fixed_rbf_admm_maxit,
      abs_tol = control$fixed_rbf_admm_abs_tol,
      rel_tol = control$fixed_rbf_admm_rel_tol
    ))
  }
  if (learner == "growing_rbf_sieve") {
    return(fit_threshold_growing_rbf_cv(
      data = data,
      y = y,
      level = level,
      weight = weight,
      x_cols = x_cols,
      bandwidth_grid = control$growing_rbf_bandwidth_grid,
      rank_multiplier = control$growing_rbf_rank_multiplier,
      rank_exponent = control$growing_rbf_rank_exponent,
      minimum_rank = control$growing_rbf_min_rank,
      center_skip = control$growing_rbf_center_skip,
      bandwidth_rule = control$growing_rbf_bandwidth_rule,
      bandwidth_sample = control$growing_rbf_bandwidth_sample,
      ridge_numerator = control$growing_rbf_ridge_numerator,
      cv_folds = control$cv_folds,
      selection_rule = control$selection_rule,
      rho = control$growing_rbf_admm_rho,
      maxit = control$growing_rbf_admm_maxit,
      abs_tol = control$growing_rbf_admm_abs_tol,
      rel_tol = control$growing_rbf_admm_rel_tol,
      plan = rkhs_plan,
      seed = seed
    ))
  }
  fit_threshold_rkhs_cv(
    data = data,
    y = y,
    level = level,
    weight = weight,
    x_cols = x_cols,
    backend = control$rkhs_backend,
    lambda_grid = control$rkhs_lambda_grid,
    bandwidth_multipliers =
      control$rkhs_bandwidth_multipliers,
    rank = control$rkhs_rank,
    cv_folds = control$cv_folds,
    selection_rule = control$selection_rule,
    bandwidth_sample = control$rkhs_bandwidth_sample,
    admm_rho = control$rkhs_admm_rho,
    admm_maxit = control$rkhs_admm_maxit,
    admm_abs_tol = control$rkhs_admm_abs_tol,
    admm_rel_tol = control$rkhs_admm_rel_tol,
    admm_check_interval = control$rkhs_admm_check_interval,
    condition_tolerance = control$rkhs_condition_tolerance,
    numerical_jitter = control$rkhs_numerical_jitter,
    plan = rkhs_plan,
    seed = seed
  )
}

# =============================================================================
# 9. Projected envelope formulas and compatibility stabilization
# =============================================================================

compute_compatible_gamma_interval <- function(
    nuisance,
    t_value,
    Lambda0,
    Lambda1
) {
  pi0 <- nuisance$pi0
  pi1 <- nuisance$pi1

  if (t_value == 1) {
    gamma_lower <- (
      (1 - pi0) / Lambda0 + pi0
    ) / (
      (1 - pi1) * Lambda1 + pi1
    )

    gamma_upper <- (
      Lambda0 * (1 - pi0) + pi0
    ) / (
      (1 - pi1) / Lambda1 + pi1
    )
  } else {
    gamma_lower <- (
      pi0 / Lambda0 + 1 - pi0
    ) / (
      pi1 * Lambda1 + 1 - pi1
    )

    gamma_upper <- (
      Lambda0 * pi0 + 1 - pi0
    ) / (
      pi1 / Lambda1 + 1 - pi1
    )
  }

  list(
    lower = gamma_lower,
    upper = gamma_upper
  )
}

project_gamma_to_compatible_interval <- function(
    gamma,
    nuisance,
    t_value,
    Lambda0,
    Lambda1,
    margin_fraction = 1e-4
) {
  interval <- compute_compatible_gamma_interval(
    nuisance = nuisance,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1
  )

  width <- interval$upper - interval$lower

  if (
    any(!is.finite(width)) ||
      any(width <= 0)
  ) {
    stop(
      sprintf(
        "The compatible gamma interval is empty for treatment %d.",
        t_value
      )
    )
  }

  interior_lower <- interval$lower +
    margin_fraction * width

  interior_upper <- interval$upper -
    margin_fraction * width

  projected <- pmin(
    interior_upper,
    pmax(
      interior_lower,
      gamma
    )
  )

  list(
    gamma = projected,
    original = gamma,
    lower = interval$lower,
    upper = interval$upper,
    projected = abs(projected - gamma) >
      10 * .Machine$double.eps
  )
}

compute_projected_envelopes <- function(
    nuisance,
    gamma,
    t_value,
    Lambda0,
    Lambda1,
    envelope_smoothing_epsilon = 0,
    gamma_compatibility_action =
      c("stop", "project"),
    gamma_projection_margin = 1e-4
) {
  gamma_compatibility_action <- match.arg(
    gamma_compatibility_action
  )

  if (
    Lambda0 < 1 ||
      Lambda1 < 1
  ) {
    stop(
      "Lambda0 and Lambda1 must both be at least one."
    )
  }

  if (
    length(envelope_smoothing_epsilon) != 1L ||
      !is.finite(envelope_smoothing_epsilon) ||
      envelope_smoothing_epsilon < 0
  ) {
    stop("envelope_smoothing_epsilon must be one nonnegative finite number.")
  }
  if (
    envelope_smoothing_epsilon > 0 &&
      envelope_smoothing_epsilon >= 1 / Lambda0
  ) {
    stop(
      paste0(
        "envelope_smoothing_epsilon must be smaller than 1/Lambda0 so ",
        "the conservative smoothed lower envelope remains positive."
      )
    )
  }

  # At Lambda0 = Lambda1 = 1 the sensitivity class is a singleton:
  # lambda_t(y, x) = 1.  The compatible gamma interval also degenerates to
  # {1}, so the usual interior projection is undefined and unnecessary for
  # evaluating the no-unmeasured-confounding plug-in benchmark.  Preserve the
  # incompatibility indicator rather than pretending the observed law is
  # compatible with this benchmark.
  if (Lambda0 == 1 && Lambda1 == 1) {
    interval <- compute_compatible_gamma_interval(
      nuisance = nuisance,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1
    )

    incompatible_original <- !is.finite(gamma) |
      abs(gamma - 1) > 10 * .Machine$double.eps

    return(list(
      A = rep(NA_real_, length(gamma)),
      B = rep(NA_real_, length(gamma)),
      L = rep(1, length(gamma)),
      U = rep(1, length(gamma)),
      D = rep(0, length(gamma)),
      hard_L = rep(1, length(gamma)),
      hard_U = rep(1, length(gamma)),
      candidate_lower = rep(1, length(gamma)),
      candidate_upper = rep(1, length(gamma)),
      lower_smoothing_gap = rep(0, length(gamma)),
      upper_smoothing_gap = rep(0, length(gamma)),
      envelope_smoothing_epsilon = 0,
      envelope_smoothing_bandwidth = 0,
      gamma_original = gamma,
      gamma_used = gamma,
      gamma_compatible_lower = interval$lower,
      gamma_compatible_upper = interval$upper,
      gamma_incompatible_original = incompatible_original,
      gamma_projected = rep(FALSE, length(gamma))
    ))
  }

  interval <- compute_compatible_gamma_interval(
    nuisance = nuisance,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1
  )

  incompatible_original <- (
    gamma < interval$lower |
      gamma > interval$upper |
      !is.finite(gamma)
  )

  if (
    gamma_compatibility_action == "project"
  ) {
    gamma_projection <-
      project_gamma_to_compatible_interval(
        gamma = gamma,
        nuisance = nuisance,
        t_value = t_value,
        Lambda0 = Lambda0,
        Lambda1 = Lambda1,
        margin_fraction =
          gamma_projection_margin
      )

    gamma_used <- gamma_projection$gamma
    gamma_projected <-
      gamma_projection$projected
  } else {
    gamma_used <- gamma
    gamma_projected <- rep(
      FALSE,
      length(gamma)
    )
  }

  pi0 <- nuisance$pi0
  pi1 <- nuisance$pi1

  if (t_value == 1) {
    A <- (
      (1 - pi1) /
        (1 - pi0)
    ) * gamma_used

    B <- (
      pi1 * gamma_used -
        pi0
    ) / (
      1 - pi0
    )
  } else {
    A <- (
      pi1 /
        pi0
    ) * gamma_used

    B <- (
      (1 - pi1) * gamma_used -
        (1 - pi0)
    ) / pi0
  }

  candidate_lower <-
    A / Lambda1 + B
  candidate_upper <-
    A * Lambda1 + B

  hard_L <- pmax(
    1 / Lambda0,
    candidate_lower
  )

  hard_U <- pmin(
    Lambda0,
    candidate_upper
  )

  smoothing_bandwidth <- if (envelope_smoothing_epsilon > 0) {
    envelope_smoothing_epsilon / log(2)
  } else {
    0
  }

  if (envelope_smoothing_epsilon > 0) {
    # Shifted log-sum-exp/log-sum-exp-of-negatives smooth the max/min in the
    # conservative direction.  Pointwise,
    #
    #   L_smooth <= hard_L <= hard_U <= U_smooth,
    #
    # and each one-sided discrepancy is at most
    # envelope_smoothing_epsilon = smoothing_bandwidth * log(2).
    smooth_lower <- function(first, second, bandwidth) {
      maximum <- pmax(first, second)
      value <- maximum + bandwidth * log(
        (
          exp((first - maximum) / bandwidth) +
            exp((second - maximum) / bandwidth)
        ) / 2
      )
      value[is.infinite(maximum) & maximum > 0] <- Inf
      value
    }
    smooth_upper <- function(first, second, bandwidth) {
      minimum <- pmin(first, second)
      value <- minimum - bandwidth * log(
        (
          exp(-(first - minimum) / bandwidth) +
            exp(-(second - minimum) / bandwidth)
        ) / 2
      )
      value[is.infinite(minimum) & minimum < 0] <- -Inf
      value
    }
    L <- smooth_lower(
      1 / Lambda0,
      candidate_lower,
      smoothing_bandwidth
    )
    U <- smooth_upper(
      Lambda0,
      candidate_upper,
      smoothing_bandwidth
    )
  } else {
    L <- hard_L
    U <- hard_U
  }

  list(
    A = A,
    B = B,
    L = L,
    U = U,
    D = U - L,
    hard_L = hard_L,
    hard_U = hard_U,
    candidate_lower = candidate_lower,
    candidate_upper = candidate_upper,
    lower_smoothing_gap = hard_L - L,
    upper_smoothing_gap = U - hard_U,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    envelope_smoothing_bandwidth = smoothing_bandwidth,
    gamma_original = gamma,
    gamma_used = gamma_used,
    gamma_compatible_lower =
      interval$lower,
    gamma_compatible_upper =
      interval$upper,
    gamma_incompatible_original =
      incompatible_original,
    gamma_projected =
      gamma_projected
  )
}

compute_target_odds <- function(
    nuisance,
    t_value
) {
  if (t_value == 1) {
    (1 - nuisance$pi0) /
      nuisance$pi0
  } else {
    nuisance$pi0 /
      (1 - nuisance$pi0)
  }
}

evaluate_threshold_bound_contributions <- function(
    data,
    nuisance,
    envelope,
    t_value,
    tau,
    q_lower,
    q_upper
) {
  n <- nrow(data)
  input_lengths <- c(
    n,
    length(nuisance$e),
    length(tau),
    length(q_lower),
    length(q_upper),
    length(envelope$L),
    length(envelope$U),
    length(envelope$D)
  )
  if (length(unique(input_lengths)) != 1L) {
    stop("Threshold-bound contribution inputs must have a common length.")
  }
  baseline <- data$T_obs == t_value & data$Z == 0
  odds <- compute_target_odds(nuisance, t_value)
  inverse_instrument <- 1 / nuisance$e

  lower_below <- as.numeric(data$Y <= q_lower)
  upper_below <- as.numeric(data$Y <= q_upper)
  lambda_lower <- envelope$L + envelope$D * lower_below
  lambda_upper <- envelope$U - envelope$D * upper_below

  lower_level <- 1 - tau
  upper_level <- tau
  interpolated_lambda <-
    tau * envelope$L + (1 - tau) * envelope$U

  direct_lower <- direct_upper <- pseudo_lower <- pseudo_upper <-
    lambda_one <- numeric(n)
  lower_score_residual <- upper_score_residual <- numeric(n)

  direct_lower[baseline] <- inverse_instrument[baseline] *
    data$Y[baseline] *
    (1 + odds[baseline] * lambda_lower[baseline])
  direct_upper[baseline] <- inverse_instrument[baseline] *
    data$Y[baseline] *
    (1 + odds[baseline] * lambda_upper[baseline])
  lambda_one[baseline] <- inverse_instrument[baseline] *
    data$Y[baseline] * (1 + odds[baseline])

  pseudo_lower[baseline] <- inverse_instrument[baseline] * (
    data$Y[baseline] * (
      1 + odds[baseline] * interpolated_lambda[baseline]
    ) -
      odds[baseline] * envelope$D[baseline] *
        check_loss(
          data$Y[baseline],
          q_lower[baseline],
          lower_level[baseline]
        )
  )
  pseudo_upper[baseline] <- inverse_instrument[baseline] * (
    data$Y[baseline] * (
      1 + odds[baseline] * interpolated_lambda[baseline]
    ) +
      odds[baseline] * envelope$D[baseline] *
        check_loss(
          data$Y[baseline],
          q_upper[baseline],
          upper_level[baseline]
        )
  )

  # These are the intercept score equations multiplied by q.  Algebraically
  # they equal direct minus pseudo row by row; their sample means therefore
  # diagnose whether a fitted restricted threshold has attained the
  # normalization identity needed to evaluate the original direct functional.
  lower_score_residual[baseline] <- inverse_instrument[baseline] *
    odds[baseline] * envelope$D[baseline] * q_lower[baseline] *
    (lower_below[baseline] - lower_level[baseline])
  upper_score_residual[baseline] <- -inverse_instrument[baseline] *
    odds[baseline] * envelope$D[baseline] * q_upper[baseline] *
    (upper_below[baseline] - upper_level[baseline])

  list(
    baseline = baseline,
    lambda_lower = lambda_lower,
    lambda_upper = lambda_upper,
    direct_lower = direct_lower,
    direct_upper = direct_upper,
    pseudo_lower = pseudo_lower,
    pseudo_upper = pseudo_upper,
    lambda_one = lambda_one,
    lower_score_residual = lower_score_residual,
    upper_score_residual = upper_score_residual
  )
}

compute_eif_envelope_derivatives <- function(
    nuisance,
    gamma,
    t_value,
    Lambda0,
    Lambda1,
    envelope_smoothing_epsilon = 0
) {
  if (Lambda0 <= 1 || Lambda1 <= 1) {
    stop(
      paste0(
        "The generic EIF envelope derivative requires ",
        "Lambda0 > 1 and Lambda1 > 1."
      )
    )
  }

  r0 <- if (t_value == 1) {
    nuisance$pi0
  } else {
    1 - nuisance$pi0
  }

  r1 <- if (t_value == 1) {
    nuisance$pi1
  } else {
    1 - nuisance$pi1
  }

  c_lower <- 1 / Lambda1
  c_upper <- Lambda1
  h_lower <- r1 + (1 - r1) * c_lower
  h_upper <- r1 + (1 - r1) * c_upper

  candidate_lower <- (
    gamma * h_lower - r0
  ) / (
    1 - r0
  )

  candidate_upper <- (
    gamma * h_upper - r0
  ) / (
    1 - r0
  )

  active_lower <- candidate_lower > 1 / Lambda0
  active_upper <- candidate_upper < Lambda0

  if (
    length(envelope_smoothing_epsilon) != 1L ||
      !is.finite(envelope_smoothing_epsilon) ||
      envelope_smoothing_epsilon < 0
  ) {
    stop("envelope_smoothing_epsilon must be one nonnegative finite number.")
  }

  if (envelope_smoothing_epsilon > 0) {
    smoothing_bandwidth <- envelope_smoothing_epsilon / log(2)
    lower_multiplier <- stats::plogis(
      (candidate_lower - 1 / Lambda0) / smoothing_bandwidth
    )
    upper_multiplier <- stats::plogis(
      (Lambda0 - candidate_upper) / smoothing_bandwidth
    )
  } else {
    smoothing_bandwidth <- 0
    lower_multiplier <- as.numeric(active_lower)
    upper_multiplier <- as.numeric(active_upper)
  }

  derivative_candidate <- function(
      c_value,
      h_value,
      multiplier
  ) {
    list(
      gamma = multiplier *
        h_value / (1 - r0),
      r0 = multiplier *
        (gamma * h_value - 1) /
        (1 - r0)^2,
      r1 = multiplier *
        gamma * (1 - c_value) /
        (1 - r0)
    )
  }

  list(
    r0 = r0,
    r1 = r1,
    candidate_lower = candidate_lower,
    candidate_upper = candidate_upper,
    active_lower = active_lower,
    active_upper = active_upper,
    lower_derivative_weight = lower_multiplier,
    upper_derivative_weight = upper_multiplier,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    envelope_smoothing_bandwidth = smoothing_bandwidth,
    lower = derivative_candidate(
      c_value = c_lower,
      h_value = h_lower,
      multiplier = lower_multiplier
    ),
    upper = derivative_candidate(
      c_value = c_upper,
      h_value = h_upper,
      multiplier = upper_multiplier
    )
  )
}

compute_eif_row_primitives <- function(
    data,
    nuisance,
    envelope,
    t_value,
    Lambda0,
    Lambda1,
    q_lower,
    q_upper
) {
  required_lengths <- c(
    nrow(data),
    length(nuisance$pi0),
    length(nuisance$pi1),
    length(envelope$L),
    length(envelope$U),
    length(envelope$gamma_used),
    length(q_lower),
    length(q_upper)
  )

  if (length(unique(required_lengths)) != 1L) {
    stop("EIF row inputs must have a common length.")
  }

  if (
    any(!is.finite(q_lower)) ||
      any(!is.finite(q_upper))
  ) {
    stop("EIF threshold predictions must be finite.")
  }

  derivative <- compute_eif_envelope_derivatives(
    nuisance = nuisance,
    gamma = envelope$gamma_used,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    envelope_smoothing_epsilon = if (
      is.null(envelope$envelope_smoothing_epsilon)
    ) {
      0
    } else {
      envelope$envelope_smoothing_epsilon
    }
  )

  make_endpoint <- function(
      q,
      direction = c("lower", "upper")
  ) {
    direction <- match.arg(direction)

    allocate_upper <- if (direction == "lower") {
      data$Y <= q
    } else {
      data$Y > q
    }

    allocate_upper <- as.numeric(
      allocate_upper
    )

    allocate_lower <- 1 - allocate_upper

    lambda <- allocate_upper * envelope$U +
      allocate_lower * envelope$L

    envelope_derivative <- function(name) {
      (data$Y - q) * (
        allocate_upper * derivative$upper[[name]] +
          allocate_lower * derivative$lower[[name]]
      )
    }

    D_gamma <- envelope_derivative("gamma")
    D_r0 <- envelope_derivative("r0")
    D_r1 <- envelope_derivative("r1")

    H0 <- derivative$r0 * data$Y +
      (1 - derivative$r0) * (
        lambda * (data$Y - q) -
          envelope$gamma_used * D_gamma
      )

    H1 <- (1 - derivative$r0) * D_gamma

    data.frame(
      q = q,
      lambda = lambda,
      D_gamma = D_gamma,
      D_r0 = D_r0,
      D_r1 = D_r1,
      H0 = H0,
      H1 = H1
    )
  }

  list(
    r0 = derivative$r0,
    r1 = derivative$r1,
    active_lower = derivative$active_lower,
    active_upper = derivative$active_upper,
    lower_derivative_weight = derivative$lower_derivative_weight,
    upper_derivative_weight = derivative$upper_derivative_weight,
    lower = make_endpoint(q_lower, "lower"),
    upper = make_endpoint(q_upper, "upper")
  )
}

compute_fixed_q_linked_components <- function(
    m0,
    Hbar0,
    Hbar1,
    q,
    r0,
    mean_D_r0
) {
  input_lengths <- lengths(list(
    m0 = m0,
    Hbar0 = Hbar0,
    Hbar1 = Hbar1,
    q = q,
    r0 = r0,
    mean_D_r0 = mean_D_r0
  ))
  if (length(unique(input_lengths)) != 1L || input_lengths[1L] == 0L) {
    stop("Fixed-q linked components must have a common positive length.")
  }
  inputs <- c(m0, Hbar0, Hbar1, q, r0, mean_D_r0)
  if (any(!is.finite(inputs)) || any(r0 <= 0) || any(r0 >= 1)) {
    stop("Fixed-q linked component inputs must be finite with 0 < r0 < 1.")
  }

  # For an arbitrary held-fixed threshold,
  #
  #   theta_q = q + E_0[lambda_q(Y - q) | X],
  #   M_q     = r0 m0 + (1-r0) theta_q
  #           = Hbar0 + Hbar1 + (1-r0) q.
  #
  # This differs from E_0(lambda_q Y | X) unless lambda_q is normalized.
  M <- Hbar0 + Hbar1 + (1 - r0) * q
  theta_fixed_q <- (M - r0 * m0) / (1 - r0)
  C0 <- m0 - theta_fixed_q + (1 - r0) * mean_D_r0

  list(
    M = M,
    theta_fixed_q = theta_fixed_q,
    C0 = C0
  )
}

fit_eif_arm_fold <- function(
    train,
    validation,
    row_id_validation,
    nuisance_train,
    nuisance_validation,
    envelope_train,
    envelope_validation,
    t_value,
    Lambda0,
    Lambda1,
    q_lower_train,
    q_upper_train,
    q_lower_validation,
    q_upper_validation,
    x_cols = "X",
    conditional_spline_df = 5,
    conditional_learner = c(
      "spline_fixed", "spline_cv", "gam", "ranger", "library_cv"
    ),
    conditional_control = NULL,
    conditional_seed = 1L,
    conditional_prediction_provider = NULL,
    instrument_probability_validation = NULL
) {
  conditional_learner <- match.arg(conditional_learner)
  if (
    !is.null(conditional_prediction_provider) &&
      !is.function(conditional_prediction_provider)
  ) {
    stop("conditional_prediction_provider must be NULL or a function.")
  }
  resolved_conditional_control <- resolve_eif_conditional_control(
    conditional_control
  )
  if (is.null(instrument_probability_validation)) {
    instrument_probability_validation <- nuisance_validation$e
  }
  if (
    length(instrument_probability_validation) != nrow(validation) ||
      any(!is.finite(instrument_probability_validation)) ||
      any(instrument_probability_validation <= 0) ||
      any(instrument_probability_validation >= 1)
  ) {
    stop(
      paste0(
        "instrument_probability_validation must contain one finite value ",
        "strictly between zero and one per validation row."
      )
    )
  }
  primitive_train <- compute_eif_row_primitives(
    data = train,
    nuisance = nuisance_train,
    envelope = envelope_train,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    q_lower = q_lower_train,
    q_upper = q_upper_train
  )

  primitive_validation <- compute_eif_row_primitives(
    data = validation,
    nuisance = nuisance_validation,
    envelope = envelope_validation,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    q_lower = q_lower_validation,
    q_upper = q_upper_validation
  )

  baseline_train <- train$T_obs == t_value &
    train$Z == 0

  comparison_train <- train$T_obs == t_value &
    train$Z == 1

  learner_complexity <- switch(
    conditional_learner,
    spline_fixed = conditional_spline_df,
    spline_cv = max(resolved_conditional_control$spline_df_grid),
    gam = resolved_conditional_control$gam_k,
    ranger = 2L,
    library_cv = max(vapply(
      resolved_conditional_control$library_candidates,
      function(candidate) switch(
        candidate,
        linear = 1L,
        ranger = 2L,
        gam = resolved_conditional_control$gam_k
      ),
      integer(1)
    ))
  )
  minimum_size <- learner_complexity *
    length(x_cols) + 3L

  if (
    sum(baseline_train) <= minimum_size ||
      sum(comparison_train) <= minimum_size
  ) {
    stop(
      sprintf(
        paste0(
          "Too few observations to fit EIF conditional means for ",
          "treatment %d. Baseline=%d; comparison=%d."
        ),
        t_value,
        sum(baseline_train),
        sum(comparison_train)
      )
    )
  }

  baseline_data <- train[
    baseline_train,
    ,
    drop = FALSE
  ]

  baseline_data$.m0 <- baseline_data$Y

  for (endpoint in c("lower", "upper")) {
    primitive <- primitive_train[[endpoint]][
      baseline_train,
      ,
      drop = FALSE
    ]

    baseline_data[[paste0(".theta_", endpoint)]] <-
      primitive$lambda * baseline_data$Y
    baseline_data[[paste0(".D_r0_", endpoint)]] <-
      primitive$D_r0
    baseline_data[[paste0(".D_r1_", endpoint)]] <-
      primitive$D_r1
    baseline_data[[paste0(".H0_", endpoint)]] <-
      primitive$H0
  }

  baseline_outcomes <- c(
    ".m0",
    unlist(lapply(
      c("lower", "upper"),
      function(endpoint) {
        paste0(
          c(".theta_", ".D_r0_", ".D_r1_", ".H0_"),
          endpoint
        )
      }
    ))
  )

  baseline_fit <- if (is.null(conditional_prediction_provider)) {
    fit_eif_conditional_means(
      data = baseline_data,
      outcomes = baseline_outcomes,
      x_cols = x_cols,
      learner = conditional_learner,
      spline_df = conditional_spline_df,
      control = resolved_conditional_control,
      seed = conditional_seed
    )
  } else {
    NULL
  }

  comparison_data <- train[
    comparison_train,
    ,
    drop = FALSE
  ]

  for (endpoint in c("lower", "upper")) {
    comparison_data[[paste0(".H1_", endpoint)]] <-
      primitive_train[[endpoint]]$H1[
        comparison_train
      ]
  }

  comparison_outcomes <- paste0(
    ".H1_",
    c("lower", "upper")
  )

  comparison_fit <- if (is.null(conditional_prediction_provider)) {
    fit_eif_conditional_means(
      data = comparison_data,
      outcomes = comparison_outcomes,
      x_cols = x_cols,
      learner = conditional_learner,
      spline_df = conditional_spline_df,
      control = resolved_conditional_control,
      seed = conditional_seed + 1L
    )
  } else {
    NULL
  }

  baseline_prediction <- if (is.null(baseline_fit)) {
    NULL
  } else {
    baseline_fit$predict(validation)
  }

  comparison_prediction <- if (is.null(comparison_fit)) {
    NULL
  } else {
    comparison_fit$predict(validation)
  }

  provider_output <- NULL
  if (!is.null(conditional_prediction_provider)) {
    provider_output <- conditional_prediction_provider(
      validation = validation,
      nuisance_validation = nuisance_validation,
      envelope_validation = envelope_validation,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      q_lower_validation = q_lower_validation,
      q_upper_validation = q_upper_validation,
      baseline_outcomes = baseline_outcomes,
      comparison_outcomes = comparison_outcomes
    )
    if (
      !is.list(provider_output) ||
        !all(
          c("baseline_prediction", "comparison_prediction") %in%
            names(provider_output)
        )
    ) {
      stop(
        paste0(
          "conditional_prediction_provider must return baseline_prediction ",
          "and comparison_prediction matrices."
        )
      )
    }
    baseline_prediction <- as.matrix(provider_output$baseline_prediction)
    comparison_prediction <- as.matrix(
      provider_output$comparison_prediction
    )
    if (
      nrow(baseline_prediction) != nrow(validation) ||
        nrow(comparison_prediction) != nrow(validation) ||
        !all(baseline_outcomes %in% colnames(baseline_prediction)) ||
        !all(comparison_outcomes %in% colnames(comparison_prediction)) ||
        any(!is.finite(baseline_prediction[, baseline_outcomes])) ||
        any(!is.finite(comparison_prediction[, comparison_outcomes]))
    ) {
      stop(
        paste0(
          "The conditional prediction provider returned invalid dimensions, ",
          "columns, or non-finite values."
        )
      )
    }
  }

  observed_arm <- as.numeric(
    validation$T_obs == t_value
  )

  instrument0 <- as.numeric(
    validation$Z == 0
  )

  instrument1 <- 1 - instrument0

  make_signal <- function(endpoint) {
    primitive <- primitive_validation[[endpoint]]
    m0 <- baseline_prediction[, ".m0"]
    theta_direct <- baseline_prediction[,
      paste0(".theta_", endpoint)
    ]

    Hbar0 <- baseline_prediction[,
      paste0(".H0_", endpoint)
    ]

    Hbar1 <- comparison_prediction[,
      paste0(".H1_", endpoint)
    ]

    # The exact fixed-q conditional identity
    #
    #   M_q = E_0(H0 | X) + E_1(H1 | X) + (1-r0) q
    #       = r0 m0 + (1-r0) theta_q
    #
    # must also be enforced for fitted conditional regressions. Estimating M
    # independently from the direct moment E_0(lambda_q Y | X) is valid only
    # at a normalized optimizer. Linking both M_q and theta_q to the Hbar
    # regressions gives the fixed-q pseudo-bound EIF even when the training-fold
    # threshold has a nonzero normalization residual.
    fixed_q <- compute_fixed_q_linked_components(
      m0 = m0,
      Hbar0 = Hbar0,
      Hbar1 = Hbar1,
      q = primitive$q,
      r0 = primitive_validation$r0,
      mean_D_r0 = baseline_prediction[,
        paste0(".D_r0_", endpoint)
      ]
    )
    M <- fixed_q$M
    C0 <- fixed_q$C0

    C1 <- (1 - primitive_validation$r0) *
      baseline_prediction[,
        paste0(".D_r1_", endpoint)
      ]

    treatment_correction <-
      instrument0 / instrument_probability_validation *
      C0 * (observed_arm - primitive_validation$r0) +
      instrument1 / (1 - instrument_probability_validation) *
      C1 * (observed_arm - primitive_validation$r1)

    outcome_correction <-
      instrument0 / instrument_probability_validation *
      observed_arm / primitive_validation$r0 *
      (primitive$H0 - Hbar0) +
      instrument1 / (1 - instrument_probability_validation) *
      observed_arm / primitive_validation$r1 *
      (primitive$H1 - Hbar1)

    data.frame(
      q = primitive$q,
      M = M,
      theta_direct = theta_direct,
      theta_fixed_q = fixed_q$theta_fixed_q,
      C0 = C0,
      C1 = C1,
      Hbar0 = Hbar0,
      Hbar1 = Hbar1,
      treatment_correction = treatment_correction,
      outcome_correction = outcome_correction,
      signal = M +
        treatment_correction +
        outcome_correction
    )
  }

  lower_signal <- make_signal("lower")
  upper_signal <- make_signal("upper")

  diagnostic <- data.frame(
    treatment = t_value,
    active_lower_rate_train = mean(
      primitive_train$active_lower
    ),
    active_upper_rate_train = mean(
      primitive_train$active_upper
    ),
    active_lower_rate_validation = mean(
      primitive_validation$active_lower
    ),
    active_upper_rate_validation = mean(
      primitive_validation$active_upper
    ),
    lower_mean_treatment_correction = mean(
      lower_signal$treatment_correction
    ),
    lower_mean_outcome_correction = mean(
      lower_signal$outcome_correction
    ),
    upper_mean_treatment_correction = mean(
      upper_signal$treatment_correction
    ),
    upper_mean_outcome_correction = mean(
      upper_signal$outcome_correction
    ),
    lower_mean_absolute_fixed_q_theta_adjustment = mean(abs(
      lower_signal$theta_fixed_q - lower_signal$theta_direct
    )),
    upper_mean_absolute_fixed_q_theta_adjustment = mean(abs(
      upper_signal$theta_fixed_q - upper_signal$theta_direct
    )),
    baseline_regression_rank = if (is.null(provider_output)) {
      baseline_fit$rank
    } else {
      NA_real_
    },
    comparison_regression_rank = if (is.null(provider_output)) {
      comparison_fit$rank
    } else {
      NA_real_
    },
    conditional_learner = if (
      !is.null(provider_output) &&
        !is.null(provider_output$learner)
    ) {
      as.character(provider_output$learner)[1L]
    } else {
      conditional_learner
    },
    baseline_selected_spline_df = if (is.null(provider_output)) {
      baseline_fit$selected_spline_df
    } else {
      NA_integer_
    },
    comparison_selected_spline_df = if (is.null(provider_output)) {
      comparison_fit$selected_spline_df
    } else {
      NA_integer_
    },
    baseline_effective_df = if (is.null(provider_output)) {
      baseline_fit$effective_df
    } else {
      NA_real_
    },
    comparison_effective_df = if (is.null(provider_output)) {
      comparison_fit$effective_df
    } else {
      NA_real_
    },
    mean_instrument_probability_validation = mean(
      instrument_probability_validation
    ),
    mean_absolute_instrument_probability_perturbation = mean(abs(
      instrument_probability_validation - nuisance_validation$e
    )),
    maximum_absolute_instrument_probability_perturbation = max(abs(
      instrument_probability_validation - nuisance_validation$e
    ))
  )
  if (
    !is.null(provider_output) &&
      !is.null(provider_output$diagnostics)
  ) {
    provider_diagnostics <- as.data.frame(provider_output$diagnostics)
    if (nrow(provider_diagnostics) != 1L) {
      stop("Provider diagnostics must contain exactly one row.")
    }
    duplicate_names <- intersect(
      names(diagnostic), names(provider_diagnostics)
    )
    if (length(duplicate_names) > 0L) {
      stop(
        "Provider diagnostics duplicate core columns: ",
        paste(duplicate_names, collapse = ", ")
      )
    }
    diagnostic <- cbind(diagnostic, provider_diagnostics)
  }

  list(
    signals = data.frame(
      row_id = row_id_validation,
      treatment = t_value,
      r0 = primitive_validation$r0,
      lower_q = lower_signal$q,
      lower_M = lower_signal$M,
      lower_theta_direct = lower_signal$theta_direct,
      lower_theta_fixed_q = lower_signal$theta_fixed_q,
      lower_C0 = lower_signal$C0,
      lower_Hbar0 = lower_signal$Hbar0,
      lower_Hbar1 = lower_signal$Hbar1,
      lower_treatment_correction =
        lower_signal$treatment_correction,
      lower_outcome_correction =
        lower_signal$outcome_correction,
      lower_signal = lower_signal$signal,
      upper_q = upper_signal$q,
      upper_M = upper_signal$M,
      upper_theta_direct = upper_signal$theta_direct,
      upper_theta_fixed_q = upper_signal$theta_fixed_q,
      upper_C0 = upper_signal$C0,
      upper_Hbar0 = upper_signal$Hbar0,
      upper_Hbar1 = upper_signal$Hbar1,
      upper_treatment_correction =
        upper_signal$treatment_correction,
      upper_outcome_correction =
        upper_signal$outcome_correction,
      upper_signal = upper_signal$signal
    ),
    diagnostics = diagnostic,
    conditional_cv = if (is.null(provider_output)) {
      bind_conditional_cv_tables(
        baseline_fit$cv_table,
        comparison_fit$cv_table
      )
    } else {
      data.frame()
    }
  )
}

# =============================================================================
# 10. One treatment-arm bound fit in one outer fold
# =============================================================================

fit_bound_arm_fold <- function(
    train,
    validation,
    row_id_validation,
    nuisance_train,
    nuisance_validation,
    gamma_train,
    gamma_validation,
    t_value,
    Lambda0,
    Lambda1,
    mode,
    x_cols = "X",
    conditional_mean_spline_df = 5,
    threshold_spline_df = 5,
    threshold_learner = c(
      "spline_fixed", "spline_cv", "rkhs", "fixed_rbf",
      "growing_rbf_sieve", "external"
    ),
    threshold_control = NULL,
    threshold_seed = 1L,
    envelope_smoothing_epsilon = 0,
    compatibility_tolerance = 1e-6,
    compatibility_action = c("stop", "project"),
    gamma_projection_margin = 1e-4,
    max_gamma_projection_rate = 0.02,
    tau_epsilon = 1e-4,
    threshold_maxit = 4000,
    compute_eif = FALSE,
    allow_projected_eif = FALSE,
    eif_conditional_spline_df =
      conditional_mean_spline_df,
    eif_conditional_learner = c(
      "spline_fixed", "spline_cv", "gam", "ranger", "library_cv"
    ),
    eif_conditional_control = NULL,
    eif_conditional_seed = threshold_seed + 20000L,
    eif_conditional_prediction_provider = NULL,
    eif_instrument_probability_validation = NULL,
    bound_component_provider = NULL,
    fixed_threshold_target_evaluator = NULL
) {
  threshold_learner <- match.arg(threshold_learner)
  eif_conditional_learner <- match.arg(eif_conditional_learner)
  compatibility_action <- match.arg(
    compatibility_action
  )
  if (
    !is.null(bound_component_provider) &&
      !is.list(bound_component_provider)
  ) {
    stop(
      "bound_component_provider must be NULL or a list."
    )
  }
  if (
    !is.null(bound_component_provider$envelope_means) &&
      !is.function(bound_component_provider$envelope_means)
  ) {
    stop("bound_component_provider$envelope_means must be a function.")
  }
  if (
    !is.null(bound_component_provider$thresholds) &&
      !is.function(bound_component_provider$thresholds)
  ) {
    stop("bound_component_provider$thresholds must be a function.")
  }
  if (
    threshold_learner == "external" &&
      (is.null(bound_component_provider) ||
        !is.function(bound_component_provider$thresholds))
  ) {
    stop(
      paste0(
        "threshold_learner='external' requires ",
        "bound_component_provider$thresholds."
      )
    )
  }
  if (
    !is.null(fixed_threshold_target_evaluator) &&
      !is.function(fixed_threshold_target_evaluator)
  ) {
    stop("fixed_threshold_target_evaluator must be NULL or a function.")
  }

  envelope_train <- compute_projected_envelopes(
    nuisance = nuisance_train,
    gamma = gamma_train,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    gamma_compatibility_action =
      compatibility_action,
    gamma_projection_margin =
      gamma_projection_margin
  )

  envelope_validation <- compute_projected_envelopes(
    nuisance = nuisance_validation,
    gamma = gamma_validation,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    gamma_compatibility_action =
      compatibility_action,
    gamma_projection_margin =
      gamma_projection_margin
  )

  baseline_train <- (
    train$T_obs == t_value &
      train$Z == 0
  )

  baseline_validation <- (
    validation$T_obs == t_value &
      validation$Z == 0
  )

  pointwise_crossing_train <- (
    envelope_train$L >
      envelope_train$U +
      compatibility_tolerance
  )

  pointwise_crossing_validation <- (
    envelope_validation$L >
      envelope_validation$U +
      compatibility_tolerance
  )

  gamma_projection_rate_train <- mean(
    envelope_train$gamma_projected[
      baseline_train
    ]
  )

  gamma_projection_rate_validation <- mean(
    envelope_validation$gamma_projected[
      baseline_validation
    ]
  )

  gamma_incompatibility_rate_train <- mean(
    envelope_train$gamma_incompatible_original[
      baseline_train
    ]
  )

  gamma_incompatibility_rate_validation <- mean(
    envelope_validation$gamma_incompatible_original[
      baseline_validation
    ]
  )

  # This is a diagnostic flag, not a formal compatibility test. Estimated
  # gamma can require projection because of finite-sample nuisance noise,
  # especially in small cross-fitting strata. Returning the estimate and its
  # projection diagnostics avoids selectively dropping difficult Monte Carlo
  # replications. Pointwise oracle envelopes remain strict below.
  gamma_projection_threshold_exceeded <- (
    compatibility_action == "project" &&
      (
        gamma_projection_rate_train >
          max_gamma_projection_rate ||
          gamma_projection_rate_validation >
            max_gamma_projection_rate
      )
  )

  if (
    any(pointwise_crossing_train[baseline_train]) ||
      any(pointwise_crossing_validation[baseline_validation])
  ) {
    stop(
      sprintf(
        paste0(
          "Pointwise projected envelopes cross for mode %s, ",
          "treatment %d after the requested compatibility action. ",
          "Training rate %.4f; validation rate %.4f."
        ),
        mode,
        t_value,
        mean(pointwise_crossing_train[baseline_train]),
        mean(pointwise_crossing_validation[baseline_validation])
      )
    )
  }

  # The Lambda0 = Lambda1 = 1 class has no threshold optimization: both
  # extremizers equal one pointwise.  Handling the singleton directly avoids
  # sending zero weights to the check-loss optimizer and makes the reported
  # endpoints exactly equal to the lambda-one plug-in benchmark.
  if (Lambda0 == 1 && Lambda1 == 1) {
    target_odds_validation <- compute_target_odds(
      nuisance = nuisance_validation,
      t_value = t_value
    )

    contribution <- rep(0, nrow(validation))
    contribution[baseline_validation] <-
      validation$Y[baseline_validation] /
      nuisance_validation$e[baseline_validation] *
      (1 + target_odds_validation[baseline_validation])

    prediction_table <- data.frame(
      row_id = row_id_validation[baseline_validation],
      treatment = t_value,
      Y = validation$Y[baseline_validation],
      X = validation[baseline_validation, x_cols[1], drop = TRUE],
      L = 1,
      U = 1,
      D = 0,
      hard_L = 1,
      hard_U = 1,
      lower_smoothing_gap = 0,
      upper_smoothing_gap = 0,
      Lbar = 1,
      Ubar = 1,
      tau_raw = NA_real_,
      tau = NA_real_,
      q_lower = NA_real_,
      q_upper = NA_real_,
      lambda_lower = 1,
      lambda_upper = 1,
      direct_lower_contribution = contribution[baseline_validation],
      direct_upper_contribution = contribution[baseline_validation],
      pseudo_lower_contribution = contribution[baseline_validation],
      pseudo_upper_contribution = contribution[baseline_validation],
      lower_q_score_residual = 0,
      upper_q_score_residual = 0,
      threshold_learner = threshold_learner,
      threshold_rkhs_backend = NA_character_,
      gamma_original = envelope_validation$gamma_original[baseline_validation],
      gamma_used = envelope_validation$gamma_used[baseline_validation],
      gamma_compatible_lower = envelope_validation$gamma_compatible_lower[baseline_validation],
      gamma_compatible_upper = envelope_validation$gamma_compatible_upper[baseline_validation],
      gamma_projected = envelope_validation$gamma_projected[baseline_validation]
    )

    diagnostics <- data.frame(
      treatment = t_value,
      n_train_baseline = sum(baseline_train),
      n_validation_baseline = sum(baseline_validation),
      pointwise_crossing_rate_train = 0,
      pointwise_crossing_rate_validation = 0,
      gamma_incompatibility_rate_train = gamma_incompatibility_rate_train,
      gamma_incompatibility_rate_validation = gamma_incompatibility_rate_validation,
      gamma_projection_rate_train = 0,
      gamma_projection_rate_validation = 0,
      gamma_projection_threshold_exceeded = FALSE,
      envelope_smoothing_epsilon = 0,
      mean_lower_smoothing_gap_validation = 0,
      mean_upper_smoothing_gap_validation = 0,
      maximum_lower_smoothing_gap_validation = 0,
      maximum_upper_smoothing_gap_validation = 0,
      conditional_projection_rate_train = 0,
      conditional_projection_rate_validation = 0,
      tau_invalid_rate_train = NA_real_,
      tau_invalid_rate_validation = NA_real_,
      lower_threshold_convergence = NA_integer_,
      upper_threshold_convergence = NA_integer_,
      lower_threshold_objective = NA_real_,
      upper_threshold_objective = NA_real_,
      mean_L_validation = 1,
      mean_U_validation = 1,
      mean_tau_validation = NA_real_,
      mean_lambda_lower_validation = 1,
      mean_lambda_upper_validation = 1,
      threshold_learner = threshold_learner,
      threshold_rkhs_backend = NA_character_,
      lower_selected_spline_df = NA_integer_,
      upper_selected_spline_df = NA_integer_,
      lower_selected_lambda = NA_real_,
      upper_selected_lambda = NA_real_,
      lower_selected_bandwidth_multiplier = NA_real_,
      upper_selected_bandwidth_multiplier = NA_real_,
      lower_selected_bandwidth = NA_real_,
      upper_selected_bandwidth = NA_real_,
      lower_base_bandwidth = NA_real_,
      upper_base_bandwidth = NA_real_,
      lower_selected_rank = NA_integer_,
      upper_selected_rank = NA_integer_,
      lower_threshold_iterations = NA_real_,
      upper_threshold_iterations = NA_real_,
      lower_threshold_initial_convergence = NA_integer_,
      upper_threshold_initial_convergence = NA_integer_,
      lower_threshold_retry_used = NA,
      upper_threshold_retry_used = NA,
      lower_threshold_retry_accepted = NA,
      upper_threshold_retry_accepted = NA,
      lower_threshold_retry_convergence = NA_integer_,
      upper_threshold_retry_convergence = NA_integer_,
      lower_threshold_initial_penalized_objective = NA_real_,
      upper_threshold_initial_penalized_objective = NA_real_,
      lower_threshold_penalized_objective = NA_real_,
      upper_threshold_penalized_objective = NA_real_,
      lower_threshold_tie_break_ridge = NA_real_,
      upper_threshold_tie_break_ridge = NA_real_,
      lower_threshold_tie_break_ridge_numerator = NA_real_,
      upper_threshold_tie_break_ridge_numerator = NA_real_,
      lower_threshold_tie_break_ridge_training_size = NA_integer_,
      upper_threshold_tie_break_ridge_training_size = NA_integer_,
      lower_threshold_primal_residual = NA_real_,
      upper_threshold_primal_residual = NA_real_,
      lower_threshold_dual_residual = NA_real_,
      upper_threshold_dual_residual = NA_real_,
      lower_threshold_final_rho = NA_real_,
      upper_threshold_final_rho = NA_real_,
      lower_threshold_reciprocal_condition = NA_real_,
      upper_threshold_reciprocal_condition = NA_real_,
      lower_threshold_numerical_jitter = NA_real_,
      upper_threshold_numerical_jitter = NA_real_,
      lower_kernel_minimum_raw_eigenvalue = NA_real_,
      upper_kernel_minimum_raw_eigenvalue = NA_real_,
      fold_lower_mean = mean(contribution),
      fold_upper_mean = mean(contribution),
      fold_lambda_one_mean = mean(contribution),
      fold_lower_pseudo_mean = mean(contribution),
      fold_upper_pseudo_mean = mean(contribution),
      lower_direct_minus_pseudo_train = 0,
      upper_direct_minus_pseudo_train = 0,
      lower_direct_minus_pseudo_validation = 0,
      upper_direct_minus_pseudo_validation = 0,
      lower_q_score_residual_train = 0,
      upper_q_score_residual_train = 0,
      lower_q_score_residual_validation = 0,
      upper_q_score_residual_validation = 0,
      lower_direct_pseudo_identity_error_validation = 0,
      upper_direct_pseudo_identity_error_validation = 0
    )

    return(list(
      lower_contribution = contribution,
      upper_contribution = contribution,
      lower_pseudo_contribution = contribution,
      upper_pseudo_contribution = contribution,
      lambda_one_contribution = contribution,
      predictions = prediction_table,
      diagnostics = diagnostics,
      threshold_cv = data.frame(),
      fixed_threshold_target_evaluation = data.frame(),
      eif = NULL
    ))
  }

  if (
    sum(baseline_train) <=
      conditional_mean_spline_df + 2
  ) {
    stop(
      sprintf(
        "Too few T=%d,Z=0 training observations.",
        t_value
      )
    )
  }

  regression_data <- train[
    baseline_train,
    ,
    drop = FALSE
  ]

  regression_data$.L <- envelope_train$L[
    baseline_train
  ]

  regression_data$.U <- envelope_train$U[
    baseline_train
  ]

  Lbar_fit <- fit_spline_mean(
    data = regression_data,
    outcome = ".L",
    x_cols = x_cols,
    spline_df = conditional_mean_spline_df
  )

  Ubar_fit <- fit_spline_mean(
    data = regression_data,
    outcome = ".U",
    x_cols = x_cols,
    spline_df = conditional_mean_spline_df
  )

  Lbar_train <- Lbar_fit$predict(
    train
  )

  Ubar_train <- Ubar_fit$predict(
    train
  )

  Lbar_validation <- Lbar_fit$predict(
    validation
  )

  Ubar_validation <- Ubar_fit$predict(
    validation
  )

  if (
    !is.null(bound_component_provider) &&
      is.function(bound_component_provider$envelope_means)
  ) {
    provider_envelope_means <- bound_component_provider$envelope_means(
      train = train,
      validation = validation,
      nuisance_train = nuisance_train,
      nuisance_validation = nuisance_validation,
      envelope_train = envelope_train,
      envelope_validation = envelope_validation,
      treatment = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      compatibility_action = compatibility_action,
      gamma_projection_margin = gamma_projection_margin,
      envelope_smoothing_epsilon = envelope_smoothing_epsilon
    )
    required_envelope_means <- c(
      "Lbar_train", "Ubar_train", "Lbar_validation", "Ubar_validation"
    )
    if (
      !is.list(provider_envelope_means) ||
        !all(required_envelope_means %in% names(provider_envelope_means))
    ) {
      stop(
        paste0(
          "bound_component_provider$envelope_means must return ",
          paste(required_envelope_means, collapse = ", "), "."
        )
      )
    }
    provider_lengths <- vapply(
      provider_envelope_means[required_envelope_means], length, integer(1)
    )
    if (
      !identical(
        unname(provider_lengths),
        c(nrow(train), nrow(train), nrow(validation), nrow(validation))
      ) ||
        any(!is.finite(unlist(provider_envelope_means[required_envelope_means])))
    ) {
      stop("The external envelope-mean provider returned invalid predictions.")
    }
    Lbar_train <- as.numeric(provider_envelope_means$Lbar_train)
    Ubar_train <- as.numeric(provider_envelope_means$Ubar_train)
    Lbar_validation <- as.numeric(provider_envelope_means$Lbar_validation)
    Ubar_validation <- as.numeric(provider_envelope_means$Ubar_validation)
  }

  invalid_train <- (
    Lbar_train > 1 + compatibility_tolerance |
      Ubar_train < 1 - compatibility_tolerance |
      Ubar_train - Lbar_train <=
        compatibility_tolerance
  )

  invalid_validation <- (
    Lbar_validation > 1 + compatibility_tolerance |
      Ubar_validation < 1 - compatibility_tolerance |
      Ubar_validation - Lbar_validation <=
        compatibility_tolerance
  )

  projection_rate_train <- mean(
    invalid_train[baseline_train]
  )

  projection_rate_validation <- mean(
    invalid_validation[baseline_validation]
  )

  # Lbar and Ubar are fitted regressions in every mode, including oracle_all.
  # Their finite-sample violations are numerical estimation error rather than
  # a formal model-compatibility test. Stabilize minimally, retain the original
  # violation rates above, and continue the replication.
  Lbar_train <- pmin(
    Lbar_train,
    1 - tau_epsilon
  )

  Ubar_train <- pmax(
    Ubar_train,
    1 + tau_epsilon
  )

  Lbar_validation <- pmin(
    Lbar_validation,
    1 - tau_epsilon
  )

  Ubar_validation <- pmax(
    Ubar_validation,
    1 + tau_epsilon
  )

  tau_train_raw <- (
    Ubar_train - 1
  ) / (
    Ubar_train - Lbar_train
  )

  tau_validation_raw <- (
    Ubar_validation - 1
  ) / (
    Ubar_validation - Lbar_validation
  )

  tau_train <- tau_train_raw
  tau_validation <- tau_validation_raw

  tau_invalid_train <- (
    !is.finite(tau_train) |
      tau_train <= 0 |
      tau_train >= 1
  )

  tau_invalid_validation <- (
    !is.finite(tau_validation) |
      tau_validation <= 0 |
      tau_validation >= 1
  )

  tau_train[!is.finite(tau_train)] <- 0.5
  tau_validation[!is.finite(tau_validation)] <- 0.5

  tau_train <- pmin(
    1 - tau_epsilon,
    pmax(
      tau_epsilon,
      tau_train
    )
  )

  tau_validation <- pmin(
    1 - tau_epsilon,
    pmax(
      tau_epsilon,
      tau_validation
    )
  )

  target_odds_train <- compute_target_odds(
    nuisance = nuisance_train,
    t_value = t_value
  )

  threshold_weight <- (
    1 / nuisance_train$e[
      baseline_train
    ]
  ) *
    target_odds_train[
      baseline_train
    ] *
    envelope_train$D[
      baseline_train
    ]

  if (
    any(threshold_weight < 0)
  ) {
    stop(
      sprintf(
        "Negative threshold weights for treatment %d.",
        t_value
      )
    )
  }

  resolved_threshold_control <- resolve_threshold_control(
    threshold_control
  )
  threshold_rkhs_plan <- if (threshold_learner == "rkhs") {
    make_threshold_rkhs_plan(
      data = train[baseline_train, , drop = FALSE],
      x_cols = x_cols,
      backend = resolved_threshold_control$rkhs_backend,
      bandwidth_multipliers =
        resolved_threshold_control$rkhs_bandwidth_multipliers,
      rank = resolved_threshold_control$rkhs_rank,
      cv_folds = resolved_threshold_control$cv_folds,
      bandwidth_sample =
        resolved_threshold_control$rkhs_bandwidth_sample,
      seed = threshold_seed
    )
  } else if (threshold_learner == "growing_rbf_sieve") {
    make_threshold_growing_rbf_plan(
      data = train[baseline_train, , drop = FALSE],
      x_cols = x_cols,
      bandwidth_grid =
        resolved_threshold_control$growing_rbf_bandwidth_grid,
      rank_multiplier =
        resolved_threshold_control$growing_rbf_rank_multiplier,
      rank_exponent =
        resolved_threshold_control$growing_rbf_rank_exponent,
      minimum_rank = resolved_threshold_control$growing_rbf_min_rank,
      center_skip = resolved_threshold_control$growing_rbf_center_skip,
      bandwidth_rule =
        resolved_threshold_control$growing_rbf_bandwidth_rule,
      bandwidth_sample =
        resolved_threshold_control$growing_rbf_bandwidth_sample,
      cv_folds = resolved_threshold_control$cv_folds,
      seed = threshold_seed
    )
  } else {
    NULL
  }

  make_external_threshold_fit <- function(predict_function, q_train, level) {
    list(
      predict = predict_function,
      convergence = 0L,
      objective = sum(
        threshold_weight * check_loss(
          train$Y[baseline_train], q_train[baseline_train],
          level[baseline_train]
        )
      ),
      penalized_objective = NA_real_,
      counts = c(iterations = 0L),
      selected_spline_df = NA_integer_,
      selected_lambda = NA_real_,
      selected_bandwidth_multiplier = NA_real_,
      selected_bandwidth = NA_real_,
      base_bandwidth = NA_real_,
      selected_rank = NA_integer_,
      cv_table = data.frame(),
      primal_residual = NA_real_,
      dual_residual = NA_real_,
      final_rho = NA_real_,
      rkhs_backend = "external_parametric",
      minimum_reciprocal_condition = NA_real_,
      maximum_numerical_jitter = NA_real_,
      kernel_minimum_raw_eigenvalue = NA_real_,
      tie_break_ridge = NA_real_,
      tie_break_ridge_numerator = NA_real_,
      tie_break_ridge_training_size = sum(baseline_train)
    )
  }

  if (threshold_learner == "external") {
    provider_thresholds <- bound_component_provider$thresholds(
      train = train,
      validation = validation,
      nuisance_train = nuisance_train,
      nuisance_validation = nuisance_validation,
      envelope_train = envelope_train,
      envelope_validation = envelope_validation,
      treatment = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      tau_train = tau_train,
      tau_validation = tau_validation,
      compatibility_action = compatibility_action,
      gamma_projection_margin = gamma_projection_margin,
      envelope_smoothing_epsilon = envelope_smoothing_epsilon,
      tau_epsilon = tau_epsilon
    )
    required_thresholds <- c(
      "q_lower_train", "q_upper_train", "q_lower_validation",
      "q_upper_validation", "lower_predict", "upper_predict"
    )
    if (
      !is.list(provider_thresholds) ||
        !all(required_thresholds %in% names(provider_thresholds)) ||
        !is.function(provider_thresholds$lower_predict) ||
        !is.function(provider_thresholds$upper_predict)
    ) {
      stop(
        paste0(
          "bound_component_provider$thresholds returned an incompatible ",
          "threshold object."
        )
      )
    }
    q_lower_train <- as.numeric(provider_thresholds$q_lower_train)
    q_upper_train <- as.numeric(provider_thresholds$q_upper_train)
    q_lower_validation <- as.numeric(provider_thresholds$q_lower_validation)
    q_upper_validation <- as.numeric(provider_thresholds$q_upper_validation)
    if (
      length(q_lower_train) != nrow(train) ||
        length(q_upper_train) != nrow(train) ||
        length(q_lower_validation) != nrow(validation) ||
        length(q_upper_validation) != nrow(validation) ||
        any(!is.finite(c(
          q_lower_train, q_upper_train,
          q_lower_validation, q_upper_validation
        )))
    ) {
      stop("The external threshold provider returned invalid predictions.")
    }
    lower_threshold_fit <- make_external_threshold_fit(
      provider_thresholds$lower_predict, q_lower_train, 1 - tau_train
    )
    upper_threshold_fit <- make_external_threshold_fit(
      provider_thresholds$upper_predict, q_upper_train, tau_train
    )
  } else {
    lower_threshold_fit <- fit_threshold_learner(
      data = train[
        baseline_train,
        ,
        drop = FALSE
      ],
      y = train$Y[
        baseline_train
      ],
      level = 1 -
        tau_train[
          baseline_train
        ],
      weight = threshold_weight,
      x_cols = x_cols,
      learner = threshold_learner,
      spline_df = threshold_spline_df,
      maxit = threshold_maxit,
      control = resolved_threshold_control,
      rkhs_plan = threshold_rkhs_plan,
      seed = threshold_seed
    )

    upper_threshold_fit <- fit_threshold_learner(
      data = train[
        baseline_train,
        ,
        drop = FALSE
      ],
      y = train$Y[
        baseline_train
      ],
      level =
        tau_train[
          baseline_train
        ],
      weight = threshold_weight,
      x_cols = x_cols,
      learner = threshold_learner,
      spline_df = threshold_spline_df,
      maxit = threshold_maxit,
      control = resolved_threshold_control,
      rkhs_plan = threshold_rkhs_plan,
      seed = threshold_seed
    )
    q_lower_train <- lower_threshold_fit$predict(train)
    q_upper_train <- upper_threshold_fit$predict(train)
    q_lower_validation <- lower_threshold_fit$predict(validation)
    q_upper_validation <- upper_threshold_fit$predict(validation)
  }

  fixed_threshold_target_evaluation <- if (
    is.null(fixed_threshold_target_evaluator)
  ) {
    data.frame()
  } else {
    value <- fixed_threshold_target_evaluator(
      lower_predict = lower_threshold_fit$predict,
      upper_predict = upper_threshold_fit$predict,
      treatment = as.integer(t_value),
      mode = mode
    )
    if (!is.data.frame(value)) {
      stop("fixed_threshold_target_evaluator must return a data frame.")
    }
    reserved_columns <- intersect(
      names(value), c("fold", "treatment", "mode")
    )
    if (length(reserved_columns) > 0L) {
      stop(
        "fixed_threshold_target_evaluator returned reserved columns: ",
        paste(reserved_columns, collapse = ", "), "."
      )
    }
    value
  }

  train_contributions <- evaluate_threshold_bound_contributions(
    data = train,
    nuisance = nuisance_train,
    envelope = envelope_train,
    t_value = t_value,
    tau = tau_train,
    q_lower = q_lower_train,
    q_upper = q_upper_train
  )
  validation_contributions <- evaluate_threshold_bound_contributions(
    data = validation,
    nuisance = nuisance_validation,
    envelope = envelope_validation,
    t_value = t_value,
    tau = tau_validation,
    q_lower = q_lower_validation,
    q_upper = q_upper_validation
  )

  lambda_lower <- validation_contributions$lambda_lower
  lambda_upper <- validation_contributions$lambda_upper
  contribution_lower <- validation_contributions$direct_lower
  contribution_upper <- validation_contributions$direct_upper
  contribution_lower_pseudo <- validation_contributions$pseudo_lower
  contribution_upper_pseudo <- validation_contributions$pseudo_upper
  contribution_lambda_one <- validation_contributions$lambda_one

  prediction_table <- data.frame(
    row_id = row_id_validation[
      baseline_validation
    ],
    treatment = t_value,
    Y = validation$Y[
      baseline_validation
    ],
    X = validation[
      baseline_validation,
      x_cols[1],
      drop = TRUE
    ],
    L = envelope_validation$L[
      baseline_validation
    ],
    U = envelope_validation$U[
      baseline_validation
    ],
    D = envelope_validation$D[
      baseline_validation
    ],
    hard_L = envelope_validation$hard_L[
      baseline_validation
    ],
    hard_U = envelope_validation$hard_U[
      baseline_validation
    ],
    lower_smoothing_gap = envelope_validation$lower_smoothing_gap[
      baseline_validation
    ],
    upper_smoothing_gap = envelope_validation$upper_smoothing_gap[
      baseline_validation
    ],
    Lbar = Lbar_validation[
      baseline_validation
    ],
    Ubar = Ubar_validation[
      baseline_validation
    ],
    tau_raw = tau_validation_raw[
      baseline_validation
    ],
    tau = tau_validation[
      baseline_validation
    ],
    q_lower = q_lower_validation[
      baseline_validation
    ],
    q_upper = q_upper_validation[
      baseline_validation
    ],
    lambda_lower = lambda_lower[
      baseline_validation
    ],
    lambda_upper = lambda_upper[
      baseline_validation
    ],
    direct_lower_contribution = contribution_lower[
      baseline_validation
    ],
    direct_upper_contribution = contribution_upper[
      baseline_validation
    ],
    pseudo_lower_contribution = contribution_lower_pseudo[
      baseline_validation
    ],
    pseudo_upper_contribution = contribution_upper_pseudo[
      baseline_validation
    ],
    lower_q_score_residual =
      validation_contributions$lower_score_residual[
        baseline_validation
      ],
    upper_q_score_residual =
      validation_contributions$upper_score_residual[
        baseline_validation
      ],
    threshold_learner = threshold_learner,
    threshold_rkhs_backend = lower_threshold_fit$rkhs_backend,
    gamma_original =
      envelope_validation$gamma_original[
        baseline_validation
      ],
    gamma_used =
      envelope_validation$gamma_used[
        baseline_validation
      ],
    gamma_compatible_lower =
      envelope_validation$gamma_compatible_lower[
        baseline_validation
      ],
    gamma_compatible_upper =
      envelope_validation$gamma_compatible_upper[
        baseline_validation
      ],
    gamma_projected =
      envelope_validation$gamma_projected[
        baseline_validation
      ]
  )

  eif_fit <- if (compute_eif) {
    eif_arm_train <- train$T_obs == t_value
    eif_arm_validation <- validation$T_obs == t_value
    eif_gamma_projection_rate_train <- mean(
      envelope_train$gamma_projected[eif_arm_train]
    )
    eif_gamma_projection_rate_validation <- mean(
      envelope_validation$gamma_projected[eif_arm_validation]
    )
    eif_projection_active <-
      eif_gamma_projection_rate_train > 0 ||
      eif_gamma_projection_rate_validation > 0
    if (eif_projection_active && !allow_projected_eif) {
      stop(
        paste0(
          "The EIF fit encountered projected gamma predictions. ",
          "Set allow_projected_eif=TRUE only when retaining the ",
          "projection rate as an inference diagnostic."
        )
      )
    }

    result <- fit_eif_arm_fold(
      train = train,
      validation = validation,
      row_id_validation = row_id_validation,
      nuisance_train = nuisance_train,
      nuisance_validation = nuisance_validation,
      envelope_train = envelope_train,
      envelope_validation = envelope_validation,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      q_lower_train = q_lower_train,
      q_upper_train = q_upper_train,
      q_lower_validation = q_lower_validation,
      q_upper_validation = q_upper_validation,
      x_cols = x_cols,
      conditional_spline_df =
        eif_conditional_spline_df,
      conditional_learner = eif_conditional_learner,
      conditional_control = eif_conditional_control,
      conditional_seed = eif_conditional_seed,
      conditional_prediction_provider =
        eif_conditional_prediction_provider,
      instrument_probability_validation =
        eif_instrument_probability_validation
    )
    result$diagnostics$gamma_projection_active <-
      eif_projection_active
    result$diagnostics$gamma_projection_rate_train <-
      eif_gamma_projection_rate_train
    result$diagnostics$gamma_projection_rate_validation <-
      eif_gamma_projection_rate_validation
    result$diagnostics$conditional_projection_rate_train <-
      projection_rate_train
    result$diagnostics$conditional_projection_rate_validation <-
      projection_rate_validation
    result$diagnostics$projected_nuisance_inference_flag <-
      eif_projection_active
    result
  } else {
    NULL
  }

  threshold_cv <- list()
  if (nrow(lower_threshold_fit$cv_table) > 0L) {
    local <- lower_threshold_fit$cv_table
    local$endpoint <- "lower"
    local$learner <- threshold_learner
    threshold_cv[[length(threshold_cv) + 1L]] <- local
  }
  if (nrow(upper_threshold_fit$cv_table) > 0L) {
    local <- upper_threshold_fit$cv_table
    local$endpoint <- "upper"
    local$learner <- threshold_learner
    threshold_cv[[length(threshold_cv) + 1L]] <- local
  }
  threshold_cv <- if (length(threshold_cv) > 0L) {
    do.call(rbind, threshold_cv)
  } else {
    data.frame()
  }

  threshold_iteration_count <- function(fit) {
    if ("iterations" %in% names(fit$counts)) {
      return(unname(fit$counts[["iterations"]]))
    }
    if ("function" %in% names(fit$counts)) {
      return(unname(fit$counts[["function"]]))
    }
    NA_real_
  }
  threshold_fit_value <- function(fit, name, default) {
    value <- fit[[name]]
    if (is.null(value) || length(value) != 1L) default else value
  }

  lower_identity_error <- max(abs(
    contribution_lower - contribution_lower_pseudo -
      validation_contributions$lower_score_residual
  ))
  upper_identity_error <- max(abs(
    contribution_upper - contribution_upper_pseudo -
      validation_contributions$upper_score_residual
  ))

  list(
    lower_contribution = contribution_lower,
    upper_contribution = contribution_upper,
    lower_pseudo_contribution = contribution_lower_pseudo,
    upper_pseudo_contribution = contribution_upper_pseudo,
    lambda_one_contribution =
      contribution_lambda_one,
    predictions = prediction_table,
    diagnostics = data.frame(
      treatment = t_value,
      n_train_baseline =
        sum(baseline_train),
      n_validation_baseline =
        sum(baseline_validation),
      pointwise_crossing_rate_train =
        mean(pointwise_crossing_train[baseline_train]),
      pointwise_crossing_rate_validation =
        mean(pointwise_crossing_validation[baseline_validation]),
      gamma_incompatibility_rate_train =
        gamma_incompatibility_rate_train,
      gamma_incompatibility_rate_validation =
        gamma_incompatibility_rate_validation,
      gamma_projection_rate_train =
        gamma_projection_rate_train,
      gamma_projection_rate_validation =
        gamma_projection_rate_validation,
      gamma_projection_threshold_exceeded =
        gamma_projection_threshold_exceeded,
      envelope_smoothing_epsilon = envelope_smoothing_epsilon,
      mean_lower_smoothing_gap_validation = mean(
        envelope_validation$lower_smoothing_gap[baseline_validation]
      ),
      mean_upper_smoothing_gap_validation = mean(
        envelope_validation$upper_smoothing_gap[baseline_validation]
      ),
      maximum_lower_smoothing_gap_validation = max(
        envelope_validation$lower_smoothing_gap[baseline_validation]
      ),
      maximum_upper_smoothing_gap_validation = max(
        envelope_validation$upper_smoothing_gap[baseline_validation]
      ),
      conditional_projection_rate_train =
        projection_rate_train,
      conditional_projection_rate_validation =
        projection_rate_validation,
      tau_invalid_rate_train =
        mean(tau_invalid_train[baseline_train]),
      tau_invalid_rate_validation =
        mean(tau_invalid_validation[baseline_validation]),
      lower_threshold_convergence =
        lower_threshold_fit$convergence,
      upper_threshold_convergence =
        upper_threshold_fit$convergence,
      lower_threshold_objective =
        lower_threshold_fit$objective,
      upper_threshold_objective =
        upper_threshold_fit$objective,
      mean_L_validation =
        mean(
          envelope_validation$L[
            baseline_validation
          ]
        ),
      mean_U_validation =
        mean(
          envelope_validation$U[
            baseline_validation
          ]
        ),
      mean_tau_validation =
        mean(
          tau_validation[
            baseline_validation
          ]
        ),
      mean_lambda_lower_validation =
        mean(
          lambda_lower[
            baseline_validation
          ]
        ),
      mean_lambda_upper_validation =
        mean(
          lambda_upper[
            baseline_validation
          ]
        ),
      threshold_learner = threshold_learner,
      threshold_rkhs_backend = lower_threshold_fit$rkhs_backend,
      lower_selected_spline_df =
        lower_threshold_fit$selected_spline_df,
      upper_selected_spline_df =
        upper_threshold_fit$selected_spline_df,
      lower_selected_lambda =
        lower_threshold_fit$selected_lambda,
      upper_selected_lambda =
        upper_threshold_fit$selected_lambda,
      lower_selected_bandwidth_multiplier =
        lower_threshold_fit$selected_bandwidth_multiplier,
      upper_selected_bandwidth_multiplier =
        upper_threshold_fit$selected_bandwidth_multiplier,
      lower_selected_bandwidth =
        lower_threshold_fit$selected_bandwidth,
      upper_selected_bandwidth =
        upper_threshold_fit$selected_bandwidth,
      lower_base_bandwidth =
        lower_threshold_fit$base_bandwidth,
      upper_base_bandwidth =
        upper_threshold_fit$base_bandwidth,
      lower_selected_rank =
        lower_threshold_fit$selected_rank,
      upper_selected_rank =
        upper_threshold_fit$selected_rank,
      lower_threshold_iterations =
        threshold_iteration_count(lower_threshold_fit),
      upper_threshold_iterations =
        threshold_iteration_count(upper_threshold_fit),
      lower_threshold_initial_convergence =
        threshold_fit_value(
          lower_threshold_fit, "initial_convergence", NA_integer_
        ),
      upper_threshold_initial_convergence =
        threshold_fit_value(
          upper_threshold_fit, "initial_convergence", NA_integer_
        ),
      lower_threshold_retry_used =
        threshold_fit_value(lower_threshold_fit, "retry_used", NA),
      upper_threshold_retry_used =
        threshold_fit_value(upper_threshold_fit, "retry_used", NA),
      lower_threshold_retry_accepted =
        threshold_fit_value(lower_threshold_fit, "retry_accepted", NA),
      upper_threshold_retry_accepted =
        threshold_fit_value(upper_threshold_fit, "retry_accepted", NA),
      lower_threshold_retry_convergence =
        threshold_fit_value(
          lower_threshold_fit, "retry_convergence", NA_integer_
        ),
      upper_threshold_retry_convergence =
        threshold_fit_value(
          upper_threshold_fit, "retry_convergence", NA_integer_
        ),
      lower_threshold_initial_penalized_objective =
        threshold_fit_value(
          lower_threshold_fit,
          "initial_penalized_objective",
          NA_real_
        ),
      upper_threshold_initial_penalized_objective =
        threshold_fit_value(
          upper_threshold_fit,
          "initial_penalized_objective",
          NA_real_
        ),
      lower_threshold_penalized_objective =
        lower_threshold_fit$penalized_objective,
      upper_threshold_penalized_objective =
        upper_threshold_fit$penalized_objective,
      lower_threshold_tie_break_ridge =
        lower_threshold_fit$tie_break_ridge,
      upper_threshold_tie_break_ridge =
        upper_threshold_fit$tie_break_ridge,
      lower_threshold_tie_break_ridge_numerator =
        lower_threshold_fit$tie_break_ridge_numerator,
      upper_threshold_tie_break_ridge_numerator =
        upper_threshold_fit$tie_break_ridge_numerator,
      lower_threshold_tie_break_ridge_training_size =
        lower_threshold_fit$tie_break_ridge_training_size,
      upper_threshold_tie_break_ridge_training_size =
        upper_threshold_fit$tie_break_ridge_training_size,
      lower_threshold_primal_residual =
        lower_threshold_fit$primal_residual,
      upper_threshold_primal_residual =
        upper_threshold_fit$primal_residual,
      lower_threshold_dual_residual =
        lower_threshold_fit$dual_residual,
      upper_threshold_dual_residual =
        upper_threshold_fit$dual_residual,
      lower_threshold_final_rho =
        lower_threshold_fit$final_rho,
      upper_threshold_final_rho =
        upper_threshold_fit$final_rho,
      lower_threshold_reciprocal_condition =
        lower_threshold_fit$minimum_reciprocal_condition,
      upper_threshold_reciprocal_condition =
        upper_threshold_fit$minimum_reciprocal_condition,
      lower_threshold_numerical_jitter =
        lower_threshold_fit$maximum_numerical_jitter,
      upper_threshold_numerical_jitter =
        upper_threshold_fit$maximum_numerical_jitter,
      lower_kernel_minimum_raw_eigenvalue =
        lower_threshold_fit$kernel_minimum_raw_eigenvalue,
      upper_kernel_minimum_raw_eigenvalue =
        upper_threshold_fit$kernel_minimum_raw_eigenvalue,
      fold_lower_mean =
        mean(contribution_lower),
      fold_upper_mean =
        mean(contribution_upper),
      fold_lambda_one_mean =
        mean(contribution_lambda_one),
      fold_lower_pseudo_mean =
        mean(contribution_lower_pseudo),
      fold_upper_pseudo_mean =
        mean(contribution_upper_pseudo),
      lower_direct_minus_pseudo_train = mean(
        train_contributions$direct_lower -
          train_contributions$pseudo_lower
      ),
      upper_direct_minus_pseudo_train = mean(
        train_contributions$direct_upper -
          train_contributions$pseudo_upper
      ),
      lower_direct_minus_pseudo_validation = mean(
        contribution_lower - contribution_lower_pseudo
      ),
      upper_direct_minus_pseudo_validation = mean(
        contribution_upper - contribution_upper_pseudo
      ),
      lower_q_score_residual_train = mean(
        train_contributions$lower_score_residual
      ),
      upper_q_score_residual_train = mean(
        train_contributions$upper_score_residual
      ),
      lower_q_score_residual_validation = mean(
        validation_contributions$lower_score_residual
      ),
      upper_q_score_residual_validation = mean(
        validation_contributions$upper_score_residual
      ),
      lower_direct_pseudo_identity_error_validation =
        lower_identity_error,
      upper_direct_pseudo_identity_error_validation =
        upper_identity_error
    ),
    threshold_cv = threshold_cv,
    fixed_threshold_target_evaluation =
      fixed_threshold_target_evaluation,
    eif = eif_fit
  )
}

# =============================================================================
# 11. Diagnostics
# =============================================================================

gamma_accuracy_metrics_simple <- function(
    estimated,
    oracle
) {
  keep <- is.finite(estimated) &
    is.finite(oracle)

  estimated <- estimated[keep]
  oracle <- oracle[keep]

  if (length(estimated) < 1L) {
    return(data.frame(
      n = 0L,
      bias = NA_real_,
      rmse = NA_real_,
      mae = NA_real_,
      correlation = NA_real_,
      mean_estimated = NA_real_,
      mean_oracle = NA_real_,
      rmse_constant_one = NA_real_,
      explained_fraction = NA_real_
    ))
  }

  data.frame(
    n = length(estimated),
    bias = mean(estimated - oracle),
    rmse = sqrt(
      mean(
        (estimated - oracle)^2
      )
    ),
    mae = mean(
      abs(estimated - oracle)
    ),
    correlation = if (
      length(estimated) > 1L &&
        is.finite(stats::sd(estimated)) &&
        is.finite(stats::sd(oracle)) &&
        stats::sd(estimated) > 0 &&
        stats::sd(oracle) > 0
    ) {
      stats::cor(
        estimated,
        oracle
      )
    } else {
      NA_real_
    },
    mean_estimated = mean(estimated),
    mean_oracle = mean(oracle),
    rmse_constant_one = sqrt(
      mean(
        (1 - oracle)^2
      )
    ),
    explained_fraction = 1 -
      mean(
        (estimated - oracle)^2
      ) /
      mean(
        (1 - oracle)^2
      )
  )
}

compute_nuisance_diagnostics_simple <- function(
    data
) {
  accuracy <- function(estimated, truth) {
    correlation <- if (
      stats::sd(estimated) > 0 && stats::sd(truth) > 0
    ) {
      stats::cor(estimated, truth)
    } else {
      NA_real_
    }
    data.frame(
      bias = mean(estimated - truth),
      rmse = sqrt(
        mean(
          (estimated - truth)^2
        )
      ),
      mae = mean(
        abs(estimated - truth)
      ),
      correlation = correlation,
      mean_estimated = mean(estimated),
      mean_true = mean(truth)
    )
  }

  estimated_nuisance <- list(
    e = data$e_hat,
    pi0 = data$pi0_hat,
    pi1 = data$pi1_hat
  )

  if (!all(c("e_true", "pi0_true", "pi1_true") %in% names(data))) {
    unavailable <- data.frame(
      bias = NA_real_, rmse = NA_real_, mae = NA_real_,
      correlation = NA_real_, mean_estimated = NA_real_,
      mean_true = NA_real_
    )
    return(list(
      e = unavailable, pi0 = unavailable, pi1 = unavailable,
      w1 = unavailable, w0 = unavailable
    ))
  }

  oracle_nuisance <- list(
    e = data$e_true,
    pi0 = data$pi0_true,
    pi1 = data$pi1_true
  )

  w1_hat <- compute_riesz_weight(
    estimated_nuisance,
    t_value = 1
  )

  w1_true <- compute_riesz_weight(
    oracle_nuisance,
    t_value = 1
  )

  w0_hat <- compute_riesz_weight(
    estimated_nuisance,
    t_value = 0
  )

  w0_true <- compute_riesz_weight(
    oracle_nuisance,
    t_value = 0
  )

  list(
    e = accuracy(
      data$e_hat,
      data$e_true
    ),
    pi0 = accuracy(
      data$pi0_hat,
      data$pi0_true
    ),
    pi1 = accuracy(
      data$pi1_hat,
      data$pi1_true
    ),
    w1 = accuracy(
      w1_hat,
      w1_true
    ),
    w0 = accuracy(
      w0_hat,
      w0_true
    )
  )
}

compute_ipw_identities_simple <- function(
    data
) {
  c(
    t1_z0_estimated = mean(
      (data$Z == 0) *
        data$T_obs /
        (
          data$e_hat *
            data$pi0_hat
        )
    ),
    t1_z0_oracle = mean(
      (data$Z == 0) *
        data$T_obs /
        (
          data$e_true *
            data$pi0_true
        )
    ),
    t0_z0_estimated = mean(
      (data$Z == 0) *
        (1 - data$T_obs) /
        (
          data$e_hat *
            (1 - data$pi0_hat)
        )
    ),
    t0_z0_oracle = mean(
      (data$Z == 0) *
        (1 - data$T_obs) /
        (
          data$e_true *
            (1 - data$pi0_true)
        )
    ),
    t1_z1_estimated = mean(
      data$Z *
        data$T_obs /
        (
          (1 - data$e_hat) *
            data$pi1_hat
        )
    ),
    t1_z1_oracle = mean(
      data$Z *
        data$T_obs /
        (
          (1 - data$e_true) *
            data$pi1_true
        )
    ),
    t0_z1_estimated = mean(
      data$Z *
        (1 - data$T_obs) /
        (
          (1 - data$e_hat) *
            (1 - data$pi1_hat)
        )
    ),
    t0_z1_oracle = mean(
      data$Z *
        (1 - data$T_obs) /
        (
          (1 - data$e_true) *
            (1 - data$pi1_true)
        )
    )
  )
}

build_normalization_by_x_bin <- function(
    predictions,
    x_bins = 10
) {
  output <- list()
  counter <- 1

  groups <- split(
    predictions,
    interaction(
      predictions$mode,
      predictions$treatment,
      drop = TRUE
    )
  )

  for (group in groups) {
    ranks <- rank(
      group$X,
      ties.method = "average"
    )

    bin <- pmin(
      x_bins,
      ceiling(
        x_bins *
          ranks /
          length(ranks)
      )
    )

    for (b in seq_len(x_bins)) {
      local <- group[
        bin == b,
        ,
        drop = FALSE
      ]

      if (nrow(local) == 0) {
        next
      }

      output[[counter]] <- data.frame(
        mode = local$mode[1],
        treatment = local$treatment[1],
        x_bin = b,
        n = nrow(local),
        x_min = min(local$X),
        x_max = max(local$X),
        x_mean = mean(local$X),
        mean_Lbar = mean(local$Lbar),
        mean_Ubar = mean(local$Ubar),
        mean_tau = mean(local$tau),
        mean_lambda_lower =
          mean(local$lambda_lower),
        mean_lambda_upper =
          mean(local$lambda_upper),
        mean_gamma =
          mean(local$gamma_used)
      )

      counter <- counter + 1
    }
  }

  do.call(
    rbind,
    output
  )
}

# =============================================================================
# 12. Complete cross-fitted plug-in sharp-bound estimator
# =============================================================================

estimate_sharp_bounds_plugin <- function(
    data,
    Lambda0_T1 = 5,
    Lambda1_T1 = 5,
    Lambda0_T0 = 5,
    Lambda1_T0 = 5,
    x_cols = "X",
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    K = 5,
    fold_id = NULL,
    positivity_epsilon = 0.05,
    nuisance_learner = c(
      "gam", "ranger", "glm_spline", "library_cv", "super_learner",
      "integrated_parametric"
    ),
    treatment_spline_df = 5,
    nuisance_gam_k = NULL,
    nuisance_gam_min_k = 5L,
    nuisance_gam_max_k = 20L,
    nuisance_gam_growth_exponent = 1 / 3,
    nuisance_gam_method = "REML",
    nuisance_library_candidates = c("glm_spline", "gam", "ranger"),
    nuisance_library_cv_folds = 3L,
    nuisance_library_selection_rule = c("one_se", "minimum"),
    nuisance_library_ranger_num_trees = 500L,
    nuisance_super_learner_candidates = NULL,
    nuisance_super_learner_cv_folds = 3L,
    nuisance_super_learner_ranger_num_trees = 500L,
    nuisance_super_learner_hal_max_degree = 2L,
    nuisance_integrated_nodes = 40L,
    nuisance_integrated_maxit = 500L,
    gamma_floor = 1e-6,
    gamma_cap = 100,
    oracle_nodes = 80,
    oracle_batch_size = 1000,
    rkhs_lambda_by_t = NULL,
    rkhs_lambda_grid = default_gamma_lambda_grid,
    rkhs_bandwidth_multipliers = 1 / sqrt(2),
    rkhs_cv_folds = 3,
    rkhs_cv_repeats = 2L,
    rkhs_selection_rule =
      c("minimum", "one_se"),
    gamma_rkhs_backend = c("exact", "nystrom"),
    gamma_rkhs_kernel = c(
      "gaussian", "sobolev_matern32", "sobolev_matern52"
    ),
    gamma_nystrom_rank = NULL,
    gamma_nystrom_eigen_tolerance = 1e-10,
    gamma_nystrom_jitter = 1e-10,
    bandwidth_sample = 1000,
    conditional_mean_spline_df = 5,
    threshold_spline_df = 5,
    threshold_learner = c(
      "spline_fixed", "spline_cv", "rkhs", "fixed_rbf",
      "growing_rbf_sieve", "external"
    ),
    threshold_control = NULL,
    envelope_smoothing_epsilon = 0,
    estimated_compatibility_action =
      c("project", "stop"),
    compatibility_tolerance = 1e-6,
    gamma_projection_margin = 1e-4,
    max_gamma_projection_rate = 0.02,
    tau_epsilon = 1e-4,
    threshold_maxit = 4000,
    compute_oracle_eif = TRUE,
    compute_estimated_eif = FALSE,
    one_step_modes = NULL,
    estimation_modes = c("oracle_all", "oracle_gamma", "rkhs"),
    eif_conditional_spline_df =
      conditional_mean_spline_df,
    eif_conditional_learner = c(
      "spline_fixed", "spline_cv", "gam", "ranger", "library_cv"
    ),
    eif_conditional_control = NULL,
    eif_conditional_provider_factory = NULL,
    eif_instrument_probability_transform = NULL,
    nuisance_fit_provider_factory = NULL,
    gamma_fit_provider_factory = NULL,
    gamma_fit_transform_factory = NULL,
    bound_component_provider_factory = NULL,
    fixed_threshold_target_evaluator = NULL,
    x_bins = 10,
    seed = 1
) {
  threshold_learner <- match.arg(threshold_learner)
  nuisance_learner <- match.arg(nuisance_learner)
  nuisance_library_selection_rule <- match.arg(
    nuisance_library_selection_rule
  )
  eif_conditional_learner <- match.arg(eif_conditional_learner)
  if (
    !is.null(eif_conditional_provider_factory) &&
      !is.function(eif_conditional_provider_factory)
  ) {
    stop("eif_conditional_provider_factory must be NULL or a function.")
  }
  if (
    !is.null(eif_instrument_probability_transform) &&
      !is.function(eif_instrument_probability_transform)
  ) {
    stop(
      "eif_instrument_probability_transform must be NULL or a function."
    )
  }
  if (
    !is.null(nuisance_fit_provider_factory) &&
      !is.function(nuisance_fit_provider_factory)
  ) {
    stop("nuisance_fit_provider_factory must be NULL or a function.")
  }
  if (
    !is.null(gamma_fit_provider_factory) &&
      !is.function(gamma_fit_provider_factory)
  ) {
    stop("gamma_fit_provider_factory must be NULL or a function.")
  }
  if (
    !is.null(gamma_fit_transform_factory) &&
      !is.function(gamma_fit_transform_factory)
  ) {
    stop("gamma_fit_transform_factory must be NULL or a function.")
  }
  if (
    !is.null(bound_component_provider_factory) &&
      !is.function(bound_component_provider_factory)
  ) {
    stop("bound_component_provider_factory must be NULL or a function.")
  }
  if (
    !is.null(fixed_threshold_target_evaluator) &&
      !is.function(fixed_threshold_target_evaluator)
  ) {
    stop("fixed_threshold_target_evaluator must be NULL or a function.")
  }
  rkhs_selection_rule <- match.arg(
    rkhs_selection_rule
  )
  gamma_rkhs_backend <- match.arg(
    gamma_rkhs_backend
  )
  gamma_rkhs_kernel <- match.arg(gamma_rkhs_kernel)

  estimated_compatibility_action <- match.arg(
    estimated_compatibility_action
  )

  if (
    length(envelope_smoothing_epsilon) != 1L ||
      !is.finite(envelope_smoothing_epsilon) ||
      envelope_smoothing_epsilon < 0
  ) {
    stop("envelope_smoothing_epsilon must be one nonnegative finite number.")
  }

  required <- unique(
    c(
      "Y",
      "T_obs",
      "Z",
      x_cols
    )
  )

  missing_columns <- setdiff(
    required,
    names(data)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  analysis_data <- data[
    stats::complete.cases(
      data[, required, drop = FALSE]
    ),
    ,
    drop = FALSE
  ]

  analysis_data <- as.data.frame(
    analysis_data
  )

  analysis_data$.row_id <-
    seq_len(nrow(analysis_data))

  oracle_columns <- c(
    "e_true",
    "pi0_true",
    "pi1_true",
    "gamma0_oracle",
    "gamma1_oracle"
  )

  if (
    any(c("oracle_all", "oracle_gamma") %in% estimation_modes) &&
      !all(
      oracle_columns %in%
        names(analysis_data)
    )
  ) {
    oracle <- compute_oracle_quantities(
      data = analysis_data,
      confounder_strength =
        confounder_strength,
      treatment_intercept =
        treatment_intercept,
      instrument_strength =
        instrument_strength,
      treatment_effect =
        treatment_effect,
      n_nodes = oracle_nodes,
      batch_size =
        oracle_batch_size
    )

    for (name in names(oracle)) {
      analysis_data[[name]] <-
        oracle[[name]]
    }
  }

  n <- nrow(analysis_data)

  if (is.null(fold_id)) {
    fold_id <- make_tz_stratified_folds(
      data = analysis_data,
      K = K,
      seed = seed
    )
  } else {
    if (length(fold_id) != n) {
      stop(
        "fold_id must have one value per analysis row."
      )
    }

    K <- max(fold_id)
  }

  analysis_data$fold_id <- fold_id

  available_modes <- c(
    "oracle_all",
    "oracle_gamma",
    "rkhs"
  )
  modes <- match.arg(
    estimation_modes,
    choices = available_modes,
    several.ok = TRUE
  )
  modes <- unique(as.character(modes))

  eif_modes <- if (is.null(one_step_modes)) {
    intersect(c(
      if (isTRUE(compute_oracle_eif)) "oracle_all",
      if (isTRUE(compute_estimated_eif)) c("oracle_gamma", "rkhs")
    ), modes)
  } else {
    invalid_eif_modes <- setdiff(one_step_modes, available_modes)
    if (length(invalid_eif_modes) > 0L) {
      stop(
        "Unknown one_step_modes: ",
        paste(invalid_eif_modes, collapse = ", ")
      )
    }
    missing_estimation_modes <- setdiff(one_step_modes, modes)
    if (length(missing_estimation_modes) > 0L) {
      stop(
        "Every one_step_modes entry must also be included in estimation_modes: ",
        paste(missing_estimation_modes, collapse = ", ")
      )
    }
    unique(as.character(one_step_modes))
  }

  contribution_storage <- list()

  for (mode in modes) {
    contribution_storage[[mode]] <- list()

    for (t_value in c(0, 1)) {
      contribution_storage[[mode]][[
        as.character(t_value)
      ]] <- list(
        lower = rep(
          NA_real_,
          n
        ),
        upper = rep(
          NA_real_,
          n
        ),
        lower_pseudo = rep(
          NA_real_,
          n
        ),
        upper_pseudo = rep(
          NA_real_,
          n
        ),
        lambda_one = rep(
          NA_real_,
          n
        )
      )
    }
  }

  e_hat <- pi0_hat <- pi1_hat <-
    rep(NA_real_, n)

  gamma0_rkhs <- gamma1_rkhs <-
    rep(NA_real_, n)

  gamma_accuracy_rows <- list()
  gamma_accuracy_counter <- 1

  fold_diagnostics <- list()
  diagnostic_counter <- 1

  prediction_tables <- list()
  prediction_counter <- 1

  selected_lambdas <- list()
  selected_lambda_counter <- 1

  gamma_fit_diagnostics <- list()
  gamma_fit_diagnostic_counter <- 1L

  gamma_transform_diagnostics <- list()
  gamma_transform_diagnostic_counter <- 1L

  cv_tables <- list()
  cv_counter <- 1

  threshold_cv_tables <- list()
  threshold_cv_counter <- 1L

  fixed_threshold_target_evaluations <- list()
  fixed_threshold_target_evaluation_counter <- 1L

  nuisance_fit_diagnostics <- list()
  nuisance_fit_diagnostic_counter <- 1L

  eif_storage <- setNames(
    lapply(eif_modes, function(mode) {
      list(
        "0" = list(
          lower = rep(NA_real_, n),
          upper = rep(NA_real_, n)
        ),
        "1" = list(
          lower = rep(NA_real_, n),
          upper = rep(NA_real_, n)
        )
      )
    }),
    eif_modes
  )

  eif_predictions <- list()
  eif_prediction_counter <- 1L

  eif_diagnostics <- list()
  eif_diagnostic_counter <- 1L

  eif_conditional_cv_tables <- list()
  eif_conditional_cv_counter <- 1L

  for (fold in seq_len(K)) {
    validation_index <- which(
      fold_id == fold
    )

    training_index <- which(
      fold_id != fold
    )

    train <- analysis_data[
      training_index,
      ,
      drop = FALSE
    ]

    validation <- analysis_data[
      validation_index,
      ,
      drop = FALSE
    ]

    nuisance_fit <- if (is.null(nuisance_fit_provider_factory)) {
      fit_nuisance_models(
        train = train,
        x_cols = x_cols,
        positivity_epsilon =
          positivity_epsilon,
        learner = nuisance_learner,
        treatment_spline_df =
          treatment_spline_df,
        gam_basis_dimension = nuisance_gam_k,
        gam_minimum_basis_dimension = nuisance_gam_min_k,
        gam_maximum_basis_dimension = nuisance_gam_max_k,
        gam_basis_growth_exponent = nuisance_gam_growth_exponent,
        gam_method = nuisance_gam_method,
        library_candidates = nuisance_library_candidates,
        library_cv_folds = nuisance_library_cv_folds,
        library_selection_rule = nuisance_library_selection_rule,
        library_ranger_num_trees = nuisance_library_ranger_num_trees,
        library_seed = seed + 90000L * fold,
        super_learner_candidates = nuisance_super_learner_candidates,
        super_learner_cv_folds = nuisance_super_learner_cv_folds,
        super_learner_ranger_num_trees =
          nuisance_super_learner_ranger_num_trees,
        super_learner_hal_max_degree =
          nuisance_super_learner_hal_max_degree,
        integrated_confounder_strength =
          confounder_strength,
        integrated_quadrature_nodes =
          nuisance_integrated_nodes,
        integrated_maxit =
          nuisance_integrated_maxit
      )
    } else {
      nuisance_provider_arguments <- list(
        fold = fold,
        train = train,
        x_cols = x_cols,
        seed = seed + 90000L * fold
      )
      nuisance_provider_formals <- names(formals(
        nuisance_fit_provider_factory
      ))
      if ("..." %in% nuisance_provider_formals) {
        nuisance_provider_arguments$positivity_epsilon <- positivity_epsilon
        nuisance_provider_arguments$learner <- nuisance_learner
      }
      do.call(
        nuisance_fit_provider_factory,
        nuisance_provider_arguments
      )
    }

    if (!is.list(nuisance_fit) ||
        !is.function(nuisance_fit$predict) ||
        !is.data.frame(nuisance_fit$diagnostics)) {
      stop(
        paste0(
          "The nuisance fit provider must return a list containing a ",
          "predict function and a diagnostics data frame."
        )
      )
    }

    local_nuisance_diagnostics <- nuisance_fit$diagnostics
    local_nuisance_diagnostics$fold <- fold
    if (!"n_training" %in% names(local_nuisance_diagnostics)) {
      local_nuisance_diagnostics$n_training <- c(
        nrow(train),
        sum(train$Z == 0),
        sum(train$Z == 1)
      )
    }
    nuisance_fit_diagnostics[[
      nuisance_fit_diagnostic_counter
    ]] <- local_nuisance_diagnostics
    nuisance_fit_diagnostic_counter <-
      nuisance_fit_diagnostic_counter + 1L

    nuisance_train_estimated <-
      nuisance_fit$predict(
        train
      )

    nuisance_validation_estimated <-
      nuisance_fit$predict(
        validation
      )

    nuisance_train_oracle <- list(
      e = train$e_true,
      pi0 = train$pi0_true,
      pi1 = train$pi1_true
    )

    nuisance_validation_oracle <- list(
      e = validation$e_true,
      pi0 = validation$pi0_true,
      pi1 = validation$pi1_true
    )

    e_hat[validation_index] <-
      nuisance_validation_estimated$e

    pi0_hat[validation_index] <-
      nuisance_validation_estimated$pi0

    pi1_hat[validation_index] <-
      nuisance_validation_estimated$pi1

    for (t_value in c(0, 1)) {
      lambda_key <- as.character(
        t_value
      )

      if (!is.null(gamma_fit_provider_factory)) {
        gamma_provider_arguments <- list(
          fold = fold,
          treatment = t_value,
          train = train,
          nuisance_train = nuisance_train_estimated,
          x_cols = x_cols,
          seed = seed + 100000 * fold + 100 * t_value
        )
        gamma_provider_formals <- names(formals(gamma_fit_provider_factory))
        if (
          "nuisance_predict" %in% gamma_provider_formals ||
            "..." %in% gamma_provider_formals
        ) {
          gamma_provider_arguments$nuisance_predict <- nuisance_fit$predict
        }
        rkhs_fit <- do.call(
          gamma_fit_provider_factory,
          gamma_provider_arguments
        )
        required_gamma_fit_fields <- c(
          "predict", "selected_lambda", "selected_bandwidth_multiplier",
          "base_bandwidth", "selected_bandwidth", "normalization_constant",
          "training_negative_rate", "rkhs_backend",
          "nystrom_rank_requested", "nystrom_landmark_count",
          "nystrom_effective_rank", "nystrom_min_retained_eigenvalue",
          "nystrom_max_eigenvalue"
        )
        if (
          !is.list(rkhs_fit) ||
            !all(required_gamma_fit_fields %in% names(rkhs_fit)) ||
            !is.function(rkhs_fit$predict)
        ) {
          stop(
            paste0(
              "gamma_fit_provider_factory returned an incompatible gamma ",
              "fit object."
            )
          )
        }
        if (!is.null(rkhs_fit$cv_table) && nrow(rkhs_fit$cv_table) > 0L) {
          cv_table <- as.data.frame(
            rkhs_fit$cv_table,
            stringsAsFactors = FALSE
          )
          cv_table$fold <- fold
          cv_table$treatment <- t_value
          cv_tables[[cv_counter]] <- cv_table
          cv_counter <- cv_counter + 1L
        }
      } else if (is.null(rkhs_lambda_by_t)) {
        rkhs_fit <- fit_gamma_rkhs(
          train = train,
          nuisance_train =
            nuisance_train_estimated,
          t_value = t_value,
          x_cols = x_cols,
          lambda_grid =
            rkhs_lambda_grid,
          bandwidth_multipliers =
            rkhs_bandwidth_multipliers,
          cv_folds =
            rkhs_cv_folds,
          cv_repeats =
            rkhs_cv_repeats,
          selection_rule =
            rkhs_selection_rule,
          gamma_floor =
            gamma_floor,
          gamma_cap =
            gamma_cap,
          bandwidth_sample =
            bandwidth_sample,
          backend = gamma_rkhs_backend,
          nystrom_rank = gamma_nystrom_rank,
          nystrom_eigen_tolerance =
            gamma_nystrom_eigen_tolerance,
          nystrom_jitter = gamma_nystrom_jitter,
          kernel_type = gamma_rkhs_kernel,
          seed = seed +
            100000 * fold +
            100 * t_value
        )

        cv_table <- rkhs_fit$cv_table
        cv_table$fold <- fold
        cv_table$treatment <- t_value

        cv_tables[[cv_counter]] <-
          cv_table

        cv_counter <- cv_counter + 1
      } else {
        if (
          !lambda_key %in%
            names(rkhs_lambda_by_t)
        ) {
          stop(
            "rkhs_lambda_by_t must be named with entries '0' and '1'."
          )
        }

        rkhs_fit <-
          fit_gamma_rkhs_fixed(
            train = train,
            nuisance_train =
              nuisance_train_estimated,
            t_value = t_value,
            x_cols = x_cols,
            lambda =
              rkhs_lambda_by_t[
                lambda_key
              ],
            gamma_floor =
              gamma_floor,
            gamma_cap =
              gamma_cap,
            bandwidth_sample =
              bandwidth_sample,
            backend = gamma_rkhs_backend,
            nystrom_rank = gamma_nystrom_rank,
            nystrom_eigen_tolerance =
              gamma_nystrom_eigen_tolerance,
            nystrom_jitter = gamma_nystrom_jitter,
            kernel_type = gamma_rkhs_kernel,
            seed = seed +
              100000 * fold +
              100 * t_value
          )
      }

      gamma_fit_pipeline <- rkhs_fit
      gamma_transform <- NULL
      if (!is.null(gamma_fit_transform_factory)) {
        gamma_transform <- gamma_fit_transform_factory(
          fold = fold,
          treatment = t_value,
          gamma_fit = rkhs_fit,
          x_cols = x_cols
        )
        if (
          !is.list(gamma_transform) ||
            !is.function(gamma_transform$transform)
        ) {
          stop(
            paste0(
              "gamma_fit_transform_factory must return a list containing ",
              "a transform function."
            )
          )
        }
        original_predict <- rkhs_fit$predict
        transform_function <- gamma_transform$transform
        gamma_fit_pipeline$predict <- local({
          local_predict <- original_predict
          local_transform <- transform_function
          function(newdata) {
            prediction <- local_predict(newdata)
            transformed <- as.numeric(local_transform(
              newdata = newdata,
              normalized_predictions = prediction$normalized
            ))
            if (
              length(transformed) != nrow(newdata) ||
                any(!is.finite(transformed)) ||
                any(transformed <= 0)
            ) {
              stop(
                "The fitted-gamma transformation returned invalid predictions."
              )
            }
            prediction$normalized <- transformed
            prediction
          }
        })
      }

      gamma_train_rkhs <-
        gamma_fit_pipeline$predict(
          train
        )$normalized

      gamma_validation_rkhs <-
        gamma_fit_pipeline$predict(
          validation
        )$normalized

      if (t_value == 1) {
        gamma1_rkhs[
          validation_index
        ] <- gamma_validation_rkhs

        gamma_oracle_train <-
          train$gamma1_oracle

        gamma_oracle_validation <-
          validation$gamma1_oracle
      } else {
        gamma0_rkhs[
          validation_index
        ] <- gamma_validation_rkhs

        gamma_oracle_train <-
          train$gamma0_oracle

        gamma_oracle_validation <-
          validation$gamma0_oracle
      }

      baseline_validation <- (
        validation$T_obs == t_value &
          validation$Z == 0
      )

      accuracy <- gamma_accuracy_metrics_simple(
        estimated =
          gamma_validation_rkhs[
            baseline_validation
          ],
        oracle =
          gamma_oracle_validation[
            baseline_validation
          ]
      )

      accuracy$fold <- fold
      accuracy$treatment <- t_value

      gamma_accuracy_rows[[
        gamma_accuracy_counter
      ]] <- accuracy

      gamma_accuracy_counter <-
        gamma_accuracy_counter + 1

      selected_lambdas[[
        selected_lambda_counter
      ]] <- data.frame(
        fold = fold,
        treatment = t_value,
        gamma_estimator = if (
          !is.null(rkhs_fit$gamma_estimator)
        ) {
          rkhs_fit$gamma_estimator
        } else {
          "direct_rkhs_riesz"
        },
        lambda =
          rkhs_fit$selected_lambda,
        bandwidth_multiplier =
          rkhs_fit$selected_bandwidth_multiplier,
        base_bandwidth =
          rkhs_fit$base_bandwidth,
        selected_bandwidth =
          rkhs_fit$selected_bandwidth,
        normalization_constant =
          rkhs_fit$normalization_constant,
        training_negative_rate =
          rkhs_fit$training_negative_rate,
        rkhs_backend = rkhs_fit$rkhs_backend,
        rkhs_kernel = if (!is.null(rkhs_fit$rkhs_kernel)) {
          rkhs_fit$rkhs_kernel
        } else {
          gamma_rkhs_kernel
        },
        nystrom_rank_requested =
          rkhs_fit$nystrom_rank_requested,
        nystrom_landmark_count =
          rkhs_fit$nystrom_landmark_count,
        nystrom_effective_rank =
          rkhs_fit$nystrom_effective_rank,
        nystrom_min_retained_eigenvalue =
          rkhs_fit$nystrom_min_retained_eigenvalue,
        nystrom_max_eigenvalue =
          rkhs_fit$nystrom_max_eigenvalue
      )

      selected_lambda_counter <-
        selected_lambda_counter + 1

      if (!is.null(rkhs_fit$fit_diagnostics)) {
        local_gamma_fit_diagnostics <- as.data.frame(
          rkhs_fit$fit_diagnostics,
          stringsAsFactors = FALSE
        )
        local_gamma_fit_diagnostics$fold <- fold
        local_gamma_fit_diagnostics$treatment <- t_value
        gamma_fit_diagnostics[[gamma_fit_diagnostic_counter]] <-
          local_gamma_fit_diagnostics
        gamma_fit_diagnostic_counter <- gamma_fit_diagnostic_counter + 1L
      }

      mode_inputs <- list(
        oracle_all = list(
          nuisance_train =
            nuisance_train_oracle,
          nuisance_validation =
            nuisance_validation_oracle,
          gamma_train =
            gamma_oracle_train,
          gamma_validation =
            gamma_oracle_validation
        ),
        oracle_gamma = list(
          nuisance_train =
            nuisance_train_estimated,
          nuisance_validation =
            nuisance_validation_estimated,
          gamma_train =
            gamma_oracle_train,
          gamma_validation =
            gamma_oracle_validation
        ),
        rkhs = list(
          nuisance_train =
            nuisance_train_estimated,
          nuisance_validation =
            nuisance_validation_estimated,
          gamma_train =
            gamma_train_rkhs,
          gamma_validation =
            gamma_validation_rkhs
        )
      )

      Lambda0 <- if (
        t_value == 1
      ) {
        Lambda0_T1
      } else {
        Lambda0_T0
      }

      Lambda1 <- if (
        t_value == 1
      ) {
        Lambda1_T1
      } else {
        Lambda1_T0
      }

      for (mode in modes) {
        bound_component_provider <- NULL
        if (!is.null(bound_component_provider_factory)) {
          bound_component_provider <- bound_component_provider_factory(
            mode = mode,
            fold = fold,
            treatment = t_value,
            gamma_fit = gamma_fit_pipeline,
            nuisance_fit = nuisance_fit,
            x_cols = x_cols
          )
        }
        eif_conditional_prediction_provider <- NULL
        if (
          mode %in% eif_modes &&
            !is.null(eif_conditional_provider_factory)
        ) {
          eif_conditional_prediction_provider <-
            eif_conditional_provider_factory(
              mode = mode,
              fold = fold,
              treatment = t_value,
              gamma_fit = gamma_fit_pipeline,
              x_cols = x_cols
            )
          if (
            !is.null(eif_conditional_prediction_provider) &&
              !is.function(eif_conditional_prediction_provider)
          ) {
            stop(
              paste0(
                "eif_conditional_provider_factory must return NULL or a ",
                "prediction-provider function."
              )
            )
          }
        }
        eif_instrument_probability_validation <- NULL
        if (mode %in% eif_modes) {
          eif_instrument_probability_validation <-
            mode_inputs[[mode]]$nuisance_validation$e
          if (!is.null(eif_instrument_probability_transform)) {
            eif_instrument_probability_validation <- as.numeric(
              eif_instrument_probability_transform(
                probability = eif_instrument_probability_validation,
                data = validation,
                mode = mode,
                fold = fold,
                treatment = t_value
              )
            )
          }
        }
        fit <- fit_bound_arm_fold(
          train = train,
          validation = validation,
          row_id_validation =
            validation$.row_id,
          nuisance_train =
            mode_inputs[[mode]]$
              nuisance_train,
          nuisance_validation =
            mode_inputs[[mode]]$
              nuisance_validation,
          gamma_train =
            mode_inputs[[mode]]$
              gamma_train,
          gamma_validation =
            mode_inputs[[mode]]$
              gamma_validation,
          t_value = t_value,
          Lambda0 = Lambda0,
          Lambda1 = Lambda1,
          mode = mode,
          x_cols = x_cols,
          conditional_mean_spline_df =
            conditional_mean_spline_df,
          threshold_spline_df =
            threshold_spline_df,
          threshold_learner =
            threshold_learner,
          threshold_control =
            threshold_control,
          threshold_seed = seed +
            500000L * fold +
            1000L * t_value +
            match(mode, modes),
          envelope_smoothing_epsilon = envelope_smoothing_epsilon,
          compatibility_tolerance =
            compatibility_tolerance,
          compatibility_action =
            if (mode == "oracle_all") {
              "stop"
            } else {
              estimated_compatibility_action
            },
          gamma_projection_margin =
            gamma_projection_margin,
          max_gamma_projection_rate =
            max_gamma_projection_rate,
          tau_epsilon =
            tau_epsilon,
          threshold_maxit =
            threshold_maxit,
          compute_eif =
            mode %in% eif_modes &&
            Lambda0 > 1 &&
            Lambda1 > 1,
          allow_projected_eif =
            mode != "oracle_all",
          eif_conditional_spline_df =
            eif_conditional_spline_df,
          eif_conditional_learner = eif_conditional_learner,
          eif_conditional_control = eif_conditional_control,
          eif_conditional_seed = seed +
            700000L * fold +
            10000L * t_value +
            100L * match(mode, modes),
          eif_conditional_prediction_provider =
            eif_conditional_prediction_provider,
          eif_instrument_probability_validation =
            eif_instrument_probability_validation,
          bound_component_provider = bound_component_provider,
          fixed_threshold_target_evaluator =
            fixed_threshold_target_evaluator
        )

        contribution_storage[[mode]][[
          lambda_key
        ]]$lower[
          validation_index
        ] <- fit$lower_contribution

        contribution_storage[[mode]][[
          lambda_key
        ]]$upper[
          validation_index
        ] <- fit$upper_contribution

        contribution_storage[[mode]][[
          lambda_key
        ]]$lambda_one[
          validation_index
        ] <- fit$lambda_one_contribution

        contribution_storage[[mode]][[
          lambda_key
        ]]$lower_pseudo[
          validation_index
        ] <- fit$lower_pseudo_contribution

        contribution_storage[[mode]][[
          lambda_key
        ]]$upper_pseudo[
          validation_index
        ] <- fit$upper_pseudo_contribution

        diagnostic <- fit$diagnostics
        diagnostic$mode <- mode
        diagnostic$fold <- fold

        fold_diagnostics[[
          diagnostic_counter
        ]] <- diagnostic

        diagnostic_counter <-
          diagnostic_counter + 1

        prediction <- fit$predictions
        prediction$mode <- mode
        prediction$fold <- fold

        prediction_tables[[
          prediction_counter
        ]] <- prediction

        prediction_counter <-
          prediction_counter + 1

        if (nrow(fit$threshold_cv) > 0L) {
          local_threshold_cv <- fit$threshold_cv
          local_threshold_cv$mode <- mode
          local_threshold_cv$fold <- fold
          local_threshold_cv$treatment <- t_value
          threshold_cv_tables[[threshold_cv_counter]] <-
            local_threshold_cv
          threshold_cv_counter <- threshold_cv_counter + 1L
        }

        if (nrow(fit$fixed_threshold_target_evaluation) > 0L) {
          local_fixed_target <- fit$fixed_threshold_target_evaluation
          local_fixed_target$mode <- mode
          local_fixed_target$fold <- fold
          local_fixed_target$treatment <- t_value
          fixed_threshold_target_evaluations[[
            fixed_threshold_target_evaluation_counter
          ]] <- local_fixed_target
          fixed_threshold_target_evaluation_counter <-
            fixed_threshold_target_evaluation_counter + 1L
        }

        if (!is.null(fit$eif)) {
          eif_storage[[mode]][[
            lambda_key
          ]]$lower[
            validation_index
          ] <- fit$eif$signals$lower_signal

          eif_storage[[mode]][[
            lambda_key
          ]]$upper[
            validation_index
          ] <- fit$eif$signals$upper_signal

          eif_prediction <- fit$eif$signals
          eif_prediction$mode <- mode
          eif_prediction$fold <- fold

          eif_predictions[[
            eif_prediction_counter
          ]] <- eif_prediction

          eif_prediction_counter <-
            eif_prediction_counter + 1L

          eif_diagnostic <- fit$eif$diagnostics
          eif_diagnostic$mode <- mode
          eif_diagnostic$fold <- fold

          eif_diagnostics[[
            eif_diagnostic_counter
          ]] <- eif_diagnostic

          eif_diagnostic_counter <-
            eif_diagnostic_counter + 1L

          if (nrow(fit$eif$conditional_cv) > 0L) {
            local_conditional_cv <- fit$eif$conditional_cv
            local_conditional_cv$mode <- mode
            local_conditional_cv$fold <- fold
            local_conditional_cv$treatment <- t_value
            eif_conditional_cv_tables[[
              eif_conditional_cv_counter
            ]] <- local_conditional_cv
            eif_conditional_cv_counter <-
              eif_conditional_cv_counter + 1L
          }
        }
      }

      if (
        !is.null(gamma_transform) &&
          is.function(gamma_transform$diagnostics)
      ) {
        local_gamma_transform_diagnostics <-
          gamma_transform$diagnostics()
        if (!is.data.frame(local_gamma_transform_diagnostics)) {
          stop("gamma transformation diagnostics must be a data frame.")
        }
        local_gamma_transform_diagnostics$fold <- fold
        local_gamma_transform_diagnostics$treatment <- t_value
        gamma_transform_diagnostics[[
          gamma_transform_diagnostic_counter
        ]] <- local_gamma_transform_diagnostics
        gamma_transform_diagnostic_counter <-
          gamma_transform_diagnostic_counter + 1L
      }
    }
  }

  analysis_data$e_hat <- e_hat
  analysis_data$pi0_hat <- pi0_hat
  analysis_data$pi1_hat <- pi1_hat
  analysis_data$gamma0_rkhs <-
    gamma0_rkhs
  analysis_data$gamma1_rkhs <-
    gamma1_rkhs

  summary_rows <- list()
  summary_counter <- 1

  pseudo_summary_rows <- list()
  pseudo_summary_counter <- 1

  unconfounded_rows <- list()
  unconfounded_counter <- 1

  for (mode in modes) {
    mu1_lower <- mean(
      contribution_storage[[mode]][[
        "1"
      ]]$lower
    )

    mu1_upper <- mean(
      contribution_storage[[mode]][[
        "1"
      ]]$upper
    )

    mu0_lower <- mean(
      contribution_storage[[mode]][[
        "0"
      ]]$lower
    )

    mu0_upper <- mean(
      contribution_storage[[mode]][[
        "0"
      ]]$upper
    )

    summary_rows[[
      summary_counter
    ]] <- data.frame(
      mode = mode,
      mu1_lower = mu1_lower,
      mu1_upper = mu1_upper,
      mu0_lower = mu0_lower,
      mu0_upper = mu0_upper,
      ate_lower =
        mu1_lower -
        mu0_upper,
      ate_upper =
        mu1_upper -
        mu0_lower
    )

    summary_counter <-
      summary_counter + 1

    pseudo_mu1_lower <- mean(
      contribution_storage[[mode]][["1"]]$lower_pseudo
    )
    pseudo_mu1_upper <- mean(
      contribution_storage[[mode]][["1"]]$upper_pseudo
    )
    pseudo_mu0_lower <- mean(
      contribution_storage[[mode]][["0"]]$lower_pseudo
    )
    pseudo_mu0_upper <- mean(
      contribution_storage[[mode]][["0"]]$upper_pseudo
    )
    pseudo_summary_rows[[pseudo_summary_counter]] <- data.frame(
      mode = mode,
      mu1_lower = pseudo_mu1_lower,
      mu1_upper = pseudo_mu1_upper,
      mu0_lower = pseudo_mu0_lower,
      mu0_upper = pseudo_mu0_upper,
      ate_lower = pseudo_mu1_lower - pseudo_mu0_upper,
      ate_upper = pseudo_mu1_upper - pseudo_mu0_lower
    )
    pseudo_summary_counter <- pseudo_summary_counter + 1L

    unconfounded_rows[[
      unconfounded_counter
    ]] <- data.frame(
      mode = mode,
      mu1_no_unmeasured_confounding_benchmark = mean(
        contribution_storage[[mode]][[
          "1"
        ]]$lambda_one
      ),
      mu0_no_unmeasured_confounding_benchmark = mean(
        contribution_storage[[mode]][[
          "0"
        ]]$lambda_one
      )
    )

    unconfounded_counter <-
      unconfounded_counter + 1
  }

  prediction_table <- do.call(
    rbind,
    prediction_tables
  )

  fold_diagnostic_table <- do.call(
    rbind,
    fold_diagnostics
  )

  gamma_accuracy_table <- do.call(
    rbind,
    gamma_accuracy_rows
  )

  selected_lambda_table <- do.call(
    rbind,
    selected_lambdas
  )

  summary_table <- do.call(
    rbind,
    summary_rows
  )

  pseudo_summary_table <- do.call(
    rbind,
    pseudo_summary_rows
  )

  endpoint_names <- c(
    "mu1_lower",
    "mu1_upper",
    "mu0_lower",
    "mu0_upper",
    "ate_lower",
    "ate_upper"
  )

  eif_complete_modes <- eif_modes[vapply(
    eif_modes,
    function(mode) {
      all(vapply(
        eif_storage[[mode]],
        function(arm) {
          all(is.finite(arm$lower)) &&
            all(is.finite(arm$upper))
        },
        logical(1)
      ))
    },
    logical(1)
  )]

  one_step_signal_tables <- lapply(
    eif_complete_modes,
    function(mode) {
      local <- data.frame(
        row_id = seq_len(n),
        mode = mode,
        mu1_lower = eif_storage[[mode]][["1"]]$lower,
        mu1_upper = eif_storage[[mode]][["1"]]$upper,
        mu0_lower = eif_storage[[mode]][["0"]]$lower,
        mu0_upper = eif_storage[[mode]][["0"]]$upper,
        stringsAsFactors = FALSE
      )
      local$ate_lower <- local$mu1_lower - local$mu0_upper
      local$ate_upper <- local$mu1_upper - local$mu0_lower
      local
    }
  )

  one_step_signal_table <- if (length(one_step_signal_tables) > 0L) {
    do.call(rbind, one_step_signal_tables)
  } else {
    data.frame()
  }

  eif_prediction_table <- if (length(eif_predictions) > 0L) {
    do.call(rbind, eif_predictions)
  } else {
    data.frame()
  }

  eif_diagnostic_table <- if (length(eif_diagnostics) > 0L) {
    do.call(rbind, eif_diagnostics)
  } else {
    data.frame()
  }

  one_step_summary <- if (length(eif_complete_modes) > 0L) {
    do.call(rbind, lapply(eif_complete_modes, function(mode) {
      signals <- one_step_signal_table[
        one_step_signal_table$mode == mode,
        ,
        drop = FALSE
      ]
      plugin <- summary_table[
        summary_table$mode == mode,
        ,
        drop = FALSE
      ]
      local_diagnostics <- fold_diagnostic_table[
        fold_diagnostic_table$mode == mode,
        ,
        drop = FALSE
      ]
      local_eif_diagnostics <- eif_diagnostic_table[
        eif_diagnostic_table$mode == mode,
        ,
        drop = FALSE
      ]
      maximum_gamma_projection_rate <- max(
        local_eif_diagnostics$gamma_projection_rate_train,
        local_eif_diagnostics$gamma_projection_rate_validation
      )
      maximum_conditional_projection_rate <- max(
        local_diagnostics$conditional_projection_rate_train,
        local_diagnostics$conditional_projection_rate_validation
      )
      projection_active <- any(
        local_eif_diagnostics$gamma_projection_active
      )

      do.call(rbind, lapply(endpoint_names, function(endpoint) {
        signal <- signals[[endpoint]]
        estimate <- mean(signal)
        standard_error <- stats::sd(signal) / sqrt(length(signal))
        data.frame(
          mode = mode,
          gamma_nuisance = if (mode == "rkhs") {
            paste0(
              "estimated_",
              unique(selected_lambda_table$gamma_estimator)[1L]
            )
          } else {
            "oracle"
          },
          threshold_learner = threshold_learner,
          endpoint = endpoint,
          plugin_estimate = plugin[[endpoint]][1L],
          one_step_estimate = estimate,
          one_step_correction = estimate - plugin[[endpoint]][1L],
          standard_error = standard_error,
          confidence_level = 0.95,
          confidence_lower = estimate -
            stats::qnorm(0.975) * standard_error,
          confidence_upper = estimate +
            stats::qnorm(0.975) * standard_error,
          maximum_gamma_projection_rate =
            maximum_gamma_projection_rate,
          maximum_conditional_projection_rate =
            maximum_conditional_projection_rate,
          projected_nuisance_inference_flag = projection_active,
          inference_status = if (projection_active) {
            "if_with_projected_gamma_diagnostic"
          } else {
            "if_under_nuisance_rate_conditions"
          },
          stringsAsFactors = FALSE
        )
      }))
    }))
  } else {
    data.frame()
  }

  oracle_eif_signal_table <- if (
    "oracle_all" %in% eif_complete_modes
  ) {
    local <- one_step_signal_table[
      one_step_signal_table$mode == "oracle_all",
      ,
      drop = FALSE
    ]
    local$mode <- NULL
    local
  } else {
    data.frame()
  }

  oracle_one_step_summary <- if (nrow(one_step_summary) > 0L) {
    one_step_summary[
      one_step_summary$mode == "oracle_all",
      ,
      drop = FALSE
    ]
  } else {
    data.frame()
  }

  list(
    summary = summary_table,
    pseudo_summary = pseudo_summary_table,
    truth = data.frame(
      sample_mu1 = if (
        "Y1" %in% names(analysis_data)
      ) {
        mean(analysis_data$Y1)
      } else {
        NA_real_
      },
      sample_mu0 = if (
        "Y0" %in% names(analysis_data)
      ) {
        mean(analysis_data$Y0)
      } else {
        NA_real_
      },
      sample_ate = if (
        all(
          c("Y1", "Y0") %in%
            names(analysis_data)
        )
      ) {
        mean(
          analysis_data$Y1 -
            analysis_data$Y0
        )
      } else {
        NA_real_
      },
      population_mu1 = treatment_effect,
      population_mu0 = 0,
      population_ate = treatment_effect
    ),
    unconfounded_checks = do.call(
      rbind,
      unconfounded_rows
    ),
    gamma_accuracy_by_fold =
      gamma_accuracy_table,
    gamma_accuracy = stats::aggregate(
      cbind(
        bias,
        rmse,
        mae,
        correlation,
        mean_estimated,
        mean_oracle,
        rmse_constant_one,
        explained_fraction
      ) ~ treatment,
      data = gamma_accuracy_table,
      FUN = mean,
      na.action = stats::na.pass
    ),
    nuisance =
      compute_nuisance_diagnostics_simple(
        analysis_data
      ),
    nuisance_learner = nuisance_learner,
    nuisance_fit_diagnostics = do.call(
      rbind,
      nuisance_fit_diagnostics
    ),
    ipw_identities =
      compute_ipw_identities_simple(
        analysis_data
      ),
    selected_lambdas =
      selected_lambda_table,
    gamma_fit_provider_active =
      !is.null(gamma_fit_provider_factory),
    gamma_fit_diagnostics = if (length(gamma_fit_diagnostics) > 0L) {
      do.call(rbind, gamma_fit_diagnostics)
    } else {
      data.frame()
    },
    rkhs_cv = if (
      length(cv_tables) > 0
    ) {
      do.call(
        rbind,
        cv_tables
      )
    } else {
      data.frame()
    },
    threshold_cv = if (
      length(threshold_cv_tables) > 0L
    ) {
      do.call(rbind, threshold_cv_tables)
    } else {
      data.frame()
    },
    fixed_threshold_target_evaluation = if (
      length(fixed_threshold_target_evaluations) > 0L
    ) {
      do.call(rbind, fixed_threshold_target_evaluations)
    } else {
      data.frame()
    },
    fold_diagnostics =
      fold_diagnostic_table,
    one_step = one_step_summary,
    fully_estimated_one_step = if (nrow(one_step_summary) > 0L) {
      one_step_summary[
        one_step_summary$mode == "rkhs",
        ,
        drop = FALSE
      ]
    } else {
      data.frame()
    },
    one_step_signals = one_step_signal_table,
    eif_predictions = eif_prediction_table,
    eif_diagnostics = eif_diagnostic_table,
    eif_conditional_cv = if (length(eif_conditional_cv_tables) > 0L) {
      do.call(rbind, eif_conditional_cv_tables)
    } else {
      data.frame()
    },
    eif_conditional_provider_active =
      !is.null(eif_conditional_provider_factory),
    gamma_fit_transform_active =
      !is.null(gamma_fit_transform_factory),
    gamma_transform_diagnostics = if (
      length(gamma_transform_diagnostics) > 0L
    ) {
      do.call(rbind, gamma_transform_diagnostics)
    } else {
      data.frame()
    },
    oracle_one_step =
      oracle_one_step_summary,
    oracle_eif_signals =
      oracle_eif_signal_table,
    oracle_eif_predictions = if (
      nrow(eif_prediction_table) > 0L
    ) {
      eif_prediction_table[
        eif_prediction_table$mode == "oracle_all",
        ,
        drop = FALSE
      ]
    } else {
      data.frame()
    },
    oracle_eif_diagnostics = if (
      nrow(eif_diagnostic_table) > 0L
    ) {
      eif_diagnostic_table[
        eif_diagnostic_table$mode == "oracle_all",
        ,
        drop = FALSE
      ]
    } else {
      data.frame()
    },
    bound_predictions =
      prediction_table,
    normalization_by_x_bin =
      build_normalization_by_x_bin(
        predictions =
          prediction_table,
        x_bins = x_bins
      ),
    data = analysis_data,
    fold_id = fold_id,
    envelope_smoothing = data.frame(
      method = if (envelope_smoothing_epsilon > 0) {
        "conservative_log_mean_exp"
      } else {
        "none"
      },
      epsilon = envelope_smoothing_epsilon,
      bandwidth = if (envelope_smoothing_epsilon > 0) {
        envelope_smoothing_epsilon / log(2)
      } else {
        0
      },
      stringsAsFactors = FALSE
    )
  )
}

# Fully estimated one-step entry point for the current validation stage.
# The `rkhs` diagnostic mode still refers to estimation of gamma_t(Y,X);
# threshold_learner independently selects the q_t^- and q_t^+ learner.
estimate_sharp_bounds_one_step <- function(
    data,
    ...,
    compute_oracle_layers = TRUE
) {
  if (
    length(compute_oracle_layers) != 1L ||
      is.na(compute_oracle_layers)
  ) {
    stop("compute_oracle_layers must be TRUE or FALSE.")
  }
  arguments <- list(...)
  threshold_learner <- "spline_cv"
  if ("threshold_learner" %in% names(arguments)) {
    threshold_learner <- match.arg(
      arguments$threshold_learner,
      c(
        "spline_fixed", "spline_cv", "rkhs", "fixed_rbf",
        "growing_rbf_sieve", "external"
      )
    )
    arguments$threshold_learner <- NULL
  }
  reserved <- intersect(
    names(arguments),
    c(
      "compute_oracle_eif", "compute_estimated_eif", "one_step_modes",
      "estimation_modes"
    )
  )
  if (length(reserved) > 0L) {
    stop(
      "Use compute_oracle_layers rather than supplying ",
      paste(reserved, collapse = ", "),
      "."
    )
  }
  result <- do.call(
    estimate_sharp_bounds_plugin,
    c(
      list(
        data = data,
        threshold_learner = threshold_learner,
        estimation_modes = if (isTRUE(compute_oracle_layers)) {
          c("oracle_all", "oracle_gamma", "rkhs")
        } else {
          "rkhs"
        },
        one_step_modes = if (isTRUE(compute_oracle_layers)) {
          c("oracle_all", "oracle_gamma", "rkhs")
        } else {
          "rkhs"
        }
      ),
      arguments
    )
  )
  result$one_step_configuration <- list(
    threshold_learner = threshold_learner,
    fully_estimated_mode = "rkhs",
    nuisance_learner = if (
      "nuisance_learner" %in% names(arguments)
    ) {
      match.arg(
        arguments$nuisance_learner,
        c(
          "gam", "ranger", "glm_spline", "library_cv", "super_learner",
          "integrated_parametric"
        )
      )
    } else {
      "gam"
    },
    nuisance_gam_basis_selection = if (
      "nuisance_learner" %in% names(arguments) &&
        !identical(arguments$nuisance_learner, "gam")
    ) {
      NA_character_
    } else if (
      "nuisance_gam_k" %in% names(arguments) &&
      !is.null(arguments$nuisance_gam_k)
    ) {
      "fixed"
    } else {
      "adaptive_training_stratum_size"
    },
    nuisance_gam_min_k = if (
      "nuisance_gam_min_k" %in% names(arguments)
    ) arguments$nuisance_gam_min_k else 5L,
    nuisance_gam_max_k = if (
      "nuisance_gam_max_k" %in% names(arguments)
    ) arguments$nuisance_gam_max_k else 20L,
    nuisance_gam_growth_exponent = if (
      "nuisance_gam_growth_exponent" %in% names(arguments)
    ) arguments$nuisance_gam_growth_exponent else 1 / 3,
    gamma_estimator = if (
      "gamma_fit_provider_factory" %in% names(arguments) &&
        !is.null(arguments$gamma_fit_provider_factory)
    ) {
      unique(result$selected_lambdas$gamma_estimator)[1L]
    } else {
      "direct_rkhs_riesz"
    },
    gamma_rkhs_backend = if (
      "gamma_rkhs_backend" %in% names(arguments)
    ) {
      match.arg(arguments$gamma_rkhs_backend, c("exact", "nystrom"))
    } else {
      "exact"
    },
    gamma_rkhs_kernel = if (
      "gamma_rkhs_kernel" %in% names(arguments)
    ) {
      match.arg(
        arguments$gamma_rkhs_kernel,
        c("gaussian", "sobolev_matern32", "sobolev_matern52")
      )
    } else {
      "gaussian"
    },
    eif_conditional_learner = if (
      "eif_conditional_learner" %in% names(arguments)
    ) {
      match.arg(
        arguments$eif_conditional_learner,
        c(
          "spline_fixed", "spline_cv", "gam", "ranger", "library_cv"
        )
      )
    } else {
      "spline_fixed"
    },
    eif_conditional_spline_df = if (
      "eif_conditional_spline_df" %in% names(arguments)
    ) {
      arguments$eif_conditional_spline_df
    } else if ("conditional_mean_spline_df" %in% names(arguments)) {
      arguments$conditional_mean_spline_df
    } else {
      5L
    },
    custom_eif_instrument_probability =
      "eif_instrument_probability_transform" %in% names(arguments) &&
      is.function(arguments$eif_instrument_probability_transform),
    envelope_smoothing_epsilon = if (
      "envelope_smoothing_epsilon" %in% names(arguments)
    ) {
      arguments$envelope_smoothing_epsilon
    } else {
      0
    },
    oracle_layers_computed = isTRUE(compute_oracle_layers)
  )
  result
}

# =============================================================================
# 13. Recommended commands
# =============================================================================
#
# source("R/sharp_bounds_plugin.R")
#
# sim_data <- generate_iv_data(
#   n = 2000,
#   confounder_strength = 1,
#   seed = 2026
# )
#
# bound_results <- estimate_sharp_bounds_plugin(
#   data = sim_data,
#
#   Lambda0_T1 = 5,
#   Lambda1_T1 = 5,
#   Lambda0_T0 = 5,
#   Lambda1_T0 = 5,
#
#   x_cols = "X",
#   confounder_strength = 1,
#   K = 5,
#
#   positivity_epsilon = 0.05,
#   treatment_spline_df = 5,
#
#   gamma_floor = 1e-6,
#   gamma_cap = 100,
#
#   # NULL selects lambda separately by nested cross-validation within
#   # each outer fold and treatment arm. A named vector can still be
#   # supplied explicitly for fast debugging.
#   rkhs_lambda_by_t = NULL,
#
#   bandwidth_sample = 1000,
#   conditional_mean_spline_df = 5,
#   threshold_spline_df = 5,
#   # "spline_fixed" preserves the original implementation. Alternatives are
#   # "spline_cv" and "rkhs"; their tuning tables are returned in threshold_cv.
#   threshold_learner = "spline_fixed",
#
#   # Oracle mode remains strict; estimated modes use bounded projection.
#   estimated_compatibility_action = "project",
#   max_gamma_projection_rate = 0.02,
#
#   seed = 2026
# )
#
# bound_results$summary
# bound_results$pseudo_summary
# bound_results$threshold_cv
# bound_results$oracle_one_step
# bound_results$oracle_eif_diagnostics
# bound_results$truth
# bound_results$unconfounded_checks
#
# bound_results$gamma_accuracy
# bound_results$nuisance$e
# bound_results$nuisance$pi0
# bound_results$nuisance$pi1
# bound_results$nuisance$w1
# bound_results$nuisance$w0
# bound_results$ipw_identities
#
# bound_results$selected_lambdas
#
# aggregate(
#   cbind(
#     gamma_projection_rate_validation,
#     conditional_projection_rate_validation,
#     tau_invalid_rate_validation,
#     mean_lambda_lower_validation,
#     mean_lambda_upper_validation,
#     fold_lower_mean,
#     fold_upper_mean,
#     fold_lambda_one_mean
#   ) ~ mode + treatment,
#   data = bound_results$fold_diagnostics,
#   FUN = mean
# )
#
# bound_results$normalization_by_x_bin
#
# subset(
#   bound_results$normalization_by_x_bin,
#   mode == "oracle_all" &
#     treatment == 1
# )
#
# subset(
#   bound_results$normalization_by_x_bin,
#   mode == "rkhs" &
#     treatment == 1
# )
