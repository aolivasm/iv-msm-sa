# Correctly specified, estimated parametric nuisance layer for the observed
# continuous-X5 simulation law.
#
# This is a simulation benchmark, not an oracle layer. Every outer training
# fold estimates the coefficients of P(Z | X), P(T | X,Z), the Z=0 outcome
# location, and the Z=1 exponential tilt. The remaining nuisance functions
# (gamma, conditional envelope means, thresholds, and EIF conditional means)
# are then induced by that fitted observed-data law using deterministic
# one-dimensional quadrature.

observed_x5_parametric_design <- function(data, component) {
  required <- paste0("X", 1:5)
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Missing X5 covariates: ", paste(missing, collapse = ", "), ".")
  }
  x1 <- data$X1
  x2 <- data$X2
  x3 <- data$X3
  x4 <- data$X4
  x5 <- data$X5
  switch(
    component,
    instrument = cbind(
      `(Intercept)` = 1, X1 = x1, X2 = x2, X3 = x3
    ),
    treatment = cbind(
      `(Intercept)` = 1, X1 = x1, X2 = x2, X3 = x3,
      X1_X2 = x1 * x2, sin_X4 = sin(x4), Z = data$Z
    ),
    outcome = cbind(
      `(Intercept)` = 1, T = data$T_obs, X1 = x1, X2 = x2, X3 = x3,
      X4_X5 = x4 * x5, X1_sq_minus_1 = x1^2 - 1
    ),
    stop("Unknown X5 parametric design component: ", component, ".")
  )
}

observed_x5_parametric_tilt_shape <- function(data) {
  sin((data$X1 + data$X2 - data$X3) / sqrt(3))^2
}

observed_x5_parametric_truncnorm_mean <- function(
    mean, lower = -2, upper = 2
) {
  mean + (
    stats::dnorm(lower - mean) - stats::dnorm(upper - mean)
  ) / (
    stats::pnorm(upper - mean) - stats::pnorm(lower - mean)
  )
}

observed_x5_parametric_qtruncnorm <- function(
    probability, mean = 0, lower = -2, upper = 2
) {
  probability <- pmin(1 - 1e-12, pmax(1e-12, probability))
  mean <- rep_len(as.numeric(mean), length(probability))
  lower_probability <- stats::pnorm(lower - mean)
  upper_probability <- stats::pnorm(upper - mean)
  mean + stats::qnorm(
    lower_probability + probability * (upper_probability - lower_probability)
  )
}

fit_observed_x5_parametric_law <- function(
    train, positivity_epsilon = 0.05, tilt_bound = 0.25
) {
  if (!all(c("Y", "T_obs", "Z", paste0("X", 1:5)) %in% names(train))) {
    stop("The continuous-X5 parametric law received incomplete training data.")
  }

  instrument_design <- observed_x5_parametric_design(train, "instrument")
  instrument_fit <- stats::glm.fit(
    x = instrument_design, y = train$Z, family = stats::binomial(),
    control = stats::glm.control(maxit = 100)
  )

  treatment_design <- observed_x5_parametric_design(train, "treatment")
  treatment_fit <- stats::glm.fit(
    x = treatment_design, y = train$T_obs, family = stats::binomial(),
    control = stats::glm.control(maxit = 100)
  )

  baseline <- train$Z == 0
  outcome_design <- observed_x5_parametric_design(train, "outcome")
  outcome_fit <- stats::lm.fit(
    x = outcome_design[baseline, , drop = FALSE],
    y = train$Y[baseline]
  )
  outcome_coefficients <- outcome_fit$coefficients
  outcome_coefficients[!is.finite(outcome_coefficients)] <- 0
  residual <- train$Y - as.vector(outcome_design %*% outcome_coefficients)

  comparison <- train$Z == 1
  comparison_residual <- residual[comparison]
  comparison_sign <- 2 * train$T_obs[comparison] - 1
  comparison_shape <- observed_x5_parametric_tilt_shape(train[comparison, ])
  tilt_design <- cbind(
    intercept = comparison_sign,
    shape = comparison_sign * comparison_shape
  )
  tilt_objective <- function(coefficient) {
    eta <- as.vector(tilt_design %*% coefficient)
    normalizer <- stats::pnorm(2 - eta) - stats::pnorm(-2 - eta)
    -sum(-0.5 * (comparison_residual - eta)^2 - log(normalizer))
  }
  tilt_fit <- stats::optim(
    par = c(intercept = 0.03, shape = 0.01),
    fn = tilt_objective,
    method = "L-BFGS-B",
    lower = rep(-tilt_bound, 2L),
    upper = rep(tilt_bound, 2L),
    control = list(factr = 1e7, pgtol = 1e-10, maxit = 500L)
  )

  instrument_coefficient <- instrument_fit$coefficients
  treatment_coefficient <- treatment_fit$coefficients
  instrument_coefficient[!is.finite(instrument_coefficient)] <- 0
  treatment_coefficient[!is.finite(treatment_coefficient)] <- 0

  predict_nuisance <- function(newdata) {
    instrument <- observed_x5_parametric_design(newdata, "instrument")
    treatment0 <- newdata
    treatment1 <- newdata
    treatment0$Z <- 0
    treatment1$Z <- 1
    p_z1 <- stats::plogis(as.vector(instrument %*% instrument_coefficient))
    pi0 <- stats::plogis(as.vector(
      observed_x5_parametric_design(treatment0, "treatment") %*%
        treatment_coefficient
    ))
    pi1 <- stats::plogis(as.vector(
      observed_x5_parametric_design(treatment1, "treatment") %*%
        treatment_coefficient
    ))
    list(
      e = clip_probability(1 - p_z1, positivity_epsilon),
      pi0 = clip_probability(pi0, positivity_epsilon),
      pi1 = clip_probability(pi1, positivity_epsilon)
    )
  }

  predict_location <- function(newdata, treatment = NULL) {
    local <- newdata
    if (!is.null(treatment)) local$T_obs <- rep_len(treatment, nrow(local))
    as.vector(
      observed_x5_parametric_design(local, "outcome") %*%
        outcome_coefficients
    )
  }

  predict_tilt <- function(newdata, treatment = NULL) {
    local_treatment <- if (is.null(treatment)) {
      newdata$T_obs
    } else {
      rep_len(treatment, nrow(newdata))
    }
    (2 * local_treatment - 1) * (
      tilt_fit$par[["intercept"]] +
        tilt_fit$par[["shape"]] * observed_x5_parametric_tilt_shape(newdata)
    )
  }

  predict_gamma <- function(newdata, treatment = NULL) {
    location <- predict_location(newdata, treatment)
    eta <- predict_tilt(newdata, treatment)
    residual_value <- newdata$Y - location
    base_normalizer <- stats::pnorm(2) - stats::pnorm(-2)
    tilted_normalizer <- stats::pnorm(2 - eta) - stats::pnorm(-2 - eta)
    exp(eta * residual_value - eta^2 / 2) *
      base_normalizer / tilted_normalizer
  }

  list(
    predict_nuisance = predict_nuisance,
    predict_location = predict_location,
    predict_tilt = predict_tilt,
    predict_gamma = predict_gamma,
    coefficients = list(
      instrument = instrument_coefficient,
      treatment = treatment_coefficient,
      outcome = outcome_coefficients,
      tilt = tilt_fit$par
    ),
    diagnostics = data.frame(
      component = c("e", "pi0", "pi1", "outcome_location", "gamma_tilt"),
      learner = "correct_parametric_observed_law",
      converged = c(
        isTRUE(instrument_fit$converged),
        rep(isTRUE(treatment_fit$converged), 2L),
        outcome_fit$rank == ncol(outcome_design),
        identical(tilt_fit$convergence, 0L)
      ),
      model_rank = c(
        instrument_fit$rank, rep(treatment_fit$rank, 2L),
        outcome_fit$rank, 2L
      ),
      effective_df = c(
        instrument_fit$rank, rep(treatment_fit$rank, 2L),
        outcome_fit$rank, 2L
      ),
      n_training = c(
        nrow(train), nrow(train), nrow(train), sum(baseline), sum(comparison)
      ),
      stringsAsFactors = FALSE
    )
  )
}

make_observed_x5_parametric_provider_bundle <- function(
    quadrature_nodes = 24L, positivity_epsilon = 0.05
) {
  quadrature_nodes <- as.integer(quadrature_nodes)
  if (quadrature_nodes < 8L) stop("quadrature_nodes must be at least eight.")
  rule <- observed_dgp_gauss_legendre(quadrature_nodes, -2, 2)
  cache <- new.env(parent = emptyenv())

  get_fit <- function(fold) {
    key <- paste0("fold_", as.integer(fold))
    if (!exists(key, envir = cache, inherits = FALSE)) {
      stop("The parametric observed-law fit is missing for fold ", fold, ".")
    }
    get(key, envir = cache, inherits = FALSE)
  }

  nuisance_provider <- function(fold, train, x_cols, seed, ...) {
    if (!identical(x_cols, paste0("X", 1:5))) {
      stop("The correctly specified X5 provider requires X1,...,X5.")
    }
    key <- paste0("fold_", as.integer(fold))
    if (!exists(key, envir = cache, inherits = FALSE)) {
      assign(
        key,
        fit_observed_x5_parametric_law(
          train, positivity_epsilon = positivity_epsilon
        ),
        envir = cache
      )
    }
    fit <- get_fit(fold)
    list(
      predict = fit$predict_nuisance,
      diagnostics = fit$diagnostics,
      learner = "correct_parametric_observed_law",
      parametric_law_fit = fit
    )
  }

  gamma_provider <- function(
      fold, treatment, train, nuisance_train, x_cols, seed, ...
  ) {
    fit <- get_fit(fold)
    predict_gamma <- function(newdata) {
      value <- as.numeric(fit$predict_gamma(newdata, treatment))
      list(raw = value, positive = value, normalized = value)
    }
    list(
      predict = predict_gamma,
      gamma_estimator = "correct_parametric_density_ratio",
      selected_lambda = NA_real_,
      selected_bandwidth_multiplier = NA_real_,
      base_bandwidth = NA_real_,
      selected_bandwidth = NA_real_,
      normalization_constant = 1,
      training_negative_rate = 0,
      rkhs_backend = "parametric_observed_law",
      rkhs_kernel = NA_character_,
      nystrom_rank_requested = NA_integer_,
      nystrom_landmark_count = NA_integer_,
      nystrom_effective_rank = NA_integer_,
      nystrom_min_retained_eigenvalue = NA_real_,
      nystrom_max_eigenvalue = NA_real_,
      cv_table = data.frame()
    )
  }

  expand_rows <- function(newdata, treatment, z, fit) {
    n <- nrow(newdata)
    index <- rep(seq_len(n), each = quadrature_nodes)
    expanded <- newdata[index, , drop = FALSE]
    expanded$T_obs <- treatment
    expanded$Z <- z
    residual <- rep(rule$nodes, times = n)
    expanded$Y <- fit$predict_location(expanded, treatment) + residual
    list(data = expanded, index = index, residual = residual, n = n)
  }

  row_integral <- function(value, weight, n) {
    value_matrix <- matrix(value, nrow = n, byrow = TRUE)
    weight_matrix <- matrix(weight, nrow = n, byrow = TRUE)
    weight_matrix <- weight_matrix / rowSums(weight_matrix)
    rowSums(value_matrix * weight_matrix)
  }

  conditional_envelope_means <- function(
      newdata, treatment, Lambda0, Lambda1,
      envelope_smoothing_epsilon, compatibility_action,
      gamma_projection_margin, fit
  ) {
    expanded <- expand_rows(newdata, treatment, 0L, fit)
    nuisance <- fit$predict_nuisance(expanded$data)
    gamma <- fit$predict_gamma(expanded$data, treatment)
    envelope <- compute_projected_envelopes(
      nuisance = nuisance,
      gamma = gamma,
      t_value = treatment,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      envelope_smoothing_epsilon = envelope_smoothing_epsilon,
      gamma_compatibility_action = compatibility_action,
      gamma_projection_margin = gamma_projection_margin
    )
    base_weight <- rule$weights * observed_dgp_dtruncnorm(rule$nodes, 0)
    weight <- rep(base_weight, times = expanded$n)
    list(
      Lbar = row_integral(envelope$L, weight, expanded$n),
      Ubar = row_integral(envelope$U, weight, expanded$n)
    )
  }

  conditional_weighted_threshold <- function(
      newdata, treatment, probability, Lambda0, Lambda1,
      envelope_smoothing_epsilon, compatibility_action,
      gamma_projection_margin, fit
  ) {
    expanded <- expand_rows(newdata, treatment, 0L, fit)
    nuisance <- fit$predict_nuisance(expanded$data)
    gamma <- fit$predict_gamma(expanded$data, treatment)
    envelope <- compute_projected_envelopes(
      nuisance = nuisance,
      gamma = gamma,
      t_value = treatment,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      envelope_smoothing_epsilon = envelope_smoothing_epsilon,
      gamma_compatibility_action = compatibility_action,
      gamma_projection_margin = gamma_projection_margin
    )
    base_weight <- rule$weights * observed_dgp_dtruncnorm(rule$nodes, 0)
    weighted_mass <- matrix(
      envelope$D * rep(base_weight, times = expanded$n),
      nrow = expanded$n,
      byrow = TRUE
    )
    weighted_mass <- weighted_mass / rowSums(weighted_mass)
    probability <- rep_len(as.numeric(probability), expanded$n)
    residual_quantile <- vapply(seq_len(expanded$n), function(index) {
      mass <- weighted_mass[index, ]
      # The half-mass locations turn the Gaussian-quadrature masses into a
      # continuous CDF approximation. Interpolation avoids the artificial
      # jumps obtained by treating quadrature nodes as a discrete outcome.
      cdf_at_node <- cumsum(mass) - mass / 2
      stats::approx(
        x = c(0, cdf_at_node, 1),
        y = c(-2, rule$nodes, 2),
        xout = probability[index],
        ties = "ordered",
        rule = 2
      )$y
    }, numeric(1))
    fit$predict_location(newdata, treatment) + residual_quantile
  }

  bound_component_factory <- function(
      mode, fold, treatment, gamma_fit, nuisance_fit, x_cols
  ) {
    fit <- get_fit(fold)
    envelope_means <- function(
        train, validation, nuisance_train, nuisance_validation,
        envelope_train, envelope_validation, treatment, Lambda0, Lambda1,
        compatibility_action, gamma_projection_margin,
        envelope_smoothing_epsilon
    ) {
      train_mean <- conditional_envelope_means(
        train, treatment, Lambda0, Lambda1, envelope_smoothing_epsilon,
        compatibility_action, gamma_projection_margin, fit
      )
      validation_mean <- conditional_envelope_means(
        validation, treatment, Lambda0, Lambda1, envelope_smoothing_epsilon,
        compatibility_action, gamma_projection_margin, fit
      )
      list(
        Lbar_train = train_mean$Lbar,
        Ubar_train = train_mean$Ubar,
        Lbar_validation = validation_mean$Lbar,
        Ubar_validation = validation_mean$Ubar
      )
    }

    thresholds <- function(
        train, validation, nuisance_train, nuisance_validation,
        envelope_train, envelope_validation, treatment, Lambda0, Lambda1,
        tau_train, tau_validation, compatibility_action,
        gamma_projection_margin, envelope_smoothing_epsilon, tau_epsilon
    ) {
      threshold_from_tau <- function(newdata, tau, endpoint) {
        probability <- if (endpoint == "lower") 1 - tau else tau
        conditional_weighted_threshold(
          newdata = newdata,
          treatment = treatment,
          probability = probability,
          Lambda0 = Lambda0,
          Lambda1 = Lambda1,
          envelope_smoothing_epsilon = envelope_smoothing_epsilon,
          compatibility_action = compatibility_action,
          gamma_projection_margin = gamma_projection_margin,
          fit = fit
        )
      }
      predict_endpoint <- function(endpoint) {
        force(endpoint)
        function(newdata) {
          means <- conditional_envelope_means(
            newdata, treatment, Lambda0, Lambda1,
            envelope_smoothing_epsilon, compatibility_action,
            gamma_projection_margin, fit
          )
          tau <- (means$Ubar - 1) / (means$Ubar - means$Lbar)
          tau <- pmin(1 - tau_epsilon, pmax(tau_epsilon, tau))
          threshold_from_tau(newdata, tau, endpoint)
        }
      }
      list(
        q_lower_train = threshold_from_tau(train, tau_train, "lower"),
        q_upper_train = threshold_from_tau(train, tau_train, "upper"),
        q_lower_validation = threshold_from_tau(
          validation, tau_validation, "lower"
        ),
        q_upper_validation = threshold_from_tau(
          validation, tau_validation, "upper"
        ),
        lower_predict = predict_endpoint("lower"),
        upper_predict = predict_endpoint("upper")
      )
    }
    list(envelope_means = envelope_means, thresholds = thresholds)
  }

  eif_provider_factory <- function(
      mode, fold, treatment, gamma_fit, x_cols
  ) {
    fit <- get_fit(fold)
    function(
        validation, nuisance_validation, envelope_validation, t_value,
        Lambda0, Lambda1, q_lower_validation, q_upper_validation,
        baseline_outcomes, comparison_outcomes
    ) {
      smoothing_epsilon <- envelope_validation$envelope_smoothing_epsilon
      if (is.null(smoothing_epsilon)) smoothing_epsilon <- 0

      integrate_primitives <- function(z) {
        expanded <- expand_rows(validation, treatment, z, fit)
        nuisance <- lapply(nuisance_validation, function(value) {
          value[expanded$index]
        })
        gamma <- fit$predict_gamma(expanded$data, treatment)
        envelope <- compute_projected_envelopes(
          nuisance = nuisance,
          gamma = gamma,
          t_value = treatment,
          Lambda0 = Lambda0,
          Lambda1 = Lambda1,
          envelope_smoothing_epsilon = smoothing_epsilon,
          gamma_compatibility_action = "project",
          gamma_projection_margin = 1e-4
        )
        primitive <- compute_eif_row_primitives(
          data = expanded$data,
          nuisance = nuisance,
          envelope = envelope,
          t_value = treatment,
          Lambda0 = Lambda0,
          Lambda1 = Lambda1,
          q_lower = q_lower_validation[expanded$index],
          q_upper = q_upper_validation[expanded$index]
        )
        eta <- if (z == 0L) 0 else fit$predict_tilt(expanded$data, treatment)
        weight <- rule$weights[rep(seq_len(quadrature_nodes), times = expanded$n)] *
          observed_dgp_dtruncnorm(expanded$residual, eta)
        integrate <- function(value) row_integral(value, weight, expanded$n)
        list(primitive = primitive, integrate = integrate, expanded = expanded)
      }

      baseline <- integrate_primitives(0L)
      comparison <- integrate_primitives(1L)
      baseline_prediction <- matrix(
        NA_real_, nrow(validation), length(baseline_outcomes),
        dimnames = list(NULL, baseline_outcomes)
      )
      comparison_prediction <- matrix(
        NA_real_, nrow(validation), length(comparison_outcomes),
        dimnames = list(NULL, comparison_outcomes)
      )
      baseline_prediction[, ".m0"] <- baseline$integrate(
        baseline$expanded$data$Y
      )
      for (endpoint in c("lower", "upper")) {
        primitive0 <- baseline$primitive[[endpoint]]
        primitive1 <- comparison$primitive[[endpoint]]
        baseline_prediction[, paste0(".theta_", endpoint)] <-
          baseline$integrate(primitive0$lambda * baseline$expanded$data$Y)
        baseline_prediction[, paste0(".D_r0_", endpoint)] <-
          baseline$integrate(primitive0$D_r0)
        baseline_prediction[, paste0(".D_r1_", endpoint)] <-
          baseline$integrate(primitive0$D_r1)
        baseline_prediction[, paste0(".H0_", endpoint)] <-
          baseline$integrate(primitive0$H0)
        comparison_prediction[, paste0(".H1_", endpoint)] <-
          comparison$integrate(primitive1$H1)
      }
      list(
        baseline_prediction = baseline_prediction,
        comparison_prediction = comparison_prediction,
        learner = "correct_parametric_quadrature",
        diagnostics = data.frame(
          parametric_quadrature_nodes = quadrature_nodes,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  list(
    nuisance_fit_provider_factory = nuisance_provider,
    gamma_fit_provider_factory = gamma_provider,
    bound_component_provider_factory = bound_component_factory,
    eif_conditional_provider_factory = eif_provider_factory,
    cache = cache,
    quadrature_nodes = quadrature_nodes
  )
}
