# Fixed-sieve compatibility test for the Dual Marginal Sensitivity Model.
#
# This file implements a constrained-versus-unrestricted minimum-Riesz-risk
# prototype.  It is deliberately separate from R/compatibility_diagnostic.R:
# the latter counts violations of a previously fitted gamma_t and remains a
# descriptive diagnostic, whereas this file re-estimates gamma_t under the
# candidate Dual-MSM restrictions.
#
# Required before sourcing this file:
#   source("R/sharp_bounds_plugin.R")
#
# Inferential scope:
#   * Pointwise gamma restrictions are imposed on the treatment-arm sample.
#   * Conditional envelope-normalization inequalities are approximated by a
#     finite collection of nonnegative box instruments in X.
#   * A null-imposed multiplier bootstrap is supplied for the resulting fixed
#     sieve problem.
#   * The bootstrap holds primitive nuisance estimates and the sieve fixed.
#
# The final item is important.  This is a calibrated fixed-sieve prototype,
# not yet the paper's final uniform growing-sieve test.  A final theorem must
# account for sieve approximation and estimated primitive nuisances, following
# the constrained conditional-moment framework described in the accompanying
# documentation.

check_formal_compatibility_dependencies <- function() {
  required_functions <- c(
    "prepare_rkhs_coordinates",
    "prepare_gamma_nystrom_landmarks",
    "make_gamma_nystrom_features",
    "compute_riesz_weight",
    "compute_compatible_gamma_interval",
    "fit_nuisance_models",
    "make_tz_stratified_folds"
  )
  missing <- required_functions[
    !vapply(required_functions, exists, logical(1), mode = "function")
  ]
  if (length(missing) > 0L) {
    stop(
      "Source R/sharp_bounds_plugin.R before this file. Missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop(
      "The formal compatibility prototype requires the quadprog package.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_compatibility_nuisance <- function(nuisance, n) {
  required <- c("e", "pi0", "pi1")
  missing <- setdiff(required, names(nuisance))
  if (length(missing) > 0L) {
    stop("nuisance omits: ", paste(missing, collapse = ", "))
  }
  invalid <- vapply(required, function(name) {
    length(nuisance[[name]]) != n ||
      any(!is.finite(nuisance[[name]])) ||
      any(nuisance[[name]] <= 0 | nuisance[[name]] >= 1)
  }, logical(1))
  if (any(invalid)) {
    stop(
      "Every nuisance component must contain n finite probabilities strictly ",
      "between zero and one. Invalid: ",
      paste(required[invalid], collapse = ", ")
    )
  }
  invisible(TRUE)
}

crossfit_compatibility_nuisance <- function(
    data,
    x_cols = "X",
    folds = 5L,
    learner = c("glm_spline", "gam", "integrated_parametric"),
    positivity_epsilon = 0.05,
    treatment_spline_df = 5L,
    gam_basis_dimension = 6L,
    integrated_confounder_strength = 1,
    integrated_quadrature_nodes = 40L,
    seed = 1L
) {
  check_formal_compatibility_dependencies()
  learner <- match.arg(learner)
  fold_id <- make_tz_stratified_folds(data, K = folds, seed = seed)
  output <- list(
    e = rep(NA_real_, nrow(data)),
    pi0 = rep(NA_real_, nrow(data)),
    pi1 = rep(NA_real_, nrow(data))
  )
  diagnostics <- vector("list", folds)
  for (fold in seq_len(folds)) {
    training <- fold_id != fold
    validation <- !training
    fit <- fit_nuisance_models(
      train = data[training, , drop = FALSE],
      x_cols = x_cols,
      positivity_epsilon = positivity_epsilon,
      learner = learner,
      treatment_spline_df = treatment_spline_df,
      gam_basis_dimension = gam_basis_dimension,
      integrated_confounder_strength = integrated_confounder_strength,
      integrated_quadrature_nodes = integrated_quadrature_nodes
    )
    prediction <- fit$predict(data[validation, , drop = FALSE])
    for (name in names(output)) {
      output[[name]][validation] <- prediction[[name]]
    }
    diagnostic <- fit$diagnostics
    diagnostic$fold <- fold
    diagnostics[[fold]] <- diagnostic
  }
  validate_compatibility_nuisance(output, nrow(data))
  output$fold_id <- fold_id
  output$diagnostics <- do.call(rbind, diagnostics)
  output$learner <- learner
  output
}

make_multiscale_box_instruments <- function(
    data,
    x_cols = "X",
    resolutions = c(2L, 4L),
    minimum_cell_count = 5L
) {
  missing <- setdiff(x_cols, names(data))
  if (length(missing) > 0L) {
    stop("Instrument data omit: ", paste(missing, collapse = ", "))
  }
  if (nrow(data) < 2L) {
    stop("At least two baseline observations are required.")
  }
  resolutions <- sort(unique(as.integer(resolutions)))
  if (
    length(resolutions) == 0L ||
      any(!is.finite(resolutions)) ||
      any(resolutions < 2L)
  ) {
    stop("resolutions must contain integers of at least two.")
  }

  columns <- list(one = rep(1, nrow(data)))
  metadata <- list(data.frame(
    instrument = "one",
    covariate = "all",
    resolution = 1L,
    cell = 1L,
    cell_count = nrow(data),
    stringsAsFactors = FALSE
  ))
  for (x_name in x_cols) {
    x <- data[[x_name]]
    if (!is.numeric(x) || any(!is.finite(x))) {
      stop("Box-instrument covariates must be finite numeric variables.")
    }
    for (resolution in resolutions) {
      breaks <- unique(as.numeric(stats::quantile(
        x,
        probs = seq(0, 1, length.out = resolution + 1L),
        type = 8,
        names = FALSE
      )))
      if (length(breaks) < 3L) {
        next
      }
      breaks[1L] <- -Inf
      breaks[length(breaks)] <- Inf
      cell <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)
      for (cell_value in sort(unique(cell))) {
        indicator <- as.numeric(cell == cell_value)
        cell_count <- sum(indicator)
        if (cell_count < minimum_cell_count) {
          next
        }
        name <- paste0(x_name, "_r", resolution, "_c", cell_value)
        columns[[name]] <- indicator
        metadata[[length(metadata) + 1L]] <- data.frame(
          instrument = name,
          covariate = x_name,
          resolution = resolution,
          cell = cell_value,
          cell_count = cell_count,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  matrix <- do.call(cbind, columns)
  storage.mode(matrix) <- "double"
  colnames(matrix) <- names(columns)
  list(
    matrix = matrix,
    metadata = do.call(rbind, metadata),
    description = paste0(
      "nonnegative marginal quantile boxes at resolutions ",
      paste(resolutions, collapse = ",")
    )
  )
}

make_compatibility_gamma_sieve <- function(
    arm_data,
    x_cols = "X",
    rank = 25L,
    bandwidth_multiplier = 4,
    bandwidth_sample = 1000L,
    seed = 1L,
    eigen_tolerance = 1e-10,
    numerical_jitter = 1e-10
) {
  check_formal_compatibility_dependencies()
  coordinate_spec <- prepare_rkhs_coordinates(
    train_t = arm_data,
    x_cols = x_cols,
    bandwidth_sample = bandwidth_sample,
    seed = seed
  )
  landmark_spec <- prepare_gamma_nystrom_landmarks(
    coordinate_spec = coordinate_spec,
    rank = rank,
    seed = seed + 1000L
  )
  feature_spec <- make_gamma_nystrom_features(
    coordinate_spec = coordinate_spec,
    landmark_spec = landmark_spec,
    bandwidth_multiplier = bandwidth_multiplier,
    eigen_tolerance = eigen_tolerance,
    numerical_jitter = numerical_jitter
  )
  raw <- cbind(intercept = 1, feature_spec$features_train)
  scale <- sqrt(colMeans(raw^2))
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  scale[1L] <- 1
  features <- sweep(raw, 2L, scale, "/")
  colnames(features) <- c(
    "intercept",
    paste0("nystrom_", seq_len(ncol(features) - 1L))
  )
  list(
    features = features,
    scale = scale,
    coordinate_spec = coordinate_spec,
    landmark_spec = landmark_spec,
    feature_spec = feature_spec,
    rank_requested = feature_spec$rank_requested,
    effective_rank = feature_spec$effective_rank,
    base_bandwidth = coordinate_spec$sigma,
    bandwidth_multiplier = bandwidth_multiplier,
    selected_bandwidth = feature_spec$sigma
  )
}

envelope_affine_coefficients <- function(
    nuisance,
    t_value,
    Lambda0,
    Lambda1
) {
  pi0 <- nuisance$pi0
  pi1 <- nuisance$pi1
  if (t_value == 1L) {
    coefficient_A <- (1 - pi1) / (1 - pi0)
    coefficient_B <- pi1 / (1 - pi0)
    intercept_B <- -pi0 / (1 - pi0)
  } else if (t_value == 0L) {
    coefficient_A <- pi1 / pi0
    coefficient_B <- (1 - pi1) / pi0
    intercept_B <- -(1 - pi0) / pi0
  } else {
    stop("t_value must equal zero or one.")
  }
  list(
    lower_slope = coefficient_A / Lambda1 + coefficient_B,
    upper_slope = coefficient_A * Lambda1 + coefficient_B,
    shared_intercept = intercept_B,
    lower_constant = 1 / Lambda0,
    upper_constant = Lambda0
  )
}

evaluate_sieve_normalization_moments <- function(problem, beta) {
  gamma <- as.vector(problem$features_baseline %*% beta)
  affine <- problem$envelope_affine
  lower_affine <- affine$lower_slope * gamma + affine$shared_intercept
  upper_affine <- affine$upper_slope * gamma + affine$shared_intercept
  lower_active <- lower_affine >= affine$lower_constant
  upper_active <- upper_affine <= affine$upper_constant
  L <- ifelse(lower_active, lower_affine, affine$lower_constant)
  U <- ifelse(upper_active, upper_affine, affine$upper_constant)
  H <- problem$moment_instruments
  lower <- colMeans(H * (L - 1))
  upper <- colMeans(H * (1 - U))
  list(
    lower = lower,
    upper = upper,
    maximum = max(c(lower, upper)),
    L = L,
    U = U,
    gamma = gamma,
    lower_active = lower_active,
    upper_active = upper_active,
    lower_affine = lower_affine,
    upper_affine = upper_affine
  )
}

make_normalization_cut <- function(problem, beta, side, instrument_index) {
  evaluation <- evaluate_sieve_normalization_moments(problem, beta)
  H <- problem$moment_instruments[, instrument_index]
  Phi <- problem$features_baseline
  affine <- problem$envelope_affine
  if (side == "lower") {
    active <- as.numeric(evaluation$lower_active)
    slope <- colMeans(Phi * (H * active * affine$lower_slope))
    intercept <- mean(H * (
      active * affine$shared_intercept +
        (1 - active) * affine$lower_constant - 1
    ))
  } else if (side == "upper") {
    active <- as.numeric(evaluation$upper_active)
    slope <- -colMeans(Phi * (H * active * affine$upper_slope))
    intercept <- mean(H * (
      1 - active * affine$shared_intercept -
        (1 - active) * affine$upper_constant
    ))
  } else {
    stop("side must be 'lower' or 'upper'.")
  }
  # The active affine piece has form slope' beta + intercept <= 0.
  # solve.QP uses C beta >= b, hence C=-slope and b=intercept.
  list(
    C = -slope,
    b = intercept,
    side = side,
    instrument_index = instrument_index,
    instrument = colnames(problem$moment_instruments)[instrument_index]
  )
}

deduplicate_linear_constraints <- function(C, b, digits = 11L) {
  if (nrow(C) == 0L) {
    return(list(C = C, b = b))
  }
  key <- apply(
    cbind(C, b),
    1L,
    function(row) paste(format(round(row, digits), scientific = TRUE), collapse = "|")
  )
  keep <- !duplicated(key)
  list(C = C[keep, , drop = FALSE], b = b[keep])
}

solve_compatibility_qp <- function(
    problem,
    dvec = problem$riesz_linear,
    initial_C = problem$pointwise_C,
    initial_b = problem$pointwise_b,
    normalization_tolerance = 1e-8,
    maximum_cutting_iterations = 60L
) {
  C <- initial_C
  b <- initial_b
  cut_log <- list()
  solve_once <- function(C, b) {
    result <- tryCatch(
      quadprog::solve.QP(
        Dmat = problem$riesz_hessian,
        dvec = dvec,
        Amat = t(C),
        bvec = b,
        meq = 0L,
        factorized = FALSE
      ),
      error = function(error) error
    )
    if (inherits(result, "error")) {
      stop(
        "The fixed-sieve compatibility quadratic program is infeasible or ",
        "numerically singular: ", conditionMessage(result),
        call. = FALSE
      )
    }
    result
  }

  for (iteration in seq_len(maximum_cutting_iterations)) {
    deduplicated <- deduplicate_linear_constraints(C, b)
    C <- deduplicated$C
    b <- deduplicated$b
    fit <- solve_once(C, b)
    beta <- as.vector(fit$solution)
    moments <- evaluate_sieve_normalization_moments(problem, beta)
    lower_violations <- which(moments$lower > normalization_tolerance)
    upper_violations <- which(moments$upper > normalization_tolerance)
    if (length(lower_violations) + length(upper_violations) == 0L) {
      criterion <- as.numeric(
        0.5 * crossprod(beta, problem$riesz_hessian %*% beta) -
          crossprod(dvec, beta)
      )
      return(list(
        beta = beta,
        criterion = criterion,
        moments = moments,
        constraint_C = C,
        constraint_b = b,
        minimum_constraint_slack = min(as.vector(C %*% beta - b)),
        cutting_iterations = iteration,
        cut_log = if (length(cut_log) == 0L) {
          data.frame()
        } else {
          do.call(rbind, cut_log)
        },
        solver = fit
      ))
    }

    new_cuts <- c(
      lapply(lower_violations, function(index) {
        make_normalization_cut(problem, beta, "lower", index)
      }),
      lapply(upper_violations, function(index) {
        make_normalization_cut(problem, beta, "upper", index)
      })
    )
    new_C <- do.call(rbind, lapply(new_cuts, `[[`, "C"))
    new_b <- vapply(new_cuts, `[[`, numeric(1), "b")
    cut_log[[length(cut_log) + 1L]] <- data.frame(
      iteration = iteration,
      side = vapply(new_cuts, `[[`, character(1), "side"),
      instrument = vapply(new_cuts, `[[`, character(1), "instrument"),
      violation = c(
        moments$lower[lower_violations],
        moments$upper[upper_violations]
      ),
      stringsAsFactors = FALSE
    )
    old_count <- nrow(deduplicate_linear_constraints(C, b)$C)
    augmented <- deduplicate_linear_constraints(rbind(C, new_C), c(b, new_b))
    C <- augmented$C
    b <- augmented$b
    if (nrow(C) == old_count) {
      stop(
        "The normalization cutting-plane algorithm stalled with maximum ",
        sprintf("violation %.6g.", moments$maximum),
        call. = FALSE
      )
    }
  }
  stop(
    "The normalization cutting-plane algorithm exceeded its iteration limit.",
    call. = FALSE
  )
}

prepare_compatibility_arm_problem <- function(
    data,
    nuisance,
    t_value,
    Lambda0,
    Lambda1,
    x_cols = "X",
    gamma_rank = 25L,
    gamma_bandwidth_multiplier = 4,
    moment_resolutions = c(2L, 4L),
    ridge = 1e-8,
    pointwise_slack = 1e-9,
    seed = 1L
) {
  check_formal_compatibility_dependencies()
  validate_compatibility_nuisance(nuisance, nrow(data))
  if (!(t_value %in% c(0L, 1L))) {
    stop("t_value must equal zero or one.")
  }
  if (
    length(Lambda0) != 1L || length(Lambda1) != 1L ||
      !is.finite(Lambda0) || !is.finite(Lambda1) ||
      Lambda0 < 1 || Lambda1 < 1
  ) {
    stop("Lambda0 and Lambda1 must be finite scalars of at least one.")
  }
  arm_index <- which(data$T_obs == t_value)
  arm_data <- data[arm_index, , drop = FALSE]
  nuisance_arm <- list(
    e = nuisance$e[arm_index],
    pi0 = nuisance$pi0[arm_index],
    pi1 = nuisance$pi1[arm_index]
  )
  baseline <- arm_data$Z == 0
  if (sum(baseline) < 20L || sum(!baseline) < 20L) {
    stop("Each treatment arm needs at least 20 observations at both Z levels.")
  }

  sieve <- make_compatibility_gamma_sieve(
    arm_data = arm_data,
    x_cols = x_cols,
    rank = gamma_rank,
    bandwidth_multiplier = gamma_bandwidth_multiplier,
    seed = seed
  )
  Phi <- sieve$features
  W0 <- as.numeric(arm_data$Z == 0)
  w1 <- as.numeric(arm_data$Z == 1) * compute_riesz_weight(
    nuisance = nuisance_arm,
    t_value = t_value
  )
  n_t <- nrow(arm_data)
  penalty <- diag(ncol(Phi))
  penalty[1L, 1L] <- 0
  hessian <- crossprod(Phi, Phi * W0) / n_t + ridge * penalty
  # quadprog requires strict positive definiteness.  The tiny intercept ridge
  # is numerical only and is reported separately from the sieve penalty.
  minimum_eigenvalue <- min(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values)
  numerical_intercept_ridge <- max(0, 1e-10 - minimum_eigenvalue)
  if (numerical_intercept_ridge > 0) {
    diag(hessian) <- diag(hessian) + numerical_intercept_ridge
  }
  linear <- as.vector(crossprod(Phi, w1) / n_t)

  interval <- compute_compatible_gamma_interval(
    nuisance = nuisance_arm,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1
  )
  if (
    any(!is.finite(interval$lower)) || any(!is.finite(interval$upper)) ||
      any(interval$lower > interval$upper)
  ) {
    stop("The nuisance-based compatible gamma interval is invalid.")
  }
  pointwise_C <- rbind(Phi, -Phi)
  pointwise_b <- c(
    interval$lower - pointwise_slack,
    -interval$upper - pointwise_slack
  )

  instruments <- make_multiscale_box_instruments(
    data = arm_data[baseline, , drop = FALSE],
    x_cols = x_cols,
    resolutions = moment_resolutions
  )
  problem <- list(
    t_value = as.integer(t_value),
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    arm_index = arm_index,
    arm_data = arm_data,
    nuisance_arm = nuisance_arm,
    baseline = baseline,
    features = Phi,
    features_baseline = Phi[baseline, , drop = FALSE],
    sieve = sieve,
    moment_instruments = instruments$matrix,
    moment_instrument_metadata = instruments$metadata,
    envelope_affine = envelope_affine_coefficients(
      nuisance = list(
        pi0 = nuisance_arm$pi0[baseline],
        pi1 = nuisance_arm$pi1[baseline]
      ),
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1
    ),
    gamma_interval = interval,
    pointwise_C = pointwise_C,
    pointwise_b = pointwise_b,
    riesz_hessian = hessian,
    riesz_linear = linear,
    W0 = W0,
    w1 = w1,
    ridge = ridge,
    numerical_intercept_ridge = numerical_intercept_ridge,
    n_arm = n_t,
    n_baseline = sum(baseline),
    seed = seed
  )
  class(problem) <- "dual_msm_compatibility_problem"
  problem
}

fit_compatibility_arm <- function(
    problem,
    normalization_tolerance = 1e-8,
    maximum_cutting_iterations = 60L
) {
  H <- problem$riesz_hessian
  r <- problem$riesz_linear
  beta_unrestricted <- tryCatch(
    solve(H, r),
    error = function(error) qr.solve(H, r, tol = 1e-10)
  )
  criterion <- function(beta, dvec = r) {
    as.numeric(
      0.5 * crossprod(beta, H %*% beta) - crossprod(dvec, beta)
    )
  }
  constrained <- solve_compatibility_qp(
    problem = problem,
    dvec = r,
    normalization_tolerance = normalization_tolerance,
    maximum_cutting_iterations = maximum_cutting_iterations
  )
  criterion_unrestricted <- criterion(beta_unrestricted)
  gap <- max(0, constrained$criterion - criterion_unrestricted)
  statistic <- 2 * problem$n_arm * gap
  unrestricted_moments <- evaluate_sieve_normalization_moments(
    problem, beta_unrestricted
  )
  gamma_unrestricted <- as.vector(problem$features %*% beta_unrestricted)
  gamma_constrained <- as.vector(problem$features %*% constrained$beta)
  interval <- problem$gamma_interval
  pointwise_excess <- pmax(interval$lower - gamma_unrestricted, 0) +
    pmax(gamma_unrestricted - interval$upper, 0)

  list(
    problem = problem,
    beta_unrestricted = as.vector(beta_unrestricted),
    beta_constrained = constrained$beta,
    gamma_unrestricted = gamma_unrestricted,
    gamma_constrained = gamma_constrained,
    criterion_unrestricted = criterion_unrestricted,
    criterion_constrained = constrained$criterion,
    criterion_gap = gap,
    statistic = statistic,
    unrestricted_moments = unrestricted_moments,
    constrained_moments = constrained$moments,
    unrestricted_pointwise_violation_rate = mean(pointwise_excess > 1e-8),
    unrestricted_pointwise_maximum_excess = max(pointwise_excess),
    minimum_constraint_slack = constrained$minimum_constraint_slack,
    cutting_iterations = constrained$cutting_iterations,
    cut_log = constrained$cut_log,
    constraint_C = constrained$constraint_C,
    constraint_b = constrained$constraint_b
  )
}

bootstrap_compatibility_arm <- function(
    fit,
    multipliers,
    normalization_tolerance = 1e-8,
    maximum_cutting_iterations = 60L
) {
  problem <- fit$problem
  if (length(multipliers) != problem$n_arm) {
    stop("multipliers must have one value per treatment-arm observation.")
  }
  Phi <- problem$features
  gamma_constrained <- fit$gamma_constrained
  score <- Phi * (
    problem$W0 * gamma_constrained - problem$w1
  )
  score <- sweep(score, 2L, colMeans(score), "-")
  perturbation <- as.vector(crossprod(score, multipliers) / problem$n_arm)
  dvec <- as.vector(
    problem$riesz_hessian %*% fit$beta_constrained - perturbation
  )
  beta_unrestricted <- tryCatch(
    solve(problem$riesz_hessian, dvec),
    error = function(error) {
      qr.solve(problem$riesz_hessian, dvec, tol = 1e-10)
    }
  )
  constrained <- solve_compatibility_qp(
    problem = problem,
    dvec = dvec,
    initial_C = fit$constraint_C,
    initial_b = fit$constraint_b,
    normalization_tolerance = normalization_tolerance,
    maximum_cutting_iterations = maximum_cutting_iterations
  )
  criterion <- function(beta) {
    as.numeric(
      0.5 * crossprod(beta, problem$riesz_hessian %*% beta) -
        crossprod(dvec, beta)
    )
  }
  gap <- max(0, constrained$criterion - criterion(beta_unrestricted))
  2 * problem$n_arm * gap
}

test_dual_msm_compatibility <- function(
    data,
    nuisance,
    Lambda0,
    Lambda1,
    x_cols = "X",
    gamma_rank = 25L,
    gamma_bandwidth_multiplier = 4,
    moment_resolutions = c(2L, 4L),
    ridge = 1e-8,
    bootstrap_replications = 199L,
    confidence_level = 0.95,
    normalization_tolerance = 1e-8,
    seed = 1L
) {
  check_formal_compatibility_dependencies()
  if (
    length(bootstrap_replications) != 1L ||
      bootstrap_replications < 19L ||
      abs(bootstrap_replications - round(bootstrap_replications)) > 0
  ) {
    stop("bootstrap_replications must be one integer of at least 19.")
  }
  if (
    length(confidence_level) != 1L ||
      !is.finite(confidence_level) ||
      confidence_level <= 0 || confidence_level >= 1
  ) {
    stop("confidence_level must lie strictly between zero and one.")
  }

  arm_fits <- lapply(c(0L, 1L), function(t_value) {
    problem <- prepare_compatibility_arm_problem(
      data = data,
      nuisance = nuisance,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      x_cols = x_cols,
      gamma_rank = gamma_rank,
      gamma_bandwidth_multiplier = gamma_bandwidth_multiplier,
      moment_resolutions = moment_resolutions,
      ridge = ridge,
      seed = seed + 10000L * t_value
    )
    fit_compatibility_arm(
      problem,
      normalization_tolerance = normalization_tolerance
    )
  })
  names(arm_fits) <- c("0", "1")
  observed_joint <- sum(vapply(arm_fits, `[[`, numeric(1), "statistic"))

  set.seed(seed + 900000L)
  bootstrap_arm <- matrix(
    NA_real_, nrow = bootstrap_replications, ncol = 2L,
    dimnames = list(NULL, c("t0", "t1"))
  )
  for (replication in seq_len(bootstrap_replications)) {
    for (arm in seq_along(arm_fits)) {
      bootstrap_arm[replication, arm] <- bootstrap_compatibility_arm(
        fit = arm_fits[[arm]],
        multipliers = stats::rnorm(arm_fits[[arm]]$problem$n_arm),
        normalization_tolerance = normalization_tolerance
      )
    }
  }
  bootstrap_joint <- rowSums(bootstrap_arm)
  alpha <- 1 - confidence_level
  joint_p_value <- (
    1 + sum(bootstrap_joint >= observed_joint - 1e-12)
  ) / (bootstrap_replications + 1)
  arm_p_values <- vapply(seq_along(arm_fits), function(index) {
    observed <- arm_fits[[index]]$statistic
    (1 + sum(bootstrap_arm[, index] >= observed - 1e-12)) /
      (bootstrap_replications + 1)
  }, numeric(1))

  arm_summary <- do.call(rbind, lapply(seq_along(arm_fits), function(index) {
    fit <- arm_fits[[index]]
    data.frame(
      treatment = index - 1L,
      n_arm = fit$problem$n_arm,
      n_baseline = fit$problem$n_baseline,
      statistic = fit$statistic,
      bootstrap_p_value = arm_p_values[index],
      criterion_gap = fit$criterion_gap,
      unrestricted_pointwise_violation_rate =
        fit$unrestricted_pointwise_violation_rate,
      unrestricted_pointwise_maximum_excess =
        fit$unrestricted_pointwise_maximum_excess,
      unrestricted_maximum_normalization_moment =
        fit$unrestricted_moments$maximum,
      constrained_maximum_normalization_moment =
        fit$constrained_moments$maximum,
      cutting_iterations = fit$cutting_iterations,
      gamma_rank_requested = fit$problem$sieve$rank_requested,
      gamma_effective_rank = fit$problem$sieve$effective_rank,
      base_bandwidth = fit$problem$sieve$base_bandwidth,
      bandwidth_multiplier = fit$problem$sieve$bandwidth_multiplier,
      selected_bandwidth = fit$problem$sieve$selected_bandwidth,
      stringsAsFactors = FALSE
    )
  }))
  joint_summary <- data.frame(
    method = "fixed_sieve_null_imposed_multiplier_prototype",
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    statistic = observed_joint,
    bootstrap_critical_value = as.numeric(stats::quantile(
      bootstrap_joint,
      probs = confidence_level,
      type = 8,
      names = FALSE
    )),
    bootstrap_p_value = joint_p_value,
    confidence_level = confidence_level,
    reject = joint_p_value <= alpha + 10 * .Machine$double.eps,
    bootstrap_replications = bootstrap_replications,
    gamma_rank = gamma_rank,
    bandwidth_multiplier = gamma_bandwidth_multiplier,
    moment_sieve = paste(moment_resolutions, collapse = ","),
    formal_status = paste(
      "Fixed-sieve pointwise prototype; primitive nuisances and sieve are",
      "held fixed in the multiplier bootstrap."
    ),
    stringsAsFactors = FALSE
  )

  list(
    joint = joint_summary,
    arms = arm_summary,
    arm_fits = arm_fits,
    bootstrap_joint = bootstrap_joint,
    bootstrap_arm = bootstrap_arm
  )
}

invert_symmetric_compatibility_tests <- function(
    data,
    nuisance,
    lambda_grid,
    seed = 1L,
    ...
) {
  lambda_grid <- sort(unique(lambda_grid))
  if (length(lambda_grid) == 0L || any(lambda_grid < 1)) {
    stop("lambda_grid must contain sensitivity values of at least one.")
  }
  fits <- lapply(seq_along(lambda_grid), function(index) {
    Lambda <- lambda_grid[index]
    test_dual_msm_compatibility(
      data = data,
      nuisance = nuisance,
      Lambda0 = Lambda,
      Lambda1 = Lambda,
      # Keep the feature sieve, landmarks, and multiplier draws common across
      # Lambda.  Otherwise numerical changes in the approximation could mask
      # the nesting of the population compatibility models.
      seed = seed,
      ...
    )
  })
  table <- do.call(rbind, lapply(seq_along(fits), function(index) {
    row <- fits[[index]]$joint
    row$Lambda <- lambda_grid[index]
    row
  }))
  accepted <- table$Lambda[!table$reject]
  list(
    tests = fits,
    summary = table,
    smallest_grid_value_not_rejected = if (length(accepted) == 0L) {
      Inf
    } else {
      min(accepted)
    },
    interpretation = paste(
      "This is an inverted confidence set for the fixed-sieve population",
      "compatibility frontier, not an estimate of the latent confounding",
      "magnitude."
    )
  )
}

make_default_compatibility_omnibus_specifications <- function(data) {
  if (!"T_obs" %in% names(data)) {
    stop("data must contain T_obs.")
  }
  arm_counts <- table(factor(data$T_obs, levels = c(0, 1)))
  if (any(arm_counts < 20L)) {
    stop("Both treatment arms need at least 20 observations.")
  }
  reference_size <- min(as.numeric(arm_counts))
  ranks <- unique(pmax(
    5L,
    as.integer(ceiling(c(2, 4) * reference_size^(1 / 4)))
  ))
  grid <- expand.grid(
    gamma_rank = ranks,
    gamma_bandwidth_multiplier = c(1, 2, 4),
    ridge_constant = c(0.1, 1),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$ridge <- grid$ridge_constant / sqrt(reference_size)
  grid$specification <- sprintf(
    "rank%d_bw%s_c%s",
    grid$gamma_rank,
    format(grid$gamma_bandwidth_multiplier, trim = TRUE, scientific = FALSE),
    format(grid$ridge_constant, trim = TRUE, scientific = FALSE)
  )
  grid[, c(
    "specification", "gamma_rank", "gamma_bandwidth_multiplier",
    "ridge_constant", "ridge"
  )]
}

make_focused_compatibility_omnibus_specifications <- function(data) {
  if (!"T_obs" %in% names(data)) {
    stop("data must contain T_obs.")
  }
  arm_counts <- table(factor(data$T_obs, levels = c(0, 1)))
  if (any(arm_counts < 20L)) {
    stop("Both treatment arms need at least 20 observations.")
  }
  reference_size <- min(as.numeric(arm_counts))
  ranks <- unique(pmax(
    5L,
    as.integer(ceiling(c(2, 4) * reference_size^(1 / 4)))
  ))
  if (length(ranks) == 1L) {
    ranks <- c(ranks, ranks + 1L)
  }
  low_rank <- min(ranks)
  high_rank <- max(ranks)

  # This collection was fixed from the oracle approximation audit, without
  # reference to bootstrap rejection rates.  The low-rank fits emphasize
  # global gamma accuracy.  The progressively less regularized high-rank fits
  # add sensitivity to incompatibilities concentrated in thin tail regions.
  specifications <- data.frame(
    specification = c(
      "stable_low_bw4_c1",
      "stable_low_bw2_c1",
      "stable_low_bw1_c1",
      "tail_high_bw4_c1",
      "tail_high_bw4_c0p3",
      "tail_high_bw4_c0p1"
    ),
    gamma_rank = c(
      low_rank, low_rank, low_rank,
      high_rank, high_rank, high_rank
    ),
    gamma_bandwidth_multiplier = c(4, 2, 1, 4, 4, 4),
    ridge_constant = c(1, 1, 1, 1, 0.3, 0.1),
    stringsAsFactors = FALSE
  )
  specifications$ridge <-
    specifications$ridge_constant / sqrt(reference_size)
  specifications
}

validate_compatibility_omnibus_specifications <- function(specifications) {
  required <- c(
    "specification", "gamma_rank", "gamma_bandwidth_multiplier"
  )
  missing <- setdiff(required, names(specifications))
  if (length(missing) > 0L) {
    stop("Omnibus specifications omit: ", paste(missing, collapse = ", "))
  }
  if (
    nrow(specifications) < 1L ||
      anyNA(specifications[, required, drop = FALSE]) ||
      anyDuplicated(as.character(specifications$specification)) ||
      any(specifications$gamma_rank < 1) ||
      any(abs(specifications$gamma_rank - round(specifications$gamma_rank)) > 0) ||
      any(!is.finite(specifications$gamma_bandwidth_multiplier)) ||
      any(specifications$gamma_bandwidth_multiplier <= 0)
  ) {
    stop("The omnibus specification table is invalid.")
  }
  specifications$specification <- as.character(specifications$specification)
  specifications$gamma_rank <- as.integer(specifications$gamma_rank)
  specifications$gamma_bandwidth_multiplier <- as.numeric(
    specifications$gamma_bandwidth_multiplier
  )
  if (!"ridge" %in% names(specifications)) {
    specifications$ridge <- NA_real_
  }
  if (
    any(!is.na(specifications$ridge) &
      (!is.finite(specifications$ridge) | specifications$ridge <= 0))
  ) {
    stop("Every specification-specific ridge must be positive and finite.")
  }
  specifications$ridge <- as.numeric(specifications$ridge)
  specifications
}

empirical_upper_tail_probability <- function(reference, value) {
  (1 + sum(reference >= value - 1e-12)) / (length(reference) + 1)
}

bootstrap_self_upper_tail_probabilities <- function(reference) {
  vapply(reference, function(value) {
    # Treat each draw as pseudo-observed and use the other B-1 draws as its
    # reference distribution.  Including its own tie in the numerator gives
    # the algebraically equivalent count/B expression.
    sum(reference >= value - 1e-12) / length(reference)
  }, numeric(1))
}

test_dual_msm_compatibility_omnibus <- function(
    data,
    nuisance,
    Lambda0,
    Lambda1,
    x_cols = "X",
    specifications = NULL,
    moment_resolutions = c(2L, 4L, 8L, 16L),
    ridge = 1e-6,
    bootstrap_replications = 199L,
    confidence_level = 0.95,
    normalization_tolerance = 1e-8,
    seed = 1L
) {
  if (is.null(specifications)) {
    specifications <- make_default_compatibility_omnibus_specifications(data)
  }
  specifications <- validate_compatibility_omnibus_specifications(
    specifications
  )
  if (length(specifications$specification) < 2L) {
    warning(
      "The omnibus test contains only one specification; use the ordinary ",
      "fixed-sieve test unless this is intentional."
    )
  }

  number_of_specs <- nrow(specifications)
  bootstrap_by_specification <- matrix(
    NA_real_,
    nrow = bootstrap_replications,
    ncol = number_of_specs,
    dimnames = list(NULL, specifications$specification)
  )
  specification_rows <- vector("list", number_of_specs)
  arm_rows <- vector("list", number_of_specs)
  for (index in seq_len(number_of_specs)) {
    specification <- specifications[index, , drop = FALSE]
    fit <- test_dual_msm_compatibility(
      data = data,
      nuisance = nuisance,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      x_cols = x_cols,
      gamma_rank = specification$gamma_rank,
      gamma_bandwidth_multiplier =
        specification$gamma_bandwidth_multiplier,
      moment_resolutions = moment_resolutions,
      ridge = if (is.finite(specification$ridge)) {
        specification$ridge
      } else {
        ridge
      },
      bootstrap_replications = bootstrap_replications,
      confidence_level = confidence_level,
      normalization_tolerance = normalization_tolerance,
      # Keeping this seed common makes the Gaussian multipliers identical
      # across specifications.  The joint bootstrap therefore calibrates the
      # search over bandwidth and rank rather than treating it as fixed.
      seed = seed
    )
    bootstrap_by_specification[, index] <- fit$bootstrap_joint
    row <- fit$joint
    row$specification <- specification$specification
    row$gamma_rank <- specification$gamma_rank
    row$bandwidth_multiplier <-
      specification$gamma_bandwidth_multiplier
    row$ridge <- if (is.finite(specification$ridge)) {
      specification$ridge
    } else {
      ridge
    }
    row$ridge_constant <- if (
      "ridge_constant" %in% names(specification) &&
        is.finite(specification$ridge_constant)
    ) {
      specification$ridge_constant
    } else {
      NA_real_
    }
    specification_rows[[index]] <- row
    arms <- fit$arms
    arms$specification <- specification$specification
    arms$ridge <- if (is.finite(specification$ridge)) {
      specification$ridge
    } else {
      ridge
    }
    arms$ridge_constant <- row$ridge_constant
    arm_rows[[index]] <- arms
    rm(fit)
  }
  specification_summary <- do.call(rbind, specification_rows)
  rownames(specification_summary) <- NULL
  arms <- do.call(rbind, arm_rows)
  rownames(arms) <- NULL

  observed_statistics <- specification_summary$statistic
  observed_marginal_p <- vapply(seq_len(number_of_specs), function(index) {
    empirical_upper_tail_probability(
      bootstrap_by_specification[, index],
      observed_statistics[index]
    )
  }, numeric(1))
  bootstrap_marginal_p <- apply(
    bootstrap_by_specification,
    2L,
    bootstrap_self_upper_tail_probabilities
  )
  if (number_of_specs == 1L) {
    bootstrap_marginal_p <- matrix(
      bootstrap_marginal_p,
      ncol = 1L,
      dimnames = list(NULL, specifications$specification)
    )
  }
  observed_score_by_specification <- -log(observed_marginal_p)
  bootstrap_score_by_specification <- -log(bootstrap_marginal_p)
  observed_omnibus_score <- max(observed_score_by_specification)
  bootstrap_omnibus_score <- apply(
    bootstrap_score_by_specification,
    1L,
    max
  )
  omnibus_p_value <- (
    1 + sum(bootstrap_omnibus_score >= observed_omnibus_score - 1e-12)
  ) / (bootstrap_replications + 1)
  alpha <- 1 - confidence_level
  selected_index <- which.min(observed_marginal_p)
  specification_summary$marginal_p_value <- observed_marginal_p
  specification_summary$marginal_score <- observed_score_by_specification
  specification_summary$selected_by_omnibus <-
    seq_len(number_of_specs) == selected_index

  omnibus <- data.frame(
    method = "fixed_sieve_joint_min_p_omnibus_prototype",
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    omnibus_score = observed_omnibus_score,
    bootstrap_critical_value = as.numeric(stats::quantile(
      bootstrap_omnibus_score,
      probs = confidence_level,
      type = 8,
      names = FALSE
    )),
    bootstrap_p_value = omnibus_p_value,
    confidence_level = confidence_level,
    reject = omnibus_p_value <= alpha + 10 * .Machine$double.eps,
    number_of_specifications = number_of_specs,
    selected_specification = specifications$specification[selected_index],
    selected_marginal_p_value = observed_marginal_p[selected_index],
    moment_sieve = paste(moment_resolutions, collapse = ","),
    bootstrap_replications = bootstrap_replications,
    formal_status = paste(
      "Finite multiscale pointwise prototype; the common-multiplier min-p",
      "bootstrap accounts for the prespecified rank/bandwidth search, while",
      "primitive nuisances and the finite specification set remain fixed."
    ),
    stringsAsFactors = FALSE
  )

  list(
    omnibus = omnibus,
    specifications = specification_summary,
    arms = arms,
    bootstrap_omnibus_score = bootstrap_omnibus_score,
    bootstrap_by_specification = bootstrap_by_specification,
    bootstrap_marginal_p = bootstrap_marginal_p
  )
}

evaluate_gamma_normalization_with_instruments <- function(
    nuisance,
    gamma,
    t_value,
    Lambda0,
    Lambda1,
    instruments
) {
  affine <- envelope_affine_coefficients(
    nuisance = nuisance,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1
  )
  lower_affine <- affine$lower_slope * gamma + affine$shared_intercept
  upper_affine <- affine$upper_slope * gamma + affine$shared_intercept
  L <- pmax(affine$lower_constant, lower_affine)
  U <- pmin(affine$upper_constant, upper_affine)
  lower <- colMeans(instruments * (L - 1))
  upper <- colMeans(instruments * (1 - U))
  list(
    lower = lower,
    upper = upper,
    maximum = max(c(lower, upper)),
    L = L,
    U = U
  )
}

audit_compatibility_gamma_sieves <- function(
    data,
    nuisance,
    lambda_grid = c(1, 2, 3),
    x_cols = "X",
    specifications = NULL,
    moment_resolutions = c(2L, 4L, 8L, 16L),
    ridge = 1e-6,
    violation_tolerance = 1e-8,
    seed = 1L
) {
  required <- c("T_obs", "Z", "gamma0_oracle", "gamma1_oracle")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Oracle sieve audit data omit: ", paste(missing, collapse = ", "))
  }
  validate_compatibility_nuisance(nuisance, nrow(data))
  lambda_grid <- sort(unique(lambda_grid))
  if (length(lambda_grid) == 0L || any(lambda_grid < 1)) {
    stop("lambda_grid must contain values of at least one.")
  }
  if (is.null(specifications)) {
    default <- make_default_compatibility_omnibus_specifications(data)
    smooth_reference <- data.frame(
      specification = "rank15_bw4_reference",
      gamma_rank = 15L,
      gamma_bandwidth_multiplier = 4,
      ridge = ridge,
      stringsAsFactors = FALSE
    )
    shared_names <- union(names(default), names(smooth_reference))
    for (name in setdiff(shared_names, names(default))) default[[name]] <- NA
    for (name in setdiff(shared_names, names(smooth_reference))) {
      smooth_reference[[name]] <- NA
    }
    specifications <- rbind(
      default[, shared_names, drop = FALSE],
      smooth_reference[, shared_names, drop = FALSE]
    )
    specifications <- specifications[
      !duplicated(specifications$specification), , drop = FALSE
    ]
  }
  specifications <- validate_compatibility_omnibus_specifications(
    specifications
  )

  accuracy_rows <- list()
  compatibility_rows <- list()
  row_accuracy <- 0L
  row_compatibility <- 0L
  for (specification_index in seq_len(nrow(specifications))) {
    specification <- specifications[specification_index, , drop = FALSE]
    for (t_value in c(0L, 1L)) {
      problem <- prepare_compatibility_arm_problem(
        data = data,
        nuisance = nuisance,
        t_value = t_value,
        Lambda0 = lambda_grid[1L],
        Lambda1 = lambda_grid[1L],
        x_cols = x_cols,
        gamma_rank = specification$gamma_rank,
        gamma_bandwidth_multiplier =
          specification$gamma_bandwidth_multiplier,
        moment_resolutions = moment_resolutions,
        ridge = if (is.finite(specification$ridge)) {
          specification$ridge
        } else {
          ridge
        },
        seed = seed
      )
      beta <- tryCatch(
        solve(problem$riesz_hessian, problem$riesz_linear),
        error = function(error) {
          qr.solve(
            problem$riesz_hessian,
            problem$riesz_linear,
            tol = 1e-10
          )
        }
      )
      gamma_fitted <- as.vector(problem$features %*% beta)
      gamma_oracle <- if (t_value == 1L) {
        data$gamma1_oracle[problem$arm_index]
      } else {
        data$gamma0_oracle[problem$arm_index]
      }
      baseline <- problem$baseline
      summarize_accuracy <- function(index, region) {
        difference <- gamma_fitted[index] - gamma_oracle[index]
        correlation <- if (
          sum(index) >= 3L &&
            stats::sd(gamma_fitted[index]) > 0 &&
            stats::sd(gamma_oracle[index]) > 0
        ) {
          stats::cor(gamma_fitted[index], gamma_oracle[index])
        } else {
          NA_real_
        }
        data.frame(
          specification = specification$specification,
          treatment = t_value,
          region = region,
          observations = sum(index),
          bias = mean(difference),
          rmse = sqrt(mean(difference^2)),
          mae = mean(abs(difference)),
          correlation = correlation,
          mean_fitted = mean(gamma_fitted[index]),
          mean_oracle = mean(gamma_oracle[index]),
          gamma_rank = specification$gamma_rank,
          gamma_effective_rank = problem$sieve$effective_rank,
          bandwidth_multiplier =
            specification$gamma_bandwidth_multiplier,
          base_bandwidth = problem$sieve$base_bandwidth,
          selected_bandwidth = problem$sieve$selected_bandwidth,
          ridge = if (is.finite(specification$ridge)) {
            specification$ridge
          } else {
            ridge
          },
          ridge_constant = if (
            "ridge_constant" %in% names(specification) &&
              is.finite(specification$ridge_constant)
          ) {
            specification$ridge_constant
          } else {
            NA_real_
          },
          stringsAsFactors = FALSE
        )
      }
      for (region in c("arm", "baseline")) {
        row_accuracy <- row_accuracy + 1L
        accuracy_rows[[row_accuracy]] <- summarize_accuracy(
          if (region == "arm") rep(TRUE, length(baseline)) else baseline,
          region
        )
      }

      nuisance_arm <- problem$nuisance_arm
      nuisance_baseline <- list(
        pi0 = nuisance_arm$pi0[baseline],
        pi1 = nuisance_arm$pi1[baseline]
      )
      instruments <- problem$moment_instruments
      for (Lambda in lambda_grid) {
        interval <- compute_compatible_gamma_interval(
          nuisance = nuisance_arm,
          t_value = t_value,
          Lambda0 = Lambda,
          Lambda1 = Lambda
        )
        pointwise_excess <- function(gamma) {
          pmax(interval$lower - gamma, 0) +
            pmax(gamma - interval$upper, 0)
        }
        oracle_excess <- pointwise_excess(gamma_oracle)
        fitted_excess <- pointwise_excess(gamma_fitted)
        oracle_violation <- oracle_excess > violation_tolerance
        fitted_violation <- fitted_excess > violation_tolerance
        oracle_normalization <- evaluate_gamma_normalization_with_instruments(
          nuisance = nuisance_baseline,
          gamma = gamma_oracle[baseline],
          t_value = t_value,
          Lambda0 = Lambda,
          Lambda1 = Lambda,
          instruments = instruments
        )
        fitted_normalization <- evaluate_gamma_normalization_with_instruments(
          nuisance = nuisance_baseline,
          gamma = gamma_fitted[baseline],
          t_value = t_value,
          Lambda0 = Lambda,
          Lambda1 = Lambda,
          instruments = instruments
        )
        row_compatibility <- row_compatibility + 1L
        compatibility_rows[[row_compatibility]] <- data.frame(
          specification = specification$specification,
          treatment = t_value,
          Lambda = Lambda,
          oracle_pointwise_violation_rate = mean(oracle_violation),
          fitted_pointwise_violation_rate = mean(fitted_violation),
          fitted_detection_rate_on_oracle_violations = if (
            any(oracle_violation)
          ) {
            mean(fitted_violation[oracle_violation])
          } else {
            NA_real_
          },
          fitted_false_positive_rate_off_oracle_violations = if (
            any(!oracle_violation)
          ) {
            mean(fitted_violation[!oracle_violation])
          } else {
            NA_real_
          },
          oracle_mean_pointwise_excess = mean(oracle_excess),
          fitted_mean_pointwise_excess = mean(fitted_excess),
          oracle_maximum_pointwise_excess = max(oracle_excess),
          fitted_maximum_pointwise_excess = max(fitted_excess),
          oracle_maximum_normalization_moment =
            oracle_normalization$maximum,
          fitted_maximum_normalization_moment =
            fitted_normalization$maximum,
          gamma_rank = specification$gamma_rank,
          gamma_effective_rank = problem$sieve$effective_rank,
          bandwidth_multiplier =
            specification$gamma_bandwidth_multiplier,
          selected_bandwidth = problem$sieve$selected_bandwidth,
          ridge = if (is.finite(specification$ridge)) {
            specification$ridge
          } else {
            ridge
          },
          ridge_constant = if (
            "ridge_constant" %in% names(specification) &&
              is.finite(specification$ridge_constant)
          ) {
            specification$ridge_constant
          } else {
            NA_real_
          },
          stringsAsFactors = FALSE
        )
      }
    }
  }
  list(
    specifications = specifications,
    accuracy = do.call(rbind, accuracy_rows),
    compatibility = do.call(rbind, compatibility_rows),
    interpretation = paste(
      "This audit selects candidate sieves by approximation of oracle gamma",
      "and recovery of oracle compatibility violations, not by bootstrap",
      "rejection rates."
    )
  )
}
