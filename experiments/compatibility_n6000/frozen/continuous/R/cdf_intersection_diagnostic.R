# Model-assisted CDF-intersection diagnostic for the Dual MSM.
#
# This file intentionally does not modify the sharp-bound estimator.  It
# implements a descriptive necessary-condition diagnostic that avoids direct
# estimation of gamma_t.  For each treatment arm, outcome threshold, and
# instrument stratum, it applies the ordinary MSM to the bounded transform
# I(Y <= y).  IV exchangeability and distributional exclusion require the two
# stratum-specific intervals for P(Y^t <= y | X) to intersect.
#
# The resulting scores are heuristic.  They are not p-values, formal tests, or
# confidence bounds for the exact compatibility frontier.

check_cdf_intersection_dependencies <- function() {
  required <- c(
    "make_spline_design",
    "make_tz_stratified_folds",
    "fit_nuisance_models",
    "gauss_legendre_rule",
    "dtnormtc"
  )
  missing <- required[
    !vapply(required, exists, logical(1), mode = "function")
  ]
  if (length(missing) > 0L) {
    stop(
      "Source R/sharp_bounds_plugin.R first. Missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

cdf_clip_probability <- function(value, epsilon = 1e-6) {
  output <- value
  output[] <- pmin(1 - epsilon, pmax(epsilon, value))
  output
}

truncated_normal_cdf <- function(
    value,
    mean = 0,
    standard_deviation = 1
) {
  standardized <- (value - mean) / standard_deviation
  denominator <- stats::pnorm(2) - stats::pnorm(-2)
  output <- (
    stats::pnorm(pmin(2, pmax(-2, standardized))) - stats::pnorm(-2)
  ) / denominator
  output[standardized <= -2] <- 0
  output[standardized >= 2] <- 1
  output
}

make_cdf_threshold_grid <- function(
    data,
    treatment_col = "T_obs",
    outcome_col = "Y",
    probabilities = seq(0.05, 0.95, by = 0.05)
) {
  required <- c(treatment_col, outcome_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Threshold-grid data omit: ", paste(missing, collapse = ", "))
  }
  if (
    length(probabilities) < 3L ||
      any(!is.finite(probabilities)) ||
      any(probabilities <= 0 | probabilities >= 1)
  ) {
    stop("probabilities must contain at least three values in (0,1).")
  }
  treatment <- data[[treatment_col]]
  outcome <- data[[outcome_col]]
  output <- lapply(c(0L, 1L), function(t_value) {
    values <- outcome[treatment == t_value]
    if (length(values) < 20L || any(!is.finite(values))) {
      stop("Each treatment arm needs at least 20 finite outcomes.")
    }
    unique(as.numeric(stats::quantile(
      values,
      probs = sort(unique(probabilities)),
      type = 8,
      names = FALSE
    )))
  })
  names(output) <- c("0", "1")
  if (any(vapply(output, length, integer(1)) < 3L)) {
    stop("The outcome has too few distinct values for a CDF threshold grid.")
  }
  output
}

msm_binary_event_interval <- function(
    observed_probability,
    arm_probability,
    Lambda
) {
  if (
    length(Lambda) != 1L || !is.finite(Lambda) || Lambda < 1
  ) {
    stop("Lambda must be one finite value of at least one.")
  }
  if (
    any(!is.finite(observed_probability)) ||
      any(observed_probability < 0 | observed_probability > 1) ||
      any(!is.finite(arm_probability)) ||
      any(arm_probability <= 0 | arm_probability >= 1)
  ) {
    stop("Event and treatment-arm probabilities are invalid.")
  }
  if (
    length(arm_probability) != 1L &&
      length(arm_probability) != length(observed_probability)
  ) {
    stop("arm_probability must be scalar or conformable with the event CDF.")
  }

  r <- observed_probability
  p <- arm_probability
  missing_lower <- pmax(r / Lambda, 1 - Lambda * (1 - r))
  missing_upper <- pmin(Lambda * r, 1 - (1 - r) / Lambda)
  lower <- p * r + (1 - p) * missing_lower
  upper <- p * r + (1 - p) * missing_upper
  list(
    lower = pmax(0, pmin(1, lower)),
    upper = pmax(0, pmin(1, upper))
  )
}

rearrange_cdf_rows <- function(probability_matrix) {
  probability_matrix <- as.matrix(probability_matrix)
  if (ncol(probability_matrix) <= 1L) return(probability_matrix)
  output <- t(apply(probability_matrix, 1L, sort))
  dimnames(output) <- dimnames(probability_matrix)
  output
}

fit_penalized_spline_distribution <- function(
    training_data,
    thresholds,
    basis_data = training_data,
    x_cols = "X",
    spline_df = 4L,
    ridge = 1e-4,
    probability_epsilon = 1e-5,
    maximum_iterations = 500L
) {
  check_cdf_intersection_dependencies()
  if (nrow(training_data) < spline_df + 10L) {
    stop("Too few observations for a spline distribution regression.")
  }
  design <- make_spline_design(
    data = basis_data,
    x_cols = x_cols,
    spline_df = spline_df
  )
  X <- design$predict_matrix(training_data)
  y <- training_data$Y
  penalty <- rep(1, ncol(X))
  penalty[1L] <- 0
  coefficients <- matrix(
    NA_real_,
    nrow = ncol(X),
    ncol = length(thresholds)
  )
  convergence <- integer(length(thresholds))
  constant_response <- logical(length(thresholds))
  previous <- rep(0, ncol(X))

  for (index in seq_along(thresholds)) {
    response <- as.numeric(y <= thresholds[index])
    event_mean <- (sum(response) + 0.5) / (length(response) + 1)
    initial <- previous
    initial[1L] <- stats::qlogis(event_mean)
    if (length(unique(response)) == 1L) {
      constant_response[index] <- TRUE
      coefficients[, index] <- 0
      coefficients[1L, index] <- stats::qlogis(event_mean)
      previous <- coefficients[, index]
      next
    }
    objective <- function(beta) {
      eta <- as.vector(X %*% beta)
      logistic_loss <- pmax(eta, 0) + log1p(exp(-abs(eta))) - response * eta
      mean(logistic_loss) + 0.5 * ridge * sum(penalty * beta^2)
    }
    gradient <- function(beta) {
      eta <- as.vector(X %*% beta)
      as.vector(crossprod(X, stats::plogis(eta) - response)) / nrow(X) +
        ridge * penalty * beta
    }
    fit <- stats::optim(
      par = initial,
      fn = objective,
      gr = gradient,
      method = "BFGS",
      control = list(maxit = maximum_iterations, reltol = 1e-8)
    )
    if (any(!is.finite(fit$par))) {
      stop("The spline distribution regression returned nonfinite values.")
    }
    coefficients[, index] <- fit$par
    convergence[index] <- fit$convergence
    previous <- fit$par
  }

  predict_cdf <- function(newdata, rearrange = TRUE) {
    prediction <- stats::plogis(
      design$predict_matrix(newdata) %*% coefficients
    )
    prediction <- cdf_clip_probability(
      prediction,
      epsilon = probability_epsilon
    )
    colnames(prediction) <- format(thresholds, digits = 12)
    if (rearrange) rearrange_cdf_rows(prediction) else prediction
  }
  training_raw <- predict_cdf(training_data, rearrange = FALSE)
  training_rearranged <- rearrange_cdf_rows(training_raw)

  list(
    predict = predict_cdf,
    coefficients = coefficients,
    thresholds = thresholds,
    diagnostics = data.frame(
      observations = nrow(training_data),
      basis_observations = nrow(basis_data),
      thresholds = length(thresholds),
      spline_df = spline_df,
      ridge = ridge,
      nonconvergence_rate = mean(convergence != 0L),
      constant_response_rate = mean(constant_response),
      rearrangement_change_rate = mean(
        abs(training_raw - training_rearranged) > 1e-10
      ),
      stringsAsFactors = FALSE
    )
  )
}

crossfit_spline_cdf_layer <- function(
    data,
    thresholds_by_t,
    x_cols = "X",
    folds = 5L,
    spline_df = 4L,
    distribution_ridge = 1e-4,
    nuisance_learner = "glm_spline",
    positivity_epsilon = 0.05,
    seed = 1L
) {
  check_cdf_intersection_dependencies()
  fold_id <- make_tz_stratified_folds(data, K = folds, seed = seed)
  n <- nrow(data)
  nuisance <- list(
    e = rep(NA_real_, n),
    pi0 = rep(NA_real_, n),
    pi1 = rep(NA_real_, n)
  )
  cdf <- list()
  for (t_value in c(0L, 1L)) {
    for (z_value in c(0L, 1L)) {
      key <- paste0("t", t_value, "_z", z_value)
      cdf[[key]] <- matrix(
        NA_real_,
        nrow = n,
        ncol = length(thresholds_by_t[[as.character(t_value)]])
      )
    }
  }
  diagnostics <- list()
  diagnostic_index <- 0L
  for (fold in seq_len(folds)) {
    training <- fold_id != fold
    validation <- !training
    nuisance_fit <- fit_nuisance_models(
      train = data[training, , drop = FALSE],
      x_cols = x_cols,
      positivity_epsilon = positivity_epsilon,
      learner = nuisance_learner,
      treatment_spline_df = spline_df
    )
    nuisance_prediction <- nuisance_fit$predict(
      data[validation, , drop = FALSE]
    )
    for (name in names(nuisance)) {
      nuisance[[name]][validation] <- nuisance_prediction[[name]]
    }

    for (t_value in c(0L, 1L)) {
      thresholds <- thresholds_by_t[[as.character(t_value)]]
      for (z_value in c(0L, 1L)) {
        cell <- training & data$T_obs == t_value & data$Z == z_value
        fit <- fit_penalized_spline_distribution(
          training_data = data[cell, , drop = FALSE],
          thresholds = thresholds,
          basis_data = data[training, , drop = FALSE],
          x_cols = x_cols,
          spline_df = spline_df,
          ridge = distribution_ridge
        )
        key <- paste0("t", t_value, "_z", z_value)
        cdf[[key]][validation, ] <- fit$predict(
          data[validation, , drop = FALSE]
        )
        diagnostic_index <- diagnostic_index + 1L
        local <- fit$diagnostics
        local$fold <- fold
        local$treatment <- t_value
        local$instrument <- z_value
        diagnostics[[diagnostic_index]] <- local
      }
    }
  }
  if (
    any(vapply(nuisance, function(value) any(!is.finite(value)), logical(1))) ||
      any(vapply(cdf, function(value) any(!is.finite(value)), logical(1)))
  ) {
    stop("Cross-fitted CDF or nuisance predictions are incomplete.")
  }
  list(
    layer = "estimated_spline",
    cdf = cdf,
    nuisance = nuisance,
    fold_id = fold_id,
    diagnostics = do.call(rbind, diagnostics),
    specification = list(
      folds = folds,
      spline_df = spline_df,
      distribution_ridge = distribution_ridge,
      nuisance_learner = nuisance_learner,
      rearrangement = "increasing row rearrangement"
    )
  )
}

compute_dgp_oracle_cdf_layer <- function(
    data,
    thresholds_by_t,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 0.75,
    quadrature_nodes = 60L,
    batch_size = 500L
) {
  check_cdf_intersection_dependencies()
  if (!"X" %in% names(data)) {
    stop("The current DGP oracle CDF calculation requires X.")
  }
  rule <- gauss_legendre_rule(quadrature_nodes, -2, 2)
  u <- rule$nodes
  weights <- rule$weights * dtnormtc(u)
  weights <- weights / sum(weights)
  n <- nrow(data)
  cdf <- list()
  for (t_value in c(0L, 1L)) {
    for (z_value in c(0L, 1L)) {
      key <- paste0("t", t_value, "_z", z_value)
      cdf[[key]] <- matrix(
        NA_real_,
        nrow = n,
        ncol = length(thresholds_by_t[[as.character(t_value)]])
      )
    }
  }
  pi0 <- pi1 <- rep(NA_real_, n)
  starts <- seq.int(1L, n, by = batch_size)
  for (start in starts) {
    finish <- min(n, start + batch_size - 1L)
    index <- start:finish
    x <- data$X[index]
    x_matrix <- matrix(x, nrow = length(index), ncol = length(u))
    u_matrix <- matrix(
      u,
      nrow = length(index),
      ncol = length(u),
      byrow = TRUE
    )
    eta0 <- treatment_intercept + x_matrix +
      confounder_strength * u_matrix
    probability_t1 <- list(
      `0` = stats::plogis(eta0),
      `1` = stats::plogis(eta0 + instrument_strength)
    )
    pi0[index] <- as.vector(probability_t1[["0"]] %*% weights)
    pi1[index] <- as.vector(probability_t1[["1"]] %*% weights)

    for (t_value in c(0L, 1L)) {
      thresholds <- thresholds_by_t[[as.character(t_value)]]
      outcome_mean <- treatment_effect * t_value +
        1.5 * x_matrix + 1.5 * u_matrix
      for (z_value in c(0L, 1L)) {
        arm_weight <- if (t_value == 1L) {
          probability_t1[[as.character(z_value)]]
        } else {
          1 - probability_t1[[as.character(z_value)]]
        }
        denominator <- as.vector(arm_weight %*% weights)
        key <- paste0("t", t_value, "_z", z_value)
        for (threshold_index in seq_along(thresholds)) {
          conditional_event <- truncated_normal_cdf(
            value = thresholds[threshold_index],
            mean = outcome_mean
          )
          numerator <- as.vector(
            (conditional_event * arm_weight) %*% weights
          )
          cdf[[key]][index, threshold_index] <- numerator / denominator
        }
      }
    }
  }
  for (name in names(cdf)) {
    cdf[[name]] <- rearrange_cdf_rows(cdf_clip_probability(cdf[[name]]))
  }
  list(
    layer = "dgp_oracle",
    cdf = cdf,
    nuisance = list(
      e = stats::plogis(-0.5 * data$X),
      pi0 = cdf_clip_probability(pi0),
      pi1 = cdf_clip_probability(pi1)
    ),
    diagnostics = data.frame(
      quadrature_nodes = quadrature_nodes,
      max_cdf_monotonicity_violation = max(vapply(
        cdf,
        function(value) max(pmax(-t(apply(value, 1L, diff)), 0)),
        numeric(1)
      )),
      stringsAsFactors = FALSE
    )
  )
}

central_covariate_index <- function(
    data,
    x_cols = "X",
    trim_probability = 0.005
) {
  if (
    length(trim_probability) != 1L ||
      !is.finite(trim_probability) ||
      trim_probability < 0 || trim_probability >= 0.5
  ) {
    stop("trim_probability must lie in [0,0.5).")
  }
  keep <- rep(TRUE, nrow(data))
  for (name in x_cols) {
    limits <- stats::quantile(
      data[[name]],
      probs = c(trim_probability, 1 - trim_probability),
      type = 8,
      names = FALSE
    )
    keep <- keep & data[[name]] >= limits[1L] & data[[name]] <= limits[2L]
  }
  keep
}

upper_tail_mean <- function(value, fraction = 0.05) {
  if (
    length(fraction) != 1L || !is.finite(fraction) ||
      fraction <= 0 || fraction > 1
  ) {
    stop("fraction must lie in (0,1].")
  }
  count <- max(1L, as.integer(ceiling(fraction * length(value))))
  mean(sort(value, decreasing = TRUE)[seq_len(count)])
}

evaluate_cdf_intersection_layer <- function(
    data,
    layer,
    thresholds_by_t,
    lambda_grid = c(1, 1.25, 1.5, 2),
    x_cols = "X",
    trim_probability = 0.005,
    x_tail_fraction = 0.05,
    numerical_tolerance = 1e-10
) {
  lambda_grid <- sort(unique(lambda_grid))
  if (length(lambda_grid) == 0L || any(lambda_grid < 1)) {
    stop("lambda_grid must contain values of at least one.")
  }
  keep <- central_covariate_index(
    data = data,
    x_cols = x_cols,
    trim_probability = trim_probability
  )
  rows <- list()
  detail <- list()
  row_index <- 0L
  for (t_value in c(0L, 1L)) {
    thresholds <- thresholds_by_t[[as.character(t_value)]]
    cdf0 <- layer$cdf[[paste0("t", t_value, "_z0")]]
    cdf1 <- layer$cdf[[paste0("t", t_value, "_z1")]]
    pi0 <- if (t_value == 1L) layer$nuisance$pi0 else 1 - layer$nuisance$pi0
    pi1 <- if (t_value == 1L) layer$nuisance$pi1 else 1 - layer$nuisance$pi1
    for (Lambda in lambda_grid) {
      interval0 <- msm_binary_event_interval(
        observed_probability = as.vector(cdf0),
        arm_probability = rep(pi0, times = ncol(cdf0)),
        Lambda = Lambda
      )
      interval1 <- msm_binary_event_interval(
        observed_probability = as.vector(cdf1),
        arm_probability = rep(pi1, times = ncol(cdf1)),
        Lambda = Lambda
      )
      lower0 <- matrix(interval0$lower, nrow = nrow(cdf0))
      upper0 <- matrix(interval0$upper, nrow = nrow(cdf0))
      lower1 <- matrix(interval1$lower, nrow = nrow(cdf1))
      upper1 <- matrix(interval1$upper, nrow = nrow(cdf1))
      gap <- pmax(lower0, lower1) - pmin(upper0, upper1)
      positive_gap <- pmax(gap, 0)
      gap_kept <- gap[keep, , drop = FALSE]
      positive_kept <- positive_gap[keep, , drop = FALSE]
      x_score <- apply(positive_kept, 1L, max)
      maximum_index <- which(gap_kept == max(gap_kept), arr.ind = TRUE)[1L, ]
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        layer = layer$layer,
        treatment = t_value,
        Lambda = Lambda,
        observations = sum(keep),
        thresholds = length(thresholds),
        trim_probability = trim_probability,
        maximum_gap = max(gap_kept),
        mean_positive_gap = mean(positive_kept),
        rms_positive_gap = sqrt(mean(positive_kept^2)),
        positive_grid_fraction = mean(gap_kept > numerical_tolerance),
        x_any_positive_fraction = mean(
          apply(gap_kept > numerical_tolerance, 1L, any)
        ),
        mean_x_maximum_gap = mean(x_score),
        q95_x_maximum_gap = as.numeric(stats::quantile(
          x_score, 0.95, type = 8, names = FALSE
        )),
        q99_x_maximum_gap = as.numeric(stats::quantile(
          x_score, 0.99, type = 8, names = FALSE
        )),
        upper_tail_mean_x_gap = upper_tail_mean(
          x_score, fraction = x_tail_fraction
        ),
        x_at_maximum = data[[x_cols[1L]]][keep][maximum_index[1L]],
        threshold_at_maximum = thresholds[maximum_index[2L]],
        stringsAsFactors = FALSE
      )
      detail[[row_index]] <- data.frame(
        layer = layer$layer,
        treatment = t_value,
        Lambda = Lambda,
        row_id = which(keep),
        x_maximum_gap = x_score,
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    summary = do.call(rbind, rows),
    x_detail = do.call(rbind, detail),
    interpretation = paste(
      "Positive values measure nonoverlap of stratum-specific MSM CDF",
      "intervals. They are descriptive scores, not hypothesis-test results."
    )
  )
}

summarize_cdf_prediction_accuracy <- function(
    estimated_layer,
    oracle_layer
) {
  rows <- lapply(names(oracle_layer$cdf), function(name) {
    estimated <- estimated_layer$cdf[[name]]
    oracle <- oracle_layer$cdf[[name]]
    difference <- as.vector(estimated - oracle)
    components <- strcapture(
      "t([01])_z([01])",
      name,
      proto = list(treatment = integer(), instrument = integer())
    )
    data.frame(
      treatment = components$treatment,
      instrument = components$instrument,
      observations = nrow(oracle),
      thresholds = ncol(oracle),
      bias = mean(difference),
      rmse = sqrt(mean(difference^2)),
      mae = mean(abs(difference)),
      correlation = stats::cor(as.vector(estimated), as.vector(oracle)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
