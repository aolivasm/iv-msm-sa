# Diagnostic-only compatibility utilities.
#
# These functions assume R/sharp_bounds_plugin.R has already been sourced.
# They refit the final gamma RKHS learner at tuning values selected in a
# previous run, construct envelopes without compatibility projection, and
# evaluate the pointwise and conditional-normalization implications.

refit_gamma_rkhs_saved_tuning <- function(
    train,
    nuisance_train,
    t_value,
    tuning,
    x_cols = "X",
    fit_seed = 1L,
    gamma_floor = 1e-6,
    gamma_cap = 100,
    bandwidth_sample = 1000L,
    nystrom_eigen_tolerance = 1e-10,
    nystrom_jitter = 1e-10,
    reconstruction_tolerance = 1e-8
) {
  if (nrow(tuning) != 1L) {
    stop("tuning must contain exactly one saved fold-arm row.")
  }
  required <- c(
    "selected_lambda", "selected_bandwidth_multiplier", "base_bandwidth",
    "selected_bandwidth", "normalization_constant", "rkhs_backend",
    "nystrom_rank_requested", "nystrom_landmark_count",
    "nystrom_effective_rank"
  )
  missing <- setdiff(required, names(tuning))
  if (length(missing) > 0L) {
    stop("Saved tuning row omits: ", paste(missing, collapse = ", "))
  }
  if (!identical(as.character(tuning$rkhs_backend), "nystrom")) {
    stop("The production diagnostic refit expects the Nyström gamma backend.")
  }

  treatment_index <- which(train$T_obs == t_value)
  train_t <- train[treatment_index, , drop = FALSE]
  nuisance_t <- subset_nuisance(nuisance_train, treatment_index)
  requested_rank <- as.integer(tuning$nystrom_rank_requested)

  # These offsets reproduce make_rkhs_cv_plan()'s final-fit construction.
  coordinate_spec <- prepare_rkhs_coordinates(
    train_t = train_t,
    x_cols = x_cols,
    bandwidth_sample = bandwidth_sample,
    seed = fit_seed + 50000L
  )
  landmark_spec <- prepare_gamma_nystrom_landmarks(
    coordinate_spec = coordinate_spec,
    rank = requested_rank,
    seed = fit_seed + 60000L
  )
  feature_spec <- make_gamma_nystrom_features(
    coordinate_spec = coordinate_spec,
    landmark_spec = landmark_spec,
    bandwidth_multiplier = tuning$selected_bandwidth_multiplier,
    eigen_tolerance = nystrom_eigen_tolerance,
    numerical_jitter = nystrom_jitter
  )
  reduced_fit <- fit_gamma_rkhs_reduced(
    train_t = train_t,
    nuisance_t = nuisance_t,
    t_value = t_value,
    x_cols = x_cols,
    lambda = tuning$selected_lambda,
    bandwidth_multiplier = tuning$selected_bandwidth_multiplier,
    coordinate_spec = coordinate_spec,
    bandwidth_sample = bandwidth_sample,
    backend = "nystrom",
    nystrom_rank = requested_rank,
    nystrom_landmark_spec = landmark_spec,
    nystrom_feature_spec = feature_spec,
    nystrom_eigen_tolerance = nystrom_eigen_tolerance,
    nystrom_jitter = nystrom_jitter
  )

  raw_train <- reduced_fit$predict_raw(train_t)
  positive_train <- clip_positive(
    raw_train, floor = gamma_floor, cap = gamma_cap
  )
  normalization_constant <- max(
    mean(positive_train[train_t$Z == 0]), gamma_floor
  )

  relative_difference <- function(value, reference) {
    abs(value - reference) / max(1, abs(reference))
  }
  reconstruction <- data.frame(
    base_bandwidth_difference =
      reduced_fit$base_sigma - tuning$base_bandwidth,
    selected_bandwidth_difference =
      reduced_fit$sigma - tuning$selected_bandwidth,
    normalization_constant_difference =
      normalization_constant - tuning$normalization_constant,
    landmark_count_difference =
      reduced_fit$nystrom_landmark_count - tuning$nystrom_landmark_count,
    effective_rank_difference =
      reduced_fit$nystrom_effective_rank - tuning$nystrom_effective_rank,
    stringsAsFactors = FALSE
  )
  numerical_match <-
    relative_difference(reduced_fit$base_sigma, tuning$base_bandwidth) <=
      reconstruction_tolerance &&
    relative_difference(reduced_fit$sigma, tuning$selected_bandwidth) <=
      reconstruction_tolerance &&
    relative_difference(
      normalization_constant, tuning$normalization_constant
    ) <= reconstruction_tolerance &&
    reconstruction$landmark_count_difference == 0 &&
    reconstruction$effective_rank_difference == 0
  reconstruction$numerical_match <- numerical_match
  if (!numerical_match) {
    stop(
      "The fixed-tuning gamma refit did not reproduce the saved final fit."
    )
  }

  predict_gamma <- function(newdata) {
    raw <- reduced_fit$predict_raw(newdata)
    positive <- clip_positive(raw, floor = gamma_floor, cap = gamma_cap)
    clip_positive(
      positive / normalization_constant,
      floor = gamma_floor,
      cap = gamma_cap
    )
  }

  list(
    predict = predict_gamma,
    reconstruction = reconstruction,
    selected_lambda = tuning$selected_lambda,
    selected_bandwidth_multiplier = tuning$selected_bandwidth_multiplier,
    normalization_constant = normalization_constant,
    nystrom_landmark_indices = reduced_fit$nystrom_landmark_indices
  )
}

unprojected_envelope_matrices <- function(
    nuisance,
    gamma,
    t_value,
    lambda_grid
) {
  pi0 <- nuisance$pi0
  pi1 <- nuisance$pi1
  if (t_value == 1) {
    A <- ((1 - pi1) / (1 - pi0)) * gamma
    B <- (pi1 * gamma - pi0) / (1 - pi0)
  } else {
    A <- (pi1 / pi0) * gamma
    B <- ((1 - pi1) * gamma - (1 - pi0)) / pi0
  }
  n <- length(gamma)
  G <- length(lambda_grid)
  L <- U <- matrix(NA_real_, nrow = n, ncol = G)
  for (index in seq_along(lambda_grid)) {
    Lambda <- lambda_grid[index]
    L[, index] <- pmax(1 / Lambda, A / Lambda + B)
    U[, index] <- pmin(Lambda, A * Lambda + B)
  }
  colnames(L) <- colnames(U) <- format(lambda_grid, scientific = FALSE)
  list(L = L, U = U, A = A, B = B)
}

evaluate_unprojected_compatibility_fold <- function(
    train,
    validation,
    nuisance_train,
    nuisance_validation,
    gamma_train,
    gamma_validation,
    t_value,
    lambda_grid,
    x_cols = "X",
    conditional_spline_df = 5L,
    compatibility_tolerance = 1e-6
) {
  baseline_train <- train$T_obs == t_value & train$Z == 0
  baseline_validation <- validation$T_obs == t_value & validation$Z == 0
  if (sum(baseline_train) <= conditional_spline_df + 2L) {
    stop("Too few baseline training observations for the normalization fit.")
  }
  if (!any(baseline_validation)) {
    stop("The validation fold has no baseline observations for one arm.")
  }

  envelope_train <- unprojected_envelope_matrices(
    nuisance_train, gamma_train, t_value, lambda_grid
  )
  envelope_validation <- unprojected_envelope_matrices(
    nuisance_validation, gamma_validation, t_value, lambda_grid
  )
  L_train <- envelope_train$L[baseline_train, , drop = FALSE]
  U_train <- envelope_train$U[baseline_train, , drop = FALSE]
  L_validation <- envelope_validation$L[
    baseline_validation, , drop = FALSE
  ]
  U_validation <- envelope_validation$U[
    baseline_validation, , drop = FALSE
  ]

  # Fit all Lambda-specific conditional means with a common spline design.
  training_data <- train[baseline_train, , drop = FALSE]
  validation_data <- validation[baseline_validation, , drop = FALSE]
  design <- make_spline_design(
    data = training_data,
    x_cols = x_cols,
    spline_df = conditional_spline_df
  )
  response <- cbind(L_train, U_train)
  conditional_fit <- stats::lm.fit(x = design$matrix, y = response)
  coefficients <- conditional_fit$coefficients
  coefficients[!is.finite(coefficients)] <- 0
  conditional_prediction <- design$predict_matrix(validation_data) %*%
    coefficients
  G <- length(lambda_grid)
  Lbar_validation <- conditional_prediction[, seq_len(G), drop = FALSE]
  Ubar_validation <- conditional_prediction[, G + seq_len(G), drop = FALSE]

  pointwise_excess <- pmax(L_validation - U_validation, 0)
  normalization_lower_excess <- pmax(Lbar_validation - 1, 0)
  normalization_upper_excess <- pmax(1 - Ubar_validation, 0)
  normalization_excess <- normalization_lower_excess +
    normalization_upper_excess
  conditional_gap <- Ubar_validation - Lbar_validation

  summarize_matrix <- function(value, indicator) {
    list(
      count = colSums(indicator),
      rate = colMeans(indicator),
      mean_magnitude = colMeans(value),
      maximum_magnitude = apply(value, 2L, max)
    )
  }
  pointwise <- summarize_matrix(
    pointwise_excess, pointwise_excess > compatibility_tolerance
  )
  normalization <- summarize_matrix(
    normalization_excess,
    normalization_excess > compatibility_tolerance
  )

  data.frame(
    treatment = t_value,
    Lambda = lambda_grid,
    n_validation_baseline = sum(baseline_validation),
    pointwise_violation_count = pointwise$count,
    pointwise_violation_rate = pointwise$rate,
    pointwise_mean_violation_magnitude = pointwise$mean_magnitude,
    pointwise_maximum_violation_magnitude = pointwise$maximum_magnitude,
    normalization_violation_count = normalization$count,
    normalization_violation_rate = normalization$rate,
    normalization_lower_violation_rate = colMeans(
      normalization_lower_excess > compatibility_tolerance
    ),
    normalization_upper_violation_rate = colMeans(
      normalization_upper_excess > compatibility_tolerance
    ),
    normalization_mean_violation_magnitude =
      normalization$mean_magnitude,
    normalization_maximum_violation_magnitude =
      normalization$maximum_magnitude,
    conditional_degeneracy_rate = colMeans(
      conditional_gap <= compatibility_tolerance
    ),
    mean_Lbar = colMeans(Lbar_validation),
    mean_Ubar = colMeans(Ubar_validation),
    stringsAsFactors = FALSE
  )
}

aggregate_compatibility_fold_rates <- function(fold_rates) {
  split_rows <- split(
    fold_rates,
    interaction(fold_rates$treatment, fold_rates$Lambda, drop = TRUE)
  )
  rows <- lapply(split_rows, function(group) {
    n_total <- sum(group$n_validation_baseline)
    weighted_mean <- function(column) {
      sum(group$n_validation_baseline * group[[column]]) / n_total
    }
    data.frame(
      treatment = group$treatment[1L],
      Lambda = group$Lambda[1L],
      n_validation_baseline = n_total,
      pointwise_violation_count = sum(group$pointwise_violation_count),
      pointwise_violation_rate =
        sum(group$pointwise_violation_count) / n_total,
      pointwise_mean_violation_magnitude = weighted_mean(
        "pointwise_mean_violation_magnitude"
      ),
      pointwise_maximum_violation_magnitude = max(
        group$pointwise_maximum_violation_magnitude
      ),
      normalization_violation_count =
        sum(group$normalization_violation_count),
      normalization_violation_rate =
        sum(group$normalization_violation_count) / n_total,
      normalization_lower_violation_rate = weighted_mean(
        "normalization_lower_violation_rate"
      ),
      normalization_upper_violation_rate = weighted_mean(
        "normalization_upper_violation_rate"
      ),
      normalization_mean_violation_magnitude = weighted_mean(
        "normalization_mean_violation_magnitude"
      ),
      normalization_maximum_violation_magnitude = max(
        group$normalization_maximum_violation_magnitude
      ),
      conditional_degeneracy_rate = weighted_mean(
        "conditional_degeneracy_rate"
      ),
      mean_Lbar = weighted_mean("mean_Lbar"),
      mean_Ubar = weighted_mean("mean_Ubar"),
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output[order(output$Lambda, output$treatment), , drop = FALSE]
}

compatibility_frontier_from_rates <- function(
    arm_rates,
    violation_tolerances = c(0.10, 0.05, 0.01)
) {
  lambda_grid <- sort(unique(arm_rates$Lambda))
  by_lambda <- lapply(lambda_grid, function(Lambda) {
    group <- arm_rates[arm_rates$Lambda == Lambda, , drop = FALSE]
    if (!setequal(group$treatment, c(0, 1))) {
      stop("Every Lambda must have diagnostics for both treatment arms.")
    }
    data.frame(
      Lambda = Lambda,
      maximum_pointwise_violation_rate = max(
        group$pointwise_violation_rate
      ),
      maximum_normalization_violation_rate = max(
        group$normalization_violation_rate
      )
    )
  })
  rate_curve <- do.call(rbind, by_lambda)

  frontier_for_pass <- function(pass) {
    all_larger_pass <- rev(cumprod(rev(as.integer(pass)))) == 1L
    candidate <- which(all_larger_pass)
    if (length(candidate) == 0L) {
      return(list(
        minimum_nonrefuted_lambda = Inf,
        largest_refuted_grid_lambda = max(lambda_grid),
        raw_largest_flagged_lambda = if (any(!pass)) {
          max(lambda_grid[!pass])
        } else NA_real_,
        status = "above_grid"
      ))
    }
    first <- candidate[1L]
    list(
      minimum_nonrefuted_lambda = lambda_grid[first],
      largest_refuted_grid_lambda = if (first == 1L) {
        NA_real_
      } else lambda_grid[first - 1L],
      raw_largest_flagged_lambda = if (any(!pass)) {
        max(lambda_grid[!pass])
      } else NA_real_,
      status = if (first == 1L) "below_grid" else "within_grid"
    )
  }

  rows <- lapply(violation_tolerances, function(tolerance) {
    pointwise_pass <-
      rate_curve$maximum_pointwise_violation_rate <= tolerance
    normalization_pass <-
      rate_curve$maximum_normalization_violation_rate <= tolerance
    combined_pass <- pointwise_pass & normalization_pass
    pointwise <- frontier_for_pass(pointwise_pass)
    normalization <- frontier_for_pass(normalization_pass)
    combined <- frontier_for_pass(combined_pass)
    data.frame(
      violation_tolerance = tolerance,
      compatibility_percentage = 100 * (1 - tolerance),
      primary = tolerance == 0.10,
      pointwise_minimum_nonrefuted_lambda =
        pointwise$minimum_nonrefuted_lambda,
      normalization_minimum_nonrefuted_lambda =
        normalization$minimum_nonrefuted_lambda,
      combined_minimum_nonrefuted_lambda =
        combined$minimum_nonrefuted_lambda,
      combined_largest_refuted_grid_lambda =
        combined$largest_refuted_grid_lambda,
      combined_raw_largest_flagged_lambda =
        combined$raw_largest_flagged_lambda,
      pointwise_frontier_status = pointwise$status,
      normalization_frontier_status = normalization$status,
      combined_frontier_status = combined$status,
      stringsAsFactors = FALSE
    )
  })
  list(frontiers = do.call(rbind, rows), rate_curve = rate_curve)
}
