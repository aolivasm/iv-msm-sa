# Diagnostics for detecting unmeasured confounding with a binary instrument.
#
# This file supplies two deliberately distinct procedures:
#
# 1. A Guo-style IV propensity-score subclassification test.  Under
#    individual exclusion, instrument exchangeability including potential
#    treatments, monotonicity, and relevance, it compares identified
#    compliance-class mean potential outcomes.  These assumptions are
#    stronger than the IV restrictions imposed by the Dual MSM.
#
# 2. The heuristic Dual-MSM compatibility diagnostic already implemented in
#    R/compatibility_diagnostic.R, evaluated with cross-fitted estimated
#    nuisances and, when available, with oracle simulation nuisances.
#
# The second procedure remains a descriptive diagnostic.  Its violation-rate
# cutoffs are not calibrated test critical values.

fit_iv_propensity_subclasses <- function(
    data,
    x_cols = "X",
    n_subclasses = 5L,
    probability_epsilon = 0.01
) {
  check_binary(data$Z, "Z")
  if (
    length(n_subclasses) != 1L ||
      !is.finite(n_subclasses) ||
      n_subclasses < 2L ||
      abs(n_subclasses - round(n_subclasses)) >
        sqrt(.Machine$double.eps)
  ) {
    stop("n_subclasses must be one integer of at least 2.")
  }
  missing_x <- setdiff(x_cols, names(data))
  if (length(missing_x) > 0L) {
    stop("Missing IV propensity predictors: ", paste(missing_x, collapse = ", "))
  }

  formula <- stats::reformulate(x_cols, response = "Z")
  fit <- stats::glm(
    formula = formula,
    data = data,
    family = stats::binomial(),
    control = stats::glm.control(maxit = 100L)
  )
  score <- clip_probability(
    stats::predict(fit, newdata = data, type = "response"),
    probability_epsilon
  )
  probabilities <- seq(0, 1, length.out = as.integer(n_subclasses) + 1L)
  breaks <- unique(as.numeric(stats::quantile(
    score,
    probs = probabilities,
    names = FALSE,
    type = 8
  )))
  if (length(breaks) != as.integer(n_subclasses) + 1L) {
    stop(
      "The estimated IV propensity score does not support the requested ",
      "number of distinct subclasses."
    )
  }
  interior <- breaks[-c(1L, length(breaks))]
  subclass <- findInterval(score, interior) + 1L

  list(
    subclass = as.integer(subclass),
    score = as.numeric(score),
    fit = fit,
    breaks = breaks,
    n_subclasses = as.integer(n_subclasses)
  )
}

estimate_compliance_contrasts_by_subclass <- function(
    data,
    subclass,
    target_distribution = c("compliers", "full_sample"),
    minimum_complier_probability = 1e-4
) {
  target_distribution <- match.arg(target_distribution)
  required <- c("Y", "T_obs", "Z")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Data omit: ", paste(missing, collapse = ", "))
  }
  check_binary(data$T_obs, "T_obs")
  check_binary(data$Z, "Z")
  if (length(subclass) != nrow(data) || anyNA(subclass)) {
    stop("subclass must be observed once for every data row.")
  }

  subclass_values <- sort(unique(subclass))
  rows <- lapply(subclass_values, function(value) {
    in_subclass <- subclass == value
    z0 <- in_subclass & data$Z == 0
    z1 <- in_subclass & data$Z == 1
    if (!any(z0) || !any(z1)) {
      stop("Every IV propensity subclass must contain both instrument levels.")
    }

    probability_always <- mean(data$T_obs[z0])
    probability_never <- mean(1 - data$T_obs[z1])
    probability_complier <- 1 - probability_always - probability_never
    if (
      !is.finite(probability_complier) ||
        probability_complier <= minimum_complier_probability
    ) {
      stop(
        sprintf(
          "Estimated complier probability %.6f is too small in subclass %s.",
          probability_complier,
          as.character(value)
        )
      )
    }

    z0_treated <- z0 & data$T_obs == 1
    z1_treated <- z1 & data$T_obs == 1
    z0_control <- z0 & data$T_obs == 0
    z1_control <- z1 & data$T_obs == 0
    if (
      !any(z0_treated) || !any(z1_treated) ||
        !any(z0_control) || !any(z1_control)
    ) {
      stop(
        "Every subclass-by-instrument group must contain treated and control observations."
      )
    }

    mean_always_y1 <- mean(data$Y[z0_treated])
    mean_never_y0 <- mean(data$Y[z1_control])
    treated_mean_z1 <- mean(data$Y[z1_treated])
    control_mean_z0 <- mean(data$Y[z0_control])
    mean_complier_y1 <- (
      (probability_always + probability_complier) * treated_mean_z1 -
        probability_always * mean_always_y1
    ) / probability_complier
    mean_complier_y0 <- (
      (probability_never + probability_complier) * control_mean_z0 -
        probability_never * mean_never_y0
    ) / probability_complier

    data.frame(
      subclass = value,
      n_subclass = sum(in_subclass),
      n_z0 = sum(z0),
      n_z1 = sum(z1),
      probability_always = probability_always,
      probability_complier = probability_complier,
      probability_never = probability_never,
      mean_always_y1 = mean_always_y1,
      mean_complier_y1 = mean_complier_y1,
      mean_never_y0 = mean_never_y0,
      mean_complier_y0 = mean_complier_y0,
      treated_contrast = mean_always_y1 - mean_complier_y1,
      control_contrast = mean_never_y0 - mean_complier_y0,
      stringsAsFactors = FALSE
    )
  })
  detail <- do.call(rbind, rows)
  rownames(detail) <- NULL
  unnormalized_weights <- if (target_distribution == "compliers") {
    detail$n_subclass * detail$probability_complier
  } else {
    detail$n_subclass
  }
  weights <- unnormalized_weights / sum(unnormalized_weights)
  estimates <- c(
    treated_contrast = sum(weights * detail$treated_contrast),
    control_contrast = sum(weights * detail$control_contrast)
  )

  list(
    estimates = estimates,
    subclass_detail = detail,
    target_distribution = if (target_distribution == "compliers") {
      "estimated distribution of X among compliers"
    } else {
      "empirical distribution of X in the full sample"
    }
  )
}

resample_pairs_within_iv_subclasses <- function(data, subclass) {
  output <- data
  cells <- interaction(subclass, data$Z, drop = TRUE)
  for (cell in levels(cells)) {
    index <- which(cells == cell)
    sampled <- index[sample.int(
      length(index), size = length(index), replace = TRUE
    )]
    output$Y[index] <- data$Y[sampled]
    output$T_obs[index] <- data$T_obs[sampled]
  }
  output
}

estimate_guo_subclassification_test <- function(
    data,
    x_cols = "X",
    n_subclasses = 5L,
    target_distribution = c("compliers", "full_sample"),
    bootstrap_replications = 199L,
    confidence_level = 0.95,
    probability_epsilon = 0.01,
    minimum_complier_probability = 1e-4,
    seed = 1L
) {
  target_distribution <- match.arg(target_distribution)
  if (
    length(bootstrap_replications) != 1L ||
      !is.finite(bootstrap_replications) ||
      bootstrap_replications < 19L ||
      abs(bootstrap_replications - round(bootstrap_replications)) >
        sqrt(.Machine$double.eps)
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

  subclass_fit <- fit_iv_propensity_subclasses(
    data = data,
    x_cols = x_cols,
    n_subclasses = n_subclasses,
    probability_epsilon = probability_epsilon
  )
  point <- estimate_compliance_contrasts_by_subclass(
    data = data,
    subclass = subclass_fit$subclass,
    target_distribution = target_distribution,
    minimum_complier_probability = minimum_complier_probability
  )

  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_random_seed) {
    old_random_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", old_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  bootstrap <- matrix(
    NA_real_,
    nrow = as.integer(bootstrap_replications),
    ncol = 2L,
    dimnames = list(NULL, names(point$estimates))
  )
  valid_count <- 0L
  attempted_count <- 0L
  maximum_attempts <- 10L * as.integer(bootstrap_replications)
  while (
    valid_count < as.integer(bootstrap_replications) &&
      attempted_count < maximum_attempts
  ) {
    attempted_count <- attempted_count + 1L
    bootstrap_data <- resample_pairs_within_iv_subclasses(
      data,
      subclass_fit$subclass
    )
    draw <- tryCatch(
      estimate_compliance_contrasts_by_subclass(
        data = bootstrap_data,
        subclass = subclass_fit$subclass,
        target_distribution = target_distribution,
        minimum_complier_probability = minimum_complier_probability
      )$estimates,
      error = function(error) rep(NA_real_, 2L)
    )
    if (all(is.finite(draw))) {
      valid_count <- valid_count + 1L
      bootstrap[valid_count, ] <- draw
    }
  }
  valid <- stats::complete.cases(bootstrap)
  if (sum(valid) < bootstrap_replications) {
    stop(
      sprintf(
        paste0(
          "Only %d valid Guo-style bootstrap draws were obtained in %d ",
          "attempts."
        ),
        sum(valid), attempted_count
      )
    )
  }
  bootstrap_valid <- bootstrap[valid, , drop = FALSE]
  alpha <- 1 - confidence_level
  intervals <- apply(
    bootstrap_valid,
    2L,
    stats::quantile,
    probs = c(alpha / 2, 1 - alpha / 2),
    names = FALSE,
    type = 8
  )
  if (is.null(dim(intervals))) {
    intervals <- matrix(intervals, nrow = 2L)
  }
  standard_errors <- apply(bootstrap_valid, 2L, stats::sd)
  p_values <- 2 * stats::pnorm(
    -abs(point$estimates / pmax(standard_errors, 1e-12))
  )
  rejected <- intervals[1L, ] > 0 | intervals[2L, ] < 0
  contrasts <- data.frame(
    contrast = names(point$estimates),
    estimate = as.numeric(point$estimates),
    bootstrap_standard_error = as.numeric(standard_errors),
    confidence_level = confidence_level,
    confidence_lower = as.numeric(intervals[1L, ]),
    confidence_upper = as.numeric(intervals[2L, ]),
    normal_approximation_p_value = as.numeric(p_values),
    rejected_by_percentile_interval = as.logical(rejected),
    stringsAsFactors = FALSE
  )

  covariance <- stats::cov(bootstrap_valid)
  eigen_covariance <- eigen(covariance, symmetric = TRUE)
  cutoff <- max(eigen_covariance$values) * 1e-10
  inverse_values <- ifelse(
    eigen_covariance$values > cutoff,
    1 / eigen_covariance$values,
    0
  )
  covariance_inverse <- eigen_covariance$vectors %*%
    diag(inverse_values, nrow = length(inverse_values)) %*%
    t(eigen_covariance$vectors)
  joint_statistic <- as.numeric(
    t(point$estimates) %*% covariance_inverse %*% point$estimates
  )
  joint_rank <- sum(eigen_covariance$values > cutoff)
  joint_p_value <- if (joint_rank > 0L) {
    stats::pchisq(joint_statistic, df = joint_rank, lower.tail = FALSE)
  } else {
    NA_real_
  }

  overall <- data.frame(
    confidence_level = confidence_level,
    rejected_any_percentile_interval = any(rejected),
    joint_wald_statistic = joint_statistic,
    joint_wald_degrees_freedom = joint_rank,
    joint_wald_p_value = joint_p_value,
    joint_wald_rejected = is.finite(joint_p_value) && joint_p_value < alpha,
    requested_bootstrap_replications = as.integer(bootstrap_replications),
    valid_bootstrap_replications = sum(valid),
    attempted_bootstrap_replications = attempted_count,
    invalid_bootstrap_rate = 1 - sum(valid) / attempted_count,
    minimum_subclass_complier_probability = min(
      point$subclass_detail$probability_complier
    ),
    minimum_iv_propensity = min(subclass_fit$score),
    maximum_iv_propensity = max(subclass_fit$score),
    stringsAsFactors = FALSE
  )

  list(
    contrasts = contrasts,
    overall = overall,
    subclass_detail = point$subclass_detail,
    target_distribution = point$target_distribution,
    iv_propensity_breaks = subclass_fit$breaks,
    bootstrap_draws = bootstrap,
    assumptions = paste(
      "Guo-style comparator: conditional IV independence including potential",
      "treatments, individual-level exclusion, monotonicity, relevance, and",
      "approximate covariate balance within IV propensity subclasses."
    ),
    inferential_note = paste(
      "Separate percentile intervals are not multiplicity adjusted; the joint",
      "Wald result is also reported. The bootstrap fixes X, Z, and subclasses",
      "and resamples observed (Y,T) pairs within subclass-by-Z cells. Draws",
      "that violate the estimated monotonic first stage are discarded and",
      "their rate is reported."
    )
  )
}

make_compatibility_refutation_rows <- function(
    arm_rates,
    violation_tolerances = c(0.10, 0.05, 0.01)
) {
  required <- c(
    "treatment", "Lambda", "pointwise_violation_rate",
    "normalization_violation_rate"
  )
  missing <- setdiff(required, names(arm_rates))
  if (length(missing) > 0L) {
    stop("arm_rates omit: ", paste(missing, collapse = ", "))
  }
  layers <- if ("nuisance_layer" %in% names(arm_rates)) {
    unique(arm_rates$nuisance_layer)
  } else {
    "unspecified"
  }
  rows <- list()
  counter <- 1L
  for (layer in layers) {
    layer_rows <- if ("nuisance_layer" %in% names(arm_rates)) {
      arm_rates[arm_rates$nuisance_layer == layer, , drop = FALSE]
    } else {
      arm_rates
    }
    for (Lambda in sort(unique(layer_rows$Lambda))) {
      group <- layer_rows[layer_rows$Lambda == Lambda, , drop = FALSE]
      if (!setequal(group$treatment, c(0, 1))) {
        stop("Every nuisance layer and Lambda must contain both treatment arms.")
      }
      for (tolerance in violation_tolerances) {
        pointwise <- any(group$pointwise_violation_rate > tolerance)
        normalization <- any(group$normalization_violation_rate > tolerance)
        rows[[counter]] <- data.frame(
          nuisance_layer = layer,
          Lambda = Lambda,
          violation_tolerance = tolerance,
          pointwise_refuted = pointwise,
          normalization_refuted = normalization,
          combined_refuted = pointwise || normalization,
          maximum_pointwise_violation_rate = max(
            group$pointwise_violation_rate
          ),
          maximum_normalization_violation_rate = max(
            group$normalization_violation_rate
          ),
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

fit_crossfitted_compatibility_calibration <- function(
    data,
    lambda_grid = c(1, 2, 3),
    violation_tolerances = c(0.10, 0.05, 0.01),
    x_cols = "X",
    K = 5L,
    nuisance_gam_k = 6L,
    conditional_spline_df = 5L,
    compatibility_tolerance = 1e-6,
    gamma_lambda_constants = 10^seq(-1.5, 1, by = 0.5),
    gamma_bandwidth_multipliers = c(0.25, 0.5, 1, 2, 4),
    gamma_cv_folds = 3L,
    gamma_selection_rule = "minimum",
    gamma_nystrom_rank = NULL,
    include_oracle_layer = TRUE,
    seed = 1L
) {
  if (include_oracle_layer) {
    oracle_columns <- c(
      "e_true", "pi0_true", "pi1_true",
      "gamma0_oracle", "gamma1_oracle"
    )
    missing_oracle <- setdiff(oracle_columns, names(data))
    if (length(missing_oracle) > 0L) {
      stop(
        "The requested oracle layer is missing: ",
        paste(missing_oracle, collapse = ", ")
      )
    }
  }
  if (!".row_id" %in% names(data)) data$.row_id <- seq_len(nrow(data))
  fold_id <- make_tz_stratified_folds(data, K = K, seed = seed)
  gamma_provider <- make_cached_gamma_rkhs_provider(
    lambda_grid = function(n_training, treatment) {
      gamma_lambda_constants / sqrt(n_training)
    },
    bandwidth_multipliers = gamma_bandwidth_multipliers,
    cv_folds = gamma_cv_folds,
    selection_rule = gamma_selection_rule,
    backend = "nystrom",
    nystrom_rank = gamma_nystrom_rank
  )

  fold_rows <- list()
  counter <- 1L
  for (fold in seq_len(K)) {
    validation_index <- which(fold_id == fold)
    training_index <- which(fold_id != fold)
    train <- data[training_index, , drop = FALSE]
    validation <- data[validation_index, , drop = FALSE]
    nuisance_fit <- fit_nuisance_models(
      train = train,
      x_cols = x_cols,
      positivity_epsilon = 0.05,
      learner = "gam",
      gam_basis_dimension = nuisance_gam_k,
      gam_method = "REML"
    )
    nuisance_train <- nuisance_fit$predict(train)
    nuisance_validation <- nuisance_fit$predict(validation)

    for (t_value in c(0, 1)) {
      fit_seed <- seed + 100000L * fold + 100L * t_value
      gamma_fit <- gamma_provider(
        fold = fold,
        treatment = t_value,
        train = train,
        nuisance_train = nuisance_train,
        x_cols = x_cols,
        seed = fit_seed
      )
      estimated <- evaluate_unprojected_compatibility_fold(
        train = train,
        validation = validation,
        nuisance_train = nuisance_train,
        nuisance_validation = nuisance_validation,
        gamma_train = gamma_fit$predict(train)$normalized,
        gamma_validation = gamma_fit$predict(validation)$normalized,
        t_value = t_value,
        lambda_grid = lambda_grid,
        x_cols = x_cols,
        conditional_spline_df = conditional_spline_df,
        compatibility_tolerance = compatibility_tolerance
      )
      estimated$fold <- fold
      estimated$nuisance_layer <- "estimated_rkhs"
      fold_rows[[counter]] <- estimated
      counter <- counter + 1L

      if (include_oracle_layer) {
        oracle_train <- list(
          e = train$e_true,
          pi0 = train$pi0_true,
          pi1 = train$pi1_true
        )
        oracle_validation <- list(
          e = validation$e_true,
          pi0 = validation$pi0_true,
          pi1 = validation$pi1_true
        )
        gamma_column <- if (t_value == 1) {
          "gamma1_oracle"
        } else {
          "gamma0_oracle"
        }
        oracle <- evaluate_unprojected_compatibility_fold(
          train = train,
          validation = validation,
          nuisance_train = oracle_train,
          nuisance_validation = oracle_validation,
          gamma_train = train[[gamma_column]],
          gamma_validation = validation[[gamma_column]],
          t_value = t_value,
          lambda_grid = lambda_grid,
          x_cols = x_cols,
          conditional_spline_df = conditional_spline_df,
          compatibility_tolerance = compatibility_tolerance
        )
        oracle$fold <- fold
        oracle$nuisance_layer <- "oracle"
        fold_rows[[counter]] <- oracle
        counter <- counter + 1L
      }
    }
  }

  fold_rates <- do.call(rbind, fold_rows)
  rownames(fold_rates) <- NULL
  layer_values <- unique(fold_rates$nuisance_layer)
  arm_rates <- do.call(rbind, lapply(layer_values, function(layer) {
    local <- aggregate_compatibility_fold_rates(
      fold_rates[fold_rates$nuisance_layer == layer, , drop = FALSE]
    )
    local$nuisance_layer <- layer
    local
  }))
  rownames(arm_rates) <- NULL
  refutations <- make_compatibility_refutation_rows(
    arm_rates,
    violation_tolerances = violation_tolerances
  )

  cached <- collect_cached_gamma_rkhs_fits(gamma_provider)
  tuning <- do.call(rbind, lapply(cached, function(entry) {
    fit <- entry$fit
    data.frame(
      fold = entry$fold,
      treatment = entry$treatment,
      selected_lambda = fit$selected_lambda,
      selected_bandwidth_multiplier = fit$selected_bandwidth_multiplier,
      base_bandwidth = fit$base_bandwidth,
      selected_bandwidth = fit$selected_bandwidth,
      normalization_constant = fit$normalization_constant,
      training_negative_rate = fit$training_negative_rate,
      nystrom_rank_requested = fit$nystrom_rank_requested,
      nystrom_landmark_count = fit$nystrom_landmark_count,
      nystrom_effective_rank = fit$nystrom_effective_rank,
      stringsAsFactors = FALSE
    )
  }))
  rownames(tuning) <- NULL

  list(
    fold_rates = fold_rates,
    arm_rates = arm_rates,
    refutations = refutations,
    gamma_tuning = tuning,
    fold_id = fold_id,
    configuration = list(
      lambda_grid = lambda_grid,
      violation_tolerances = violation_tolerances,
      nuisance_learner = "gam",
      nuisance_gam_k = nuisance_gam_k,
      gamma_lambda_grid =
        "10^seq(-1.5,1,by=0.5)/sqrt(n_t_training)",
      gamma_bandwidth_multipliers = gamma_bandwidth_multipliers,
      gamma_selection_rule = gamma_selection_rule,
      gamma_backend = "nystrom",
      compatibility_projection_used = FALSE
    )
  )
}
