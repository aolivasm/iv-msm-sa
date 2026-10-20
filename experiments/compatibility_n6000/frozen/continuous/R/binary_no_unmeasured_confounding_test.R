# Direct test of no unmeasured confounding for binary outcomes.
#
# Under the maintained distributional IV assumptions and positivity, the
# Dual-MSM null Lambda0 = Lambda1 = 1 implies
#
#   P(Y = 1 | T = t, X, Z = 1) = P(Y = 1 | T = t, X, Z = 0),  t = 0, 1.
#
# Because Y is binary, these two conditional-mean equalities are equivalent to
# equality of the corresponding conditional outcome distributions.  This file
# tests a finite collection of implications of those equalities.  It does not
# require monotonicity, compliance classes, or individual-level exclusion.
#
# Required first:
#   source("R/binary_outcome_bounds.R")
#
# The test uses cross-fitted outcome regressions m_tz(X) and the conditional
# instrument propensity g_t(X) = P(Z = 1 | T = t, X).  The latter is obtained
# by Bayes' rule from cross-fitted estimates of e(X) and pi_z(X).  For each arm,
# the doubly robust pseudo-outcome
#
#   D_t = m_t1 - m_t0
#       + Z / g_t       * (Y - m_t1)
#       - (1-Z)/(1-g_t) * (Y - m_t0)
#
# satisfies E(D_t | T=t, X) = m_t1,true(X)-m_t0,true(X) when either g_t or both
# outcome regressions are consistently estimated.  The primary test applies a
# two-sided Gaussian multiplier maximum to prespecified series moments of D_t.
# A separately reported augmentation adds the cross-fitted learned direction
# m_t1(X)-m_t0(X) to the witness dictionary.

check_binary_no_unmeasured_confounding_dependencies <- function() {
  required <- c(
    "binary_check_indicator", "binary_clip_probability",
    "binary_make_tz_folds", "fit_binary_cross_fitted_nuisances"
  )
  missing <- required[
    !vapply(required, exists, logical(1), mode = "function")
  ]
  if (length(missing) > 0L) {
    stop(
      "Source R/binary_outcome_bounds.R first. Missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

make_binary_conditional_independence_basis <- function(
    data,
    x_cols,
    spline_df = 3L,
    include_pairwise = TRUE,
    max_columns = 40L
) {
  missing <- setdiff(x_cols, names(data))
  if (length(missing) > 0L) {
    stop("Missing test-basis covariates: ", paste(missing, collapse = ", "))
  }
  non_numeric <- x_cols[
    !vapply(data[, x_cols, drop = FALSE], is.numeric, logical(1))
  ]
  if (length(non_numeric) > 0L) {
    stop(
      "Test-basis covariates must be numeric: ",
      paste(non_numeric, collapse = ", ")
    )
  }
  spline_df <- as.integer(spline_df)
  max_columns <- as.integer(max_columns)
  if (length(spline_df) != 1L || is.na(spline_df) || spline_df < 1L) {
    stop("spline_df must be one positive integer.")
  }
  if (length(max_columns) != 1L || is.na(max_columns) || max_columns < 1L) {
    stop("max_columns must be one positive integer.")
  }

  n <- nrow(data)
  columns <- list(intercept = rep(1, n))
  first_order <- list()
  for (variable in x_cols) {
    value <- as.numeric(data[[variable]])
    if (any(!is.finite(value))) {
      stop("Non-finite values in test-basis covariate ", variable, ".")
    }
    unique_count <- length(unique(value))
    if (unique_count <= 1L) next
    centered <- value - mean(value)
    scale <- sqrt(mean(centered^2))
    if (!is.finite(scale) || scale <= 1e-12) next
    standardized <- centered / scale
    first_order[[variable]] <- standardized

    if (unique_count <= 2L || spline_df == 1L) {
      columns[[paste0("main_", variable)]] <- standardized
    } else {
      local_df <- min(spline_df, unique_count - 1L)
      spline <- splines::ns(value, df = local_df)
      for (index in seq_len(ncol(spline))) {
        local <- as.numeric(spline[, index])
        local <- local - mean(local)
        local_scale <- sqrt(mean(local^2))
        if (is.finite(local_scale) && local_scale > 1e-12) {
          columns[[paste0("spline_", variable, "_", index)]] <-
            local / local_scale
        }
      }
    }
  }

  if (isTRUE(include_pairwise) && length(first_order) >= 2L) {
    pairs <- utils::combn(names(first_order), 2L, simplify = FALSE)
    for (pair in pairs) {
      if (length(columns) >= max_columns) break
      product <- first_order[[pair[1L]]] * first_order[[pair[2L]]]
      product <- product - mean(product)
      product_scale <- sqrt(mean(product^2))
      if (is.finite(product_scale) && product_scale > 1e-12) {
        columns[[paste0("interaction_", pair[1L], "_", pair[2L])]] <-
          product / product_scale
      }
    }
  }

  if (length(columns) > max_columns) {
    columns <- columns[seq_len(max_columns)]
  }
  basis <- do.call(cbind, columns)
  if (is.null(dim(basis))) basis <- matrix(basis, ncol = 1L)
  colnames(basis) <- names(columns)

  # Remove exact linear dependencies without using Y, T, or Z.  The first
  # max_columns pivots are deterministic functions of the observed X design.
  decomposition <- qr(basis, tol = 1e-10)
  keep <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  basis <- basis[, keep, drop = FALSE]
  attr(basis, "configuration") <- list(
    x_cols = x_cols,
    spline_df = spline_df,
    include_pairwise = isTRUE(include_pairwise),
    max_columns = max_columns
  )
  basis
}

binary_conditional_instrument_probability <- function(
    predictions,
    treatment,
    epsilon = 0.025
) {
  required <- c("e", "pi0", "pi1")
  missing <- setdiff(required, names(predictions))
  if (length(missing) > 0L) {
    stop("Nuisance predictions omit: ", paste(missing, collapse = ", "))
  }
  if (!treatment %in% c(0, 1)) stop("treatment must equal zero or one.")
  r0 <- if (treatment == 1L) predictions$pi0 else 1 - predictions$pi0
  r1 <- if (treatment == 1L) predictions$pi1 else 1 - predictions$pi1
  probability_t <- predictions$e * r0 + (1 - predictions$e) * r1
  raw <- (1 - predictions$e) * r1 / pmax(probability_t, 1e-12)
  list(
    probability = binary_clip_probability(raw, epsilon),
    raw = raw,
    treatment_probability = probability_t,
    lower_truncation_rate = mean(raw < epsilon),
    upper_truncation_rate = mean(raw > 1 - epsilon)
  )
}

make_binary_conditional_independence_signals <- function(
    data,
    nuisance_fit,
    x_cols,
    instrument_overlap_epsilon = 0.025,
    series_spline_df = 3L,
    include_pairwise = TRUE,
    max_basis_columns = 40L,
    include_learned_direction = TRUE,
    learned_direction_floor = 0.01
) {
  check_binary_no_unmeasured_confounding_dependencies()
  required <- c("Y", "T_obs", "Z", x_cols)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Data omit: ", paste(missing, collapse = ", "))
  }
  binary_check_indicator(data$Y, "Y")
  binary_check_indicator(data$T_obs, "T_obs")
  binary_check_indicator(data$Z, "Z")
  if (is.null(nuisance_fit$predictions)) {
    stop("nuisance_fit must contain cross-fitted predictions.")
  }
  predictions <- nuisance_fit$predictions
  required_predictions <- c(
    "e", "pi0", "pi1",
    "a_t0_z0", "a_t0_z1", "a_t1_z0", "a_t1_z1"
  )
  missing_predictions <- setdiff(required_predictions, names(predictions))
  if (length(missing_predictions) > 0L) {
    stop(
      "Nuisance fit omits: ", paste(missing_predictions, collapse = ", ")
    )
  }
  if (nrow(predictions) != nrow(data) ||
      any(!is.finite(as.matrix(predictions[, required_predictions])))) {
    stop("Cross-fitted nuisance predictions are invalid.")
  }

  basis <- make_binary_conditional_independence_basis(
    data = data,
    x_cols = x_cols,
    spline_df = series_spline_df,
    include_pairwise = include_pairwise,
    max_columns = max_basis_columns
  )
  primary_signals <- augmented_signals <- list()
  metadata <- diagnostics <- list()
  learned_metadata <- list()

  for (treatment in c(0, 1)) {
    arm <- data$T_obs == treatment
    arm_probability <- mean(arm)
    if (!is.finite(arm_probability) || arm_probability <= 0) {
      stop("A treatment arm is empty.")
    }
    instrument <- binary_conditional_instrument_probability(
      predictions = predictions,
      treatment = treatment,
      epsilon = instrument_overlap_epsilon
    )
    m0_name <- sprintf("a_t%d_z0", treatment)
    m1_name <- sprintf("a_t%d_z1", treatment)
    m0 <- as.numeric(predictions[[m0_name]])
    m1 <- as.numeric(predictions[[m1_name]])
    g <- instrument$probability
    pseudo_outcome <-
      m1 - m0 +
      data$Z / g * (data$Y - m1) -
      (1 - data$Z) / (1 - g) * (data$Y - m0)

    arm_series <- basis * as.numeric(arm) * pseudo_outcome / arm_probability
    colnames(arm_series) <- paste0(
      "t", treatment, "_", colnames(basis)
    )
    primary_signals[[as.character(treatment)]] <- arm_series
    augmented_signals[[as.character(treatment)]] <- arm_series
    metadata[[as.character(treatment)]] <- data.frame(
      moment = colnames(arm_series),
      treatment = treatment,
      witness = colnames(basis),
      witness_type = "prespecified_series",
      stringsAsFactors = FALSE
    )

    learned_difference <- m1 - m0
    learned_scale <- sqrt(mean(learned_difference[arm]^2))
    learned_active <- isTRUE(include_learned_direction) &&
      is.finite(learned_scale) && learned_scale > learned_direction_floor
    if (learned_active) {
      learned_witness <- pmax(
        -1, pmin(1, learned_difference / learned_scale)
      )
      learned_signal <- matrix(
        as.numeric(arm) * learned_witness * pseudo_outcome / arm_probability,
        ncol = 1L,
        dimnames = list(NULL, paste0("t", treatment, "_learned_direction"))
      )
      augmented_signals[[paste0(treatment, "_learned")]] <- learned_signal
      learned_metadata[[as.character(treatment)]] <- data.frame(
        moment = colnames(learned_signal),
        treatment = treatment,
        witness = "cross_fitted_m1_minus_m0",
        witness_type = "learned_direction",
        stringsAsFactors = FALSE
      )
    }
    diagnostics[[as.character(treatment)]] <- data.frame(
      treatment = treatment,
      observations = sum(arm),
      treatment_fraction = arm_probability,
      raw_g_minimum = min(instrument$raw[arm]),
      raw_g_maximum = max(instrument$raw[arm]),
      g_lower_truncation_rate = mean(
        instrument$raw[arm] < instrument_overlap_epsilon
      ),
      g_upper_truncation_rate = mean(
        instrument$raw[arm] > 1 - instrument_overlap_epsilon
      ),
      learned_direction_rms = learned_scale,
      learned_direction_included = learned_active,
      stringsAsFactors = FALSE
    )
  }

  primary <- do.call(cbind, primary_signals)
  augmented <- do.call(cbind, augmented_signals)
  primary_metadata <- do.call(rbind, metadata)
  augmented_metadata <- rbind(
    primary_metadata,
    if (length(learned_metadata) > 0L) {
      do.call(rbind, learned_metadata)
    } else {
      primary_metadata[FALSE, , drop = FALSE]
    }
  )
  augmented <- augmented[, augmented_metadata$moment, drop = FALSE]
  rownames(primary_metadata) <- rownames(augmented_metadata) <- NULL

  list(
    primary = primary,
    primary_metadata = primary_metadata,
    augmented = augmented,
    augmented_metadata = augmented_metadata,
    arm_diagnostics = do.call(rbind, diagnostics),
    basis = basis,
    basis_configuration = attr(basis, "configuration")
  )
}

calibrate_binary_conditional_independence_maximum <- function(
    signals,
    metadata,
    bootstrap_replications = 999L,
    confidence_level = 0.95,
    seed = 1L,
    bootstrap_batch_size = 100L
) {
  signals <- as.matrix(signals)
  n <- nrow(signals)
  if (n < 2L || ncol(signals) < 1L || nrow(metadata) != ncol(signals)) {
    stop("Invalid signal matrix or metadata.")
  }
  bootstrap_replications <- as.integer(bootstrap_replications)
  bootstrap_batch_size <- as.integer(bootstrap_batch_size)
  if (bootstrap_replications < 99L) {
    stop("Use at least 99 multiplier-bootstrap replications.")
  }
  if (bootstrap_batch_size < 1L) {
    stop("bootstrap_batch_size must be positive.")
  }
  if (confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must lie strictly between zero and one.")
  }

  estimate <- colMeans(signals)
  standard_deviation <- apply(signals, 2L, stats::sd)
  keep <- is.finite(estimate) & is.finite(standard_deviation) &
    standard_deviation > 1e-10
  if (!any(keep)) stop("Every conditional-independence moment is degenerate.")
  signals <- signals[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  estimate <- estimate[keep]
  standard_deviation <- standard_deviation[keep]
  z_value <- sqrt(n) * estimate / standard_deviation
  observed_joint <- max(abs(z_value))
  arm_values <- sort(unique(metadata$treatment))
  observed_arm <- stats::setNames(vapply(arm_values, function(treatment) {
    max(abs(z_value[metadata$treatment == treatment]))
  }, numeric(1)), arm_values)

  centered_standardized <- sweep(signals, 2L, estimate, "-")
  centered_standardized <- sweep(
    centered_standardized, 2L, standard_deviation, "/"
  )
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))

  bootstrap_joint <- numeric(bootstrap_replications)
  bootstrap_arm <- matrix(
    NA_real_, nrow = bootstrap_replications, ncol = length(arm_values),
    dimnames = list(NULL, as.character(arm_values))
  )
  completed <- 0L
  while (completed < bootstrap_replications) {
    batch <- min(bootstrap_batch_size, bootstrap_replications - completed)
    multipliers <- matrix(stats::rnorm(batch * n), nrow = batch, ncol = n)
    draw <- multipliers %*% centered_standardized / sqrt(n)
    index <- completed + seq_len(batch)
    bootstrap_joint[index] <- apply(abs(draw), 1L, max)
    for (arm_index in seq_along(arm_values)) {
      columns <- metadata$treatment == arm_values[arm_index]
      bootstrap_arm[index, arm_index] <- apply(
        abs(draw[, columns, drop = FALSE]), 1L, max
      )
    }
    completed <- completed + batch
  }

  alpha <- 1 - confidence_level
  joint <- data.frame(
    statistic = observed_joint,
    critical_value = as.numeric(stats::quantile(
      bootstrap_joint, confidence_level, names = FALSE, type = 8
    )),
    p_value = (1 + sum(bootstrap_joint >= observed_joint)) /
      (bootstrap_replications + 1),
    reject = FALSE,
    moments_tested = ncol(signals),
    stringsAsFactors = FALSE
  )
  joint$reject <- joint$p_value <= alpha

  arm_summary <- do.call(rbind, lapply(seq_along(arm_values), function(index) {
    treatment <- arm_values[index]
    local_bootstrap <- bootstrap_arm[, index]
    data.frame(
      treatment = treatment,
      statistic = observed_arm[[as.character(treatment)]],
      critical_value = as.numeric(stats::quantile(
        local_bootstrap, confidence_level, names = FALSE, type = 8
      )),
      p_value = (1 + sum(
        local_bootstrap >= observed_arm[[as.character(treatment)]]
      )) / (bootstrap_replications + 1),
      reject = FALSE,
      moments_tested = sum(metadata$treatment == treatment),
      stringsAsFactors = FALSE
    )
  }))
  arm_summary$reject <- arm_summary$p_value <= alpha

  moment_detail <- cbind(
    metadata,
    estimate = estimate,
    standard_error = standard_deviation / sqrt(n),
    z_value = z_value,
    two_sided_normal_p_value = 2 * stats::pnorm(-abs(z_value))
  )
  moment_detail <- moment_detail[
    order(moment_detail$treatment, -abs(moment_detail$z_value)),
    , drop = FALSE
  ]
  rownames(moment_detail) <- NULL

  list(
    joint = joint,
    arm_summary = arm_summary,
    moment_detail = moment_detail,
    bootstrap_joint = bootstrap_joint,
    bootstrap_arm = bootstrap_arm
  )
}

test_binary_no_unmeasured_confounding <- function(
    data,
    x_cols = "X",
    folds = 5L,
    fold_id = NULL,
    nuisance_fit = NULL,
    nuisance_learner = c(
      "library_cv", "super_learner", "gam", "glm_spline", "ranger",
      "glmnet"
    ),
    nuisance_spline_df = 5L,
    nuisance_gam_basis_dimension = NULL,
    nuisance_library_candidates = c(
      "glm_spline", "gam", "glmnet", "ranger"
    ),
    nuisance_library_cv_folds = 3L,
    nuisance_library_selection_rule = c("one_se", "minimum"),
    nuisance_ranger_num_trees = 500L,
    nuisance_super_learner_candidates =
      binary_default_super_learner_candidates(),
    nuisance_super_learner_cv_folds = 3L,
    nuisance_super_learner_ranger_num_trees = 500L,
    nuisance_super_learner_hal_max_degree = 2L,
    nuisance_overlap_epsilon = 0.05,
    outcome_probability_epsilon = 1e-4,
    instrument_overlap_epsilon = 0.025,
    series_spline_df = 3L,
    include_pairwise = TRUE,
    max_basis_columns = 40L,
    include_learned_direction = TRUE,
    learned_direction_floor = 0.01,
    bootstrap_replications = 999L,
    confidence_level = 0.95,
    seed = 1L
) {
  check_binary_no_unmeasured_confounding_dependencies()
  nuisance_learner <- match.arg(nuisance_learner)
  nuisance_library_selection_rule <- match.arg(
    nuisance_library_selection_rule
  )
  required <- c("Y", "T_obs", "Z", x_cols)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Data omit: ", paste(missing, collapse = ", "))
  }
  binary_check_indicator(data$Y, "Y")
  binary_check_indicator(data$T_obs, "T_obs")
  binary_check_indicator(data$Z, "Z")
  if (is.null(fold_id)) {
    fold_id <- binary_make_tz_folds(data, K = folds, seed = seed)
  }
  if (length(fold_id) != nrow(data)) stop("fold_id has the wrong length.")

  if (is.null(nuisance_fit)) {
    nuisance_fit <- fit_binary_cross_fitted_nuisances(
      data = data,
      x_cols = x_cols,
      K = folds,
      fold_id = fold_id,
      spline_df = nuisance_spline_df,
      learner = nuisance_learner,
      overlap_epsilon = nuisance_overlap_epsilon,
      outcome_probability_epsilon = outcome_probability_epsilon,
      gam_basis_dimension = nuisance_gam_basis_dimension,
      library_candidates = nuisance_library_candidates,
      library_cv_folds = nuisance_library_cv_folds,
      library_selection_rule = nuisance_library_selection_rule,
      library_ranger_num_trees = nuisance_ranger_num_trees,
      super_learner_candidates = nuisance_super_learner_candidates,
      super_learner_cv_folds = nuisance_super_learner_cv_folds,
      super_learner_ranger_num_trees =
        nuisance_super_learner_ranger_num_trees,
      super_learner_hal_max_degree =
        nuisance_super_learner_hal_max_degree,
      seed = seed + 10000L
    )
  }

  signals <- make_binary_conditional_independence_signals(
    data = data,
    nuisance_fit = nuisance_fit,
    x_cols = x_cols,
    instrument_overlap_epsilon = instrument_overlap_epsilon,
    series_spline_df = series_spline_df,
    include_pairwise = include_pairwise,
    max_basis_columns = max_basis_columns,
    include_learned_direction = include_learned_direction,
    learned_direction_floor = learned_direction_floor
  )
  primary <- calibrate_binary_conditional_independence_maximum(
    signals = signals$primary,
    metadata = signals$primary_metadata,
    bootstrap_replications = bootstrap_replications,
    confidence_level = confidence_level,
    seed = seed + 20000L
  )
  augmented <- calibrate_binary_conditional_independence_maximum(
    signals = signals$augmented,
    metadata = signals$augmented_metadata,
    bootstrap_replications = bootstrap_replications,
    confidence_level = confidence_level,
    seed = seed + 30000L
  )

  primary$joint$test <- "prespecified_series"
  augmented$joint$test <- "series_plus_cross_fitted_learned_direction"
  primary$arm_summary$test <- "prespecified_series"
  augmented$arm_summary$test <-
    "series_plus_cross_fitted_learned_direction"
  primary$moment_detail$test <- "prespecified_series"
  augmented$moment_detail$test <-
    "series_plus_cross_fitted_learned_direction"

  list(
    summary = rbind(primary$joint, augmented$joint),
    arm_summary = rbind(primary$arm_summary, augmented$arm_summary),
    moment_detail = rbind(primary$moment_detail, augmented$moment_detail),
    arm_diagnostics = signals$arm_diagnostics,
    nuisance_fit = nuisance_fit,
    basis = signals$basis,
    fold_id = fold_id,
    primary = primary,
    augmented = augmented,
    null = paste(
      "Under the maintained distributional IV assumptions:",
      "P(Y=1|T=t,X,Z=1)=P(Y=1|T=t,X,Z=0) for t=0,1."
    ),
    interpretation = paste(
      "Rejection falsifies the joint model containing the maintained IV",
      "assumptions and Lambda0=Lambda1=1; it does not identify which",
      "maintained assumption failed."
    ),
    configuration = list(
      x_cols = x_cols,
      folds = length(unique(fold_id)),
      nuisance_learner = nuisance_learner,
      nuisance_library_candidates = nuisance_library_candidates,
      nuisance_library_cv_folds = nuisance_library_cv_folds,
      nuisance_library_selection_rule = nuisance_library_selection_rule,
      nuisance_super_learner_candidates =
        nuisance_super_learner_candidates,
      nuisance_super_learner_cv_folds =
        nuisance_super_learner_cv_folds,
      nuisance_super_learner_ranger_num_trees =
        nuisance_super_learner_ranger_num_trees,
      nuisance_super_learner_hal_max_degree =
        nuisance_super_learner_hal_max_degree,
      nuisance_overlap_epsilon = nuisance_overlap_epsilon,
      outcome_probability_epsilon = outcome_probability_epsilon,
      instrument_overlap_epsilon = instrument_overlap_epsilon,
      series_spline_df = series_spline_df,
      include_pairwise = include_pairwise,
      max_basis_columns = max_basis_columns,
      include_learned_direction = include_learned_direction,
      bootstrap_replications = bootstrap_replications,
      confidence_level = confidence_level,
      seed = seed
    )
  )
}
