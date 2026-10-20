# Predictable online one-step estimation for binary Dual-MSM bounds.
#
# This module supplies the nuisance-refitting layer that is intentionally
# absent from R/online_one_step_bounds.R.  For every evaluation block, all
# nuisance regressions are fit using observations that occur strictly earlier
# in a randomized analysis order.  A held-out subset of that past estimates
# the predictable endpoint variances.  The evaluation block is never used to
# choose a learner, fit a nuisance regression, or estimate its variance weight.
#
# Internal cross-validation may be used *within the past training sample* to
# choose a nuisance learner.  That is tuning, not ordinary outer cross-fitting:
# independence from the current evaluation block follows from temporal sample
# splitting.

binary_fit_past_nuisances <- function(
    training_data,
    prediction_data,
    x_cols,
    spline_df = 3L,
    learner = c("library_cv", "gam", "glm_spline", "ranger", "glmnet"),
    overlap_epsilon = 0.05,
    outcome_probability_epsilon = 1e-4,
    gam_basis_dimension = NULL,
    gam_min_basis = 5L,
    gam_max_basis = 20L,
    gam_growth_exponent = 1 / 3,
    library_candidates = c("glm_spline", "gam", "glmnet", "ranger"),
    library_cv_folds = 3L,
    library_selection_rule = c("one_se", "minimum"),
    library_ranger_num_trees = 500L,
    seed = 1L
) {
  learner <- match.arg(learner)
  library_selection_rule <- match.arg(library_selection_rule)
  required_functions <- c(
    "binary_fit_regression_with_learner", "binary_clip_probability",
    "binary_check_indicator"
  )
  absent <- required_functions[!vapply(
    required_functions, exists, logical(1), mode = "function"
  )]
  if (length(absent) > 0L) {
    stop(
      "Source R/binary_outcome_bounds.R before this module; missing: ",
      paste(absent, collapse = ", "), "."
    )
  }
  required_columns <- c("Y", "T_obs", "Z", x_cols)
  for (label in c("training_data", "prediction_data")) {
    value <- get(label)
    missing <- setdiff(required_columns, names(value))
    if (length(missing) > 0L) {
      stop(label, " omits: ", paste(missing, collapse = ", "), ".")
    }
    if (nrow(value) < 1L) stop(label, " must contain observations.")
    binary_check_indicator(value$Y, paste0(label, "$Y"))
    binary_check_indicator(value$T_obs, paste0(label, "$T_obs"))
    binary_check_indicator(value$Z, paste0(label, "$Z"))
  }
  validate_epsilon <- function(value, label) {
    if (length(value) != 1L || !is.finite(value) ||
        value < 0 || value >= 0.5) {
      stop(label, " must be one number in [0,0.5).")
    }
  }
  validate_epsilon(overlap_epsilon, "overlap_epsilon")
  validate_epsilon(
    outcome_probability_epsilon, "outcome_probability_epsilon"
  )

  training_sets <- list(
    e = training_data,
    pi0 = training_data[training_data$Z == 0, , drop = FALSE],
    pi1 = training_data[training_data$Z == 1, , drop = FALSE],
    a_t0_z0 = training_data[
      training_data$T_obs == 0 & training_data$Z == 0, , drop = FALSE
    ],
    a_t0_z1 = training_data[
      training_data$T_obs == 0 & training_data$Z == 1, , drop = FALSE
    ],
    a_t1_z0 = training_data[
      training_data$T_obs == 1 & training_data$Z == 0, , drop = FALSE
    ],
    a_t1_z1 = training_data[
      training_data$T_obs == 1 & training_data$Z == 1, , drop = FALSE
    ]
  )
  if (any(vapply(training_sets, nrow, integer(1)) == 0L)) {
    empty <- names(training_sets)[vapply(training_sets, nrow, integer(1)) == 0L]
    stop(
      "The past training sample has empty nuisance cells: ",
      paste(empty, collapse = ", "), "."
    )
  }
  responses <- c(
    e = "Z", pi0 = "T_obs", pi1 = "T_obs",
    a_t0_z0 = "Y", a_t0_z1 = "Y", a_t1_z0 = "Y", a_t1_z1 = "Y"
  )
  component_seed <- c(
    e = 1L, pi0 = 2L, pi1 = 3L,
    a_t0_z0 = 4L, a_t0_z1 = 5L, a_t1_z0 = 6L, a_t1_z1 = 7L
  )
  fits <- lapply(names(training_sets), function(component) {
    binary_fit_regression_with_learner(
      data = training_sets[[component]],
      response = responses[[component]],
      x_cols = x_cols,
      spline_df = spline_df,
      learner = learner,
      gam_basis_dimension = gam_basis_dimension,
      gam_min_basis = gam_min_basis,
      gam_max_basis = gam_max_basis,
      gam_growth_exponent = gam_growth_exponent,
      library_candidates = library_candidates,
      library_cv_folds = library_cv_folds,
      library_selection_rule = library_selection_rule,
      ranger_num_trees = library_ranger_num_trees,
      seed = as.integer(seed + component_seed[[component]])
    )
  })
  names(fits) <- names(training_sets)

  raw_predictions <- lapply(fits, function(fit) fit$predict(prediction_data))
  raw_predictions$e <- 1 - raw_predictions$e
  predictions <- as.data.frame(raw_predictions, stringsAsFactors = FALSE)
  predictions$e <- binary_clip_probability(
    predictions$e, overlap_epsilon
  )
  predictions$pi0 <- binary_clip_probability(
    predictions$pi0, overlap_epsilon
  )
  predictions$pi1 <- binary_clip_probability(
    predictions$pi1, overlap_epsilon
  )
  outcome_components <- grep("^a_", names(predictions), value = TRUE)
  for (component in outcome_components) {
    predictions[[component]] <- binary_clip_probability(
      predictions[[component]], outcome_probability_epsilon
    )
  }
  if (any(!is.finite(as.matrix(predictions)))) {
    stop("A past-trained nuisance prediction is non-finite.")
  }

  scalar_or <- function(value, name, default) {
    candidate <- value[[name]]
    if (is.null(candidate) || length(candidate) == 0L) default else candidate[1L]
  }
  diagnostics <- do.call(rbind, lapply(names(fits), function(component) {
    fit <- fits[[component]]
    training <- training_sets[[component]]
    response <- responses[[component]]
    epsilon <- if (component %in% c("e", "pi0", "pi1")) {
      overlap_epsilon
    } else {
      outcome_probability_epsilon
    }
    raw <- raw_predictions[[component]]
    data.frame(
      component = component,
      response = response,
      training_observations = nrow(training),
      training_events = sum(training[[response]]),
      learner = as.character(scalar_or(fit, "learner", learner)),
      selected_learner = as.character(scalar_or(
        fit, "selected_learner", scalar_or(fit, "learner", learner)
      )),
      converged = isTRUE(scalar_or(fit, "converged", FALSE)),
      rank = as.numeric(scalar_or(fit, "rank", NA_real_)),
      constant_fit = isTRUE(scalar_or(fit, "constant", FALSE)),
      cv_mean_risk = as.numeric(scalar_or(fit, "cv_mean_risk", NA_real_)),
      cv_standard_error = as.numeric(scalar_or(
        fit, "cv_standard_error", NA_real_
      )),
      truncation_epsilon = epsilon,
      prediction_lower_truncation_rate = mean(raw < epsilon),
      prediction_upper_truncation_rate = mean(raw > 1 - epsilon),
      stringsAsFactors = FALSE
    )
  }))
  cv_tables <- lapply(names(fits), function(component) {
    table <- fits[[component]]$cv_table
    if (is.null(table) || nrow(table) == 0L) return(NULL)
    table$component <- component
    table
  })
  cv_tables <- cv_tables[!vapply(cv_tables, is.null, logical(1))]

  list(
    predictions = predictions,
    diagnostics = diagnostics,
    cv_table = if (length(cv_tables) > 0L) {
      do.call(rbind, cv_tables)
    } else {
      data.frame()
    },
    training_observations = nrow(training_data),
    prediction_observations = nrow(prediction_data),
    configuration = list(
      learner = learner,
      spline_df = as.integer(spline_df),
      overlap_epsilon = overlap_epsilon,
      outcome_probability_epsilon = outcome_probability_epsilon,
      library_candidates = library_candidates,
      library_cv_folds = as.integer(library_cv_folds),
      library_selection_rule = library_selection_rule,
      library_ranger_num_trees = as.integer(library_ranger_num_trees)
    )
  )
}

binary_online_endpoint_signals <- function(bound_result) {
  if (!is.list(bound_result) || !all(c("0", "1") %in% names(
    bound_result$signals
  ))) {
    stop("bound_result does not contain both arm-specific signal lists.")
  }
  arm0 <- bound_result$signals[["0"]]
  arm1 <- bound_result$signals[["1"]]
  output <- data.frame(
    mu1_lower = arm1$lower,
    mu1_upper = arm1$upper,
    mu0_lower = arm0$lower,
    mu0_upper = arm0$upper,
    ate_lower = arm1$lower - arm0$upper,
    ate_upper = arm1$upper - arm0$lower,
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(as.matrix(output)))) {
    stop("A binary online endpoint signal is non-finite.")
  }
  output
}

binary_online_bind_rows <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  completed <- lapply(rows, function(row) {
    row <- as.data.frame(row, stringsAsFactors = FALSE)
    for (column in setdiff(columns, names(row))) row[[column]] <- NA
    row[, columns, drop = FALSE]
  })
  output <- do.call(rbind, completed)
  rownames(output) <- NULL
  output
}

estimate_binary_online_dual_bounds <- function(
    data,
    x_cols,
    Lambda_grid,
    epsilon = 0.01,
    confidence_level = 0.95,
    initial_size = NULL,
    initial_exponent = 0.8,
    minimum_initial_size = 400L,
    chunks = 3L,
    variance_calibration_folds = 5L,
    standard_deviation_floor = 0.05,
    seed = 1L,
    permute = TRUE,
    spline_df = 3L,
    nuisance_learner = c(
      "library_cv", "gam", "glm_spline", "ranger", "glmnet"
    ),
    overlap_epsilon = 0.05,
    outcome_probability_epsilon = 1e-4,
    gam_basis_dimension = NULL,
    gam_min_basis = 5L,
    gam_max_basis = 20L,
    gam_growth_exponent = 1 / 3,
    nuisance_library_candidates = c(
      "glm_spline", "gam", "glmnet", "ranger"
    ),
    nuisance_library_cv_folds = 3L,
    nuisance_library_selection_rule = c("one_se", "minimum"),
    nuisance_library_ranger_num_trees = 500L,
    compatibility_action = c("project", "none"),
    compatibility_tolerance = 1e-12
) {
  nuisance_learner <- match.arg(nuisance_learner)
  nuisance_library_selection_rule <- match.arg(
    nuisance_library_selection_rule
  )
  compatibility_action <- match.arg(compatibility_action)
  required_functions <- c(
    "make_online_chunk_plan", "make_online_signal_chunk",
    "estimate_predictable_online_sd", "aggregate_online_one_step_bounds",
    "evaluate_binary_bounds_from_nuisances", "binary_make_tz_folds"
  )
  absent <- required_functions[!vapply(
    required_functions, exists, logical(1), mode = "function"
  )]
  if (length(absent) > 0L) {
    stop(
      "Source R/binary_outcome_bounds.R and R/online_one_step_bounds.R ",
      "before this module; missing: ", paste(absent, collapse = ", "), "."
    )
  }
  Lambda_grid <- sort(unique(as.numeric(Lambda_grid)))
  if (length(Lambda_grid) == 0L || any(!is.finite(Lambda_grid)) ||
      any(Lambda_grid < 1)) {
    stop("Lambda_grid must contain finite values of at least one.")
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon < 0 ||
      any(epsilon >= 1 / Lambda_grid)) {
    stop("epsilon must be nonnegative and smaller than every 1/Lambda.")
  }
  n <- nrow(data)
  chunks <- as.integer(chunks)
  minimum_initial_size <- as.integer(minimum_initial_size)
  variance_calibration_folds <- as.integer(variance_calibration_folds)
  if (is.null(initial_size)) {
    if (length(initial_exponent) != 1L || !is.finite(initial_exponent) ||
        initial_exponent <= 0.5 || initial_exponent >= 1) {
      stop("initial_exponent must lie strictly between 0.5 and one.")
    }
    initial_size <- max(
      minimum_initial_size, as.integer(ceiling(n^initial_exponent))
    )
  }
  initial_size <- min(as.integer(initial_size), n - chunks)
  if (variance_calibration_folds < 2L ||
      length(standard_deviation_floor) != 1L ||
      !is.finite(standard_deviation_floor) ||
      standard_deviation_floor <= 0) {
    stop("Invalid variance-calibration controls.")
  }

  plan <- make_online_chunk_plan(
    n = n,
    initial_size = initial_size,
    chunks = chunks,
    seed = seed,
    permute = permute
  )
  ordered_data <- data[plan$order, , drop = FALSE]
  rownames(ordered_data) <- NULL
  ordered_data$.online_original_row_id <- plan$order
  lambda_keys <- format(Lambda_grid, scientific = FALSE, trim = TRUE)
  online_chunks <- stats::setNames(
    lapply(Lambda_grid, function(value) vector("list", chunks)),
    lambda_keys
  )
  fit_diagnostics <- compatibility_diagnostics <- variance_diagnostics <-
    cv_tables <- list()

  for (chunk in seq_len(chunks)) {
    local_plan <- plan$plan[chunk, , drop = FALSE]
    training_end <- local_plan$training_end
    evaluation_position <- seq.int(
      local_plan$evaluation_start, local_plan$evaluation_end
    )
    past_data <- ordered_data[seq_len(training_end), , drop = FALSE]
    past_strata <- interaction(past_data$T_obs, past_data$Z, drop = TRUE)
    minimum_stratum <- min(table(past_strata))
    local_calibration_folds <- min(
      variance_calibration_folds, as.integer(minimum_stratum)
    )
    if (local_calibration_folds < 2L) {
      stop("Too few past observations for predictable variance calibration.")
    }
    calibration_assignment <- binary_make_tz_folds(
      past_data,
      K = local_calibration_folds,
      seed = as.integer(seed + 100L + chunk)
    )
    calibration_position <- which(calibration_assignment == 1L)
    fitting_position <- setdiff(seq_len(training_end), calibration_position)
    if (length(calibration_position) < 2L || length(fitting_position) < 2L) {
      stop("The past fitting/calibration split is too small.")
    }
    prediction_position <- c(calibration_position, evaluation_position)
    prediction_data <- ordered_data[prediction_position, , drop = FALSE]
    nuisance_fit <- binary_fit_past_nuisances(
      training_data = ordered_data[fitting_position, , drop = FALSE],
      prediction_data = prediction_data,
      x_cols = x_cols,
      spline_df = spline_df,
      learner = nuisance_learner,
      overlap_epsilon = overlap_epsilon,
      outcome_probability_epsilon = outcome_probability_epsilon,
      gam_basis_dimension = gam_basis_dimension,
      gam_min_basis = gam_min_basis,
      gam_max_basis = gam_max_basis,
      gam_growth_exponent = gam_growth_exponent,
      library_candidates = nuisance_library_candidates,
      library_cv_folds = nuisance_library_cv_folds,
      library_selection_rule = nuisance_library_selection_rule,
      library_ranger_num_trees = nuisance_library_ranger_num_trees,
      seed = as.integer(seed + 1000L * chunk)
    )
    local_fit_diagnostics <- nuisance_fit$diagnostics
    local_fit_diagnostics$chunk <- chunk
    local_fit_diagnostics$training_end <- training_end
    local_fit_diagnostics$fitting_observations <- length(fitting_position)
    local_fit_diagnostics$calibration_observations <-
      length(calibration_position)
    local_fit_diagnostics$evaluation_observations <-
      length(evaluation_position)
    local_fit_diagnostics$maximum_fit_source_position <-
      max(fitting_position)
    local_fit_diagnostics$minimum_evaluation_position <-
      min(evaluation_position)
    fit_diagnostics[[length(fit_diagnostics) + 1L]] <-
      local_fit_diagnostics
    if (nrow(nuisance_fit$cv_table) > 0L) {
      local_cv <- nuisance_fit$cv_table
      local_cv$chunk <- chunk
      local_cv$training_end <- training_end
      cv_tables[[length(cv_tables) + 1L]] <- local_cv
    }

    calibration_rows <- seq_along(calibration_position)
    evaluation_rows <- length(calibration_position) +
      seq_along(evaluation_position)
    for (lambda_index in seq_along(Lambda_grid)) {
      Lambda <- Lambda_grid[lambda_index]
      key <- lambda_keys[lambda_index]
      bound_result <- evaluate_binary_bounds_from_nuisances(
        data = prediction_data,
        nuisance_fit = nuisance_fit,
        Lambda0 = Lambda,
        Lambda1 = Lambda,
        epsilon = epsilon,
        confidence_level = confidence_level,
        compatibility_action = compatibility_action,
        compatibility_tolerance = compatibility_tolerance
      )
      signals <- binary_online_endpoint_signals(bound_result)
      calibration_signals <- signals[calibration_rows, , drop = FALSE]
      calibration_signals$mode <- "application"
      predictable_sd <- estimate_predictable_online_sd(
        calibration_signals,
        mode = "application",
        standard_deviation_floor = standard_deviation_floor,
        return_diagnostics = TRUE
      )
      online_chunks[[key]][[chunk]] <- make_online_signal_chunk(
        chunk = chunk,
        training_end = training_end,
        sequence_position = evaluation_position,
        signals = signals[evaluation_rows, , drop = FALSE],
        predictable_standard_deviation =
          predictable_sd$standard_deviation,
        standard_deviation_source_end = training_end,
        fit_source_end = training_end,
        predictable_standard_deviation_raw =
          predictable_sd$raw_standard_deviation,
        standard_deviation_floor =
          predictable_sd$standard_deviation_floor,
        variance_calibration_observations =
          predictable_sd$calibration_observations
      )
      local_variance <- data.frame(
        Lambda = Lambda,
        chunk = chunk,
        endpoint = names(predictable_sd$standard_deviation),
        training_end = training_end,
        fitting_observations = length(fitting_position),
        calibration_observations = length(calibration_position),
        evaluation_observations = length(evaluation_position),
        raw_standard_deviation = as.numeric(
          predictable_sd$raw_standard_deviation
        ),
        standard_deviation = as.numeric(
          predictable_sd$standard_deviation
        ),
        floor_active = as.logical(predictable_sd$floor_active),
        stringsAsFactors = FALSE
      )
      variance_diagnostics[[length(variance_diagnostics) + 1L]] <-
        local_variance
      local_compatibility <- bound_result$diagnostics
      local_compatibility$Lambda <- Lambda
      local_compatibility$chunk <- chunk
      local_compatibility$training_end <- training_end
      local_compatibility$evaluation_observations <-
        length(evaluation_position)
      compatibility_diagnostics[[length(compatibility_diagnostics) + 1L]] <-
        local_compatibility
    }
  }

  results <- lapply(seq_along(Lambda_grid), function(lambda_index) {
    Lambda <- Lambda_grid[lambda_index]
    key <- lambda_keys[lambda_index]
    value <- aggregate_online_one_step_bounds(
      online_chunks[[key]], confidence_level = confidence_level
    )
    value$summary$Lambda <- Lambda
    value$ate_outward_confidence_set$Lambda <- Lambda
    value$signals$Lambda <- Lambda
    value
  })
  names(results) <- lambda_keys
  endpoint_summary <- binary_online_bind_rows(lapply(results, `[[`, "summary"))
  ate_confidence_sets <- binary_online_bind_rows(lapply(
    results, `[[`, "ate_outward_confidence_set"
  ))
  point_rows <- lapply(results, function(value) {
    lower <- value$summary[
      value$summary$endpoint == "ate_lower", , drop = FALSE
    ]
    upper <- value$summary[
      value$summary$endpoint == "ate_upper", , drop = FALSE
    ]
    data.frame(
      Lambda = lower$Lambda,
      method = "dual_smooth_online_one_step",
      point_lower = lower$estimate,
      point_upper = upper$estimate,
      confidence_lower = lower$outward_confidence_limit,
      confidence_upper = upper$outward_confidence_limit,
      confidence_level = confidence_level,
      stringsAsFactors = FALSE
    )
  })
  point_and_confidence_sets <- binary_online_bind_rows(point_rows)
  signals <- binary_online_bind_rows(lapply(results, `[[`, "signals"))
  plan_table <- plan$plan
  plan_table$minimum_fit_source_position <- 1L
  plan_table$predictable <-
    plan_table$training_end < plan_table$evaluation_start

  list(
    point_and_confidence_sets = point_and_confidence_sets,
    endpoint_summary = endpoint_summary,
    ate_outward_confidence_sets = ate_confidence_sets,
    signals = signals,
    online_plan = plan_table,
    fit_diagnostics = binary_online_bind_rows(fit_diagnostics),
    nuisance_cv = binary_online_bind_rows(cv_tables),
    variance_diagnostics = binary_online_bind_rows(variance_diagnostics),
    compatibility_diagnostics = binary_online_bind_rows(
      compatibility_diagnostics
    ),
    results_by_lambda = results,
    settings = list(
      observations = n,
      evaluated_observations = n - initial_size,
      initial_size = initial_size,
      initial_exponent = initial_exponent,
      chunks = chunks,
      variance_calibration_folds = variance_calibration_folds,
      standard_deviation_floor = standard_deviation_floor,
      seed = seed,
      permuted = isTRUE(permute),
      Lambda_grid = Lambda_grid,
      epsilon = epsilon,
      confidence_level = confidence_level,
      x_cols = x_cols,
      nuisance_learner = nuisance_learner,
      construction = paste(
        "past-only nuisance fitting, held-past variance calibration,",
        "and future-block evaluation"
      )
    )
  )
}
