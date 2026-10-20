# Online aggregation for fixed-threshold one-step bound signals.
#
# This file implements the martingale aggregation layer adapted from
# Luedtke and van der Laan (2018), building on their 2016 optimal-treatment-
# rule construction.  It does not refit the sharp-bound
# nuisances.  Each signal chunk must have been produced by fitting q and all
# nuisances using observations preceding that chunk.  The existing
# fixed-threshold signal implementation in R/sharp_bounds_plugin.R supplies
# the six required row-level signals.
#
# The resulting outward limits are justified only under the high-level
# predictable variance and averaged nuisance-drift conditions stated in
# docs/estimation_section_revised.tex.  They are not made valid merely by
# passing the mechanical checks in this file.

online_bound_endpoint_names <- function() {
  c(
    "mu1_lower", "mu1_upper", "mu0_lower", "mu0_upper",
    "ate_lower", "ate_upper"
  )
}

make_online_chunk_plan <- function(
    n,
    initial_size,
    chunks = 5L,
    seed = 1L,
    permute = TRUE
) {
  n <- as.integer(n)
  initial_size <- as.integer(initial_size)
  chunks <- as.integer(chunks)
  if (
    length(n) != 1L || is.na(n) || n < 3L ||
      length(initial_size) != 1L || is.na(initial_size) ||
      initial_size < 2L || initial_size >= n ||
      length(chunks) != 1L || is.na(chunks) || chunks < 1L ||
      chunks > n - initial_size
  ) {
    stop("Invalid n, initial_size, or chunks for the online plan.")
  }
  sequence_rows <- seq_len(n)
  if (isTRUE(permute)) {
    set.seed(seed)
    sequence_rows <- sample(sequence_rows)
  }
  evaluation_positions <- seq.int(initial_size + 1L, n)
  chunk_id <- rep(seq_len(chunks), length.out = length(evaluation_positions))
  chunk_id <- sort(chunk_id)
  plan <- lapply(seq_len(chunks), function(chunk) {
    positions <- evaluation_positions[chunk_id == chunk]
    training_end <- min(positions) - 1L
    data.frame(
      chunk = chunk,
      training_end = training_end,
      evaluation_start = min(positions),
      evaluation_end = max(positions),
      evaluation_observations = length(positions),
      stringsAsFactors = FALSE
    )
  })
  list(
    order = sequence_rows,
    plan = do.call(rbind, plan),
    initial_size = initial_size,
    chunks = chunks,
    permuted = isTRUE(permute),
    seed = if (isTRUE(permute)) seed else NA_integer_
  )
}

make_online_signal_chunk <- function(
    chunk,
    training_end,
    sequence_position,
    signals,
    predictable_standard_deviation,
    standard_deviation_source_end = training_end,
    fit_source_end = training_end,
    predictable_standard_deviation_raw =
      predictable_standard_deviation,
    standard_deviation_floor = NA_real_,
    variance_calibration_observations = NA_integer_
) {
  endpoints <- online_bound_endpoint_names()
  signals <- as.data.frame(signals, stringsAsFactors = FALSE)
  missing <- setdiff(endpoints, names(signals))
  if (length(missing) > 0L) {
    stop("Online signals omit: ", paste(missing, collapse = ", "))
  }
  sequence_position <- as.integer(sequence_position)
  training_end <- as.integer(training_end)
  standard_deviation_source_end <- as.integer(
    standard_deviation_source_end
  )
  fit_source_end <- as.integer(fit_source_end)
  if (
    nrow(signals) == 0L ||
      length(sequence_position) != nrow(signals) ||
      anyNA(sequence_position) ||
      anyDuplicated(sequence_position) ||
      any(sequence_position <= training_end) ||
      length(training_end) != 1L || is.na(training_end) || training_end < 1L ||
      length(standard_deviation_source_end) != 1L ||
      is.na(standard_deviation_source_end) ||
      standard_deviation_source_end > training_end ||
      length(fit_source_end) != 1L || is.na(fit_source_end) ||
      fit_source_end > training_end ||
      any(!is.finite(as.matrix(signals[, endpoints, drop = FALSE])))
  ) {
    stop("The online chunk violates ordering, dimensions, or finiteness.")
  }
  predictable_standard_deviation <- unlist(
    predictable_standard_deviation,
    use.names = TRUE
  )
  if (is.null(names(predictable_standard_deviation))) {
    stop("Predictable standard deviations must be named by endpoint.")
  }
  missing_sd <- setdiff(endpoints, names(predictable_standard_deviation))
  if (
    length(missing_sd) > 0L ||
      any(!is.finite(predictable_standard_deviation[endpoints])) ||
      any(predictable_standard_deviation[endpoints] <= 0)
  ) {
    stop("Predictable standard deviations are missing or nonpositive.")
  }
  predictable_standard_deviation_raw <- unlist(
    predictable_standard_deviation_raw,
    use.names = TRUE
  )
  if (is.null(names(predictable_standard_deviation_raw))) {
    stop("Raw predictable standard deviations must be named by endpoint.")
  }
  missing_raw_sd <- setdiff(
    endpoints,
    names(predictable_standard_deviation_raw)
  )
  if (
    length(missing_raw_sd) > 0L ||
      any(!is.finite(predictable_standard_deviation_raw[endpoints])) ||
      any(predictable_standard_deviation_raw[endpoints] < 0) ||
      length(standard_deviation_floor) != 1L ||
      (!is.na(standard_deviation_floor) &&
        (!is.finite(standard_deviation_floor) ||
          standard_deviation_floor <= 0)) ||
      length(variance_calibration_observations) != 1L ||
      (!is.na(variance_calibration_observations) &&
        variance_calibration_observations < 2L)
  ) {
    stop("Raw variance-calibration diagnostics are invalid.")
  }
  signals$sequence_position <- sequence_position
  signals <- signals[order(signals$sequence_position), , drop = FALSE]
  list(
    chunk = as.integer(chunk),
    training_end = training_end,
    standard_deviation_source_end = standard_deviation_source_end,
    fit_source_end = fit_source_end,
    signals = signals,
    predictable_standard_deviation =
      predictable_standard_deviation[endpoints],
    predictable_standard_deviation_raw =
      predictable_standard_deviation_raw[endpoints],
    standard_deviation_floor = standard_deviation_floor,
    variance_calibration_observations =
      as.integer(variance_calibration_observations)
  )
}

estimate_predictable_online_sd <- function(
    past_one_step_signals,
    mode = "rkhs",
    standard_deviation_floor = 1e-3,
    return_diagnostics = FALSE
) {
  endpoints <- online_bound_endpoint_names()
  values <- as.data.frame(past_one_step_signals, stringsAsFactors = FALSE)
  if (
    !all(c("mode", endpoints) %in% names(values)) ||
      length(standard_deviation_floor) != 1L ||
      !is.finite(standard_deviation_floor) ||
      standard_deviation_floor <= 0 ||
      length(return_diagnostics) != 1L ||
      is.na(return_diagnostics)
  ) {
    stop("Past signals or the standard-deviation floor are invalid.")
  }
  values <- values[values$mode == mode, , drop = FALSE]
  if (
    nrow(values) < 2L ||
      any(!is.finite(as.matrix(values[, endpoints, drop = FALSE])))
  ) {
    stop("Too few finite past signals for the requested mode.")
  }
  raw_standard_deviation <- vapply(
    values[, endpoints, drop = FALSE],
    stats::sd,
    numeric(1)
  )
  standard_deviation <- pmax(
    raw_standard_deviation,
    standard_deviation_floor
  )
  if (!isTRUE(return_diagnostics)) {
    return(standard_deviation)
  }
  list(
    standard_deviation = standard_deviation,
    raw_standard_deviation = raw_standard_deviation,
    floor_active = raw_standard_deviation < standard_deviation_floor,
    standard_deviation_floor = standard_deviation_floor,
    calibration_observations = nrow(values),
    mode = mode
  )
}

make_online_chunk_from_one_step_result <- function(
    result,
    mode = "rkhs",
    evaluation_row_id,
    sequence_position,
    chunk,
    training_end,
    predictable_standard_deviation,
    fit_source_end = training_end,
    standard_deviation_source_end = training_end,
    predictable_standard_deviation_raw =
      predictable_standard_deviation,
    standard_deviation_floor = NA_real_,
    variance_calibration_observations = NA_integer_
) {
  if (
    !is.list(result) ||
      !"one_step_signals" %in% names(result) ||
      !is.data.frame(result$one_step_signals)
  ) {
    stop("result must contain the one_step_signals data frame.")
  }
  endpoints <- online_bound_endpoint_names()
  table <- result$one_step_signals
  if (!all(c("row_id", "mode", endpoints) %in% names(table))) {
    stop("The one-step result omits row identifiers, modes, or endpoints.")
  }
  fit_source_end <- as.integer(fit_source_end)
  if (
    length(fit_source_end) != 1L || is.na(fit_source_end) ||
      fit_source_end > training_end
  ) {
    stop("The fixed-q fit is not measurable using the declared training past.")
  }
  evaluation_row_id <- as.integer(evaluation_row_id)
  local <- table[
    table$mode == mode & table$row_id %in% evaluation_row_id,
    c("row_id", endpoints),
    drop = FALSE
  ]
  if (
    nrow(local) != length(evaluation_row_id) ||
      anyDuplicated(local$row_id) ||
      !setequal(local$row_id, evaluation_row_id)
  ) {
    stop("The requested evaluation rows are absent or duplicated.")
  }
  local <- local[match(evaluation_row_id, local$row_id), , drop = FALSE]
  local$row_id <- NULL
  make_online_signal_chunk(
    chunk = chunk,
    training_end = training_end,
    sequence_position = sequence_position,
    signals = local,
    predictable_standard_deviation = predictable_standard_deviation,
    standard_deviation_source_end = standard_deviation_source_end,
    fit_source_end = fit_source_end,
    predictable_standard_deviation_raw =
      predictable_standard_deviation_raw,
    standard_deviation_floor = standard_deviation_floor,
    variance_calibration_observations =
      variance_calibration_observations
  )
}

validate_online_signal_chunks <- function(chunks) {
  if (!is.list(chunks) || length(chunks) == 0L) {
    stop("chunks must be a nonempty list of online signal chunks.")
  }
  endpoints <- online_bound_endpoint_names()
  positions <- integer()
  previous_end <- 0L
  for (index in seq_along(chunks)) {
    value <- chunks[[index]]
    required <- c(
      "chunk", "training_end", "standard_deviation_source_end",
      "fit_source_end", "signals", "predictable_standard_deviation",
      "predictable_standard_deviation_raw", "standard_deviation_floor",
      "variance_calibration_observations"
    )
    if (!is.list(value) || !all(required %in% names(value))) {
      stop("Online chunk ", index, " is incomplete.")
    }
    local_positions <- value$signals$sequence_position
    if (
      !all(endpoints %in% names(value$signals)) ||
        any(!is.finite(as.matrix(
          value$signals[, endpoints, drop = FALSE]
        ))) ||
        any(local_positions <= value$training_end) ||
        value$standard_deviation_source_end > value$training_end ||
        value$fit_source_end > value$training_end ||
        min(local_positions) <= previous_end ||
        anyDuplicated(c(positions, local_positions))
    ) {
      stop("Online chunk ", index, " is not predictable or ordered.")
    }
    standard_deviation <- value$predictable_standard_deviation[endpoints]
    raw_standard_deviation <-
      value$predictable_standard_deviation_raw[endpoints]
    if (
      any(!is.finite(standard_deviation)) ||
        any(standard_deviation <= 0) ||
        any(!is.finite(raw_standard_deviation)) ||
        any(raw_standard_deviation < 0)
    ) {
      stop("Online chunk ", index, " has invalid variance weights.")
    }
    positions <- c(positions, local_positions)
    previous_end <- max(local_positions)
  }
  invisible(TRUE)
}

summarize_online_signal_chunks <- function(
    chunks,
    population_target,
    endpoints = online_bound_endpoint_names()
) {
  validate_online_signal_chunks(chunks)
  endpoints <- as.character(endpoints)
  if (
    !length(endpoints) || anyDuplicated(endpoints) ||
      !all(endpoints %in% online_bound_endpoint_names())
  ) {
    stop("endpoints must be distinct online bound endpoint names.")
  }
  if (
    is.null(names(population_target)) ||
      !all(endpoints %in% names(population_target)) ||
      any(!is.finite(population_target[endpoints]))
  ) {
    stop("population_target must be a finite named vector for endpoints.")
  }
  denominator <- stats::setNames(vapply(endpoints, function(endpoint) {
    sum(vapply(chunks, function(value) {
      nrow(value$signals) /
        value$predictable_standard_deviation[[endpoint]]
    }, numeric(1)))
  }, numeric(1)), endpoints)
  rows <- lapply(seq_along(chunks), function(index) {
    value <- chunks[[index]]
    observations <- nrow(value$signals)
    do.call(rbind, lapply(endpoints, function(endpoint) {
      signal <- value$signals[[endpoint]]
      predictable_sd <-
        value$predictable_standard_deviation[[endpoint]]
      raw_predictable_sd <-
        value$predictable_standard_deviation_raw[[endpoint]]
      aggregate_weight <- observations / predictable_sd /
        denominator[[endpoint]]
      score_bias <- mean(signal) - population_target[[endpoint]]
      data.frame(
        chunk = value$chunk,
        endpoint = endpoint,
        training_end = value$training_end,
        evaluation_start = min(value$signals$sequence_position),
        evaluation_end = max(value$signals$sequence_position),
        evaluation_observations = observations,
        mean_signal = mean(signal),
        population_target = population_target[[endpoint]],
        score_bias = score_bias,
        observed_signal_sd = if (observations > 1L) {
          stats::sd(signal)
        } else {
          NA_real_
        },
        predictable_standard_deviation = predictable_sd,
        raw_predictable_standard_deviation = raw_predictable_sd,
        variance_floor_active = predictable_sd > raw_predictable_sd,
        variance_calibration_observations =
          value$variance_calibration_observations,
        aggregate_weight = aggregate_weight,
        weighted_score_bias = aggregate_weight * score_bias,
        stringsAsFactors = FALSE
      )
    }))
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

aggregate_online_one_step_bounds <- function(
    chunks,
    confidence_level = 0.95
) {
  validate_online_signal_chunks(chunks)
  if (
    length(confidence_level) != 1L || !is.finite(confidence_level) ||
      confidence_level <= 0 || confidence_level >= 1
  ) {
    stop("confidence_level must lie strictly between zero and one.")
  }
  endpoints <- online_bound_endpoint_names()
  alpha <- 1 - confidence_level
  marginal_one_sided_quantile <- stats::qnorm(confidence_level)
  simultaneous_one_sided_quantile <- stats::qnorm(1 - alpha / 2)
  rows <- lapply(endpoints, function(endpoint) {
    signal <- unlist(lapply(
      chunks,
      function(value) value$signals[[endpoint]]
    ), use.names = FALSE)
    standard_deviation <- unlist(lapply(
      chunks,
      function(value) {
        rep(
          value$predictable_standard_deviation[[endpoint]],
          nrow(value$signals)
        )
      }
    ), use.names = FALSE)
    raw_standard_deviation <- unlist(lapply(
      chunks,
      function(value) {
        rep(
          value$predictable_standard_deviation_raw[[endpoint]],
          nrow(value$signals)
        )
      }
    ), use.names = FALSE)
    floor_active <- unlist(lapply(
      chunks,
      function(value) {
        rep(
          value$predictable_standard_deviation[[endpoint]] >
            value$predictable_standard_deviation_raw[[endpoint]],
          nrow(value$signals)
        )
      }
    ), use.names = FALSE)
    inverse_standard_deviation <- 1 / standard_deviation
    gamma <- mean(inverse_standard_deviation)
    estimate <- sum(inverse_standard_deviation * signal) /
      sum(inverse_standard_deviation)
    standard_error <- 1 / (gamma * sqrt(length(signal)))
    effective_sample_size <- sum(inverse_standard_deviation)^2 /
      sum(inverse_standard_deviation^2)
    standardized_residual <- unlist(lapply(chunks, function(value) {
      local_signal <- value$signals[[endpoint]]
      local_center <- mean(local_signal)
      (local_signal - local_center) /
        value$predictable_standard_deviation[[endpoint]]
    }), use.names = FALSE)
    direction <- if (grepl("_lower$", endpoint)) "lower" else "upper"
    data.frame(
      endpoint = endpoint,
      direction = direction,
      observations = length(signal),
      chunks = length(chunks),
      estimate = estimate,
      Gamma = gamma,
      standard_error = standard_error,
      effective_sample_size = effective_sample_size,
      minimum_predictable_standard_deviation = min(standard_deviation),
      maximum_predictable_standard_deviation = max(standard_deviation),
      minimum_raw_standard_deviation = min(raw_standard_deviation),
      variance_floor_active_rate = mean(floor_active),
      maximum_scaled_standardized_residual =
        max(abs(standardized_residual)) / sqrt(length(signal)),
      one_sided_outward_confidence_limit = if (direction == "lower") {
        estimate - marginal_one_sided_quantile * standard_error
      } else {
        estimate + marginal_one_sided_quantile * standard_error
      },
      outward_confidence_limit = if (direction == "lower") {
        estimate - simultaneous_one_sided_quantile * standard_error
      } else {
        estimate + simultaneous_one_sided_quantile * standard_error
      },
      confidence_level = confidence_level,
      marginal_one_sided_quantile = marginal_one_sided_quantile,
      quantile = simultaneous_one_sided_quantile,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  ate_lower <- summary[summary$endpoint == "ate_lower", , drop = FALSE]
  ate_upper <- summary[summary$endpoint == "ate_upper", , drop = FALSE]
  identified_set <- data.frame(
    confidence_level = confidence_level,
    lower = ate_lower$outward_confidence_limit,
    upper = ate_upper$outward_confidence_limit,
    ordered = ate_lower$outward_confidence_limit <=
      ate_upper$outward_confidence_limit,
    interpretation = paste(
      "Outward one-sided online confidence set; validity requires",
      "predictable variance consistency and averaged nuisance drift."
    ),
    stringsAsFactors = FALSE
  )
  signal_table <- do.call(rbind, lapply(chunks, function(value) {
    local <- value$signals
    local$chunk <- value$chunk
    local$training_end <- value$training_end
    local$standard_deviation_source_end <-
      value$standard_deviation_source_end
    local$fit_source_end <- value$fit_source_end
    local$variance_calibration_observations <-
      value$variance_calibration_observations
    local
  }))
  rownames(signal_table) <- NULL
  list(
    summary = summary,
    ate_outward_confidence_set = identified_set,
    signals = signal_table,
    specification = list(
      construction = "inverse-predictable-SD online one-step",
      confidence_level = confidence_level,
      multiplicity = "Bonferroni across the lower and upper ATE endpoints",
      target = "data-adaptive outer bounds"
    )
  )
}
