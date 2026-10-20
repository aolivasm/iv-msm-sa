# Oracle-allocation diagnostics for the learned-witness compatibility test.
#
# This file does not change the production test.  It supplies an alternative
# fold routine used only by the focused calibrated-DGP decomposition.  The
# routine can replace the estimated additive-spline conditional allocation
# means by their DGP-known values and can optionally localize the pointwise
# inequalities in X.  Witness construction and all selection remain honest:
# cutoffs, scales, and localization rules are learned only on the outer
# training fold and are then frozen on validation.

check_oracle_witness_diagnostic_dependencies <- function() {
  required <- c(
    "compute_learned_witness_scores", "make_witness_cutoffs",
    "fit_witness_allocation_means", "compute_dual_msm_moment_primitives",
    "compute_learned_family_base", "make_arm_lambda_signal_blocks",
    "observed_dgp_nuisances", "observed_dgp_outcome_location",
    "observed_dgp_tilt", "observed_dgp_gamma",
    "observed_binary_probabilities", "observed_binary_gamma"
  )
  missing <- required[
    !vapply(required, exists, logical(1L), mode = "function")
  ]
  if (length(missing)) {
    stop(
      "Source the learned-witness test and observed-data DGP files first. ",
      "Missing: ", paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

oracle_witness_gauss_legendre_rule <- function(nodes = 21L) {
  nodes <- as.integer(nodes)
  if (length(nodes) != 1L || is.na(nodes) || nodes < 3L) {
    stop("nodes must be one integer of at least three.")
  }
  index <- seq_len(nodes - 1L)
  off_diagonal <- index / sqrt(4 * index^2 - 1)
  jacobi <- matrix(0, nodes, nodes)
  jacobi[cbind(index, index + 1L)] <- off_diagonal
  jacobi[cbind(index + 1L, index)] <- off_diagonal
  decomposition <- eigen(jacobi, symmetric = TRUE)
  ordering <- order(decomposition$values)
  list(
    nodes = decomposition$values[ordering],
    weights = 2 * decomposition$vectors[1L, ordering]^2
  )
}

make_oracle_witness_response_specifications <- function(
    score_training,
    training,
    treatment,
    witness_tail_fractions,
    witness_shapes
) {
  rows <- list()
  counter <- 0L
  for (family in names(score_training)) {
    family_counter <- counter
    score <- score_training[[family]]
    active <- training$T_obs == treatment & is.finite(score)
    cutoffs <- make_witness_cutoffs(
      score,
      active = active,
      upper_tail_fractions = witness_tail_fractions
    )
    if ("indicator" %in% witness_shapes) {
      for (cutoff_name in names(cutoffs)) {
        counter <- counter + 1L
        rows[[counter]] <- data.frame(
          key = paste(family, "indicator", cutoff_name, sep = "__"),
          family = family,
          shape = "indicator",
          cutoff = unname(cutoffs[cutoff_name]),
          scale = NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
    positive <- pmax(score, 0)
    scale <- max(positive[active], na.rm = TRUE)
    if (is.finite(scale) && scale > 1e-10) {
      for (shape in intersect(witness_shapes, c("linear", "sqrt"))) {
        counter <- counter + 1L
        rows[[counter]] <- data.frame(
          key = paste(family, shape, "positive", sep = "__"),
          family = family,
          shape = shape,
          cutoff = 0,
          scale = scale,
          stringsAsFactors = FALSE
        )
      }
    }
    if (counter == family_counter) {
      counter <- counter + 1L
      rows[[counter]] <- data.frame(
        key = paste(family, "fallback", "zero", sep = "__"),
        family = family,
        shape = "indicator",
        cutoff = 0,
        scale = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

evaluate_oracle_witness_responses <- function(scores, specifications) {
  output <- matrix(
    NA_real_, nrow(scores), nrow(specifications),
    dimnames = list(NULL, specifications$key)
  )
  for (index in seq_len(nrow(specifications))) {
    specification <- specifications[index, , drop = FALSE]
    score <- scores[[specification$family]]
    output[, index] <- switch(
      specification$shape,
      indicator = as.numeric(score > specification$cutoff),
      linear = pmin(pmax(score, 0) / specification$scale, 1),
      sqrt = pmin(sqrt(pmax(score, 0) / specification$scale), 1),
      stop("Unknown witness-response shape: ", specification$shape, ".")
    )
  }
  storage.mode(output) <- "double"
  output
}

oracle_witness_response_at_y <- function(
    newdata,
    y,
    treatment,
    Lambda,
    specification,
    response_specifications,
    probability
) {
  x <- newdata[, observed_dgp_x_columns(specification$dimension), drop = FALSE]
  gamma <- if (identical(specification$outcome_type, "binary")) {
    observed_binary_gamma(y, x, treatment, specification)
  } else {
    observed_dgp_gamma(y, x, treatment, specification)
  }
  scores <- compute_learned_witness_scores(
    gamma,
    probability$pi0,
    probability$pi1,
    treatment,
    Lambda
  )
  evaluate_oracle_witness_responses(scores, response_specifications)
}

oracle_witness_indicator_probabilities <- function(
    newdata,
    treatment,
    z_value,
    Lambda,
    specification,
    response_specifications,
    probability
) {
  indicator <- which(response_specifications$shape == "indicator")
  if (!length(indicator)) {
    return(matrix(numeric(), nrow(newdata), 0L))
  }
  x <- newdata[, observed_dgp_x_columns(specification$dimension), drop = FALSE]
  n <- nrow(newdata)
  lower <- specification$residual_lower
  upper <- specification$residual_upper
  location <- observed_dgp_outcome_location(x, treatment, specification)
  tilt <- observed_dgp_tilt(x, treatment, specification)
  residual_mean <- if (z_value == 0L) rep(0, n) else tilt
  normalizer <- stats::pnorm(upper, mean = residual_mean) -
    stats::pnorm(lower, mean = residual_mean)
  p0 <- if (treatment == 1L) probability$pi0 else 1 - probability$pi0
  p1 <- if (treatment == 1L) probability$pi1 else 1 - probability$pi1
  primitive <- compute_dual_msm_moment_primitives(p0, p1, Lambda, Lambda)
  gamma_mid <- observed_dgp_gamma(
    location, x, treatment, specification
  )
  gamma_lower_residual <- observed_dgp_gamma(
    location + lower, x, treatment, specification
  )
  gamma_upper_residual <- observed_dgp_gamma(
    location + upper, x, treatment, specification
  )
  family_coefficients <- list(
    pointwise_lower = list(A = primitive$gamma_lower, B = rep(-1, n)),
    pointwise_upper = list(A = -primitive$gamma_upper, B = rep(1, n)),
    normalization_lower = list(
      A = primitive$b - 1 / Lambda,
      B = primitive$c_lower
    ),
    normalization_upper = list(
      A = Lambda - primitive$b,
      B = -primitive$c_upper
    )
  )
  output <- matrix(
    NA_real_, n, length(indicator),
    dimnames = list(NULL, response_specifications$key[indicator])
  )
  for (output_index in seq_along(indicator)) {
    response_index <- indicator[output_index]
    response_specification <-
      response_specifications[response_index, , drop = FALSE]
    coefficient <- family_coefficients[[response_specification$family]]
    score_lower <- coefficient$A + coefficient$B * gamma_lower_residual
    score_upper <- coefficient$A + coefficient$B * gamma_upper_residual
    cutoff <- response_specification$cutoff
    minimum <- pmin(score_lower, score_upper)
    maximum <- pmax(score_lower, score_upper)
    value <- numeric(n)
    value[cutoff < minimum] <- 1
    interior <- cutoff >= minimum & cutoff < maximum
    if (any(interior)) {
      ratio <- (cutoff - coefficient$A[interior]) /
        (coefficient$B[interior] * gamma_mid[interior])
      if (any(!is.finite(ratio) | ratio <= 0)) {
        stop("Could not invert an interior continuous indicator witness.")
      }
      threshold <- log(ratio) / tilt[interior]
      conditional_cdf <- (
        stats::pnorm(threshold, mean = residual_mean[interior]) -
          stats::pnorm(lower, mean = residual_mean[interior])
      ) / normalizer[interior]
      conditional_cdf <- pmin(1, pmax(0, conditional_cdf))
      increasing <- coefficient$B[interior] * tilt[interior] > 0
      value[interior] <- ifelse(
        increasing, 1 - conditional_cdf, conditional_cdf
      )
    }
    output[, output_index] <- value
  }
  output
}

oracle_witness_conditional_allocation_means <- function(
    newdata,
    treatment,
    z_value,
    Lambda,
    specification,
    response_specifications,
    quadrature_nodes = 21L
) {
  x_columns <- observed_dgp_x_columns(specification$dimension)
  x <- newdata[, x_columns, drop = FALSE]
  probability <- observed_dgp_nuisances(x, specification)
  n <- nrow(newdata)

  if (identical(specification$outcome_type, "binary")) {
    outcome <- observed_binary_probabilities(x, treatment, specification)
    event_probability <- if (z_value == 0L) outcome$a0 else outcome$a1
    response_zero <- oracle_witness_response_at_y(
      newdata, rep(0, n), treatment, Lambda, specification,
      response_specifications, probability
    )
    response_one <- oracle_witness_response_at_y(
      newdata, rep(1, n), treatment, Lambda, specification,
      response_specifications, probability
    )
    output <- response_zero * (1 - event_probability) +
      response_one * event_probability
  } else {
    rule <- oracle_witness_gauss_legendre_rule(quadrature_nodes)
    lower <- specification$residual_lower
    upper <- specification$residual_upper
    residual <- (upper - lower) * rule$nodes / 2 + (upper + lower) / 2
    integration_weight <- (upper - lower) * rule$weights / 2
    location <- observed_dgp_outcome_location(x, treatment, specification)
    tilt <- observed_dgp_tilt(x, treatment, specification)
    residual_mean <- if (z_value == 0L) rep(0, n) else tilt
    normalizer <- stats::pnorm(upper, mean = residual_mean) -
      stats::pnorm(lower, mean = residual_mean)
    output <- matrix(
      0, n, nrow(response_specifications),
      dimnames = list(NULL, response_specifications$key)
    )
    indicator <- which(response_specifications$shape == "indicator")
    if (length(indicator)) {
      output[, indicator] <- oracle_witness_indicator_probabilities(
        newdata,
        treatment,
        z_value,
        Lambda,
        specification,
        response_specifications,
        probability
      )
    }
    smooth <- which(response_specifications$shape != "indicator")
    if (length(smooth)) {
      for (node_index in seq_along(residual)) {
        density_weight <- integration_weight[node_index] *
          stats::dnorm(residual[node_index], mean = residual_mean) / normalizer
        response <- oracle_witness_response_at_y(
          newdata,
          location + residual[node_index],
          treatment,
          Lambda,
          specification,
          response_specifications[smooth, , drop = FALSE],
          probability
        )
        output[, smooth] <- output[, smooth, drop = FALSE] +
          response * density_weight
      }
    }
  }
  output[] <- pmin(1 - 1e-8, pmax(1e-8, output))
  output
}

fit_oracle_witness_allocation_means <- function(
    training,
    validation,
    treatment,
    Lambda,
    specification,
    response_specifications,
    quadrature_nodes = 21L
) {
  output <- list()
  for (z_value in 0:1) {
    output[[paste0("train_z", z_value)]] <-
      oracle_witness_conditional_allocation_means(
        training, treatment, z_value, Lambda, specification,
        response_specifications, quadrature_nodes
      )
    output[[paste0("validation_z", z_value)]] <-
      oracle_witness_conditional_allocation_means(
        validation, treatment, z_value, Lambda, specification,
        response_specifications, quadrature_nodes
      )
  }
  output
}

make_oracle_diagnostic_localization_cutoffs <- function(
    base,
    family,
    upper_tail_fractions,
    localize_pointwise
) {
  if (!localize_pointwise && !grepl("^normalization_", family)) {
    return(c(all = -Inf))
  }
  local <- base[is.finite(base)]
  cutoffs <- c(
    positive = 0,
    stats::setNames(
      as.numeric(stats::quantile(
        local,
        probs = 1 - upper_tail_fractions,
        type = 8,
        names = FALSE
      )),
      paste0("top", format(100 * upper_tail_fractions, trim = TRUE))
    ),
    all = -Inf
  )
  cutoffs[!duplicated(round(cutoffs, 12))]
}

select_oracle_diagnostic_witness <- function(
    family,
    response_training,
    response_validation,
    mean_fit,
    primitive_training,
    primitive_validation,
    Lambda,
    localization_tail_fractions,
    localize_pointwise
) {
  candidates <- list()
  counter <- 0L
  for (column in seq_len(ncol(response_training))) {
    base_training <- compute_learned_family_base(
      family,
      primitive_training,
      mean_fit$train_z0[, column],
      mean_fit$train_z1[, column],
      Lambda
    )
    base_validation <- compute_learned_family_base(
      family,
      primitive_validation,
      mean_fit$validation_z0[, column],
      mean_fit$validation_z1[, column],
      Lambda
    )
    h_cutoffs <- make_oracle_diagnostic_localization_cutoffs(
      base_training,
      family,
      localization_tail_fractions,
      localize_pointwise
    )
    for (h_name in names(h_cutoffs)) {
      h_training <- as.numeric(base_training > h_cutoffs[h_name])
      active_fraction <- mean(h_training)
      if (active_fraction < 0.025) next
      plug_in <- h_training * base_training
      standard_deviation <- stats::sd(plug_in)
      criterion <- if (
        is.finite(standard_deviation) && standard_deviation > 1e-10
      ) {
        sqrt(length(plug_in)) * mean(plug_in) / standard_deviation
      } else {
        -Inf
      }
      counter <- counter + 1L
      candidates[[counter]] <- list(
        column = column,
        r_name = colnames(response_training)[column],
        h_name = h_name,
        h_cutoff = unname(h_cutoffs[h_name]),
        criterion = criterion,
        h_validation = as.numeric(base_validation > h_cutoffs[h_name]),
        base_training_mean = mean(plug_in),
        base_training_standard_deviation = standard_deviation,
        h_training_fraction = active_fraction
      )
    }
  }
  if (!length(candidates)) stop("No diagnostic witness candidate remained.")
  criterion <- vapply(candidates, `[[`, numeric(1L), "criterion")
  selected <- candidates[[which.max(criterion)]]
  selected$R_validation <- response_validation[, selected$column, drop = FALSE]
  selected$m0_validation <-
    mean_fit$validation_z0[, selected$column, drop = FALSE]
  selected$m1_validation <-
    mean_fit$validation_z1[, selected$column, drop = FALSE]
  selected$H_validation <- matrix(
    selected$h_validation,
    ncol = 1L,
    dimnames = list(NULL, paste0("learned_", family))
  )
  selected
}

# Drop-in replacement for fit_fold_learned_witness_signals().  The calling
# simulation must define .oracle_witness_diagnostic_control in the global
# environment with allocation_method, localize_pointwise, specification, and
# quadrature_nodes fields.
fit_fold_oracle_witness_diagnostic_signals <- function(
    training,
    validation,
    probability_training,
    probability_validation,
    gamma_training,
    gamma_validation,
    treatment,
    Lambda,
    x_cols,
    allocation_spline_df,
    witness_tail_fractions,
    localization_tail_fractions,
    witness_shapes
) {
  control <- get(
    ".oracle_witness_diagnostic_control",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  if (!control$allocation_method %in% c("additive_series", "oracle")) {
    stop("Unknown diagnostic allocation method.")
  }
  score_training <- compute_learned_witness_scores(
    gamma_training,
    probability_training$pi0,
    probability_training$pi1,
    treatment,
    Lambda
  )
  score_validation <- compute_learned_witness_scores(
    gamma_validation,
    probability_validation$pi0,
    probability_validation$pi1,
    treatment,
    Lambda
  )
  response_specifications <- make_oracle_witness_response_specifications(
    score_training,
    training,
    treatment,
    witness_tail_fractions,
    witness_shapes
  )
  response_training <- evaluate_oracle_witness_responses(
    score_training, response_specifications
  )
  response_validation <- evaluate_oracle_witness_responses(
    score_validation, response_specifications
  )
  response_family <- stats::setNames(
    response_specifications$family,
    response_specifications$key
  )

  mean_fit <- if (control$allocation_method == "oracle") {
    fit_oracle_witness_allocation_means(
      training,
      validation,
      treatment,
      Lambda,
      control$specification,
      response_specifications,
      control$quadrature_nodes
    )
  } else {
    fit_witness_allocation_means(
      training,
      validation,
      response_training,
      treatment,
      x_cols,
      allocation_spline_df
    )
  }

  p0_training <- if (treatment == 1L) {
    probability_training$pi0
  } else 1 - probability_training$pi0
  p1_training <- if (treatment == 1L) {
    probability_training$pi1
  } else 1 - probability_training$pi1
  p0_validation <- if (treatment == 1L) {
    probability_validation$pi0
  } else 1 - probability_validation$pi0
  p1_validation <- if (treatment == 1L) {
    probability_validation$pi1
  } else 1 - probability_validation$pi1
  primitive_training <- compute_dual_msm_moment_primitives(
    p0_training, p1_training, Lambda, Lambda
  )
  primitive_validation <- compute_dual_msm_moment_primitives(
    p0_validation, p1_validation, Lambda, Lambda
  )

  output_signals <- output_metadata <- diagnostics <- list()
  for (family in names(score_training)) {
    columns <- which(response_family == family)
    local_mean_fit <- lapply(mean_fit, function(value) {
      value[, columns, drop = FALSE]
    })
    selected <- select_oracle_diagnostic_witness(
      family,
      response_training[, columns, drop = FALSE],
      response_validation[, columns, drop = FALSE],
      local_mean_fit,
      primitive_training,
      primitive_validation,
      Lambda,
      localization_tail_fractions,
      control$localize_pointwise
    )
    block <- make_arm_lambda_signal_blocks(
      data = validation,
      H = selected$H_validation,
      R = selected$R_validation,
      m0 = selected$m0_validation,
      m1 = selected$m1_validation,
      e = probability_validation$e,
      pi0 = probability_validation$pi0,
      pi1 = probability_validation$pi1,
      t_value = treatment,
      Lambda0 = Lambda,
      Lambda1 = Lambda
    )
    keep <- block$metadata$family == family
    if (sum(keep) != 1L) stop("A diagnostic witness did not map to one moment.")
    output_signals[[family]] <- block$signals[, keep, drop = FALSE]
    local_metadata <- block$metadata[keep, , drop = FALSE]
    local_metadata$r <- selected$r_name
    local_metadata$h <- paste0("learned_", selected$h_name)
    output_metadata[[family]] <- local_metadata
    diagnostics[[family]] <- data.frame(
      treatment = treatment,
      Lambda = Lambda,
      family = family,
      selected_r = selected$r_name,
      selected_h = selected$h_name,
      selection_criterion = selected$criterion,
      training_plugin_mean = selected$base_training_mean,
      training_plugin_standard_deviation =
        selected$base_training_standard_deviation,
      training_h_fraction = selected$h_training_fraction,
      validation_r_fraction = mean(selected$R_validation),
      validation_h_fraction = mean(selected$H_validation),
      diagnostic_allocation_method = control$allocation_method,
      diagnostic_localize_pointwise = control$localize_pointwise,
      stringsAsFactors = FALSE
    )
  }
  list(
    signals = do.call(cbind, output_signals),
    metadata = do.call(rbind, output_metadata),
    diagnostics = do.call(rbind, diagnostics)
  )
}
