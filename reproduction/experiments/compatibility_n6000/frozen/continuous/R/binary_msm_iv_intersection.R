# Sharp binary-outcome bounds under an ordinary marginal sensitivity model
# intersected with either a valid instrumental-variable restriction or the
# stronger treatment-proxy bridge restriction obtained by using Z as a proxy.
#
# The MSM is imposed after marginalizing over the instrument.  This differs
# from the stratum-specific dual IV-MSM implemented in binary_outcome_bounds.R.
# For a binary outcome, the counterfactual law is characterized by
# h_t(x) = P(Y^t = 1 | X = x).  The estimator intersects the exact MSM interval
# for h_t(x) with the exact IV interval and evaluates the resulting endpoints.
#
# The implementation accepts cross-fitted predictions of
#   e(x)       = P(Z = 0 | X = x),
#   pi_z(x)    = P(T = 1 | X = x, Z = z), and
#   a_tz(x)    = P(Y = 1 | T = t, X = x, Z = z).
# It reports plug-in bounds, cross-fitted one-step estimates, outward Wald
# confidence sets, active-branch diagnostics, tie rates, and estimated
# crossing rates.  Crossed envelopes are flagged and never repaired.  These
# crossing diagnostics are not formal compatibility tests.

msm_iv_binary_required_prediction_names <- function() {
  c(
    "e", "pi0", "pi1",
    "a_t0_z0", "a_t0_z1", "a_t1_z0", "a_t1_z1"
  )
}

msm_iv_binary_validate_inputs <- function(data, predictions) {
  required_data <- c("Y", "T_obs", "Z")
  missing_data <- setdiff(required_data, names(data))
  if (length(missing_data)) {
    stop("Missing data columns: ", paste(missing_data, collapse = ", "))
  }
  required_predictions <- msm_iv_binary_required_prediction_names()
  missing_predictions <- setdiff(required_predictions, names(predictions))
  if (length(missing_predictions)) {
    stop(
      "Missing prediction columns: ",
      paste(missing_predictions, collapse = ", ")
    )
  }
  if (nrow(data) != nrow(predictions)) {
    stop("data and predictions must have the same number of rows.")
  }
  for (name in required_data) {
    observed <- sort(unique(stats::na.omit(data[[name]])))
    if (!length(observed) || !all(observed %in% c(0, 1))) {
      stop(name, " must be coded 0/1.")
    }
  }
  prediction_matrix <- as.matrix(predictions[, required_predictions])
  if (any(!is.finite(prediction_matrix))) {
    stop("All nuisance predictions must be finite.")
  }
  if (any(prediction_matrix <= 0 | prediction_matrix >= 1)) {
    stop("All nuisance predictions must lie strictly between zero and one.")
  }
  invisible(TRUE)
}

msm_iv_binary_arm_nuisance <- function(predictions, t_value) {
  if (!t_value %in% c(0, 1)) stop("t_value must be zero or one.")
  if (t_value == 1) {
    return(list(
      e = predictions$e,
      r0 = predictions$pi0,
      r1 = predictions$pi1,
      a0 = predictions$a_t1_z0,
      a1 = predictions$a_t1_z1
    ))
  }
  list(
    e = predictions$e,
    r0 = 1 - predictions$pi0,
    r1 = 1 - predictions$pi1,
    a0 = predictions$a_t0_z0,
    a1 = predictions$a_t0_z1
  )
}

msm_iv_binary_zero_gradient <- function(n) {
  matrix(
    0,
    nrow = n,
    ncol = 5L,
    dimnames = list(NULL, c("e", "r0", "r1", "a0", "a1"))
  )
}

msm_iv_binary_choose <- function(
    first,
    second,
    gradient_first,
    gradient_second,
    operation = c("max", "min"),
    tie_tolerance = 1e-8
) {
  operation <- match.arg(operation)
  if (operation == "max") {
    first_active <- first > second + tie_tolerance
    second_active <- second > first + tie_tolerance
    value <- pmax(first, second)
  } else {
    first_active <- first < second - tie_tolerance
    second_active <- second < first - tie_tolerance
    value <- pmin(first, second)
  }
  tie <- !(first_active | second_active)
  weight_first <- as.numeric(first_active) + 0.5 * as.numeric(tie)
  weight_second <- as.numeric(second_active) + 0.5 * as.numeric(tie)
  gradient <- gradient_first * weight_first + gradient_second * weight_second
  list(
    value = value,
    gradient = gradient,
    tie = tie,
    first_active = first_active,
    second_active = second_active,
    gap = abs(first - second)
  )
}

msm_iv_binary_arm_endpoint <- function(
    predictions,
    t_value,
    Lambda,
    tie_tolerance = 1e-8
) {
  if (length(Lambda) != 1L || !is.finite(Lambda) || Lambda < 1) {
    stop("Lambda must be one finite number of at least one.")
  }
  if (
    length(tie_tolerance) != 1L || !is.finite(tie_tolerance) ||
      tie_tolerance < 0
  ) {
    stop("tie_tolerance must be one nonnegative finite number.")
  }

  nuisance <- msm_iv_binary_arm_nuisance(predictions, t_value)
  e <- nuisance$e
  r0 <- nuisance$r0
  r1 <- nuisance$r1
  a0 <- nuisance$a0
  a1 <- nuisance$a1
  n_observations <- length(e)
  zero <- msm_iv_binary_zero_gradient(n_observations)

  # Marginalize over Z to obtain the ordinary-MSM nuisances.
  r <- e * r0 + (1 - e) * r1
  observed_joint_one <- e * r0 * a0 + (1 - e) * r1 * a1
  a <- observed_joint_one / r

  gradient_r <- zero
  gradient_r[, "e"] <- r0 - r1
  gradient_r[, "r0"] <- e
  gradient_r[, "r1"] <- 1 - e

  gradient_joint <- zero
  gradient_joint[, "e"] <- r0 * a0 - r1 * a1
  gradient_joint[, "r0"] <- e * a0
  gradient_joint[, "r1"] <- (1 - e) * a1
  gradient_joint[, "a0"] <- e * r0
  gradient_joint[, "a1"] <- (1 - e) * r1
  gradient_a <- (gradient_joint - gradient_r * a) / r

  # The ordinary MSM bounds the missing-arm Bernoulli probability theta_t.
  theta_lower_first <- a / Lambda
  theta_lower_second <- 1 - Lambda * (1 - a)
  gradient_theta_lower_first <- gradient_a / Lambda
  gradient_theta_lower_second <- gradient_a * Lambda
  theta_lower <- msm_iv_binary_choose(
    theta_lower_first,
    theta_lower_second,
    gradient_theta_lower_first,
    gradient_theta_lower_second,
    operation = "max",
    tie_tolerance = tie_tolerance
  )

  theta_upper_first <- Lambda * a
  theta_upper_second <- 1 - (1 - a) / Lambda
  gradient_theta_upper_first <- gradient_a * Lambda
  gradient_theta_upper_second <- gradient_a / Lambda
  theta_upper <- msm_iv_binary_choose(
    theta_upper_first,
    theta_upper_second,
    gradient_theta_upper_first,
    gradient_theta_upper_second,
    operation = "min",
    tie_tolerance = tie_tolerance
  )

  msm_lower <- observed_joint_one + (1 - r) * theta_lower$value
  msm_upper <- observed_joint_one + (1 - r) * theta_upper$value
  gradient_msm_lower <- gradient_joint -
    gradient_r * theta_lower$value +
    theta_lower$gradient * (1 - r)
  gradient_msm_upper <- gradient_joint -
    gradient_r * theta_upper$value +
    theta_upper$gradient * (1 - r)

  # IV-only binary-outcome bounds on h_t(x) = P(Y^t = 1 | X = x).
  iv_lower_zero <- r0 * a0
  iv_lower_one <- r1 * a1
  gradient_iv_lower_zero <- zero
  gradient_iv_lower_zero[, "r0"] <- a0
  gradient_iv_lower_zero[, "a0"] <- r0
  gradient_iv_lower_one <- zero
  gradient_iv_lower_one[, "r1"] <- a1
  gradient_iv_lower_one[, "a1"] <- r1
  iv_lower <- msm_iv_binary_choose(
    iv_lower_zero,
    iv_lower_one,
    gradient_iv_lower_zero,
    gradient_iv_lower_one,
    operation = "max",
    tie_tolerance = tie_tolerance
  )

  iv_upper_zero <- 1 - r0 * (1 - a0)
  iv_upper_one <- 1 - r1 * (1 - a1)
  gradient_iv_upper_zero <- zero
  gradient_iv_upper_zero[, "r0"] <- a0 - 1
  gradient_iv_upper_zero[, "a0"] <- r0
  gradient_iv_upper_one <- zero
  gradient_iv_upper_one[, "r1"] <- a1 - 1
  gradient_iv_upper_one[, "a1"] <- r1
  iv_upper <- msm_iv_binary_choose(
    iv_upper_zero,
    iv_upper_one,
    gradient_iv_upper_zero,
    gradient_iv_upper_one,
    operation = "min",
    tie_tolerance = tie_tolerance
  )

  intersection_lower <- msm_iv_binary_choose(
    msm_lower,
    iv_lower$value,
    gradient_msm_lower,
    iv_lower$gradient,
    operation = "max",
    tie_tolerance = tie_tolerance
  )
  intersection_upper <- msm_iv_binary_choose(
    msm_upper,
    iv_upper$value,
    gradient_msm_upper,
    iv_upper$gradient,
    operation = "min",
    tie_tolerance = tie_tolerance
  )

  # Treatment-proxy bridge restriction.  For a binary outcome, the proxy
  # density-ratio envelope is equivalent to requiring the missing-arm risk
  # theta_t(x) to lie in the convex hull of the two observed Z-specific risks:
  #   min_z a_tz(x) <= theta_t(x) <= max_z a_tz(x).
  # This uses more than IV validity: it also requires existence of the
  # relevant nonnegative treatment-confounding bridge.
  gradient_a0 <- zero
  gradient_a0[, "a0"] <- 1
  gradient_a1 <- zero
  gradient_a1[, "a1"] <- 1
  proxy_theta_lower <- msm_iv_binary_choose(
    a0,
    a1,
    gradient_a0,
    gradient_a1,
    operation = "min",
    tie_tolerance = tie_tolerance
  )
  proxy_theta_upper <- msm_iv_binary_choose(
    a0,
    a1,
    gradient_a0,
    gradient_a1,
    operation = "max",
    tie_tolerance = tie_tolerance
  )
  treatment_proxy_theta_lower <- msm_iv_binary_choose(
    theta_lower$value,
    proxy_theta_lower$value,
    theta_lower$gradient,
    proxy_theta_lower$gradient,
    operation = "max",
    tie_tolerance = tie_tolerance
  )
  treatment_proxy_theta_upper <- msm_iv_binary_choose(
    theta_upper$value,
    proxy_theta_upper$value,
    theta_upper$gradient,
    proxy_theta_upper$gradient,
    operation = "min",
    tie_tolerance = tie_tolerance
  )
  treatment_proxy_lower <- observed_joint_one +
    (1 - r) * treatment_proxy_theta_lower$value
  treatment_proxy_upper <- observed_joint_one +
    (1 - r) * treatment_proxy_theta_upper$value
  gradient_treatment_proxy_lower <- gradient_joint -
    gradient_r * treatment_proxy_theta_lower$value +
    treatment_proxy_theta_lower$gradient * (1 - r)
  gradient_treatment_proxy_upper <- gradient_joint -
    gradient_r * treatment_proxy_theta_upper$value +
    treatment_proxy_theta_upper$gradient * (1 - r)

  methods <- list(
    msm_only = list(
      lower = msm_lower,
      upper = msm_upper,
      gradient_lower = gradient_msm_lower,
      gradient_upper = gradient_msm_upper
    ),
    iv_only = list(
      lower = iv_lower$value,
      upper = iv_upper$value,
      gradient_lower = iv_lower$gradient,
      gradient_upper = iv_upper$gradient
    ),
    msm_iv_intersection = list(
      lower = intersection_lower$value,
      upper = intersection_upper$value,
      gradient_lower = intersection_lower$gradient,
      gradient_upper = intersection_upper$gradient
    ),
    msm_treatment_proxy_intersection = list(
      lower = treatment_proxy_lower,
      upper = treatment_proxy_upper,
      gradient_lower = gradient_treatment_proxy_lower,
      gradient_upper = gradient_treatment_proxy_upper
    )
  )

  list(
    treatment = t_value,
    Lambda = Lambda,
    nuisance = nuisance,
    marginal = list(r = r, a = a),
    methods = methods,
    candidates = list(
      msm_lower = msm_lower,
      msm_upper = msm_upper,
      iv_lower = iv_lower$value,
      iv_upper = iv_upper$value,
      proxy_theta_lower = proxy_theta_lower$value,
      proxy_theta_upper = proxy_theta_upper$value,
      treatment_proxy_lower = treatment_proxy_lower,
      treatment_proxy_upper = treatment_proxy_upper
    ),
    diagnostics = data.frame(
      treatment = t_value,
      Lambda = Lambda,
      msm_lower_branch_tie_rate = mean(theta_lower$tie),
      msm_upper_branch_tie_rate = mean(theta_upper$tie),
      iv_lower_instrument_tie_rate = mean(iv_lower$tie),
      iv_upper_instrument_tie_rate = mean(iv_upper$tie),
      intersection_lower_tie_rate = mean(intersection_lower$tie),
      intersection_upper_tie_rate = mean(intersection_upper$tie),
      iv_active_lower_rate = mean(
        iv_lower$value > msm_lower + tie_tolerance
      ),
      iv_active_upper_rate = mean(
        iv_upper$value < msm_upper - tie_tolerance
      ),
      incompatibility_rate = mean(
        intersection_lower$value >
          intersection_upper$value + tie_tolerance
      ),
      maximum_crossing = max(
        intersection_lower$value - intersection_upper$value
      ),
      minimum_intersection_width = min(
        intersection_upper$value - intersection_lower$value
      ),
      mean_intersection_width = mean(
        intersection_upper$value - intersection_lower$value
      ),
      proxy_lower_instrument_tie_rate = mean(proxy_theta_lower$tie),
      proxy_upper_instrument_tie_rate = mean(proxy_theta_upper$tie),
      treatment_proxy_lower_tie_rate = mean(
        treatment_proxy_theta_lower$tie
      ),
      treatment_proxy_upper_tie_rate = mean(
        treatment_proxy_theta_upper$tie
      ),
      proxy_active_lower_rate = mean(
        proxy_theta_lower$value > theta_lower$value + tie_tolerance
      ),
      proxy_active_upper_rate = mean(
        proxy_theta_upper$value < theta_upper$value - tie_tolerance
      ),
      treatment_proxy_crossing_rate = mean(
        treatment_proxy_lower >
          treatment_proxy_upper + tie_tolerance
      ),
      treatment_proxy_maximum_crossing = max(
        treatment_proxy_lower - treatment_proxy_upper
      ),
      treatment_proxy_minimum_width = min(
        treatment_proxy_upper - treatment_proxy_lower
      ),
      treatment_proxy_mean_width = mean(
        treatment_proxy_upper - treatment_proxy_lower
      ),
      stringsAsFactors = FALSE
    )
  )
}

msm_iv_binary_one_step_signal <- function(
    data,
    arm_endpoint,
    method = c(
      "msm_iv_intersection", "msm_treatment_proxy_intersection",
      "msm_only", "iv_only"
    ),
    direction = c("lower", "upper")
) {
  method <- match.arg(method)
  direction <- match.arg(direction)
  nuisance <- arm_endpoint$nuisance
  component <- arm_endpoint$methods[[method]]
  value <- component[[direction]]
  gradient <- component[[paste0("gradient_", direction)]]
  treatment <- arm_endpoint$treatment
  observed_arm <- as.numeric(data$T_obs == treatment)
  z_zero <- as.numeric(data$Z == 0)
  z_one <- 1 - z_zero

  signal <- value + gradient[, "e"] * (z_zero - nuisance$e)
  signal <- signal + gradient[, "r0"] * z_zero / nuisance$e *
    (observed_arm - nuisance$r0)
  signal <- signal + gradient[, "r1"] * z_one / (1 - nuisance$e) *
    (observed_arm - nuisance$r1)
  signal <- signal + gradient[, "a0"] * z_zero * observed_arm /
    (nuisance$e * nuisance$r0) * (data$Y - nuisance$a0)
  signal <- signal + gradient[, "a1"] * z_one * observed_arm /
    ((1 - nuisance$e) * nuisance$r1) * (data$Y - nuisance$a1)
  signal
}

msm_iv_binary_evaluate <- function(
    data,
    nuisance_fit,
    Lambda,
    confidence_level = 0.95,
    tie_tolerance = 1e-8,
    incompatibility_tolerance = 1e-8
) {
  predictions <- nuisance_fit$predictions
  msm_iv_binary_validate_inputs(data, predictions)
  if (
    length(confidence_level) != 1L || !is.finite(confidence_level) ||
      confidence_level <= 0 || confidence_level >= 1
  ) {
    stop("confidence_level must lie in (0,1).")
  }
  if (
    length(incompatibility_tolerance) != 1L ||
      !is.finite(incompatibility_tolerance) ||
      incompatibility_tolerance < 0
  ) {
    stop("incompatibility_tolerance must be nonnegative and finite.")
  }

  arms <- list(
    `0` = msm_iv_binary_arm_endpoint(
      predictions, 0, Lambda, tie_tolerance
    ),
    `1` = msm_iv_binary_arm_endpoint(
      predictions, 1, Lambda, tie_tolerance
    )
  )
  diagnostics <- do.call(rbind, lapply(arms, `[[`, "diagnostics"))
  rownames(diagnostics) <- NULL
  no_pointwise_crossing <- all(
    diagnostics$maximum_crossing <= incompatibility_tolerance
  )
  treatment_proxy_no_pointwise_crossing <- all(
    diagnostics$treatment_proxy_maximum_crossing <=
      incompatibility_tolerance
  )
  maximum_incompatibility_rate <- max(diagnostics$incompatibility_rate)
  maximum_crossing <- max(diagnostics$maximum_crossing)
  treatment_proxy_maximum_crossing_rate <- max(
    diagnostics$treatment_proxy_crossing_rate
  )
  treatment_proxy_maximum_crossing <- max(
    diagnostics$treatment_proxy_maximum_crossing
  )

  methods <- c(
    "msm_iv_intersection", "msm_treatment_proxy_intersection",
    "msm_only", "iv_only"
  )
  critical_value <- stats::qnorm(1 - (1 - confidence_level) / 2)
  endpoint_rows <- list()
  summary_rows <- list()
  signal_store <- list()

  for (method in methods) {
    method_no_pointwise_crossing <- switch(
      method,
      msm_iv_intersection = no_pointwise_crossing,
      msm_treatment_proxy_intersection =
        treatment_proxy_no_pointwise_crossing,
      TRUE
    )
    # A finite-sample crossing diagnostic is not a formal compatibility test.
    # Always report the estimator of the combined max/min bound functional and
    # retain any local crossing separately in the diagnostic columns.
    method_is_available <- TRUE
    signals <- list()
    plugin <- list()
    for (treatment in c(0, 1)) {
      arm <- arms[[as.character(treatment)]]
      for (direction in c("lower", "upper")) {
        key <- paste0("mu", treatment, "_", direction)
        plugin[[key]] <- mean(arm$methods[[method]][[direction]])
        signals[[key]] <- msm_iv_binary_one_step_signal(
          data, arm, method, direction
        )
      }
    }
    signals$ate_lower <- signals$mu1_lower - signals$mu0_upper
    signals$ate_upper <- signals$mu1_upper - signals$mu0_lower
    plugin$ate_lower <- plugin$mu1_lower - plugin$mu0_upper
    plugin$ate_upper <- plugin$mu1_upper - plugin$mu0_lower
    signal_store[[method]] <- signals

    endpoint_names <- c(
      "mu1_lower", "mu1_upper", "mu0_lower", "mu0_upper",
      "ate_lower", "ate_upper"
    )
    raw_one_step <- vapply(signals[endpoint_names], mean, numeric(1))
    raw_standard_error <- vapply(
      signals[endpoint_names], stats::sd, numeric(1)
    ) / sqrt(nrow(data))
    raw_plugin <- unlist(plugin[endpoint_names], use.names = TRUE)
    plugin_values <- raw_plugin
    one_step <- raw_one_step
    standard_error <- raw_standard_error
    endpoint_rows[[method]] <- data.frame(
      method = method,
      endpoint = endpoint_names,
      plugin = unname(plugin_values),
      one_step = unname(one_step),
      standard_error = unname(standard_error),
      raw_plugin = unname(raw_plugin),
      raw_one_step = unname(raw_one_step),
      raw_standard_error = unname(raw_standard_error),
      available = method_is_available,
      stringsAsFactors = FALSE
    )

    confidence_lower <- one_step["ate_lower"] -
      critical_value * standard_error["ate_lower"]
    confidence_upper <- one_step["ate_upper"] +
      critical_value * standard_error["ate_upper"]
    raw_confidence_lower <- raw_one_step["ate_lower"] -
      critical_value * raw_standard_error["ate_lower"]
    raw_confidence_upper <- raw_one_step["ate_upper"] +
      critical_value * raw_standard_error["ate_upper"]
    summary_rows[[method]] <- data.frame(
      method = method,
      Lambda = Lambda,
      pointwise_no_crossing = method_no_pointwise_crossing,
      available = method_is_available,
      availability_reason = "reported; crossings retained as diagnostics",
      plugin_lower = unname(plugin_values["ate_lower"]),
      plugin_upper = unname(plugin_values["ate_upper"]),
      plugin_width = unname(
        plugin_values["ate_upper"] - plugin_values["ate_lower"]
      ),
      point_lower = unname(one_step["ate_lower"]),
      point_upper = unname(one_step["ate_upper"]),
      confidence_level = confidence_level,
      confidence_lower = unname(confidence_lower),
      confidence_upper = unname(confidence_upper),
      point_width = unname(
        one_step["ate_upper"] - one_step["ate_lower"]
      ),
      confidence_width = unname(confidence_upper - confidence_lower),
      raw_plugin_lower = unname(raw_plugin["ate_lower"]),
      raw_plugin_upper = unname(raw_plugin["ate_upper"]),
      raw_plugin_width = unname(
        raw_plugin["ate_upper"] - raw_plugin["ate_lower"]
      ),
      raw_point_lower = unname(raw_one_step["ate_lower"]),
      raw_point_upper = unname(raw_one_step["ate_upper"]),
      raw_confidence_lower = unname(raw_confidence_lower),
      raw_confidence_upper = unname(raw_confidence_upper),
      raw_confidence_width = unname(
        raw_confidence_upper - raw_confidence_lower
      ),
      maximum_incompatibility_rate = switch(
        method,
        msm_iv_intersection = maximum_incompatibility_rate,
        msm_treatment_proxy_intersection =
          treatment_proxy_maximum_crossing_rate,
        0
      ),
      maximum_crossing = switch(
        method,
        msm_iv_intersection = maximum_crossing,
        msm_treatment_proxy_intersection =
          treatment_proxy_maximum_crossing,
        0
      ),
      stringsAsFactors = FALSE
    )
  }

  list(
    summary = do.call(rbind, summary_rows),
    endpoints = do.call(rbind, endpoint_rows),
    diagnostics = diagnostics,
    arms = arms,
    signals = signal_store,
    no_pointwise_crossing = no_pointwise_crossing,
    no_pointwise_crossing_by_method = c(
      msm_iv_intersection = no_pointwise_crossing,
      msm_treatment_proxy_intersection =
        treatment_proxy_no_pointwise_crossing
    ),
    configuration = list(
      Lambda = Lambda,
      confidence_level = confidence_level,
      tie_tolerance = tie_tolerance,
      incompatibility_tolerance = incompatibility_tolerance
    )
  )
}
