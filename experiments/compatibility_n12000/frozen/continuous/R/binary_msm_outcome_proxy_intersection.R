# Sharp binary-outcome bounds under an ordinary marginal sensitivity model
# intersected with a binary outcome-confounding-proxy restriction.
#
# For treatment arm t, write
#   p_t(x) = P(T=t | X=x),
#   a_t(x) = P(Y=1 | T=t, X=x), and
#   s_t(x) = P(W=1 | T=t, X=x).
# The outcome-proxy density-ratio envelope is the range of
#   P(W=w | T=1-t,X) / P(W=w | T=t,X), w in {0,1}.
# Intersecting that range with [1/Lambda,Lambda] produces density-ratio
# limits L_t and M_t.  The exact binary missing-arm risk interval is
#   max{a_t L_t, 1-(1-a_t)M_t} <= theta_t <=
#   min{a_t M_t, 1-(1-a_t)L_t}.
#
# This file evaluates the resulting plug-in and cross-fitted one-step bounds.
# Hard min/max branches are differentiated analytically; ties are averaged and
# reported.  No crossed envelope is repaired.

msm_op_binary_validate_inputs <- function(data, predictions) {
  required_data <- c("Y", "T_obs", "W")
  required_predictions <- c("p", "a0", "a1", "w0", "w1")
  missing_data <- setdiff(required_data, names(data))
  missing_predictions <- setdiff(required_predictions, names(predictions))
  if (length(missing_data)) {
    stop("Missing data columns: ", paste(missing_data, collapse = ", "))
  }
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

msm_op_binary_zero_gradient <- function(n) {
  matrix(
    0,
    nrow = n,
    ncol = 5L,
    dimnames = list(NULL, c("p", "a0", "a1", "w0", "w1"))
  )
}

msm_op_binary_choose <- function(
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
  first_weight <- as.numeric(first_active) + 0.5 * as.numeric(tie)
  second_weight <- as.numeric(second_active) + 0.5 * as.numeric(tie)
  list(
    value = value,
    gradient = gradient_first * first_weight +
      gradient_second * second_weight,
    tie = tie,
    first_active = first_active,
    second_active = second_active
  )
}

msm_op_binary_arm_endpoint <- function(
    predictions,
    t_value,
    Lambda,
    tie_tolerance = 1e-8
) {
  if (!t_value %in% c(0, 1)) stop("t_value must be zero or one.")
  if (length(Lambda) != 1L || !is.finite(Lambda) || Lambda < 1) {
    stop("Lambda must be one finite number of at least one.")
  }
  n <- nrow(predictions)
  zero <- msm_op_binary_zero_gradient(n)
  p <- predictions$p
  r <- if (t_value == 1) p else 1 - p
  a <- if (t_value == 1) predictions$a1 else predictions$a0
  s_t <- if (t_value == 1) predictions$w1 else predictions$w0
  s_other <- if (t_value == 1) predictions$w0 else predictions$w1
  a_name <- paste0("a", t_value)
  w_t_name <- paste0("w", t_value)
  w_other_name <- paste0("w", 1 - t_value)

  gradient_a <- zero
  gradient_a[, a_name] <- 1
  gradient_r <- zero
  gradient_r[, "p"] <- if (t_value == 1) 1 else -1

  ratio_one <- s_other / s_t
  gradient_ratio_one <- zero
  gradient_ratio_one[, w_t_name] <- -s_other / s_t^2
  gradient_ratio_one[, w_other_name] <- 1 / s_t

  ratio_zero <- (1 - s_other) / (1 - s_t)
  gradient_ratio_zero <- zero
  gradient_ratio_zero[, w_t_name] <-
    (1 - s_other) / (1 - s_t)^2
  gradient_ratio_zero[, w_other_name] <- -1 / (1 - s_t)

  proxy_lower <- msm_op_binary_choose(
    ratio_zero,
    ratio_one,
    gradient_ratio_zero,
    gradient_ratio_one,
    operation = "min",
    tie_tolerance = tie_tolerance
  )
  proxy_upper <- msm_op_binary_choose(
    ratio_zero,
    ratio_one,
    gradient_ratio_zero,
    gradient_ratio_one,
    operation = "max",
    tie_tolerance = tie_tolerance
  )

  constant_lower_gradient <- zero
  constant_upper_gradient <- zero
  combined_lower <- msm_op_binary_choose(
    rep(1 / Lambda, n),
    proxy_lower$value,
    constant_lower_gradient,
    proxy_lower$gradient,
    operation = "max",
    tie_tolerance = tie_tolerance
  )
  combined_upper <- msm_op_binary_choose(
    rep(Lambda, n),
    proxy_upper$value,
    constant_upper_gradient,
    proxy_upper$gradient,
    operation = "min",
    tie_tolerance = tie_tolerance
  )

  theta_lower_first <- a * combined_lower$value
  gradient_theta_lower_first <-
    gradient_a * combined_lower$value +
    combined_lower$gradient * a
  theta_lower_second <- 1 - (1 - a) * combined_upper$value
  gradient_theta_lower_second <-
    gradient_a * combined_upper$value -
    combined_upper$gradient * (1 - a)
  theta_lower <- msm_op_binary_choose(
    theta_lower_first,
    theta_lower_second,
    gradient_theta_lower_first,
    gradient_theta_lower_second,
    operation = "max",
    tie_tolerance = tie_tolerance
  )

  theta_upper_first <- a * combined_upper$value
  gradient_theta_upper_first <-
    gradient_a * combined_upper$value +
    combined_upper$gradient * a
  theta_upper_second <- 1 - (1 - a) * combined_lower$value
  gradient_theta_upper_second <-
    gradient_a * combined_lower$value -
    combined_lower$gradient * (1 - a)
  theta_upper <- msm_op_binary_choose(
    theta_upper_first,
    theta_upper_second,
    gradient_theta_upper_first,
    gradient_theta_upper_second,
    operation = "min",
    tie_tolerance = tie_tolerance
  )

  lower <- r * a + (1 - r) * theta_lower$value
  upper <- r * a + (1 - r) * theta_upper$value
  gradient_lower <- gradient_r * (a - theta_lower$value) +
    gradient_a * r + theta_lower$gradient * (1 - r)
  gradient_upper <- gradient_r * (a - theta_upper$value) +
    gradient_a * r + theta_upper$gradient * (1 - r)

  # Ordinary MSM comparator, evaluated from the same nuisances.
  msm_l <- 1 / Lambda
  msm_m <- Lambda
  msm_theta_lower_first <- a * msm_l
  msm_theta_lower_second <- 1 - (1 - a) * msm_m
  msm_theta_lower <- msm_op_binary_choose(
    msm_theta_lower_first,
    msm_theta_lower_second,
    gradient_a * msm_l,
    gradient_a * msm_m,
    operation = "max",
    tie_tolerance = tie_tolerance
  )
  msm_theta_upper_first <- a * msm_m
  msm_theta_upper_second <- 1 - (1 - a) * msm_l
  msm_theta_upper <- msm_op_binary_choose(
    msm_theta_upper_first,
    msm_theta_upper_second,
    gradient_a * msm_m,
    gradient_a * msm_l,
    operation = "min",
    tie_tolerance = tie_tolerance
  )
  msm_lower <- r * a + (1 - r) * msm_theta_lower$value
  msm_upper <- r * a + (1 - r) * msm_theta_upper$value
  gradient_msm_lower <- gradient_r * (a - msm_theta_lower$value) +
    gradient_a * r + msm_theta_lower$gradient * (1 - r)
  gradient_msm_upper <- gradient_r * (a - msm_theta_upper$value) +
    gradient_a * r + msm_theta_upper$gradient * (1 - r)

  list(
    treatment = t_value,
    Lambda = Lambda,
    nuisance = list(
      p = p,
      a0 = predictions$a0,
      a1 = predictions$a1,
      w0 = predictions$w0,
      w1 = predictions$w1
    ),
    methods = list(
      msm_outcome_proxy_intersection = list(
        lower = lower,
        upper = upper,
        gradient_lower = gradient_lower,
        gradient_upper = gradient_upper
      ),
      msm_only = list(
        lower = msm_lower,
        upper = msm_upper,
        gradient_lower = gradient_msm_lower,
        gradient_upper = gradient_msm_upper
      )
    ),
    candidates = list(
      ratio_zero = ratio_zero,
      ratio_one = ratio_one,
      proxy_lower = proxy_lower$value,
      proxy_upper = proxy_upper$value,
      combined_lower = combined_lower$value,
      combined_upper = combined_upper$value,
      theta_lower = theta_lower$value,
      theta_upper = theta_upper$value
    ),
    diagnostics = data.frame(
      treatment = t_value,
      Lambda = Lambda,
      proxy_ratio_tie_rate = mean(proxy_lower$tie),
      proxy_active_lower_rate = mean(
        proxy_lower$value > 1 / Lambda + tie_tolerance
      ),
      proxy_active_upper_rate = mean(
        proxy_upper$value < Lambda - tie_tolerance
      ),
      lower_normalization_branch_tie_rate = mean(theta_lower$tie),
      upper_normalization_branch_tie_rate = mean(theta_upper$tie),
      crossing_rate = mean(lower > upper + tie_tolerance),
      maximum_crossing = max(lower - upper),
      minimum_width = min(upper - lower),
      mean_width = mean(upper - lower),
      proxy_lower_minimum = min(proxy_lower$value),
      proxy_upper_maximum = max(proxy_upper$value),
      stringsAsFactors = FALSE
    )
  )
}

msm_op_binary_one_step_signal <- function(
    data,
    arm_endpoint,
    method = c("msm_outcome_proxy_intersection", "msm_only"),
    direction = c("lower", "upper")
) {
  method <- match.arg(method)
  direction <- match.arg(direction)
  component <- arm_endpoint$methods[[method]]
  gradient <- component[[paste0("gradient_", direction)]]
  nuisance <- arm_endpoint$nuisance
  t_one <- as.numeric(data$T_obs == 1)
  t_zero <- 1 - t_one

  signal <- component[[direction]] +
    gradient[, "p"] * (t_one - nuisance$p)
  signal <- signal + gradient[, "a0"] * t_zero / (1 - nuisance$p) *
    (data$Y - nuisance$a0)
  signal <- signal + gradient[, "a1"] * t_one / nuisance$p *
    (data$Y - nuisance$a1)
  signal <- signal + gradient[, "w0"] * t_zero / (1 - nuisance$p) *
    (data$W - nuisance$w0)
  signal <- signal + gradient[, "w1"] * t_one / nuisance$p *
    (data$W - nuisance$w1)
  signal
}

msm_op_binary_evaluate <- function(
    data,
    nuisance_fit,
    Lambda,
    confidence_level = 0.95,
    tie_tolerance = 1e-8
) {
  predictions <- nuisance_fit$predictions
  msm_op_binary_validate_inputs(data, predictions)
  arms <- list(
    `0` = msm_op_binary_arm_endpoint(
      predictions, 0, Lambda, tie_tolerance
    ),
    `1` = msm_op_binary_arm_endpoint(
      predictions, 1, Lambda, tie_tolerance
    )
  )
  diagnostics <- do.call(rbind, lapply(arms, `[[`, "diagnostics"))
  rownames(diagnostics) <- NULL
  methods <- c("msm_outcome_proxy_intersection", "msm_only")
  critical_value <- stats::qnorm(1 - (1 - confidence_level) / 2)
  endpoint_rows <- list()
  summary_rows <- list()
  signals_store <- list()

  for (method in methods) {
    plugin <- list()
    signals <- list()
    for (treatment in c(0, 1)) {
      arm <- arms[[as.character(treatment)]]
      for (direction in c("lower", "upper")) {
        key <- paste0("mu", treatment, "_", direction)
        plugin[[key]] <- mean(arm$methods[[method]][[direction]])
        signals[[key]] <- msm_op_binary_one_step_signal(
          data, arm, method, direction
        )
      }
    }
    plugin$ate_lower <- plugin$mu1_lower - plugin$mu0_upper
    plugin$ate_upper <- plugin$mu1_upper - plugin$mu0_lower
    signals$ate_lower <- signals$mu1_lower - signals$mu0_upper
    signals$ate_upper <- signals$mu1_upper - signals$mu0_lower
    signals_store[[method]] <- signals
    endpoint_names <- c(
      "mu1_lower", "mu1_upper", "mu0_lower", "mu0_upper",
      "ate_lower", "ate_upper"
    )
    plugin_values <- unlist(plugin[endpoint_names], use.names = TRUE)
    one_step <- vapply(signals[endpoint_names], mean, numeric(1))
    standard_error <- vapply(
      signals[endpoint_names], stats::sd, numeric(1)
    ) / sqrt(nrow(data))
    endpoint_rows[[method]] <- data.frame(
      method = method,
      endpoint = endpoint_names,
      plugin = unname(plugin_values),
      one_step = unname(one_step),
      standard_error = unname(standard_error),
      stringsAsFactors = FALSE
    )
    confidence_lower <- one_step["ate_lower"] -
      critical_value * standard_error["ate_lower"]
    confidence_upper <- one_step["ate_upper"] +
      critical_value * standard_error["ate_upper"]
    summary_rows[[method]] <- data.frame(
      method = method,
      Lambda = Lambda,
      plugin_lower = unname(plugin_values["ate_lower"]),
      plugin_upper = unname(plugin_values["ate_upper"]),
      plugin_width = unname(
        plugin_values["ate_upper"] - plugin_values["ate_lower"]
      ),
      point_lower = unname(one_step["ate_lower"]),
      point_upper = unname(one_step["ate_upper"]),
      point_width = unname(
        one_step["ate_upper"] - one_step["ate_lower"]
      ),
      confidence_level = confidence_level,
      confidence_lower = unname(confidence_lower),
      confidence_upper = unname(confidence_upper),
      confidence_width = unname(confidence_upper - confidence_lower),
      stringsAsFactors = FALSE
    )
  }

  list(
    summary = do.call(rbind, summary_rows),
    endpoints = do.call(rbind, endpoint_rows),
    diagnostics = diagnostics,
    arms = arms,
    signals = signals_store,
    configuration = list(
      Lambda = Lambda,
      confidence_level = confidence_level,
      tie_tolerance = tie_tolerance
    )
  )
}
