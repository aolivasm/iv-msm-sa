# Observed-data law for branch-variance instability at exact gamma-envelope
# ties.  The construction keeps the threshold optimizer unique.
#
# For T = 1, the baseline Z = 0 outcome residual has three separated,
# continuous components.  Their probabilities are
#
#   regular = 0.20, lower-envelope tie = 0.70, upper-envelope tie = 0.10.
#
# The lower-tie component is placed below the regular component and the
# upper-tie component above it.  The density ratio gamma_1 is exactly equal
# to the corresponding envelope switch on the two tie components.  Its value
# on the regular component enforces E_0{gamma_1(Y)} = 1.  At Lambda = 4,
# pi0 = 0.10, and pi1 = 0.50, the upper-bound threshold lies strictly inside
# the regular component, so the envelope is nonregular while the weighted
# check-loss optimizer remains unique.

make_gamma_plateau_variance_spec <- function(
    sensitivity_lambda = 4,
    probability_z1 = 0.8,
    treatment_probability_z0 = 0.1,
    treatment_probability_z1 = 0.5,
    component_probabilities = c(
      regular = 0.20,
      lower_tie = 0.70,
      upper_tie = 0.10
    ),
    component_order = c("lower_tie", "regular", "upper_tie"),
    component_width = 0.04,
    treatment_effect = 0.75
) {
  Lambda <- as.numeric(sensitivity_lambda)
  probability_z1 <- as.numeric(probability_z1)
  pi0 <- as.numeric(treatment_probability_z0)
  pi1 <- as.numeric(treatment_probability_z1)
  component_probability <- as.numeric(component_probabilities)
  names(component_probability) <- names(component_probabilities)
  labels <- c("regular", "lower_tie", "upper_tie")
  component_order <- as.character(component_order)
  component_width <- as.numeric(component_width)
  treatment_effect <- as.numeric(treatment_effect)
  if (
    any(!is.finite(c(
      Lambda, probability_z1, pi0, pi1,
      component_probability, component_width, treatment_effect
    ))) || Lambda <= 1 || probability_z1 <= 0 || probability_z1 >= 1 ||
      pi0 <= 0 || pi0 >= pi1 || pi1 >= 1 ||
      !identical(sort(names(component_probability)), sort(labels)) ||
      any(component_probability <= 0) ||
      abs(sum(component_probability) - 1) > 1e-12 ||
      length(component_order) != 3L ||
      !identical(sort(component_order), sort(labels)) ||
      component_width <= 0 || component_width >= 0.25
  ) {
    stop("Invalid gamma-plateau variance specification.")
  }

  lower_switch <- (
    pi0 + (1 - pi0) / Lambda
  ) / (
    pi1 + (1 - pi1) / Lambda
  )
  upper_switch <- (
    pi0 + (1 - pi0) * Lambda
  ) / (
    pi1 + (1 - pi1) * Lambda
  )
  regular_gamma <- (
    1 - component_probability["lower_tie"] * lower_switch -
      component_probability["upper_tie"] * upper_switch
  ) / component_probability["regular"]
  gamma_value <- c(
    regular = unname(regular_gamma),
    lower_tie = unname(lower_switch),
    upper_tie = unname(upper_switch)
  )

  compatible_lower_t1 <- (
    (1 - pi0) / Lambda + pi0
  ) / (
    (1 - pi1) * Lambda + pi1
  )
  compatible_upper_t1 <- (
    Lambda * (1 - pi0) + pi0
  ) / (
    (1 - pi1) / Lambda + pi1
  )
  if (
    any(gamma_value <= 0) ||
      any(gamma_value < compatible_lower_t1) ||
      any(gamma_value > compatible_upper_t1)
  ) {
    stop("The gamma plateaus are not compatible with the sensitivity model.")
  }

  centers <- c(0, 0.5, 1)
  left_by_position <- centers - component_width / 2
  right_by_position <- centers + component_width / 2
  component_left <- component_right <- setNames(numeric(3L), labels)
  for (position in seq_along(component_order)) {
    label <- component_order[position]
    component_left[label] <- left_by_position[position]
    component_right[label] <- right_by_position[position]
  }

  specification <- list(
    name = "unique_threshold_gamma_plateau_branch_variance",
    sensitivity_lambda = Lambda,
    probability_z1 = probability_z1,
    pi0 = pi0,
    pi1 = pi1,
    overlap_constant = min(
      probability_z1, 1 - probability_z1,
      pi0, 1 - pi0, pi1, 1 - pi1
    ),
    marginal_first_stage = pi1 - pi0,
    component_labels = labels,
    component_probabilities = component_probability,
    component_order = component_order,
    component_left = component_left,
    component_right = component_right,
    component_width = component_width,
    interval_probabilities = component_probability,
    covariate_design = "constant",
    outcome_geometry = "separated_clusters",
    separated_cluster_width = component_width,
    gamma1_values = gamma_value,
    gamma1_lower_switch = lower_switch,
    gamma1_upper_switch = upper_switch,
    gamma1_compatible_interval = c(
      lower = compatible_lower_t1,
      upper = compatible_upper_t1
    ),
    gamma0_value = 1,
    treatment_effect = treatment_effect,
    tie_probability = sum(component_probability[c(
      "lower_tie", "upper_tie"
    )])
  )
  class(specification) <- c(
    "gamma_plateau_variance_specification", "list"
  )
  specification
}

gamma_plateau_variance_region <- function(
    residual,
    specification = make_gamma_plateau_variance_spec()
) {
  residual <- as.numeric(residual)
  if (any(!is.finite(residual))) stop("residual must be finite.")
  distance <- vapply(
    specification$component_labels,
    function(label) {
      midpoint <- (
        specification$component_left[label] +
          specification$component_right[label]
      ) / 2
      abs(residual - midpoint)
    },
    numeric(length(residual))
  )
  if (is.vector(distance)) {
    distance <- matrix(
      distance,
      nrow = length(residual),
      ncol = length(specification$component_labels)
    )
  }
  specification$component_labels[max.col(-distance, ties.method = "first")]
}

generate_gamma_plateau_variance_data <- function(
    n,
    specification = make_gamma_plateau_variance_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  z <- stats::rbinom(n, 1L, specification$probability_z1)
  treatment_probability <- ifelse(
    z == 1L, specification$pi1, specification$pi0
  )
  treatment <- stats::rbinom(n, 1L, treatment_probability)

  region <- sample(
    specification$component_labels,
    n,
    replace = TRUE,
    prob = specification$component_probabilities
  )
  treated_z1 <- treatment == 1L & z == 1L
  if (any(treated_z1)) {
    comparison_probability <-
      specification$component_probabilities * specification$gamma1_values
    comparison_probability <- comparison_probability /
      sum(comparison_probability)
    region[treated_z1] <- sample(
      specification$component_labels,
      sum(treated_z1),
      replace = TRUE,
      prob = comparison_probability
    )
  }
  residual <- stats::runif(
    n,
    specification$component_left[region],
    specification$component_right[region]
  )
  data.frame(
    Y = specification$treatment_effect * treatment + residual,
    T_obs = treatment,
    Z = z,
    X = 0,
    outcome_region = factor(
      region,
      levels = specification$component_labels
    ),
    e_true = rep(1 - specification$probability_z1, n),
    pi0_true = rep(specification$pi0, n),
    pi1_true = rep(specification$pi1, n),
    gamma0_oracle = rep(1, n),
    gamma1_oracle = unname(specification$gamma1_values[region]),
    stringsAsFactors = FALSE
  )
}

gamma_plateau_variance_rule <- function(
    specification = make_gamma_plateau_variance_spec(),
    nodes_per_component = 10000L
) {
  nodes_per_component <- as.integer(nodes_per_component)
  if (is.na(nodes_per_component) || nodes_per_component < 100L) {
    stop("nodes_per_component must be an integer of at least 100.")
  }
  u <- (seq_len(nodes_per_component) - 0.5) / nodes_per_component
  y <- weight <- numeric()
  region <- character()
  for (label in specification$component_labels) {
    y <- c(
      y,
      specification$component_left[label] +
        specification$component_width * u
    )
    weight <- c(
      weight,
      rep(
        specification$component_probabilities[label] /
          nodes_per_component,
        nodes_per_component
      )
    )
    region <- c(region, rep(label, nodes_per_component))
  }
  list(y = y, weight = weight, region = region)
}

gamma_plateau_variance_population <- function(
    specification = make_gamma_plateau_variance_spec(),
    envelope_smoothing_epsilon = 0,
    nodes_per_component = 10000L
) {
  if (!exists("compute_projected_envelopes", mode = "function")) {
    stop("Source R/sharp_bounds_plugin.R before computing population bounds.")
  }
  rule <- gamma_plateau_variance_rule(
    specification, nodes_per_component
  )
  n <- length(rule$y)
  nuisance <- list(
    e = rep(1 - specification$probability_z1, n),
    pi0 = rep(specification$pi0, n),
    pi1 = rep(specification$pi1, n)
  )
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = unname(specification$gamma1_values[rule$region]),
    t_value = 1L,
    Lambda0 = specification$sensitivity_lambda,
    Lambda1 = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    gamma_compatibility_action = "stop"
  )
  Lbar <- sum(rule$weight * envelope$L)
  Ubar <- sum(rule$weight * envelope$U)
  tau <- (Ubar - 1) / (Ubar - Lbar)
  weighted_quantile <- function(probability) {
    ordering <- order(rule$y)
    value <- rule$y[ordering]
    weight <- (rule$weight * envelope$D)[ordering]
    weight <- weight / sum(weight)
    value[which(cumsum(weight) >= probability)[1L]]
  }
  q_lower <- weighted_quantile(1 - tau)
  q_upper <- weighted_quantile(tau)
  lambda_lower <- envelope$L + envelope$D * (rule$y <= q_lower)
  lambda_upper <- envelope$U - envelope$D * (rule$y <= q_upper)
  baseline_mean <- sum(rule$weight * rule$y)
  lower <- specification$treatment_effect +
    specification$pi0 * baseline_mean +
    (1 - specification$pi0) *
      sum(rule$weight * lambda_lower * rule$y)
  upper <- specification$treatment_effect +
    specification$pi0 * baseline_mean +
    (1 - specification$pi0) *
      sum(rule$weight * lambda_upper * rule$y)
  threshold_region <- function(q) {
    specification$component_labels[
      q >= specification$component_left &
        q <= specification$component_right
    ]
  }
  data.frame(
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    mu1_lower = lower,
    mu1_upper = upper,
    lower_threshold = specification$treatment_effect + q_lower,
    upper_threshold = specification$treatment_effect + q_upper,
    lower_threshold_region = threshold_region(q_lower),
    upper_threshold_region = threshold_region(q_upper),
    mean_lower_smoothing_gap = sum(
      rule$weight * envelope$lower_smoothing_gap
    ),
    mean_upper_smoothing_gap = sum(
      rule$weight * envelope$upper_smoothing_gap
    ),
    stringsAsFactors = FALSE
  )
}

gamma_plateau_variance_population_grid <- function(
    specification = make_gamma_plateau_variance_spec(),
    epsilon_grid = c(0, 0.01),
    nodes_per_component = 10000L
) {
  epsilon_grid <- sort(unique(as.numeric(epsilon_grid)))
  if (
    !length(epsilon_grid) || any(!is.finite(epsilon_grid)) ||
      any(epsilon_grid < 0) || !any(epsilon_grid == 0)
  ) {
    stop("epsilon_grid must include zero and valid nonnegative values.")
  }
  output <- do.call(rbind, lapply(epsilon_grid, function(epsilon) {
    gamma_plateau_variance_population(
      specification,
      envelope_smoothing_epsilon = epsilon,
      nodes_per_component = nodes_per_component
    )
  }))
  rownames(output) <- NULL
  hard <- output[output$envelope_smoothing_epsilon == 0, ]
  output$lower_displacement_from_sharp <-
    output$mu1_lower - hard$mu1_lower
  output$upper_displacement_from_sharp <-
    output$mu1_upper - hard$mu1_upper
  output
}

# Full treatment-arm and ATE population targets. These functions use the
# observed-data law only for Monte Carlo scoring; no population quantity is
# supplied to the estimator.
gamma_plateau_variance_arm_population <- function(
    treatment,
    endpoint = c("lower", "upper"),
    specification = make_gamma_plateau_variance_spec(),
    envelope_smoothing_epsilon = 0,
    nodes_per_component = 10000L
) {
  endpoint <- match.arg(endpoint)
  treatment <- as.integer(treatment)
  if (length(treatment) != 1L || is.na(treatment) ||
      !treatment %in% 0:1) {
    stop("treatment must equal zero or one.")
  }
  if (!exists("compute_projected_envelopes", mode = "function")) {
    stop("Source R/sharp_bounds_plugin.R before computing population bounds.")
  }
  rule <- gamma_plateau_variance_rule(specification, nodes_per_component)
  n <- length(rule$y)
  nuisance <- list(
    e = rep(1 - specification$probability_z1, n),
    pi0 = rep(specification$pi0, n),
    pi1 = rep(specification$pi1, n)
  )
  gamma <- if (treatment == 1L) {
    unname(specification$gamma1_values[rule$region])
  } else {
    rep(specification$gamma0_value, n)
  }
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = gamma,
    t_value = treatment,
    Lambda0 = specification$sensitivity_lambda,
    Lambda1 = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    gamma_compatibility_action = "stop"
  )
  if (any(envelope$gamma_incompatible_original) ||
      any(envelope$L >= envelope$U)) {
    stop("The gamma-plateau law is incompatible at the requested Lambda.")
  }
  Lbar <- sum(rule$weight * envelope$L)
  Ubar <- sum(rule$weight * envelope$U)
  tau <- (Ubar - 1) / (Ubar - Lbar)
  weighted_quantile <- function(probability) {
    ordering <- order(rule$y)
    value <- rule$y[ordering]
    weight <- (rule$weight * envelope$D)[ordering]
    weight <- weight / sum(weight)
    value[which(cumsum(weight) >= probability)[1L]]
  }
  q <- if (endpoint == "lower") {
    weighted_quantile(1 - tau)
  } else {
    weighted_quantile(tau)
  }
  lambda <- if (endpoint == "lower") {
    envelope$L + envelope$D * (rule$y <= q)
  } else {
    envelope$U - envelope$D * (rule$y <= q)
  }
  location <- treatment * specification$treatment_effect
  y <- location + rule$y
  r0 <- if (treatment == 1L) specification$pi0 else 1 - specification$pi0
  value <- r0 * sum(rule$weight * y) +
    (1 - r0) * sum(rule$weight * lambda * y)
  list(
    value = value,
    threshold = location + q,
    Lbar = Lbar,
    Ubar = Ubar,
    tau = tau,
    mean_lower_smoothing_gap = sum(
      rule$weight * envelope$lower_smoothing_gap
    ),
    mean_upper_smoothing_gap = sum(
      rule$weight * envelope$upper_smoothing_gap
    )
  )
}

gamma_plateau_variance_ate_population <- function(
    specification = make_gamma_plateau_variance_spec(),
    envelope_smoothing_epsilon = 0,
    nodes_per_component = 10000L
) {
  evaluated <- list()
  for (treatment in 0:1) {
    for (endpoint in c("lower", "upper")) {
      name <- paste0("mu", treatment, "_", endpoint)
      evaluated[[name]] <- gamma_plateau_variance_arm_population(
        treatment = treatment,
        endpoint = endpoint,
        specification = specification,
        envelope_smoothing_epsilon = envelope_smoothing_epsilon,
        nodes_per_component = nodes_per_component
      )
    }
  }
  output <- data.frame(
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    mu1_lower = evaluated$mu1_lower$value,
    mu1_upper = evaluated$mu1_upper$value,
    mu0_lower = evaluated$mu0_lower$value,
    mu0_upper = evaluated$mu0_upper$value,
    stringsAsFactors = FALSE
  )
  output$ate_lower <- output$mu1_lower - output$mu0_upper
  output$ate_upper <- output$mu1_upper - output$mu0_lower
  output$ate_width <- output$ate_upper - output$ate_lower
  output$reference_ate <- specification$treatment_effect
  output
}

gamma_plateau_variance_ate_population_grid <- function(
    specification = make_gamma_plateau_variance_spec(),
    epsilon_grid = c(0, 0.01, 0.05),
    nodes_per_component = 10000L
) {
  epsilon_grid <- sort(unique(as.numeric(epsilon_grid)))
  if (length(epsilon_grid) < 2L || any(!is.finite(epsilon_grid)) ||
      any(epsilon_grid < 0) || !any(epsilon_grid == 0)) {
    stop("epsilon_grid must include zero and another nonnegative value.")
  }
  output <- do.call(rbind, lapply(epsilon_grid, function(epsilon) {
    gamma_plateau_variance_ate_population(
      specification = specification,
      envelope_smoothing_epsilon = epsilon,
      nodes_per_component = nodes_per_component
    )
  }))
  hard <- output[output$envelope_smoothing_epsilon == 0, , drop = FALSE]
  output$ate_lower_shift_from_hard <- output$ate_lower - hard$ate_lower
  output$ate_upper_shift_from_hard <- output$ate_upper - hard$ate_upper
  output$ate_width_increase_from_hard <- output$ate_width - hard$ate_width
  rownames(output) <- NULL
  output
}
