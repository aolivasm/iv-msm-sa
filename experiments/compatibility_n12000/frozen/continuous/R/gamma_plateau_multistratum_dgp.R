# Multi-stratum version of the exact gamma-envelope plateau design.
#
# Each covariate stratum has the same three-component conditional residual law
# and two positive-mass sets on which gamma_1 is exactly at an envelope switch.
# The strata are represented by a saturated set of numeric indicators.  Their
# outcome locations differ, so cross-fitted gamma fits can make separate branch
# decisions while the threshold optimizer remains unique within every stratum.

make_gamma_plateau_multistratum_spec <- function(
    sensitivity_lambda = 4,
    probability_z1 = 0.7,
    treatment_probability_z0 = 0.1,
    treatment_probability_z1 = 0.5,
    number_of_strata = 4L,
    stratum_probabilities = NULL,
    component_probabilities = c(
      regular = 0.15,
      lower_tie = 0.75,
      upper_tie = 0.10
    ),
    component_order = c("lower_tie", "regular", "upper_tie"),
    component_width = 0.04,
    stratum_location_span = 0.6,
    treatment_effect = 0.75
) {
  if (!exists("make_gamma_plateau_variance_spec", mode = "function")) {
    stop("Source R/gamma_plateau_variance_dgp.R first.")
  }
  number_of_strata <- as.integer(number_of_strata)
  if (
    length(number_of_strata) != 1L || is.na(number_of_strata) ||
      number_of_strata < 2L || number_of_strata > 12L
  ) {
    stop("number_of_strata must be one integer from 2 through 12.")
  }
  if (is.null(stratum_probabilities)) {
    stratum_probabilities <- rep(1 / number_of_strata, number_of_strata)
  }
  stratum_probabilities <- as.numeric(stratum_probabilities)
  if (
    length(stratum_probabilities) != number_of_strata ||
      any(!is.finite(stratum_probabilities)) ||
      any(stratum_probabilities <= 0) ||
      abs(sum(stratum_probabilities) - 1) > 1e-12
  ) {
    stop("stratum_probabilities must be positive and sum to one.")
  }
  stratum_location_span <- as.numeric(stratum_location_span)
  if (
    length(stratum_location_span) != 1L ||
      !is.finite(stratum_location_span) || stratum_location_span < 0
  ) {
    stop("stratum_location_span must be one nonnegative finite number.")
  }

  base_specification <- make_gamma_plateau_variance_spec(
    sensitivity_lambda = sensitivity_lambda,
    probability_z1 = probability_z1,
    treatment_probability_z0 = treatment_probability_z0,
    treatment_probability_z1 = treatment_probability_z1,
    component_probabilities = component_probabilities,
    component_order = component_order,
    component_width = component_width,
    treatment_effect = treatment_effect
  )
  stratum_labels <- paste0("stratum_", seq_len(number_of_strata))
  locations <- if (number_of_strata == 2L) {
    c(-0.5, 0.5) * stratum_location_span
  } else {
    seq(
      -stratum_location_span / 2,
      stratum_location_span / 2,
      length.out = number_of_strata
    )
  }
  indicator_columns <- paste0("X_stratum_", seq.int(2L, number_of_strata))

  specification <- base_specification
  specification$name <- "multi_stratum_unique_threshold_gamma_plateau"
  specification$base_specification <- base_specification
  specification$number_of_strata <- number_of_strata
  specification$stratum_labels <- stratum_labels
  specification$stratum_probabilities <- setNames(
    stratum_probabilities, stratum_labels
  )
  specification$stratum_locations <- setNames(locations, stratum_labels)
  specification$stratum_location_span <- stratum_location_span
  specification$x_columns <- indicator_columns
  specification$covariate_design <- "saturated_categorical_strata"
  specification$tie_probability_by_stratum <- rep(
    base_specification$tie_probability, number_of_strata
  )
  specification$minimum_stratum_probability <- min(stratum_probabilities)
  class(specification) <- c(
    "gamma_plateau_multistratum_specification", "list"
  )
  specification
}

gamma_plateau_multistratum_x_columns <- function(specification) {
  if (!inherits(specification, "gamma_plateau_multistratum_specification")) {
    stop("specification must be a multi-stratum gamma-plateau specification.")
  }
  specification$x_columns
}

generate_gamma_plateau_multistratum_data <- function(
    n,
    specification = make_gamma_plateau_multistratum_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  stratum <- sample(
    specification$stratum_labels,
    n,
    replace = TRUE,
    prob = specification$stratum_probabilities
  )
  data <- generate_gamma_plateau_variance_data(
    n = n,
    specification = specification$base_specification,
    seed = NULL
  )
  data$Y <- data$Y + unname(specification$stratum_locations[stratum])
  data$X <- match(stratum, specification$stratum_labels)
  data$covariate_stratum <- factor(
    stratum, levels = specification$stratum_labels
  )
  for (index in seq_along(specification$x_columns)) {
    stratum_index <- index + 1L
    data[[specification$x_columns[index]]] <- as.numeric(
      data$X == stratum_index
    )
  }
  data
}

gamma_plateau_multistratum_threshold_margin <- function(
    specification,
    envelope_smoothing_epsilon = 0,
    nodes_per_component = 10000L
) {
  base <- specification$base_specification
  rule <- gamma_plateau_variance_rule(base, nodes_per_component)
  n <- length(rule$y)
  nuisance <- list(
    e = rep(1 - base$probability_z1, n),
    pi0 = rep(base$pi0, n),
    pi1 = rep(base$pi1, n)
  )
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = unname(base$gamma1_values[rule$region]),
    t_value = 1L,
    Lambda0 = base$sensitivity_lambda,
    Lambda1 = base$sensitivity_lambda,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    gamma_compatibility_action = "stop"
  )
  Lbar <- sum(rule$weight * envelope$L)
  Ubar <- sum(rule$weight * envelope$U)
  tau <- (Ubar - 1) / (Ubar - Lbar)
  ordered_region <- unique(rule$region[order(rule$y)])
  weighted_mass <- vapply(
    ordered_region,
    function(label) {
      index <- rule$region == label
      sum(rule$weight[index] * envelope$D[index])
    },
    numeric(1)
  )
  boundaries <- cumsum(weighted_mass) / sum(weighted_mass)
  boundaries <- boundaries[-length(boundaries)]
  c(
    lower = min(abs((1 - tau) - boundaries)),
    upper = min(abs(tau - boundaries))
  )
}

gamma_plateau_multistratum_population <- function(
    specification = make_gamma_plateau_multistratum_spec(),
    envelope_smoothing_epsilon = 0,
    nodes_per_component = 10000L
) {
  base_population <- gamma_plateau_variance_population(
    specification$base_specification,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    nodes_per_component = nodes_per_component
  )
  average_location <- sum(
    specification$stratum_probabilities * specification$stratum_locations
  )
  margin <- gamma_plateau_multistratum_threshold_margin(
    specification,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    nodes_per_component = nodes_per_component
  )
  base_population$mu1_lower <- base_population$mu1_lower + average_location
  base_population$mu1_upper <- base_population$mu1_upper + average_location
  base_population$lower_threshold <- NA_real_
  base_population$upper_threshold <- NA_real_
  base_population$number_of_strata <- specification$number_of_strata
  base_population$minimum_lower_threshold_margin <- unname(margin["lower"])
  base_population$minimum_upper_threshold_margin <- unname(margin["upper"])
  base_population
}

gamma_plateau_multistratum_population_grid <- function(
    specification = make_gamma_plateau_multistratum_spec(),
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
    gamma_plateau_multistratum_population(
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
