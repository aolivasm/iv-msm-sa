# Two-stratum observed-data law for estimated density-ratio envelope ties.
#
# Eighty percent of the population lies in a boundary stratum with
# pi0 = 0.10 and pi1 = 0.12.  Within this stratum gamma_1(Y, X) has
# positive-mass plateaus at the exact lower and upper envelope switching
# values.  The remaining twenty percent is a regular strong-first-stage
# stratum with pi0 = 0.10, pi1 = 0.90, and gamma_1 = 1.  Thus the default
# marginal first stage is 0.176 even though the density-ratio switches are
# close enough to one to be learnable by the direct RKHS estimator.

make_gamma_boundary_two_stratum_spec <- function(
    sensitivity_lambda = 4,
    probability_z1 = 0.5,
    boundary_stratum_probability = 0.8,
    treatment_probability_z0 = 0.1,
    boundary_treatment_probability_z1 = 0.12,
    regular_treatment_probability_z1 = 0.9,
    baseline_interval_probabilities = c(0.2, 0.6, 0.2),
    treatment_effect = 0.75,
    outcome_geometry = c(
      "contiguous", "separated_clusters", "near_discrete_clusters"
    ),
    separated_cluster_width = 0.1
) {
  outcome_geometry <- match.arg(outcome_geometry)
  separated_cluster_width <- as.numeric(separated_cluster_width)
  Lambda <- as.numeric(sensitivity_lambda)
  probability_z1 <- as.numeric(probability_z1)
  boundary_probability <- as.numeric(boundary_stratum_probability)
  pi0 <- as.numeric(treatment_probability_z0)
  pi1_boundary <- as.numeric(boundary_treatment_probability_z1)
  pi1_regular <- as.numeric(regular_treatment_probability_z1)
  interval_probability <- as.numeric(baseline_interval_probabilities)
  treatment_effect <- as.numeric(treatment_effect)
  values <- c(
    Lambda, probability_z1, boundary_probability, pi0,
    pi1_boundary, pi1_regular, interval_probability, treatment_effect,
    separated_cluster_width
  )
  if (
    any(!is.finite(values)) || Lambda <= 1 ||
      probability_z1 <= 0 || probability_z1 >= 1 ||
      boundary_probability <= 0 || boundary_probability >= 1 ||
      pi0 <= 0 || pi0 >= 1 ||
      pi1_boundary <= pi0 || pi1_boundary >= pi1_regular ||
      pi1_regular >= 1 || length(interval_probability) != 3L ||
      any(interval_probability <= 0) ||
      abs(sum(interval_probability) - 1) > 1e-12 ||
      separated_cluster_width <= 0 || separated_cluster_width >= 0.4
  ) {
    stop("Invalid two-stratum gamma-boundary specification.")
  }

  lower_switch <- (
    pi0 + (1 - pi0) / Lambda
  ) / (
    pi1_boundary + (1 - pi1_boundary) / Lambda
  )
  upper_switch <- (
    pi0 + (1 - pi0) * Lambda
  ) / (
    pi1_boundary + (1 - pi1_boundary) * Lambda
  )
  regular_piece <- (
    1 - interval_probability[2L] * lower_switch -
      interval_probability[3L] * upper_switch
  ) / interval_probability[1L]
  boundary_gamma <- c(regular_piece, lower_switch, upper_switch)

  compatible_interval <- function(pi1, treatment) {
    if (treatment == 1L) {
      return(c(
        lower = ((1 - pi0) / Lambda + pi0) /
          ((1 - pi1) * Lambda + pi1),
        upper = (Lambda * (1 - pi0) + pi0) /
          ((1 - pi1) / Lambda + pi1)
      ))
    }
    c(
      lower = (pi0 / Lambda + 1 - pi0) /
        (pi1 * Lambda + 1 - pi1),
      upper = (Lambda * pi0 + 1 - pi0) /
        (pi1 / Lambda + 1 - pi1)
    )
  }
  boundary_t1_interval <- compatible_interval(pi1_boundary, 1L)
  regular_t1_interval <- compatible_interval(pi1_regular, 1L)
  boundary_t0_interval <- compatible_interval(pi1_boundary, 0L)
  regular_t0_interval <- compatible_interval(pi1_regular, 0L)
  if (
    any(boundary_gamma <= 0) ||
      any(boundary_gamma < boundary_t1_interval["lower"]) ||
      any(boundary_gamma > boundary_t1_interval["upper"]) ||
      1 < regular_t1_interval["lower"] ||
      1 > regular_t1_interval["upper"] ||
      1 < boundary_t0_interval["lower"] ||
      1 > boundary_t0_interval["upper"] ||
      1 < regular_t0_interval["lower"] ||
      1 > regular_t0_interval["upper"]
  ) {
    stop("The two-stratum density ratios are not model-compatible.")
  }

  x_probability <- c(
    boundary = boundary_probability,
    regular = 1 - boundary_probability
  )
  pi0_by_x <- c(boundary = pi0, regular = pi0)
  pi1_by_x <- c(boundary = pi1_boundary, regular = pi1_regular)
  interval_left <- switch(
    outcome_geometry,
    contiguous = c(0, cumsum(interval_probability)[-3L]),
    separated_clusters = c(
      0,
      0.5 - separated_cluster_width / 2,
      1 - separated_cluster_width
    ),
    near_discrete_clusters = c(0, 0.499, 0.998)
  )
  interval_right <- switch(
    outcome_geometry,
    contiguous = cumsum(interval_probability),
    separated_clusters = c(
      separated_cluster_width,
      0.5 + separated_cluster_width / 2,
      1
    ),
    near_discrete_clusters = c(0.002, 0.501, 1)
  )
  specification <- list(
    name = "two_stratum_positive_mass_gamma_boundary_plateaus",
    sensitivity_lambda = Lambda,
    probability_z1 = probability_z1,
    x_values = c(boundary = 0, regular = 1),
    x_probabilities = x_probability,
    pi0 = pi0_by_x,
    pi1 = pi1_by_x,
    overlap_constant = min(
      probability_z1, 1 - probability_z1,
      pi0_by_x, 1 - pi0_by_x, pi1_by_x, 1 - pi1_by_x
    ),
    marginal_first_stage = sum(
      x_probability * (pi1_by_x - pi0_by_x)
    ),
    interval_labels = c("regular_piece", "lower_tie", "upper_tie"),
    interval_probabilities = interval_probability,
    interval_left = interval_left,
    interval_right = interval_right,
    interval_boundaries = c(0, cumsum(interval_probability)),
    gamma1_values = boundary_gamma,
    gamma1_boundary_values = boundary_gamma,
    gamma1_regular_values = rep(1, 3L),
    gamma1_lower_switch = lower_switch,
    gamma1_upper_switch = upper_switch,
    gamma1_compatible_boundary = boundary_t1_interval,
    gamma1_compatible_regular = regular_t1_interval,
    gamma0_value = 1,
    treatment_effect = treatment_effect,
    covariate_design = "two_stratum",
    outcome_geometry = outcome_geometry,
    separated_cluster_width = if (
      outcome_geometry == "near_discrete_clusters"
    ) 0.002 else if (outcome_geometry == "separated_clusters") {
      separated_cluster_width
    } else {
      NA_real_
    }
  )
  class(specification) <- c(
    "gamma_boundary_two_stratum_specification", "list"
  )
  specification
}

gamma_boundary_two_stratum_region <- function(
    residual,
    specification = make_gamma_boundary_two_stratum_spec()
) {
  residual <- as.numeric(residual)
  if (any(!is.finite(residual))) stop("residual must be finite.")
  region <- rep(NA_integer_, length(residual))
  for (index in seq_along(specification$interval_left)) {
    inside <- residual >= specification$interval_left[index] &
      residual <= specification$interval_right[index]
    region[inside & is.na(region)] <- index
  }
  if (anyNA(region)) {
    midpoint <- (
      specification$interval_left + specification$interval_right
    ) / 2
    missing <- which(is.na(region))
    region[missing] <- vapply(residual[missing], function(value) {
      which.min(abs(value - midpoint))
    }, integer(1))
  }
  region
}

generate_gamma_boundary_two_stratum_data <- function(
    n,
    specification = make_gamma_boundary_two_stratum_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  boundary_x <- stats::rbinom(
    n, 1L, specification$x_probabilities["regular"]
  )
  x_index <- boundary_x + 1L
  z <- stats::rbinom(n, 1L, specification$probability_z1)
  pi0 <- specification$pi0[x_index]
  pi1 <- specification$pi1[x_index]
  treatment_probability <- ifelse(z == 1L, pi1, pi0)
  treatment <- stats::rbinom(n, 1L, treatment_probability)
  baseline_region <- sample(
    seq_along(specification$interval_probabilities),
    size = n,
    replace = TRUE,
    prob = specification$interval_probabilities
  )
  residual <- stats::runif(
    n,
    specification$interval_left[baseline_region],
    specification$interval_right[baseline_region]
  )

  boundary_treated_z1 <- boundary_x == 0L & treatment == 1L & z == 1L
  if (any(boundary_treated_z1)) {
    sampling_probability <- specification$interval_probabilities *
      specification$gamma1_boundary_values
    sampling_probability <- sampling_probability / sum(sampling_probability)
    region <- sample(
      seq_along(sampling_probability),
      size = sum(boundary_treated_z1),
      replace = TRUE,
      prob = sampling_probability
    )
    residual[boundary_treated_z1] <- stats::runif(
      sum(boundary_treated_z1),
      specification$interval_left[region],
      specification$interval_right[region]
    )
  }

  region <- gamma_boundary_two_stratum_region(residual, specification)
  gamma1 <- ifelse(
    boundary_x == 0L,
    specification$gamma1_boundary_values[region],
    1
  )
  data.frame(
    Y = specification$treatment_effect * treatment + residual,
    T_obs = treatment,
    Z = z,
    X = as.numeric(boundary_x),
    X_boundary = as.numeric(boundary_x == 0L),
    X_regular = as.numeric(boundary_x == 1L),
    outcome_region = factor(
      specification$interval_labels[region],
      levels = specification$interval_labels
    ),
    e_true = rep(1 - specification$probability_z1, n),
    pi0_true = pi0,
    pi1_true = pi1,
    gamma0_oracle = rep(1, n),
    gamma1_oracle = gamma1,
    stringsAsFactors = FALSE
  )
}

gamma_boundary_two_stratum_envelopes <- function(
    pi0,
    pi1,
    gamma,
    sensitivity_lambda,
    envelope_smoothing_epsilon = 0
) {
  Lambda <- sensitivity_lambda
  epsilon <- envelope_smoothing_epsilon
  candidate_lower <- (
    gamma * (pi1 + (1 - pi1) / Lambda) - pi0
  ) / (1 - pi0)
  candidate_upper <- (
    gamma * (pi1 + (1 - pi1) * Lambda) - pi0
  ) / (1 - pi0)
  hard_lower <- pmax(1 / Lambda, candidate_lower)
  hard_upper <- pmin(Lambda, candidate_upper)
  if (epsilon == 0) {
    lower <- hard_lower
    upper <- hard_upper
  } else {
    bandwidth <- epsilon / log(2)
    maximum <- pmax(1 / Lambda, candidate_lower)
    lower <- maximum + bandwidth * log(
      (
        exp((1 / Lambda - maximum) / bandwidth) +
          exp((candidate_lower - maximum) / bandwidth)
      ) / 2
    )
    minimum <- pmin(Lambda, candidate_upper)
    upper <- minimum - bandwidth * log(
      (
        exp(-(Lambda - minimum) / bandwidth) +
          exp(-(candidate_upper - minimum) / bandwidth)
      ) / 2
    )
  }
  list(
    lower = lower,
    upper = upper,
    hard_lower = hard_lower,
    hard_upper = hard_upper
  )
}

gamma_boundary_two_stratum_local_population <- function(
    pi0,
    pi1,
    gamma,
    specification,
    envelope_smoothing_epsilon
) {
  envelope <- gamma_boundary_two_stratum_envelopes(
    pi0, pi1, gamma, specification$sensitivity_lambda,
    envelope_smoothing_epsilon
  )
  left <- specification$interval_left
  right <- specification$interval_right
  probability <- specification$interval_probabilities
  density <- probability / (right - left)
  first_moment <- density * (right^2 - left^2) / 2
  difference <- envelope$upper - envelope$lower
  cumulative_difference <- cumsum(probability * difference)

  integrate_difference_to <- function(required_mass) {
    segment <- which(cumulative_difference >= required_mass - 1e-14)[1L]
    if (is.na(segment)) stop("No threshold solves the normalization equation.")
    mass_before <- if (segment == 1L) 0 else cumulative_difference[segment - 1L]
    local_probability <- (
      required_mass - mass_before
    ) / difference[segment]
    threshold <- left[segment] +
      local_probability / density[segment]
    moment_before <- if (segment == 1L) {
      0
    } else {
      sum(difference[seq_len(segment - 1L)] *
        first_moment[seq_len(segment - 1L)])
    }
    list(
      threshold = threshold,
      moment = moment_before + difference[segment] * density[segment] *
        (threshold^2 - left[segment]^2) / 2
    )
  }

  lower_integral <- integrate_difference_to(
    1 - sum(probability * envelope$lower)
  )
  upper_integral <- integrate_difference_to(
    sum(probability * envelope$upper) - 1
  )
  tilted_lower <- sum(envelope$lower * first_moment) + lower_integral$moment
  tilted_upper <- sum(envelope$upper * first_moment) - upper_integral$moment
  list(
    lower = specification$treatment_effect +
      pi0 * 0.5 + (1 - pi0) * tilted_lower,
    upper = specification$treatment_effect +
      pi0 * 0.5 + (1 - pi0) * tilted_upper,
    lower_threshold = specification$treatment_effect + lower_integral$threshold,
    upper_threshold = specification$treatment_effect + upper_integral$threshold,
    mean_lower_smoothing_gap = sum(
      probability * (envelope$hard_lower - envelope$lower)
    ),
    mean_upper_smoothing_gap = sum(
      probability * (envelope$upper - envelope$hard_upper)
    )
  )
}

gamma_boundary_two_stratum_population <- function(
    specification = make_gamma_boundary_two_stratum_spec(),
    envelope_smoothing_epsilon = 0
) {
  epsilon <- as.numeric(envelope_smoothing_epsilon)
  if (
    length(epsilon) != 1L || !is.finite(epsilon) || epsilon < 0 ||
      epsilon >= 1 / specification$sensitivity_lambda
  ) {
    stop("Invalid envelope_smoothing_epsilon.")
  }
  boundary <- gamma_boundary_two_stratum_local_population(
    specification$pi0["boundary"], specification$pi1["boundary"],
    specification$gamma1_boundary_values, specification, epsilon
  )
  regular <- gamma_boundary_two_stratum_local_population(
    specification$pi0["regular"], specification$pi1["regular"],
    specification$gamma1_regular_values, specification, epsilon
  )
  weight <- specification$x_probabilities
  data.frame(
    envelope_smoothing_epsilon = epsilon,
    mu1_lower = weight["boundary"] * boundary$lower +
      weight["regular"] * regular$lower,
    mu1_upper = weight["boundary"] * boundary$upper +
      weight["regular"] * regular$upper,
    boundary_lower_threshold = boundary$lower_threshold,
    boundary_upper_threshold = boundary$upper_threshold,
    regular_lower_threshold = regular$lower_threshold,
    regular_upper_threshold = regular$upper_threshold,
    mean_lower_smoothing_gap = weight["boundary"] *
      boundary$mean_lower_smoothing_gap + weight["regular"] *
      regular$mean_lower_smoothing_gap,
    mean_upper_smoothing_gap = weight["boundary"] *
      boundary$mean_upper_smoothing_gap + weight["regular"] *
      regular$mean_upper_smoothing_gap,
    stringsAsFactors = FALSE
  )
}

gamma_boundary_two_stratum_population_grid <- function(
    specification = make_gamma_boundary_two_stratum_spec(),
    epsilon_grid = c(0, 0.01)
) {
  epsilon_grid <- sort(unique(as.numeric(epsilon_grid)))
  if (
    length(epsilon_grid) < 1L || any(!is.finite(epsilon_grid)) ||
      any(epsilon_grid < 0) || !any(epsilon_grid == 0) ||
      any(epsilon_grid >= 1 / specification$sensitivity_lambda)
  ) {
    stop("epsilon_grid must include zero and contain valid values.")
  }
  output <- do.call(rbind, lapply(epsilon_grid, function(epsilon) {
    gamma_boundary_two_stratum_population(specification, epsilon)
  }))
  rownames(output) <- NULL
  hard <- output[output$envelope_smoothing_epsilon == 0, , drop = FALSE]
  output$lower_displacement_from_sharp <- output$mu1_lower - hard$mu1_lower
  output$upper_displacement_from_sharp <- output$mu1_upper - hard$mu1_upper
  output
}
