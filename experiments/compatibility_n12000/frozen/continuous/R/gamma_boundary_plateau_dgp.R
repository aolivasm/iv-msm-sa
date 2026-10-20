# Observed-data law with positive-mass density-ratio plateaus at the exact
# switching boundaries of the sharp marginal-sensitivity envelopes.
#
# The design targets nonregularity in the estimated density ratio rather than
# in the treatment regressions.  There are no treatment-probability ties:
#
#   P(T = 1 | Z = 0, X) = pi0,
#   P(T = 1 | Z = 1, X) = pi1,
#
# with defaults pi0 = 0.10 and pi1 = 0.30.  X is independent noise, so all
# nuisance functions are nevertheless estimated as functions of X in the
# fully estimated experiment.
#
# For treated observations, the baseline residual law under Z = 0 is
# Uniform(0, 1).  Its density ratio under Z = 1 is piecewise constant.  On a
# middle interval it equals the exact lower-envelope switching value and on
# the upper interval it equals the exact upper-envelope switching value.  The
# value on the remaining lower interval is determined by normalization.  For
# the defaults, the three baseline interval probabilities are 0.20, 0.60,
# and 0.20.  Untreated outcomes have gamma_0 = 1.

make_gamma_boundary_plateau_spec <- function(
    sensitivity_lambda = 4,
    probability_z1 = 0.5,
    treatment_probability_z0 = 0.1,
    treatment_probability_z1 = 0.3,
    baseline_interval_probabilities = c(0.2, 0.6, 0.2),
    treatment_effect = 0.75,
    covariate_design = c("independent_uniform", "constant")
) {
  covariate_design <- match.arg(covariate_design)
  Lambda <- as.numeric(sensitivity_lambda)
  probability_z1 <- as.numeric(probability_z1)
  pi0 <- as.numeric(treatment_probability_z0)
  pi1 <- as.numeric(treatment_probability_z1)
  interval_probabilities <- as.numeric(baseline_interval_probabilities)
  treatment_effect <- as.numeric(treatment_effect)

  if (
    length(Lambda) != 1L || !is.finite(Lambda) || Lambda <= 1 ||
      length(probability_z1) != 1L || !is.finite(probability_z1) ||
      probability_z1 <= 0 || probability_z1 >= 1 ||
      length(pi0) != 1L || length(pi1) != 1L ||
      any(!is.finite(c(pi0, pi1))) || pi0 <= 0 || pi1 >= 1 ||
      pi1 <= pi0 ||
      length(interval_probabilities) != 3L ||
      any(!is.finite(interval_probabilities)) ||
      any(interval_probabilities <= 0) ||
      abs(sum(interval_probabilities) - 1) > 1e-12 ||
      length(treatment_effect) != 1L || !is.finite(treatment_effect)
  ) {
    stop("Invalid gamma-boundary plateau specification.")
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

  # Interval order: regular lower tail, lower-envelope tie plateau,
  # upper-envelope tie plateau.
  regular_gamma <- (
    1 - interval_probabilities[2L] * lower_switch -
      interval_probabilities[3L] * upper_switch
  ) / interval_probabilities[1L]
  gamma_values <- c(regular_gamma, lower_switch, upper_switch)

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
  compatible_lower_t0 <- (
    pi0 / Lambda + 1 - pi0
  ) / (
    pi1 * Lambda + 1 - pi1
  )
  compatible_upper_t0 <- (
    Lambda * pi0 + 1 - pi0
  ) / (
    pi1 / Lambda + 1 - pi1
  )

  if (
    any(!is.finite(gamma_values)) || any(gamma_values <= 0) ||
      any(gamma_values < compatible_lower_t1) ||
      any(gamma_values > compatible_upper_t1) ||
      1 < compatible_lower_t0 || 1 > compatible_upper_t0
  ) {
    stop(
      paste0(
        "The requested plateau density ratios are not compatible with ",
        "the sensitivity model."
      )
    )
  }

  interval_boundaries <- c(0, cumsum(interval_probabilities))
  specification <- list(
    name = "positive_mass_gamma_envelope_boundary_plateaus",
    sensitivity_lambda = Lambda,
    probability_z1 = probability_z1,
    pi0 = pi0,
    pi1 = pi1,
    marginal_first_stage = pi1 - pi0,
    overlap_constant = min(
      probability_z1, 1 - probability_z1,
      pi0, 1 - pi0, pi1, 1 - pi1
    ),
    interval_labels = c("regular", "lower_tie", "upper_tie"),
    interval_probabilities = interval_probabilities,
    interval_boundaries = interval_boundaries,
    gamma1_values = gamma_values,
    gamma1_lower_switch = lower_switch,
    gamma1_upper_switch = upper_switch,
    gamma1_compatible_interval = c(
      lower = compatible_lower_t1,
      upper = compatible_upper_t1
    ),
    gamma0_value = 1,
    gamma0_compatible_interval = c(
      lower = compatible_lower_t0,
      upper = compatible_upper_t0
    ),
    treatment_effect = treatment_effect,
    covariate_design = covariate_design
  )
  class(specification) <- c(
    "gamma_boundary_plateau_specification", "list"
  )
  specification
}

gamma_boundary_plateau_region <- function(
    residual,
    specification = make_gamma_boundary_plateau_spec()
) {
  residual <- as.numeric(residual)
  if (any(!is.finite(residual))) {
    stop("residual must be finite.")
  }
  boundaries <- specification$interval_boundaries
  # Include the right endpoint in the final interval.  Exact boundary draws
  # have probability zero under the generator, but this convention keeps the
  # function deterministic for hand-constructed test values.
  region <- findInterval(
    pmin(1, pmax(0, residual)),
    boundaries,
    rightmost.closed = TRUE,
    all.inside = TRUE
  )
  pmin(length(specification$interval_probabilities), pmax(1L, region))
}

gamma_boundary_plateau_gamma1 <- function(
    residual,
    specification = make_gamma_boundary_plateau_spec()
) {
  region <- gamma_boundary_plateau_region(residual, specification)
  specification$gamma1_values[region]
}

generate_gamma_boundary_plateau_data <- function(
    n,
    specification = make_gamma_boundary_plateau_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  x <- if (identical(specification$covariate_design, "constant")) {
    rep(0, n)
  } else {
    stats::runif(n, -1, 1)
  }
  z <- stats::rbinom(n, 1L, specification$probability_z1)
  treatment_probability <- ifelse(
    z == 1L, specification$pi1, specification$pi0
  )
  treatment <- stats::rbinom(n, 1L, treatment_probability)
  residual <- stats::runif(n)

  treated_z1 <- treatment == 1L & z == 1L
  if (any(treated_z1)) {
    sampling_probability <- specification$interval_probabilities *
      specification$gamma1_values
    sampling_probability <- sampling_probability / sum(sampling_probability)
    region <- sample(
      seq_along(sampling_probability),
      size = sum(treated_z1),
      replace = TRUE,
      prob = sampling_probability
    )
    lower <- specification$interval_boundaries[region]
    upper <- specification$interval_boundaries[region + 1L]
    residual[treated_z1] <- stats::runif(sum(treated_z1), lower, upper)
  }

  region <- gamma_boundary_plateau_region(residual, specification)
  data.frame(
    Y = specification$treatment_effect * treatment + residual,
    T_obs = treatment,
    Z = z,
    X = x,
    outcome_region = factor(
      specification$interval_labels[region],
      levels = specification$interval_labels
    ),
    e_true = rep(1 - specification$probability_z1, n),
    pi0_true = rep(specification$pi0, n),
    pi1_true = rep(specification$pi1, n),
    gamma0_oracle = rep(specification$gamma0_value, n),
    gamma1_oracle = specification$gamma1_values[region],
    stringsAsFactors = FALSE
  )
}

gamma_boundary_plateau_envelopes <- function(
    specification = make_gamma_boundary_plateau_spec(),
    envelope_smoothing_epsilon = 0
) {
  epsilon <- as.numeric(envelope_smoothing_epsilon)
  if (
    length(epsilon) != 1L || !is.finite(epsilon) || epsilon < 0 ||
      epsilon >= 1 / specification$sensitivity_lambda
  ) {
    stop("Invalid envelope_smoothing_epsilon.")
  }

  Lambda <- specification$sensitivity_lambda
  pi0 <- specification$pi0
  pi1 <- specification$pi1
  gamma <- specification$gamma1_values
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
    bandwidth <- 0
  } else {
    bandwidth <- epsilon / log(2)
    stable_shifted_softmax <- function(first, second) {
      maximum <- pmax(first, second)
      maximum + bandwidth * log(
        (
          exp((first - maximum) / bandwidth) +
            exp((second - maximum) / bandwidth)
        ) / 2
      )
    }
    stable_shifted_softmin <- function(first, second) {
      minimum <- pmin(first, second)
      minimum - bandwidth * log(
        (
          exp(-(first - minimum) / bandwidth) +
            exp(-(second - minimum) / bandwidth)
        ) / 2
      )
    }
    lower <- stable_shifted_softmax(1 / Lambda, candidate_lower)
    upper <- stable_shifted_softmin(Lambda, candidate_upper)
  }

  if (any(!is.finite(lower)) || any(!is.finite(upper)) || any(lower >= upper)) {
    stop("The plateau specification produced invalid envelopes.")
  }
  data.frame(
    region = specification$interval_labels,
    gamma1 = gamma,
    candidate_lower = candidate_lower,
    candidate_upper = candidate_upper,
    hard_lower = hard_lower,
    hard_upper = hard_upper,
    lower = lower,
    upper = upper,
    lower_smoothing_gap = hard_lower - lower,
    upper_smoothing_gap = upper - hard_upper,
    envelope_smoothing_epsilon = epsilon,
    envelope_smoothing_bandwidth = bandwidth,
    stringsAsFactors = FALSE
  )
}

gamma_boundary_plateau_population <- function(
    specification = make_gamma_boundary_plateau_spec(),
    envelope_smoothing_epsilon = 0
) {
  envelope <- gamma_boundary_plateau_envelopes(
    specification, envelope_smoothing_epsilon
  )
  boundaries <- specification$interval_boundaries
  left <- boundaries[-length(boundaries)]
  right <- boundaries[-1L]
  probability <- right - left
  first_moment <- (right^2 - left^2) / 2
  difference <- envelope$upper - envelope$lower
  cumulative_difference <- cumsum(probability * difference)

  integrate_difference_to <- function(required_mass) {
    segment <- which(cumulative_difference >= required_mass - 1e-14)[1L]
    if (is.na(segment)) {
      stop("The normalization equation has no solution.")
    }
    mass_before <- if (segment == 1L) 0 else cumulative_difference[segment - 1L]
    threshold <- left[segment] +
      (required_mass - mass_before) / difference[segment]
    moment_before <- if (segment == 1L) {
      0
    } else {
      sum(difference[seq_len(segment - 1L)] *
        first_moment[seq_len(segment - 1L)])
    }
    moment <- moment_before + difference[segment] *
      (threshold^2 - left[segment]^2) / 2
    list(threshold = threshold, moment = moment, segment = segment)
  }

  lower_required <- 1 - sum(probability * envelope$lower)
  upper_required <- sum(probability * envelope$upper) - 1
  lower_integral <- integrate_difference_to(lower_required)
  upper_integral <- integrate_difference_to(upper_required)
  tilted_lower <- sum(envelope$lower * first_moment) + lower_integral$moment
  tilted_upper <- sum(envelope$upper * first_moment) - upper_integral$moment
  baseline_mean <- 0.5
  pi0 <- specification$pi0

  data.frame(
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    mu1_lower = specification$treatment_effect +
      pi0 * baseline_mean + (1 - pi0) * tilted_lower,
    mu1_upper = specification$treatment_effect +
      pi0 * baseline_mean + (1 - pi0) * tilted_upper,
    lower_threshold = specification$treatment_effect +
      lower_integral$threshold,
    upper_threshold = specification$treatment_effect +
      upper_integral$threshold,
    mean_lower_smoothing_gap = sum(
      probability * envelope$lower_smoothing_gap
    ),
    mean_upper_smoothing_gap = sum(
      probability * envelope$upper_smoothing_gap
    ),
    stringsAsFactors = FALSE
  )
}

gamma_boundary_plateau_population_grid <- function(
    specification = make_gamma_boundary_plateau_spec(),
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
    gamma_boundary_plateau_population(specification, epsilon)
  }))
  rownames(output) <- NULL
  hard <- output[output$envelope_smoothing_epsilon == 0, , drop = FALSE]
  output$lower_displacement_from_sharp <- output$mu1_lower - hard$mu1_lower
  output$upper_displacement_from_sharp <- output$mu1_upper - hard$mu1_upper
  output
}
