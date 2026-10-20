# Continuous exceptional-law DGP for envelope smoothing.
#
# X is uniform on [-1, 1].  On the interval [-1, tie_cutoff], the two
# instrument-specific treatment probabilities are exactly equal, so both
# treatment-one envelope candidates meet their Lambda caps.  Above the
# cutoff, pi1(X) rises smoothly while pi0(X) remains constant.  The default
# exceptional set therefore has probability 0.75, while the marginal first
# stage is 0.10 and every treatment probability is bounded below by 0.10.
#
# Outcomes follow
#
#   Y = treatment_effect * T + E,  E ~ Uniform(0, 1),
#
# independently of Z conditional on (T, X).  Hence gamma_0 = gamma_1 = 1.

make_continuous_envelope_exceptional_spec <- function(
    sensitivity_lambda = 4,
    probability_z1 = 0.5,
    treatment_probability_z0 = 0.1,
    maximum_treatment_probability_z1 = 0.9,
    tie_cutoff = 0.5,
    treatment_effect = 0.75
) {
  values <- as.numeric(c(
    sensitivity_lambda, probability_z1, treatment_probability_z0,
    maximum_treatment_probability_z1, tie_cutoff, treatment_effect
  ))
  if (
    any(!is.finite(values)) || sensitivity_lambda <= 1 ||
      probability_z1 <= 0 || probability_z1 >= 1 ||
      treatment_probability_z0 <= 0 || treatment_probability_z0 >= 0.5 ||
      maximum_treatment_probability_z1 <= treatment_probability_z0 ||
      maximum_treatment_probability_z1 >= 1 ||
      tie_cutoff <= -1 || tie_cutoff >= 1
  ) {
    stop("Invalid continuous exceptional-law specification.")
  }

  tie_mass <- (tie_cutoff + 1) / 2
  regular_mass <- 1 - tie_mass
  marginal_first_stage <- regular_mass *
    (maximum_treatment_probability_z1 - treatment_probability_z0) / 2
  specification <- list(
    name = "continuous_positive_mass_envelope_ties",
    sensitivity_lambda = as.numeric(sensitivity_lambda),
    probability_z1 = as.numeric(probability_z1),
    pi0 = as.numeric(treatment_probability_z0),
    pi1 = c(
      tie = as.numeric(treatment_probability_z0),
      maximum = as.numeric(maximum_treatment_probability_z1)
    ),
    tie_cutoff = as.numeric(tie_cutoff),
    tie_mass = tie_mass,
    overlap_constant = min(
      probability_z1, 1 - probability_z1,
      treatment_probability_z0, 1 - treatment_probability_z0,
      maximum_treatment_probability_z1,
      1 - maximum_treatment_probability_z1
    ),
    marginal_first_stage = marginal_first_stage,
    interval_probabilities = c(uniform_outcome = 1),
    gamma1_values = c(constant = 1),
    gamma0_value = 1,
    treatment_effect = as.numeric(treatment_effect),
    covariate_design = "continuous_exceptional",
    outcome_geometry = "uniform",
    separated_cluster_width = NA_real_
  )
  class(specification) <- c(
    "continuous_envelope_exceptional_specification", "list"
  )
  specification
}

continuous_envelope_exceptional_pi1 <- function(
    x,
    specification = make_continuous_envelope_exceptional_spec()
) {
  x <- as.numeric(x)
  if (any(!is.finite(x)) || any(x < -1) || any(x > 1)) {
    stop("x must be finite and lie in [-1, 1].")
  }
  pi0 <- unname(specification$pi0)
  pi1_maximum <- unname(specification$pi1["maximum"])
  scaled <- pmax(
    0,
    (x - specification$tie_cutoff) / (1 - specification$tie_cutoff)
  )
  smoothstep <- 3 * scaled^2 - 2 * scaled^3
  pi0 + (pi1_maximum - pi0) * smoothstep
}

generate_continuous_envelope_exceptional_data <- function(
    n,
    specification = make_continuous_envelope_exceptional_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  x <- stats::runif(n, -1, 1)
  pi0 <- rep(unname(specification$pi0), n)
  pi1 <- continuous_envelope_exceptional_pi1(x, specification)
  z <- stats::rbinom(n, 1L, specification$probability_z1)
  probability_treatment <- ifelse(z == 1L, pi1, pi0)
  treatment <- stats::rbinom(n, 1L, probability_treatment)
  residual <- stats::runif(n)

  data.frame(
    Y = specification$treatment_effect * treatment + residual,
    T_obs = treatment,
    Z = z,
    X = x,
    X_exceptional = as.numeric(x <= specification$tie_cutoff),
    e_true = rep(1 - specification$probability_z1, n),
    pi0_true = pi0,
    pi1_true = pi1,
    gamma0_oracle = rep(1, n),
    gamma1_oracle = rep(1, n),
    stringsAsFactors = FALSE
  )
}

continuous_envelope_exceptional_local_mu1 <- function(
    pi0,
    pi1,
    endpoint = c("lower", "upper"),
    specification = make_continuous_envelope_exceptional_spec(),
    envelope_smoothing_epsilon = 0
) {
  endpoint <- match.arg(endpoint)
  pi0 <- rep_len(as.numeric(pi0), length(pi1))
  pi1 <- as.numeric(pi1)
  Lambda <- specification$sensitivity_lambda
  epsilon <- as.numeric(envelope_smoothing_epsilon)
  if (
    any(!is.finite(c(pi0, pi1, epsilon))) || any(pi0 <= 0) ||
      any(pi0 >= 1) || any(pi1 <= 0) || any(pi1 >= 1) ||
      length(epsilon) != 1L || epsilon < 0 || epsilon >= 1 / Lambda
  ) {
    stop("Invalid local population-evaluation inputs.")
  }

  candidate_lower <- (
    pi1 + (1 - pi1) / Lambda - pi0
  ) / (1 - pi0)
  candidate_upper <- (
    pi1 + (1 - pi1) * Lambda - pi0
  ) / (1 - pi0)
  if (epsilon == 0) {
    lower <- pmax(1 / Lambda, candidate_lower)
    upper <- pmin(Lambda, candidate_upper)
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
  if (any(lower >= upper)) stop("The local envelopes cross.")

  treatment_effect <- specification$treatment_effect
  outcome_mean <- treatment_effect + 0.5
  if (endpoint == "lower") {
    quantile_residual <- (1 - lower) / (upper - lower)
    lower_partial_moment <- treatment_effect * quantile_residual +
      quantile_residual^2 / 2
    tilted_mean <- lower * outcome_mean +
      (upper - lower) * lower_partial_moment
  } else {
    quantile_residual <- (upper - 1) / (upper - lower)
    lower_partial_moment <- treatment_effect * quantile_residual +
      quantile_residual^2 / 2
    tilted_mean <- upper * outcome_mean +
      (lower - upper) * lower_partial_moment
  }
  pi0 * outcome_mean + (1 - pi0) * tilted_mean
}

continuous_envelope_exceptional_population <- function(
    specification = make_continuous_envelope_exceptional_spec(),
    envelope_smoothing_epsilon = 0,
    integration_points = 200000L
) {
  integration_points <- as.integer(integration_points)
  if (
    length(integration_points) != 1L || is.na(integration_points) ||
      integration_points < 1000L
  ) {
    stop("integration_points must be an integer of at least 1000.")
  }
  x <- -1 + (seq_len(integration_points) - 0.5) * 2 / integration_points
  pi1 <- continuous_envelope_exceptional_pi1(x, specification)
  lower <- continuous_envelope_exceptional_local_mu1(
    specification$pi0, pi1, "lower", specification,
    envelope_smoothing_epsilon
  )
  upper <- continuous_envelope_exceptional_local_mu1(
    specification$pi0, pi1, "upper", specification,
    envelope_smoothing_epsilon
  )
  data.frame(
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    mu1_lower = mean(lower),
    mu1_upper = mean(upper),
    mu1_width = mean(upper) - mean(lower),
    tie_mass = specification$tie_mass,
    marginal_first_stage = specification$marginal_first_stage,
    stringsAsFactors = FALSE
  )
}

continuous_envelope_exceptional_population_grid <- function(
    specification = make_continuous_envelope_exceptional_spec(),
    epsilon_grid = c(0, 0.01),
    integration_points = 200000L
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
    continuous_envelope_exceptional_population(
      specification, epsilon, integration_points
    )
  }))
  rownames(output) <- NULL
  hard <- output[output$envelope_smoothing_epsilon == 0, , drop = FALSE]
  output$lower_displacement_from_sharp <- output$mu1_lower - hard$mu1_lower
  output$upper_displacement_from_sharp <- output$mu1_upper - hard$mu1_upper
  output
}
