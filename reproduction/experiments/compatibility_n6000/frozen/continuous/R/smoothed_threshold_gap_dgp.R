# Observed-data IV law where conservative envelope smoothing resolves an exact
# nonunique threshold.
#
# X is uniform on four strata, Z is randomized, and
#
#   P(T = 1 | Z = 0) = 0.60,
#   P(T = 1 | Z = 1) = 0.80.
#
# For treated observations under Z = 0, the outcome is a 32:68 mixture of two
# narrow uniform intervals separated by a gap.  Under Z = 1 the corresponding
# mixture is 20:80.  Thus gamma_1 is 5/8 on the lower interval and 20/17 on
# the upper interval.  Controls have outcome zero and gamma_0 = 1.
#
# At Lambda = 4, the lower treated extremizer is normalized by assigning U to
# the lower interval and L to the upper interval.  Every threshold in the gap
# is therefore a population optimizer.  With four strata, all 2^4 choices of
# the two gap endpoints are distinct optimal threshold functions.  The
# conservative epsilon = 0.01 envelopes perturb the normalization equation,
# moving the optimizer into the upper interval and making it unique.

if (!exists("compute_projected_envelopes", mode = "function")) {
  if (!file.exists("R/sharp_bounds_plugin.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/sharp_bounds_plugin.R")
}

make_smoothed_threshold_gap_spec <- function(
    sensitivity_lambda = 4,
    strata = 4L,
    probability_z0 = 0.5,
    treatment_probability_z0 = 0.6,
    treatment_probability_z1 = 0.8,
    baseline_interval_probabilities = c(0.32, 0.68),
    gamma1_values = c(5 / 8, 20 / 17),
    lower_interval = c(-0.001, 0),
    upper_interval = c(1, 1.001),
    control_outcome = 0
) {
  strata <- as.integer(strata)
  interval_probability <- as.numeric(baseline_interval_probabilities)
  gamma1 <- as.numeric(gamma1_values)
  lower_interval <- as.numeric(lower_interval)
  upper_interval <- as.numeric(upper_interval)
  values <- c(
    sensitivity_lambda, probability_z0, treatment_probability_z0,
    treatment_probability_z1, interval_probability, gamma1,
    lower_interval, upper_interval, control_outcome
  )
  if (
    any(!is.finite(values)) || sensitivity_lambda <= 1 ||
      is.na(strata) || strata < 2L ||
      probability_z0 <= 0 || probability_z0 >= 1 ||
      treatment_probability_z0 <= 0 || treatment_probability_z0 >= 1 ||
      treatment_probability_z1 <= treatment_probability_z0 ||
      treatment_probability_z1 >= 1 ||
      length(interval_probability) != 2L ||
      any(interval_probability <= 0) ||
      abs(sum(interval_probability) - 1) > 1e-12 ||
      length(gamma1) != 2L || any(gamma1 <= 0) ||
      abs(sum(interval_probability * gamma1) - 1) > 1e-12 ||
      length(lower_interval) != 2L || length(upper_interval) != 2L ||
      lower_interval[1L] >= lower_interval[2L] ||
      upper_interval[1L] >= upper_interval[2L] ||
      lower_interval[2L] >= upper_interval[1L]
  ) {
    stop("Invalid smoothed-threshold-gap specification.")
  }

  nuisance <- list(
    e = rep(probability_z0, 2L),
    pi0 = rep(treatment_probability_z0, 2L),
    pi1 = rep(treatment_probability_z1, 2L)
  )
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = gamma1,
    t_value = 1L,
    Lambda0 = sensitivity_lambda,
    Lambda1 = sensitivity_lambda,
    envelope_smoothing_epsilon = 0,
    gamma_compatibility_action = "stop"
  )
  hard_normalization <- interval_probability[1L] * envelope$U[1L] +
    interval_probability[2L] * envelope$L[2L]
  if (
    any(envelope$gamma_incompatible_original) ||
      any(envelope$L >= envelope$U) ||
      abs(hard_normalization - 1) > 1e-10
  ) {
    stop("The sharp lower endpoint does not have the requested gap tie.")
  }

  specification <- list(
    name = "smoothed_threshold_gap_complete_binary",
    sensitivity_lambda = as.numeric(sensitivity_lambda),
    strata = strata,
    candidates = 2^strata,
    probability_z0 = as.numeric(probability_z0),
    probability_z1 = 1 - as.numeric(probability_z0),
    pi0 = as.numeric(treatment_probability_z0),
    pi1 = as.numeric(treatment_probability_z1),
    overlap_constant = min(
      probability_z0, 1 - probability_z0,
      treatment_probability_z0, 1 - treatment_probability_z0,
      treatment_probability_z1, 1 - treatment_probability_z1
    ),
    marginal_first_stage = treatment_probability_z1 -
      treatment_probability_z0,
    interval_labels = c("lower", "upper"),
    interval_probabilities = interval_probability,
    interval_left = c(lower_interval[1L], upper_interval[1L]),
    interval_right = c(lower_interval[2L], upper_interval[2L]),
    gamma1_values = gamma1,
    gamma0_value = 1,
    control_outcome = as.numeric(control_outcome),
    treatment_effect = NA_real_,
    covariate_design = "four_stratum_complete_binary",
    outcome_geometry = "two_interval_gap",
    separated_cluster_width = max(diff(lower_interval), diff(upper_interval)),
    hard_gap_normalization = hard_normalization
  )
  class(specification) <- c(
    "smoothed_threshold_gap_specification", "list"
  )
  specification
}

smoothed_threshold_gap_x_columns <- function(
    specification = make_smoothed_threshold_gap_spec()
) {
  paste0("X", seq_len(specification$strata - 1L))
}

smoothed_threshold_gap_add_x_columns <- function(data, specification) {
  for (stratum in seq_len(specification$strata - 1L)) {
    data[[paste0("X", stratum)]] <- as.numeric(data$X == stratum)
  }
  data
}

generate_smoothed_threshold_gap_data <- function(
    n,
    specification = make_smoothed_threshold_gap_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  x <- sample.int(specification$strata, n, replace = TRUE)
  z <- stats::rbinom(n, 1L, specification$probability_z1)
  pi0 <- rep(specification$pi0, n)
  pi1 <- rep(specification$pi1, n)
  probability_treatment <- ifelse(z == 1L, pi1, pi0)
  treatment <- stats::rbinom(n, 1L, probability_treatment)

  y <- rep(specification$control_outcome, n)
  region <- rep(NA_integer_, n)
  treated <- treatment == 1L
  lower_probability <- ifelse(
    z == 0L,
    specification$interval_probabilities[1L],
    specification$interval_probabilities[1L] *
      specification$gamma1_values[1L]
  )
  region[treated] <- ifelse(
    stats::runif(sum(treated)) < lower_probability[treated], 1L, 2L
  )
  y[treated] <- stats::runif(
    sum(treated),
    specification$interval_left[region[treated]],
    specification$interval_right[region[treated]]
  )

  gamma1 <- rep(specification$gamma1_values[1L], n)
  gamma1[treated] <- specification$gamma1_values[region[treated]]
  output <- data.frame(
    Y = y,
    T_obs = treatment,
    Z = z,
    X = x,
    outcome_region = factor(
      ifelse(treated, specification$interval_labels[region], "control"),
      levels = c("control", specification$interval_labels)
    ),
    e_true = rep(specification$probability_z0, n),
    pi0_true = pi0,
    pi1_true = pi1,
    gamma0_oracle = rep(1, n),
    gamma1_oracle = gamma1,
    stringsAsFactors = FALSE
  )
  smoothed_threshold_gap_add_x_columns(output, specification)
}

smoothed_threshold_gap_arm_population <- function(
    specification = make_smoothed_threshold_gap_spec(),
    envelope_smoothing_epsilon = 0
) {
  epsilon <- as.numeric(envelope_smoothing_epsilon)
  if (
    length(epsilon) != 1L || !is.finite(epsilon) || epsilon < 0 ||
      epsilon >= 1 / specification$sensitivity_lambda
  ) {
    stop("Invalid envelope_smoothing_epsilon.")
  }
  probability <- specification$interval_probabilities
  left <- specification$interval_left
  right <- specification$interval_right
  density <- probability / (right - left)
  first_moment <- density * (right^2 - left^2) / 2
  nuisance <- list(
    e = rep(specification$probability_z0, 2L),
    pi0 = rep(specification$pi0, 2L),
    pi1 = rep(specification$pi1, 2L)
  )
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = specification$gamma1_values,
    t_value = 1L,
    Lambda0 = specification$sensitivity_lambda,
    Lambda1 = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = epsilon,
    gamma_compatibility_action = "stop"
  )
  difference <- envelope$D
  cumulative_difference <- cumsum(probability * difference)

  integrate_difference_to <- function(required_mass) {
    segment <- which(cumulative_difference >= required_mass - 1e-14)[1L]
    if (is.na(segment)) stop("No threshold solves the normalization equation.")
    mass_before <- if (segment == 1L) 0 else {
      cumulative_difference[segment - 1L]
    }
    local_probability <- (required_mass - mass_before) / difference[segment]
    local_probability <- pmin(pmax(local_probability, 0), probability[segment])
    threshold <- left[segment] + local_probability / density[segment]
    moment_before <- if (segment == 1L) 0 else {
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
    1 - sum(probability * envelope$L)
  )
  upper_integral <- integrate_difference_to(
    sum(probability * envelope$U) - 1
  )
  tilted_lower <- sum(envelope$L * first_moment) + lower_integral$moment
  tilted_upper <- sum(envelope$U * first_moment) - upper_integral$moment
  observed_mean <- sum(first_moment)
  r0 <- specification$pi0

  data.frame(
    envelope_smoothing_epsilon = epsilon,
    mu1_lower = r0 * observed_mean + (1 - r0) * tilted_lower,
    mu1_upper = r0 * observed_mean + (1 - r0) * tilted_upper,
    mu0_lower = specification$control_outcome,
    mu0_upper = specification$control_outcome,
    ate_lower = r0 * observed_mean + (1 - r0) * tilted_lower -
      specification$control_outcome,
    ate_upper = r0 * observed_mean + (1 - r0) * tilted_upper -
      specification$control_outcome,
    lower_threshold = lower_integral$threshold,
    upper_threshold = upper_integral$threshold,
    gap_lower = right[1L],
    gap_upper = left[2L],
    gap_lower_normalization_residual =
      probability[1L] * envelope$U[1L] +
      probability[2L] * envelope$L[2L] - 1,
    stringsAsFactors = FALSE
  )
}

smoothed_threshold_gap_population_grid <- function(
    specification = make_smoothed_threshold_gap_spec(),
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
    smoothed_threshold_gap_arm_population(specification, epsilon)
  }))
  rownames(output) <- NULL
  hard <- output[output$envelope_smoothing_epsilon == 0, , drop = FALSE]
  output$lower_displacement_from_sharp <- output$mu1_lower - hard$mu1_lower
  output$upper_displacement_from_sharp <- output$mu1_upper - hard$mu1_upper
  output
}
