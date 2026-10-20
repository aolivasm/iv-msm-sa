# Continuous exceptional-law analogue for the dual-MSM envelope kink.
#
# This design adapts the core geometry of the continuous exceptional-law
# simulation in Luedtke and van der Laan (2016).  X is uniform on [-1, 1].
# On the interval [-1/2, 1/3], pi1(X) = pi0(X) and gamma_t(Y, X) = 1, so both
# max/min branches defining the dual-MSM envelope tie exactly on a set of
# positive probability.  Outside that interval, the first stage and a mild
# outcome-density tilt increase linearly with distance from the tie region.
# Outcomes remain continuous, so threshold optimizers are unique: the only
# intended nonregularity is the envelope branch choice.

if (!exists("observed_dgp_rtruncnorm", mode = "function")) {
  if (!file.exists("R/observed_data_dgp.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/observed_data_dgp.R")
}

make_continuous_envelope_tie_spec <- function(
    sensitivity_lambda = 3,
    tie_interval = c(-0.5, 1 / 3),
    probability_z1 = 0.5,
    pi0_intercept = 0.45,
    pi0_slope = 0.05,
    maximum_first_stage = 0.35,
    local_tie_first_stage = 0,
    outcome_treatment_shift = 0.75,
    outcome_x_slope = 0.3,
    comparison_tilt_t0 = -0.075,
    comparison_tilt_t1 = 0.075,
    residual_lower = -2,
    residual_upper = 2
) {
  values <- c(
    sensitivity_lambda, tie_interval, probability_z1,
    pi0_intercept, pi0_slope, maximum_first_stage,
    local_tie_first_stage,
    outcome_treatment_shift, outcome_x_slope,
    comparison_tilt_t0, comparison_tilt_t1,
    residual_lower, residual_upper
  )
  if (any(!is.finite(values)) || length(sensitivity_lambda) != 1L ||
      sensitivity_lambda <= 1) {
    stop("All design parameters must be finite and Lambda must exceed one.")
  }
  if (length(tie_interval) != 2L || tie_interval[1L] <= -1 ||
      tie_interval[2L] >= 1 || tie_interval[1L] >= tie_interval[2L]) {
    stop("tie_interval must be a nonempty strict subinterval of [-1, 1].")
  }
  if (probability_z1 <= 0 || probability_z1 >= 1 ||
      maximum_first_stage <= 0 || local_tie_first_stage < 0) {
    stop("Instrument probabilities and strength are invalid.")
  }
  endpoint_pi0 <- pi0_intercept + pi0_slope * c(-1, 1)
  endpoint_pi1 <- endpoint_pi0 +
    max(maximum_first_stage, local_tie_first_stage)
  if (min(endpoint_pi0) <= 0.05 || max(endpoint_pi1) >= 0.95) {
    stop("The requested treatment probabilities violate the 0.05 overlap margin.")
  }
  if (residual_lower >= residual_upper) {
    stop("The residual support must be nonempty.")
  }

  left_probability <- (tie_interval[1L] + 1) / 2
  right_probability <- (1 - tie_interval[2L]) / 2
  tie_probability <- diff(tie_interval) / 2
  marginal_first_stage <- maximum_first_stage *
    (left_probability + right_probability) / 2 +
    local_tie_first_stage * tie_probability

  list(
    name = "continuous_envelope_tie_exceptional",
    dimension = 1L,
    sensitivity_lambda = sensitivity_lambda,
    tie_interval = tie_interval,
    probability_z1 = probability_z1,
    pi0_intercept = pi0_intercept,
    pi0_slope = pi0_slope,
    maximum_first_stage = maximum_first_stage,
    local_tie_first_stage = local_tie_first_stage,
    outcome_treatment_shift = outcome_treatment_shift,
    outcome_x_slope = outcome_x_slope,
    comparison_tilt_t0 = comparison_tilt_t0,
    comparison_tilt_t1 = comparison_tilt_t1,
    residual_lower = residual_lower,
    residual_upper = residual_upper,
    population_tie_probability = if (local_tie_first_stage == 0) {
      tie_probability
    } else 0,
    population_local_to_tie_probability = tie_probability,
    marginal_first_stage = marginal_first_stage
  )
}

make_local_to_tie_envelope_spec <- function(
    sample_size,
    local_coefficient = 0.15,
    sensitivity_lambda = 3,
    ...
) {
  sample_size <- as.integer(sample_size)
  local_coefficient <- as.numeric(local_coefficient)
  if (length(sample_size) != 1L || is.na(sample_size) || sample_size < 1L ||
      length(local_coefficient) != 1L || !is.finite(local_coefficient) ||
      local_coefficient <= 0) {
    stop("sample_size and local_coefficient must be positive finite scalars.")
  }
  specification <- make_continuous_envelope_tie_spec(
    sensitivity_lambda = sensitivity_lambda,
    local_tie_first_stage = local_coefficient / sqrt(sample_size),
    ...
  )
  specification$name <- "continuous_envelope_local_to_tie"
  specification$nominal_sample_size <- sample_size
  specification$local_coefficient <- local_coefficient
  specification
}

continuous_envelope_tie_validate_x <- function(x) {
  if (is.data.frame(x)) {
    if (!"X" %in% names(x)) stop("x data frame must contain X.")
    x <- x$X
  }
  x <- as.numeric(x)
  if (any(!is.finite(x)) || any(x < -1) || any(x > 1)) {
    stop("X must be finite and lie in [-1, 1].")
  }
  x
}

continuous_envelope_tie_shape <- function(x, specification) {
  x <- continuous_envelope_tie_validate_x(x)
  lower <- specification$tie_interval[1L]
  upper <- specification$tie_interval[2L]
  ifelse(
    x < lower,
    (lower - x) / (lower + 1),
    ifelse(x > upper, (x - upper) / (1 - upper), 0)
  )
}

continuous_envelope_tie_region <- function(x, specification) {
  x <- continuous_envelope_tie_validate_x(x)
  x >= specification$tie_interval[1L] &
    x <= specification$tie_interval[2L]
}

continuous_envelope_tie_nuisances <- function(x, specification) {
  x <- continuous_envelope_tie_validate_x(x)
  pi0 <- specification$pi0_intercept + specification$pi0_slope * x
  pi1 <- pi0 + specification$maximum_first_stage *
    continuous_envelope_tie_shape(x, specification) +
    specification$local_tie_first_stage * as.numeric(
      continuous_envelope_tie_region(x, specification)
    )
  list(
    e = rep(1 - specification$probability_z1, length(x)),
    pi0 = pi0,
    pi1 = pi1
  )
}

continuous_envelope_tie_location <- function(x, treatment, specification) {
  x <- continuous_envelope_tie_validate_x(x)
  treatment <- rep_len(as.numeric(treatment), length(x))
  if (any(!treatment %in% 0:1)) {
    stop("treatment must contain only zeroes and ones.")
  }
  treatment * specification$outcome_treatment_shift +
    specification$outcome_x_slope * x
}

continuous_envelope_tie_tilt <- function(x, treatment, specification) {
  x <- continuous_envelope_tie_validate_x(x)
  treatment <- rep_len(as.numeric(treatment), length(x))
  if (any(!treatment %in% 0:1)) {
    stop("treatment must contain only zeroes and ones.")
  }
  arm_tilt <- ifelse(
    treatment == 1,
    specification$comparison_tilt_t1,
    specification$comparison_tilt_t0
  )
  arm_tilt * continuous_envelope_tie_shape(x, specification)
}

continuous_envelope_tie_gamma <- function(y, x, treatment, specification) {
  x <- continuous_envelope_tie_validate_x(x)
  y <- rep_len(as.numeric(y), length(x))
  treatment <- rep_len(as.numeric(treatment), length(x))
  location <- continuous_envelope_tie_location(x, treatment, specification)
  tilt <- continuous_envelope_tie_tilt(x, treatment, specification)
  residual <- y - location
  lower <- specification$residual_lower
  upper <- specification$residual_upper
  base_normalizer <- stats::pnorm(upper) - stats::pnorm(lower)
  tilted_normalizer <- stats::pnorm(upper, mean = tilt) -
    stats::pnorm(lower, mean = tilt)
  exp(tilt * residual - tilt^2 / 2) *
    base_normalizer / tilted_normalizer
}

generate_continuous_envelope_tie_data <- function(
    n,
    specification = make_continuous_envelope_tie_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  x <- stats::runif(n, -1, 1)
  nuisance <- continuous_envelope_tie_nuisances(x, specification)
  z <- stats::rbinom(n, 1L, specification$probability_z1)
  probability_t <- ifelse(z == 1L, nuisance$pi1, nuisance$pi0)
  treatment <- stats::rbinom(n, 1L, probability_t)
  location <- continuous_envelope_tie_location(x, treatment, specification)
  tilt <- ifelse(
    z == 1L,
    continuous_envelope_tie_tilt(x, treatment, specification),
    0
  )
  residual <- observed_dgp_rtruncnorm(
    n,
    mean = tilt,
    lower = specification$residual_lower,
    upper = specification$residual_upper
  )
  y <- location + residual
  output <- data.frame(Y = y, T_obs = treatment, Z = z, X = x)
  output$e_true <- nuisance$e
  output$pi0_true <- nuisance$pi0
  output$pi1_true <- nuisance$pi1
  output$gamma0_oracle <- continuous_envelope_tie_gamma(
    y, x, 0L, specification
  )
  output$gamma1_oracle <- continuous_envelope_tie_gamma(
    y, x, 1L, specification
  )
  output$prob_Z <- rep(specification$probability_z1, n)
  output$prob_T <- probability_t
  attr(output, "continuous_envelope_tie_specification") <- specification
  output
}

continuous_envelope_tie_weighted_quantile <- function(
    y,
    probability_weight,
    level
) {
  ordering <- order(y)
  y <- y[ordering]
  probability_weight <- probability_weight[ordering]
  probability_weight <- probability_weight / sum(probability_weight)
  y[which(cumsum(probability_weight) >= level)[1L]]
}

continuous_envelope_tie_x_rule <- function(specification, nodes = 120L) {
  nodes <- as.integer(nodes)
  if (length(nodes) != 1L || is.na(nodes) || nodes < 30L) {
    stop("x_nodes must be an integer of at least 30.")
  }
  breaks <- c(-1, specification$tie_interval, 1)
  widths <- diff(breaks)
  counts <- pmax(10L, round(nodes * widths / 2))
  values <- weights <- region <- list()
  for (index in seq_along(widths)) {
    step <- widths[index] / counts[index]
    values[[index]] <- breaks[index] +
      (seq_len(counts[index]) - 0.5) * step
    weights[[index]] <- rep(step / 2, counts[index])
    region[[index]] <- rep(index, counts[index])
  }
  list(
    x = unlist(values),
    weight = unlist(weights),
    tie = unlist(region) == 2L
  )
}

compute_continuous_envelope_tie_population <- function(
    specification = make_continuous_envelope_tie_spec(),
    Lambda = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = 0,
    x_nodes = 120L,
    y_nodes = 4000L
) {
  if (!exists("compute_projected_envelopes", mode = "function")) {
    stop("Source R/sharp_bounds_plugin.R before computing population bounds.")
  }
  x_rule <- continuous_envelope_tie_x_rule(specification, x_nodes)
  y_nodes <- as.integer(y_nodes)
  if (length(y_nodes) != 1L || is.na(y_nodes) || y_nodes < 1000L) {
    stop("y_nodes must be an integer of at least 1000.")
  }
  residual_step <- (specification$residual_upper -
    specification$residual_lower) / y_nodes
  residual <- specification$residual_lower +
    (seq_len(y_nodes) - 0.5) * residual_step
  residual_weight <- residual_step * observed_dgp_dtruncnorm(
    residual,
    mean = 0,
    lower = specification$residual_lower,
    upper = specification$residual_upper
  )
  residual_weight <- residual_weight / sum(residual_weight)

  arm_bounds <- vector("list", 2L)
  diagnostics <- list()
  counter <- 1L
  for (treatment in 0:1) {
    lower_total <- upper_total <- 0
    for (x_index in seq_along(x_rule$x)) {
      x_value <- x_rule$x[x_index]
      y <- continuous_envelope_tie_location(
        x_value, treatment, specification
      ) + residual
      x <- rep(x_value, y_nodes)
      nuisance <- continuous_envelope_tie_nuisances(x, specification)
      gamma <- continuous_envelope_tie_gamma(
        y, x, treatment, specification
      )
      envelope <- compute_projected_envelopes(
        nuisance = nuisance,
        gamma = gamma,
        t_value = treatment,
        Lambda0 = Lambda,
        Lambda1 = Lambda,
        envelope_smoothing_epsilon = envelope_smoothing_epsilon,
        gamma_compatibility_action = "stop"
      )
      Lbar <- sum(residual_weight * envelope$L)
      Ubar <- sum(residual_weight * envelope$U)
      compatible <- !any(envelope$gamma_incompatible_original) &&
        all(envelope$L <= envelope$U + 1e-10) &&
        Lbar <= 1 + 1e-10 && Ubar >= 1 - 1e-10
      if (!compatible) {
        stop("The continuous population law is incompatible at this Lambda.")
      }
      tau <- (Ubar - 1) / (Ubar - Lbar)
      weighted_mass <- residual_weight * envelope$D
      q_lower <- continuous_envelope_tie_weighted_quantile(
        y, weighted_mass, 1 - tau
      )
      q_upper <- continuous_envelope_tie_weighted_quantile(
        y, weighted_mass, tau
      )
      lambda_lower <- envelope$L + envelope$D * (y <= q_lower)
      lambda_upper <- envelope$U - envelope$D * (y <= q_upper)
      r0 <- if (treatment == 1L) {
        nuisance$pi0[1L]
      } else {
        1 - nuisance$pi0[1L]
      }
      baseline_mean <- sum(residual_weight * y)
      conditional_lower <- r0 * baseline_mean +
        (1 - r0) * sum(residual_weight * lambda_lower * y)
      conditional_upper <- r0 * baseline_mean +
        (1 - r0) * sum(residual_weight * lambda_upper * y)
      lower_total <- lower_total + x_rule$weight[x_index] * conditional_lower
      upper_total <- upper_total + x_rule$weight[x_index] * conditional_upper
      diagnostics[[counter]] <- data.frame(
        treatment = treatment,
        X = x_value,
        X_weight = x_rule$weight[x_index],
        exact_tie_region = x_rule$tie[x_index],
        lower_branch_gap = sum(
          residual_weight * (envelope$candidate_lower - 1 / Lambda)
        ),
        upper_branch_gap = sum(
          residual_weight * (envelope$candidate_upper - Lambda)
        ),
        gamma_violation_probability = sum(
          residual_weight * envelope$gamma_incompatible_original
        ),
        envelope_crossing_probability = sum(
          residual_weight * (envelope$L > envelope$U + 1e-10)
        ),
        lower_normalization_residual =
          sum(residual_weight * lambda_lower) - 1,
        upper_normalization_residual =
          sum(residual_weight * lambda_upper) - 1,
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
    arm_bounds[[treatment + 1L]] <- c(
      lower = lower_total, upper = upper_total
    )
  }
  endpoints <- data.frame(
    design = specification$name,
    Lambda = Lambda,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    mu1_lower = arm_bounds[[2L]][["lower"]],
    mu1_upper = arm_bounds[[2L]][["upper"]],
    mu0_lower = arm_bounds[[1L]][["lower"]],
    mu0_upper = arm_bounds[[1L]][["upper"]],
    stringsAsFactors = FALSE
  )
  endpoints$ate_lower <- endpoints$mu1_lower - endpoints$mu0_upper
  endpoints$ate_upper <- endpoints$mu1_upper - endpoints$mu0_lower
  endpoints$ate_width <- endpoints$ate_upper - endpoints$ate_lower
  diagnostics <- do.call(rbind, diagnostics)
  rownames(diagnostics) <- NULL
  list(
    endpoints = endpoints,
    diagnostics = diagnostics,
    compatible = all(
      diagnostics$gamma_violation_probability <= 1e-12 &
        diagnostics$envelope_crossing_probability <= 1e-12
    ),
    population_tie_probability = specification$population_tie_probability,
    population_local_to_tie_probability =
      specification$population_local_to_tie_probability,
    marginal_first_stage = specification$marginal_first_stage,
    reference_ate = specification$outcome_treatment_shift
  )
}
