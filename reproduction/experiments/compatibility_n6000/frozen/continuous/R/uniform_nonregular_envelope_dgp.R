# Fixed and local observed-data laws for nonuniform exact-envelope inference.
#
# The center law is deliberately more demanding than a finite-stratum tie
# design.  Let X ~ Uniform(-1, 1), I0 = [-rho, rho], and
#
#   pi_1(x) - pi_0(x) = d_max g(x),
#
# where g(x) is zero on I0 and, outside I0, is infinitely flat at the
# boundary of I0:
#
#   g(x) = exp{a - a / u(x)^2},
#
# with u(x) the distance from I0 rescaled to [0, 1].  Outcomes are independent
# of Z conditional on (T, X), so gamma_t(y, x) = 1.  Consequently both
# max/min branches defining the exact dual-MSM envelopes tie throughout I0.
# The center law therefore has P(Delta = 0) = P(X in I0) > 0.
#
# The positive gaps immediately outside I0 also violate every polynomial
# margin condition.  In fact, for 0 < r < d_max,
#
#   P{0 < |pi_1(X)-pi_0(X)| <= r}
#     = P(X outside I0) / sqrt{1 + log(d_max/r)/a},
#
# which decreases more slowly than r^alpha for every alpha > 0.
#
# For a uniform local experiment, the exact-tie region is perturbed by
#
#   h s(X) / sqrt(n),
#
# where s(X) is one smooth sine wave across I0.  The perturbation is zero at
# the two boundaries, integrates to zero, preserves the marginal first stage,
# and approaches the fixed center law at the root-n rate.  A symmetric grid
# of h values therefore probes nonuniform coverage in shrinking neighborhoods
# of the exceptional law without changing the estimator.

if (!exists("observed_dgp_rtruncnorm", mode = "function")) {
  if (!file.exists("R/observed_data_dgp.R")) {
    stop("Run this file from the repository root.", call. = FALSE)
  }
  source("R/observed_data_dgp.R")
}

make_uniform_nonregular_envelope_spec <- function(
    sensitivity_lambda = 4,
    tie_interval = c(-0.8, 0.8),
    probability_z1 = 0.5,
    pi0_intercept = 0.3,
    pi0_slope = 0,
    maximum_first_stage = 0.6,
    outside_flatness = 0.1,
    local_coefficient = 0,
    nominal_sample_size = NULL,
    outcome_treatment_shift = 0.75,
    outcome_x_slope = 0.4,
    residual_lower = -2,
    residual_upper = 2
) {
  numeric_values <- c(
    sensitivity_lambda, tie_interval, probability_z1,
    pi0_intercept, pi0_slope, maximum_first_stage, outside_flatness,
    local_coefficient, outcome_treatment_shift, outcome_x_slope,
    residual_lower, residual_upper
  )
  if (any(!is.finite(numeric_values)) || sensitivity_lambda <= 1 ||
      length(tie_interval) != 2L || tie_interval[1L] <= -1 ||
      tie_interval[2L] >= 1 || tie_interval[1L] >= tie_interval[2L] ||
      probability_z1 <= 0 || probability_z1 >= 1 ||
      maximum_first_stage <= 0 || outside_flatness <= 0 ||
      residual_lower >= residual_upper) {
    stop("Invalid uniform nonregular envelope specification.")
  }

  local_coefficient <- as.numeric(local_coefficient)
  if (is.null(nominal_sample_size)) {
    if (abs(local_coefficient) > 0) {
      stop("nominal_sample_size is required when local_coefficient is nonzero.")
    }
    nominal_sample_size <- NA_integer_
    local_gap <- 0
  } else {
    nominal_sample_size <- as.integer(nominal_sample_size)
    if (length(nominal_sample_size) != 1L || is.na(nominal_sample_size) ||
        nominal_sample_size < 1L) {
      stop("nominal_sample_size must be one positive integer.")
    }
    local_gap <- local_coefficient / sqrt(nominal_sample_size)
  }

  tie_mass <- diff(tie_interval) / 2
  outside_mass <- 1 - tie_mass
  outside_shape_mean <- stats::integrate(
    function(u) exp(outside_flatness * (1 - 1 / u^2)),
    lower = 0, upper = 1, rel.tol = 1e-11
  )$value
  marginal_first_stage <-
    maximum_first_stage * outside_mass * outside_shape_mean

  specification <- list(
    name = if (local_gap == 0) {
      "fixed_uniform_nonregular_envelope"
    } else {
      "local_uniform_nonregular_envelope"
    },
    dimension = 1L,
    sensitivity_lambda = sensitivity_lambda,
    tie_interval = tie_interval,
    tie_midpoint = mean(tie_interval),
    tie_mass = tie_mass,
    outside_mass = outside_mass,
    probability_z1 = probability_z1,
    pi0_intercept = pi0_intercept,
    pi0_slope = pi0_slope,
    maximum_first_stage = maximum_first_stage,
    outside_flatness = outside_flatness,
    outside_shape_mean = outside_shape_mean,
    local_coefficient = local_coefficient,
    nominal_sample_size = nominal_sample_size,
    local_gap = local_gap,
    outcome_treatment_shift = outcome_treatment_shift,
    outcome_x_slope = outcome_x_slope,
    residual_lower = residual_lower,
    residual_upper = residual_upper,
    gamma_t = 1,
    population_exact_tie_probability = if (local_gap == 0) tie_mass else 0,
    population_limit_tie_probability = tie_mass,
    marginal_first_stage = marginal_first_stage,
    local_perturbation_integral = 0,
    purpose = paste0(
      "fixed positive-mass envelope tie with an infinitely flat near-tie ",
      "continuum and a symmetric root-n local perturbation submodel"
    )
  )
  class(specification) <- c(
    "uniform_nonregular_envelope_specification", "list"
  )

  check_x <- seq(-1, 1, length.out = 20001L)
  nuisance <- uniform_nonregular_envelope_nuisances(check_x, specification)
  overlap_constant <- min(
    probability_z1, 1 - probability_z1,
    nuisance$pi0, 1 - nuisance$pi0,
    nuisance$pi1, 1 - nuisance$pi1
  )
  if (!is.finite(overlap_constant) || overlap_constant <= 0.05) {
    stop(
      "The requested treatment probabilities violate the 0.05 overlap margin."
    )
  }
  specification$overlap_constant <- overlap_constant
  specification
}

uniform_nonregular_envelope_validate_x <- function(x) {
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

uniform_nonregular_envelope_limit_tie_region <- function(x, specification) {
  x <- uniform_nonregular_envelope_validate_x(x)
  x >= specification$tie_interval[1L] &
    x <= specification$tie_interval[2L]
}

uniform_nonregular_envelope_local_direction <- function(x, specification) {
  x <- uniform_nonregular_envelope_validate_x(x)
  tie <- uniform_nonregular_envelope_limit_tie_region(x, specification)
  output <- numeric(length(x))
  output[tie] <- sin(
    2 * pi *
      (x[tie] - specification$tie_interval[1L]) /
      diff(specification$tie_interval)
  )
  output
}

uniform_nonregular_envelope_outside_coordinate <- function(x, specification) {
  x <- uniform_nonregular_envelope_validate_x(x)
  lower <- specification$tie_interval[1L]
  upper <- specification$tie_interval[2L]
  ifelse(
    x < lower,
    (lower - x) / (lower + 1),
    ifelse(x > upper, (x - upper) / (1 - upper), 0)
  )
}

uniform_nonregular_envelope_outside_shape <- function(x, specification) {
  coordinate <- uniform_nonregular_envelope_outside_coordinate(
    x, specification
  )
  output <- numeric(length(coordinate))
  positive <- coordinate > 0
  output[positive] <- exp(
    specification$outside_flatness *
      (1 - 1 / coordinate[positive]^2)
  )
  output
}

uniform_nonregular_envelope_propensity_gap <- function(x, specification) {
  x <- uniform_nonregular_envelope_validate_x(x)
  specification$maximum_first_stage *
    uniform_nonregular_envelope_outside_shape(x, specification) +
    specification$local_gap *
      uniform_nonregular_envelope_local_direction(x, specification)
}

uniform_nonregular_envelope_nuisances <- function(x, specification) {
  x <- uniform_nonregular_envelope_validate_x(x)
  pi0 <- specification$pi0_intercept + specification$pi0_slope * x
  pi1 <- pi0 + uniform_nonregular_envelope_propensity_gap(x, specification)
  list(
    e = rep(1 - specification$probability_z1, length(x)),
    pi0 = pi0,
    pi1 = pi1
  )
}

uniform_nonregular_envelope_location <- function(x, treatment, specification) {
  x <- uniform_nonregular_envelope_validate_x(x)
  treatment <- rep_len(as.numeric(treatment), length(x))
  if (any(!treatment %in% 0:1)) {
    stop("treatment must contain only zeroes and ones.")
  }
  specification$outcome_treatment_shift * treatment +
    specification$outcome_x_slope * x
}

uniform_nonregular_envelope_gamma <- function(y, x, treatment, specification) {
  x <- uniform_nonregular_envelope_validate_x(x)
  rep(1, max(length(y), length(x), length(treatment)))
}

generate_uniform_nonregular_envelope_data <- function(
    n,
    specification = make_uniform_nonregular_envelope_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  x <- stats::runif(n, -1, 1)
  nuisance <- uniform_nonregular_envelope_nuisances(x, specification)
  z <- stats::rbinom(n, 1L, specification$probability_z1)
  treatment_probability <- ifelse(z == 1L, nuisance$pi1, nuisance$pi0)
  treatment <- stats::rbinom(n, 1L, treatment_probability)
  location <- uniform_nonregular_envelope_location(
    x, treatment, specification
  )
  residual <- observed_dgp_rtruncnorm(
    n, mean = 0,
    lower = specification$residual_lower,
    upper = specification$residual_upper
  )
  output <- data.frame(
    Y = location + residual,
    T_obs = treatment,
    Z = z,
    X = x,
    X_limit_tie = as.numeric(
      uniform_nonregular_envelope_limit_tie_region(x, specification)
    ),
    stringsAsFactors = FALSE
  )
  output$e_true <- nuisance$e
  output$pi0_true <- nuisance$pi0
  output$pi1_true <- nuisance$pi1
  output$gamma0_oracle <- rep(1, n)
  output$gamma1_oracle <- rep(1, n)
  output$prob_Z <- rep(specification$probability_z1, n)
  output$prob_T <- treatment_probability
  attr(output, "uniform_nonregular_envelope_specification") <- specification
  output
}

uniform_nonregular_envelope_positive_gap_cdf <- function(
    threshold,
    specification = make_uniform_nonregular_envelope_spec()
) {
  threshold <- as.numeric(threshold)
  if (any(!is.finite(threshold)) || any(threshold < 0)) {
    stop("threshold must contain finite nonnegative values.")
  }
  maximum <- specification$maximum_first_stage
  flatness <- specification$outside_flatness
  outside_probability <- ifelse(
    threshold <= 0,
    0,
    ifelse(
      threshold >= maximum,
      specification$outside_mass,
      specification$outside_mass /
        sqrt(1 + log(maximum / threshold) / flatness)
    )
  )
  local_amplitude <- abs(specification$local_gap)
  local_probability <- if (local_amplitude == 0) {
    rep(0, length(threshold))
  } else {
    specification$tie_mass * ifelse(
      threshold <= 0,
      0,
      ifelse(
        threshold >= local_amplitude,
        1,
        2 / pi * asin(pmin(1, threshold / local_amplitude))
      )
    )
  }
  outside_probability + local_probability
}

uniform_nonregular_envelope_margin_audit <- function(
    specification = make_uniform_nonregular_envelope_spec(),
    thresholds = 10^seq(-1, -12, by = -1),
    exponents = c(0.25, 0.5, 1)
) {
  thresholds <- sort(unique(as.numeric(thresholds)), decreasing = TRUE)
  exponents <- sort(unique(as.numeric(exponents)))
  if (!length(thresholds) || any(!is.finite(thresholds)) ||
      any(thresholds <= 0) || !length(exponents) ||
      any(!is.finite(exponents)) || any(exponents <= 0)) {
    stop("thresholds and exponents must be positive finite values.")
  }
  positive_probability <- uniform_nonregular_envelope_positive_gap_cdf(
    thresholds, specification
  )
  output <- expand.grid(
    threshold = thresholds,
    exponent = exponents,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  output$exact_tie_probability <-
    specification$population_exact_tie_probability
  output$positive_small_gap_probability <- rep(
    positive_probability, times = length(exponents)
  )
  output$total_small_gap_probability <-
    output$exact_tie_probability + output$positive_small_gap_probability
  output$polynomial_margin_ratio <-
    output$positive_small_gap_probability / output$threshold^output$exponent
  output
}

uniform_nonregular_envelope_x_rule <- function(
    specification,
    nodes = 240L
) {
  nodes <- as.integer(nodes)
  if (length(nodes) != 1L || is.na(nodes) || nodes < 60L) {
    stop("x_nodes must be one integer of at least 60.")
  }
  breaks <- c(-1, specification$tie_interval, 1)
  widths <- diff(breaks)
  counts <- pmax(20L, round(nodes * widths / 2))
  values <- weights <- regions <- list()
  for (index in seq_along(widths)) {
    step <- widths[index] / counts[index]
    values[[index]] <- breaks[index] +
      (seq_len(counts[index]) - 0.5) * step
    weights[[index]] <- rep(step / 2, counts[index])
    regions[[index]] <- rep(index, counts[index])
  }
  list(
    x = unlist(values),
    weight = unlist(weights),
    limit_tie_region = unlist(regions) == 2L
  )
}

uniform_nonregular_envelope_weighted_quantile <- function(
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

compute_uniform_nonregular_envelope_population <- function(
    specification = make_uniform_nonregular_envelope_spec(),
    Lambda = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = 0,
    x_nodes = 240L,
    y_nodes = 4000L
) {
  if (!exists("compute_projected_envelopes", mode = "function")) {
    stop("Source R/sharp_bounds_plugin.R before computing population bounds.")
  }
  x_rule <- uniform_nonregular_envelope_x_rule(specification, x_nodes)
  y_nodes <- as.integer(y_nodes)
  if (length(y_nodes) != 1L || is.na(y_nodes) || y_nodes < 1000L) {
    stop("y_nodes must be one integer of at least 1000.")
  }
  residual_step <- (specification$residual_upper -
    specification$residual_lower) / y_nodes
  residual <- specification$residual_lower +
    (seq_len(y_nodes) - 0.5) * residual_step
  residual_weight <- residual_step * observed_dgp_dtruncnorm(
    residual, mean = 0,
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
      y <- uniform_nonregular_envelope_location(
        x_value, treatment, specification
      ) + residual
      nuisance <- uniform_nonregular_envelope_nuisances(
        rep(x_value, y_nodes), specification
      )
      envelope <- compute_projected_envelopes(
        nuisance = nuisance,
        gamma = rep(1, y_nodes),
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
        stop("The uniform nonregular law is incompatible at this Lambda.")
      }
      tau <- (Ubar - 1) / (Ubar - Lbar)
      weighted_mass <- residual_weight * envelope$D
      q_lower <- uniform_nonregular_envelope_weighted_quantile(
        y, weighted_mass, 1 - tau
      )
      q_upper <- uniform_nonregular_envelope_weighted_quantile(
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
      lower_total <- lower_total +
        x_rule$weight[x_index] * conditional_lower
      upper_total <- upper_total +
        x_rule$weight[x_index] * conditional_upper
      diagnostics[[counter]] <- data.frame(
        treatment = treatment,
        X = x_value,
        X_weight = x_rule$weight[x_index],
        limit_tie_region = x_rule$limit_tie_region[x_index],
        propensity_gap = nuisance$pi1[1L] - nuisance$pi0[1L],
        lower_branch_gap = sum(
          residual_weight *
            (envelope$candidate_lower - 1 / Lambda)
        ),
        upper_branch_gap = sum(
          residual_weight *
            (envelope$candidate_upper - Lambda)
        ),
        exact_lower_branch_tie = max(abs(
          envelope$candidate_lower - 1 / Lambda
        )) < 1e-12,
        exact_upper_branch_tie = max(abs(
          envelope$candidate_upper - Lambda
        )) < 1e-12,
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
    nominal_sample_size = specification$nominal_sample_size,
    local_coefficient = specification$local_coefficient,
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
    compatible = TRUE,
    population_exact_tie_probability =
      specification$population_exact_tie_probability,
    population_limit_tie_probability =
      specification$population_limit_tie_probability,
    marginal_first_stage = specification$marginal_first_stage,
    reference_ate = specification$outcome_treatment_shift
  )
}
