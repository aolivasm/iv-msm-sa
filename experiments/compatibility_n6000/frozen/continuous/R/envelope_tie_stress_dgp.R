# Observed-data stress law for envelope smoothing at moderate sensitivity.
#
# The design isolates the max/min envelope kink rather than threshold
# nonuniqueness.  X has three strata.  In the first two strata the instrument
# does not shift treatment and f(Y | T, X, Z = 1) equals
# f(Y | T, X, Z = 0), so gamma_t = 1 and both envelope branches tie exactly.
# Treatment is deliberately common in the first stratum and rare in the
# second.  This makes the unsmoothed envelope derivative large for t = 1 and
# t = 0, respectively.  The third stratum retains a regular, nonzero first
# stage, so the overall instrument remains relevant.
#
# At an exact branch tie, the conservative shifted log-mean-exp envelope has
# exactly the same value as the hard envelope but assigns derivative weight
# 1/2 instead of making a discontinuous 0/1 branch choice.  This is the
# intended mechanism of the simulation.

if (!exists("observed_dgp_rtruncnorm", mode = "function")) {
  if (!file.exists("R/observed_data_dgp.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/observed_data_dgp.R")
}

make_envelope_tie_stress_spec <- function(
    sensitivity_lambda = 3,
    x_probabilities = c(0.4, 0.4, 0.2),
    probability_z1 = 0.5,
    pi0 = c(0.9, 0.1, 0.1),
    pi1 = c(0.9, 0.1, 0.9),
    outcome_treatment_shift = 0.75,
    outcome_x_shift = c(-0.3, 0, 0.3),
    comparison_tilt_t0 = -0.075,
    comparison_tilt_t1 = 0.075,
    residual_lower = -2,
    residual_upper = 2
) {
  sensitivity_lambda <- as.numeric(sensitivity_lambda)
  x_probabilities <- as.numeric(x_probabilities)
  probability_z1 <- as.numeric(probability_z1)
  pi0 <- as.numeric(pi0)
  pi1 <- as.numeric(pi1)
  outcome_x_shift <- as.numeric(outcome_x_shift)
  comparison_tilt_t0 <- as.numeric(comparison_tilt_t0)
  comparison_tilt_t1 <- as.numeric(comparison_tilt_t1)

  if (length(sensitivity_lambda) != 1L ||
      !is.finite(sensitivity_lambda) || sensitivity_lambda <= 1) {
    stop("sensitivity_lambda must be one finite number greater than one.")
  }
  if (length(x_probabilities) != 3L ||
      any(!is.finite(x_probabilities)) ||
      any(x_probabilities <= 0) ||
      abs(sum(x_probabilities) - 1) > 1e-10) {
    stop("x_probabilities must contain three positive values summing to one.")
  }
  if (length(probability_z1) != 1L ||
      !is.finite(probability_z1) ||
      probability_z1 <= 0 || probability_z1 >= 1) {
    stop("probability_z1 must be strictly between zero and one.")
  }
  if (length(pi0) != 3L || length(pi1) != 3L ||
      any(!is.finite(c(pi0, pi1))) ||
      any(c(pi0, pi1) <= 0) || any(c(pi0, pi1) >= 1)) {
    stop("pi0 and pi1 must each contain three probabilities in (0, 1).")
  }
  if (!isTRUE(all.equal(pi0[1:2], pi1[1:2], tolerance = 1e-12))) {
    stop("The first two strata must be exact envelope-tie strata (pi0 = pi1).")
  }
  if (abs(sum(x_probabilities * (pi1 - pi0))) <= 1e-10) {
    stop("The marginal first stage must be nonzero.")
  }
  if (length(outcome_treatment_shift) != 1L ||
      !is.finite(outcome_treatment_shift) ||
      length(outcome_x_shift) != 3L ||
      any(!is.finite(outcome_x_shift)) ||
      length(comparison_tilt_t0) != 1L ||
      length(comparison_tilt_t1) != 1L ||
      !is.finite(comparison_tilt_t0) ||
      !is.finite(comparison_tilt_t1)) {
    stop("Outcome shifts must be finite and have the requested dimensions.")
  }
  if (length(residual_lower) != 1L || length(residual_upper) != 1L ||
      !is.finite(residual_lower) || !is.finite(residual_upper) ||
      residual_lower >= residual_upper) {
    stop("The residual support must be a nonempty finite interval.")
  }

  list(
    name = "envelope_tie_high_low_propensity",
    dimension = 2L,
    sensitivity_lambda = sensitivity_lambda,
    x_values = 0:2,
    x_probabilities = x_probabilities,
    probability_z1 = probability_z1,
    pi0 = pi0,
    pi1 = pi1,
    outcome_treatment_shift = outcome_treatment_shift,
    outcome_x_shift = outcome_x_shift,
    comparison_tilt_t0 = comparison_tilt_t0,
    comparison_tilt_t1 = comparison_tilt_t1,
    residual_lower = residual_lower,
    residual_upper = residual_upper,
    tie_strata = c(TRUE, TRUE, FALSE),
    marginal_first_stage = sum(x_probabilities * (pi1 - pi0)),
    population_tie_probability = sum(x_probabilities[1:2])
  )
}

envelope_tie_stress_x_index <- function(x, specification) {
  if (is.data.frame(x)) {
    if (!"X" %in% names(x)) stop("x data frame must contain X.")
    x <- x$X
  }
  x <- as.numeric(x)
  index <- match(x, specification$x_values)
  if (any(!is.finite(x)) || anyNA(index)) {
    stop("X must be coded numerically as 0, 1, or 2.")
  }
  index
}

envelope_tie_stress_nuisances <- function(x, specification) {
  index <- envelope_tie_stress_x_index(x, specification)
  list(
    e = rep(1 - specification$probability_z1, length(index)),
    pi0 = specification$pi0[index],
    pi1 = specification$pi1[index]
  )
}

envelope_tie_stress_location <- function(x, treatment, specification) {
  index <- envelope_tie_stress_x_index(x, specification)
  treatment <- rep_len(as.numeric(treatment), length(index))
  if (any(!is.finite(treatment)) || any(!treatment %in% 0:1)) {
    stop("treatment must contain only zeroes and ones.")
  }
  treatment * specification$outcome_treatment_shift +
    specification$outcome_x_shift[index]
}

envelope_tie_stress_gamma <- function(y, x, treatment, specification) {
  index <- envelope_tie_stress_x_index(x, specification)
  y <- rep_len(as.numeric(y), length(index))
  treatment <- rep_len(as.numeric(treatment), length(index))
  if (any(!is.finite(y)) || any(!treatment %in% 0:1)) {
    stop("y must be finite and treatment must contain only zeroes and ones.")
  }
  location <- envelope_tie_stress_location(x, treatment, specification)
  tilt <- ifelse(
    index == 3L,
    ifelse(
      treatment == 1,
      specification$comparison_tilt_t1,
      specification$comparison_tilt_t0
    ),
    0
  )
  residual <- y - location
  lower <- specification$residual_lower
  upper <- specification$residual_upper
  base_normalizer <- stats::pnorm(upper) - stats::pnorm(lower)
  tilted_normalizer <- stats::pnorm(upper, mean = tilt) -
    stats::pnorm(lower, mean = tilt)
  exp(tilt * residual - tilt^2 / 2) *
    base_normalizer / tilted_normalizer
}

generate_envelope_tie_stress_data <- function(
    n,
    specification = make_envelope_tie_stress_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  x <- sample(
    specification$x_values,
    size = n,
    replace = TRUE,
    prob = specification$x_probabilities
  )
  nuisance <- envelope_tie_stress_nuisances(x, specification)
  z <- stats::rbinom(n, 1L, specification$probability_z1)
  probability_t <- ifelse(z == 1L, nuisance$pi1, nuisance$pi0)
  treatment <- stats::rbinom(n, 1L, probability_t)
  location <- envelope_tie_stress_location(x, treatment, specification)
  x_index <- envelope_tie_stress_x_index(x, specification)
  comparison_tilt <- ifelse(
    treatment == 1L,
    specification$comparison_tilt_t1,
    specification$comparison_tilt_t0
  )
  tilt <- ifelse(z == 1L & x_index == 3L, comparison_tilt, 0)
  residual <- observed_dgp_rtruncnorm(
    n,
    mean = tilt,
    lower = specification$residual_lower,
    upper = specification$residual_upper
  )
  y <- location + residual

  output <- data.frame(
    Y = y,
    T_obs = treatment,
    Z = z,
    X = x,
    X1 = as.numeric(x == 0),
    X2 = as.numeric(x == 1)
  )
  output$e_true <- nuisance$e
  output$pi0_true <- nuisance$pi0
  output$pi1_true <- nuisance$pi1
  output$gamma0_oracle <- envelope_tie_stress_gamma(
    y, x, 0L, specification
  )
  output$gamma1_oracle <- envelope_tie_stress_gamma(
    y, x, 1L, specification
  )
  output$prob_Z <- rep(specification$probability_z1, n)
  output$prob_T <- probability_t
  attr(output, "envelope_tie_stress_specification") <- specification
  output
}

envelope_tie_stress_baseline_rule <- function(
    x,
    treatment,
    specification,
    nodes = 20000L
) {
  nodes <- as.integer(nodes)
  if (length(nodes) != 1L || is.na(nodes) || nodes < 1000L) {
    stop("nodes must be one integer of at least 1000.")
  }
  width <- (specification$residual_upper -
    specification$residual_lower) / nodes
  residual <- specification$residual_lower +
    (seq_len(nodes) - 0.5) * width
  weight <- width * observed_dgp_dtruncnorm(
    residual,
    mean = 0,
    lower = specification$residual_lower,
    upper = specification$residual_upper
  )
  weight <- weight / sum(weight)
  location <- envelope_tie_stress_location(x, treatment, specification)
  list(y = location + residual, residual = residual, weight = weight)
}

envelope_tie_stress_weighted_quantile <- function(
    y,
    probability_weight,
    level
) {
  if (!is.finite(level) || level < 0 || level > 1 ||
      length(y) != length(probability_weight) ||
      any(!is.finite(y)) || any(!is.finite(probability_weight)) ||
      any(probability_weight < 0) || sum(probability_weight) <= 0) {
    stop("Invalid weighted-quantile inputs.")
  }
  ordering <- order(y)
  y <- y[ordering]
  probability_weight <- probability_weight[ordering]
  probability_weight <- probability_weight / sum(probability_weight)
  y[which(cumsum(probability_weight) >= level)[1L]]
}

compute_envelope_tie_stress_population <- function(
    specification = make_envelope_tie_stress_spec(),
    Lambda = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = 0,
    nodes = 20000L
) {
  if (!exists("compute_projected_envelopes", mode = "function")) {
    stop("Source R/sharp_bounds_plugin.R before computing population bounds.")
  }
  Lambda <- as.numeric(Lambda)
  if (length(Lambda) != 1L || !is.finite(Lambda) || Lambda <= 1) {
    stop("Lambda must be one finite number greater than one.")
  }

  arm_bounds <- vector("list", 2L)
  diagnostic_rows <- list()
  counter <- 1L

  for (treatment in 0:1) {
    lower_total <- 0
    upper_total <- 0
    for (x_index in seq_along(specification$x_values)) {
      x_value <- specification$x_values[x_index]
      rule <- envelope_tie_stress_baseline_rule(
        x_value, treatment, specification, nodes
      )
      x <- rep(x_value, length(rule$y))
      nuisance <- envelope_tie_stress_nuisances(x, specification)
      gamma <- envelope_tie_stress_gamma(
        rule$y, x, treatment, specification
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

      Lbar <- sum(rule$weight * envelope$L)
      Ubar <- sum(rule$weight * envelope$U)
      compatible <- !any(envelope$gamma_incompatible_original) &&
        all(envelope$L <= envelope$U + 1e-10) &&
        Lbar <= 1 + 1e-10 && Ubar >= 1 - 1e-10
      if (!compatible) {
        stop("The population law is incompatible at the requested Lambda.")
      }
      tau <- (Ubar - 1) / (Ubar - Lbar)
      weighted_mass <- rule$weight * envelope$D
      q_lower <- envelope_tie_stress_weighted_quantile(
        rule$y, weighted_mass, 1 - tau
      )
      q_upper <- envelope_tie_stress_weighted_quantile(
        rule$y, weighted_mass, tau
      )
      lambda_lower <- envelope$L + envelope$D * (rule$y <= q_lower)
      lambda_upper <- envelope$U - envelope$D * (rule$y <= q_upper)

      r0 <- if (treatment == 1L) {
        nuisance$pi0[1L]
      } else {
        1 - nuisance$pi0[1L]
      }
      baseline_mean <- sum(rule$weight * rule$y)
      conditional_lower <- r0 * baseline_mean +
        (1 - r0) * sum(rule$weight * lambda_lower * rule$y)
      conditional_upper <- r0 * baseline_mean +
        (1 - r0) * sum(rule$weight * lambda_upper * rule$y)
      lower_total <- lower_total +
        specification$x_probabilities[x_index] * conditional_lower
      upper_total <- upper_total +
        specification$x_probabilities[x_index] * conditional_upper

      lower_branch_gap <- envelope$candidate_lower - 1 / Lambda
      upper_branch_gap <- envelope$candidate_upper - Lambda
      diagnostic_rows[[counter]] <- data.frame(
        design = specification$name,
        treatment = treatment,
        X = x_value,
        X_probability = specification$x_probabilities[x_index],
        Lambda = Lambda,
        envelope_smoothing_epsilon = envelope_smoothing_epsilon,
        exact_tie_stratum = specification$tie_strata[x_index],
        compatible = compatible,
        gamma_violation_probability = sum(
          rule$weight * envelope$gamma_incompatible_original
        ),
        envelope_crossing_probability = sum(
          rule$weight * (envelope$L > envelope$U + 1e-10)
        ),
        lower_branch_gap = sum(rule$weight * lower_branch_gap),
        upper_branch_gap = sum(rule$weight * upper_branch_gap),
        lower_smoothing_gap = sum(
          rule$weight * envelope$lower_smoothing_gap
        ),
        upper_smoothing_gap = sum(
          rule$weight * envelope$upper_smoothing_gap
        ),
        Lbar = Lbar,
        Ubar = Ubar,
        tau = tau,
        q_lower = q_lower,
        q_upper = q_upper,
        lower_normalization_residual =
          sum(rule$weight * lambda_lower) - 1,
        upper_normalization_residual =
          sum(rule$weight * lambda_upper) - 1,
        conditional_lower = conditional_lower,
        conditional_upper = conditional_upper,
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
    arm_bounds[[treatment + 1L]] <- c(
      lower = lower_total,
      upper = upper_total
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

  diagnostics <- do.call(rbind, diagnostic_rows)
  rownames(diagnostics) <- NULL
  list(
    endpoints = endpoints,
    diagnostics = diagnostics,
    compatible = all(diagnostics$compatible),
    population_tie_probability = specification$population_tie_probability,
    marginal_first_stage = specification$marginal_first_stage,
    reference_ate = specification$outcome_treatment_shift
  )
}

compute_envelope_tie_stress_derivatives <- function(
    specification = make_envelope_tie_stress_spec(),
    Lambda = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = 0.01
) {
  if (!exists("compute_eif_envelope_derivatives", mode = "function")) {
    stop("Source R/sharp_bounds_plugin.R before computing derivatives.")
  }
  rows <- list()
  counter <- 1L
  for (treatment in 0:1) {
    for (x_index in seq_along(specification$x_values)) {
      x_value <- specification$x_values[x_index]
      nuisance <- envelope_tie_stress_nuisances(x_value, specification)
      hard <- compute_eif_envelope_derivatives(
        nuisance, 1, treatment, Lambda, Lambda,
        envelope_smoothing_epsilon = 0
      )
      smooth <- compute_eif_envelope_derivatives(
        nuisance, 1, treatment, Lambda, Lambda,
        envelope_smoothing_epsilon = envelope_smoothing_epsilon
      )
      rows[[counter]] <- data.frame(
        treatment = treatment,
        X = x_value,
        exact_tie_stratum = specification$tie_strata[x_index],
        r0 = hard$r0,
        r1 = hard$r1,
        lower_candidate = hard$candidate_lower,
        upper_candidate = hard$candidate_upper,
        hard_lower_weight = hard$lower_derivative_weight,
        smooth_lower_weight = smooth$lower_derivative_weight,
        hard_upper_weight = hard$upper_derivative_weight,
        smooth_upper_weight = smooth$upper_derivative_weight,
        hard_lower_gamma_derivative = hard$lower$gamma,
        smooth_lower_gamma_derivative = smooth$lower$gamma,
        hard_upper_gamma_derivative = hard$upper$gamma,
        smooth_upper_gamma_derivative = smooth$upper$gamma,
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}
