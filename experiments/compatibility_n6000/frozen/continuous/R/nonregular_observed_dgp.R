# Observed-data laws for targeted nonregularity experiments.
#
# The laws are specified directly through
#
#   P(X) P(Z | X) P(T | X, Z) P(Y | T, X, Z).
#
# X is binary with equal stratum probabilities.  In stratum X = 0 the
# instrument does not shift treatment and the comparison outcome law equals
# the baseline outcome law.  In stratum X = 1 the instrument has a nonzero
# first stage and the comparison outcome density is an exponential tilt of
# the baseline density.  Hence the complete observed law is incompatible with
# Lambda = 1 but compatible with the study value Lambda = 3.
#
# The `branch_tie` law has ordinary truncated-normal outcomes in X = 0.  Both
# branches defining L and U tie there with positive probability.  The
# `nonunique_threshold` law replaces that stratum's residual distribution by
# a two-interval mixture.  Its lower-component probability is 1/(Lambda+1),
# making the lower check-loss optimizer equal to the entire gap between the
# two intervals.  The population endpoint itself remains uniquely defined.
#
# The `atomic_nonunique_threshold` law is a stronger stress design. In X = 0,
# its residual distribution places probabilities
#
#   {1/(Lambda+1), (Lambda-1)/(Lambda+1), 1/(Lambda+1)}
#
# on three separated atoms. Both relevant quantile levels then have
# set-valued optimizer intervals, and alternative endpoints classify an atom
# of positive probability differently. This remains useful as a diagnostic of
# the check-loss pseudo-functional. It is not the official direct-functional
# nonregularity design: an endpoint that reallocates an atom can have a nonzero
# extremizer-normalization residual. The official continuous-gap construction,
# whose alternative thresholds preserve the same normalized direct extremizer,
# is implemented in R/multibranch_nonunique_threshold_pilot.R.

if (!exists("observed_dgp_gauss_legendre", mode = "function")) {
  if (!file.exists("R/observed_data_dgp.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/observed_data_dgp.R")
}

make_nonregular_observed_dgp_spec <- function(
    design = c(
      "branch_tie", "nonunique_threshold", "atomic_nonunique_threshold"
    ),
    sensitivity_lambda = 3
) {
  design <- match.arg(design)
  sensitivity_lambda <- as.numeric(sensitivity_lambda)
  if (length(sensitivity_lambda) != 1L ||
      !is.finite(sensitivity_lambda) || sensitivity_lambda <= 1) {
    stop("sensitivity_lambda must be one finite number greater than one.")
  }
  atomic_design <- identical(design, "atomic_nonunique_threshold")
  lower_mixture_probability <- 1 / (sensitivity_lambda + 1)
  list(
    name = design,
    dimension = 1L,
    sensitivity_lambda = sensitivity_lambda,
# Put 80% probability on the atomic pseudo-functional diagnostic stratum while
    # retaining a nonzero observed first stage through X = 1.
    probability_x1 = if (atomic_design) 0.2 else 0.5,
    probability_z1 = 0.5,
    pi0_x0 = 0.5,
    pi1_x0 = 0.5,
    pi0_x1 = 0.3,
    pi1_x1 = 0.7,
    outcome_treatment_shift = 0.75,
    outcome_x_shift = 0.35,
    residual_lower = -2,
    residual_upper = 2,
    comparison_tilt_t0 = -0.075,
    comparison_tilt_t1 = 0.075,
    lower_mixture_probability = lower_mixture_probability,
    low_interval = c(-2, -1),
    high_interval = c(1, 2),
    atomic_residual_values = c(-4, 0, 4),
    atomic_residual_probabilities = c(
      lower_mixture_probability,
      1 - 2 * lower_mixture_probability,
      lower_mixture_probability
    )
  )
}

nonregular_observed_is_nonunique <- function(specification) {
  specification$name %in% c(
    "nonunique_threshold", "atomic_nonunique_threshold"
  )
}

nonregular_observed_is_atomic <- function(specification) {
  identical(specification$name, "atomic_nonunique_threshold")
}

nonregular_observed_validate_x <- function(x) {
  if (is.data.frame(x)) {
    if (!"X" %in% names(x)) stop("x data frame must contain X.")
    x <- x$X
  }
  x <- as.numeric(x)
  if (any(!is.finite(x)) || any(!x %in% c(0, 1))) {
    stop("X must be coded numerically as zero or one.")
  }
  x
}

nonregular_observed_nuisances <- function(x, specification) {
  x <- nonregular_observed_validate_x(x)
  list(
    e = rep(1 - specification$probability_z1, length(x)),
    pi0 = ifelse(x == 0, specification$pi0_x0, specification$pi0_x1),
    pi1 = ifelse(x == 0, specification$pi1_x0, specification$pi1_x1)
  )
}

nonregular_observed_location <- function(x, treatment, specification) {
  x <- nonregular_observed_validate_x(x)
  treatment <- as.integer(treatment)
  if (length(treatment) != 1L || !treatment %in% 0:1) {
    stop("treatment must be zero or one.")
  }
  (treatment - 0.5) * specification$outcome_treatment_shift +
    specification$outcome_x_shift * x
}

nonregular_observed_tilt <- function(x, treatment, specification) {
  x <- nonregular_observed_validate_x(x)
  treatment <- as.integer(treatment)
  if (length(treatment) != 1L || !treatment %in% 0:1) {
    stop("treatment must be zero or one.")
  }
  value <- if (treatment == 1L) {
    specification$comparison_tilt_t1
  } else {
    specification$comparison_tilt_t0
  }
  ifelse(x == 0, 0, value)
}

nonregular_observed_gamma <- function(y, x, treatment, specification) {
  x <- nonregular_observed_validate_x(x)
  y <- rep_len(as.numeric(y), length(x))
  location <- nonregular_observed_location(x, treatment, specification)
  tilt <- nonregular_observed_tilt(x, treatment, specification)
  residual <- y - location
  lower <- specification$residual_lower
  upper <- specification$residual_upper
  base_normalizer <- stats::pnorm(upper) - stats::pnorm(lower)
  tilted_normalizer <- stats::pnorm(upper, mean = tilt) -
    stats::pnorm(lower, mean = tilt)
  value <- exp(tilt * residual - tilt^2 / 2) *
    base_normalizer / tilted_normalizer
  value[x == 0] <- 1
  value
}

nonregular_observed_sample_gap_residual <- function(n, specification) {
  low <- stats::rbinom(
    n, 1L, specification$lower_mixture_probability
  ) == 1L
  output <- numeric(n)
  output[low] <- stats::runif(
    sum(low), specification$low_interval[1L],
    specification$low_interval[2L]
  )
  output[!low] <- stats::runif(
    sum(!low), specification$high_interval[1L],
    specification$high_interval[2L]
  )
  output
}

nonregular_observed_sample_atomic_residual <- function(n, specification) {
  specification$atomic_residual_values[base::sample.int(
    length(specification$atomic_residual_values),
    size = n,
    replace = TRUE,
    prob = specification$atomic_residual_probabilities
  )]
}

generate_nonregular_observed_iv_data <- function(
    n,
    specification = make_nonregular_observed_dgp_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  x <- stats::rbinom(n, 1L, specification$probability_x1)
  nuisance <- nonregular_observed_nuisances(x, specification)
  z <- stats::rbinom(n, 1L, specification$probability_z1)
  probability_t <- ifelse(z == 1L, nuisance$pi1, nuisance$pi0)
  treatment <- stats::rbinom(n, 1L, probability_t)
  location <- numeric(n)
  tilt <- numeric(n)
  for (arm in 0:1) {
    index <- treatment == arm
    if (any(index)) {
      location[index] <- nonregular_observed_location(
        x[index], arm, specification
      )
      tilt[index] <- ifelse(
        z[index] == 1L,
        nonregular_observed_tilt(x[index], arm, specification),
        0
      )
    }
  }
  residual <- numeric(n)
  gap <- x == 0 & identical(specification$name, "nonunique_threshold")
  if (any(gap)) {
    residual[gap] <- nonregular_observed_sample_gap_residual(
      sum(gap), specification
    )
  }
  atomic <- x == 0 & nonregular_observed_is_atomic(specification)
  if (any(atomic)) {
    residual[atomic] <- nonregular_observed_sample_atomic_residual(
      sum(atomic), specification
    )
  }
  regular <- !(gap | atomic)
  if (any(regular)) {
    residual[regular] <- observed_dgp_rtruncnorm(
      sum(regular), mean = tilt[regular],
      lower = specification$residual_lower,
      upper = specification$residual_upper
    )
  }
  y <- location + residual
  output <- data.frame(Y = y, T_obs = treatment, Z = z, X = x)
  output$e_true <- nuisance$e
  output$pi0_true <- nuisance$pi0
  output$pi1_true <- nuisance$pi1
  output$gamma0_oracle <- nonregular_observed_gamma(
    y, x, 0L, specification
  )
  output$gamma1_oracle <- nonregular_observed_gamma(
    y, x, 1L, specification
  )
  output$prob_Z <- rep(specification$probability_z1, n)
  output$prob_T <- probability_t
  attr(output, "nonregular_observed_dgp_specification") <- specification
  output
}

nonregular_observed_baseline_rule <- function(
    x,
    treatment,
    specification,
    nodes = 10000L
) {
  nodes <- as.integer(nodes)
  if (length(nodes) != 1L || is.na(nodes) || nodes < 100L) {
    stop("nodes must be one integer of at least 100.")
  }
  midpoint_rule <- function(lower, upper, count) {
    width <- (upper - lower) / count
    list(
      nodes = lower + (seq_len(count) - 0.5) * width,
      weights = rep(width, count)
    )
  }
  location <- nonregular_observed_location(x, treatment, specification)
  if (x == 0 && identical(specification$name, "nonunique_threshold")) {
    low <- midpoint_rule(
      specification$low_interval[1L],
      specification$low_interval[2L], nodes
    )
    high <- midpoint_rule(
      specification$high_interval[1L],
      specification$high_interval[2L], nodes
    )
    residual <- c(low$nodes, high$nodes)
    weight <- c(
      specification$lower_mixture_probability * low$weights,
      (1 - specification$lower_mixture_probability) * high$weights
    )
  } else if (x == 0 && nonregular_observed_is_atomic(specification)) {
    residual <- specification$atomic_residual_values
    weight <- specification$atomic_residual_probabilities
  } else {
    rule <- midpoint_rule(
      specification$residual_lower, specification$residual_upper, nodes
    )
    residual <- rule$nodes
    weight <- rule$weights * observed_dgp_dtruncnorm(
      residual, mean = 0,
      lower = specification$residual_lower,
      upper = specification$residual_upper
    )
  }
  weight <- weight / sum(weight)
  list(y = location + residual, residual = residual, weight = weight)
}

nonregular_weighted_quantile <- function(y, probability_weight, level) {
  ordering <- order(y)
  y <- y[ordering]
  probability_weight <- probability_weight[ordering]
  probability_weight <- probability_weight / sum(probability_weight)
  y[which(cumsum(probability_weight) >= level)[1L]]
}

compute_nonregular_tan_population_bounds <- function(
    specification = make_nonregular_observed_dgp_spec(),
    Lambda = specification$sensitivity_lambda,
    nodes = 100000L
) {
  Lambda <- as.numeric(Lambda)
  if (length(Lambda) != 1L || !is.finite(Lambda) || Lambda < 1) {
    stop("Lambda must be one finite value of at least one.")
  }
  nodes <- as.integer(nodes)
  if (length(nodes) != 1L || is.na(nodes) || nodes < 1000L) {
    stop("nodes must be one integer of at least 1000.")
  }

  x_values <- 0:1
  x_weights <- c(
    1 - specification$probability_x1,
    specification$probability_x1
  )
  tan_x <- array(NA_real_, c(2L, 2L, 2L))
  tan_xz <- array(NA_real_, c(2L, 2L, 2L))
  maximum_mass_error <- 0

  for (x_index in seq_along(x_values)) {
    x <- x_values[x_index]
    nuisance <- nonregular_observed_nuisances(x, specification)
    e <- nuisance$e
    pi0 <- nuisance$pi0
    pi1 <- nuisance$pi1
    marginal_treatment_probability <- e * pi0 + (1 - e) * pi1

    for (treatment in 0:1) {
      baseline <- nonregular_observed_baseline_rule(
        x, treatment, specification, nodes = nodes
      )
      y <- baseline$y
      mass0 <- baseline$weight / sum(baseline$weight)
      gamma <- nonregular_observed_gamma(
        y, rep(x, length(y)), treatment, specification
      )
      mass1 <- mass0 * gamma
      mass1 <- mass1 / sum(mass1)
      maximum_mass_error <- max(
        maximum_mass_error,
        abs(sum(mass0) - 1), abs(sum(mass1) - 1)
      )

      probability_t0 <- if (treatment == 1L) pi0 else 1 - pi0
      probability_t1 <- if (treatment == 1L) pi1 else 1 - pi1
      probability_t <- if (treatment == 1L) {
        marginal_treatment_probability
      } else {
        1 - marginal_treatment_probability
      }
      mixture_weight0 <- e * probability_t0 / probability_t
      mixture_mass <- mixture_weight0 * mass0 +
        (1 - mixture_weight0) * mass1
      mixture_mass <- mixture_mass / sum(mixture_mass)

      tan_x[x_index, treatment + 1L, ] <-
        observed_dgp_standard_msm_arm(
          y, mixture_mass, marginal_treatment_probability,
          treatment, Lambda
        )
      tan_z0 <- observed_dgp_standard_msm_arm(
        y, mass0, pi0, treatment, Lambda
      )
      tan_z1 <- observed_dgp_standard_msm_arm(
        y, mass1, pi1, treatment, Lambda
      )
      tan_xz[x_index, treatment + 1L, ] <-
        e * tan_z0 + (1 - e) * tan_z1
    }
  }

  summarize <- function(values, method) {
    mu0 <- colSums(values[, 1L, , drop = FALSE][, 1L, ] * x_weights)
    mu1 <- colSums(values[, 2L, , drop = FALSE][, 1L, ] * x_weights)
    data.frame(
      method = method,
      Lambda = Lambda,
      epsilon = 0,
      mu1_lower = mu1[1L],
      mu1_upper = mu1[2L],
      mu0_lower = mu0[1L],
      mu0_upper = mu0[2L],
      ate_lower = mu1[1L] - mu0[2L],
      ate_upper = mu1[2L] - mu0[1L],
      stringsAsFactors = FALSE
    )
  }
  bounds <- rbind(
    summarize(tan_x, "Tan MSM, adjust X"),
    summarize(tan_xz, "Tan MSM, adjust X and Z")
  )
  bounds$ate_width <- bounds$ate_upper - bounds$ate_lower
  list(
    bounds = bounds,
    diagnostics = data.frame(
      design = specification$name,
      Lambda = Lambda,
      integration_nodes_per_stratum_arm = nodes,
      maximum_density_mass_error = maximum_mass_error,
      stringsAsFactors = FALSE
    )
  )
}

assess_nonregular_observed_compatibility <- function(
    specification = make_nonregular_observed_dgp_spec(),
    Lambda,
    nodes = 4000L
) {
  if (!exists("compute_projected_envelopes", mode = "function")) {
    stop("Source R/sharp_bounds_plugin.R before assessing compatibility.")
  }
  rows <- list()
  counter <- 1L
  for (treatment in 0:1) {
    for (x_value in 0:1) {
      rule <- nonregular_observed_baseline_rule(
        x_value, treatment, specification, nodes
      )
      x <- rep(x_value, length(rule$y))
      nuisance <- nonregular_observed_nuisances(x, specification)
      gamma <- nonregular_observed_gamma(
        rule$y, x, treatment, specification
      )
      envelope <- compute_projected_envelopes(
        nuisance = nuisance,
        gamma = gamma,
        t_value = treatment,
        Lambda0 = Lambda,
        Lambda1 = Lambda,
        envelope_smoothing_epsilon = 0,
        gamma_compatibility_action = "stop"
      )
      Lbar <- sum(rule$weight * envelope$L)
      Ubar <- sum(rule$weight * envelope$U)
      rows[[counter]] <- data.frame(
        design = specification$name,
        treatment = treatment,
        X = x_value,
        Lambda = Lambda,
        gamma_violation_probability = sum(
          rule$weight * envelope$gamma_incompatible_original
        ),
        envelope_crossing_probability = sum(
          rule$weight * (envelope$L > envelope$U + 1e-10)
        ),
        Lbar = Lbar,
        Ubar = Ubar,
        conditional_mean_compatible =
          Lbar <= 1 + 1e-10 && Ubar >= 1 - 1e-10,
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }
  detail <- do.call(rbind, rows)
  rownames(detail) <- NULL
  list(
    compatible = all(
      detail$gamma_violation_probability <= 1e-10 &
        detail$envelope_crossing_probability <= 1e-10 &
        detail$conditional_mean_compatible
    ),
    detail = detail
  )
}

nonregular_observed_optimizer_intervals <- function(
    x,
    treatment,
    specification
) {
  location <- nonregular_observed_location(x, treatment, specification)
  empty <- c(lower = NA_real_, upper = NA_real_)
  if (x != 0 || !nonregular_observed_is_nonunique(specification)) {
    return(list(lower = empty, upper = empty))
  }
  if (nonregular_observed_is_atomic(specification)) {
    atoms <- location + specification$atomic_residual_values
    return(list(
      lower = c(lower = atoms[1L], upper = atoms[2L]),
      upper = c(lower = atoms[2L], upper = atoms[3L])
    ))
  }
  list(
    lower = c(
      lower = location + specification$low_interval[2L],
      upper = location + specification$high_interval[1L]
    ),
    upper = empty
  )
}

compute_nonregular_observed_population <- function(
    specification = make_nonregular_observed_dgp_spec(),
    Lambda = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = 0.01,
    nodes = 100000L
) {
  if (!exists("compute_projected_envelopes", mode = "function")) {
    stop("Source R/sharp_bounds_plugin.R before computing the population law.")
  }
  Lambda <- as.numeric(Lambda)
  arm_rows <- list()
  diagnostic_rows <- list()
  counter <- 1L
  for (treatment in 0:1) {
    lower_total <- upper_total <- 0
    for (x_value in 0:1) {
      rule <- nonregular_observed_baseline_rule(
        x_value, treatment, specification, nodes
      )
      x <- rep(x_value, length(rule$y))
      nuisance <- nonregular_observed_nuisances(x, specification)
      gamma <- nonregular_observed_gamma(
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
      tau <- if (Ubar > Lbar) (Ubar - 1) / (Ubar - Lbar) else NA_real_
      lower_level <- 1 - tau
      upper_level <- tau
      weighted_mass <- rule$weight * envelope$D
      optimizer_intervals <- nonregular_observed_optimizer_intervals(
        x_value, treatment, specification
      )
      if (all(is.finite(optimizer_intervals$lower))) {
        q_lower <- mean(optimizer_intervals$lower)
      } else {
        q_lower <- nonregular_weighted_quantile(
          rule$y, weighted_mass, lower_level
        )
      }
      q_upper <- if (all(is.finite(optimizer_intervals$upper))) {
        mean(optimizer_intervals$upper)
      } else {
        nonregular_weighted_quantile(
          rule$y, weighted_mass, upper_level
        )
      }
      lambda_lower <- envelope$L + envelope$D * (rule$y <= q_lower)
      lambda_upper <- envelope$U - envelope$D * (rule$y <= q_upper)
      pi0 <- nuisance$pi0[1L]
      baseline_treatment_probability <- if (treatment == 1L) {
        pi0
      } else {
        1 - pi0
      }
      missing_treatment_probability <- 1 - baseline_treatment_probability
      conditional_lower <- baseline_treatment_probability *
        sum(rule$weight * rule$y) + missing_treatment_probability *
        sum(rule$weight * lambda_lower * rule$y)
      conditional_upper <- baseline_treatment_probability *
        sum(rule$weight * rule$y) + missing_treatment_probability *
        sum(rule$weight * lambda_upper * rule$y)
      x_probability <- if (x_value == 1L) {
        specification$probability_x1
      } else {
        1 - specification$probability_x1
      }
      lower_total <- lower_total + x_probability * conditional_lower
      upper_total <- upper_total + x_probability * conditional_upper
      diagnostic_rows[[counter]] <- data.frame(
        design = specification$name,
        treatment = treatment,
        X = x_value,
        Lambda = Lambda,
        compatible = compatible,
        gamma_violation_probability = sum(
          rule$weight * envelope$gamma_incompatible_original
        ),
        Lbar = Lbar,
        Ubar = Ubar,
        tau = tau,
        lower_level = lower_level,
        upper_level = upper_level,
        q_lower = q_lower,
        q_upper = q_upper,
        branch_tie_probability = if (x_value == 0) 1 else 0,
        nonunique_lower_threshold = all(is.finite(
          optimizer_intervals$lower
        )),
        nonunique_upper_threshold = all(is.finite(
          optimizer_intervals$upper
        )),
        # Retain the original lower-interval names for backward compatibility.
        nonunique_interval_lower = optimizer_intervals$lower[["lower"]],
        nonunique_interval_upper = optimizer_intervals$lower[["upper"]],
        lower_optimizer_interval_lower =
          optimizer_intervals$lower[["lower"]],
        lower_optimizer_interval_upper =
          optimizer_intervals$lower[["upper"]],
        upper_optimizer_interval_lower =
          optimizer_intervals$upper[["lower"]],
        upper_optimizer_interval_upper =
          optimizer_intervals$upper[["upper"]],
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
    arm_rows[[as.character(treatment)]] <- c(
      lower = lower_total, upper = upper_total
    )
  }
  endpoints <- data.frame(
    design = specification$name,
    Lambda = Lambda,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    mu1_lower = arm_rows[["1"]][["lower"]],
    mu1_upper = arm_rows[["1"]][["upper"]],
    mu0_lower = arm_rows[["0"]][["lower"]],
    mu0_upper = arm_rows[["0"]][["upper"]],
    stringsAsFactors = FALSE
  )
  endpoints$ate_lower <- endpoints$mu1_lower - endpoints$mu0_upper
  endpoints$ate_upper <- endpoints$mu1_upper - endpoints$mu0_lower
  diagnostics <- do.call(rbind, diagnostic_rows)
  rownames(diagnostics) <- NULL
  list(
    endpoints = endpoints,
    diagnostics = diagnostics,
    compatible = all(diagnostics$compatible),
    population_nonunique_threshold_probability =
      (1 - specification$probability_x1) *
      as.numeric(nonregular_observed_is_nonunique(specification))
  )
}
