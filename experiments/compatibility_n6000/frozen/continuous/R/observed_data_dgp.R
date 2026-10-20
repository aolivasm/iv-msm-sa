# Observed-data-first simulation laws for the Dual IV--MSM.
#
# These generators specify the observed law
#
#   P(X) P(Z | X) P(T | X, Z) P(Y | T, X, Z)
#
# directly. They deliberately do not generate U, Y0, Y1, or a distinguished
# full-data extension. Population simulation targets are the sharp bound
# endpoints implied by the observed law at a candidate sensitivity parameter.
#
# For continuous outcomes, the Z = 0 conditional outcome density is a normal
# residual truncated to [-2, 2]. The Z = 1 density is its normalized
# exponential tilt. Consequently gamma_t(y, x) is known exactly, bounded on
# the observed support, and normalized by construction.

observed_dgp_clip_probability <- function(value, epsilon = 1e-8) {
  pmin(1 - epsilon, pmax(epsilon, value))
}

observed_dgp_rtruncnorm <- function(n, mean, lower = -2, upper = 2) {
  n <- as.integer(n)
  mean <- rep_len(as.numeric(mean), n)
  if (n < 1L || any(!is.finite(mean)) || !is.finite(lower) ||
      !is.finite(upper) || lower >= upper) {
    stop("Invalid truncated-normal sampling arguments.")
  }
  lower_probability <- stats::pnorm(lower, mean = mean)
  upper_probability <- stats::pnorm(upper, mean = mean)
  probability <- stats::runif(n, lower_probability, upper_probability)
  stats::qnorm(probability, mean = mean)
}

observed_dgp_dtruncnorm <- function(value, mean = 0, lower = -2, upper = 2) {
  value <- as.numeric(value)
  mean <- rep_len(as.numeric(mean), length(value))
  normalizer <- stats::pnorm(upper, mean = mean) -
    stats::pnorm(lower, mean = mean)
  ifelse(
    value >= lower & value <= upper,
    stats::dnorm(value, mean = mean) / normalizer,
    0
  )
}

observed_dgp_rtruncated_standard_normal <- function(n) {
  observed_dgp_rtruncnorm(n, mean = 0, lower = -2, upper = 2)
}

observed_dgp_x_columns <- function(dimension) {
  if (dimension == 1L) "X" else paste0("X", seq_len(dimension))
}

make_observed_data_dgp_spec <- function(
    design = c("continuous_scalar_main", "continuous_x5_nonlinear",
               "continuous_x5_instrument_strength",
               "continuous_scalar_strong_tilt",
               "continuous_scalar_compatibility_separated",
               "continuous_x5_compatibility_separated")
) {
  design <- match.arg(design)
  switch(
    design,
    continuous_scalar_main = list(
      name = design,
      outcome_type = "continuous",
      dimension = 1L,
      nonlinear = FALSE,
      tilt_scale = 0.075,
      observed_treatment_shift = 0.75,
      instrument_log_odds_strength = 1.10,
      residual_lower = -2,
      residual_upper = 2,
      purpose = "primary scalar-X endpoint-inference design"
    ),
    continuous_x5_nonlinear = list(
      name = design,
      outcome_type = "continuous",
      dimension = 5L,
      nonlinear = TRUE,
      # The nonlinear treatment index reaches more extreme values than in the
      # scalar design. A slightly smaller tilt keeps the observable
      # compatibility frontier just below Lambda = 1.5.
      tilt_scale = 0.045,
      observed_treatment_shift = 0.75,
      instrument_log_odds_strength = 1.10,
      residual_lower = -2,
      residual_upper = 2,
      purpose = "five-covariate flexible-nuisance design"
    ),
    continuous_x5_instrument_strength = list(
      name = design,
      outcome_type = "continuous",
      dimension = 5L,
      nonlinear = TRUE,
      tilt_scale = 0.16,
      observed_treatment_shift = 0.75,
      instrument_log_odds_strength = 1.10,
      treatment_intercept = -0.15,
      treatment_feature_scale = 0.40,
      x5_treatment_coefficient = 0.18,
      residual_lower = -2,
      residual_upper = 2,
      purpose = paste0(
        "five-covariate instrument-strength design in which X5 is also ",
        "a treatment--outcome confounder"
      )
    ),
    continuous_scalar_strong_tilt = list(
      name = design,
      outcome_type = "continuous",
      dimension = 1L,
      nonlinear = FALSE,
      tilt_scale = 0.12,
      observed_treatment_shift = 0.75,
      instrument_log_odds_strength = 1.10,
      residual_lower = -2,
      residual_upper = 2,
      purpose = "separated observable alternative for compatibility testing"
    ),
    continuous_scalar_compatibility_separated = list(
      name = design,
      outcome_type = "continuous",
      dimension = 1L,
      nonlinear = FALSE,
      flat_compatibility_design = TRUE,
      tilt_scale = 0.30,
      observed_treatment_shift = 0.75,
      instrument_log_odds_strength = 1.10,
      residual_lower = -2,
      residual_upper = 2,
      purpose = paste0(
        "observable compatibility alternative with violations spread over ",
        "the covariate distribution"
      )
    ),
    continuous_x5_compatibility_separated = list(
      name = design,
      outcome_type = "continuous",
      dimension = 5L,
      nonlinear = TRUE,
      flat_compatibility_treatment = TRUE,
      constant_compatibility_tilt_shape = TRUE,
      tilt_scale = 0.30,
      observed_treatment_shift = 0.75,
      instrument_log_odds_strength = 1.10,
      residual_lower = -2,
      residual_upper = 2,
      purpose = paste0(
        "five-covariate nonlinear continuous compatibility alternative ",
        "with violations spread over the covariate distribution"
      )
    )
  )
}

configure_observed_x5_instrument_strength <- function(
    specification,
    instrument_log_odds_strength,
    x5_treatment_coefficient = 0.18,
    treatment_feature_scale = NULL,
    reference_instrument_strength = 1.10,
    reference_treatment_intercept = -0.15
) {
  if (specification$dimension != 5L) {
    stop("Instrument-strength configuration requires a five-covariate DGP.")
  }
  if (is.null(treatment_feature_scale)) {
    treatment_feature_scale <- if (
      is.null(specification$treatment_feature_scale)
    ) 1 else specification$treatment_feature_scale
  }
  values <- c(
    instrument_log_odds_strength, x5_treatment_coefficient,
    treatment_feature_scale, reference_instrument_strength,
    reference_treatment_intercept
  )
  if (any(!is.finite(values)) || instrument_log_odds_strength <= 0) {
    stop("Instrument strength must be positive and all settings finite.")
  }
  if (treatment_feature_scale <= 0) {
    stop("treatment_feature_scale must be positive.")
  }
  specification$instrument_log_odds_strength <-
    instrument_log_odds_strength
  specification$treatment_intercept <- reference_treatment_intercept +
    (reference_instrument_strength - instrument_log_odds_strength) / 2
  specification$x5_treatment_coefficient <- x5_treatment_coefficient
  specification$treatment_feature_scale <- treatment_feature_scale
  specification$instrument_strength_reference <-
    reference_instrument_strength
  specification
}

observed_dgp_validate_x <- function(x, specification) {
  columns <- observed_dgp_x_columns(specification$dimension)
  if (is.data.frame(x)) {
    missing <- setdiff(columns, names(x))
    if (length(missing)) {
      stop("Missing observed-DGP covariates: ", paste(missing, collapse = ", "))
    }
    matrix_value <- as.matrix(x[, columns, drop = FALSE])
  } else {
    matrix_value <- as.matrix(x)
    if (specification$dimension == 1L && ncol(matrix_value) != 1L) {
      matrix_value <- matrix(as.numeric(x), ncol = 1L)
    }
  }
  storage.mode(matrix_value) <- "double"
  if (ncol(matrix_value) != specification$dimension ||
      any(!is.finite(matrix_value))) {
    stop("x has the wrong dimension or contains non-finite values.")
  }
  colnames(matrix_value) <- columns
  matrix_value
}

observed_dgp_features <- function(x, specification) {
  x <- observed_dgp_validate_x(x, specification)
  if (isTRUE(specification$flat_compatibility_design)) {
    x1 <- x[, 1L]
    return(list(
      instrument = 0.20 * x1,
      treatment = rep(0, length(x1)),
      outcome = 0.65 * x1,
      tilt_shape = rep(1, length(x1))
    ))
  }
  if (specification$dimension == 1L) {
    x1 <- x[, 1L]
    return(list(
      instrument = 0.35 * x1,
      treatment = 0.70 * x1,
      outcome = 0.65 * x1 + 0.15 * sin(pi * x1 / 2),
      tilt_shape = 0.75 + 0.25 * sin(pi * x1 / 2)^2
    ))
  }

  x1 <- x[, 1L]
  x2 <- x[, 2L]
  x3 <- x[, 3L]
  x4 <- x[, 4L]
  x5 <- x[, 5L]
  x5_treatment_coefficient <- if (
    is.null(specification$x5_treatment_coefficient)
  ) 0 else specification$x5_treatment_coefficient
  treatment_feature <- 0.42 * x1 - 0.30 * x2 + 0.22 * x3 +
    0.16 * x1 * x2 + 0.14 * sin(x4) +
    x5_treatment_coefficient * x5
  tilt_shape <- 0.75 +
    0.25 * sin((x1 + x2 - x3) / sqrt(3))^2
  if (isTRUE(specification$flat_compatibility_treatment)) {
    treatment_feature <- rep(0, length(x1))
  }
  treatment_feature_scale <- if (
    is.null(specification$treatment_feature_scale)
  ) 1 else specification$treatment_feature_scale
  treatment_feature <- treatment_feature_scale * treatment_feature
  if (isTRUE(specification$constant_compatibility_tilt_shape)) {
    tilt_shape <- rep(1, length(x1))
  }
  list(
    instrument = 0.28 * x1 - 0.20 * x2 + 0.12 * x3,
    treatment = treatment_feature,
    outcome = 0.45 * x1 + 0.30 * x2 - 0.22 * x3 +
      0.16 * x4 * x5 + 0.12 * (x1^2 - 1),
    tilt_shape = tilt_shape
  )
}

observed_dgp_nuisances <- function(x, specification) {
  features <- observed_dgp_features(x, specification)
  treatment_intercept <- if (
    is.null(specification$treatment_intercept)
  ) -0.15 else specification$treatment_intercept
  probability_z1 <- stats::plogis(features$instrument)
  pi0 <- stats::plogis(treatment_intercept + features$treatment)
  pi1 <- stats::plogis(
    treatment_intercept + features$treatment +
      specification$instrument_log_odds_strength
  )
  list(
    e = 1 - probability_z1,
    probability_z1 = probability_z1,
    pi0 = pi0,
    pi1 = pi1,
    outcome_feature = features$outcome,
    tilt_shape = features$tilt_shape
  )
}

observed_dgp_outcome_location <- function(x, treatment, specification) {
  treatment <- as.integer(treatment)
  if (!treatment %in% 0:1) stop("treatment must be zero or one.")
  features <- observed_dgp_features(x, specification)
  (treatment - 0.5) * specification$observed_treatment_shift +
    features$outcome
}

observed_dgp_tilt <- function(x, treatment, specification) {
  treatment <- as.integer(treatment)
  if (!treatment %in% 0:1) stop("treatment must be zero or one.")
  features <- observed_dgp_features(x, specification)
  ifelse(treatment == 1L, 1, -1) * specification$tilt_scale *
    features$tilt_shape
}

observed_dgp_gamma <- function(y, x, treatment, specification) {
  x <- observed_dgp_validate_x(x, specification)
  y <- rep_len(as.numeric(y), nrow(x))
  location <- observed_dgp_outcome_location(x, treatment, specification)
  tilt <- observed_dgp_tilt(x, treatment, specification)
  residual <- y - location
  lower <- specification$residual_lower
  upper <- specification$residual_upper
  base_normalizer <- stats::pnorm(upper) - stats::pnorm(lower)
  tilted_normalizer <- stats::pnorm(upper, mean = tilt) -
    stats::pnorm(lower, mean = tilt)
  exp(tilt * residual - tilt^2 / 2) *
    base_normalizer / tilted_normalizer
}

generate_observed_iv_data <- function(
    n,
    specification = make_observed_data_dgp_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  columns <- observed_dgp_x_columns(specification$dimension)
  x <- matrix(
    observed_dgp_rtruncated_standard_normal(n * specification$dimension),
    nrow = n,
    ncol = specification$dimension,
    dimnames = list(NULL, columns)
  )
  nuisance <- observed_dgp_nuisances(x, specification)
  z <- stats::rbinom(n, 1L, nuisance$probability_z1)
  probability_t <- ifelse(z == 1L, nuisance$pi1, nuisance$pi0)
  treatment <- stats::rbinom(n, 1L, probability_t)
  location <- numeric(n)
  tilt <- numeric(n)
  for (arm in 0:1) {
    index <- treatment == arm
    if (any(index)) {
      location[index] <- observed_dgp_outcome_location(
        x[index, , drop = FALSE], arm, specification
      )
      tilt[index] <- ifelse(
        z[index] == 1L,
        observed_dgp_tilt(x[index, , drop = FALSE], arm, specification),
        0
      )
    }
  }
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
    stringsAsFactors = FALSE
  )
  output <- cbind(output, as.data.frame(x))
  output$e_true <- nuisance$e
  output$pi0_true <- nuisance$pi0
  output$pi1_true <- nuisance$pi1
  output$gamma0_oracle <- observed_dgp_gamma(y, x, 0L, specification)
  output$gamma1_oracle <- observed_dgp_gamma(y, x, 1L, specification)
  output$prob_Z <- nuisance$probability_z1
  output$prob_T <- probability_t
  attr(output, "observed_data_dgp_specification") <- specification
  output
}

observed_dgp_gauss_legendre <- function(n_nodes, lower = -2, upper = 2) {
  n_nodes <- as.integer(n_nodes)
  if (n_nodes < 2L) stop("n_nodes must be at least two.")
  index <- seq_len(n_nodes - 1L)
  beta <- index / sqrt(4 * index^2 - 1)
  jacobi <- matrix(0, n_nodes, n_nodes)
  jacobi[cbind(index, index + 1L)] <- beta
  jacobi[cbind(index + 1L, index)] <- beta
  decomposition <- eigen(jacobi, symmetric = TRUE)
  ordering <- order(decomposition$values)
  list(
    nodes = (upper - lower) * decomposition$values[ordering] / 2 +
      (upper + lower) / 2,
    weights = (upper - lower) * decomposition$vectors[1L, ordering]^2
  )
}

observed_dgp_trapezoid_rule <- function(lower, upper, n_nodes) {
  nodes <- seq(lower, upper, length.out = as.integer(n_nodes))
  step <- nodes[2L] - nodes[1L]
  weights <- rep(step, length(nodes))
  weights[c(1L, length(nodes))] <- step / 2
  list(nodes = nodes, weights = weights)
}

observed_dgp_smooth_max_below <- function(first, second, epsilon) {
  if (epsilon == 0) return(pmax(first, second))
  bandwidth <- epsilon / log(2)
  maximum <- pmax(first, second)
  maximum + bandwidth * log(
    (exp((first - maximum) / bandwidth) +
       exp((second - maximum) / bandwidth)) / 2
  )
}

observed_dgp_smooth_min_above <- function(first, second, epsilon) {
  if (epsilon == 0) return(pmin(first, second))
  bandwidth <- epsilon / log(2)
  minimum <- pmin(first, second)
  minimum - bandwidth * log(
    (exp(-(first - minimum) / bandwidth) +
       exp(-(second - minimum) / bandwidth)) / 2
  )
}

observed_dgp_envelopes <- function(
    gamma, pi0, pi1, treatment, Lambda0, Lambda1, epsilon = 0
) {
  if (treatment == 1L) {
    gamma_lower <- ((1 - pi0) / Lambda0 + pi0) /
      ((1 - pi1) * Lambda1 + pi1)
    gamma_upper <- (Lambda0 * (1 - pi0) + pi0) /
      ((1 - pi1) / Lambda1 + pi1)
    A <- (1 - pi1) / (1 - pi0) * gamma
    B <- (pi1 * gamma - pi0) / (1 - pi0)
  } else {
    gamma_lower <- (pi0 / Lambda0 + 1 - pi0) /
      (pi1 * Lambda1 + 1 - pi1)
    gamma_upper <- (Lambda0 * pi0 + 1 - pi0) /
      (pi1 / Lambda1 + 1 - pi1)
    A <- pi1 / pi0 * gamma
    B <- ((1 - pi1) * gamma - (1 - pi0)) / pi0
  }
  candidate_lower <- A / Lambda1 + B
  candidate_upper <- A * Lambda1 + B
  list(
    L = observed_dgp_smooth_max_below(
      rep(1 / Lambda0, length(gamma)), candidate_lower, epsilon
    ),
    U = observed_dgp_smooth_min_above(
      rep(Lambda0, length(gamma)), candidate_upper, epsilon
    ),
    hard_L = pmax(1 / Lambda0, candidate_lower),
    hard_U = pmin(Lambda0, candidate_upper),
    gamma_lower = gamma_lower,
    gamma_upper = gamma_upper,
    candidate_lower = candidate_lower,
    candidate_upper = candidate_upper
  )
}

observed_dgp_tail_moment <- function(y, mass, target, upper_tail = FALSE) {
  if (upper_tail) {
    result <- observed_dgp_tail_moment(
      rev(-y), rev(mass), target, upper_tail = FALSE
    )
    return(list(mass = result$mass, y_moment = -result$y_moment))
  }
  target <- pmax(0, pmin(sum(mass), target))
  if (target <= 0) return(list(mass = 0, y_moment = 0))
  cumulative <- cumsum(mass)
  boundary <- which(cumulative >= target)[1L]
  before <- if (boundary == 1L) 0 else cumulative[boundary - 1L]
  remainder <- target - before
  fraction <- if (mass[boundary] > 0) remainder / mass[boundary] else 0
  preceding <- if (boundary == 1L) numeric() else seq_len(boundary - 1L)
  list(
    mass = before + fraction * mass[boundary],
    y_moment = sum(y[preceding] * mass[preceding]) +
      fraction * y[boundary] * mass[boundary]
  )
}

observed_dgp_extreme_weighted_mean <- function(y, density_mass, L, U) {
  D <- U - L
  lower_base_mass <- sum(density_mass * L)
  target <- 1 - lower_base_mass
  d_mass <- density_mass * D
  if (target < -1e-8 || target > sum(d_mass) + 1e-8 || any(D < -1e-8)) {
    return(list(lower = NA_real_, upper = NA_real_, normalization_error = Inf))
  }
  lower_tail <- observed_dgp_tail_moment(y, d_mass, target)
  upper_tail <- observed_dgp_tail_moment(y, d_mass, target, upper_tail = TRUE)
  base_y <- sum(y * density_mass * L)
  list(
    lower = base_y + lower_tail$y_moment,
    upper = base_y + upper_tail$y_moment,
    normalization_error = max(
      abs(lower_base_mass + lower_tail$mass - 1),
      abs(lower_base_mass + upper_tail$mass - 1)
    )
  )
}

observed_dgp_standard_msm_arm <- function(
    y, density_mass, propensity, treatment, Lambda
) {
  L <- rep(1 / Lambda, length(y))
  U <- rep(Lambda, length(y))
  extreme <- observed_dgp_extreme_weighted_mean(y, density_mass, L, U)
  observed_mean <- sum(y * density_mass)
  if (treatment == 1L) {
    c(
      lower = propensity * observed_mean + (1 - propensity) * extreme$lower,
      upper = propensity * observed_mean + (1 - propensity) * extreme$upper
    )
  } else {
    c(
      lower = (1 - propensity) * observed_mean + propensity * extreme$lower,
      upper = (1 - propensity) * observed_mean + propensity * extreme$upper
    )
  }
}

observed_dgp_covariate_rule <- function(
    specification,
    x_nodes = 120L,
    x_nodes_per_dimension = 7L
) {
  if (specification$dimension == 1L) {
    marginal_rule <- observed_dgp_gauss_legendre(x_nodes, -2, 2)
    marginal_weight <- marginal_rule$weights * observed_dgp_dtruncnorm(
      marginal_rule$nodes, mean = 0
    )
    x_matrix <- matrix(marginal_rule$nodes, ncol = 1L)
    x_weight <- marginal_weight / sum(marginal_weight)
  } else {
    marginal_rule <- observed_dgp_gauss_legendre(
      x_nodes_per_dimension, -2, 2
    )
    marginal_weight <- marginal_rule$weights * observed_dgp_dtruncnorm(
      marginal_rule$nodes, mean = 0
    )
    marginal_weight <- marginal_weight / sum(marginal_weight)
    index_grid <- expand.grid(rep(
      list(seq_along(marginal_rule$nodes)), specification$dimension
    ))
    x_matrix <- matrix(
      marginal_rule$nodes[as.matrix(index_grid)],
      nrow = nrow(index_grid),
      ncol = specification$dimension
    )
    x_weight <- apply(
      matrix(
        marginal_weight[as.matrix(index_grid)],
        nrow = nrow(index_grid),
        ncol = specification$dimension
      ),
      1L,
      prod
    )
    x_weight <- x_weight / sum(x_weight)
  }
  colnames(x_matrix) <- observed_dgp_x_columns(specification$dimension)
  list(x = x_matrix, weights = x_weight)
}

observed_dgp_instrument_strength_metrics <- function(
    specification,
    x_nodes = 160L,
    x_nodes_per_dimension = 9L
) {
  x_rule <- observed_dgp_covariate_rule(
    specification, x_nodes, x_nodes_per_dimension
  )
  x <- x_rule$x
  weight <- x_rule$weights
  nuisance <- observed_dgp_nuisances(x, specification)
  probability_z1 <- nuisance$probability_z1
  probability_z0 <- 1 - probability_z1
  pi0 <- nuisance$pi0
  pi1 <- nuisance$pi1
  treatment_given_x <- probability_z0 * pi0 + probability_z1 * pi1

  incremental_variance <- sum(
    weight * probability_z0 * probability_z1 * (pi1 - pi0)^2
  )
  residual_variance_after_x <- sum(
    weight * treatment_given_x * (1 - treatment_given_x)
  )
  residual_variance_after_xz <- sum(weight * (
    probability_z0 * pi0 * (1 - pi0) +
      probability_z1 * pi1 * (1 - pi1)
  ))
  partial_r_squared <- incremental_variance / residual_variance_after_x

  probability_z <- c(
    z0 = sum(weight * probability_z0),
    z1 = sum(weight * probability_z1)
  )
  treatment_by_observed_z <- c(
    z0 = unname(
      sum(weight * probability_z0 * pi0) / probability_z["z0"]
    ),
    z1 = unname(
      sum(weight * probability_z1 * pi1) / probability_z["z1"]
    )
  )
  cell_probability <- c(
    z0_t0 = sum(weight * probability_z0 * (1 - pi0)),
    z0_t1 = sum(weight * probability_z0 * pi0),
    z1_t0 = sum(weight * probability_z1 * (1 - pi1)),
    z1_t1 = sum(weight * probability_z1 * pi1)
  )
  treatment_intercept <- if (
    is.null(specification$treatment_intercept)
  ) -0.15 else specification$treatment_intercept
  x5_treatment_coefficient <- if (
    is.null(specification$x5_treatment_coefficient)
  ) 0 else specification$x5_treatment_coefficient
  treatment_feature_scale <- if (
    is.null(specification$treatment_feature_scale)
  ) 1 else specification$treatment_feature_scale

  data.frame(
    design = specification$name,
    instrument_log_odds_coefficient =
      specification$instrument_log_odds_strength,
    conditional_instrument_odds_ratio =
      exp(specification$instrument_log_odds_strength),
    treatment_intercept = treatment_intercept,
    treatment_feature_scale = treatment_feature_scale,
    x5_treatment_coefficient = x5_treatment_coefficient,
    standardized_first_stage = sum(weight * (pi1 - pi0)),
    observed_marginal_first_stage =
      unname(treatment_by_observed_z["z1"] - treatment_by_observed_z["z0"]),
    treatment_probability_z0 = unname(treatment_by_observed_z["z0"]),
    treatment_probability_z1 = unname(treatment_by_observed_z["z1"]),
    population_partial_r_squared = partial_r_squared,
    incremental_first_stage_variance = incremental_variance,
    residual_variance_after_x = residual_variance_after_x,
    residual_variance_after_xz = residual_variance_after_xz,
    marginal_treatment_probability = sum(weight * treatment_given_x),
    minimum_treatment_instrument_cell_probability = min(cell_probability),
    pi0_tail_mass = sum(weight * (pi0 < 0.05 | pi0 > 0.95)),
    pi1_tail_mass = sum(weight * (pi1 < 0.05 | pi1 > 0.95)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

observed_dgp_sample_instrument_strength_metrics <- function(
    data,
    specification
) {
  treatment_name <- if ("T" %in% names(data)) {
    "T"
  } else if ("T_obs" %in% names(data)) {
    "T_obs"
  } else {
    "T"
  }
  required <- c(
    treatment_name, "Z", observed_dgp_x_columns(specification$dimension)
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Missing variables for the sample first-stage diagnostic: ",
      paste(missing, collapse = ", "), "."
    )
  }
  treatment <- as.numeric(data[[treatment_name]])
  instrument <- as.numeric(data$Z)
  if (
    any(!is.finite(treatment)) || any(!treatment %in% c(0, 1)) ||
      any(!is.finite(instrument)) || any(!instrument %in% c(0, 1))
  ) {
    stop("T and Z must be finite binary variables.")
  }
  x <- observed_dgp_validate_x(data, specification)
  controls <- cbind(`(Intercept)` = 1, x)
  if (specification$dimension == 5L) {
    controls <- cbind(
      controls,
      X1_X2 = x[, "X1"] * x[, "X2"],
      sin_X4 = sin(x[, "X4"])
    )
  }
  full_design <- cbind(controls, Z = instrument)
  full_fit <- stats::lm.fit(full_design, treatment)
  restricted_fit <- stats::lm.fit(controls, treatment)
  if (full_fit$rank != ncol(full_design)) {
    stop("The sample first-stage design matrix is rank deficient.")
  }

  sample_size <- nrow(full_design)
  parameter_count <- ncol(full_design)
  inverse_cross_product <- chol2inv(chol(crossprod(full_design)))
  score_cross_product <- crossprod(
    full_design,
    full_design * as.numeric(full_fit$residuals)^2
  )
  robust_vcov <- sample_size / (sample_size - parameter_count) *
    inverse_cross_product %*% score_cross_product %*%
    inverse_cross_product
  instrument_index <- ncol(full_design)
  instrument_coefficient <- unname(full_fit$coefficients[instrument_index])
  instrument_standard_error <- sqrt(
    robust_vcov[instrument_index, instrument_index]
  )
  robust_wald_statistic <-
    (instrument_coefficient / instrument_standard_error)^2
  restricted_rss <- sum(restricted_fit$residuals^2)
  full_rss <- sum(full_fit$residuals^2)
  cell_counts <- table(
    factor(instrument, levels = c(0, 1)),
    factor(treatment, levels = c(0, 1))
  )
  probability_treatment_by_z <- vapply(c(0, 1), function(z) {
    mean(treatment[instrument == z])
  }, numeric(1))

  data.frame(
    sample_first_stage_coefficient = instrument_coefficient,
    sample_first_stage_robust_standard_error = instrument_standard_error,
    sample_first_stage_robust_wald_statistic = robust_wald_statistic,
    sample_first_stage_partial_r_squared =
      max(0, (restricted_rss - full_rss) / restricted_rss),
    sample_unadjusted_first_stage =
      probability_treatment_by_z[2L] - probability_treatment_by_z[1L],
    sample_minimum_treatment_instrument_cell_count = min(cell_counts),
    sample_minimum_treatment_instrument_cell_proportion =
      min(cell_counts) / sample_size,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

compute_observed_data_population_bounds <- function(
    specification = make_observed_data_dgp_spec(),
    Lambda,
    epsilon = 0,
    x_nodes = 120L,
    x_nodes_per_dimension = 7L,
    residual_nodes = 4001L,
    tolerance = 1e-7
) {
  if (length(Lambda) != 1L || !is.finite(Lambda) || Lambda < 1) {
    stop("Lambda must be one finite value of at least one.")
  }
  if (epsilon < 0 || (epsilon > 0 && epsilon >= 1 / Lambda)) {
    stop("epsilon must be nonnegative and smaller than 1/Lambda.")
  }

  x_rule <- observed_dgp_covariate_rule(
    specification, x_nodes, x_nodes_per_dimension
  )
  x_matrix <- x_rule$x
  x_weight <- x_rule$weights
  residual_rule <- observed_dgp_trapezoid_rule(
    specification$residual_lower,
    specification$residual_upper,
    residual_nodes
  )
  residual <- residual_rule$nodes

  dual_conditional <- array(NA_real_, c(nrow(x_matrix), 2L, 2L))
  tan_x_conditional <- array(NA_real_, c(nrow(x_matrix), 2L, 2L))
  tan_xz_conditional <- array(NA_real_, c(nrow(x_matrix), 2L, 2L))
  maximum_pointwise_violation <- 0
  maximum_envelope_mean_violation <- 0
  maximum_normalization_error <- 0
  maximum_violation_mass <- 0
  violation_mass_by_x_arm <- matrix(0, nrow(x_matrix), 2L)
  envelope_mean_violation_by_x_arm <- matrix(0, nrow(x_matrix), 2L)

  for (x_index in seq_len(nrow(x_matrix))) {
    x <- x_matrix[x_index, , drop = FALSE]
    nuisance <- observed_dgp_nuisances(x, specification)
    e <- nuisance$e
    pi0 <- nuisance$pi0
    pi1 <- nuisance$pi1
    p_treatment <- e * pi0 + (1 - e) * pi1

    for (treatment in 0:1) {
      location <- observed_dgp_outcome_location(x, treatment, specification)
      tilt <- observed_dgp_tilt(x, treatment, specification)
      y <- location + residual
      f0 <- observed_dgp_dtruncnorm(
        residual, 0, specification$residual_lower,
        specification$residual_upper
      )
      f1 <- observed_dgp_dtruncnorm(
        residual, tilt, specification$residual_lower,
        specification$residual_upper
      )
      mass0 <- f0 * residual_rule$weights
      mass1 <- f1 * residual_rule$weights
      mass0 <- mass0 / sum(mass0)
      mass1 <- mass1 / sum(mass1)
      gamma <- observed_dgp_gamma(y, x[rep(1L, length(y)), , drop = FALSE],
                                  treatment, specification)
      envelope <- observed_dgp_envelopes(
        gamma, pi0, pi1, treatment, Lambda, Lambda, epsilon
      )
      pointwise_violation <- pmax(
        envelope$gamma_lower - gamma,
        gamma - envelope$gamma_upper,
        envelope$hard_L - envelope$hard_U,
        0
      )
      violation_indicator <- pointwise_violation > tolerance
      maximum_pointwise_violation <- max(
        maximum_pointwise_violation, pointwise_violation
      )
      maximum_violation_mass <- max(
        maximum_violation_mass, sum(mass0[violation_indicator])
      )
      violation_mass_by_x_arm[x_index, treatment + 1L] <-
        sum(mass0[violation_indicator])
      mean_L <- sum(mass0 * envelope$hard_L)
      mean_U <- sum(mass0 * envelope$hard_U)
      maximum_envelope_mean_violation <- max(
        maximum_envelope_mean_violation, mean_L - 1, 1 - mean_U, 0
      )
      envelope_mean_violation_by_x_arm[x_index, treatment + 1L] <-
        max(mean_L - 1, 1 - mean_U, 0)
      extreme <- observed_dgp_extreme_weighted_mean(
        y, mass0, envelope$L, envelope$U
      )
      maximum_normalization_error <- max(
        maximum_normalization_error, extreme$normalization_error
      )
      observed_mean0 <- sum(y * mass0)
      if (treatment == 1L) {
        dual_conditional[x_index, treatment + 1L, ] <- c(
          pi0 * observed_mean0 + (1 - pi0) * extreme$lower,
          pi0 * observed_mean0 + (1 - pi0) * extreme$upper
        )
      } else {
        dual_conditional[x_index, treatment + 1L, ] <- c(
          (1 - pi0) * observed_mean0 + pi0 * extreme$lower,
          (1 - pi0) * observed_mean0 + pi0 * extreme$upper
        )
      }

      probability_t0 <- if (treatment == 1L) pi0 else 1 - pi0
      probability_t1 <- if (treatment == 1L) pi1 else 1 - pi1
      probability_t <- if (treatment == 1L) p_treatment else 1 - p_treatment
      mixture_weight0 <- e * probability_t0 / probability_t
      mixture_mass <- mixture_weight0 * mass0 + (1 - mixture_weight0) * mass1
      tan_x_conditional[x_index, treatment + 1L, ] <-
        observed_dgp_standard_msm_arm(
          y, mixture_mass, p_treatment, treatment, Lambda
        )
      tan_z0 <- observed_dgp_standard_msm_arm(
        y, mass0, pi0, treatment, Lambda
      )
      tan_z1 <- observed_dgp_standard_msm_arm(
        y, mass1, pi1, treatment, Lambda
      )
      tan_xz_conditional[x_index, treatment + 1L, ] <-
        e * tan_z0 + (1 - e) * tan_z1
    }
  }

  summarize_method <- function(values, method) {
    mu0 <- colSums(values[, 1L, , drop = FALSE][, 1L, ] * x_weight)
    mu1 <- colSums(values[, 2L, , drop = FALSE][, 1L, ] * x_weight)
    data.frame(
      method = method,
      Lambda = Lambda,
      epsilon = epsilon,
      mu1_lower = mu1[1L],
      mu1_upper = mu1[2L],
      mu0_lower = mu0[1L],
      mu0_upper = mu0[2L],
      ate_lower = mu1[1L] - mu0[2L],
      ate_upper = mu1[2L] - mu0[1L],
      stringsAsFactors = FALSE
    )
  }

  compatible <- maximum_pointwise_violation <= tolerance &&
    maximum_envelope_mean_violation <= tolerance
  dual <- summarize_method(dual_conditional, "Dual IV--MSM")
  if (!compatible) dual[, 4:9] <- NA_real_
  bounds <- rbind(
    dual,
    summarize_method(tan_x_conditional, "Tan MSM, adjust X"),
    summarize_method(tan_xz_conditional, "Tan MSM, adjust X and Z")
  )
  bounds$ate_width <- bounds$ate_upper - bounds$ate_lower

  list(
    bounds = bounds,
    diagnostics = data.frame(
      design = specification$name,
      Lambda = Lambda,
      epsilon = epsilon,
      compatible = compatible,
      maximum_pointwise_violation = maximum_pointwise_violation,
      maximum_pointwise_violation_mass = maximum_violation_mass,
      maximum_population_weighted_pointwise_violation_mass = max(
        colSums(violation_mass_by_x_arm * x_weight)
      ),
      maximum_envelope_mean_violation = maximum_envelope_mean_violation,
      maximum_population_weighted_envelope_mean_violation = max(
        colSums(envelope_mean_violation_by_x_arm * x_weight)
      ),
      maximum_extremizer_normalization_error = maximum_normalization_error,
      stringsAsFactors = FALSE
    )
  )
}
