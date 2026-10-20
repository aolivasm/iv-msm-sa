# Observed-data-first binary-outcome DGP and exact population bounds.

if (!exists("make_observed_data_dgp_spec", mode = "function")) {
  if (!file.exists("R/observed_data_dgp.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/observed_data_dgp.R")
}
if (!exists("binary_arm_endpoint", mode = "function")) {
  if (!file.exists("R/binary_outcome_bounds.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/binary_outcome_bounds.R")
}

make_observed_binary_dgp_spec <- function(
    design = c(
      "binary_scalar_main", "binary_scalar_compatibility_separated",
      "binary_x5_nonlinear", "binary_x5_instrument_strength",
      "binary_x5_compatibility_separated"
    )
) {
  design <- match.arg(design)
  base_design <- switch(
    design,
    binary_scalar_main = "continuous_scalar_main",
    binary_scalar_compatibility_separated =
      "continuous_scalar_compatibility_separated",
    binary_x5_nonlinear = "continuous_x5_nonlinear",
    binary_x5_instrument_strength =
      "continuous_x5_instrument_strength",
    binary_x5_compatibility_separated = "continuous_x5_nonlinear"
  )
  base <- make_observed_data_dgp_spec(base_design)
  base$name <- design
  base$outcome_type <- "binary"
  base$binary_outcome_intercept <- -0.40
  base$binary_observed_treatment_log_odds_shift <- 1.50
  base$binary_tilt_scale <- switch(
    design,
    binary_scalar_main = 0.15,
    binary_scalar_compatibility_separated = 0.75,
    binary_x5_nonlinear = 0.10,
    binary_x5_instrument_strength = 0.35,
    binary_x5_compatibility_separated = 0.75
  )
  if (design == "binary_x5_compatibility_separated") {
    base$flat_compatibility_treatment <- TRUE
    base$constant_compatibility_tilt_shape <- TRUE
  }
  base$purpose <- switch(
    design,
    binary_scalar_main = "primary scalar-X binary endpoint-inference design",
    binary_scalar_compatibility_separated = paste0(
      "binary observable compatibility alternative with violations spread ",
      "over the covariate distribution"
    ),
    binary_x5_nonlinear = "five-covariate nonlinear binary design",
    binary_x5_instrument_strength = paste0(
      "five-covariate binary instrument-strength design in which X5 is ",
      "also a treatment--outcome confounder"
    ),
    binary_x5_compatibility_separated = paste0(
      "five-covariate nonlinear binary compatibility alternative with ",
      "violations spread over the covariate distribution"
    )
  )
  base
}

observed_binary_probabilities <- function(x, treatment, specification) {
  treatment <- as.integer(treatment)
  if (!treatment %in% 0:1) stop("treatment must be zero or one.")
  features <- observed_dgp_features(x, specification)
  baseline <- stats::plogis(
    specification$binary_outcome_intercept +
      specification$binary_observed_treatment_log_odds_shift * treatment +
      features$outcome
  )
  comparison <- stats::plogis(
    stats::qlogis(baseline) +
      ifelse(treatment == 1L, 1, -1) *
      specification$binary_tilt_scale * features$tilt_shape
  )
  list(a0 = baseline, a1 = comparison)
}

observed_binary_gamma <- function(y, x, treatment, specification) {
  probability <- observed_binary_probabilities(x, treatment, specification)
  y <- rep_len(as.integer(y), length(probability$a0))
  ifelse(
    y == 1L,
    probability$a1 / probability$a0,
    (1 - probability$a1) / (1 - probability$a0)
  )
}

generate_observed_binary_iv_data <- function(
    n,
    specification = make_observed_binary_dgp_spec(),
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
  probability_y <- numeric(n)
  for (arm in 0:1) {
    index <- treatment == arm
    if (any(index)) {
      arm_probability <- observed_binary_probabilities(
        x[index, , drop = FALSE], arm, specification
      )
      probability_y[index] <- ifelse(
        z[index] == 1L, arm_probability$a1, arm_probability$a0
      )
    }
  }
  y <- stats::rbinom(n, 1L, probability_y)
  output <- data.frame(Y = y, T_obs = treatment, Z = z)
  output <- cbind(output, as.data.frame(x))
  output$e_true <- nuisance$e
  output$pi0_true <- nuisance$pi0
  output$pi1_true <- nuisance$pi1
  output$gamma0_oracle <- observed_binary_gamma(
    y, x, 0L, specification
  )
  output$gamma1_oracle <- observed_binary_gamma(
    y, x, 1L, specification
  )
  output$prob_Z <- nuisance$probability_z1
  output$prob_T <- probability_t
  output$prob_Y <- probability_y
  attr(output, "observed_binary_dgp_specification") <- specification
  output
}

compute_observed_binary_population_bounds <- function(
    specification = make_observed_binary_dgp_spec(),
    Lambda,
    epsilon = 0.01,
    x_nodes = 160L,
    x_nodes_per_dimension = 7L
) {
  if (length(Lambda) != 1L || !is.finite(Lambda) || Lambda < 1) {
    stop("Lambda must be one finite number of at least one.")
  }
  x_rule <- observed_dgp_covariate_rule(
    specification, x_nodes, x_nodes_per_dimension
  )
  x <- x_rule$x
  weight <- x_rule$weights
  nuisance <- observed_dgp_nuisances(x, specification)
  hard <- smooth <- tan_x <- tan_xz <- iv_only <-
    iv_only_smooth <- list()

  for (treatment in 0:1) {
    key <- as.character(treatment)
    outcome <- observed_binary_probabilities(x, treatment, specification)
    arm <- list(
      r0 = if (treatment == 1L) nuisance$pi0 else 1 - nuisance$pi0,
      r1 = if (treatment == 1L) nuisance$pi1 else 1 - nuisance$pi1,
      a0 = outcome$a0,
      a1 = outcome$a1
    )
    hard[[key]] <- binary_arm_endpoint(arm, Lambda, Lambda, 0, FALSE)
    smooth[[key]] <- binary_arm_endpoint(
      arm, Lambda, Lambda, epsilon, FALSE
    )
    marginal_r <- nuisance$e * arm$r0 + (1 - nuisance$e) * arm$r1
    marginal_a <- (
      nuisance$e * arm$r0 * arm$a0 +
        (1 - nuisance$e) * arm$r1 * arm$a1
    ) / marginal_r
    tan_x[[key]] <- binary_standard_msm_arm(
      marginal_r, marginal_a, Lambda
    )
    tan_z0 <- binary_standard_msm_arm(arm$r0, arm$a0, Lambda)
    tan_z1 <- binary_standard_msm_arm(arm$r1, arm$a1, Lambda)
    tan_xz[[key]] <- list(
      lower = nuisance$e * tan_z0$lower + (1 - nuisance$e) * tan_z1$lower,
      upper = nuisance$e * tan_z0$upper + (1 - nuisance$e) * tan_z1$upper
    )
    iv_only[[key]] <- list(
      lower = pmax(arm$r0 * arm$a0, arm$r1 * arm$a1),
      upper = pmin(
        arm$r0 * arm$a0 + 1 - arm$r0,
        arm$r1 * arm$a1 + 1 - arm$r1
      )
    )
    iv_only_smooth[[key]] <- binary_iv_arm_endpoint(
      arm, epsilon = epsilon, compute_gradient = FALSE
    )
  }

  compatible <- !any(c(
    hard[["0"]]$gamma_violation,
    hard[["1"]]$gamma_violation,
    hard[["0"]]$envelope_crossing,
    hard[["1"]]$envelope_crossing,
    hard[["0"]]$theta_crossing,
    hard[["1"]]$theta_crossing
  ))
  weighted <- function(value) sum(weight * value)
  make_row <- function(method, value, epsilon_value = 0) {
    mu1_lower <- weighted(value[["1"]]$lower)
    mu1_upper <- weighted(value[["1"]]$upper)
    mu0_lower <- weighted(value[["0"]]$lower)
    mu0_upper <- weighted(value[["0"]]$upper)
    data.frame(
      method = method,
      Lambda = Lambda,
      epsilon = epsilon_value,
      mu1_lower = mu1_lower,
      mu1_upper = mu1_upper,
      mu0_lower = mu0_lower,
      mu0_upper = mu0_upper,
      ate_lower = mu1_lower - mu0_upper,
      ate_upper = mu1_upper - mu0_lower,
      stringsAsFactors = FALSE
    )
  }
  hard_row <- make_row("Dual IV--MSM, hard", hard)
  smooth_row <- make_row("Dual IV--MSM, smoothed", smooth, epsilon)
  if (!compatible) {
    hard_row[, 4:9] <- NA_real_
    smooth_row[, 4:9] <- NA_real_
  }
  bounds <- rbind(
    hard_row,
    smooth_row,
    make_row("Tan MSM, adjust X", tan_x),
    make_row("Tan MSM, adjust X and Z", tan_xz),
    make_row("IV only", iv_only),
    make_row("IV only, smoothed outer", iv_only_smooth, epsilon)
  )
  bounds$ate_width <- bounds$ate_upper - bounds$ate_lower

  arm_diagnostic <- function(treatment) {
    fit <- hard[[as.character(treatment)]]
    data.frame(
      treatment = treatment,
      gamma_violation_probability = weighted(fit$gamma_violation),
      envelope_crossing_probability = weighted(fit$envelope_crossing),
      theta_crossing_probability = weighted(fit$theta_crossing),
      stringsAsFactors = FALSE
    )
  }
  list(
    bounds = bounds,
    diagnostics = cbind(
      data.frame(
        design = specification$name,
        Lambda = Lambda,
        compatible = compatible,
        stringsAsFactors = FALSE
      ),
      rbind(arm_diagnostic(0L), arm_diagnostic(1L))
    )
  )
}
