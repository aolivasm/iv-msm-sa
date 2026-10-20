# Observed-data DGP for fully estimated exact and local-to-envelope-tie
# experiments.
#
# The law is designed to isolate nonregularity of the hard max/min envelopes
# while retaining a nonzero marginal first stage. There are six tie-design
# strata and one regular stratum. At the limiting law of a tie-design stratum,
#
#   pi_0(x) = P(T = 1 | X = x, Z = 0)
#           = P(T = 1 | X = x, Z = 1) = pi_1(x)
#
# and the conditional outcome law is the same for Z = 0 and Z = 1. Hence
# gamma_t(y, x) = 1, and at Lambda > 1 both envelope branches meet exactly.
# The default exact-tie design alternates between high and low treatment
# propensity so that the kink matters for both treatment arms. The remaining
# stratum has a strong first stage and gives the instrument a nonzero marginal
# first stage.
#
# The optional instrument-specific tie probabilities also permit a triangular
# array in which pi_1 - pi_0 is of order n^{-1/2}. All nuisance quantities are
# estimated in the simulation. The population functions below are used only
# to compute Monte Carlo coverage targets.

make_fully_estimated_ate_envelope_tie_spec <- function(
    sensitivity_lambda = 3,
    tie_strata = 6L,
    tie_mass = 0.9,
    probability_z1 = 0.5,
    low_treatment_probability = 0.1,
    high_treatment_probability = 0.9,
    tie_treatment_probabilities = NULL,
    tie_treatment_probabilities_z1 = NULL,
    treatment_effect = 0.75,
    outcome_x_shifts = NULL,
    outcome_distribution = c("uniform", "two_cluster_mixture"),
    lower_cluster_probability = 0.2,
    cluster_width = 0.01,
    noise_covariates = 0L
) {
  outcome_distribution <- match.arg(outcome_distribution)
  sensitivity_lambda <- as.numeric(sensitivity_lambda)
  tie_strata <- as.integer(tie_strata)
  noise_covariates <- as.integer(noise_covariates)
  numeric_values <- c(
    sensitivity_lambda, tie_mass, probability_z1,
    low_treatment_probability, high_treatment_probability,
    treatment_effect, lower_cluster_probability, cluster_width
  )
  if (
    length(tie_strata) != 1L || is.na(tie_strata) || tie_strata < 2L ||
      length(noise_covariates) != 1L || is.na(noise_covariates) ||
      noise_covariates < 0L ||
      any(!is.finite(numeric_values)) || sensitivity_lambda <= 1 ||
      tie_mass <= 0 || tie_mass >= 1 ||
      probability_z1 <= 0 || probability_z1 >= 1 ||
      low_treatment_probability <= 0 || low_treatment_probability >= 0.5 ||
      high_treatment_probability <= 0.5 || high_treatment_probability >= 1 ||
      low_treatment_probability >= high_treatment_probability ||
      lower_cluster_probability <= 0 || lower_cluster_probability >= 1 ||
      cluster_width <= 0 || cluster_width >= 0.5
  ) {
    stop("Invalid fully estimated ATE envelope-tie specification.")
  }

  stratum_count <- tie_strata + 1L
  if (is.null(outcome_x_shifts)) {
    outcome_x_shifts <- seq(-0.3, 0.3, length.out = stratum_count)
  }
  outcome_x_shifts <- as.numeric(outcome_x_shifts)
  if (length(outcome_x_shifts) != stratum_count ||
      any(!is.finite(outcome_x_shifts))) {
    stop("outcome_x_shifts must have one finite value per stratum.")
  }

  if (is.null(tie_treatment_probabilities)) {
    tie_probabilities_z0 <- rep(
      c(high_treatment_probability, low_treatment_probability),
      length.out = tie_strata
    )
  } else {
    tie_probabilities_z0 <- as.numeric(tie_treatment_probabilities)
    if (length(tie_probabilities_z0) != tie_strata ||
        any(!is.finite(tie_probabilities_z0)) ||
        any(tie_probabilities_z0 < low_treatment_probability) ||
        any(tie_probabilities_z0 > high_treatment_probability)) {
      stop(paste0(
        "tie_treatment_probabilities must contain one valid probability ",
        "per tie stratum."
      ))
    }
  }
  tie_probabilities_z1 <- if (is.null(tie_treatment_probabilities_z1)) {
    tie_probabilities_z0
  } else {
    as.numeric(tie_treatment_probabilities_z1)
  }
  if (
    length(tie_probabilities_z1) != tie_strata ||
      any(!is.finite(tie_probabilities_z1)) ||
      any(tie_probabilities_z1 <= 0) ||
      any(tie_probabilities_z1 >= 1)
  ) {
    stop(paste0(
      "tie_treatment_probabilities_z1 must contain one probability ",
      "strictly between zero and one per tie stratum."
    ))
  }
  x_probabilities <- c(
    rep(tie_mass / tie_strata, tie_strata),
    1 - tie_mass
  )
  pi0 <- c(tie_probabilities_z0, low_treatment_probability)
  pi1 <- c(tie_probabilities_z1, high_treatment_probability)
  overlap_constant <- min(
    probability_z1, 1 - probability_z1,
    pi0, 1 - pi0, pi1, 1 - pi1
  )

  specification <- list(
    name = if (all(abs(tie_probabilities_z1 - tie_probabilities_z0) < 1e-14)) {
      "fully_estimated_ate_exact_envelope_ties"
    } else {
      "fully_estimated_ate_local_to_envelope_ties"
    },
    sensitivity_lambda = sensitivity_lambda,
    tie_strata = tie_strata,
    regular_stratum = stratum_count,
    x_values = seq_len(stratum_count),
    x_probabilities = x_probabilities,
    tie_mass = tie_mass,
    probability_z1 = probability_z1,
    pi0 = pi0,
    pi1 = pi1,
    tie_treatment_probabilities_z0 = tie_probabilities_z0,
    tie_treatment_probabilities_z1 = tie_probabilities_z1,
    exact_tie_strata = abs(
      tie_probabilities_z1 - tie_probabilities_z0
    ) < 1e-14,
    treatment_effect = treatment_effect,
    outcome_x_shifts = outcome_x_shifts,
    outcome_distribution = outcome_distribution,
    lower_cluster_probability = lower_cluster_probability,
    cluster_width = cluster_width,
    noise_covariates = noise_covariates,
    overlap_constant = overlap_constant,
    marginal_treatment_probability_z0 = sum(x_probabilities * pi0),
    marginal_treatment_probability_z1 = sum(x_probabilities * pi1),
    marginal_first_stage = sum(x_probabilities * (pi1 - pi0)),
    gamma_t = 1
  )
  class(specification) <- c(
    "fully_estimated_ate_envelope_tie_specification", "list"
  )
  specification
}

fully_estimated_ate_envelope_tie_x_index <- function(x, specification) {
  index <- match(as.integer(x), specification$x_values)
  if (anyNA(index)) {
    stop("X contains a value outside the envelope-tie support.")
  }
  index
}

fully_estimated_ate_envelope_tie_x_columns <- function(specification) {
  c(
    paste0("X_tie_", seq_len(specification$tie_strata)),
    if (specification$noise_covariates > 0L) {
      paste0("X_noise_", seq_len(specification$noise_covariates))
    } else {
      character()
    }
  )
}

fully_estimated_ate_envelope_tie_residual_mean <- function(specification) {
  if (identical(specification$outcome_distribution, "uniform")) return(0.5)
  probability <- specification$lower_cluster_probability
  width <- specification$cluster_width
  probability * width / 2 + (1 - probability) * (1 - width / 2)
}

fully_estimated_ate_envelope_tie_residual_quantile <- function(
    probability,
    specification
) {
  probability <- pmin(pmax(as.numeric(probability), 0), 1)
  if (identical(specification$outcome_distribution, "uniform")) {
    return(probability)
  }
  lower_probability <- specification$lower_cluster_probability
  width <- specification$cluster_width
  ifelse(
    probability <= lower_probability,
    width * probability / lower_probability,
    1 - width + width *
      (probability - lower_probability) / (1 - lower_probability)
  )
}

fully_estimated_ate_envelope_tie_residual_partial_mean <- function(
    probability,
    specification
) {
  probability <- pmin(pmax(as.numeric(probability), 0), 1)
  if (identical(specification$outcome_distribution, "uniform")) {
    return(probability^2 / 2)
  }
  lower_probability <- specification$lower_cluster_probability
  width <- specification$cluster_width
  quantile <- fully_estimated_ate_envelope_tie_residual_quantile(
    probability, specification
  )
  ifelse(
    probability <= lower_probability,
    lower_probability * quantile^2 / (2 * width),
    lower_probability * width / 2 +
      (1 - lower_probability) *
        (quantile^2 - (1 - width)^2) / (2 * width)
  )
}

generate_fully_estimated_ate_envelope_tie_data <- function(
    n,
    specification = make_fully_estimated_ate_envelope_tie_spec(),
    seed = NULL,
    include_truth_for_diagnostics = TRUE
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
  x_index <- fully_estimated_ate_envelope_tie_x_index(x, specification)
  z <- stats::rbinom(n, 1L, specification$probability_z1)
  pi0 <- specification$pi0[x_index]
  pi1 <- specification$pi1[x_index]
  treatment_probability <- ifelse(z == 1L, pi1, pi0)
  treatment <- stats::rbinom(n, 1L, treatment_probability)
  location <- specification$treatment_effect * treatment +
    specification$outcome_x_shifts[x_index]
  residual <- if (identical(specification$outcome_distribution, "uniform")) {
    stats::runif(n)
  } else {
    in_lower_cluster <- stats::rbinom(
      n, 1L, specification$lower_cluster_probability
    ) == 1L
    ifelse(
      in_lower_cluster,
      stats::runif(n, 0, specification$cluster_width),
      stats::runif(n, 1 - specification$cluster_width, 1)
    )
  }
  y <- location + residual

  output <- data.frame(
    Y = y,
    T_obs = treatment,
    Z = z,
    X_stratum = x,
    stringsAsFactors = FALSE
  )
  for (level in seq_len(specification$tie_strata)) {
    output[[paste0("X_tie_", level)]] <- as.numeric(x == level)
  }
  if (specification$noise_covariates > 0L) {
    for (index in seq_len(specification$noise_covariates)) {
      output[[paste0("X_noise_", index)]] <- stats::rnorm(n)
    }
  }

  # estimate_sharp_bounds_one_step() currently uses these columns only for
  # Monte Carlo diagnostics when compute_oracle_layers = FALSE. The sole
  # estimated mode is rkhs; these values never enter its nuisance fits,
  # point estimates, influence-function signals, or standard errors.
  if (isTRUE(include_truth_for_diagnostics)) {
    output$e_true <- rep(1 - specification$probability_z1, n)
    output$pi0_true <- pi0
    output$pi1_true <- pi1
    output$gamma0_oracle <- rep(1, n)
    output$gamma1_oracle <- rep(1, n)
  }
  output
}

fully_estimated_ate_envelope_tie_arm_endpoint <- function(
    treatment,
    endpoint = c("lower", "upper"),
    specification = make_fully_estimated_ate_envelope_tie_spec(),
    envelope_smoothing_epsilon = 0
) {
  endpoint <- match.arg(endpoint)
  treatment <- as.integer(treatment)
  if (length(treatment) != 1L || is.na(treatment) ||
      !treatment %in% 0:1) {
    stop("treatment must equal zero or one.")
  }
  if (!exists("compute_projected_envelopes", mode = "function")) {
    stop("Source R/sharp_bounds_plugin.R before computing population bounds.")
  }

  nuisance <- list(
    e = rep(1 - specification$probability_z1,
            length(specification$x_values)),
    pi0 = specification$pi0,
    pi1 = specification$pi1
  )
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = rep(1, length(specification$x_values)),
    t_value = treatment,
    Lambda0 = specification$sensitivity_lambda,
    Lambda1 = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    gamma_compatibility_action = "stop"
  )
  if (any(envelope$gamma_incompatible_original) ||
      any(envelope$L >= envelope$U)) {
    stop("The population law is incompatible at the requested Lambda.")
  }

  location <- specification$outcome_x_shifts +
    specification$treatment_effect * treatment
  outcome_mean <- location +
    fully_estimated_ate_envelope_tie_residual_mean(specification)
  L <- envelope$L
  U <- envelope$U
  if (endpoint == "lower") {
    cutoff_probability <- (1 - L) / (U - L)
    cutoff_residual <- fully_estimated_ate_envelope_tie_residual_quantile(
      cutoff_probability, specification
    )
    lower_partial_moment <- location * cutoff_probability +
      fully_estimated_ate_envelope_tie_residual_partial_mean(
        cutoff_probability, specification
      )
    tilted_mean <- L * outcome_mean +
      (U - L) * lower_partial_moment
  } else {
    cutoff_probability <- (U - 1) / (U - L)
    cutoff_residual <- fully_estimated_ate_envelope_tie_residual_quantile(
      cutoff_probability, specification
    )
    lower_partial_moment <- location * cutoff_probability +
      fully_estimated_ate_envelope_tie_residual_partial_mean(
        cutoff_probability, specification
      )
    tilted_mean <- U * outcome_mean +
      (L - U) * lower_partial_moment
  }
  r0 <- if (treatment == 1L) {
    specification$pi0
  } else {
    1 - specification$pi0
  }
  conditional_endpoint <- r0 * outcome_mean + (1 - r0) * tilted_mean

  list(
    value = sum(specification$x_probabilities * conditional_endpoint),
    conditional_endpoint = conditional_endpoint,
    L = L,
    U = U,
    cutoff = location + cutoff_residual,
    candidate_lower = envelope$candidate_lower,
    candidate_upper = envelope$candidate_upper,
    exact_lower_branch_tie = abs(
      envelope$candidate_lower - 1 / specification$sensitivity_lambda
    ) < 1e-12,
    exact_upper_branch_tie = abs(
      envelope$candidate_upper - specification$sensitivity_lambda
    ) < 1e-12
  )
}

fully_estimated_ate_envelope_tie_population <- function(
    specification = make_fully_estimated_ate_envelope_tie_spec(),
    envelope_smoothing_epsilon = 0
) {
  arm <- lapply(0:1, function(treatment) {
    lapply(c("lower", "upper"), function(endpoint) {
      fully_estimated_ate_envelope_tie_arm_endpoint(
        treatment = treatment,
        endpoint = endpoint,
        specification = specification,
        envelope_smoothing_epsilon = envelope_smoothing_epsilon
      )
    })
  })
  mu0_lower <- arm[[1L]][[1L]]$value
  mu0_upper <- arm[[1L]][[2L]]$value
  mu1_lower <- arm[[2L]][[1L]]$value
  mu1_upper <- arm[[2L]][[2L]]$value
  data.frame(
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    mu1_lower = mu1_lower,
    mu1_upper = mu1_upper,
    mu0_lower = mu0_lower,
    mu0_upper = mu0_upper,
    ate_lower = mu1_lower - mu0_upper,
    ate_upper = mu1_upper - mu0_lower,
    ate_width = (mu1_upper - mu0_lower) - (mu1_lower - mu0_upper),
    reference_ate = specification$treatment_effect,
    stringsAsFactors = FALSE
  )
}

fully_estimated_ate_envelope_tie_population_grid <- function(
    specification = make_fully_estimated_ate_envelope_tie_spec(),
    epsilon_grid = c(0, 0.01)
) {
  epsilon_grid <- sort(unique(as.numeric(epsilon_grid)))
  if (length(epsilon_grid) < 2L || any(!is.finite(epsilon_grid)) ||
      any(epsilon_grid < 0) || !any(epsilon_grid == 0)) {
    stop("epsilon_grid must include zero and another nonnegative value.")
  }
  output <- do.call(rbind, lapply(epsilon_grid, function(epsilon) {
    fully_estimated_ate_envelope_tie_population(
      specification = specification,
      envelope_smoothing_epsilon = epsilon
    )
  }))
  hard <- output[output$envelope_smoothing_epsilon == 0, , drop = FALSE]
  output$ate_lower_shift_from_hard <- output$ate_lower - hard$ate_lower
  output$ate_upper_shift_from_hard <- output$ate_upper - hard$ate_upper
  output$ate_width_increase_from_hard <- output$ate_width - hard$ate_width
  rownames(output) <- NULL
  output
}
