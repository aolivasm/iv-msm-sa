# Selected complete observed-data law from the Lambda <= 4 DGP search.
#
# The law has an exact nonunique upper threshold for the treated-arm endpoint
# at Lambda = 3.5.  It was selected subject to P(T=1 | Z=0) >= 0.10 to avoid
# turning the experiment into a near-positivity counterexample.

if (!exists("make_multibranch_nonunique_threshold_spec", mode = "function")) {
  if (!file.exists("R/multibranch_nonunique_threshold_pilot.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/multibranch_nonunique_threshold_pilot.R")
}

make_adverse_fivefold_nonunique_threshold_spec <- function() {
  make_multibranch_nonunique_threshold_spec(
    sensitivity_lambda = 3.5,
    strata = 3L,
    candidate_design = "complete_binary",
    endpoint_direction = "upper",
    treatment = 1L,
    probability_z0 = 0.2725077,
    pi0 = 0.10845878,
    pi1 = 0.6041288,
    lower_atom = 0,
    boundary_atom = 1,
    tail_atom = 4.866573,
    tail_probability = 0.05,
    threshold_gap_fractions = c(0.0001, 0.9999),
    gamma_atoms = c(1.4717060, 0.7117858, 1.772809)
  )
}

build_adverse_fivefold_nonunique_threshold_oracle <- function(
    specification = make_adverse_fivefold_nonunique_threshold_spec()
) {
  build_multibranch_nonunique_threshold_oracle(specification)
}

adverse_fivefold_x_columns <- function(specification) {
  paste0("X", seq_len(specification$strata - 1L))
}

generate_adverse_fivefold_nonunique_threshold_data <- function(
    n,
    specification = make_adverse_fivefold_nonunique_threshold_spec(),
    seed = NULL,
    oracle = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (is.null(oracle)) {
    oracle <- build_adverse_fivefold_nonunique_threshold_oracle(
      specification
    )
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  category_index <- sample.int(
    nrow(oracle$categories), n, replace = TRUE,
    prob = oracle$categories$probability
  )
  sampled <- oracle$categories[category_index, , drop = FALSE]
  output <- data.frame(
    Y = sampled$Y,
    T_obs = as.integer(sampled$T_obs),
    Z = as.integer(sampled$Z),
    X = as.integer(sampled$X),
    atom = as.integer(sampled$atom),
    category_index = category_index,
    e_true = specification$probability_z0,
    pi0_true = specification$pi0,
    pi1_true = specification$pi1,
    gamma0_oracle = specification$gamma_atoms[sampled$atom],
    gamma1_oracle = specification$gamma_atoms[sampled$atom],
    prob_Z = 1 - specification$probability_z0,
    prob_T = ifelse(
      sampled$Z == 0L, specification$pi0, specification$pi1
    ),
    stringsAsFactors = FALSE
  )
  for (stratum in seq_len(specification$strata - 1L)) {
    output[[paste0("X", stratum)]] <- as.numeric(output$X == stratum)
  }
  attr(output, "adverse_fivefold_nonunique_threshold_specification") <-
    specification
  output
}

audit_adverse_fivefold_nonunique_threshold_dgp <- function(
    specification = make_adverse_fivefold_nonunique_threshold_spec(),
    oracle = NULL
) {
  if (is.null(oracle)) {
    oracle <- build_adverse_fivefold_nonunique_threshold_oracle(
      specification
    )
  }
  baseline_probability <- oracle$atom_probability
  comparison_probability <-
    baseline_probability * specification$gamma_atoms
  covariance <- stats::cov.wt(
    oracle$score,
    wt = oracle$categories$probability,
    center = oracle$population_mean,
    method = "ML"
  )$cov
  correlation <- stats::cov2cor(covariance)
  data.frame(
    Lambda = specification$sensitivity_lambda,
    threshold_level = oracle$threshold_level,
    target = oracle$target,
    strata = specification$strata,
    candidates = specification$candidates,
    probability_z0 = specification$probability_z0,
    pi0 = specification$pi0,
    pi1 = specification$pi1,
    first_stage = specification$pi1 - specification$pi0,
    minimum_treatment_probability = min(
      specification$pi0, specification$pi1
    ),
    maximum_direct_normalization_residual = max(
      oracle$direct_normalization_residual
    ),
    maximum_direct_target_range =
      max(oracle$direct_target) - min(oracle$direct_target),
    maximum_selection_risk_range =
      max(oracle$population_risk) - min(oracle$population_risk),
    maximum_gamma_normalization_residual = abs(
      sum(comparison_probability) - 1
    ),
    minimum_score_correlation = min(
      correlation[upper.tri(correlation)]
    ),
    stringsAsFactors = FALSE
  )
}
