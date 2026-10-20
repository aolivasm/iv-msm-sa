# Local-to-envelope-tie binary DGP.
#
# This file defines a triangular-array stress experiment for nonuniform
# inference near the positive-mass envelope ties in
# R/binary_envelope_tie_dgp.R.  Within the limit-tie strata,
#
#   pi_1(x) - pi_0(x) = s(x) * kappa / sqrt(n * p_x),
#
# where p_x is the probability of one limit-tie stratum and s(x) alternates
# between -1 and 1.  The alternating signs preserve the aggregate first stage,
# while kappa measures the local separation in units of the elementary
# within-stratum propensity-difference standard error.  At kappa = 0 the law
# is the fixed exact-tie DGP.  For every fixed kappa, the sequence approaches
# that tied law at the root-n rate.

if (!exists("make_binary_envelope_tie_spec", mode = "function")) {
  if (!file.exists("R/binary_envelope_tie_dgp.R")) {
    stop("Run this file from the repository root.", call. = FALSE)
  }
  source("R/binary_envelope_tie_dgp.R")
}

make_binary_local_envelope_tie_spec <- function(
    sample_size,
    kappa,
    Lambda = 4,
    tie_strata = 8L,
    tie_mass = 0.85,
    probability_z1 = 0.5,
    regular_treatment_probabilities = c(z0 = 0.2, z1 = 0.8),
    tie_outcome_t0 = NULL,
    tie_outcome_t1 = NULL,
    regular_outcome_probabilities = c(t0 = 0.375, t1 = 0.625)
) {
  sample_size <- as.integer(sample_size)
  kappa <- as.numeric(kappa)
  tie_strata <- as.integer(tie_strata)
  if (length(sample_size) != 1L || is.na(sample_size) || sample_size < 200L) {
    stop("sample_size must be one integer of at least 200.")
  }
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("kappa must be one nonnegative finite number.")
  }
  if (length(tie_strata) != 1L || is.na(tie_strata) ||
      tie_strata < 2L || tie_strata %% 2L != 0L) {
    stop("tie_strata must be one positive even integer of at least two.")
  }

  specification <- make_binary_envelope_tie_spec(
    Lambda = Lambda,
    tie_strata = tie_strata,
    tie_mass = tie_mass,
    probability_z1 = probability_z1,
    tie_treatment_probability = 0.5,
    regular_treatment_probabilities = regular_treatment_probabilities,
    tie_outcome_t0 = tie_outcome_t0,
    tie_outcome_t1 = tie_outcome_t1,
    regular_outcome_probabilities = regular_outcome_probabilities
  )

  tie_stratum_probability <- tie_mass / tie_strata
  local_probability_difference <-
    kappa / sqrt(sample_size * tie_stratum_probability)
  local_sign <- rep(c(-1, 1), length.out = tie_strata)
  pi0 <- 0.5 - local_sign * local_probability_difference / 2
  pi1 <- 0.5 + local_sign * local_probability_difference / 2
  if (any(pi0 <= 0.05) || any(pi0 >= 0.95) ||
      any(pi1 <= 0.05) || any(pi1 >= 0.95)) {
    stop(
      paste0(
        "The requested local perturbation violates the [0.05,0.95] ",
        "propensity range; reduce kappa or increase sample_size."
      )
    )
  }

  specification$name <- "binary_local_to_envelope_tie"
  specification$sample_size <- sample_size
  specification$local_kappa <- kappa
  specification$local_sign <- local_sign
  specification$tie_stratum_probability <- tie_stratum_probability
  specification$local_probability_difference <-
    local_probability_difference
  specification$tie_treatment_probability_z0 <- pi0
  specification$tie_treatment_probability_z1 <- pi1
  specification$local_gap_definition <- paste0(
    "pi1-pi0 = local_sign*kappa/sqrt(n*p_tie_stratum)"
  )
  specification$purpose <- paste0(
    "root-n local-to-envelope-tie binary stress experiment with no learned ",
    "threshold q"
  )
  specification
}

binary_local_envelope_tie_population_audit <- function(specification) {
  if (!is.list(specification) ||
      !identical(specification$name, "binary_local_to_envelope_tie")) {
    stop("Invalid local-to-envelope-tie specification.")
  }
  nuisance <- binary_envelope_tie_nuisance_table(specification)
  tie <- nuisance$tie_stratum
  aggregate_first_stage <- sum(
    nuisance$probability * (nuisance$pi1 - nuisance$pi0)
  )
  tie_first_stage <- sum(
    nuisance$probability[tie] *
      (nuisance$pi1[tie] - nuisance$pi0[tie])
  )
  data.frame(
    sample_size = specification$sample_size,
    kappa = specification$local_kappa,
    tie_strata = specification$tie_strata,
    tie_mass = specification$tie_mass,
    tie_stratum_probability = specification$tie_stratum_probability,
    local_probability_difference =
      specification$local_probability_difference,
    root_n_probability_difference = sqrt(specification$sample_size) *
      specification$local_probability_difference,
    tie_first_stage = tie_first_stage,
    aggregate_first_stage = aggregate_first_stage,
    minimum_propensity = min(nuisance$pi0, nuisance$pi1),
    maximum_propensity = max(nuisance$pi0, nuisance$pi1),
    stringsAsFactors = FALSE
  )
}
