# Full observed-data IV embedding of the local lower-envelope kink.
#
# The first covariate stratum carries the local nonregularity. Its treatment
# probability under Z = 0 is small but positive, and the Z = 1 increment is
# calibrated so that, for t = 1 and gamma_1 = 1,
#
#   candidate_lower - 1 / Lambda = local_gap_coefficient / sqrt(n).
#
# The remaining strata provide a regular first stage. Outcomes are continuous
# truncated normals, so threshold optimizers are unique. A mild density tilt
# in the final regular stratum keeps oracle gamma nonconstant for diagnostics.

if (!exists("make_envelope_tie_stress_spec", mode = "function")) {
  if (!file.exists("R/envelope_tie_stress_dgp.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/envelope_tie_stress_dgp.R")
}

make_full_iv_envelope_kink_spec <- function(
    sample_size,
    sensitivity_lambda = 3,
    local_gap_coefficient = 0.232379000772445,
    stress_treatment_probability = 0.05,
    x_probabilities = c(0.9, 0.05, 0.05),
    regular_pi0 = 0.3,
    regular_pi1 = 0.7
) {
  sample_size <- as.integer(sample_size)
  values <- c(
    sensitivity_lambda, local_gap_coefficient,
    stress_treatment_probability, x_probabilities,
    regular_pi0, regular_pi1
  )
  if (length(sample_size) != 1L || is.na(sample_size) || sample_size < 100L ||
      any(!is.finite(values)) || sensitivity_lambda <= 1 ||
      local_gap_coefficient <= 0 ||
      stress_treatment_probability <= 0 ||
      stress_treatment_probability >= 0.5 ||
      length(x_probabilities) != 3L || any(x_probabilities <= 0) ||
      abs(sum(x_probabilities) - 1) > 1e-12 ||
      regular_pi0 <= 0 || regular_pi1 >= 1 ||
      regular_pi0 >= regular_pi1) {
    stop("Invalid full-IV envelope-kink specification.")
  }

  lower_branch_slope <- (1 - 1 / sensitivity_lambda) /
    (1 - stress_treatment_probability)
  local_first_stage <- local_gap_coefficient /
    (sqrt(sample_size) * lower_branch_slope)
  if (stress_treatment_probability + local_first_stage >= 0.5) {
    stop("The local first stage is too large for the stress stratum.")
  }

  specification <- list(
    name = "full_iv_local_envelope_kink",
    dimension = 2L,
    sensitivity_lambda = sensitivity_lambda,
    x_values = 0:2,
    x_probabilities = as.numeric(x_probabilities),
    probability_z1 = 0.5,
    pi0 = c(
      stress_treatment_probability,
      regular_pi0,
      regular_pi0
    ),
    pi1 = c(
      stress_treatment_probability + local_first_stage,
      regular_pi1,
      regular_pi1
    ),
    outcome_treatment_shift = 0.75,
    outcome_x_shift = c(-0.3, 0, 0.3),
    comparison_tilt_t0 = -0.075,
    comparison_tilt_t1 = 0.075,
    residual_lower = -2,
    residual_upper = 2,
    tie_strata = c(FALSE, FALSE, FALSE),
    nominal_sample_size = sample_size,
    local_gap_coefficient = local_gap_coefficient,
    local_first_stage = local_first_stage,
    stress_treatment_probability = stress_treatment_probability,
    population_tie_probability = 0,
    population_local_to_tie_probability = x_probabilities[1L],
    marginal_first_stage = sum(
      x_probabilities * c(
        local_first_stage,
        regular_pi1 - regular_pi0,
        regular_pi1 - regular_pi0
      )
    )
  )
  class(specification) <- c(
    "full_iv_envelope_kink_specification", "list"
  )
  specification
}

full_iv_envelope_kink_branch_diagnostics <- function(specification) {
  if (!exists("compute_eif_envelope_derivatives", mode = "function")) {
    stop("Source R/sharp_bounds_plugin.R before computing branch diagnostics.")
  }
  nuisance <- envelope_tie_stress_nuisances(0, specification)
  hard <- compute_eif_envelope_derivatives(
    nuisance = nuisance,
    gamma = 1,
    t_value = 1,
    Lambda0 = specification$sensitivity_lambda,
    Lambda1 = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = 0
  )
  smooth <- compute_eif_envelope_derivatives(
    nuisance = nuisance,
    gamma = 1,
    t_value = 1,
    Lambda0 = specification$sensitivity_lambda,
    Lambda1 = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = 0.01
  )
  data.frame(
    sample_size = specification$nominal_sample_size,
    local_first_stage = specification$local_first_stage,
    requested_lower_gap = specification$local_gap_coefficient /
      sqrt(specification$nominal_sample_size),
    actual_lower_gap = hard$candidate_lower -
      1 / specification$sensitivity_lambda,
    actual_upper_gap = hard$candidate_upper -
      specification$sensitivity_lambda,
    hard_lower_weight = hard$lower_derivative_weight,
    smooth_lower_weight = smooth$lower_derivative_weight,
    hard_upper_weight = hard$upper_derivative_weight,
    smooth_upper_weight = smooth$upper_derivative_weight,
    stringsAsFactors = FALSE
  )
}
