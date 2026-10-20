# Fixed observed-data binary DGP with positive-mass envelope branch ties.
#
# The tie strata satisfy, for each treatment arm,
#   r_0(x) = r_1(x),  a_0(x) = a_1(x),  gamma_t(y,x) = 1.
# Consequently, both branches defining each hard Dual-IV--MSM envelope agree:
#   L_t(y,x) = max(1 / Lambda, 1 / Lambda),
#   U_t(y,x) = min(Lambda, Lambda).
# The Bernoulli success probabilities are kept away from 1/(Lambda+1) and
# Lambda/(Lambda+1), so the subsequent binary theta max/min operations have
# unique active branches.  This isolates the envelope kink without invoking a
# conditional threshold q or a second binary-outcome tie.

if (!exists("binary_arm_endpoint", mode = "function")) {
  if (!file.exists("R/binary_outcome_bounds.R")) {
    stop("Run this file from the repository root.", call. = FALSE)
  }
  source("R/binary_outcome_bounds.R")
}

make_binary_envelope_tie_spec <- function(
    Lambda = 4,
    tie_strata = 8L,
    tie_mass = 0.85,
    probability_z1 = 0.5,
    tie_treatment_probability = 0.5,
    regular_treatment_probabilities = c(z0 = 0.2, z1 = 0.8),
    tie_outcome_t0 = NULL,
    tie_outcome_t1 = NULL,
    regular_outcome_probabilities = c(t0 = 0.375, t1 = 0.625)
) {
  Lambda <- as.numeric(Lambda)
  tie_strata <- as.integer(tie_strata)
  if (length(Lambda) != 1L || !is.finite(Lambda) || Lambda <= 1) {
    stop("Lambda must be one finite number greater than one.")
  }
  if (length(tie_strata) != 1L || is.na(tie_strata) || tie_strata < 2L) {
    stop("tie_strata must be one integer of at least two.")
  }
  scalar_probability <- function(value, name) {
    if (length(value) != 1L || !is.finite(value) || value <= 0 || value >= 1) {
      stop(name, " must lie strictly between zero and one.")
    }
    as.numeric(value)
  }
  tie_mass <- scalar_probability(tie_mass, "tie_mass")
  probability_z1 <- scalar_probability(probability_z1, "probability_z1")
  if (length(tie_treatment_probability) == 1L) {
    tie_treatment_probability <- rep(
      tie_treatment_probability, tie_strata
    )
  }
  if (length(tie_treatment_probability) != tie_strata ||
      any(!is.finite(tie_treatment_probability)) ||
      any(tie_treatment_probability <= 0) ||
      any(tie_treatment_probability >= 1)) {
    stop(
      paste0(
        "tie_treatment_probability must be one interior probability or ",
        "one interior probability per tie stratum."
      )
    )
  }
  tie_treatment_probability <- as.numeric(tie_treatment_probability)
  if (length(regular_treatment_probabilities) != 2L ||
      any(!is.finite(regular_treatment_probabilities)) ||
      any(regular_treatment_probabilities <= 0) ||
      any(regular_treatment_probabilities >= 1)) {
    stop("regular_treatment_probabilities must contain two interior values.")
  }
  regular_treatment_probabilities <- stats::setNames(
    as.numeric(regular_treatment_probabilities), c("z0", "z1")
  )
  if (regular_treatment_probabilities[["z0"]] ==
      regular_treatment_probabilities[["z1"]]) {
    stop("The regular stratum must have a nonzero conditional first stage.")
  }

  if (is.null(tie_outcome_t1)) {
    tie_outcome_t1 <- seq(0.55, 0.70, length.out = tie_strata)
  }
  if (is.null(tie_outcome_t0)) {
    tie_outcome_t0 <- seq(0.45, 0.30, length.out = tie_strata)
  }
  if (length(tie_outcome_t0) != tie_strata ||
      length(tie_outcome_t1) != tie_strata ||
      any(!is.finite(c(tie_outcome_t0, tie_outcome_t1))) ||
      any(c(tie_outcome_t0, tie_outcome_t1) <= 0) ||
      any(c(tie_outcome_t0, tie_outcome_t1) >= 1)) {
    stop("Tie-stratum outcome probabilities are invalid.")
  }
  if (length(regular_outcome_probabilities) != 2L ||
      any(!is.finite(regular_outcome_probabilities)) ||
      any(regular_outcome_probabilities <= 0) ||
      any(regular_outcome_probabilities >= 1)) {
    stop("regular_outcome_probabilities must contain two interior values.")
  }
  regular_outcome_probabilities <- stats::setNames(
    as.numeric(regular_outcome_probabilities), c("t0", "t1")
  )

  # These are exactly the success probabilities at which the final binary
  # theta lower and upper branches tie.  Staying strictly between them makes
  # the final branches unique while retaining the earlier envelope ties.
  theta_upper_tie <- 1 / (Lambda + 1)
  theta_lower_tie <- Lambda / (Lambda + 1)
  all_outcome_probabilities <- c(
    tie_outcome_t0, tie_outcome_t1, regular_outcome_probabilities
  )
  if (any(all_outcome_probabilities <= theta_upper_tie) ||
      any(all_outcome_probabilities >= theta_lower_tie)) {
    stop(
      paste0(
        "Outcome probabilities must lie strictly between 1/(Lambda+1) and ",
        "Lambda/(Lambda+1) to exclude final binary theta ties."
      )
    )
  }

  total_strata <- tie_strata + 1L
  x_cols <- paste0("X", 2:total_strata)
  list(
    name = "binary_positive_mass_envelope_tie",
    Lambda = Lambda,
    tie_strata = tie_strata,
    total_strata = total_strata,
    tie_mass = tie_mass,
    stratum_probabilities = c(
      rep(tie_mass / tie_strata, tie_strata), 1 - tie_mass
    ),
    probability_z1 = probability_z1,
    tie_treatment_probability = tie_treatment_probability,
    regular_treatment_probabilities = regular_treatment_probabilities,
    outcome_t0 = c(
      as.numeric(tie_outcome_t0),
      regular_outcome_probabilities[["t0"]]
    ),
    outcome_t1 = c(
      as.numeric(tie_outcome_t1),
      regular_outcome_probabilities[["t1"]]
    ),
    theta_upper_tie = theta_upper_tie,
    theta_lower_tie = theta_lower_tie,
    x_cols = x_cols,
    purpose = paste0(
      "fixed-law binary envelope-kink experiment with no learned threshold q"
    )
  )
}

binary_envelope_tie_nuisance_table <- function(specification) {
  if (!is.list(specification) ||
      !specification$name %in% c(
        "binary_positive_mass_envelope_tie",
        "binary_local_to_envelope_tie"
      )) {
    stop("Invalid binary envelope-tie specification.")
  }
  stratum <- seq_len(specification$total_strata)
  tie <- stratum <= specification$tie_strata
  tie_pi0 <- if (!is.null(specification$tie_treatment_probability_z0)) {
    specification$tie_treatment_probability_z0
  } else {
    specification$tie_treatment_probability
  }
  tie_pi1 <- if (!is.null(specification$tie_treatment_probability_z1)) {
    specification$tie_treatment_probability_z1
  } else {
    specification$tie_treatment_probability
  }
  if (length(tie_pi0) != specification$tie_strata ||
      length(tie_pi1) != specification$tie_strata ||
      any(!is.finite(c(tie_pi0, tie_pi1))) ||
      any(c(tie_pi0, tie_pi1) <= 0) ||
      any(c(tie_pi0, tie_pi1) >= 1)) {
    stop("The tie-stratum treatment propensities are invalid.")
  }
  pi0 <- c(
    tie_pi0,
    specification$regular_treatment_probabilities[["z0"]]
  )
  pi1 <- c(
    tie_pi1,
    specification$regular_treatment_probabilities[["z1"]]
  )
  data.frame(
    stratum = stratum,
    probability = specification$stratum_probabilities,
    tie_stratum = tie,
    e = 1 - specification$probability_z1,
    pi0 = pi0,
    pi1 = pi1,
    outcome_t0_z0 = specification$outcome_t0,
    outcome_t0_z1 = specification$outcome_t0,
    outcome_t1_z0 = specification$outcome_t1,
    outcome_t1_z1 = specification$outcome_t1,
    stringsAsFactors = FALSE
  )
}

binary_envelope_tie_design_matrix <- function(stratum, specification) {
  stratum <- as.integer(stratum)
  if (anyNA(stratum) || any(stratum < 1L) ||
      any(stratum > specification$total_strata)) {
    stop("Stratum indices are invalid.")
  }
  matrix_value <- vapply(
    2:specification$total_strata,
    function(value) as.numeric(stratum == value),
    numeric(length(stratum))
  )
  if (is.null(dim(matrix_value))) {
    matrix_value <- matrix(matrix_value, ncol = 1L)
  }
  colnames(matrix_value) <- specification$x_cols
  matrix_value
}

generate_binary_envelope_tie_data <- function(
    n,
    specification = make_binary_envelope_tie_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  stratum <- sample.int(
    specification$total_strata,
    size = n,
    replace = TRUE,
    prob = specification$stratum_probabilities
  )
  nuisance <- binary_envelope_tie_nuisance_table(specification)
  z <- stats::rbinom(n, 1L, specification$probability_z1)
  pi0 <- nuisance$pi0[stratum]
  pi1 <- nuisance$pi1[stratum]
  probability_t <- ifelse(z == 1L, pi1, pi0)
  treatment <- stats::rbinom(n, 1L, probability_t)
  outcome_t0 <- specification$outcome_t0[stratum]
  outcome_t1 <- specification$outcome_t1[stratum]
  probability_y <- ifelse(treatment == 1L, outcome_t1, outcome_t0)
  y <- stats::rbinom(n, 1L, probability_y)
  output <- data.frame(
    Y = y,
    T_obs = treatment,
    Z = z,
    stratum = stratum,
    tie_stratum = stratum <= specification$tie_strata,
    stringsAsFactors = FALSE
  )
  output <- cbind(
    output,
    as.data.frame(binary_envelope_tie_design_matrix(stratum, specification))
  )
  output$e_true <- 1 - specification$probability_z1
  output$pi0_true <- pi0
  output$pi1_true <- pi1
  output$a_t0_z0_true <- outcome_t0
  output$a_t0_z1_true <- outcome_t0
  output$a_t1_z0_true <- outcome_t1
  output$a_t1_z1_true <- outcome_t1
  output$gamma0_oracle <- 1
  output$gamma1_oracle <- 1
  output$prob_Z <- specification$probability_z1
  output$prob_T <- probability_t
  output$prob_Y <- probability_y
  attr(output, "binary_envelope_tie_specification") <- specification
  output
}

compute_binary_envelope_tie_population <- function(
    specification = make_binary_envelope_tie_spec(),
    epsilon_grid = c(0, 0.01)
) {
  epsilon_grid <- sort(unique(as.numeric(epsilon_grid)))
  if (!length(epsilon_grid) || any(!is.finite(epsilon_grid)) ||
      any(epsilon_grid < 0) ||
      any(epsilon_grid >= 1 / specification$Lambda)) {
    stop("epsilon_grid contains invalid smoothing values.")
  }
  nuisance <- binary_envelope_tie_nuisance_table(specification)
  weights <- nuisance$probability
  Lambda <- specification$Lambda
  endpoint_names <- c(
    "mu1_lower", "mu1_upper", "mu0_lower", "mu0_upper",
    "ate_lower", "ate_upper"
  )
  rows <- vector("list", length(epsilon_grid))
  arm_storage <- vector("list", length(epsilon_grid))
  for (epsilon_index in seq_along(epsilon_grid)) {
    epsilon <- epsilon_grid[epsilon_index]
    arms <- list()
    for (treatment in 0:1) {
      arms[[as.character(treatment)]] <- binary_arm_endpoint(
        nuisance = list(
          r0 = if (treatment == 1L) nuisance$pi0 else 1 - nuisance$pi0,
          r1 = if (treatment == 1L) nuisance$pi1 else 1 - nuisance$pi1,
          a0 = if (treatment == 1L) {
            nuisance$outcome_t1_z0
          } else {
            nuisance$outcome_t0_z0
          },
          a1 = if (treatment == 1L) {
            nuisance$outcome_t1_z1
          } else {
            nuisance$outcome_t0_z1
          }
        ),
        Lambda0 = Lambda,
        Lambda1 = Lambda,
        epsilon = epsilon,
        compute_gradient = FALSE
      )
    }
    weighted <- function(value) sum(weights * value)
    values <- c(
      mu1_lower = weighted(arms[["1"]]$lower),
      mu1_upper = weighted(arms[["1"]]$upper),
      mu0_lower = weighted(arms[["0"]]$lower),
      mu0_upper = weighted(arms[["0"]]$upper)
    )
    values <- c(
      values,
      ate_lower = values[["mu1_lower"]] - values[["mu0_upper"]],
      ate_upper = values[["mu1_upper"]] - values[["mu0_lower"]]
    )
    rows[[epsilon_index]] <- data.frame(
      Lambda = Lambda,
      epsilon = epsilon,
      as.list(values[endpoint_names]),
      ate_width = values[["ate_upper"]] - values[["ate_lower"]],
      stringsAsFactors = FALSE
    )
    arm_storage[[epsilon_index]] <- arms
  }
  bounds <- do.call(rbind, rows)

  hard_arms <- arm_storage[[which.min(abs(epsilon_grid))]]
  audit_rows <- list()
  for (treatment in 0:1) {
    arm <- hard_arms[[as.character(treatment)]]
    r0 <- if (treatment == 1L) nuisance$pi0 else 1 - nuisance$pi0
    r1 <- if (treatment == 1L) nuisance$pi1 else 1 - nuisance$pi1
    a0 <- if (treatment == 1L) {
      nuisance$outcome_t1_z0
    } else {
      nuisance$outcome_t0_z0
    }
    a1 <- if (treatment == 1L) {
      nuisance$outcome_t1_z1
    } else {
      nuisance$outcome_t0_z1
    }
    gamma_one <- a1 / a0
    gamma_zero <- (1 - a1) / (1 - a0)
    candidate_lower <- function(gamma) {
      (gamma * (r1 + (1 - r1) / Lambda) - r0) / (1 - r0)
    }
    candidate_upper <- function(gamma) {
      (gamma * (r1 + (1 - r1) * Lambda) - r0) / (1 - r0)
    }
    lower_first <- a0 * arm$L_one
    lower_second <- 1 - (1 - a0) * arm$U_zero
    upper_first <- a0 * arm$U_one
    upper_second <- 1 - (1 - a0) * arm$L_zero
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      stratum = nuisance$stratum,
      probability = nuisance$probability,
      tie_stratum = nuisance$tie_stratum,
      treatment = treatment,
      lower_one_signed_gap = candidate_lower(gamma_one) - 1 / Lambda,
      upper_one_signed_gap = candidate_upper(gamma_one) - Lambda,
      lower_zero_signed_gap = candidate_lower(gamma_zero) - 1 / Lambda,
      upper_zero_signed_gap = candidate_upper(gamma_zero) - Lambda,
      lower_one_branch_gap = abs(candidate_lower(gamma_one) - 1 / Lambda),
      upper_one_branch_gap = abs(candidate_upper(gamma_one) - Lambda),
      lower_zero_branch_gap = abs(candidate_lower(gamma_zero) - 1 / Lambda),
      upper_zero_branch_gap = abs(candidate_upper(gamma_zero) - Lambda),
      theta_lower_branch_gap = abs(lower_first - lower_second),
      theta_upper_branch_gap = abs(upper_first - upper_second),
      gamma_violation = arm$gamma_violation,
      envelope_crossing = arm$envelope_crossing,
      theta_crossing = arm$theta_crossing,
      stringsAsFactors = FALSE
    )
  }
  audit <- do.call(rbind, audit_rows)
  list(
    bounds = bounds,
    audit = audit,
    nuisance = nuisance,
    arms = stats::setNames(arm_storage, as.character(epsilon_grid)),
    compatible = !any(
      audit$gamma_violation | audit$envelope_crossing | audit$theta_crossing
    )
  )
}
