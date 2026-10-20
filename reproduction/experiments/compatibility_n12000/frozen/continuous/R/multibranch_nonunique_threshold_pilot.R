# Oracle pilot for score-distinct, nonunique sharp-bound thresholds.
#
# This file constructs a finite library of threshold functions for one arm of
# the lower sharp-bound endpoint.  Every library member minimizes the same
# population weighted check loss.  The alternative thresholds lie strictly
# inside the same outcome-support gap, so they induce the same normalized
# direct extremizer at the population law.  Their pathwise scores nevertheless
# differ because perturbations of the envelope enter through (Y - q).  This
# isolates genuine nonuniqueness of the original direct sharp-bound functional,
# rather than a boundary-atom normalization residual in its check-loss
# representation.
#
# The construction isolates the nonregular selection problem considered by
# Luedtke and van der Laan: nuisance functions and candidate-specific score
# variances are known exactly.  It compares
#
#   1. ordinary K-fold cross-fitting, which selects a threshold on the
#      overlapping training folds and uses the pooled empirical standard
#      error; and
#   2. predictable online selection, which selects only from past data and
#      uses the exact candidate-specific predictable variance.
#
# It is a focused gate before integrating a new law into the expensive RKHS
# simulation driver.  It does not replace the observed-data simulation DGP.

if (!exists("compute_eif_row_primitives", mode = "function")) {
  if (!file.exists("R/sharp_bounds_plugin.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/sharp_bounds_plugin.R")
}

multibranch_hadamard <- function(order) {
  order <- as.integer(order)
  if (length(order) != 1L || is.na(order) || order < 2L ||
      bitwAnd(order, order - 1L) != 0L) {
    stop("order must be a power of two of at least two.")
  }
  output <- matrix(1, nrow = 1L, ncol = 1L)
  while (nrow(output) < order) {
    output <- rbind(
      cbind(output, output),
      cbind(output, -output)
    )
  }
  output
}

make_multibranch_nonunique_threshold_spec <- function(
    sensitivity_lambda = 3,
    strata = 16L,
    candidates = 16L,
    candidate_design = c("hadamard", "complete_binary"),
    endpoint_direction = c("lower", "upper"),
    treatment = 1L,
    probability_z0 = 0.5,
    pi0 = 0.3,
    pi1 = 0.7,
    lower_atom = 0,
    boundary_atom = 1,
    tail_atom = 8,
    tail_probability = 0.05,
    threshold_gap_fractions = c(0.25, 0.75),
    gamma_atoms = c(1, 1, 1),
    atoms = NULL,
    gap_index = NULL,
    probability_null_fraction = 0.5,
    atom_probability_override = NULL,
    atom_shift_by_stratum = NULL,
    envelope_smoothing_epsilon = 0
) {
  candidate_design <- match.arg(candidate_design)
  endpoint_direction <- match.arg(endpoint_direction)
  sensitivity_lambda <- as.numeric(sensitivity_lambda)
  strata <- as.integer(strata)
  candidates <- as.integer(candidates)
  treatment <- as.integer(treatment)
  if (is.null(atoms)) {
    atoms <- c(lower_atom, boundary_atom, tail_atom)
  } else {
    atoms <- as.numeric(atoms)
  }
  if (is.null(gap_index)) {
    gap_index <- if (endpoint_direction == "lower") 1L else length(atoms) - 1L
  }
  gap_index <- as.integer(gap_index)
  probability_null_fraction <- as.numeric(probability_null_fraction)
  if (!is.null(atom_probability_override)) {
    atom_probability_override <- as.numeric(atom_probability_override)
  }
  if (is.null(atom_shift_by_stratum)) {
    atom_shift_by_stratum <- rep(0, strata)
  } else {
    atom_shift_by_stratum <- as.numeric(atom_shift_by_stratum)
  }
  numeric_values <- c(
    sensitivity_lambda, probability_z0, pi0, pi1,
    atoms, tail_probability, probability_null_fraction,
    envelope_smoothing_epsilon
  )
  gamma_atoms <- as.numeric(gamma_atoms)
  threshold_gap_fractions <- as.numeric(threshold_gap_fractions)
  if (
    any(!is.finite(numeric_values)) ||
      sensitivity_lambda <= 1 ||
      strata < 2L ||
      candidates < 2L ||
      !treatment %in% 0:1 ||
      probability_z0 <= 0 || probability_z0 >= 1 ||
      pi0 <= 0 || pi0 >= 1 || pi1 <= 0 || pi1 >= 1 ||
      length(atoms) < 3L || any(diff(atoms) <= 0) ||
      length(gap_index) != 1L || is.na(gap_index) ||
      gap_index < 1L || gap_index >= length(atoms) ||
      tail_probability <= 0 || tail_probability >= 1 ||
      length(gamma_atoms) != length(atoms) ||
      any(!is.finite(gamma_atoms)) ||
      any(gamma_atoms <= 0) ||
      length(threshold_gap_fractions) != 2L ||
      any(!is.finite(threshold_gap_fractions)) ||
      any(threshold_gap_fractions <= 0) ||
      any(threshold_gap_fractions >= 1) ||
      threshold_gap_fractions[1L] >= threshold_gap_fractions[2L] ||
      probability_null_fraction <= 0 || probability_null_fraction >= 1 ||
      (!is.null(atom_probability_override) && (
        length(atom_probability_override) != length(atoms) ||
          any(!is.finite(atom_probability_override)) ||
          any(atom_probability_override <= 0) ||
          abs(sum(atom_probability_override) - 1) > 1e-10
      )) ||
      length(atom_shift_by_stratum) != strata ||
      any(!is.finite(atom_shift_by_stratum)) ||
      envelope_smoothing_epsilon < 0
  ) {
    stop("Invalid multi-branch nonunique-threshold specification.")
  }
  if (candidate_design == "hadamard") {
    if (bitwAnd(strata, strata - 1L) != 0L || candidates > strata) {
      stop(
        "The Hadamard design requires power-of-two strata and candidates <= strata."
      )
    }
    code <- multibranch_hadamard(strata)[
      seq_len(candidates), , drop = FALSE
    ] > 0
  } else {
    if (strata > 10L) {
      stop("The complete-binary design is limited to at most ten strata.")
    }
    code <- as.matrix(expand.grid(
      rep(list(c(FALSE, TRUE)), strata),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ))
    candidates <- nrow(code)
  }
  list(
    sensitivity_lambda = sensitivity_lambda,
    strata = strata,
    candidates = candidates,
    candidate_design = candidate_design,
    endpoint_direction = endpoint_direction,
    treatment = treatment,
    probability_z0 = probability_z0,
    pi0 = pi0,
    pi1 = pi1,
    atoms = atoms,
    gap_index = gap_index,
    probability_null_fraction = probability_null_fraction,
    atom_probability_override = atom_probability_override,
    atom_shift_by_stratum = atom_shift_by_stratum,
    tail_probability = tail_probability,
    threshold_gap_fractions = threshold_gap_fractions,
    gamma_atoms = gamma_atoms,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    # TRUE selects the right endpoint of the lower-quantile optimizer
    # interval; FALSE selects the left endpoint.
    candidate_right_endpoint = code
  )
}

multibranch_weighted_group_mean <- function(
    value,
    probability,
    group,
    levels
) {
  output <- numeric(length(levels))
  for (index in seq_along(levels)) {
    local <- group == levels[index]
    mass <- sum(probability[local])
    if (!is.finite(mass) || mass <= 0) {
      stop("An oracle conditional-mean group has zero probability.")
    }
    output[index] <- sum(probability[local] * value[local]) / mass
  }
  output
}

build_multibranch_nonunique_threshold_oracle <- function(
    specification = make_multibranch_nonunique_threshold_spec()
) {
  required <- c(
    "sensitivity_lambda", "strata", "candidates", "treatment",
    "endpoint_direction",
    "probability_z0", "pi0", "pi1", "atoms", "gap_index",
    "probability_null_fraction", "tail_probability",
    "threshold_gap_fractions",
    "gamma_atoms",
    "envelope_smoothing_epsilon", "candidate_right_endpoint"
  )
  if (!all(required %in% names(specification))) {
    stop("The multi-branch specification is incomplete.")
  }
  Lambda <- specification$sensitivity_lambda
  treatment <- specification$treatment
  x_levels <- seq_len(specification$strata)
  nuisance_atoms <- list(
    e = rep(specification$probability_z0, length(specification$atoms)),
    pi0 = rep(specification$pi0, length(specification$atoms)),
    pi1 = rep(specification$pi1, length(specification$atoms))
  )
  envelope_atoms <- compute_projected_envelopes(
    nuisance = nuisance_atoms,
    gamma = specification$gamma_atoms,
    t_value = treatment,
    Lambda0 = Lambda,
    Lambda1 = Lambda,
    envelope_smoothing_epsilon =
      specification$envelope_smoothing_epsilon,
    gamma_compatibility_action = "stop"
  )
  if (
    any(envelope_atoms$gamma_incompatible_original) ||
      any(!is.finite(envelope_atoms$L)) ||
      any(!is.finite(envelope_atoms$U)) ||
      any(envelope_atoms$L > envelope_atoms$U + 1e-12)
  ) {
    stop(
      paste0(
        "The all-oracle gamma values are outside the compatible interval ",
        "or produce crossed envelopes."
      )
    )
  }
  atom_probability_override <- specification$atom_probability_override
  if (!is.null(atom_probability_override)) {
    atom_probability <- atom_probability_override
    lower_group <- seq_along(specification$atoms) <=
      specification$gap_index
    normalizing_lambda <- if (specification$endpoint_direction == "lower") {
      ifelse(lower_group, envelope_atoms$U, envelope_atoms$L)
    } else {
      ifelse(lower_group, envelope_atoms$L, envelope_atoms$U)
    }
    if (
      abs(sum(atom_probability * specification$gamma_atoms) - 1) >
        1e-10 ||
        abs(sum(atom_probability * normalizing_lambda) - 1) > 1e-10
    ) {
      stop(
        "The supplied atom probabilities violate a population normalization."
      )
    }
  } else {
    constant_gamma <- max(abs(specification$gamma_atoms - 1)) < 1e-12
    if (constant_gamma) {
    provisional_tau <- (envelope_atoms$U[1L] - 1) /
      (envelope_atoms$U[1L] - envelope_atoms$L[1L])
    provisional_lower_level <- 1 - provisional_tau
    if (length(specification$atoms) == 3L &&
        specification$gap_index == if (
          specification$endpoint_direction == "lower"
        ) 1L else 2L) {
      atom_probability <- if (specification$endpoint_direction == "lower") {
        c(
          provisional_lower_level,
          1 - provisional_lower_level - specification$tail_probability,
          specification$tail_probability
        )
      } else {
        c(
          specification$tail_probability,
          provisional_tau - specification$tail_probability,
          1 - provisional_tau
        )
      }
    } else {
      lower_mass <- if (
        specification$endpoint_direction == "lower"
      ) provisional_lower_level else provisional_tau
      atom_probability <- c(
        rep(lower_mass / specification$gap_index,
            specification$gap_index),
        rep(
          (1 - lower_mass) /
            (length(specification$atoms) - specification$gap_index),
          length(specification$atoms) - specification$gap_index
        )
      )
    }
    } else {
      lower_group <- seq_along(specification$atoms) <=
        specification$gap_index
      normalizing_lambda <- if (specification$endpoint_direction == "lower") {
        ifelse(lower_group, envelope_atoms$U, envelope_atoms$L)
      } else {
        ifelse(lower_group, envelope_atoms$L, envelope_atoms$U)
      }
      equation <- rbind(
        probability = rep(1, length(specification$atoms)),
        density_ratio = specification$gamma_atoms,
        extremizer_normalization = normalizing_lambda
      )
      atom_probability <- tryCatch({
      target <- rep(1, 3L)
      if (ncol(equation) == 3L) {
        as.vector(solve(equation, target))
      } else if (ncol(equation) == 4L) {
        decomposition <- svd(equation, nu = 0L, nv = 4L)
        rank <- sum(decomposition$d > max(decomposition$d) * 1e-10)
        if (rank != 3L) stop("The probability constraints lack full rank.")
        particular <- as.vector(
          t(equation) %*% solve(equation %*% t(equation), target)
        )
        null_direction <- decomposition$v[, 4L]
        lower_bound <- -Inf
        upper_bound <- Inf
        for (index in seq_along(particular)) {
          if (null_direction[index] > 0) {
            lower_bound <- max(
              lower_bound, -particular[index] / null_direction[index]
            )
          } else if (null_direction[index] < 0) {
            upper_bound <- min(
              upper_bound, -particular[index] / null_direction[index]
            )
          } else if (particular[index] <= 0) {
            stop("No positive atom-probability solution exists.")
          }
        }
        if (!is.finite(lower_bound) || !is.finite(upper_bound) ||
            lower_bound >= upper_bound) {
          stop("No positive atom-probability interval exists.")
        }
        fraction <- specification$probability_null_fraction
        particular + (
          lower_bound + fraction * (upper_bound - lower_bound)
        ) * null_direction
      } else {
        stop("Only three- and four-atom probability systems are supported.")
      }
      }, error = function(error) rep(NA_real_, length(specification$atoms)))
    }
  }
  if (any(!is.finite(atom_probability)) || any(atom_probability <= 0)) {
    stop(
      paste0(
        "The gamma/tail specification does not produce three positive ",
        "baseline atom probabilities."
      )
    )
  }
  Lbar <- sum(atom_probability * envelope_atoms$L)
  Ubar <- sum(atom_probability * envelope_atoms$U)
  tau <- (Ubar - 1) / (Ubar - Lbar)
  lower_level <- 1 - tau

  categories <- expand.grid(
    X = x_levels,
    Z = 0:1,
    T_obs = 0:1,
    atom = seq_along(specification$atoms),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  categories$Y <- specification$atoms[categories$atom] +
    specification$atom_shift_by_stratum[categories$X]
  probability_z <- ifelse(
    categories$Z == 0,
    specification$probability_z0,
    1 - specification$probability_z0
  )
  local_pi <- ifelse(
    categories$Z == 0,
    specification$pi0,
    specification$pi1
  )
  probability_t <- ifelse(
    categories$T_obs == 1,
    local_pi,
    1 - local_pi
  )
  categories$probability <-
    (1 / specification$strata) * probability_z * probability_t *
    atom_probability[categories$atom] * ifelse(
      categories$Z == 0,
      1,
      specification$gamma_atoms[categories$atom]
    )
  categories$probability <- categories$probability /
    sum(categories$probability)

  row_count <- nrow(categories)
  nuisance <- list(
    e = rep(specification$probability_z0, row_count),
    pi0 = rep(specification$pi0, row_count),
    pi1 = rep(specification$pi1, row_count)
  )
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = specification$gamma_atoms[categories$atom],
    t_value = treatment,
    Lambda0 = Lambda,
    Lambda1 = Lambda,
    envelope_smoothing_epsilon =
      specification$envelope_smoothing_epsilon,
    gamma_compatibility_action = "stop"
  )
  baseline <- categories$Z == 0 & categories$T_obs == treatment
  comparison <- categories$Z == 1 & categories$T_obs == treatment
  observed_arm <- as.numeric(categories$T_obs == treatment)
  instrument0 <- as.numeric(categories$Z == 0)
  instrument1 <- 1 - instrument0
  r0 <- if (treatment == 1L) specification$pi0 else 1 - specification$pi0
  r1 <- if (treatment == 1L) specification$pi1 else 1 - specification$pi1

  score <- loss <- matrix(
    NA_real_, nrow = row_count, ncol = specification$candidates
  )
  q_matrix <- matrix(
    NA_real_, nrow = row_count, ncol = specification$candidates
  )
  direct_target <- normalization_residual <- numeric(
    specification$candidates
  )
  for (candidate in seq_len(specification$candidates)) {
    right <- specification$candidate_right_endpoint[
      candidate, categories$X
    ]
    gap_left <- specification$atoms[specification$gap_index]
    gap_right <- specification$atoms[specification$gap_index + 1L]
    q <- specification$atom_shift_by_stratum[categories$X] +
      gap_left + (gap_right - gap_left) * ifelse(
      right,
      specification$threshold_gap_fractions[2L],
      specification$threshold_gap_fractions[1L]
    )
    q_matrix[, candidate] <- q
    primitive_all <- compute_eif_row_primitives(
      data = categories,
      nuisance = nuisance,
      envelope = envelope,
      t_value = treatment,
      Lambda0 = Lambda,
      Lambda1 = Lambda,
      q_lower = q,
      q_upper = q
    )
    primitive <- primitive_all[[specification$endpoint_direction]]

    m0_x <- multibranch_weighted_group_mean(
      categories$Y[baseline], categories$probability[baseline],
      categories$X[baseline], x_levels
    )
    Hbar0_x <- multibranch_weighted_group_mean(
      primitive$H0[baseline], categories$probability[baseline],
      categories$X[baseline], x_levels
    )
    Hbar1_x <- multibranch_weighted_group_mean(
      primitive$H1[comparison], categories$probability[comparison],
      categories$X[comparison], x_levels
    )
    D_r0_x <- multibranch_weighted_group_mean(
      primitive$D_r0[baseline], categories$probability[baseline],
      categories$X[baseline], x_levels
    )
    D_r1_x <- multibranch_weighted_group_mean(
      primitive$D_r1[baseline], categories$probability[baseline],
      categories$X[baseline], x_levels
    )
    m0 <- m0_x[categories$X]
    Hbar0 <- Hbar0_x[categories$X]
    Hbar1 <- Hbar1_x[categories$X]
    linked <- compute_fixed_q_linked_components(
      m0 = m0,
      Hbar0 = Hbar0,
      Hbar1 = Hbar1,
      q = q,
      r0 = rep(r0, row_count),
      mean_D_r0 = D_r0_x[categories$X]
    )
    C1 <- (1 - r0) * D_r1_x[categories$X]
    treatment_correction <-
      instrument0 / specification$probability_z0 * linked$C0 *
      (observed_arm - r0) +
      instrument1 / (1 - specification$probability_z0) * C1 *
      (observed_arm - r1)
    outcome_correction <-
      instrument0 / specification$probability_z0 * observed_arm / r0 *
      (primitive$H0 - Hbar0) +
      instrument1 / (1 - specification$probability_z0) * observed_arm / r1 *
      (primitive$H1 - Hbar1)
    score[, candidate] <- linked$M + treatment_correction +
      outcome_correction

    # Candidate selection uses exactly the baseline weighted check loss used
    # by the threshold learner.  The omitted positive constant is common to
    # all candidates and therefore cannot affect the selected index.
    loss[, candidate] <- ifelse(
      baseline,
      envelope$D * check_loss(
        categories$Y,
        q,
        rep(
          if (specification$endpoint_direction == "lower") {
            lower_level
          } else {
            tau
          },
          row_count
        )
      ),
      0
    )
    lambda <- if (specification$endpoint_direction == "lower") {
      ifelse(categories$Y <= q, envelope$U, envelope$L)
    } else {
      ifelse(categories$Y <= q, envelope$L, envelope$U)
    }
    lambda_x <- multibranch_weighted_group_mean(
      lambda[baseline], categories$probability[baseline],
      categories$X[baseline], x_levels
    )
    theta_x <- multibranch_weighted_group_mean(
      (lambda * categories$Y)[baseline],
      categories$probability[baseline],
      categories$X[baseline], x_levels
    )
    normalization_residual[candidate] <- max(abs(lambda_x - 1))
    direct_target[candidate] <- mean(r0 * m0_x + (1 - r0) * theta_x)
  }

  population_mean <- colSums(categories$probability * score)
  population_risk <- colSums(categories$probability * loss)
  centered <- sweep(score, 2L, population_mean, "-")
  population_covariance <- crossprod(
    sqrt(categories$probability) * centered
  )
  population_sd <- sqrt(diag(population_covariance))
  population_correlation <- population_covariance /
    tcrossprod(population_sd)
  if (
    max(population_mean) - min(population_mean) > 1e-10 ||
      max(population_risk) - min(population_risk) > 1e-10 ||
      max(direct_target) - min(direct_target) > 1e-10 ||
      max(normalization_residual) > 1e-10 ||
      max(abs(population_mean - direct_target)) > 1e-10 ||
      any(!is.finite(population_sd)) || any(population_sd <= 0)
  ) {
    stop("The oracle library is not exactly tied with positive variance.")
  }

  list(
    specification = specification,
    categories = categories,
    q = q_matrix,
    score = score,
    selection_loss = loss,
    target = mean(direct_target),
    direct_target = direct_target,
    direct_normalization_residual = normalization_residual,
    population_mean = population_mean,
    population_risk = population_risk,
    population_sd = population_sd,
    population_correlation = population_correlation,
    tau = tau,
    lower_level = lower_level,
    threshold_level = if (
      specification$endpoint_direction == "lower"
    ) lower_level else tau,
    atom_probability = atom_probability
  )
}

# Clean continuous-gap law used for the five-fold stress experiment.  Unlike
# the atomic construction above, the candidate thresholds are the two support
# endpoints of a genuine outcome gap.  Those endpoints have probability zero,
# so every candidate induces the same normalized direct extremizer at the
# population law while retaining a candidate-specific pathwise score.
make_fivefold_gap_nonunique_threshold_spec <- function(
    sensitivity_lambda = 3,
    strata = 4L,
    candidate_design = c("complete_binary", "hadamard"),
    endpoint_direction = c("lower", "upper"),
    treatment = 0L,
    probability_z0 = 0.5,
    pi0 = 0.7,
    pi1 = 0.9,
    lower_interval = c(-0.05, 0),
    upper_interval = c(1, 1.05),
    quadrature_nodes_per_interval = 24L,
    baseline_lower_probability = NULL,
    gamma_components = c(1, 1),
    envelope_smoothing_epsilon = 0
) {
  candidate_design <- match.arg(candidate_design)
  endpoint_direction <- match.arg(endpoint_direction)
  treatment <- as.integer(treatment)
  strata <- as.integer(strata)
  quadrature_nodes_per_interval <- as.integer(
    quadrature_nodes_per_interval
  )
  lower_interval <- as.numeric(lower_interval)
  upper_interval <- as.numeric(upper_interval)
  gamma_components <- as.numeric(gamma_components)
  if (!is.null(baseline_lower_probability)) {
    baseline_lower_probability <- as.numeric(baseline_lower_probability)
  }
  numeric_values <- c(
    sensitivity_lambda, probability_z0, pi0, pi1,
    lower_interval, upper_interval, gamma_components,
    envelope_smoothing_epsilon
  )
  if (
    any(!is.finite(numeric_values)) || sensitivity_lambda <= 1 ||
      strata < 2L || !treatment %in% 0:1 ||
      probability_z0 <= 0 || probability_z0 >= 1 ||
      pi0 <= 0 || pi0 >= 1 || pi1 <= 0 || pi1 >= 1 ||
      length(lower_interval) != 2L || length(upper_interval) != 2L ||
      lower_interval[1L] >= lower_interval[2L] ||
      upper_interval[1L] >= upper_interval[2L] ||
      lower_interval[2L] >= upper_interval[1L] ||
      quadrature_nodes_per_interval < 2L ||
      length(gamma_components) != 2L || any(gamma_components <= 0) ||
      (!is.null(baseline_lower_probability) && (
        length(baseline_lower_probability) != 1L ||
          !is.finite(baseline_lower_probability) ||
          baseline_lower_probability <= 0 ||
          baseline_lower_probability >= 1 ||
          abs(
            baseline_lower_probability * gamma_components[1L] +
              (1 - baseline_lower_probability) * gamma_components[2L] - 1
          ) > 1e-10
      )) ||
      envelope_smoothing_epsilon < 0
  ) {
    stop("Invalid five-fold continuous-gap specification.")
  }
  library_specification <- make_multibranch_nonunique_threshold_spec(
    sensitivity_lambda = sensitivity_lambda,
    strata = strata,
    candidates = strata,
    candidate_design = candidate_design,
    endpoint_direction = endpoint_direction,
    treatment = treatment,
    probability_z0 = probability_z0,
    pi0 = pi0,
    pi1 = pi1,
    lower_atom = lower_interval[1L],
    boundary_atom = lower_interval[2L],
    tail_atom = upper_interval[2L],
    tail_probability = 0.05,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon
  )
  list(
    name = "fivefold_continuous_gap_nonunique_threshold",
    sensitivity_lambda = as.numeric(sensitivity_lambda),
    strata = strata,
    candidates = library_specification$candidates,
    candidate_design = candidate_design,
    endpoint_direction = endpoint_direction,
    treatment = treatment,
    probability_z0 = as.numeric(probability_z0),
    pi0 = as.numeric(pi0),
    pi1 = as.numeric(pi1),
    lower_interval = lower_interval,
    upper_interval = upper_interval,
    quadrature_nodes_per_interval = quadrature_nodes_per_interval,
    baseline_lower_probability = baseline_lower_probability,
    gamma_components = gamma_components,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    candidate_right_endpoint =
      library_specification$candidate_right_endpoint
  )
}

fivefold_gap_midpoint_rule <- function(interval, nodes) {
  width <- diff(interval) / nodes
  list(
    nodes = interval[1L] + (seq_len(nodes) - 0.5) * width,
    weights = rep(1 / nodes, nodes)
  )
}

build_fivefold_gap_nonunique_threshold_oracle <- function(
    specification = make_fivefold_gap_nonunique_threshold_spec()
) {
  required <- c(
    "sensitivity_lambda", "strata", "candidates", "treatment",
    "endpoint_direction",
    "probability_z0", "pi0", "pi1", "lower_interval",
    "upper_interval", "quadrature_nodes_per_interval",
    "baseline_lower_probability", "gamma_components",
    "envelope_smoothing_epsilon", "candidate_right_endpoint"
  )
  if (!all(required %in% names(specification))) {
    stop("The five-fold continuous-gap specification is incomplete.")
  }
  Lambda <- specification$sensitivity_lambda
  treatment <- specification$treatment
  x_levels <- seq_len(specification$strata)
  component_nuisance <- list(
    e = rep(specification$probability_z0, 2L),
    pi0 = rep(specification$pi0, 2L),
    pi1 = rep(specification$pi1, 2L)
  )
  component_envelope <- compute_projected_envelopes(
    nuisance = component_nuisance,
    gamma = specification$gamma_components,
    t_value = treatment,
    Lambda0 = Lambda,
    Lambda1 = Lambda,
    envelope_smoothing_epsilon =
      specification$envelope_smoothing_epsilon,
    gamma_compatibility_action = "stop"
  )
  lower_component_probability <- specification$baseline_lower_probability
  if (is.null(lower_component_probability)) {
    if (max(abs(specification$gamma_components - 1)) > 1e-12) {
      stop(
        "A nonconstant gamma_components vector requires a fixed baseline probability."
      )
    }
    provisional_tau <- (component_envelope$U[1L] - 1) /
      (component_envelope$U[1L] - component_envelope$L[1L])
    lower_component_probability <- if (
      specification$endpoint_direction == "lower"
    ) 1 - provisional_tau else provisional_tau
  }
  component_probability <- c(
    lower_component_probability, 1 - lower_component_probability
  )
  Lbar <- sum(component_probability * component_envelope$L)
  Ubar <- sum(component_probability * component_envelope$U)
  tau <- (Ubar - 1) / (Ubar - Lbar)
  lower_level <- 1 - tau
  if (!is.finite(lower_level) || lower_level <= 0 || lower_level >= 1) {
    stop("The continuous-gap lower component has invalid probability.")
  }

  low_rule <- fivefold_gap_midpoint_rule(
    specification$lower_interval,
    specification$quadrature_nodes_per_interval
  )
  high_rule <- fivefold_gap_midpoint_rule(
    specification$upper_interval,
    specification$quadrature_nodes_per_interval
  )
  outcome_nodes <- c(low_rule$nodes, high_rule$nodes)
  outcome_probability <- c(
    lower_component_probability * low_rule$weights,
    (1 - lower_component_probability) * high_rule$weights
  )
  categories <- expand.grid(
    X = x_levels,
    Z = 0:1,
    T_obs = 0:1,
    outcome_index = seq_along(outcome_nodes),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  categories$Y <- outcome_nodes[categories$outcome_index]
  categories$outcome_component <- ifelse(
    categories$outcome_index <= length(low_rule$nodes), 1L, 2L
  )
  probability_z <- ifelse(
    categories$Z == 0,
    specification$probability_z0,
    1 - specification$probability_z0
  )
  local_pi <- ifelse(
    categories$Z == 0, specification$pi0, specification$pi1
  )
  probability_t <- ifelse(
    categories$T_obs == 1, local_pi, 1 - local_pi
  )
  outcome_ratio <- ifelse(
    categories$Z == 1L & categories$T_obs == treatment,
    specification$gamma_components[categories$outcome_component],
    1
  )
  categories$probability <-
    (1 / specification$strata) * probability_z * probability_t *
    outcome_probability[categories$outcome_index] * outcome_ratio

  row_count <- nrow(categories)
  nuisance <- list(
    e = rep(specification$probability_z0, row_count),
    pi0 = rep(specification$pi0, row_count),
    pi1 = rep(specification$pi1, row_count)
  )
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = specification$gamma_components[categories$outcome_component],
    t_value = treatment,
    Lambda0 = Lambda,
    Lambda1 = Lambda,
    envelope_smoothing_epsilon =
      specification$envelope_smoothing_epsilon,
    gamma_compatibility_action = "stop"
  )
  baseline <- categories$Z == 0 & categories$T_obs == treatment
  comparison <- categories$Z == 1 & categories$T_obs == treatment
  observed_arm <- as.numeric(categories$T_obs == treatment)
  instrument0 <- as.numeric(categories$Z == 0)
  instrument1 <- 1 - instrument0
  r0 <- if (treatment == 1L) specification$pi0 else 1 - specification$pi0
  r1 <- if (treatment == 1L) specification$pi1 else 1 - specification$pi1

  score <- loss <- q_matrix <- matrix(
    NA_real_, nrow = row_count, ncol = specification$candidates
  )
  conditional_names <- c(
    "m0", "Hbar0", "Hbar1", "D_r0", "D_r1", "M", "C0", "C1"
  )
  conditional_components <- array(
    NA_real_,
    dim = c(
      specification$strata, specification$candidates,
      length(conditional_names)
    ),
    dimnames = list(
      X = x_levels, candidate = seq_len(specification$candidates),
      component = conditional_names
    )
  )
  normalization_residual <- numeric(specification$candidates)
  direct_target <- numeric(specification$candidates)
  for (candidate in seq_len(specification$candidates)) {
    right <- specification$candidate_right_endpoint[
      candidate, categories$X
    ]
    q <- ifelse(
      right,
      specification$upper_interval[1L],
      specification$lower_interval[2L]
    )
    q_matrix[, candidate] <- q
    primitive <- compute_eif_row_primitives(
      data = categories,
      nuisance = nuisance,
      envelope = envelope,
      t_value = treatment,
      Lambda0 = Lambda,
      Lambda1 = Lambda,
      q_lower = q,
      q_upper = q
    )[[specification$endpoint_direction]]

    m0_x <- multibranch_weighted_group_mean(
      categories$Y[baseline], categories$probability[baseline],
      categories$X[baseline], x_levels
    )
    Hbar0_x <- multibranch_weighted_group_mean(
      primitive$H0[baseline], categories$probability[baseline],
      categories$X[baseline], x_levels
    )
    Hbar1_x <- multibranch_weighted_group_mean(
      primitive$H1[comparison], categories$probability[comparison],
      categories$X[comparison], x_levels
    )
    D_r0_x <- multibranch_weighted_group_mean(
      primitive$D_r0[baseline], categories$probability[baseline],
      categories$X[baseline], x_levels
    )
    D_r1_x <- multibranch_weighted_group_mean(
      primitive$D_r1[baseline], categories$probability[baseline],
      categories$X[baseline], x_levels
    )
    m0 <- m0_x[categories$X]
    Hbar0 <- Hbar0_x[categories$X]
    Hbar1 <- Hbar1_x[categories$X]
    linked <- compute_fixed_q_linked_components(
      m0 = m0,
      Hbar0 = Hbar0,
      Hbar1 = Hbar1,
      q = q,
      r0 = rep(r0, row_count),
      mean_D_r0 = D_r0_x[categories$X]
    )
    C1 <- (1 - r0) * D_r1_x[categories$X]
    conditional_components[, candidate, "m0"] <- m0_x
    conditional_components[, candidate, "Hbar0"] <- Hbar0_x
    conditional_components[, candidate, "Hbar1"] <- Hbar1_x
    conditional_components[, candidate, "D_r0"] <- D_r0_x
    conditional_components[, candidate, "D_r1"] <- D_r1_x
    conditional_components[, candidate, "M"] <-
      linked$M[match(x_levels, categories$X)]
    conditional_components[, candidate, "C0"] <-
      linked$C0[match(x_levels, categories$X)]
    conditional_components[, candidate, "C1"] <-
      C1[match(x_levels, categories$X)]
    treatment_correction <-
      instrument0 / specification$probability_z0 * linked$C0 *
      (observed_arm - r0) +
      instrument1 / (1 - specification$probability_z0) * C1 *
      (observed_arm - r1)
    outcome_correction <-
      instrument0 / specification$probability_z0 * observed_arm / r0 *
      (primitive$H0 - Hbar0) +
      instrument1 / (1 - specification$probability_z0) *
      observed_arm / r1 * (primitive$H1 - Hbar1)
    score[, candidate] <- linked$M + treatment_correction +
      outcome_correction
    threshold_level <- if (
      specification$endpoint_direction == "lower"
    ) lower_level else tau
    loss[, candidate] <- ifelse(
      baseline,
      envelope$D * check_loss(
        categories$Y, q, rep(threshold_level, row_count)
      ),
      0
    )

    lambda <- if (specification$endpoint_direction == "lower") {
      ifelse(categories$Y <= q, envelope$U, envelope$L)
    } else {
      ifelse(categories$Y <= q, envelope$L, envelope$U)
    }
    lambda_x <- multibranch_weighted_group_mean(
      lambda[baseline], categories$probability[baseline],
      categories$X[baseline], x_levels
    )
    theta_x <- multibranch_weighted_group_mean(
      (lambda * categories$Y)[baseline],
      categories$probability[baseline],
      categories$X[baseline], x_levels
    )
    normalization_residual[candidate] <- max(abs(lambda_x - 1))
    direct_target[candidate] <- mean(r0 * m0_x + (1 - r0) * theta_x)
  }

  population_mean <- colSums(categories$probability * score)
  population_risk <- colSums(categories$probability * loss)
  centered <- sweep(score, 2L, population_mean, "-")
  population_covariance <- crossprod(
    sqrt(categories$probability) * centered
  )
  population_sd <- sqrt(diag(population_covariance))
  population_correlation <- population_covariance /
    tcrossprod(population_sd)
  if (
    max(population_mean) - min(population_mean) > 1e-10 ||
      max(population_risk) - min(population_risk) > 1e-10 ||
      max(direct_target) - min(direct_target) > 1e-10 ||
      max(normalization_residual) > 1e-10 ||
      max(abs(population_mean - direct_target)) > 1e-10 ||
      any(!is.finite(population_sd)) || any(population_sd <= 0)
  ) {
    stop("The continuous-gap oracle failed an exact population identity.")
  }
  list(
    specification = specification,
    categories = categories,
    q = q_matrix,
    score = score,
    selection_loss = loss,
    target = mean(direct_target),
    direct_target = direct_target,
    direct_normalization_residual = normalization_residual,
    population_mean = population_mean,
    population_risk = population_risk,
    population_sd = population_sd,
    population_correlation = population_correlation,
    conditional_components = conditional_components,
    tau = as.numeric(tau),
    lower_level = as.numeric(lower_level),
    threshold_level = as.numeric(if (
      specification$endpoint_direction == "lower"
    ) lower_level else tau),
    lower_component_probability = as.numeric(lower_component_probability),
    component_probability = component_probability
  )
}

evaluate_fivefold_gap_nonunique_threshold_rows <- function(data, oracle) {
  specification <- oracle$specification
  required <- c("Y", "T_obs", "Z", "X")
  if (!is.data.frame(data) || !all(required %in% names(data))) {
    stop("Continuous-gap evaluation data are incomplete.")
  }
  n <- nrow(data)
  x <- as.integer(data$X)
  if (any(!x %in% seq_len(specification$strata))) {
    stop("Continuous-gap X values are outside the oracle strata.")
  }
  Lambda <- specification$sensitivity_lambda
  treatment <- specification$treatment
  nuisance <- list(
    e = rep(specification$probability_z0, n),
    pi0 = rep(specification$pi0, n),
    pi1 = rep(specification$pi1, n)
  )
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = specification$gamma_components[ifelse(
      data$Y <= specification$lower_interval[2L], 1L, 2L
    )],
    t_value = treatment,
    Lambda0 = Lambda,
    Lambda1 = Lambda,
    envelope_smoothing_epsilon =
      specification$envelope_smoothing_epsilon,
    gamma_compatibility_action = "stop"
  )
  observed_arm <- as.numeric(data$T_obs == treatment)
  instrument0 <- as.numeric(data$Z == 0)
  instrument1 <- 1 - instrument0
  baseline <- data$Z == 0 & data$T_obs == treatment
  r0 <- if (treatment == 1L) specification$pi0 else 1 - specification$pi0
  r1 <- if (treatment == 1L) specification$pi1 else 1 - specification$pi1
  score <- loss <- matrix(
    NA_real_, nrow = n, ncol = specification$candidates
  )
  for (candidate in seq_len(specification$candidates)) {
    right <- specification$candidate_right_endpoint[candidate, x]
    q <- ifelse(
      right,
      specification$upper_interval[1L],
      specification$lower_interval[2L]
    )
    primitive <- compute_eif_row_primitives(
      data = data,
      nuisance = nuisance,
      envelope = envelope,
      t_value = treatment,
      Lambda0 = Lambda,
      Lambda1 = Lambda,
      q_lower = q,
      q_upper = q
    )[[specification$endpoint_direction]]
    component <- oracle$conditional_components[x, candidate, , drop = FALSE]
    dim(component) <- c(n, dim(oracle$conditional_components)[3L])
    colnames(component) <- dimnames(oracle$conditional_components)$component
    treatment_correction <-
      instrument0 / specification$probability_z0 * component[, "C0"] *
      (observed_arm - r0) +
      instrument1 / (1 - specification$probability_z0) *
      component[, "C1"] * (observed_arm - r1)
    outcome_correction <-
      instrument0 / specification$probability_z0 * observed_arm / r0 *
      (primitive$H0 - component[, "Hbar0"]) +
      instrument1 / (1 - specification$probability_z0) *
      observed_arm / r1 * (primitive$H1 - component[, "Hbar1"])
    score[, candidate] <- component[, "M"] + treatment_correction +
      outcome_correction
    loss[, candidate] <- ifelse(
      baseline,
      envelope$D * check_loss(
        data$Y, q, rep(oracle$threshold_level, n)
      ),
      0
    )
  }
  list(score = score, selection_loss = loss)
}

generate_fivefold_gap_nonunique_threshold_data <- function(
    n,
    specification = make_fivefold_gap_nonunique_threshold_spec(),
    seed = NULL,
    oracle = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  if (is.null(oracle)) {
    oracle <- build_fivefold_gap_nonunique_threshold_oracle(specification)
  }
  x <- sample.int(specification$strata, n, replace = TRUE)
  z <- stats::rbinom(n, 1L, 1 - specification$probability_z0)
  treatment_probability <- ifelse(
    z == 0, specification$pi0, specification$pi1
  )
  treatment <- stats::rbinom(n, 1L, treatment_probability)
  lower_probability <- rep(oracle$lower_component_probability, n)
  comparison_arm <- z == 1L & treatment == specification$treatment
  lower_probability[comparison_arm] <-
    oracle$lower_component_probability *
      specification$gamma_components[1L]
  lower <- stats::rbinom(n, 1L, lower_probability) == 1L
  y <- numeric(n)
  y[lower] <- stats::runif(
    sum(lower), specification$lower_interval[1L],
    specification$lower_interval[2L]
  )
  y[!lower] <- stats::runif(
    sum(!lower), specification$upper_interval[1L],
    specification$upper_interval[2L]
  )
  data.frame(
    Y = y,
    T_obs = treatment,
    Z = z,
    X = x,
    e_true = specification$probability_z0,
    pi0_true = specification$pi0,
    pi1_true = specification$pi1,
    gamma0_oracle = 1,
    gamma1_oracle = 1,
    prob_Z = 1 - specification$probability_z0,
    prob_T = treatment_probability
  )
}

multibranch_select_candidate <- function(loss_matrix, rows) {
  if (length(rows) == 0L) stop("Candidate selection requires training rows.")
  risk <- colMeans(loss_matrix[rows, , drop = FALSE])
  which.min(risk)
}

multibranch_crossfit_label <- function(folds) {
  folds <- as.integer(folds)
  if (folds == 5L) "five_fold_crossfit" else sprintf("%d_fold_crossfit", folds)
}

multibranch_endpoint_direction <- function(oracle) {
  direction <- oracle$specification$endpoint_direction
  if (is.null(direction)) direction <- "lower"
  match.arg(direction, c("lower", "upper"))
}

multibranch_outward_covered <- function(
    z,
    endpoint_direction = c("lower", "upper"),
    confidence_level = 0.95
) {
  endpoint_direction <- match.arg(endpoint_direction)
  quantile <- stats::qnorm(confidence_level)
  if (endpoint_direction == "lower") {
    z <= quantile
  } else {
    z >= -quantile
  }
}

multibranch_crossfit_estimate <- function(
    category_index,
    oracle,
    folds = 5L,
    fold_id = NULL,
    confidence_level = 0.95
) {
  n <- length(category_index)
  folds <- as.integer(folds)
  if (is.null(fold_id)) {
    fold_id <- sample(rep(seq_len(folds), length.out = n))
  }
  fold_id <- as.integer(fold_id)
  if (length(fold_id) != n || any(!fold_id %in% seq_len(folds))) {
    stop("Invalid cross-fitting fold assignment.")
  }
  local_loss <- oracle$selection_loss[category_index, , drop = FALSE]
  local_score <- oracle$score[category_index, , drop = FALSE]
  signal <- numeric(n)
  selected <- integer(folds)
  for (fold in seq_len(folds)) {
    validation <- fold_id == fold
    selected[fold] <- multibranch_select_candidate(
      local_loss, which(!validation)
    )
    signal[validation] <- local_score[
      cbind(which(validation), selected[fold])
    ]
  }
  estimate <- mean(signal)
  standard_error <- stats::sd(signal) / sqrt(n)
  alpha <- 1 - confidence_level
  two_sided_quantile <- stats::qnorm(1 - alpha / 2)
  z_value <- (estimate - oracle$target) / standard_error
  list(
    estimate = estimate,
    standard_error = standard_error,
    z = z_value,
    centered_coverage = abs(estimate - oracle$target) <=
      two_sided_quantile * standard_error,
    outward_coverage = multibranch_outward_covered(
      z_value, multibranch_endpoint_direction(oracle), confidence_level
    ),
    selected = selected,
    signal = signal
  )
}

multibranch_same_sample_estimate <- function(
    category_index,
    oracle,
    confidence_level = 0.95
) {
  n <- length(category_index)
  local_loss <- oracle$selection_loss[category_index, , drop = FALSE]
  local_score <- oracle$score[category_index, , drop = FALSE]
  selected <- multibranch_select_candidate(local_loss, seq_len(n))
  signal <- local_score[, selected]
  estimate <- mean(signal)
  standard_error <- stats::sd(signal) / sqrt(n)
  alpha <- 1 - confidence_level
  two_sided_quantile <- stats::qnorm(1 - alpha / 2)
  z_value <- (estimate - oracle$target) / standard_error
  list(
    estimate = estimate,
    standard_error = standard_error,
    z = z_value,
    centered_coverage = abs(estimate - oracle$target) <=
      two_sided_quantile * standard_error,
    outward_coverage = multibranch_outward_covered(
      z_value, multibranch_endpoint_direction(oracle), confidence_level
    ),
    selected = selected,
    signal = signal
  )
}

multibranch_online_estimate <- function(
    category_index,
    oracle,
    initial_size = NULL,
    chunks = 10L,
    confidence_level = 0.95
) {
  n <- length(category_index)
  if (is.null(initial_size)) initial_size <- ceiling(sqrt(n))
  initial_size <- as.integer(initial_size)
  chunks <- as.integer(chunks)
  if (
    initial_size < 2L || initial_size >= n || chunks < 1L ||
      chunks > n - initial_size
  ) {
    stop("Invalid online initial size or chunk count.")
  }
  evaluation_rows <- seq.int(initial_size + 1L, n)
  chunk_id <- sort(rep(seq_len(chunks), length.out = length(evaluation_rows)))
  local_loss <- oracle$selection_loss[category_index, , drop = FALSE]
  local_score <- oracle$score[category_index, , drop = FALSE]
  numerator <- denominator <- 0
  selected <- integer(chunks)
  signal <- numeric(length(evaluation_rows))
  predictable_sd <- numeric(length(evaluation_rows))
  signal_index <- 1L
  for (chunk in seq_len(chunks)) {
    rows <- evaluation_rows[chunk_id == chunk]
    training_end <- min(rows) - 1L
    selected[chunk] <- multibranch_select_candidate(
      local_loss, seq_len(training_end)
    )
    values <- local_score[rows, selected[chunk]]
    sd_value <- oracle$population_sd[selected[chunk]]
    count <- length(rows)
    indices <- signal_index:(signal_index + count - 1L)
    signal[indices] <- values
    predictable_sd[indices] <- sd_value
    signal_index <- signal_index + count
    numerator <- numerator + sum(values / sd_value)
    denominator <- denominator + count / sd_value
  }
  estimate <- numerator / denominator
  evaluation_count <- length(evaluation_rows)
  harmonic_sd <- evaluation_count / denominator
  standard_error <- harmonic_sd / sqrt(evaluation_count)
  alpha <- 1 - confidence_level
  two_sided_quantile <- stats::qnorm(1 - alpha / 2)
  z_value <- (estimate - oracle$target) / standard_error
  list(
    estimate = estimate,
    standard_error = standard_error,
    z = z_value,
    centered_coverage = abs(estimate - oracle$target) <=
      two_sided_quantile * standard_error,
    outward_coverage = multibranch_outward_covered(
      z_value, multibranch_endpoint_direction(oracle), confidence_level
    ),
    selected = selected,
    signal = signal,
    predictable_sd = predictable_sd,
    initial_size = initial_size
  )
}

run_multibranch_nonunique_threshold_replication <- function(
    n,
    oracle,
    seed,
    folds = 5L,
    initial_size = NULL,
    chunks = 10L
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 20L) {
    stop("n must be one integer of at least 20.")
  }
  set.seed(as.integer(seed))
  category_index <- sample.int(
    nrow(oracle$categories), n, replace = TRUE,
    prob = oracle$categories$probability
  )
  same_sample <- multibranch_same_sample_estimate(category_index, oracle)
  crossfit <- multibranch_crossfit_estimate(
    category_index, oracle, folds = folds
  )
  online <- multibranch_online_estimate(
    category_index, oracle, initial_size = initial_size, chunks = chunks
  )
  data.frame(
    sample_size = n,
    seed = as.integer(seed),
    estimator = c(
      "same_sample_one_step",
      multibranch_crossfit_label(folds),
      "online_predictable"
    ),
    target = oracle$target,
    estimate = c(
      same_sample$estimate, crossfit$estimate, online$estimate
    ),
    standard_error = c(
      same_sample$standard_error,
      crossfit$standard_error,
      online$standard_error
    ),
    z = c(same_sample$z, crossfit$z, online$z),
    centered_coverage = c(
      same_sample$centered_coverage,
      crossfit$centered_coverage,
      online$centered_coverage
    ),
    outward_coverage = c(
      same_sample$outward_coverage,
      crossfit$outward_coverage,
      online$outward_coverage
    ),
    distinct_selected = c(
      1L,
      length(unique(crossfit$selected)),
      length(unique(online$selected))
    ),
    stringsAsFactors = FALSE
  )
}

run_fivefold_gap_nonunique_threshold_replication <- function(
    n,
    oracle,
    seed,
    folds = 5L,
    initial_size = NULL,
    chunks = 10L
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 20L) {
    stop("n must be one integer of at least 20.")
  }
  data <- generate_fivefold_gap_nonunique_threshold_data(
    n = n,
    specification = oracle$specification,
    seed = seed,
    oracle = oracle
  )
  evaluated <- evaluate_fivefold_gap_nonunique_threshold_rows(data, oracle)
  sample_oracle <- oracle
  sample_oracle$score <- evaluated$score
  sample_oracle$selection_loss <- evaluated$selection_loss
  rows <- seq_len(n)
  same_sample <- multibranch_same_sample_estimate(rows, sample_oracle)
  crossfit <- multibranch_crossfit_estimate(
    rows, sample_oracle, folds = folds
  )
  online <- multibranch_online_estimate(
    rows, sample_oracle, initial_size = initial_size, chunks = chunks
  )
  data.frame(
    sample_size = n,
    seed = as.integer(seed),
    estimator = c(
      "same_sample_one_step",
      multibranch_crossfit_label(folds),
      "online_predictable"
    ),
    target = oracle$target,
    estimate = c(
      same_sample$estimate, crossfit$estimate, online$estimate
    ),
    standard_error = c(
      same_sample$standard_error,
      crossfit$standard_error,
      online$standard_error
    ),
    z = c(same_sample$z, crossfit$z, online$z),
    centered_coverage = c(
      same_sample$centered_coverage,
      crossfit$centered_coverage,
      online$centered_coverage
    ),
    outward_coverage = c(
      same_sample$outward_coverage,
      crossfit$outward_coverage,
      online$outward_coverage
    ),
    distinct_selected = c(
      1L,
      length(unique(crossfit$selected)),
      length(unique(online$selected))
    ),
    stringsAsFactors = FALSE
  )
}

summarize_multibranch_nonunique_threshold_pilot <- function(detail) {
  required <- c(
    "sample_size", "estimator", "target", "estimate", "standard_error", "z",
    "centered_coverage", "outward_coverage", "distinct_selected"
  )
  if (!all(required %in% names(detail))) {
    stop("The pilot detail table is incomplete.")
  }
  groups <- split(detail, interaction(
    detail$sample_size, detail$estimator, drop = TRUE
  ))
  output <- do.call(rbind, lapply(groups, function(local) {
    data.frame(
      sample_size = local$sample_size[1L],
      estimator = local$estimator[1L],
      replications = nrow(local),
      target = local$target[1L],
      mean_estimate = mean(local$estimate),
      empirical_bias = mean(local$estimate - local$target),
      mean_standard_error = mean(local$standard_error),
      empirical_standard_deviation = stats::sd(local$estimate),
      z_mean = mean(local$z),
      z_standard_deviation = stats::sd(local$z),
      centered_coverage = mean(local$centered_coverage),
      outward_coverage = mean(local$outward_coverage),
      mean_distinct_selected = mean(local$distinct_selected),
      stringsAsFactors = FALSE
    )
  }))
  rownames(output) <- NULL
  output[order(output$sample_size, output$estimator), , drop = FALSE]
}

simulate_multibranch_crossfit_gaussian_limit <- function(
    oracle,
    replications = 20000L,
    folds = 5L,
    seed = 1L,
    confidence_level = 0.95,
    batch_size = 1000L
) {
  replications <- as.integer(replications)
  folds <- as.integer(folds)
  batch_size <- as.integer(batch_size)
  if (
    replications < 1L || folds < 2L || batch_size < 1L ||
      !all(c(
        "selection_loss", "score", "categories", "population_risk",
        "population_mean", "population_sd"
      ) %in% names(oracle))
  ) {
    stop("Invalid Gaussian-limit simulation inputs.")
  }
  centered_loss <- sweep(
    oracle$selection_loss, 2L, oracle$population_risk, "-"
  )
  centered_score <- sweep(
    oracle$score, 2L, oracle$population_mean, "-"
  )
  process <- cbind(centered_loss, centered_score)
  covariance <- crossprod(
    sqrt(oracle$categories$probability) * process
  )
  decomposition <- eigen(covariance, symmetric = TRUE)
  tolerance <- max(decomposition$values) * 1e-11
  keep <- decomposition$values > tolerance
  if (!any(keep)) stop("The Gaussian empirical process has zero rank.")
  factor <- sweep(
    decomposition$vectors[, keep, drop = FALSE],
    2L,
    sqrt(decomposition$values[keep]),
    "*"
  )
  candidate_count <- ncol(oracle$score)
  process_dimension <- nrow(factor)
  process_rank <- ncol(factor)
  z <- numeric(replications)
  set.seed(as.integer(seed))
  completed <- 0L
  while (completed < replications) {
    current <- min(batch_size, replications - completed)
    normal <- matrix(
      stats::rnorm(current * folds * process_rank),
      nrow = current * folds,
      ncol = process_rank
    )
    draw <- normal %*% t(factor)
    draw <- array(
      draw,
      dim = c(current, folds, process_dimension)
    )
    for (replication in seq_len(current)) {
      loss_process <- draw[
        replication, , seq_len(candidate_count), drop = FALSE
      ]
      dim(loss_process) <- c(folds, candidate_count)
      score_process <- draw[
        replication, , candidate_count + seq_len(candidate_count),
        drop = FALSE
      ]
      dim(score_process) <- c(folds, candidate_count)
      total_loss <- colSums(loss_process)
      selected <- vapply(seq_len(folds), function(fold) {
        which.min(total_loss - loss_process[fold, ])
      }, integer(1))
      selected_score <- score_process[cbind(seq_len(folds), selected)]
      z[completed + replication] <- sum(selected_score) /
        sqrt(sum(oracle$population_sd[selected]^2))
    }
    completed <- completed + current
  }
  alpha <- 1 - confidence_level
  data.frame(
    estimator = paste0(
      multibranch_crossfit_label(folds), "_gaussian_limit"
    ),
    replications = replications,
    folds = folds,
    process_rank = process_rank,
    z_mean = mean(z),
    z_standard_deviation = stats::sd(z),
    centered_coverage = mean(abs(z) <= stats::qnorm(1 - alpha / 2)),
    outward_coverage = mean(multibranch_outward_covered(
      z, multibranch_endpoint_direction(oracle), confidence_level
    )),
    stringsAsFactors = FALSE
  )
}

simulate_multibranch_interval_crossfit_gaussian_limit <- function(
    lower_oracle,
    upper_oracle,
    replications = 20000L,
    folds = 5L,
    seed = 1L,
    confidence_level = 0.95,
    batch_size = 1000L
) {
  replications <- as.integer(replications)
  folds <- as.integer(folds)
  batch_size <- as.integer(batch_size)
  candidate_count <- ncol(lower_oracle$score)
  if (
    replications < 1L || folds < 2L || batch_size < 1L ||
      candidate_count < 2L ||
      ncol(upper_oracle$score) != candidate_count ||
      nrow(lower_oracle$categories) != nrow(upper_oracle$categories) ||
      max(abs(
        lower_oracle$categories$probability -
          upper_oracle$categories$probability
      )) > 1e-12 ||
      multibranch_endpoint_direction(lower_oracle) != "lower" ||
      multibranch_endpoint_direction(upper_oracle) != "upper"
  ) {
    stop("Invalid paired-interval Gaussian-limit inputs.")
  }
  probability <- lower_oracle$categories$probability
  process <- cbind(
    sweep(
      lower_oracle$selection_loss, 2L,
      lower_oracle$population_risk, "-"
    ),
    sweep(
      upper_oracle$selection_loss, 2L,
      upper_oracle$population_risk, "-"
    ),
    sweep(
      lower_oracle$score, 2L,
      lower_oracle$population_mean, "-"
    ),
    sweep(
      upper_oracle$score, 2L,
      upper_oracle$population_mean, "-"
    )
  )
  covariance <- crossprod(sqrt(probability) * process)
  decomposition <- eigen(covariance, symmetric = TRUE)
  tolerance <- max(decomposition$values) * 1e-11
  keep <- decomposition$values > tolerance
  if (!any(keep)) stop("The paired-interval process has zero rank.")
  factor <- sweep(
    decomposition$vectors[, keep, drop = FALSE],
    2L, sqrt(decomposition$values[keep]), "*"
  )
  process_dimension <- nrow(factor)
  process_rank <- ncol(factor)
  z_lower <- z_upper <- numeric(replications)
  set.seed(as.integer(seed))
  completed <- 0L
  while (completed < replications) {
    current <- min(batch_size, replications - completed)
    normal <- matrix(
      stats::rnorm(current * folds * process_rank),
      nrow = current * folds, ncol = process_rank
    )
    draw <- normal %*% t(factor)
    draw <- array(draw, dim = c(current, folds, process_dimension))
    for (replication in seq_len(current)) {
      local <- draw[replication, , , drop = FALSE]
      dim(local) <- c(folds, process_dimension)
      lower_loss <- local[
        , seq_len(candidate_count), drop = FALSE
      ]
      upper_loss <- local[
        , candidate_count + seq_len(candidate_count), drop = FALSE
      ]
      lower_score <- local[
        , 2L * candidate_count + seq_len(candidate_count), drop = FALSE
      ]
      upper_score <- local[
        , 3L * candidate_count + seq_len(candidate_count), drop = FALSE
      ]
      total_lower_loss <- colSums(lower_loss)
      total_upper_loss <- colSums(upper_loss)
      selected_lower <- vapply(seq_len(folds), function(fold) {
        which.min(total_lower_loss - lower_loss[fold, ])
      }, integer(1))
      selected_upper <- vapply(seq_len(folds), function(fold) {
        which.min(total_upper_loss - upper_loss[fold, ])
      }, integer(1))
      row <- completed + replication
      z_lower[row] <- sum(lower_score[
        cbind(seq_len(folds), selected_lower)
      ]) / sqrt(sum(lower_oracle$population_sd[selected_lower]^2))
      z_upper[row] <- sum(upper_score[
        cbind(seq_len(folds), selected_upper)
      ]) / sqrt(sum(upper_oracle$population_sd[selected_upper]^2))
    }
    completed <- completed + current
  }
  alpha <- 1 - confidence_level
  pointwise_quantile <- stats::qnorm(confidence_level)
  simultaneous_quantile <- stats::qnorm(1 - alpha / 2)
  lower_outward <- z_lower <= pointwise_quantile
  upper_outward <- z_upper >= -pointwise_quantile
  data.frame(
    estimator = paste0(
      "sharp_interval_", multibranch_crossfit_label(folds),
      "_gaussian_limit"
    ),
    replications = replications,
    folds = folds,
    candidates = candidate_count,
    process_rank = process_rank,
    lower_z_mean = mean(z_lower),
    lower_z_standard_deviation = stats::sd(z_lower),
    upper_z_mean = mean(z_upper),
    upper_z_standard_deviation = stats::sd(z_upper),
    lower_centered_coverage = mean(
      abs(z_lower) <= stats::qnorm(1 - alpha / 2)
    ),
    upper_centered_coverage = mean(
      abs(z_upper) <= stats::qnorm(1 - alpha / 2)
    ),
    lower_outward_coverage = mean(lower_outward),
    upper_outward_coverage = mean(upper_outward),
    outward_set_coverage = mean(lower_outward & upper_outward),
    outward_set_coverage_bonferroni = mean(
      z_lower <= simultaneous_quantile & z_upper >= -simultaneous_quantile
    ),
    endpoint_z_correlation = stats::cor(z_lower, z_upper),
    stringsAsFactors = FALSE
  )
}

simulate_multibranch_ate_lower_crossfit_gaussian_limit <- function(
    oracle,
    replications = 20000L,
    folds = 5L,
    seed = 1L,
    confidence_level = 0.95,
    batch_size = 1000L
) {
  replications <- as.integer(replications)
  folds <- as.integer(folds)
  batch_size <- as.integer(batch_size)
  required <- c(
    "categories", "treated_lower", "untreated_upper", "pair_sd",
    "candidate_count"
  )
  if (
    replications < 1L || folds < 2L || batch_size < 1L ||
      !all(required %in% names(oracle))
  ) {
    stop("Invalid ATE-lower Gaussian-limit simulation inputs.")
  }
  treated <- oracle$treated_lower
  untreated <- oracle$untreated_upper
  candidate_count <- oracle$candidate_count
  centered_treated_loss <- sweep(
    treated$selection_loss, 2L, treated$population_risk, "-"
  )
  centered_untreated_loss <- sweep(
    untreated$selection_loss, 2L, untreated$population_risk, "-"
  )
  centered_treated_score <- sweep(
    treated$score, 2L, treated$population_mean, "-"
  )
  centered_untreated_score <- sweep(
    untreated$score, 2L, untreated$population_mean, "-"
  )
  process <- cbind(
    centered_treated_loss,
    centered_untreated_loss,
    centered_treated_score,
    centered_untreated_score
  )
  covariance <- crossprod(
    sqrt(oracle$categories$probability) * process
  )
  decomposition <- eigen(covariance, symmetric = TRUE)
  tolerance <- max(decomposition$values) * 1e-11
  keep <- decomposition$values > tolerance
  if (!any(keep)) stop("The ATE-lower Gaussian process has zero rank.")
  factor <- sweep(
    decomposition$vectors[, keep, drop = FALSE],
    2L,
    sqrt(decomposition$values[keep]),
    "*"
  )
  process_dimension <- nrow(factor)
  process_rank <- ncol(factor)
  z <- numeric(replications)
  selected_treated_count <- numeric(replications)
  selected_untreated_count <- numeric(replications)
  set.seed(as.integer(seed))
  completed <- 0L
  while (completed < replications) {
    current <- min(batch_size, replications - completed)
    normal <- matrix(
      stats::rnorm(current * folds * process_rank),
      nrow = current * folds,
      ncol = process_rank
    )
    draw <- normal %*% t(factor)
    draw <- array(draw, dim = c(current, folds, process_dimension))
    for (replication in seq_len(current)) {
      local <- draw[replication, , , drop = FALSE]
      dim(local) <- c(folds, process_dimension)
      treated_loss <- local[, seq_len(candidate_count), drop = FALSE]
      untreated_loss <- local[
        , candidate_count + seq_len(candidate_count), drop = FALSE
      ]
      treated_score <- local[
        , 2L * candidate_count + seq_len(candidate_count), drop = FALSE
      ]
      untreated_score <- local[
        , 3L * candidate_count + seq_len(candidate_count), drop = FALSE
      ]
      total_treated_loss <- colSums(treated_loss)
      total_untreated_loss <- colSums(untreated_loss)
      selected_treated <- vapply(seq_len(folds), function(fold) {
        which.min(total_treated_loss - treated_loss[fold, ])
      }, integer(1))
      selected_untreated <- vapply(seq_len(folds), function(fold) {
        which.min(total_untreated_loss - untreated_loss[fold, ])
      }, integer(1))
      selected_score <-
        treated_score[cbind(seq_len(folds), selected_treated)] -
        untreated_score[cbind(seq_len(folds), selected_untreated)]
      selected_sd <- oracle$pair_sd[
        cbind(selected_treated, selected_untreated)
      ]
      z[completed + replication] <- sum(selected_score) /
        sqrt(sum(selected_sd^2))
      selected_treated_count[completed + replication] <-
        length(unique(selected_treated))
      selected_untreated_count[completed + replication] <-
        length(unique(selected_untreated))
    }
    completed <- completed + current
  }
  alpha <- 1 - confidence_level
  data.frame(
    estimator = paste0(
      "ate_lower_", multibranch_crossfit_label(folds), "_gaussian_limit"
    ),
    replications = replications,
    folds = folds,
    process_rank = process_rank,
    z_mean = mean(z),
    z_standard_deviation = stats::sd(z),
    centered_coverage = mean(abs(z) <= stats::qnorm(1 - alpha / 2)),
    outward_coverage = mean(z <= stats::qnorm(confidence_level)),
    mean_distinct_treated_selected = mean(selected_treated_count),
    mean_distinct_untreated_selected = mean(selected_untreated_count),
    stringsAsFactors = FALSE
  )
}

build_multibranch_ate_lower_oracle <- function(
    sensitivity_lambda = 4,
    strata = 4L,
    candidate_design = "complete_binary",
    probability_z0 = 0.5,
    pi0 = 0.3,
    pi1 = 0.1,
    atoms = c(0, 1, 1.5),
    envelope_smoothing_epsilon = 0
) {
  sensitivity_lambda <- as.numeric(sensitivity_lambda)
  symmetric_tail_probability <- 1 / (sensitivity_lambda + 1)
  common <- list(
    sensitivity_lambda = sensitivity_lambda,
    strata = strata,
    candidates = max(2L, as.integer(strata)),
    candidate_design = candidate_design,
    probability_z0 = probability_z0,
    pi0 = pi0,
    pi1 = pi1,
    lower_atom = atoms[1L],
    boundary_atom = atoms[2L],
    tail_atom = atoms[3L],
    tail_probability = symmetric_tail_probability,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon
  )
  treated_lower_spec <- do.call(
    make_multibranch_nonunique_threshold_spec,
    c(common, list(treatment = 1L, endpoint_direction = "lower"))
  )
  untreated_upper_spec <- do.call(
    make_multibranch_nonunique_threshold_spec,
    c(common, list(treatment = 0L, endpoint_direction = "upper"))
  )
  treated_lower <- build_multibranch_nonunique_threshold_oracle(
    treated_lower_spec
  )
  untreated_upper <- build_multibranch_nonunique_threshold_oracle(
    untreated_upper_spec
  )
  if (
    nrow(treated_lower$categories) != nrow(untreated_upper$categories) ||
      max(abs(
        treated_lower$categories$probability -
          untreated_upper$categories$probability
      )) > 1e-12 ||
      !identical(
        treated_lower$specification$candidate_right_endpoint,
        untreated_upper$specification$candidate_right_endpoint
      )
  ) {
    stop("The two ATE-lower oracle arms do not share one observed-data law.")
  }
  candidate_count <- ncol(treated_lower$score)
  pair_sd <- matrix(
    NA_real_, nrow = candidate_count, ncol = candidate_count
  )
  probability <- treated_lower$categories$probability
  for (treated in seq_len(candidate_count)) {
    for (untreated in seq_len(candidate_count)) {
      value <- treated_lower$score[, treated] -
        untreated_upper$score[, untreated]
      pair_sd[treated, untreated] <- sqrt(
        sum(probability * (value - sum(probability * value))^2)
      )
    }
  }
  if (any(!is.finite(pair_sd)) || any(pair_sd <= 0)) {
    stop("The ATE-lower candidate pairs must have positive finite variance.")
  }
  list(
    categories = treated_lower$categories,
    treated_lower = treated_lower,
    untreated_upper = untreated_upper,
    target = treated_lower$target - untreated_upper$target,
    pair_sd = pair_sd,
    candidate_count = candidate_count
  )
}

build_fivefold_gap_ate_oracle <- function(
    sensitivity_lambda = 4,
    strata = 3L,
    candidate_design = "complete_binary",
    endpoint = c("lower", "upper"),
    probability_z0 = 0.5,
    pi0 = 0.7,
    pi1 = 0.9,
    treated_lower_interval = c(-0.001, 0),
    treated_upper_interval = c(1, 1.001),
    untreated_lower_interval = c(-0.001, 0),
    untreated_upper_interval = c(1, 1.001),
    quadrature_nodes_per_interval = 24L,
    envelope_smoothing_epsilon = 0
) {
  endpoint <- match.arg(endpoint)
  treated_direction <- if (endpoint == "lower") "lower" else "upper"
  untreated_direction <- if (endpoint == "lower") "upper" else "lower"
  common <- list(
    sensitivity_lambda = sensitivity_lambda,
    strata = strata,
    candidate_design = candidate_design,
    probability_z0 = probability_z0,
    pi0 = pi0,
    pi1 = pi1,
    quadrature_nodes_per_interval = quadrature_nodes_per_interval,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon
  )
  treated_specification <- do.call(
    make_fivefold_gap_nonunique_threshold_spec,
    c(common, list(
      endpoint_direction = treated_direction,
      treatment = 1L,
      lower_interval = treated_lower_interval,
      upper_interval = treated_upper_interval
    ))
  )
  untreated_specification <- do.call(
    make_fivefold_gap_nonunique_threshold_spec,
    c(common, list(
      endpoint_direction = untreated_direction,
      treatment = 0L,
      lower_interval = untreated_lower_interval,
      upper_interval = untreated_upper_interval
    ))
  )
  treated_original <- build_fivefold_gap_nonunique_threshold_oracle(
    treated_specification
  )
  untreated_original <- build_fivefold_gap_nonunique_threshold_oracle(
    untreated_specification
  )
  if (
    treated_original$specification$candidates !=
      untreated_original$specification$candidates ||
      !identical(
        treated_original$specification$candidate_right_endpoint,
        untreated_original$specification$candidate_right_endpoint
      )
  ) {
    stop("The two continuous-gap ATE arms have incompatible libraries.")
  }

  treated_rows <- which(treated_original$categories$T_obs == 1L)
  untreated_rows <- which(untreated_original$categories$T_obs == 0L)
  joint_categories <- rbind(
    transform(
      treated_original$categories[treated_rows, , drop = FALSE],
      outcome_arm = "treated"
    ),
    transform(
      untreated_original$categories[untreated_rows, , drop = FALSE],
      outcome_arm = "untreated"
    )
  )
  rownames(joint_categories) <- NULL
  if (abs(sum(joint_categories$probability) - 1) > 1e-10) {
    stop("The continuous-gap ATE joint law does not integrate to one.")
  }
  candidate_count <- treated_original$specification$candidates
  row_count <- nrow(joint_categories)
  treated_score <- untreated_score <- matrix(
    NA_real_, nrow = row_count, ncol = candidate_count
  )
  treated_loss <- untreated_loss <- matrix(
    0, nrow = row_count, ncol = candidate_count
  )
  joint_treated <- joint_categories$T_obs == 1L
  joint_untreated <- !joint_treated
  treated_score[joint_treated, ] <-
    treated_original$score[treated_rows, , drop = FALSE]
  treated_loss[joint_treated, ] <-
    treated_original$selection_loss[treated_rows, , drop = FALSE]
  untreated_score[joint_untreated, ] <-
    untreated_original$score[untreated_rows, , drop = FALSE]
  untreated_loss[joint_untreated, ] <-
    untreated_original$selection_loss[untreated_rows, , drop = FALSE]

  make_key <- function(data) paste(data$X, data$Z, data$T_obs, sep = ":")
  treated_cross_reference <- treated_original$categories[
    treated_original$categories$T_obs == 0L &
      treated_original$categories$outcome_index == 1L,
    , drop = FALSE
  ]
  untreated_cross_reference <- untreated_original$categories[
    untreated_original$categories$T_obs == 1L &
      untreated_original$categories$outcome_index == 1L,
    , drop = FALSE
  ]
  treated_cross_rows <- match(
    make_key(joint_categories[joint_untreated, , drop = FALSE]),
    make_key(treated_cross_reference)
  )
  untreated_cross_rows <- match(
    make_key(joint_categories[joint_treated, , drop = FALSE]),
    make_key(untreated_cross_reference)
  )
  if (anyNA(treated_cross_rows) || anyNA(untreated_cross_rows)) {
    stop("Failed to align the two continuous-gap ATE score laws.")
  }
  treated_reference_indices <- as.integer(rownames(
    treated_cross_reference[treated_cross_rows, , drop = FALSE]
  ))
  untreated_reference_indices <- as.integer(rownames(
    untreated_cross_reference[untreated_cross_rows, , drop = FALSE]
  ))
  treated_score[joint_untreated, ] <- treated_original$score[
    treated_reference_indices, , drop = FALSE
  ]
  untreated_score[joint_treated, ] <- untreated_original$score[
    untreated_reference_indices, , drop = FALSE
  ]

  make_joint_arm <- function(score, loss, original) {
    probability <- joint_categories$probability
    population_mean <- colSums(probability * score)
    population_risk <- colSums(probability * loss)
    centered <- sweep(score, 2L, population_mean, "-")
    covariance <- crossprod(sqrt(probability) * centered)
    standard_deviation <- sqrt(diag(covariance))
    if (
      max(population_mean) - min(population_mean) > 1e-10 ||
        max(population_risk) - min(population_risk) > 1e-10 ||
        max(abs(population_mean - original$target)) > 1e-10
    ) {
      stop("A continuous-gap ATE arm lost its exact population tie.")
    }
    original$categories <- joint_categories
    original$score <- score
    original$selection_loss <- loss
    original$population_mean <- population_mean
    original$population_risk <- population_risk
    original$population_sd <- standard_deviation
    original$population_correlation <- stats::cov2cor(covariance)
    original
  }
  treated <- make_joint_arm(
    treated_score, treated_loss, treated_original
  )
  untreated <- make_joint_arm(
    untreated_score, untreated_loss, untreated_original
  )
  pair_sd <- matrix(
    NA_real_, nrow = candidate_count, ncol = candidate_count
  )
  probability <- joint_categories$probability
  for (treated_candidate in seq_len(candidate_count)) {
    for (untreated_candidate in seq_len(candidate_count)) {
      value <- treated$score[, treated_candidate] -
        untreated$score[, untreated_candidate]
      pair_sd[treated_candidate, untreated_candidate] <- sqrt(
        sum(probability * (value - sum(probability * value))^2)
      )
    }
  }
  if (any(!is.finite(pair_sd)) || any(pair_sd <= 0)) {
    stop("The continuous-gap ATE pairs need positive finite variances.")
  }
  list(
    name = paste0("fivefold_continuous_gap_ate_", endpoint),
    endpoint = endpoint,
    categories = joint_categories,
    treated_lower = treated,
    untreated_upper = untreated,
    treated = treated,
    untreated = untreated,
    treated_direction = treated_direction,
    untreated_direction = untreated_direction,
    target = treated$target - untreated$target,
    pair_sd = pair_sd,
    candidate_count = candidate_count
  )
}

multibranch_ate_lower_select <- function(local_loss, training_rows) {
  c(
    treated = multibranch_select_candidate(
      local_loss$treated, training_rows
    ),
    untreated = multibranch_select_candidate(
      local_loss$untreated, training_rows
    )
  )
}

multibranch_ate_lower_crossfit_estimate <- function(
    category_index,
    oracle,
    folds = 5L,
    fold_id = NULL,
    confidence_level = 0.95
) {
  n <- length(category_index)
  folds <- as.integer(folds)
  if (is.null(fold_id)) {
    fold_id <- sample(rep(seq_len(folds), length.out = n))
  }
  local_loss <- list(
    treated = oracle$treated_lower$selection_loss[
      category_index, , drop = FALSE
    ],
    untreated = oracle$untreated_upper$selection_loss[
      category_index, , drop = FALSE
    ]
  )
  treated_score <- oracle$treated_lower$score[
    category_index, , drop = FALSE
  ]
  untreated_score <- oracle$untreated_upper$score[
    category_index, , drop = FALSE
  ]
  signal <- numeric(n)
  selected <- matrix(
    NA_integer_, nrow = folds, ncol = 2L,
    dimnames = list(NULL, c("treated", "untreated"))
  )
  for (fold in seq_len(folds)) {
    validation <- fold_id == fold
    selected[fold, ] <- multibranch_ate_lower_select(
      local_loss, which(!validation)
    )
    rows <- which(validation)
    signal[rows] <-
      treated_score[cbind(rows, selected[fold, "treated"])] -
      untreated_score[cbind(rows, selected[fold, "untreated"])]
  }
  estimate <- mean(signal)
  standard_error <- stats::sd(signal) / sqrt(n)
  alpha <- 1 - confidence_level
  list(
    estimate = estimate,
    standard_error = standard_error,
    z = (estimate - oracle$target) / standard_error,
    centered_coverage = abs(estimate - oracle$target) <=
      stats::qnorm(1 - alpha / 2) * standard_error,
    outward_coverage = estimate - stats::qnorm(confidence_level) *
      standard_error <= oracle$target,
    selected = selected,
    signal = signal
  )
}

multibranch_ate_lower_online_estimate <- function(
    category_index,
    oracle,
    initial_size = NULL,
    chunks = 10L,
    confidence_level = 0.95
) {
  n <- length(category_index)
  if (is.null(initial_size)) initial_size <- ceiling(sqrt(n))
  initial_size <- as.integer(initial_size)
  chunks <- as.integer(chunks)
  if (
    initial_size < 2L || initial_size >= n || chunks < 1L ||
      chunks > n - initial_size
  ) {
    stop("Invalid ATE-lower online initial size or chunk count.")
  }
  evaluation_rows <- seq.int(initial_size + 1L, n)
  chunk_id <- sort(rep(seq_len(chunks), length.out = length(evaluation_rows)))
  local_loss <- list(
    treated = oracle$treated_lower$selection_loss[
      category_index, , drop = FALSE
    ],
    untreated = oracle$untreated_upper$selection_loss[
      category_index, , drop = FALSE
    ]
  )
  treated_score <- oracle$treated_lower$score[
    category_index, , drop = FALSE
  ]
  untreated_score <- oracle$untreated_upper$score[
    category_index, , drop = FALSE
  ]
  numerator <- denominator <- 0
  selected <- matrix(
    NA_integer_, nrow = chunks, ncol = 2L,
    dimnames = list(NULL, c("treated", "untreated"))
  )
  signal <- numeric(length(evaluation_rows))
  signal_index <- 1L
  for (chunk in seq_len(chunks)) {
    rows <- evaluation_rows[chunk_id == chunk]
    training_end <- min(rows) - 1L
    selected[chunk, ] <- multibranch_ate_lower_select(
      local_loss, seq_len(training_end)
    )
    values <-
      treated_score[rows, selected[chunk, "treated"]] -
      untreated_score[rows, selected[chunk, "untreated"]]
    sd_value <- oracle$pair_sd[
      selected[chunk, "treated"], selected[chunk, "untreated"]
    ]
    count <- length(rows)
    indices <- signal_index:(signal_index + count - 1L)
    signal[indices] <- values
    signal_index <- signal_index + count
    numerator <- numerator + sum(values / sd_value)
    denominator <- denominator + count / sd_value
  }
  estimate <- numerator / denominator
  evaluation_count <- length(evaluation_rows)
  standard_error <- (evaluation_count / denominator) /
    sqrt(evaluation_count)
  alpha <- 1 - confidence_level
  list(
    estimate = estimate,
    standard_error = standard_error,
    z = (estimate - oracle$target) / standard_error,
    centered_coverage = abs(estimate - oracle$target) <=
      stats::qnorm(1 - alpha / 2) * standard_error,
    outward_coverage = estimate - stats::qnorm(confidence_level) *
      standard_error <= oracle$target,
    selected = selected,
    signal = signal,
    initial_size = initial_size
  )
}

run_multibranch_ate_lower_replication <- function(
    n,
    oracle,
    seed,
    folds = 5L,
    initial_size = NULL,
    chunks = 10L
) {
  set.seed(as.integer(seed))
  category_index <- sample.int(
    nrow(oracle$categories), as.integer(n), replace = TRUE,
    prob = oracle$categories$probability
  )
  crossfit <- multibranch_ate_lower_crossfit_estimate(
    category_index, oracle, folds = folds
  )
  online <- multibranch_ate_lower_online_estimate(
    category_index, oracle, initial_size = initial_size, chunks = chunks
  )
  data.frame(
    sample_size = as.integer(n),
    seed = as.integer(seed),
    estimator = c("five_fold_crossfit", "online_predictable"),
    target = oracle$target,
    estimate = c(crossfit$estimate, online$estimate),
    standard_error = c(crossfit$standard_error, online$standard_error),
    z = c(crossfit$z, online$z),
    centered_coverage = c(
      crossfit$centered_coverage, online$centered_coverage
    ),
    outward_coverage = c(
      crossfit$outward_coverage, online$outward_coverage
    ),
    stringsAsFactors = FALSE
  )
}
