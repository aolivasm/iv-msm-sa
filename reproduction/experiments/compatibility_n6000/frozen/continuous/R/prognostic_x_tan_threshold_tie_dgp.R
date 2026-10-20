# Prognostic-X threshold-tie DGP for Tan CAL and dual-IV inference.
#
# Four equally likely X strata have shifted outcome supports.  Within every
# (X,Z) cell, the treated-outcome CDF equals 0.20 after atom 2 and 0.80 after
# atom 4.  The dual IV-MSM normalization equations are also exact, while all
# pointwise envelopes remain separated from their max/min branch kinks.

if (!exists("build_selected_two_gap_nonunique_threshold_oracle", mode = "function")) {
  if (!file.exists("R/two_gap_nonunique_threshold_dgp.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/two_gap_nonunique_threshold_dgp.R")
}

make_prognostic_x_tan_threshold_tie_specs <- function(
    shifts = c(-3, -1, 1, 3),
    candidate_design = c("complete_binary", "hadamard")
) {
  candidate_design <- match.arg(candidate_design)
  shifts <- as.numeric(shifts)
  strata <- length(shifts)
  if (strata < 2L || any(!is.finite(shifts))) {
    stop("shifts must contain at least two finite values.")
  }
  if (candidate_design == "hadamard" &&
      bitwAnd(strata, strata - 1L) != 0L) {
    stop("Hadamard prognostic-X designs require power-of-two strata.")
  }
  selected <- make_selected_two_gap_nonunique_threshold_specs(4L)$lower
  common <- list(
    sensitivity_lambda = selected$sensitivity_lambda,
    strata = strata,
    candidates = strata,
    candidate_design = candidate_design,
    treatment = selected$treatment,
    probability_z0 = selected$probability_z0,
    pi0 = selected$pi0,
    pi1 = selected$pi1,
    atoms = selected$atoms,
    threshold_gap_fractions = selected$threshold_gap_fractions,
    gamma_atoms = selected$gamma_atoms,
    atom_probability_override = selected$atom_probability_override,
    atom_shift_by_stratum = shifts,
    envelope_smoothing_epsilon = 0
  )
  list(
    lower = do.call(
      make_multibranch_nonunique_threshold_spec,
      c(common, list(endpoint_direction = "lower", gap_index = 2L))
    ),
    upper = do.call(
      make_multibranch_nonunique_threshold_spec,
      c(common, list(endpoint_direction = "upper", gap_index = 4L))
    )
  )
}

build_prognostic_x_tan_threshold_tie_oracle <- function(
    shifts = c(-3, -1, 1, 3),
    candidate_design = "complete_binary"
) {
  specifications <- make_prognostic_x_tan_threshold_tie_specs(
    shifts, candidate_design
  )
  lower <- build_multibranch_nonunique_threshold_oracle(
    specifications$lower
  )
  upper <- build_multibranch_nonunique_threshold_oracle(
    specifications$upper
  )
  if (max(abs(lower$categories$probability -
      upper$categories$probability)) > 1e-12) {
    stop("The prognostic-X endpoints do not share one observed law.")
  }
  list(
    name = "prognostic_x_joint_dual_tan_threshold_ties",
    specifications = specifications,
    categories = lower$categories,
    lower = lower,
    upper = upper
  )
}

compute_prognostic_x_tan_population <- function(
    oracle = build_prognostic_x_tan_threshold_tie_oracle()
) {
  if (!exists("observed_dgp_standard_msm_arm", mode = "function")) {
    if (!file.exists("R/observed_data_dgp.R")) {
      stop("Run this function from the repository root.")
    }
    source("R/observed_data_dgp.R")
  }
  specification <- oracle$lower$specification
  Lambda <- specification$sensitivity_lambda
  probability0 <- specification$atom_probability_override
  probability1 <- probability0 * specification$gamma_atoms
  probability1 <- probability1 / sum(probability1)
  e <- specification$probability_z0
  pi0 <- specification$pi0
  pi1 <- specification$pi1
  marginal_pi <- e * pi0 + (1 - e) * pi1
  probability_z0_given_treated <- e * pi0 / marginal_pi
  marginal_probability <-
    probability_z0_given_treated * probability0 +
    (1 - probability_z0_given_treated) * probability1
  conditional <- lapply(
    specification$atom_shift_by_stratum,
    function(shift) {
      y <- specification$atoms + shift
      tan_x_t1 <- observed_dgp_standard_msm_arm(
        y, marginal_probability, marginal_pi, 1L, Lambda
      )
      tan_x_t0 <- observed_dgp_standard_msm_arm(
        y, probability0, marginal_pi, 0L, Lambda
      )
      tan_xz_t1 <-
        e * observed_dgp_standard_msm_arm(
          y, probability0, pi0, 1L, Lambda
        ) +
        (1 - e) * observed_dgp_standard_msm_arm(
          y, probability1, pi1, 1L, Lambda
        )
      tan_xz_t0 <-
        e * observed_dgp_standard_msm_arm(
          y, probability0, pi0, 0L, Lambda
        ) +
        (1 - e) * observed_dgp_standard_msm_arm(
          y, probability0, pi1, 0L, Lambda
        )
      c(
        x_mu1_lower = tan_x_t1[["lower"]],
        x_mu1_upper = tan_x_t1[["upper"]],
        x_mu0_lower = tan_x_t0[["lower"]],
        x_mu0_upper = tan_x_t0[["upper"]],
        xz_mu1_lower = tan_xz_t1[["lower"]],
        xz_mu1_upper = tan_xz_t1[["upper"]],
        xz_mu0_lower = tan_xz_t0[["lower"]],
        xz_mu0_upper = tan_xz_t0[["upper"]]
      )
    }
  )
  conditional <- do.call(rbind, conditional)
  output <- data.frame(
    estimator = c("tan_cal_x_saturated", "tan_cal_xz_saturated"),
    Lambda = Lambda,
    mu1_lower = c(
      mean(conditional[, "x_mu1_lower"]),
      mean(conditional[, "xz_mu1_lower"])
    ),
    mu1_upper = c(
      mean(conditional[, "x_mu1_upper"]),
      mean(conditional[, "xz_mu1_upper"])
    ),
    mu0_lower = c(
      mean(conditional[, "x_mu0_lower"]),
      mean(conditional[, "xz_mu0_lower"])
    ),
    mu0_upper = c(
      mean(conditional[, "x_mu0_upper"]),
      mean(conditional[, "xz_mu0_upper"])
    ),
    stringsAsFactors = FALSE
  )
  output$ate_lower <- output$mu1_lower - output$mu0_upper
  output$ate_upper <- output$mu1_upper - output$mu0_lower
  # Backward-compatible aliases used by the treated-arm calibration pilot.
  output$lower_target <- output$mu1_lower
  output$upper_target <- output$mu1_upper
  output
}

audit_prognostic_x_tan_threshold_tie <- function(
    oracle = build_prognostic_x_tan_threshold_tie_oracle()
) {
  specification <- oracle$lower$specification
  p0 <- specification$atom_probability_override
  p1 <- p0 * specification$gamma_atoms
  base_audit <- audit_selected_two_gap_nonunique_threshold(list(
    lower = oracle$lower,
    upper = oracle$upper
  ))
  data.frame(
    Lambda = specification$sensitivity_lambda,
    strata = specification$strata,
    candidates = specification$candidates,
    candidate_design = specification$candidate_design,
    probability_z0 = specification$probability_z0,
    pi0 = specification$pi0,
    pi1 = specification$pi1,
    minimum_xzt_probability = min(c(
      specification$probability_z0 * specification$pi0,
      specification$probability_z0 * (1 - specification$pi0),
      (1 - specification$probability_z0) * specification$pi1,
      (1 - specification$probability_z0) * (1 - specification$pi1)
    )) / specification$strata,
    shift_mean = mean(specification$atom_shift_by_stratum),
    shift_minimum = min(specification$atom_shift_by_stratum),
    shift_maximum = max(specification$atom_shift_by_stratum),
    lower_target = oracle$lower$target,
    upper_target = oracle$upper$target,
    lower_risk_range = diff(range(oracle$lower$population_risk)),
    upper_risk_range = diff(range(oracle$upper$population_risk)),
    maximum_normalization_residual = max(c(
      abs(oracle$lower$direct_normalization_residual),
      abs(oracle$upper$direct_normalization_residual)
    )),
    tan_z0_lower_cdf_residual = sum(p0[1:2]) - 0.2,
    tan_z0_upper_cdf_residual = sum(p0[1:4]) - 0.8,
    tan_z1_lower_cdf_residual = sum(p1[1:2]) - 0.2,
    tan_z1_upper_cdf_residual = sum(p1[1:4]) - 0.8,
    minimum_envelope_width = base_audit$minimum_envelope_width,
    minimum_envelope_branch_gap = base_audit$minimum_envelope_branch_gap,
    stringsAsFactors = FALSE
  )
}

generate_prognostic_x_tan_threshold_tie_data <- function(
    n,
    oracle = build_prognostic_x_tan_threshold_tie_oracle(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  index <- sample.int(
    nrow(oracle$categories), n, replace = TRUE,
    prob = oracle$categories$probability
  )
  sampled <- oracle$categories[index, , drop = FALSE]
  data.frame(
    Y = sampled$Y,
    T_obs = as.integer(sampled$T_obs),
    Z = as.integer(sampled$Z),
    X = as.integer(sampled$X),
    atom = as.integer(sampled$atom),
    category_index = index,
    stringsAsFactors = FALSE
  )
}

prognostic_x_tan_x_columns <- function(oracle) {
  strata <- oracle$lower$specification$strata
  paste0("X", seq_len(strata - 1L))
}

prognostic_x_tan_add_x_columns <- function(data, oracle) {
  strata <- oracle$lower$specification$strata
  for (stratum in seq_len(strata - 1L)) {
    data[[paste0("X", stratum)]] <- as.numeric(data$X == stratum)
  }
  data
}

# Generate the complete observed-data law used by the fully estimated study.
# The oracle columns are retained only so that the general estimator interface
# is complete; the validation driver requests the fully estimated RKHS layer
# and does not use these columns to fit the reported estimator.
generate_prognostic_x_tan_full_data <- function(
    n,
    oracle = build_prognostic_x_tan_threshold_tie_oracle(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  specification <- oracle$lower$specification
  x <- sample.int(specification$strata, n, replace = TRUE)
  z <- stats::rbinom(n, 1L, 1 - specification$probability_z0)
  local_pi <- ifelse(z == 0L, specification$pi0, specification$pi1)
  treatment <- stats::rbinom(n, 1L, local_pi)
  probability0 <- specification$atom_probability_override
  probability1_treated <- probability0 * specification$gamma_atoms
  atom <- integer(n)
  group <- interaction(treatment, z, drop = TRUE)
  for (level in levels(group)) {
    rows <- which(group == level)
    treated_z1 <- treatment[rows[1L]] == 1L && z[rows[1L]] == 1L
    atom[rows] <- sample.int(
      length(specification$atoms), length(rows), replace = TRUE,
      prob = if (treated_z1) probability1_treated else probability0
    )
  }
  output <- data.frame(
    Y = specification$atoms[atom] +
      specification$atom_shift_by_stratum[x],
    T_obs = treatment,
    Z = z,
    X = x,
    atom = atom,
    stringsAsFactors = FALSE
  )
  output$e_true <- specification$probability_z0
  output$pi0_true <- specification$pi0
  output$pi1_true <- specification$pi1
  output$gamma0_oracle <- 1
  output$gamma1_oracle <- specification$gamma_atoms[output$atom]
  output$prob_Z <- 1 - specification$probability_z0
  output$prob_T <- local_pi
  output <- prognostic_x_tan_add_x_columns(output, oracle)
  attr(output, "prognostic_x_tan_threshold_tie_oracle") <- oracle
  output
}

prognostic_x_tan_extreme_weights <- function(
    y,
    probability,
    lower_envelope,
    upper_envelope,
    direction = c("lower", "upper")
) {
  direction <- match.arg(direction)
  y <- as.numeric(y)
  probability <- as.numeric(probability)
  lower_envelope <- as.numeric(lower_envelope)
  upper_envelope <- as.numeric(upper_envelope)
  if (
    length(unique(c(
      length(y), length(probability), length(lower_envelope),
      length(upper_envelope)
    ))) != 1L ||
      any(!is.finite(c(
        y, probability, lower_envelope, upper_envelope
      ))) ||
      any(probability <= 0) ||
      any(lower_envelope > upper_envelope)
  ) {
    stop("Invalid finite-support extremizer inputs.")
  }
  probability <- probability / sum(probability)
  weight <- lower_envelope
  remaining <- 1 - sum(probability * weight)
  capacity <- probability * (upper_envelope - lower_envelope)
  if (remaining < -1e-10 || remaining > sum(capacity) + 1e-10) {
    stop("The envelope does not contain a normalized density ratio.")
  }
  ordering <- if (direction == "lower") order(y) else order(-y)
  for (index in ordering) {
    if (remaining <= 1e-14) break
    allocation <- min(remaining, capacity[index])
    weight[index] <- weight[index] + allocation / probability[index]
    remaining <- remaining - allocation
  }
  if (abs(sum(probability * weight) - 1) > 1e-10) {
    stop("The finite-support extremizer is not normalized.")
  }
  list(
    weight = weight,
    mean = sum(probability * weight * y),
    normalization_residual = sum(probability * weight) - 1
  )
}

# Exact sharp targets for the complete finite observed law.  Each treatment
# arm uses the original dual functional: its Z=0 observed mean plus the
# model-implied extremized counterfactual component.  The linear programs are
# solved directly on the six atoms, including fractional mass at an atom when
# an arm does not have an exact threshold gap.
compute_prognostic_x_tan_dual_population <- function(
    oracle = build_prognostic_x_tan_threshold_tie_oracle()
) {
  specification <- oracle$lower$specification
  Lambda <- specification$sensitivity_lambda
  probability0 <- specification$atom_probability_override
  probability1 <- probability0 * specification$gamma_atoms
  if (abs(sum(probability0) - 1) > 1e-10 ||
      abs(sum(probability1) - 1) > 1e-10) {
    stop("The prognostic-X conditional outcome laws are not normalized.")
  }
  arm_results <- matrix(
    NA_real_, nrow = 2L, ncol = 2L,
    dimnames = list(treatment = c("0", "1"), bound = c("lower", "upper"))
  )
  diagnostics <- list()
  row_index <- 0L
  for (treatment in 0:1) {
    gamma <- if (treatment == 1L) {
      specification$gamma_atoms
    } else {
      rep(1, length(probability0))
    }
    nuisance <- list(
      e = rep(specification$probability_z0, length(probability0)),
      pi0 = rep(specification$pi0, length(probability0)),
      pi1 = rep(specification$pi1, length(probability0))
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
    if (any(envelope$gamma_incompatible_original) ||
        any(envelope$L > envelope$U + 1e-12)) {
      stop("The prognostic-X oracle is incompatible for one treatment arm.")
    }
    r0 <- if (treatment == 1L) {
      specification$pi0
    } else {
      1 - specification$pi0
    }
    conditional_results <- matrix(
      NA_real_, nrow = specification$strata, ncol = 2L,
      dimnames = list(NULL, c("lower", "upper"))
    )
    for (stratum in seq_len(specification$strata)) {
      y <- specification$atoms +
        specification$atom_shift_by_stratum[stratum]
      observed_mean <- sum(probability0 * y)
      for (direction in c("lower", "upper")) {
        extremizer <- prognostic_x_tan_extreme_weights(
          y = y,
          probability = probability0,
          lower_envelope = envelope$L,
          upper_envelope = envelope$U,
          direction = direction
        )
        conditional_results[stratum, direction] <-
          r0 * observed_mean + (1 - r0) * extremizer$mean
        row_index <- row_index + 1L
        row <- data.frame(
          treatment = treatment,
          stratum = stratum,
          direction = direction,
          observed_mean_z0 = observed_mean,
          extremized_mean = extremizer$mean,
          arm_target = conditional_results[stratum, direction],
          normalization_residual = extremizer$normalization_residual,
          minimum_envelope_width = min(envelope$U - envelope$L),
          gamma_projection_rate = mean(envelope$gamma_projected),
          stringsAsFactors = FALSE
        )
        for (atom in seq_along(extremizer$weight)) {
          row[[paste0("lambda_atom", atom)]] <- extremizer$weight[atom]
        }
        diagnostics[[row_index]] <- row
      }
    }
    arm_results[as.character(treatment), ] <- colMeans(conditional_results)
  }
  endpoints <- data.frame(
    design = oracle$name,
    Lambda = Lambda,
    mu1_lower = arm_results["1", "lower"],
    mu1_upper = arm_results["1", "upper"],
    mu0_lower = arm_results["0", "lower"],
    mu0_upper = arm_results["0", "upper"],
    stringsAsFactors = FALSE
  )
  endpoints$ate_lower <- endpoints$mu1_lower - endpoints$mu0_upper
  endpoints$ate_upper <- endpoints$mu1_upper - endpoints$mu0_lower
  list(
    endpoints = endpoints,
    diagnostics = do.call(rbind, diagnostics),
    compatible = TRUE,
    threshold_tie_audit = audit_prognostic_x_tan_threshold_tie(oracle)
  )
}

prognostic_x_tan_design_matrix <- function(
    data,
    adjustment = c("x_saturated", "xz_saturated")
) {
  adjustment <- match.arg(adjustment)
  if (!all(c("X", "Z") %in% names(data))) {
    stop("The saturated Tan design requires X and Z.")
  }
  formula <- if (adjustment == "x_saturated") {
    ~ factor(X)
  } else {
    ~ factor(X) * factor(Z)
  }
  design <- stats::model.matrix(formula, data = data)
  storage.mode(design) <- "double"
  design
}

estimate_prognostic_x_tan_saturated <- function(
    data,
    Lambda = 4,
    probability_epsilon = 0.01
) {
  if (!exists("estimate_tan_cal_arm_endpoint", mode = "function")) {
    if (!file.exists("R/comparator_estimators.R")) {
      stop("Run this function from the repository root.")
    }
    source("R/comparator_estimators.R")
  }
  if (!all(c("Y", "T_obs", "X", "Z") %in% names(data))) {
    stop("The prognostic-X Tan fit requires Y, T_obs, X, and Z.")
  }
  arm <- as.numeric(data$T_obs == 1L)
  configurations <- c(
    tan_cal_x_saturated = "x_saturated",
    tan_cal_xz_saturated = "xz_saturated"
  )
  rows <- list()
  row_index <- 0L
  fits <- list()
  for (estimator in names(configurations)) {
    design <- prognostic_x_tan_design_matrix(
      data, configurations[[estimator]]
    )
    fits[[estimator]] <- list()
    for (direction in c("lower", "upper")) {
      fit <- estimate_tan_cal_arm_endpoint(
        data$Y,
        arm,
        design,
        Lambda,
        direction,
        probability_epsilon
      )
      fits[[estimator]][[direction]] <- fit
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        estimator = estimator,
        direction = direction,
        estimate = fit$estimate,
        standard_error = fit$standard_error,
        propensity_truncation_rate = fit$propensity$truncation_rate,
        quantile_solver = fit$quantile$solver,
        quantile_convergence = fit$quantile$convergence,
        stringsAsFactors = FALSE
      )
    }
  }
  list(summary = two_gap_bind_rows(rows), fits = fits)
}

estimate_prognostic_x_tan_full_bounds <- function(
    data,
    Lambda = 4,
    probability_epsilon = 0.01,
    confidence_level = 0.95
) {
  if (!exists("estimate_tan_cal_arm_endpoint", mode = "function")) {
    if (!file.exists("R/comparator_estimators.R")) {
      stop("Run this function from the repository root.")
    }
    source("R/comparator_estimators.R")
  }
  if (!all(c("Y", "T_obs", "X", "Z") %in% names(data))) {
    stop("The full prognostic-X Tan fit requires Y, T_obs, X, and Z.")
  }
  if (confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must lie strictly between zero and one.")
  }
  configurations <- c(
    tan_cal_x_saturated = "x_saturated",
    tan_cal_xz_saturated = "xz_saturated"
  )
  alpha <- 1 - confidence_level
  marginal_critical <- stats::qnorm(confidence_level)
  simultaneous_critical <- stats::qnorm(1 - alpha / 2)
  endpoint_rows <- list()
  interval_rows <- list()
  diagnostic_rows <- list()
  fits <- list()
  endpoint_index <- interval_index <- diagnostic_index <- 0L
  for (estimator in names(configurations)) {
    design <- prognostic_x_tan_design_matrix(
      data, configurations[[estimator]]
    )
    local <- list()
    for (treatment in 0:1) {
      arm <- as.numeric(data$T_obs == treatment)
      for (direction in c("lower", "upper")) {
        key <- paste0("t", treatment, "_", direction)
        local[[key]] <- estimate_tan_cal_arm_endpoint(
          data$Y, arm, design, Lambda, direction, probability_epsilon
        )
        diagnostic_index <- diagnostic_index + 1L
        diagnostic_rows[[diagnostic_index]] <- data.frame(
          estimator = estimator,
          treatment = treatment,
          direction = direction,
          propensity_truncation_rate =
            local[[key]]$propensity$truncation_rate,
          quantile_solver = if (
            !is.null(local[[key]]$quantile$solver)
          ) {
            local[[key]]$quantile$solver
          } else {
            "nelder_mead_legacy"
          },
          quantile_convergence = local[[key]]$quantile$convergence,
          stringsAsFactors = FALSE
        )
      }
    }
    estimates <- c(
      mu1_lower = local$t1_lower$estimate,
      mu1_upper = local$t1_upper$estimate,
      mu0_lower = local$t0_lower$estimate,
      mu0_upper = local$t0_upper$estimate
    )
    influence <- list(
      mu1_lower = local$t1_lower$influence_value,
      mu1_upper = local$t1_upper$influence_value,
      mu0_lower = local$t0_lower$influence_value,
      mu0_upper = local$t0_upper$influence_value
    )
    estimates <- c(
      estimates,
      ate_lower = estimates[["mu1_lower"]] - estimates[["mu0_upper"]],
      ate_upper = estimates[["mu1_upper"]] - estimates[["mu0_lower"]]
    )
    influence$ate_lower <- influence$mu1_lower - influence$mu0_upper
    influence$ate_upper <- influence$mu1_upper - influence$mu0_lower
    standard_error <- vapply(
      influence, stats::sd, numeric(1)
    ) / sqrt(nrow(data))
    for (endpoint in names(estimates)) {
      direction <- if (grepl("_lower$", endpoint)) "lower" else "upper"
      outward_sign <- if (direction == "lower") -1 else 1
      endpoint_index <- endpoint_index + 1L
      endpoint_rows[[endpoint_index]] <- data.frame(
        estimator = estimator,
        endpoint = endpoint,
        direction = direction,
        estimate = unname(estimates[[endpoint]]),
        standard_error = unname(standard_error[[endpoint]]),
        confidence_lower = estimates[[endpoint]] -
          simultaneous_critical * standard_error[[endpoint]],
        confidence_upper = estimates[[endpoint]] +
          simultaneous_critical * standard_error[[endpoint]],
        one_sided_outward_limit = estimates[[endpoint]] +
          outward_sign * marginal_critical * standard_error[[endpoint]],
        simultaneous_outward_limit = estimates[[endpoint]] +
          outward_sign * simultaneous_critical * standard_error[[endpoint]],
        stringsAsFactors = FALSE
      )
    }
    interval_index <- interval_index + 1L
    interval_rows[[interval_index]] <- data.frame(
      estimator = estimator,
      point_lower = estimates[["ate_lower"]],
      point_upper = estimates[["ate_upper"]],
      point_width = estimates[["ate_upper"]] - estimates[["ate_lower"]],
      confidence_lower = estimates[["ate_lower"]] -
        simultaneous_critical * standard_error[["ate_lower"]],
      confidence_upper = estimates[["ate_upper"]] +
        simultaneous_critical * standard_error[["ate_upper"]],
      stringsAsFactors = FALSE
    )
    interval_rows[[interval_index]]$confidence_width <- with(
      interval_rows[[interval_index]], confidence_upper - confidence_lower
    )
    fits[[estimator]] <- local
  }
  list(
    endpoint = two_gap_bind_rows(endpoint_rows),
    interval = two_gap_bind_rows(interval_rows),
    diagnostics = two_gap_bind_rows(diagnostic_rows),
    fits = fits,
    critical_values = c(
      marginal_one_sided = marginal_critical,
      simultaneous_two_endpoint = simultaneous_critical
    )
  )
}

simulate_prognostic_x_tan_multinomial <- function(
    oracle = build_prognostic_x_tan_threshold_tie_oracle(),
    sample_sizes = c(1500L, 3000L, 6000L),
    replications = 10000L,
    seed = 1600000L
) {
  sample_sizes <- as.integer(sample_sizes)
  replications <- as.integer(replications)
  if (!length(sample_sizes) || any(sample_sizes < 10L) || replications < 1L) {
    stop("Invalid prognostic-X Tan simulation inputs.")
  }
  specification <- oracle$lower$specification
  Lambda <- specification$sensitivity_lambda
  multiplier <- Lambda - 1 / Lambda
  population <- compute_prognostic_x_tan_population(oracle)
  collapsed <- stats::aggregate(
    probability ~ X + Z + T_obs + atom + Y,
    oracle$categories,
    sum
  )
  collapsed <- collapsed[order(
    collapsed$X, collapsed$Z, collapsed$T_obs, collapsed$atom
  ), , drop = FALSE]

  endpoint_fit <- function(count, adjustment, direction) {
    sign <- if (direction == "lower") -1 else 1
    tau <- if (direction == "lower") {
      1 / (Lambda + 1)
    } else {
      Lambda / (Lambda + 1)
    }
    group <- if (adjustment == "tan_cal_x_saturated") {
      collapsed$X
    } else {
      (collapsed$X - 1L) * 2L + collapsed$Z + 1L
    }
    phi <- numeric(nrow(collapsed))
    for (local_group in sort(unique(group))) {
      rows <- group == local_group
      treated <- rows & collapsed$T_obs == 1L
      group_count <- sum(count[rows])
      treated_count <- sum(count[treated])
      if (group_count == 0L || treated_count == 0L ||
          treated_count == group_count) {
        return(c(estimate = NA_real_, standard_error = NA_real_))
      }
      propensity <- treated_count / group_count
      odds <- (1 - propensity) / propensity
      atom_count <- tapply(
        count[treated], collapsed$atom[treated], sum
      )
      atom_count <- atom_count[
        as.character(seq_along(specification$atoms))
      ]
      atom_count[is.na(atom_count)] <- 0
      quantile_index <- which(
        cumsum(atom_count) / treated_count >= tau
      )[1L]
      local_x <- collapsed$X[which(rows)[1L]]
      q <- specification$atoms[quantile_index] +
        specification$atom_shift_by_stratum[local_x]
      residual <- collapsed$Y[rows] - q
      rho <- tau * pmax(residual, 0) +
        (1 - tau) * pmax(-residual, 0)
      transformed <- collapsed$Y[rows] + sign * multiplier * rho
      local_treated <- collapsed$T_obs[rows] == 1L
      outcome_mean <- sum(
        count[rows][local_treated] * transformed[local_treated]
      ) / treated_count
      arm <- as.numeric(collapsed$T_obs[rows] == 1L)
      phi[rows] <- arm / propensity * collapsed$Y[rows] +
        sign * multiplier * arm * odds * rho -
        (arm / propensity - 1) * outcome_mean
    }
    n <- sum(count)
    estimate <- sum(count * phi) / n
    variance <- (sum(count * phi^2) - n * estimate^2) / (n - 1)
    c(
      estimate = estimate,
      standard_error = sqrt(pmax(variance, 0) / n)
    )
  }

  output <- list()
  output_index <- 0L
  for (sample_index in seq_along(sample_sizes)) {
    n <- sample_sizes[sample_index]
    storage <- lapply(population$estimator, function(ignored) {
      list(
        lower_estimate = rep(NA_real_, replications),
        upper_estimate = rep(NA_real_, replications),
        lower_se = rep(NA_real_, replications),
        upper_se = rep(NA_real_, replications)
      )
    })
    names(storage) <- population$estimator
    set.seed(seed + n + sample_index)
    counts <- stats::rmultinom(
      replications, n, collapsed$probability
    )
    for (replication in seq_len(replications)) {
      for (estimator in population$estimator) {
        lower <- endpoint_fit(counts[, replication], estimator, "lower")
        upper <- endpoint_fit(counts[, replication], estimator, "upper")
        storage[[estimator]]$lower_estimate[replication] <-
          lower[["estimate"]]
        storage[[estimator]]$upper_estimate[replication] <-
          upper[["estimate"]]
        storage[[estimator]]$lower_se[replication] <-
          lower[["standard_error"]]
        storage[[estimator]]$upper_se[replication] <-
          upper[["standard_error"]]
      }
    }
    for (estimator in population$estimator) {
      target <- population[population$estimator == estimator, , drop = FALSE]
      values <- storage[[estimator]]
      complete <- Reduce(
        `&`, lapply(values, is.finite)
      )
      values <- lapply(values, function(value) value[complete])
      lower_z <- (values$lower_estimate - target$lower_target) /
        values$lower_se
      upper_z <- (values$upper_estimate - target$upper_target) /
        values$upper_se
      pointwise <- stats::qnorm(0.95)
      simultaneous <- stats::qnorm(0.975)
      lower_component_failure <- lower_z > simultaneous
      upper_component_failure <- upper_z < -simultaneous
      output_index <- output_index + 1L
      output[[output_index]] <- data.frame(
        sample_size = n,
        requested_replications = replications,
        successful_replications = sum(complete),
        estimator = estimator,
        lower_target = target$lower_target,
        upper_target = target$upper_target,
        lower_bias = mean(values$lower_estimate - target$lower_target),
        upper_bias = mean(values$upper_estimate - target$upper_target),
        lower_inward_bias = mean(values$lower_estimate - target$lower_target),
        upper_inward_bias = mean(target$upper_target - values$upper_estimate),
        lower_monte_carlo_sd = stats::sd(values$lower_estimate),
        upper_monte_carlo_sd = stats::sd(values$upper_estimate),
        lower_mean_standard_error = mean(values$lower_se),
        upper_mean_standard_error = mean(values$upper_se),
        lower_se_to_sd_ratio = mean(values$lower_se) /
          stats::sd(values$lower_estimate),
        upper_se_to_sd_ratio = mean(values$upper_se) /
          stats::sd(values$upper_estimate),
        lower_z_mean = mean(lower_z),
        lower_z_standard_deviation = stats::sd(lower_z),
        upper_z_mean = mean(upper_z),
        upper_z_standard_deviation = stats::sd(upper_z),
        lower_one_sided_coverage_95 = mean(lower_z <= pointwise),
        upper_one_sided_coverage_95 = mean(upper_z >= -pointwise),
        lower_component_coverage_97_5 = mean(!lower_component_failure),
        upper_component_coverage_97_5 = mean(!upper_component_failure),
        outward_set_coverage_pointwise_95 = mean(
          lower_z <= pointwise & upper_z >= -pointwise
        ),
        outward_set_coverage_bonferroni_95 = mean(
          !lower_component_failure & !upper_component_failure
        ),
        lower_only_failure_bonferroni_95 = mean(
          lower_component_failure & !upper_component_failure
        ),
        upper_only_failure_bonferroni_95 = mean(
          !lower_component_failure & upper_component_failure
        ),
        both_fail_bonferroni_95 = mean(
          lower_component_failure & upper_component_failure
        ),
        mean_point_set_width = mean(
          values$upper_estimate - values$lower_estimate
        ),
        mean_outward_set_width = mean(
          values$upper_estimate + simultaneous * values$upper_se -
            values$lower_estimate + simultaneous * values$lower_se
        ),
        endpoint_z_correlation = stats::cor(lower_z, upper_z),
        stringsAsFactors = FALSE
      )
    }
  }
  two_gap_bind_rows(output)
}
