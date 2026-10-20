# Mixed-endpoint exact threshold-tie law for simultaneous ATE inference.
#
# One observed X subgroup uses configuration 9148, whose treated lower
# endpoint has an exact threshold interval.  The other subgroup uses
# configuration 10492, whose treated upper endpoint has an exact threshold
# interval.  The subgroup and the local binary stratum are both observed, so
# sharp population bounds separate across the resulting four X strata.  The
# control outcome is fixed at zero; hence the treated lower and upper bounds
# are also the lower and upper ATE bounds.

if (!exists("make_adverse_covariance_threshold_tie_spec", mode = "function")) {
  if (!file.exists("R/adverse_covariance_threshold_tie_dgp.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/adverse_covariance_threshold_tie_dgp.R")
}

make_mixed_endpoint_threshold_tie_spec <- function(
    lower_configuration_id = 9148L,
    upper_configuration_id = 11655L,
    lower_group_probability = 0.5
) {
  lower_group_probability <- as.numeric(lower_group_probability)
  if (
    length(lower_group_probability) != 1L ||
      !is.finite(lower_group_probability) ||
      lower_group_probability <= 0 || lower_group_probability >= 1
  ) {
    stop("lower_group_probability must lie strictly between zero and one.")
  }
  lower <- make_adverse_covariance_threshold_tie_spec(
    configuration_id = lower_configuration_id
  )
  upper <- make_adverse_covariance_threshold_tie_spec(
    configuration_id = upper_configuration_id
  )
  if (
    lower$endpoint_direction != "lower" ||
      upper$endpoint_direction != "upper" ||
      lower$treatment != 1L || upper$treatment != 1L ||
      lower$strata != 2L || upper$strata != 2L
  ) {
    stop("The mixed law requires binary treated-arm lower and upper ties.")
  }
  if (lower$sensitivity_lambda != upper$sensitivity_lambda) {
    stop("The mixed lower and upper laws must use the same Lambda.")
  }
  list(
    name = "mixed_lower_upper_ate_threshold_tie",
    sensitivity_lambda = lower$sensitivity_lambda,
    lower_group_probability = lower_group_probability,
    lower = lower,
    upper = upper,
    strata = 4L,
    lower_configuration_id = as.integer(lower_configuration_id),
    upper_configuration_id = as.integer(upper_configuration_id)
  )
}

build_mixed_endpoint_threshold_tie_oracle <- function(
    specification = make_mixed_endpoint_threshold_tie_spec()
) {
  list(
    lower = build_adverse_covariance_threshold_tie_oracle(
      specification$lower
    ),
    upper = build_adverse_covariance_threshold_tie_oracle(
      specification$upper
    )
  )
}

mixed_endpoint_x_columns <- function(specification) {
  paste0("X", seq_len(specification$strata - 1L))
}

mixed_endpoint_add_x_columns <- function(data, specification) {
  for (stratum in seq_len(specification$strata - 1L)) {
    data[[paste0("X", stratum)]] <- as.numeric(data$X == stratum)
  }
  data
}

generate_mixed_endpoint_threshold_tie_data <- function(
    n,
    specification = make_mixed_endpoint_threshold_tie_spec(),
    oracle = NULL,
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (is.null(oracle)) {
    oracle <- build_mixed_endpoint_threshold_tie_oracle(specification)
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  lower_group <- stats::rbinom(
    n, 1L, specification$lower_group_probability
  ) == 1L
  local_stratum <- sample.int(2L, n, replace = TRUE)
  group_specification <- ifelse(
    lower_group, "lower", "upper"
  )
  e <- ifelse(
    lower_group,
    specification$lower$probability_z0,
    specification$upper$probability_z0
  )
  pi0 <- ifelse(lower_group, specification$lower$pi0, specification$upper$pi0)
  pi1 <- ifelse(lower_group, specification$lower$pi1, specification$upper$pi1)
  z <- stats::rbinom(n, 1L, 1 - e)
  treatment_probability <- ifelse(z == 0L, pi0, pi1)
  treatment <- stats::rbinom(n, 1L, treatment_probability)
  atom <- rep(1L, n)
  outcome <- numeric(n)
  gamma1 <- rep(1, n)
  for (group_name in c("lower", "upper")) {
    group_rows <- which(group_specification == group_name & treatment == 1L)
    if (!length(group_rows)) next
    local_specification <- specification[[group_name]]
    local_oracle <- oracle[[group_name]]
    probability0 <- local_oracle$atom_probability
    probability1 <- probability0 * local_specification$gamma_atoms
    for (z_value in 0:1) {
      rows <- group_rows[z[group_rows] == z_value]
      if (!length(rows)) next
      atom[rows] <- sample.int(
        length(local_specification$atoms), length(rows), replace = TRUE,
        prob = if (z_value == 1L) probability1 else probability0
      )
      outcome[rows] <- local_specification$atoms[atom[rows]]
      gamma1[rows] <- local_specification$gamma_atoms[atom[rows]]
    }
  }
  output <- data.frame(
    Y = outcome,
    T_obs = treatment,
    Z = z,
    X = local_stratum + ifelse(lower_group, 0L, 2L),
    subgroup = ifelse(lower_group, "lower_tie", "upper_tie"),
    local_stratum = local_stratum,
    atom = atom,
    e_true = e,
    pi0_true = pi0,
    pi1_true = pi1,
    gamma0_oracle = 1,
    gamma1_oracle = gamma1,
    prob_Z = 1 - e,
    prob_T = treatment_probability,
    stringsAsFactors = FALSE
  )
  output <- mixed_endpoint_add_x_columns(output, specification)
  attr(output, "mixed_endpoint_threshold_tie_specification") <-
    specification
  output
}

compute_mixed_endpoint_dual_population <- function(
    specification = make_mixed_endpoint_threshold_tie_spec(),
    oracle = NULL
) {
  if (is.null(oracle)) {
    oracle <- build_mixed_endpoint_threshold_tie_oracle(specification)
  }
  lower <- compute_adverse_covariance_dual_population(
    specification$lower, oracle$lower
  )
  upper <- compute_adverse_covariance_dual_population(
    specification$upper, oracle$upper
  )
  weight <- specification$lower_group_probability
  endpoint_names <- c(
    "mu1_lower", "mu1_upper", "mu0_lower", "mu0_upper",
    "ate_lower", "ate_upper"
  )
  endpoints <- data.frame(
    design = specification$name,
    Lambda = specification$sensitivity_lambda,
    stringsAsFactors = FALSE
  )
  for (endpoint in endpoint_names) {
    endpoints[[endpoint]] <-
      weight * lower$endpoints[[endpoint]] +
      (1 - weight) * upper$endpoints[[endpoint]]
  }
  list(
    endpoints = endpoints,
    diagnostics = rbind(
      transform(lower$diagnostics, subgroup = "lower_tie"),
      transform(upper$diagnostics, subgroup = "upper_tie")
    ),
    compatible = isTRUE(lower$compatible) && isTRUE(upper$compatible),
    subgroup_population = list(lower = lower, upper = upper)
  )
}

compute_mixed_endpoint_tan_population <- function(
    specification = make_mixed_endpoint_threshold_tie_spec(),
    oracle = NULL
) {
  if (is.null(oracle)) {
    oracle <- build_mixed_endpoint_threshold_tie_oracle(specification)
  }
  lower <- compute_adverse_covariance_tan_population(
    specification$lower, oracle$lower
  )
  upper <- compute_adverse_covariance_tan_population(
    specification$upper, oracle$upper
  )
  weight <- specification$lower_group_probability
  endpoint_names <- c(
    "mu1_lower", "mu1_upper", "mu0_lower", "mu0_upper",
    "ate_lower", "ate_upper"
  )
  output <- lower
  for (endpoint in endpoint_names) {
    output[[endpoint]] <-
      weight * lower[[endpoint]] + (1 - weight) * upper[[endpoint]]
  }
  output
}

audit_mixed_endpoint_threshold_tie <- function(
    specification = make_mixed_endpoint_threshold_tie_spec(),
    oracle = NULL
) {
  if (is.null(oracle)) {
    oracle <- build_mixed_endpoint_threshold_tie_oracle(specification)
  }
  lower <- audit_adverse_covariance_threshold_tie(
    specification$lower, oracle$lower
  )
  upper <- audit_adverse_covariance_threshold_tie(
    specification$upper, oracle$upper
  )
  data.frame(
    lower_configuration_id = specification$lower_configuration_id,
    upper_configuration_id = specification$upper_configuration_id,
    lower_Lambda = lower$Lambda,
    upper_Lambda = upper$Lambda,
    lower_direction = lower$endpoint_direction,
    upper_direction = upper$endpoint_direction,
    strata = specification$strata,
    lower_group_probability = specification$lower_group_probability,
    minimum_xzt_probability = min(
      specification$lower_group_probability *
        lower$minimum_xzt_probability,
      (1 - specification$lower_group_probability) *
        upper$minimum_xzt_probability
    ),
    maximum_target_range = max(
      lower$maximum_target_range, upper$maximum_target_range
    ),
    maximum_risk_range = max(
      lower$maximum_risk_range, upper$maximum_risk_range
    ),
    maximum_normalization_residual = max(
      lower$maximum_normalization_residual,
      upper$maximum_normalization_residual
    ),
    stringsAsFactors = FALSE
  )
}
