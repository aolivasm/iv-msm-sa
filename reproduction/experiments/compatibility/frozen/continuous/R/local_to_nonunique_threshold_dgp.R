# Local-to-nonunique and fixed-limit IV laws for five-fold inference.
#
# For Lambda > 1, let alpha = 1 / (Lambda + 1) and tau = 1 - alpha.
# Conditional on T = 0, the outcome is a mixture of two separated uniforms.
# At sample size n,
#
#   P_n(Y in the lower component | T = 0, X, Z) =
#     threshold_level + sign * local_coefficient / sqrt(n),
#
# where threshold_level is alpha for a lower endpoint and tau for an upper
# endpoint.  A positive local_coefficient gives a triangular array: every
# finite-n law has a unique relevant quantile inside one component, while P_n
# approaches a law whose relevant quantile is the entire gap.  Setting
# local_coefficient = 0 gives that fixed limiting law at every sample size.
# The observed treatment and instrument probabilities are fixed, and gamma_0
# = gamma_1 = 1 exactly.
#
# The oracle below uses a finite quadrature representation of the continuous
# outcome law.  The quadrature grid is split at both local quantiles, so the
# population normalization of the true lower-bound extremizer is exact up to
# floating-point error.  A complete-binary threshold library lets each X
# stratum choose the lower- or upper-component local branch independently.

if (!exists("compute_eif_row_primitives", mode = "function")) {
  if (!file.exists("R/sharp_bounds_plugin.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/sharp_bounds_plugin.R")
}
if (!exists("multibranch_weighted_group_mean", mode = "function")) {
  if (!file.exists("R/multibranch_nonunique_threshold_pilot.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/multibranch_nonunique_threshold_pilot.R")
}

make_local_to_nonunique_threshold_spec <- function(
    sample_size,
    local_coefficient = 1,
    local_sign = c("above", "below"),
    endpoint_direction = c("lower", "upper"),
    sensitivity_lambda = 4,
    strata = 4L,
    probability_z0 = 0.5,
    pi0 = 0.7,
    pi1 = 0.9,
    lower_interval = c(-1, 0),
    upper_interval = c(1, 2),
    treated_outcome = 2,
    quadrature_nodes_per_segment = 12L
) {
  local_sign <- match.arg(local_sign)
  endpoint_direction <- match.arg(endpoint_direction)
  sample_size <- as.integer(sample_size)
  local_coefficient <- as.numeric(local_coefficient)
  sensitivity_lambda <- as.numeric(sensitivity_lambda)
  strata <- as.integer(strata)
  quadrature_nodes_per_segment <- as.integer(
    quadrature_nodes_per_segment
  )
  numeric_values <- c(
    sample_size, local_coefficient, sensitivity_lambda,
    probability_z0, pi0, pi1, lower_interval, upper_interval,
    treated_outcome, quadrature_nodes_per_segment
  )
  if (
    any(!is.finite(numeric_values)) || sample_size < 100L ||
      local_coefficient < 0 || sensitivity_lambda <= 1 ||
      strata < 2L || strata > 10L ||
      probability_z0 <= 0 || probability_z0 >= 1 ||
      pi0 <= 0 || pi0 >= 1 || pi1 <= 0 || pi1 >= 1 ||
      pi0 >= pi1 ||
      length(lower_interval) != 2L || length(upper_interval) != 2L ||
      lower_interval[1L] >= lower_interval[2L] ||
      upper_interval[1L] >= upper_interval[2L] ||
      lower_interval[2L] >= upper_interval[1L] ||
      quadrature_nodes_per_segment < 2L
  ) {
    stop("Invalid local-to-nonunique threshold specification.")
  }

  alpha <- 1 / (sensitivity_lambda + 1)
  tau <- 1 - alpha
  threshold_level <- if (endpoint_direction == "lower") alpha else tau
  local_shift <- local_coefficient / sqrt(sample_size)
  probability_below <- threshold_level - local_shift
  probability_above <- threshold_level + local_shift
  if (probability_below <= 0 || probability_above >= 1) {
    stop("The local coefficient puts a mixture probability outside (0,1).")
  }
  lower_width <- diff(lower_interval)
  upper_width <- diff(upper_interval)
  lower_branch <- lower_interval[1L] +
    threshold_level / probability_above * lower_width
  upper_branch <- upper_interval[1L] +
    (threshold_level - probability_below) /
    (1 - probability_below) * upper_width
  lower_probability <- if (local_sign == "above") {
    probability_above
  } else {
    probability_below
  }

  code <- as.matrix(expand.grid(
    rep(list(c(FALSE, TRUE)), strata),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ))
  threshold_matrix <- ifelse(code, upper_branch, lower_branch)
  dim(threshold_matrix) <- dim(code)
  dimnames(threshold_matrix) <- list(
    candidate = seq_len(nrow(code)), stratum = seq_len(strata)
  )
  exact_nonunique <- local_coefficient == 0
  if (exact_nonunique) {
    true_candidates <- seq_len(nrow(code))
    true_candidate <- true_candidates[1L]
  } else {
    true_uses_upper <- local_sign == "below"
    true_candidate <- which(
      rowSums(code == true_uses_upper) == strata
    )
    if (length(true_candidate) != 1L) {
      stop("Could not identify the unique population branch.")
    }
    true_candidates <- true_candidate
  }

  list(
    name = if (exact_nonunique) {
      "fixed_limit_nonunique_threshold"
    } else {
      "local_to_nonunique_threshold_triangular_array"
    },
    sample_size = sample_size,
    local_coefficient = local_coefficient,
    local_sign = local_sign,
    local_shift = local_shift,
    sensitivity_lambda = sensitivity_lambda,
    threshold_level = threshold_level,
    strata = strata,
    candidates = nrow(code),
    endpoint_direction = endpoint_direction,
    treatment = 0L,
    probability_z0 = probability_z0,
    pi0 = pi0,
    pi1 = pi1,
    lower_interval = as.numeric(lower_interval),
    upper_interval = as.numeric(upper_interval),
    treated_outcome = treated_outcome,
    quadrature_nodes_per_segment = quadrature_nodes_per_segment,
    lower_mixture_probability = lower_probability,
    probability_below = probability_below,
    probability_above = probability_above,
    threshold_values = c(
      lower_branch = lower_branch,
      upper_branch = upper_branch
    ),
    candidate_right_endpoint = code,
    threshold_matrix = threshold_matrix,
    exact_nonunique = exact_nonunique,
    true_candidate = true_candidate,
    true_candidates = true_candidates
  )
}

make_selected_local_to_nonunique_threshold_spec <- function(
    sample_size,
    local_coefficient = 1.5,
    sensitivity_lambda = 4,
    strata = 4L
) {
  make_local_to_nonunique_threshold_spec(
    sample_size = sample_size,
    local_coefficient = local_coefficient,
    local_sign = "below",
    endpoint_direction = "lower",
    sensitivity_lambda = sensitivity_lambda,
    strata = strata,
    probability_z0 = 0.5,
    pi0 = 0.7,
    pi1 = 0.9,
    lower_interval = c(-0.0001, 0),
    upper_interval = c(1, 1.0001),
    treated_outcome = 2,
    quadrature_nodes_per_segment = 12L
  )
}

make_fixed_limit_local_to_nonunique_threshold_spec <- function(
    sample_size,
    sensitivity_lambda = 4,
    strata = 4L
) {
  make_local_to_nonunique_threshold_spec(
    sample_size = sample_size,
    local_coefficient = 0,
    local_sign = "below",
    endpoint_direction = "lower",
    sensitivity_lambda = sensitivity_lambda,
    strata = strata,
    probability_z0 = 0.5,
    pi0 = 0.7,
    pi1 = 0.9,
    lower_interval = c(-0.0001, 0),
    upper_interval = c(1, 1.0001),
    treated_outcome = 2,
    quadrature_nodes_per_segment = 12L
  )
}

make_selected_local_to_nonunique_upper_threshold_spec <- function(
    sample_size,
    local_coefficient = 1,
    sensitivity_lambda = 4,
    strata = 4L
) {
  make_local_to_nonunique_threshold_spec(
    sample_size = sample_size,
    local_coefficient = local_coefficient,
    local_sign = "below",
    endpoint_direction = "upper",
    sensitivity_lambda = sensitivity_lambda,
    strata = strata,
    probability_z0 = 0.5,
    pi0 = 0.7,
    pi1 = 0.9,
    lower_interval = c(-0.0001, 0),
    upper_interval = c(1, 1.0001),
    treated_outcome = 2,
    quadrature_nodes_per_segment = 12L
  )
}

local_to_nonunique_x_columns <- function(specification) {
  paste0("X", seq_len(specification$strata - 1L))
}

local_to_nonunique_add_x_columns <- function(data, specification) {
  for (stratum in seq_len(specification$strata - 1L)) {
    data[[paste0("X", stratum)]] <- as.numeric(data$X == stratum)
  }
  data
}

generate_local_to_nonunique_threshold_data <- function(
    n,
    specification = make_selected_local_to_nonunique_threshold_spec(n),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L ||
      specification$sample_size != n) {
    stop("n must equal the positive sample size defining the local law.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  x <- sample.int(specification$strata, n, replace = TRUE)
  z <- stats::rbinom(n, 1L, 1 - specification$probability_z0)
  treatment_probability <- ifelse(
    z == 0L, specification$pi0, specification$pi1
  )
  treatment <- stats::rbinom(n, 1L, treatment_probability)
  y <- rep(specification$treated_outcome, n)
  untreated <- treatment == 0L
  lower_component <- rep(FALSE, n)
  lower_component[untreated] <- stats::rbinom(
    sum(untreated), 1L,
    specification$lower_mixture_probability
  ) == 1L
  lower_rows <- untreated & lower_component
  upper_rows <- untreated & !lower_component
  y[lower_rows] <- stats::runif(
    sum(lower_rows),
    specification$lower_interval[1L],
    specification$lower_interval[2L]
  )
  y[upper_rows] <- stats::runif(
    sum(upper_rows),
    specification$upper_interval[1L],
    specification$upper_interval[2L]
  )
  output <- data.frame(
    Y = y,
    T_obs = treatment,
    Z = z,
    X = x,
    lower_component = lower_component,
    e_true = specification$probability_z0,
    pi0_true = specification$pi0,
    pi1_true = specification$pi1,
    gamma0_oracle = 1,
    gamma1_oracle = 1,
    prob_Z = 1 - specification$probability_z0,
    prob_T = treatment_probability,
    stringsAsFactors = FALSE
  )
  output <- local_to_nonunique_add_x_columns(output, specification)
  attr(output, "local_to_nonunique_threshold_specification") <-
    specification
  output
}

local_to_nonunique_uniform_moment_below <- function(interval, threshold) {
  boundary <- min(max(threshold, interval[1L]), interval[2L])
  if (boundary <= interval[1L]) return(0)
  (boundary^2 - interval[1L]^2) / (2 * diff(interval))
}

local_to_nonunique_quantile <- function(level, specification) {
  p <- specification$lower_mixture_probability
  lower <- specification$lower_interval
  upper <- specification$upper_interval
  if (level <= p) {
    lower[1L] + level / p * diff(lower)
  } else {
    upper[1L] + (level - p) / (1 - p) * diff(upper)
  }
}

compute_local_to_nonunique_dual_population <- function(
    specification = make_selected_local_to_nonunique_threshold_spec(1500L)
) {
  Lambda <- specification$sensitivity_lambda
  p <- specification$lower_mixture_probability
  lower <- specification$lower_interval
  upper <- specification$upper_interval
  mean_y <- p * mean(lower) + (1 - p) * mean(upper)
  nuisance <- list(
    e = specification$probability_z0,
    pi0 = specification$pi0,
    pi1 = specification$pi1
  )
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = 1,
    t_value = 0L,
    Lambda0 = Lambda,
    Lambda1 = Lambda,
    envelope_smoothing_epsilon = 0,
    gamma_compatibility_action = "stop"
  )
  tau <- (envelope$U - 1) / envelope$D
  lower_level <- 1 - tau
  upper_level <- tau
  q_lower <- local_to_nonunique_quantile(lower_level, specification)
  q_upper <- local_to_nonunique_quantile(upper_level, specification)
  first_moment_below <- function(q) {
    p * local_to_nonunique_uniform_moment_below(lower, q) +
      (1 - p) * local_to_nonunique_uniform_moment_below(upper, q)
  }
  probability_below <- function(q) {
    p * pmin(pmax((q - lower[1L]) / diff(lower), 0), 1) +
      (1 - p) * pmin(pmax((q - upper[1L]) / diff(upper), 0), 1)
  }
  lower_weighted_mean <- envelope$L * mean_y +
    envelope$D * first_moment_below(q_lower)
  upper_weighted_mean <- envelope$U * mean_y -
    envelope$D * first_moment_below(q_upper)
  lower_normalization <- envelope$L +
    envelope$D * probability_below(q_lower)
  upper_normalization <- envelope$U -
    envelope$D * probability_below(q_upper)
  r0 <- 1 - specification$pi0
  mu0_lower <- r0 * mean_y + (1 - r0) * lower_weighted_mean
  mu0_upper <- r0 * mean_y + (1 - r0) * upper_weighted_mean
  endpoints <- data.frame(
    design = specification$name,
    sample_size = specification$sample_size,
    local_coefficient = specification$local_coefficient,
    Lambda = Lambda,
    mu1_lower = specification$treated_outcome,
    mu1_upper = specification$treated_outcome,
    mu0_lower = mu0_lower,
    mu0_upper = mu0_upper,
    ate_lower = specification$treated_outcome - mu0_upper,
    ate_upper = specification$treated_outcome - mu0_lower,
    stringsAsFactors = FALSE
  )
  diagnostics <- data.frame(
    treatment = 0L,
    lower_mixture_probability = p,
    limiting_threshold_level = specification$threshold_level,
    q_lower = q_lower,
    q_upper = q_upper,
    lower_normalization_residual = lower_normalization - 1,
    upper_normalization_residual = upper_normalization - 1,
    stringsAsFactors = FALSE
  )
  if (max(abs(c(
    diagnostics$lower_normalization_residual,
    diagnostics$upper_normalization_residual
  ))) > 1e-10) {
    stop("The local population extremizers are not normalized.")
  }
  list(endpoints = endpoints, diagnostics = diagnostics, compatible = TRUE)
}

compute_local_to_nonunique_tan_population <- function(
    specification = make_selected_local_to_nonunique_threshold_spec(1500L)
) {
  if (!exists("observed_dgp_standard_msm_arm", mode = "function")) {
    if (!file.exists("R/observed_data_dgp.R")) {
      stop("Run this function from the repository root.")
    }
    source("R/observed_data_dgp.R")
  }
  rule <- local_to_nonunique_outcome_rule(specification)
  Lambda <- specification$sensitivity_lambda
  arm_value <- function(treatment_probability, treatment) {
    if (treatment == 1L) {
      observed_dgp_standard_msm_arm(
        specification$treated_outcome, 1,
        treatment_probability, treatment, Lambda
      )
    } else {
      observed_dgp_standard_msm_arm(
        rule$y, rule$probability,
        treatment_probability, treatment, Lambda
      )
    }
  }
  marginal <- specification$probability_z0 * specification$pi0 +
    (1 - specification$probability_z0) * specification$pi1
  x0 <- arm_value(marginal, 0L)
  x1 <- arm_value(marginal, 1L)
  z00 <- arm_value(specification$pi0, 0L)
  z01 <- arm_value(specification$pi0, 1L)
  z10 <- arm_value(specification$pi1, 0L)
  z11 <- arm_value(specification$pi1, 1L)
  e <- specification$probability_z0
  make_row <- function(method, arm0, arm1) {
    data.frame(
      method = method,
      sample_size = specification$sample_size,
      Lambda = Lambda,
      mu1_lower = arm1[1L],
      mu1_upper = arm1[2L],
      mu0_lower = arm0[1L],
      mu0_upper = arm0[2L],
      ate_lower = arm1[1L] - arm0[2L],
      ate_upper = arm1[2L] - arm0[1L],
      stringsAsFactors = FALSE
    )
  }
  rbind(
    make_row("Tan MSM, adjust X", x0, x1),
    make_row(
      "Tan MSM, adjust X and Z",
      e * z00 + (1 - e) * z10,
      e * z01 + (1 - e) * z11
    )
  )
}

local_to_nonunique_midpoint_segment <- function(
    interval,
    total_probability,
    nodes
) {
  if (total_probability <= 0 || diff(interval) <= 0) {
    return(data.frame(y = numeric(), probability = numeric()))
  }
  width <- diff(interval) / nodes
  data.frame(
    y = interval[1L] + (seq_len(nodes) - 0.5) * width,
    probability = rep(total_probability / nodes, nodes)
  )
}

local_to_nonunique_outcome_rule <- function(specification) {
  p <- specification$lower_mixture_probability
  q_low <- specification$threshold_values[["lower_branch"]]
  q_high <- specification$threshold_values[["upper_branch"]]
  lower <- specification$lower_interval
  upper <- specification$upper_interval
  nodes <- specification$quadrature_nodes_per_segment
  segments <- list(
    local_to_nonunique_midpoint_segment(
      c(lower[1L], q_low),
      p * (q_low - lower[1L]) / diff(lower), nodes
    ),
    local_to_nonunique_midpoint_segment(
      c(q_low, lower[2L]),
      p * (lower[2L] - q_low) / diff(lower), nodes
    ),
    local_to_nonunique_midpoint_segment(
      c(upper[1L], q_high),
      (1 - p) * (q_high - upper[1L]) / diff(upper), nodes
    ),
    local_to_nonunique_midpoint_segment(
      c(q_high, upper[2L]),
      (1 - p) * (upper[2L] - q_high) / diff(upper), nodes
    )
  )
  output <- do.call(rbind, segments)
  rownames(output) <- NULL
  if (abs(sum(output$probability) - 1) > 1e-12) {
    stop("The local-to-nonunique outcome rule does not integrate to one.")
  }
  output
}

build_local_to_nonunique_threshold_oracle <- function(
    specification = make_local_to_nonunique_threshold_spec(1500L)
) {
  Lambda <- specification$sensitivity_lambda
  outcome_rule <- local_to_nonunique_outcome_rule(specification)
  outcome_rows <- rbind(
    data.frame(
      T_obs = 0L,
      Y = outcome_rule$y,
      outcome_probability = outcome_rule$probability
    ),
    data.frame(
      T_obs = 1L,
      Y = specification$treated_outcome,
      outcome_probability = 1
    )
  )
  categories <- merge(
    expand.grid(
      X = seq_len(specification$strata),
      Z = 0:1,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ),
    outcome_rows,
    by = NULL
  )
  local_pi <- ifelse(
    categories$Z == 0L, specification$pi0, specification$pi1
  )
  probability_z <- ifelse(
    categories$Z == 0L,
    specification$probability_z0,
    1 - specification$probability_z0
  )
  probability_t <- ifelse(
    categories$T_obs == 1L, local_pi, 1 - local_pi
  )
  categories$probability <-
    probability_z * probability_t * categories$outcome_probability /
    specification$strata
  if (abs(sum(categories$probability) - 1) > 1e-12) {
    stop("The local-to-nonunique observed law does not integrate to one.")
  }

  row_count <- nrow(categories)
  nuisance <- list(
    e = rep(specification$probability_z0, row_count),
    pi0 = rep(specification$pi0, row_count),
    pi1 = rep(specification$pi1, row_count)
  )
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = rep(1, row_count),
    t_value = 0L,
    Lambda0 = Lambda,
    Lambda1 = Lambda,
    envelope_smoothing_epsilon = 0,
    gamma_compatibility_action = "stop"
  )
  baseline <- categories$Z == 0L & categories$T_obs == 0L
  comparison <- categories$Z == 1L & categories$T_obs == 0L
  observed_arm <- as.numeric(categories$T_obs == 0L)
  instrument0 <- as.numeric(categories$Z == 0L)
  instrument1 <- 1 - instrument0
  r0 <- 1 - specification$pi0
  r1 <- 1 - specification$pi1
  x_levels <- seq_len(specification$strata)
  candidate_count <- specification$candidates
  score <- loss <- q_matrix <- matrix(
    NA_real_, nrow = row_count, ncol = candidate_count
  )
  population_direct_target <- normalization_residual <- numeric(
    candidate_count
  )

  group_mean <- function(value, subset) {
    multibranch_weighted_group_mean(
      value[subset], categories$probability[subset],
      categories$X[subset], x_levels
    )
  }

  for (candidate in seq_len(candidate_count)) {
    q <- specification$threshold_matrix[
      candidate, categories$X
    ]
    q_matrix[, candidate] <- q
    primitive <- compute_eif_row_primitives(
      data = categories,
      nuisance = nuisance,
      envelope = envelope,
      t_value = 0L,
      Lambda0 = Lambda,
      Lambda1 = Lambda,
      q_lower = q,
      q_upper = q
    )[[specification$endpoint_direction]]
    m0_x <- group_mean(categories$Y, baseline)
    Hbar0_x <- group_mean(primitive$H0, baseline)
    Hbar1_x <- group_mean(primitive$H1, comparison)
    D_r0_x <- group_mean(primitive$D_r0, baseline)
    D_r1_x <- group_mean(primitive$D_r1, baseline)
    linked <- compute_fixed_q_linked_components(
      m0 = m0_x[categories$X],
      Hbar0 = Hbar0_x[categories$X],
      Hbar1 = Hbar1_x[categories$X],
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
      (primitive$H0 - Hbar0_x[categories$X]) +
      instrument1 / (1 - specification$probability_z0) *
      observed_arm / r1 *
      (primitive$H1 - Hbar1_x[categories$X])
    score[, candidate] <- linked$M + treatment_correction +
      outcome_correction
    loss[, candidate] <- ifelse(
      baseline,
      envelope$D * check_loss(
        categories$Y, q,
        rep(specification$threshold_level, row_count)
      ),
      0
    )
    lambda <- if (specification$endpoint_direction == "lower") {
      ifelse(categories$Y <= q, envelope$U, envelope$L)
    } else {
      ifelse(categories$Y <= q, envelope$L, envelope$U)
    }
    lambda_x <- group_mean(lambda, baseline)
    theta_x <- group_mean(lambda * categories$Y, baseline)
    normalization_residual[candidate] <- max(abs(lambda_x - 1))
    population_direct_target[candidate] <- mean(
      r0 * m0_x + (1 - r0) * theta_x
    )
  }

  population_mean <- colSums(categories$probability * score)
  population_risk <- colSums(categories$probability * loss)
  centered <- sweep(score, 2L, population_mean, "-")
  population_covariance <- crossprod(
    sqrt(categories$probability) * centered
  )
  population_sd <- sqrt(diag(population_covariance))
  true_candidate <- specification$true_candidate
  expected_minimizers <- specification$true_candidates
  risk_tolerance <- 1e-12
  population_minimizers <- which(
    population_risk <= min(population_risk) + risk_tolerance
  )
  target <- population_direct_target[true_candidate]
  if (
    abs(normalization_residual[true_candidate]) > 1e-10 ||
      abs(population_mean[true_candidate] - target) > 1e-10 ||
      !setequal(population_minimizers, expected_minimizers) ||
      max(abs(
        population_direct_target[population_minimizers] - target
      )) > 1e-10 ||
      any(!is.finite(population_sd)) || any(population_sd <= 0)
  ) {
    stop("The local-to-nonunique oracle failed a population identity.")
  }
  nonminimizer_risk <- population_risk[
    !seq_along(population_risk) %in% population_minimizers
  ]
  local_risk_gap <- if (length(nonminimizer_risk)) {
    min(nonminimizer_risk) - min(population_risk)
  } else {
    0
  }

  list(
    specification = specification,
    categories = categories,
    q = q_matrix,
    score = score,
    selection_loss = loss,
    target = target,
    direct_target = population_direct_target,
    direct_normalization_residual = normalization_residual,
    population_mean = population_mean,
    population_risk = population_risk,
    population_sd = population_sd,
    population_covariance = population_covariance,
    population_correlation = stats::cov2cor(population_covariance),
    true_candidate = true_candidate,
    true_candidates = expected_minimizers,
    population_minimizers = population_minimizers,
    local_risk_gap = local_risk_gap,
    envelope = envelope
  )
}

simulate_local_to_nonunique_crossfit <- function(
    oracle,
    replications = 1000L,
    folds = 5L,
    seed = 1L,
    confidence_level = 0.95,
    batch_size = 1000L
) {
  replications <- as.integer(replications)
  folds <- as.integer(folds)
  batch_size <- as.integer(batch_size)
  n <- oracle$specification$sample_size
  if (replications < 1L || folds < 2L || n < folds || batch_size < 1L) {
    stop("Invalid local-to-nonunique simulation inputs.")
  }
  fold_sizes <- tabulate(
    rep(seq_len(folds), length.out = n), nbins = folds
  )
  probability <- oracle$categories$probability
  loss <- oracle$selection_loss
  score <- oracle$score
  score_squared <- score^2
  target <- oracle$target
  z_value <- stats::qnorm(1 - (1 - confidence_level) / 2)
  estimate <- standard_error <- z <- true_branch_rate <-
    distinct_selected <- numeric(replications)
  set.seed(as.integer(seed))
  completed <- 0L
  while (completed < replications) {
    current <- min(batch_size, replications - completed)
    counts <- if (length(unique(fold_sizes)) == 1L) {
      stats::rmultinom(
        current * folds, fold_sizes[1L], probability
      )
    } else {
      do.call(cbind, lapply(seq_len(current), function(replication) {
        do.call(cbind, lapply(
          fold_sizes,
          function(size) stats::rmultinom(1L, size, probability)
        ))
      }))
    }
    fold_loss <- array(
      t(counts) %*% loss,
      dim = c(folds, current, ncol(loss))
    )
    fold_score <- array(
      t(counts) %*% score,
      dim = c(folds, current, ncol(score))
    )
    fold_score_squared <- array(
      t(counts) %*% score_squared,
      dim = c(folds, current, ncol(score_squared))
    )
    total_loss <- apply(fold_loss, c(2L, 3L), sum)
    selected <- matrix(NA_integer_, nrow = folds, ncol = current)
    for (fold in seq_len(folds)) {
      local_loss <- total_loss - matrix(
        fold_loss[fold, , ], nrow = current
      )
      selected[fold, ] <- max.col(
        -local_loss / (n - fold_sizes[fold]), ties.method = "first"
      )
    }
    signal_sum <- signal_square_sum <- numeric(current)
    replication_index <- seq_len(current)
    for (fold in seq_len(folds)) {
      local_score <- matrix(
        fold_score[fold, , ], nrow = current
      )
      local_score_squared <- matrix(
        fold_score_squared[fold, , ], nrow = current
      )
      indices <- cbind(replication_index, selected[fold, ])
      signal_sum <- signal_sum + local_score[indices]
      signal_square_sum <- signal_square_sum +
        local_score_squared[indices]
    }
    local_estimate <- signal_sum / n
    signal_variance <- pmax(
      (signal_square_sum - n * local_estimate^2) / (n - 1), 0
    )
    local_standard_error <- sqrt(signal_variance / n)
    rows <- completed + seq_len(current)
    estimate[rows] <- local_estimate
    standard_error[rows] <- local_standard_error
    z[rows] <- (local_estimate - target) / local_standard_error
    selected_minimizer <- matrix(
      selected %in% oracle$true_candidates,
      nrow = nrow(selected), ncol = ncol(selected)
    )
    true_branch_rate[rows] <- colMeans(selected_minimizer)
    distinct_selected[rows] <- apply(selected, 2L, function(value) {
      length(unique(value))
    })
    completed <- completed + current
  }
  output <- data.frame(
    sample_size = n,
    replication = seq_len(replications),
    estimate = estimate,
    standard_error = standard_error,
    z = z,
    centered_coverage = abs(z) <= z_value,
    outward_coverage = multibranch_outward_covered(
      z, oracle$specification$endpoint_direction, confidence_level
    ),
    true_branch_selection_rate = true_branch_rate,
    distinct_selected = distinct_selected,
    stringsAsFactors = FALSE
  )
  attr(output, "target") <- target
  output
}

summarize_local_to_nonunique_crossfit <- function(detail) {
  data.frame(
    sample_size = detail$sample_size[1L],
    replications = nrow(detail),
    mean_estimate = mean(detail$estimate),
    bias = mean(detail$estimate) - attr(detail, "target"),
    empirical_standard_deviation = stats::sd(detail$estimate),
    mean_standard_error = mean(detail$standard_error),
    z_mean = mean(detail$z),
    z_standard_deviation = stats::sd(detail$z),
    centered_coverage = mean(detail$centered_coverage),
    outward_coverage = mean(detail$outward_coverage),
    mean_true_branch_selection_rate = mean(
      detail$true_branch_selection_rate
    ),
    mean_distinct_selected = mean(detail$distinct_selected),
    stringsAsFactors = FALSE
  )
}
