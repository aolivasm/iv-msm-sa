# Complete observed-data IV law embedding an exact nonunique threshold.
#
# The observed law factors as
#
#   P(X) P(Z | X) P(T | X,Z) P(Y | T,X,Z).
#
# X is uniform on 16 strata and is exported through 15 reference-cell dummy
# variables.  Z is randomized and has a nonzero first stage for T.  Conditional
# outcome laws do not depend on Z, so gamma_t(y,x) = 1 exactly.
#
# For T = 0, 1 / (Lambda + 1) of the outcome distribution is uniform on a
# narrow lower interval and the remainder is uniform on a separated upper
# interval.  At Lambda = 4, the lower quantile level is 0.2, and therefore the
# entire gap is a set of population check-loss minimizers.  Sixteen candidate
# threshold functions choose two interior points of that gap according to a
# Hadamard code.  For T = 1, the outcome is constant.  This isolates inference
# for mu0_lower and, equivalently, ate_upper = mu1_upper - mu0_lower.
#
# This file provides the DGP, exact/quadrature population references, and the
# oracle fixed-threshold EIF library needed for a focused cross-fit-versus-
# online comparison.  It does not replace the general DGPs used elsewhere.

if (!exists("compute_eif_row_primitives", mode = "function")) {
  if (!file.exists("R/sharp_bounds_plugin.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/sharp_bounds_plugin.R")
}
if (!exists("observed_dgp_standard_msm_arm", mode = "function")) {
  if (!file.exists("R/observed_data_dgp.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/observed_data_dgp.R")
}

hadamard_observed_matrix <- function(order) {
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

make_hadamard_observed_iv_spec <- function(
    sensitivity_lambda = 4,
    strata = 16L,
    candidates = 16L,
    candidate_design = c("hadamard", "complete_binary"),
    probability_z0 = 0.5,
    pi0 = 0.7,
    pi1 = 0.9,
    lower_interval = c(-0.001, 0),
    upper_interval = c(1, 1.001),
    treated_outcome = 2,
    threshold_separation_fraction = 0.8,
    threshold_values = NULL,
    quadrature_nodes_per_interval = 24L
) {
  candidate_design <- match.arg(candidate_design)
  sensitivity_lambda <- as.numeric(sensitivity_lambda)
  strata <- as.integer(strata)
  candidates <- as.integer(candidates)
  quadrature_nodes_per_interval <- as.integer(
    quadrature_nodes_per_interval
  )
  lower_interval <- as.numeric(lower_interval)
  upper_interval <- as.numeric(upper_interval)
  numeric_values <- c(
    sensitivity_lambda, probability_z0, pi0, pi1,
    lower_interval, upper_interval, treated_outcome,
    threshold_separation_fraction
  )
  if (
    any(!is.finite(numeric_values)) || sensitivity_lambda <= 1 ||
      strata < 2L || candidates < 2L ||
      probability_z0 <= 0 || probability_z0 >= 1 ||
      pi0 <= 0 || pi0 >= 1 || pi1 <= 0 || pi1 >= 1 ||
      length(lower_interval) != 2L || length(upper_interval) != 2L ||
      lower_interval[1L] >= lower_interval[2L] ||
      upper_interval[1L] >= upper_interval[2L] ||
      lower_interval[2L] >= upper_interval[1L] ||
      threshold_separation_fraction <= 0 ||
      threshold_separation_fraction >= 1 ||
      quadrature_nodes_per_interval < 2L
  ) {
    stop("Invalid observed Hadamard IV specification.")
  }
  if (candidate_design == "hadamard" &&
      (bitwAnd(strata, strata - 1L) != 0L || candidates > strata)) {
    stop(
      "The Hadamard design requires power-of-two strata and candidates <= strata."
    )
  }
  if (candidate_design == "complete_binary" && strata > 10L) {
    stop("The complete-binary design is limited to at most ten strata.")
  }
  alpha <- 1 / (sensitivity_lambda + 1)
  lower_mean <- mean(lower_interval)
  upper_mean <- mean(upper_interval)
  loss_center <- (1 - alpha) * lower_mean + alpha * upper_mean
  maximum_half_separation <- min(
    loss_center - lower_interval[2L],
    upper_interval[1L] - loss_center
  )
  threshold_half_separation <-
    threshold_separation_fraction * maximum_half_separation
  if (is.null(threshold_values)) {
    threshold_values <- c(
      minus = loss_center - threshold_half_separation,
      plus = loss_center + threshold_half_separation
    )
  } else {
    threshold_values <- as.numeric(threshold_values)
    if (length(threshold_values) != 2L ||
        any(!is.finite(threshold_values)) ||
        threshold_values[1L] >= threshold_values[2L]) {
      stop("threshold_values must contain two increasing finite values.")
    }
    names(threshold_values) <- c("minus", "plus")
    threshold_half_separation <- diff(threshold_values) / 2
  }
  if (threshold_values[["minus"]] <= lower_interval[2L] ||
      threshold_values[["plus"]] >= upper_interval[1L]) {
    stop("Candidate thresholds must lie strictly inside the outcome gap.")
  }
  code <- if (candidate_design == "hadamard") {
    hadamard_observed_matrix(strata)[seq_len(candidates), , drop = FALSE]
  } else {
    value <- as.matrix(expand.grid(
      rep(list(c(-1, 1)), strata),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ))
    candidates <- nrow(value)
    value
  }
  threshold_matrix <- ifelse(
    code > 0,
    threshold_values[["plus"]],
    threshold_values[["minus"]]
  )
  dim(threshold_matrix) <- dim(code)
  dimnames(threshold_matrix) <- list(
    candidate = seq_len(candidates),
    stratum = seq_len(strata)
  )
  list(
    name = paste0(
      candidate_design, "_observed_iv_exact_threshold_tie"
    ),
    dimension = strata - 1L,
    sensitivity_lambda = sensitivity_lambda,
    threshold_level = alpha,
    strata = strata,
    candidates = candidates,
    candidate_design = candidate_design,
    probability_z0 = probability_z0,
    pi0 = pi0,
    pi1 = pi1,
    lower_interval = lower_interval,
    upper_interval = upper_interval,
    treated_outcome = treated_outcome,
    threshold_separation_fraction = threshold_separation_fraction,
    threshold_values = threshold_values,
    threshold_matrix = threshold_matrix,
    code = code,
    quadrature_nodes_per_interval = quadrature_nodes_per_interval,
    lower_mixture_probability = alpha,
    marginal_treatment_probability =
      probability_z0 * pi0 + (1 - probability_z0) * pi1,
    marginal_first_stage = pi1 - pi0
  )
}

hadamard_observed_x_columns <- function(specification) {
  paste0("X", seq_len(specification$strata - 1L))
}

hadamard_observed_add_x_columns <- function(data, specification) {
  for (stratum in seq_len(specification$strata - 1L)) {
    data[[paste0("X", stratum)]] <- as.numeric(data$X == stratum)
  }
  data
}

hadamard_observed_nuisances <- function(x, specification) {
  x <- as.integer(x)
  if (any(!x %in% seq_len(specification$strata))) {
    stop("X is outside the specified Hadamard strata.")
  }
  list(
    e = rep(specification$probability_z0, length(x)),
    pi0 = rep(specification$pi0, length(x)),
    pi1 = rep(specification$pi1, length(x))
  )
}

hadamard_observed_gamma <- function(y, x, treatment, specification) {
  rep(1, max(length(y), length(x), length(treatment)))
}

generate_hadamard_observed_iv_data <- function(
    n,
    specification = make_hadamard_observed_iv_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
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
  low <- rep(FALSE, n)
  low[untreated] <- stats::rbinom(
    sum(untreated), 1L, specification$lower_mixture_probability
  ) == 1L
  low_untreated <- untreated & low
  high_untreated <- untreated & !low
  y[low_untreated] <- stats::runif(
    sum(low_untreated),
    specification$lower_interval[1L],
    specification$lower_interval[2L]
  )
  y[high_untreated] <- stats::runif(
    sum(high_untreated),
    specification$upper_interval[1L],
    specification$upper_interval[2L]
  )
  output <- data.frame(
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
  output <- hadamard_observed_add_x_columns(output, specification)
  attr(output, "hadamard_observed_iv_specification") <- specification
  output
}

hadamard_observed_midpoint_rule <- function(interval, nodes) {
  width <- diff(interval) / nodes
  list(
    nodes = interval[1L] + (seq_len(nodes) - 0.5) * width,
    probability = rep(1 / nodes, nodes)
  )
}

hadamard_observed_untreated_rule <- function(specification) {
  low <- hadamard_observed_midpoint_rule(
    specification$lower_interval,
    specification$quadrature_nodes_per_interval
  )
  high <- hadamard_observed_midpoint_rule(
    specification$upper_interval,
    specification$quadrature_nodes_per_interval
  )
  list(
    y = c(low$nodes, high$nodes),
    probability = c(
      specification$lower_mixture_probability * low$probability,
      (1 - specification$lower_mixture_probability) * high$probability
    )
  )
}

hadamard_observed_weighted_quantile <- function(y, weight, level) {
  ordering <- order(y)
  y <- y[ordering]
  weight <- weight[ordering] / sum(weight)
  y[which(cumsum(weight) >= level)[1L]]
}

compute_hadamard_observed_dual_population <- function(
    specification = make_hadamard_observed_iv_spec()
) {
  Lambda <- specification$sensitivity_lambda
  rule0 <- hadamard_observed_untreated_rule(specification)
  arm_results <- matrix(
    NA_real_, nrow = 2L, ncol = 2L,
    dimnames = list(treatment = c("0", "1"), bound = c("lower", "upper"))
  )
  diagnostics <- list()
  for (treatment in 0:1) {
    y <- if (treatment == 0L) rule0$y else specification$treated_outcome
    probability <- if (treatment == 0L) rule0$probability else 1
    if (treatment == 1L) {
      arm_results["1", ] <- specification$treated_outcome
      diagnostics[[2L]] <- data.frame(
        treatment = 1L,
        Lbar = NA_real_,
        Ubar = NA_real_,
        lower_level = NA_real_,
        upper_level = NA_real_,
        q_lower = specification$treated_outcome,
        q_upper = specification$treated_outcome,
        lower_normalization_residual = 0,
        upper_normalization_residual = 0,
        stringsAsFactors = FALSE
      )
      next
    }
    nuisance <- list(
      e = rep(specification$probability_z0, length(y)),
      pi0 = rep(specification$pi0, length(y)),
      pi1 = rep(specification$pi1, length(y))
    )
    envelope <- compute_projected_envelopes(
      nuisance = nuisance,
      gamma = rep(1, length(y)),
      t_value = treatment,
      Lambda0 = Lambda,
      Lambda1 = Lambda,
      envelope_smoothing_epsilon = 0,
      gamma_compatibility_action = "stop"
    )
    Lbar <- sum(probability * envelope$L)
    Ubar <- sum(probability * envelope$U)
    tau <- (Ubar - 1) / (Ubar - Lbar)
    lower_level <- 1 - tau
    upper_level <- tau
    q_lower <- if (treatment == 0L) {
      mean(specification$threshold_values)
    } else {
      specification$treated_outcome
    }
    q_upper <- hadamard_observed_weighted_quantile(
      y, probability * envelope$D, upper_level
    )
    lambda_lower <- envelope$L + envelope$D * as.numeric(y <= q_lower)
    lambda_upper <- envelope$U - envelope$D * as.numeric(y <= q_upper)
    r0 <- if (treatment == 1L) {
      specification$pi0
    } else {
      1 - specification$pi0
    }
    observed_mean <- sum(probability * y)
    arm_results[as.character(treatment), "lower"] <-
      r0 * observed_mean + (1 - r0) * sum(
        probability * lambda_lower * y
      )
    arm_results[as.character(treatment), "upper"] <-
      r0 * observed_mean + (1 - r0) * sum(
        probability * lambda_upper * y
      )
    diagnostics[[treatment + 1L]] <- data.frame(
      treatment = treatment,
      Lbar = Lbar,
      Ubar = Ubar,
      lower_level = lower_level,
      upper_level = upper_level,
      q_lower = q_lower,
      q_upper = q_upper,
      lower_normalization_residual =
        sum(probability * lambda_lower) - 1,
      upper_normalization_residual =
        sum(probability * lambda_upper) - 1,
      stringsAsFactors = FALSE
    )
  }
  endpoints <- data.frame(
    design = specification$name,
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
    compatible = TRUE
  )
}

compute_hadamard_observed_tan_population <- function(
    specification = make_hadamard_observed_iv_spec()
) {
  Lambda <- specification$sensitivity_lambda
  rule0 <- hadamard_observed_untreated_rule(specification)
  arm_value <- function(treatment_probability, treatment) {
    if (treatment == 1L) {
      observed_dgp_standard_msm_arm(
        specification$treated_outcome, 1,
        treatment_probability, treatment, Lambda
      )
    } else {
      observed_dgp_standard_msm_arm(
        rule0$y, rule0$probability,
        treatment_probability, treatment, Lambda
      )
    }
  }
  marginal <- specification$marginal_treatment_probability
  x_arm0 <- arm_value(marginal, 0L)
  x_arm1 <- arm_value(marginal, 1L)
  z0_arm0 <- arm_value(specification$pi0, 0L)
  z0_arm1 <- arm_value(specification$pi0, 1L)
  z1_arm0 <- arm_value(specification$pi1, 0L)
  z1_arm1 <- arm_value(specification$pi1, 1L)
  e <- specification$probability_z0
  make_row <- function(method, arm0, arm1) {
    data.frame(
      method = method,
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
    make_row("Tan MSM, adjust X", x_arm0, x_arm1),
    make_row(
      "Tan MSM, adjust X and Z",
      e * z0_arm0 + (1 - e) * z1_arm0,
      e * z0_arm1 + (1 - e) * z1_arm1
    )
  )
}

hadamard_observed_group_mean <- function(value, probability, group, levels) {
  vapply(levels, function(level) {
    local <- group == level
    sum(probability[local] * value[local]) / sum(probability[local])
  }, numeric(1))
}

build_hadamard_observed_oracle <- function(
    specification = make_hadamard_observed_iv_spec()
) {
  Lambda <- specification$sensitivity_lambda
  rule0 <- hadamard_observed_untreated_rule(specification)
  outcome_rows <- rbind(
    data.frame(T_obs = 0L, Y = rule0$y, outcome_probability = rule0$probability),
    data.frame(
      T_obs = 1L, Y = specification$treated_outcome,
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
    (1 / specification$strata) * probability_z * probability_t *
    categories$outcome_probability
  nuisance <- hadamard_observed_nuisances(categories$X, specification)
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = rep(1, nrow(categories)),
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
  row_count <- nrow(categories)
  candidate_count <- specification$candidates
  score <- loss <- q_matrix <- matrix(
    NA_real_, nrow = row_count, ncol = candidate_count
  )
  component_names <- c("m0", "Hbar0", "Hbar1", "D_r0", "D_r1", "M", "C0", "C1")
  conditional_components <- array(
    NA_real_,
    dim = c(specification$strata, candidate_count, length(component_names)),
    dimnames = list(
      X = x_levels, candidate = seq_len(candidate_count),
      component = component_names
    )
  )
  target <- normalization_residual <- numeric(candidate_count)
  for (candidate in seq_len(candidate_count)) {
    q <- specification$threshold_matrix[candidate, categories$X]
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
    )$lower
    group_mean <- function(value, subset) {
      hadamard_observed_group_mean(
        value[subset], categories$probability[subset],
        categories$X[subset], x_levels
      )
    }
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
    lambda <- ifelse(categories$Y <= q, envelope$U, envelope$L)
    lambda_x <- group_mean(lambda, baseline)
    theta_x <- group_mean(lambda * categories$Y, baseline)
    normalization_residual[candidate] <- max(abs(lambda_x - 1))
    target[candidate] <- mean(r0 * m0_x + (1 - r0) * theta_x)
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
  }
  population_mean <- colSums(categories$probability * score)
  population_risk <- colSums(categories$probability * loss)
  centered <- sweep(score, 2L, population_mean, "-")
  covariance <- crossprod(sqrt(categories$probability) * centered)
  population_sd <- sqrt(diag(covariance))
  population_correlation <- covariance / tcrossprod(
    population_sd, population_sd
  )
  if (
    max(target) - min(target) > 1e-10 ||
      max(population_mean) - min(population_mean) > 1e-10 ||
      max(population_risk) - min(population_risk) > 1e-10 ||
      max(normalization_residual) > 1e-10 ||
      max(abs(population_mean - target)) > 1e-10
  ) {
    stop("The observed Hadamard oracle failed an exact population identity.")
  }
  list(
    specification = specification,
    categories = categories,
    q = q_matrix,
    score = score,
    selection_loss = loss,
    target = mean(target),
    direct_target = target,
    direct_normalization_residual = normalization_residual,
    population_mean = population_mean,
    population_risk = population_risk,
    population_sd = population_sd,
    population_correlation = population_correlation,
    conditional_components = conditional_components,
    envelope = envelope
  )
}

# Audit the optimizer-gradient condition for pathwise differentiability.
#
# Every column of oracle$score is the uncentered fixed-threshold EIF associated
# with one population-optimal threshold. Because all candidates have the same
# target, centering each column gives its canonical gradient. Under the
# active-set directional expansion used in the manuscript, the optimized
# functional is pathwise differentiable if and only if all of these gradients
# agree P-almost surely.
audit_hadamard_observed_optimizer_gradients <- function(
    oracle,
    equality_tolerance = 1e-10
) {
  required <- c(
    "categories", "score", "selection_loss", "population_mean",
    "population_risk", "direct_target"
  )
  if (
    !all(required %in% names(oracle)) ||
      !is.data.frame(oracle$categories) ||
      !"probability" %in% names(oracle$categories) ||
      !is.matrix(oracle$score) || ncol(oracle$score) < 2L ||
      length(equality_tolerance) != 1L ||
      !is.finite(equality_tolerance) || equality_tolerance < 0
  ) {
    stop("Invalid Hadamard optimizer-gradient audit inputs.")
  }
  probability <- oracle$categories$probability
  if (
    length(probability) != nrow(oracle$score) ||
      any(!is.finite(probability)) || any(probability < 0) ||
      abs(sum(probability) - 1) > 1e-10
  ) {
    stop("The optimizer-gradient audit requires valid row probabilities.")
  }
  gradient <- sweep(oracle$score, 2L, oracle$population_mean, "-")
  common_gradient <- rowMeans(gradient)
  varying_gradient <- sweep(gradient, 1L, common_gradient, "-")
  weighted_gradient <- sweep(gradient, 1L, sqrt(probability), "*")
  pairwise_distance <- as.matrix(stats::dist(t(weighted_gradient)))
  pairwise_values <- pairwise_distance[upper.tri(pairwise_distance)]

  total_variance <- mean(colSums(probability * gradient^2))
  common_variance <- sum(probability * common_gradient^2)
  varying_variance <- mean(colSums(probability * varying_gradient^2))
  if (!is.finite(total_variance) || total_variance <= 0) {
    stop("The optimizer-gradient process has zero or invalid variance.")
  }

  centered_loss <- sweep(
    oracle$selection_loss, 2L, oracle$population_risk, "-"
  )
  varying_loss <- sweep(centered_loss, 1L, rowMeans(centered_loss), "-")
  alignment <- vapply(seq_len(ncol(gradient)), function(candidate) {
    covariance <- sum(
      probability * varying_loss[, candidate] *
        varying_gradient[, candidate]
    )
    denominator <- sqrt(
      sum(probability * varying_loss[, candidate]^2) *
        sum(probability * varying_gradient[, candidate]^2)
    )
    if (denominator <= 0) NA_real_ else covariance / denominator
  }, numeric(1))
  weighted_varying <- sweep(
    varying_gradient, 1L, sqrt(probability), "*"
  )
  varying_rank <- qr(weighted_varying, tol = 1e-10)$rank
  equality_holds <- max(pairwise_values) <= equality_tolerance

  data.frame(
    candidates = ncol(gradient),
    maximum_candidate_target_range = diff(range(oracle$direct_target)),
    maximum_candidate_risk_range = diff(range(oracle$population_risk)),
    minimum_pairwise_gradient_l2_distance = min(pairwise_values),
    maximum_pairwise_gradient_l2_distance = max(pairwise_values),
    optimizer_varying_gradient_rank = varying_rank,
    common_gradient_variance_fraction = common_variance / total_variance,
    optimizer_varying_gradient_variance_fraction =
      varying_variance / total_variance,
    minimum_loss_gradient_alignment = min(alignment, na.rm = TRUE),
    maximum_loss_gradient_alignment = max(alignment, na.rm = TRUE),
    optimizer_gradient_equality_holds = equality_holds,
    pathwise_differentiable_by_optimizer_gradient_criterion = equality_holds,
    equality_tolerance = equality_tolerance,
    stringsAsFactors = FALSE
  )
}

evaluate_hadamard_observed_oracle_rows <- function(data, oracle) {
  specification <- oracle$specification
  n <- nrow(data)
  nuisance <- hadamard_observed_nuisances(data$X, specification)
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = rep(1, n),
    t_value = 0L,
    Lambda0 = specification$sensitivity_lambda,
    Lambda1 = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = 0,
    gamma_compatibility_action = "stop"
  )
  observed_arm <- as.numeric(data$T_obs == 0L)
  instrument0 <- as.numeric(data$Z == 0L)
  instrument1 <- 1 - instrument0
  baseline <- data$Z == 0L & data$T_obs == 0L
  r0 <- 1 - specification$pi0
  r1 <- 1 - specification$pi1
  score <- loss <- matrix(
    NA_real_, nrow = n, ncol = specification$candidates
  )
  for (candidate in seq_len(specification$candidates)) {
    q <- specification$threshold_matrix[candidate, data$X]
    primitive <- compute_eif_row_primitives(
      data = data,
      nuisance = nuisance,
      envelope = envelope,
      t_value = 0L,
      Lambda0 = specification$sensitivity_lambda,
      Lambda1 = specification$sensitivity_lambda,
      q_lower = q,
      q_upper = q
    )$lower
    component <- oracle$conditional_components[
      data$X, candidate, , drop = FALSE
    ]
    dim(component) <- c(n, dim(oracle$conditional_components)[3L])
    colnames(component) <- dimnames(
      oracle$conditional_components
    )$component
    treatment_correction <-
      instrument0 / specification$probability_z0 * component[, "C0"] *
      (observed_arm - r0) +
      instrument1 / (1 - specification$probability_z0) *
      component[, "C1"] * (observed_arm - r1)
    outcome_correction <-
      instrument0 / specification$probability_z0 * observed_arm / r0 *
      (primitive$H0 - component[, "Hbar0"]) +
      instrument1 / (1 - specification$probability_z0) *
      observed_arm / r1 *
      (primitive$H1 - component[, "Hbar1"])
    score[, candidate] <- component[, "M"] + treatment_correction +
      outcome_correction
    loss[, candidate] <- ifelse(
      baseline,
      envelope$D * check_loss(
        data$Y, q,
        rep(specification$threshold_level, n)
      ),
      0
    )
  }
  list(score = score, selection_loss = loss)
}
