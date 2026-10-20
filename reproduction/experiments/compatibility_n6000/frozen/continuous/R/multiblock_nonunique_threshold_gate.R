# Fixed multi-block nonunique-threshold gates.
#
# The first pair of simulators below is an idealized endpoint-process gate.  It
# is useful for screening covariance geometry but is not labeled as a complete
# ATE experiment.  The `mirrored_ate_*` construction is the genuine observed-
# data ATE law: four covariate blocks contain base and reflected outcomes in
# both treatment arms, so every term defining both ATE endpoints is retained.
# The DGP is fixed; neither its probabilities nor its block count depends on n.

if (!exists("build_selected_two_gap_nonunique_threshold_oracle", mode = "function")) {
  if (!file.exists("R/two_gap_nonunique_threshold_dgp.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/two_gap_nonunique_threshold_dgp.R")
}

build_multiblock_nonunique_threshold_gate <- function(
    block_pairs = 4L,
    local_strata = 16L,
    atom_layout = "joint_adverse"
) {
  block_pairs <- as.integer(block_pairs)
  local_strata <- as.integer(local_strata)
  if (length(block_pairs) != 1L || is.na(block_pairs) || block_pairs < 1L) {
    stop("block_pairs must be one positive integer.")
  }
  selected <- build_selected_two_gap_nonunique_threshold_oracle(
    strata = local_strata,
    atom_layout = atom_layout
  )
  unit <- selected$upper
  groups <- 2L * block_pairs
  structure(
    list(
      block_pairs = block_pairs,
      groups = groups,
      local_strata = local_strata,
      atom_layout = atom_layout,
      unit = unit,
      category_probability = rep(
        unit$categories$probability / groups,
        groups
      ),
      lower_groups = seq_len(block_pairs),
      upper_groups = block_pairs + seq_len(block_pairs),
      lower_target = -unit$target / 2,
      upper_target = unit$target / 2
    ),
    class = "multiblock_nonunique_threshold_gate"
  )
}

multiblock_gate_summarize <- function(
    sample_size,
    estimator,
    lower_estimate,
    upper_estimate,
    lower_standard_error,
    upper_standard_error,
    gate,
    initial_size = NA_integer_,
    chunks = NA_integer_
) {
  lower_z <- (lower_estimate - gate$lower_target) / lower_standard_error
  upper_z <- (upper_estimate - gate$upper_target) / upper_standard_error
  one_sided <- stats::qnorm(0.95)
  bonferroni <- stats::qnorm(0.975)
  data.frame(
    sample_size = as.integer(sample_size),
    block_pairs = gate$block_pairs,
    local_strata = gate$local_strata,
    estimator = estimator,
    replications = length(lower_estimate),
    initial_size = initial_size,
    chunks = chunks,
    lower_target = gate$lower_target,
    upper_target = gate$upper_target,
    lower_bias = mean(lower_estimate - gate$lower_target),
    upper_bias = mean(upper_estimate - gate$upper_target),
    lower_inward_bias = mean(lower_estimate - gate$lower_target),
    upper_inward_bias = mean(gate$upper_target - upper_estimate),
    lower_monte_carlo_sd = stats::sd(lower_estimate),
    upper_monte_carlo_sd = stats::sd(upper_estimate),
    lower_mean_standard_error = mean(lower_standard_error),
    upper_mean_standard_error = mean(upper_standard_error),
    lower_se_to_sd_ratio = mean(lower_standard_error) /
      stats::sd(lower_estimate),
    upper_se_to_sd_ratio = mean(upper_standard_error) /
      stats::sd(upper_estimate),
    lower_z_mean = mean(lower_z),
    upper_z_mean = mean(upper_z),
    lower_z_standard_deviation = stats::sd(lower_z),
    upper_z_standard_deviation = stats::sd(upper_z),
    lower_two_sided_coverage_95 = mean(abs(lower_z) <= bonferroni),
    upper_two_sided_coverage_95 = mean(abs(upper_z) <= bonferroni),
    lower_one_sided_coverage_95 = mean(lower_z <= one_sided),
    upper_one_sided_coverage_95 = mean(upper_z >= -one_sided),
    lower_component_coverage_97_5 = mean(lower_z <= bonferroni),
    upper_component_coverage_97_5 = mean(upper_z >= -bonferroni),
    outward_set_coverage_pointwise_95 = mean(
      lower_z <= one_sided & upper_z >= -one_sided
    ),
    outward_set_coverage_bonferroni_95 = mean(
      lower_z <= bonferroni & upper_z >= -bonferroni
    ),
    mean_point_set_width = mean(upper_estimate - lower_estimate),
    mean_outward_set_width = mean(
      upper_estimate + bonferroni * upper_standard_error -
        lower_estimate + bonferroni * lower_standard_error
    ),
    endpoint_z_correlation = stats::cor(lower_z, upper_z),
    stringsAsFactors = FALSE
  )
}

simulate_multiblock_gate_crossfit <- function(
    gate = build_multiblock_nonunique_threshold_gate(),
    sample_size = 3000L,
    replications = 1000L,
    folds = 5L,
    seed = 310901L,
    batch_size = 20L
) {
  sample_size <- as.integer(sample_size)
  replications <- as.integer(replications)
  folds <- as.integer(folds)
  batch_size <- as.integer(batch_size)
  if (
    sample_size < folds || sample_size %% folds != 0L ||
      replications < 1L || folds < 2L || batch_size < 1L
  ) {
    stop("Invalid multi-block cross-fit simulation inputs.")
  }
  unit <- gate$unit
  category_count <- nrow(unit$categories)
  candidate_count <- ncol(unit$score)
  lower_estimate <- upper_estimate <- numeric(replications)
  lower_standard_error <- upper_standard_error <- numeric(replications)
  set.seed(as.integer(seed))
  completed <- 0L
  while (completed < replications) {
    current <- min(batch_size, replications - completed)
    raw_counts <- stats::rmultinom(
      current * folds,
      sample_size / folds,
      gate$category_probability
    )
    counts <- array(
      raw_counts,
      dim = c(category_count, gate$groups, folds, current)
    )
    calculate_endpoint <- function(active_groups, sign) {
      signal_sum <- signal_square_sum <- numeric(current)
      replication_index <- seq_len(current)
      for (group in active_groups) {
        local_count <- matrix(
          counts[, group, , ],
          nrow = category_count,
          ncol = folds * current
        )
        fold_loss <- array(
          t(local_count) %*% unit$selection_loss,
          dim = c(folds, current, candidate_count)
        )
        fold_score <- array(
          t(local_count) %*% unit$score,
          dim = c(folds, current, candidate_count)
        )
        fold_score_squared <- array(
          t(local_count) %*% unit$score^2,
          dim = c(folds, current, candidate_count)
        )
        total_loss <- apply(fold_loss, c(2L, 3L), sum)
        for (fold in seq_len(folds)) {
          training_loss <- total_loss - matrix(
            fold_loss[fold, , ], nrow = current
          )
          selected <- max.col(-training_loss, ties.method = "first")
          score_matrix <- matrix(
            fold_score[fold, , ], nrow = current
          )
          square_matrix <- matrix(
            fold_score_squared[fold, , ], nrow = current
          )
          index <- cbind(replication_index, selected)
          signal_sum <- signal_sum + sign * score_matrix[index]
          signal_square_sum <- signal_square_sum + square_matrix[index]
        }
      }
      estimate <- signal_sum / sample_size
      standard_error <- sqrt(pmax(
        (signal_square_sum - sample_size * estimate^2) /
          (sample_size - 1),
        0
      ) / sample_size)
      list(estimate = estimate, standard_error = standard_error)
    }
    lower <- calculate_endpoint(gate$lower_groups, -1)
    upper <- calculate_endpoint(gate$upper_groups, 1)
    rows <- completed + seq_len(current)
    lower_estimate[rows] <- lower$estimate
    upper_estimate[rows] <- upper$estimate
    lower_standard_error[rows] <- lower$standard_error
    upper_standard_error[rows] <- upper$standard_error
    completed <- completed + current
  }
  multiblock_gate_summarize(
    sample_size = sample_size,
    estimator = paste0("oracle_", folds, "_fold_crossfit"),
    lower_estimate = lower_estimate,
    upper_estimate = upper_estimate,
    lower_standard_error = lower_standard_error,
    upper_standard_error = upper_standard_error,
    gate = gate
  )
}

simulate_multiblock_gate_online <- function(
    gate = build_multiblock_nonunique_threshold_gate(),
    sample_size = 3000L,
    replications = 1000L,
    initial_exponent = 0.7,
    chunks = 10L,
    seed = 410901L
) {
  sample_size <- as.integer(sample_size)
  replications <- as.integer(replications)
  chunks <- as.integer(chunks)
  initial_size <- max(200L, ceiling(sample_size^initial_exponent))
  evaluation_count <- sample_size - initial_size
  if (evaluation_count < chunks || replications < 1L || chunks < 1L) {
    stop("Invalid multi-block online simulation inputs.")
  }
  chunk_sizes <- tabulate(
    rep(seq_len(chunks), length.out = evaluation_count),
    nbins = chunks
  )
  unit <- gate$unit
  category_count <- nrow(unit$categories)
  lower_estimate <- upper_estimate <- numeric(replications)
  lower_standard_error <- upper_standard_error <- numeric(replications)
  set.seed(as.integer(seed))
  for (replication in seq_len(replications)) {
    initial <- array(
      stats::rmultinom(1L, initial_size, gate$category_probability),
      dim = c(category_count, gate$groups)
    )
    evaluation <- lapply(chunk_sizes, function(size) {
      array(
        stats::rmultinom(1L, size, gate$category_probability),
        dim = c(category_count, gate$groups)
      )
    })
    calculate_endpoint <- function(active_groups, sign) {
      training <- initial
      numerator <- denominator <- 0
      for (chunk in seq_len(chunks)) {
        chunk_signal_sum <- 0
        selected_second_moment <- 0
        for (group in active_groups) {
          selected <- which.min(as.vector(
            training[, group] %*% unit$selection_loss
          ))
          count <- evaluation[[chunk]][, group]
          chunk_signal_sum <- chunk_signal_sum + sign * sum(
            count * unit$score[, selected]
          )
          selected_second_moment <- selected_second_moment + sum(
            unit$categories$probability * unit$score[, selected]^2
          )
        }
        selected_second_moment <- selected_second_moment / gate$groups
        endpoint_target <- if (sign < 0) {
          gate$lower_target
        } else {
          gate$upper_target
        }
        selected_standard_deviation <- sqrt(pmax(
          selected_second_moment - endpoint_target^2,
          .Machine$double.eps
        ))
        numerator <- numerator + chunk_signal_sum /
          selected_standard_deviation
        denominator <- denominator + chunk_sizes[chunk] /
          selected_standard_deviation
        training <- training + evaluation[[chunk]]
      }
      estimate <- numerator / denominator
      standard_error <- (evaluation_count / denominator) /
        sqrt(evaluation_count)
      c(estimate = estimate, standard_error = standard_error)
    }
    lower <- calculate_endpoint(gate$lower_groups, -1)
    upper <- calculate_endpoint(gate$upper_groups, 1)
    lower_estimate[replication] <- lower[["estimate"]]
    upper_estimate[replication] <- upper[["estimate"]]
    lower_standard_error[replication] <- lower[["standard_error"]]
    upper_standard_error[replication] <- upper[["standard_error"]]
  }
  multiblock_gate_summarize(
    sample_size = sample_size,
    estimator = "oracle_online_predictable",
    lower_estimate = lower_estimate,
    upper_estimate = upper_estimate,
    lower_standard_error = lower_standard_error,
    upper_standard_error = upper_standard_error,
    gate = gate,
    initial_size = initial_size,
    chunks = chunks
  )
}

# A genuine observed-data ATE construction.  Four equally likely covariate
# blocks contain, respectively, a reflected treated arm, a base control arm,
# a base treated arm, and a reflected control arm.  For a reflected outcome
# -Y, its lower endpoint and score are the negatives of the base upper
# endpoint and score, and conversely.  Thus both terms defining each ATE
# endpoint carry the adverse upper-threshold selection error in the inward
# direction.  The other two terms are the corresponding base lower endpoints.
build_mirrored_ate_nonunique_threshold_gate <- function(
    block_quartets = 1L,
    local_strata = 16L,
    atom_layout = "joint_adverse",
    outcome_atoms = c(0, 0.01, 0.02, 0.03, 1, 1.01)
) {
  block_quartets <- as.integer(block_quartets)
  local_strata <- as.integer(local_strata)
  if (
    length(block_quartets) != 1L || is.na(block_quartets) ||
      block_quartets < 1L
  ) {
    stop("block_quartets must be one positive integer.")
  }
  selected <- build_selected_two_gap_nonunique_threshold_oracle(
    strata = local_strata,
    atom_layout = atom_layout
  )
  outcome_atoms <- as.numeric(outcome_atoms)
  if (
    length(outcome_atoms) != length(selected$upper$specification$atoms) ||
      any(!is.finite(outcome_atoms)) || any(diff(outcome_atoms) <= 0)
  ) {
    stop("outcome_atoms must be a strictly increasing six-vector.")
  }
  lower_specification <- selected$lower$specification
  upper_specification <- selected$upper$specification
  lower_specification$atoms <- outcome_atoms
  upper_specification$atoms <- outcome_atoms
  selected$lower <- build_multibranch_nonunique_threshold_oracle(
    lower_specification
  )
  selected$upper <- build_multibranch_nonunique_threshold_oracle(
    upper_specification
  )
  groups <- 4L * block_quartets
  group_type <- rep(seq_len(4L), each = block_quartets)
  structure(
    list(
      block_quartets = block_quartets,
      groups = groups,
      group_type = group_type,
      local_strata = local_strata,
      atom_layout = atom_layout,
      outcome_atoms = outcome_atoms,
      lower = selected$lower,
      upper = selected$upper,
      category_probability = rep(
        selected$upper$categories$probability / groups,
        groups
      ),
      lower_target = (selected$lower$target - selected$upper$target) / 2,
      upper_target = (selected$upper$target - selected$lower$target) / 2
    ),
    class = "mirrored_ate_nonunique_threshold_gate"
  )
}

mirrored_ate_gate_x_columns <- function(gate) {
  c(
    paste0("block", seq_len(gate$groups - 1L)),
    paste0("local_x", seq_len(gate$local_strata - 1L))
  )
}

mirrored_ate_gate_cell_columns <- function(gate) {
  paste0("cell", seq_len(gate$groups * gate$local_strata - 1L))
}

generate_mirrored_ate_nonunique_threshold_data <- function(
    n,
    gate = build_mirrored_ate_nonunique_threshold_gate(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  specification <- gate$upper$specification
  block <- sample.int(gate$groups, n, replace = TRUE)
  block_type <- gate$group_type[block]
  local_x <- sample.int(gate$local_strata, n, replace = TRUE)
  z <- stats::rbinom(n, 1L, 1 - specification$probability_z0)
  active_treatment <- ifelse(block_type %in% c(1L, 3L), 1L, 0L)
  active_probability <- ifelse(z == 0L, specification$pi0, specification$pi1)
  treatment_probability <- ifelse(
    active_treatment == 1L,
    active_probability,
    1 - active_probability
  )
  treatment <- stats::rbinom(n, 1L, treatment_probability)
  active <- treatment == active_treatment
  atom <- rep(1L, n)
  probability0 <- specification$atom_probability_override
  probability1 <- probability0 * specification$gamma_atoms
  for (z_value in 0:1) {
    rows <- which(active & z == z_value)
    if (length(rows)) {
      atom[rows] <- sample.int(
        length(specification$atoms), length(rows), replace = TRUE,
        prob = if (z_value == 1L) probability1 else probability0
      )
    }
  }
  reflected <- block_type %in% c(1L, 4L)
  outcome <- rep(0, n)
  outcome[active] <- specification$atoms[atom[active]]
  outcome[active & reflected] <- -outcome[active & reflected]
  gamma_active <- specification$gamma_atoms[atom]
  output <- data.frame(
    Y = outcome,
    T_obs = treatment,
    Z = z,
    X = (block - 1L) * gate$local_strata + local_x,
    block = block,
    block_type = block_type,
    local_X = local_x,
    atom = atom,
    e_true = specification$probability_z0,
    pi0_true = ifelse(
      active_treatment == 1L,
      specification$pi0,
      1 - specification$pi0
    ),
    pi1_true = ifelse(
      active_treatment == 1L,
      specification$pi1,
      1 - specification$pi1
    ),
    gamma0_oracle = ifelse(treatment == 0L & active, gamma_active, 1),
    gamma1_oracle = ifelse(treatment == 1L & active, gamma_active, 1),
    prob_Z = 1 - specification$probability_z0,
    prob_T = treatment_probability,
    stringsAsFactors = FALSE
  )
  for (index in seq_len(gate$groups - 1L)) {
    output[[paste0("block", index)]] <- as.numeric(block == index)
  }
  for (index in seq_len(gate$local_strata - 1L)) {
    output[[paste0("local_x", index)]] <- as.numeric(local_x == index)
  }
  cell <- (block - 1L) * gate$local_strata + local_x
  for (index in seq_len(gate$groups * gate$local_strata - 1L)) {
    output[[paste0("cell", index)]] <- as.numeric(cell == index)
  }
  attr(output, "mirrored_ate_gate") <- gate
  output
}

make_mirrored_ate_parametric_nuisance_provider_factory <- function(
    gate = build_mirrored_ate_nonunique_threshold_gate(),
    pseudocount = 0.5
) {
  pseudocount <- as.numeric(pseudocount)
  if (length(pseudocount) != 1L || !is.finite(pseudocount) ||
      pseudocount <= 0) {
    stop("pseudocount must be one positive number.")
  }
  function(
      fold, train, x_cols, seed, positivity_epsilon = 0.05, ...
  ) {
    force(fold)
    force(x_cols)
    force(seed)
    probability_z0 <- (sum(train$Z == 0L) + pseudocount) /
      (nrow(train) + 2 * pseudocount)
    treatment_probability <- array(
      NA_real_, dim = c(4L, 2L),
      dimnames = list(paste0("type", 1:4), c("z0", "z1"))
    )
    for (type in 1:4) {
      for (z_value in 0:1) {
        rows <- train$block_type == type & train$Z == z_value
        treatment_probability[type, z_value + 1L] <-
          (sum(train$T_obs[rows] == 1L) + pseudocount) /
          (sum(rows) + 2 * pseudocount)
      }
    }
    predict_nuisance <- function(newdata) {
      type <- as.integer(newdata$block_type)
      list(
        e = rep(
          min(max(probability_z0, positivity_epsilon),
              1 - positivity_epsilon),
          nrow(newdata)
        ),
        pi0 = pmin(
          pmax(treatment_probability[cbind(type, 1L)], positivity_epsilon),
          1 - positivity_epsilon
        ),
        pi1 = pmin(
          pmax(treatment_probability[cbind(type, 2L)], positivity_epsilon),
          1 - positivity_epsilon
        )
      )
    }
    list(
      predict = predict_nuisance,
      learner = "correct_saturated_parametric",
      diagnostics = data.frame(
        component = c("e", "pi0", "pi1"),
        learner = "correct_saturated_parametric",
        converged = TRUE,
        model_rank = c(1L, 8L, 8L),
        effective_df = c(1L, 8L, 8L),
        n_training = c(nrow(train), sum(train$Z == 0L), sum(train$Z == 1L)),
        stringsAsFactors = FALSE
      )
    )
  }
}

make_mirrored_ate_parametric_gamma_provider_factory <- function(
    gate = build_mirrored_ate_nonunique_threshold_gate(),
    pseudocount = 0.5
) {
  pseudocount <- as.numeric(pseudocount)
  if (length(pseudocount) != 1L || !is.finite(pseudocount) ||
      pseudocount <= 0) {
    stop("pseudocount must be one positive number.")
  }
  function(fold, treatment, train, nuisance_train, x_cols, seed, ...) {
    force(fold)
    force(nuisance_train)
    force(x_cols)
    force(seed)
    lookup <- vector("list", gate$groups)
    for (group in seq_len(gate$groups)) {
      type <- gate$group_type[group]
      active_treatment <- if (type %in% c(1L, 3L)) 1L else 0L
      reflected <- type %in% c(1L, 4L)
      if (treatment != active_treatment) {
        lookup[[group]] <- list(active = FALSE)
        next
      }
      support <- if (reflected) -gate$outcome_atoms else gate$outcome_atoms
      probability <- matrix(
        NA_real_, nrow = length(support), ncol = 2L,
        dimnames = list(NULL, c("z0", "z1"))
      )
      for (z_value in 0:1) {
        rows <- train$block == group & train$T_obs == treatment &
          train$Z == z_value
        outcome_index <- match(train$Y[rows], support)
        count <- tabulate(outcome_index, nbins = length(support))
        probability[, z_value + 1L] <- (count + pseudocount) /
          (sum(count) + pseudocount * length(support))
      }
      lookup[[group]] <- list(
        active = TRUE,
        support = support,
        gamma = probability[, "z1"] / probability[, "z0"]
      )
    }
    predict_gamma <- function(newdata) {
      value <- rep(1, nrow(newdata))
      for (group in seq_len(gate$groups)) {
        information <- lookup[[group]]
        if (!isTRUE(information$active)) next
        rows <- which(newdata$block == group)
        if (!length(rows)) next
        index <- match(newdata$Y[rows], information$support)
        supported <- !is.na(index)
        value[rows[supported]] <- information$gamma[index[supported]]
      }
      list(raw = value, positive = value, normalized = value)
    }
    list(
      predict = predict_gamma,
      gamma_estimator = "correct_saturated_parametric",
      selected_lambda = NA_real_,
      selected_bandwidth_multiplier = NA_real_,
      base_bandwidth = NA_real_,
      selected_bandwidth = NA_real_,
      cv_table = data.frame(),
      normalization_constant = 1,
      training_negative_rate = 0,
      rkhs_backend = "correct_saturated_parametric",
      nystrom_rank_requested = NA_integer_,
      nystrom_landmark_count = NA_integer_,
      nystrom_effective_rank = NA_integer_,
      nystrom_min_retained_eigenvalue = NA_real_,
      nystrom_max_eigenvalue = NA_real_
    )
  }
}

# Correctly specified direct density-ratio model for the mirrored gate.  The
# active treatment arms share one six-atom ratio after reflected outcomes are
# mapped back to the base support.  We therefore estimate that common ratio by
# a pooled exponential tilt rather than by separately estimating f1 and f0 in
# every block and dividing them.  The final ratio is normalized against the
# empirical Z=0 distribution, so its training-sample baseline mean is exactly
# one.  A ridge of order 1 / n is used only for numerical stabilization.
make_mirrored_ate_direct_ratio_gamma_provider_factory <- function(
    gate = build_mirrored_ate_nonunique_threshold_gate(),
    pseudocount = 0.5,
    ridge_scale = 1
) {
  pseudocount <- as.numeric(pseudocount)
  ridge_scale <- as.numeric(ridge_scale)
  if (length(pseudocount) != 1L || !is.finite(pseudocount) ||
      pseudocount <= 0) {
    stop("pseudocount must be one positive number.")
  }
  if (length(ridge_scale) != 1L || !is.finite(ridge_scale) ||
      ridge_scale < 0) {
    stop("ridge_scale must be one nonnegative number.")
  }
  function(fold, treatment, train, nuisance_train, x_cols, seed, ...) {
    force(fold)
    force(nuisance_train)
    force(x_cols)
    force(seed)
    treatment <- as.integer(treatment)
    if (!treatment %in% 0:1) stop("treatment must be zero or one.")

    active_treatment <- ifelse(
      train$block_type %in% c(1L, 3L), 1L, 0L
    )
    reflected <- train$block_type %in% c(1L, 4L)
    canonical_outcome <- ifelse(reflected, -train$Y, train$Y)
    atom_index <- match(canonical_outcome, gate$outcome_atoms)
    active <- train$T_obs == active_treatment & !is.na(atom_index)
    if (!any(active & train$Z == 0L) || !any(active & train$Z == 1L)) {
      stop("Both instrument strata are required for direct ratio fitting.")
    }

    atom_count <- length(gate$outcome_atoms)
    count0 <- tabulate(
      atom_index[active & train$Z == 0L], nbins = atom_count
    )
    count1 <- tabulate(
      atom_index[active & train$Z == 1L], nbins = atom_count
    )
    probability0 <- (count0 + pseudocount) /
      (sum(count0) + pseudocount * atom_count)
    probability1 <- (count1 + pseudocount) /
      (sum(count1) + pseudocount * atom_count)
    ridge <- ridge_scale / sum(active)

    objective <- function(theta) {
      eta <- c(0, theta)
      maximum <- max(eta)
      log_normalizer <- maximum + log(sum(
        probability0 * exp(eta - maximum)
      ))
      log_normalizer - sum(probability1 * eta) +
        0.5 * ridge * sum(theta^2)
    }
    gradient <- function(theta) {
      eta <- c(0, theta)
      tilted <- probability0 * exp(eta - max(eta))
      tilted <- tilted / sum(tilted)
      tilted[-1L] - probability1[-1L] + ridge * theta
    }
    initial_eta <- log(probability1) - log(probability0)
    initial_theta <- initial_eta[-1L] - initial_eta[1L]
    optimization <- stats::optim(
      par = initial_theta,
      fn = objective,
      gr = gradient,
      method = "BFGS",
      control = list(maxit = 500L, reltol = 1e-12)
    )
    if (optimization$convergence != 0L ||
        any(!is.finite(optimization$par))) {
      stop(
        "Direct exponential-tilt gamma fit did not converge; code ",
        optimization$convergence, "."
      )
    }

    eta <- c(0, optimization$par)
    raw_ratio <- exp(eta - max(eta))
    empirical_probability0 <- count0 / sum(count0)
    normalization_constant <- sum(empirical_probability0 * raw_ratio)
    ratio_by_atom <- raw_ratio / normalization_constant

    predict_gamma <- function(newdata) {
      value <- rep(1, nrow(newdata))
      new_active_treatment <- ifelse(
        newdata$block_type %in% c(1L, 3L), 1L, 0L
      )
      relevant <- new_active_treatment == treatment
      new_reflected <- newdata$block_type %in% c(1L, 4L)
      new_canonical_outcome <- ifelse(
        new_reflected, -newdata$Y, newdata$Y
      )
      new_atom_index <- match(new_canonical_outcome, gate$outcome_atoms)
      supported <- relevant & !is.na(new_atom_index)
      value[supported] <- ratio_by_atom[new_atom_index[supported]]
      list(raw = value, positive = value, normalized = value)
    }
    list(
      predict = predict_gamma,
      gamma_estimator = "correct_pooled_direct_exponential_tilt",
      selected_lambda = ridge,
      selected_bandwidth_multiplier = NA_real_,
      base_bandwidth = NA_real_,
      selected_bandwidth = NA_real_,
      cv_table = data.frame(),
      normalization_constant = normalization_constant,
      training_negative_rate = 0,
      rkhs_backend = "correct_pooled_direct_exponential_tilt",
      nystrom_rank_requested = NA_integer_,
      nystrom_landmark_count = NA_integer_,
      nystrom_effective_rank = NA_integer_,
      nystrom_min_retained_eigenvalue = NA_real_,
      nystrom_max_eigenvalue = NA_real_,
      ratio_by_atom = ratio_by_atom,
      baseline_atom_count = count0,
      comparison_atom_count = count1,
      optimization_convergence = optimization$convergence,
      maximum_absolute_score = max(abs(gradient(optimization$par)))
    )
  }
}

make_mirrored_ate_parametric_variance_grid <- function(
    train,
    gate = build_mirrored_ate_nonunique_threshold_gate(),
    pseudocount = 0.5
) {
  required <- c("Y", "T_obs", "Z", "block", "block_type", "local_X")
  if (!all(required %in% names(train))) {
    stop("The parametric variance grid received incomplete training data.")
  }
  pseudocount <- as.numeric(pseudocount)
  if (length(pseudocount) != 1L || !is.finite(pseudocount) ||
      pseudocount <= 0) {
    stop("pseudocount must be one positive number.")
  }
  specification <- gate$upper$specification
  cells <- gate$groups * gate$local_strata
  observed_cell <- (train$block - 1L) * gate$local_strata + train$local_X
  cell_count <- tabulate(observed_cell, nbins = cells)
  cell_probability <- (cell_count + pseudocount) /
    (sum(cell_count) + pseudocount * cells)
  probability_z0 <- (sum(train$Z == 0L) + pseudocount) /
    (nrow(train) + 2 * pseudocount)
  treatment_probability <- array(NA_real_, dim = c(4L, 2L))
  for (type in 1:4) {
    for (z_value in 0:1) {
      rows <- train$block_type == type & train$Z == z_value
      treatment_probability[type, z_value + 1L] <-
        (sum(train$T_obs[rows] == 1L) + pseudocount) /
        (sum(rows) + 2 * pseudocount)
    }
  }
  outcome_probability <- vector("list", gate$groups)
  for (group in seq_len(gate$groups)) {
    type <- gate$group_type[group]
    active_treatment <- if (type %in% c(1L, 3L)) 1L else 0L
    reflected <- type %in% c(1L, 4L)
    support <- if (reflected) -gate$outcome_atoms else gate$outcome_atoms
    probability <- matrix(NA_real_, nrow = length(support), ncol = 2L)
    for (z_value in 0:1) {
      rows <- train$block == group & train$T_obs == active_treatment &
        train$Z == z_value
      index <- match(train$Y[rows], support)
      count <- tabulate(index, nbins = length(support))
      probability[, z_value + 1L] <- (count + pseudocount) /
        (sum(count) + pseudocount * length(support))
    }
    outcome_probability[[group]] <- list(
      active_treatment = active_treatment,
      reflected = reflected,
      support = support,
      probability = probability
    )
  }
  rows <- list()
  probabilities <- numeric()
  counter <- 1L
  for (cell in seq_len(cells)) {
    group <- (cell - 1L) %/% gate$local_strata + 1L
    local_x <- (cell - 1L) %% gate$local_strata + 1L
    type <- gate$group_type[group]
    information <- outcome_probability[[group]]
    for (z_value in 0:1) {
      probability_z <- if (z_value == 0L) probability_z0 else 1 - probability_z0
      pi <- treatment_probability[type, z_value + 1L]
      for (treatment in 0:1) {
        probability_t <- if (treatment == 1L) pi else 1 - pi
        support <- if (treatment == information$active_treatment) {
          information$support
        } else {
          0
        }
        probability_y <- if (treatment == information$active_treatment) {
          information$probability[, z_value + 1L]
        } else {
          1
        }
        rows[[counter]] <- data.frame(
          Y = support,
          T_obs = treatment,
          Z = z_value,
          X = cell,
          block = group,
          block_type = type,
          local_X = local_x,
          atom = if (treatment == information$active_treatment) {
            seq_along(support)
          } else {
            1L
          },
          stringsAsFactors = FALSE
        )
        probabilities <- c(
          probabilities,
          cell_probability[cell] * probability_z * probability_t *
            probability_y
        )
        counter <- counter + 1L
      }
    }
  }
  grid <- do.call(rbind, rows)
  rownames(grid) <- NULL
  active_treatment <- ifelse(grid$block_type %in% c(1L, 3L), 1L, 0L)
  active <- grid$T_obs == active_treatment
  reflected <- grid$block_type %in% c(1L, 4L)
  oracle_atom <- match(
    ifelse(reflected, -grid$Y, grid$Y),
    gate$outcome_atoms
  )
  gamma <- rep(1, nrow(grid))
  gamma[active] <- specification$gamma_atoms[oracle_atom[active]]
  grid$e_true <- specification$probability_z0
  grid$pi0_true <- ifelse(
    active_treatment == 1L, specification$pi0, 1 - specification$pi0
  )
  grid$pi1_true <- ifelse(
    active_treatment == 1L, specification$pi1, 1 - specification$pi1
  )
  grid$gamma0_oracle <- ifelse(grid$T_obs == 0L & active, gamma, 1)
  grid$gamma1_oracle <- ifelse(grid$T_obs == 1L & active, gamma, 1)
  grid$prob_Z <- 1 - specification$probability_z0
  grid$prob_T <- ifelse(
    grid$Z == 0L, grid$pi0_true, grid$pi1_true
  )
  for (index in seq_len(gate$groups - 1L)) {
    grid[[paste0("block", index)]] <- as.numeric(grid$block == index)
  }
  for (index in seq_len(gate$local_strata - 1L)) {
    grid[[paste0("local_x", index)]] <- as.numeric(grid$local_X == index)
  }
  cell <- (grid$block - 1L) * gate$local_strata + grid$local_X
  for (index in seq_len(cells - 1L)) {
    grid[[paste0("cell", index)]] <- as.numeric(cell == index)
  }
  probabilities <- probabilities / sum(probabilities)
  list(
    data = grid,
    probability = probabilities,
    states = nrow(grid),
    training_observations = nrow(train)
  )
}

make_mirrored_ate_finite_library_provider_factory <- function(
    gate = build_mirrored_ate_nonunique_threshold_gate()
) {
  function(mode, fold, treatment, gamma_fit, nuisance_fit, x_cols) {
    force(mode)
    force(fold)
    force(gamma_fit)
    force(nuisance_fit)
    force(x_cols)
    list(
      thresholds = function(
          train,
          validation,
          nuisance_train,
          nuisance_validation,
          envelope_train,
          envelope_validation,
          treatment,
          Lambda0,
          Lambda1,
          tau_train,
          tau_validation,
          compatibility_action,
          gamma_projection_margin,
          envelope_smoothing_epsilon,
          tau_epsilon
      ) {
        force(nuisance_validation)
        force(envelope_validation)
        force(tau_validation)
        force(compatibility_action)
        force(gamma_projection_margin)
        force(envelope_smoothing_epsilon)
        force(tau_epsilon)
        if (
          abs(Lambda0 - 4) > 1e-12 || abs(Lambda1 - 4) > 1e-12 ||
            !all(c("block", "block_type", "local_X") %in% names(train))
        ) {
          stop("The mirrored finite library received incompatible inputs.")
        }
        active_types <- if (treatment == 1L) c(1L, 3L) else c(2L, 4L)
        selected <- vector("list", gate$groups)
        for (group in seq_len(gate$groups)) {
          type <- gate$group_type[group]
          if (!type %in% active_types) next
          rows <- which(train$block == group)
          if (length(rows) < 2L) {
            stop("Too few observations in a mirrored threshold block.")
          }
          reflected <- type %in% c(1L, 4L)
          local_train <- train[rows, , drop = FALSE]
          local_train$X <- local_train$local_X
          if (reflected) local_train$Y <- -local_train$Y
          choose <- function(component) {
            baseline_count <- sum(
              local_train$Z == 0L & local_train$T_obs == treatment
            )
            # Early online prefixes can contain fewer than two observations
            # in a rare block/arm baseline cell.  Candidate one is a fixed,
            # predictable fallback and is itself an exact population optimizer.
            if (baseline_count < 2L) return(1L)
            select_two_gap_threshold_candidate(
              train = local_train,
              nuisance_train = lapply(nuisance_train, function(value) value[rows]),
              envelope_train = lapply(envelope_train, function(value) value[rows]),
              tau_train = tau_train[rows],
              treatment = treatment,
              specification = component$specification
            )
          }
          selected[[group]] <- list(
            reflected = reflected,
            lower_component = if (reflected) gate$upper else gate$lower,
            upper_component = if (reflected) gate$lower else gate$upper,
            lower_candidate = choose(if (reflected) gate$upper else gate$lower),
            upper_candidate = choose(if (reflected) gate$lower else gate$upper)
          )
        }
        predict_endpoint <- function(new_data, endpoint) {
          prediction <- rep(0, nrow(new_data))
          for (group in seq_len(gate$groups)) {
            information <- selected[[group]]
            if (is.null(information)) next
            rows <- which(new_data$block == group)
            if (!length(rows)) next
            local <- new_data[rows, , drop = FALSE]
            local$X <- local$local_X
            component <- information[[paste0(endpoint, "_component")]]
            candidate <- information[[paste0(endpoint, "_candidate")]]
            value <- selected_two_gap_candidate_thresholds(
              local, component$specification
            )[, candidate]
            if (information$reflected) value <- -value
            prediction[rows] <- value
          }
          prediction
        }
        lower_predict <- function(new_data) predict_endpoint(new_data, "lower")
        upper_predict <- function(new_data) predict_endpoint(new_data, "upper")
        list(
          q_lower_train = lower_predict(train),
          q_upper_train = upper_predict(train),
          q_lower_validation = lower_predict(validation),
          q_upper_validation = upper_predict(validation),
          lower_predict = lower_predict,
          upper_predict = upper_predict
        )
      }
    )
  }
}

make_mirrored_ate_saturated_gamma_provider_factory <- function(
    gate = build_mirrored_ate_nonunique_threshold_gate(),
    pseudocount = 0.5
) {
  pseudocount <- as.numeric(pseudocount)
  if (length(pseudocount) != 1L || !is.finite(pseudocount) || pseudocount <= 0) {
    stop("pseudocount must be one positive number.")
  }
  function(fold, treatment, train, nuisance_train, x_cols, seed, ...) {
    force(fold)
    force(nuisance_train)
    force(x_cols)
    force(seed)
    lookup <- vector("list", gate$groups)
    for (group in seq_len(gate$groups)) {
      type <- gate$group_type[group]
      active_treatment <- if (type %in% c(1L, 3L)) 1L else 0L
      reflected <- type %in% c(1L, 4L)
      # The estimator also asks for predictions on an outcome grid crossed
      # with every covariate block.  Include the full known support; values
      # impossible in a given block receive the same pseudocount under both Z
      # levels and therefore gamma one.
      support <- sort(unique(c(
        0, gate$outcome_atoms, -gate$outcome_atoms
      )))
      probability <- matrix(
        NA_real_, nrow = length(support), ncol = 2L,
        dimnames = list(as.character(seq_along(support)), c("z0", "z1"))
      )
      for (z_value in 0:1) {
        rows <- train$block == group & train$T_obs == treatment &
          train$Z == z_value
        outcome_index <- match(train$Y[rows], support)
        count <- tabulate(outcome_index, nbins = length(support))
        probability[, z_value + 1L] <- (count + pseudocount) /
          (sum(count) + pseudocount * length(support))
      }
      lookup[[group]] <- list(
        support = support,
        gamma = probability[, "z1"] / probability[, "z0"]
      )
    }
    predict_gamma <- function(newdata) {
      value <- rep(NA_real_, nrow(newdata))
      for (group in seq_len(gate$groups)) {
        rows <- which(newdata$block == group)
        if (!length(rows)) next
        index <- match(newdata$Y[rows], lookup[[group]]$support)
        if (anyNA(index)) {
          stop("A validation outcome is outside the known finite support.")
        }
        value[rows] <- lookup[[group]]$gamma[index]
      }
      list(raw = value, positive = value, normalized = value)
    }
    list(
      predict = predict_gamma,
      gamma_estimator = "cross_fitted_saturated_multinomial",
      selected_lambda = NA_real_,
      selected_bandwidth_multiplier = NA_real_,
      base_bandwidth = NA_real_,
      selected_bandwidth = NA_real_,
      cv_table = data.frame(),
      normalization_constant = 1,
      training_negative_rate = 0,
      rkhs_backend = "saturated_multinomial",
      nystrom_rank_requested = NA_integer_,
      nystrom_landmark_count = NA_integer_,
      nystrom_effective_rank = NA_integer_,
      nystrom_min_retained_eigenvalue = NA_real_,
      nystrom_max_eigenvalue = NA_real_
    )
  }
}

mirrored_ate_gate_endpoint_map <- function(gate, endpoint) {
  endpoint <- match.arg(endpoint, c("lower", "upper"))
  type <- gate$group_type
  if (endpoint == "lower") {
    use_upper <- type %in% c(1L, 2L)
    sign <- ifelse(use_upper, -1, 1)
  } else {
    use_upper <- type %in% c(3L, 4L)
    sign <- ifelse(use_upper, 1, -1)
  }
  list(
    component = lapply(use_upper, function(value) {
      if (value) gate$upper else gate$lower
    }),
    sign = sign,
    target = if (endpoint == "lower") gate$lower_target else gate$upper_target
  )
}

mirrored_ate_gate_summarize <- function(
    sample_size,
    estimator,
    lower_estimate,
    upper_estimate,
    lower_standard_error,
    upper_standard_error,
    gate,
    initial_size = NA_integer_,
    chunks = NA_integer_
) {
  proxy <- structure(
    list(
      block_pairs = gate$block_quartets,
      local_strata = gate$local_strata,
      lower_target = gate$lower_target,
      upper_target = gate$upper_target
    ),
    class = "multiblock_nonunique_threshold_gate"
  )
  output <- multiblock_gate_summarize(
    sample_size = sample_size,
    estimator = estimator,
    lower_estimate = lower_estimate,
    upper_estimate = upper_estimate,
    lower_standard_error = lower_standard_error,
    upper_standard_error = upper_standard_error,
    gate = proxy,
    initial_size = initial_size,
    chunks = chunks
  )
  names(output)[names(output) == "block_pairs"] <- "block_quartets"
  output
}

simulate_mirrored_ate_gate_crossfit <- function(
    gate = build_mirrored_ate_nonunique_threshold_gate(),
    sample_size = 3000L,
    replications = 1000L,
    folds = 5L,
    seed = 510901L,
    batch_size = 20L
) {
  sample_size <- as.integer(sample_size)
  replications <- as.integer(replications)
  folds <- as.integer(folds)
  batch_size <- as.integer(batch_size)
  if (
    sample_size < folds || sample_size %% folds != 0L ||
      replications < 1L || folds < 2L || batch_size < 1L
  ) {
    stop("Invalid mirrored ATE cross-fit simulation inputs.")
  }
  category_count <- nrow(gate$upper$categories)
  candidate_count <- ncol(gate$upper$score)
  lower_estimate <- upper_estimate <- numeric(replications)
  lower_standard_error <- upper_standard_error <- numeric(replications)
  set.seed(as.integer(seed))
  completed <- 0L
  while (completed < replications) {
    current <- min(batch_size, replications - completed)
    counts <- array(
      stats::rmultinom(
        current * folds,
        sample_size / folds,
        gate$category_probability
      ),
      dim = c(category_count, gate$groups, folds, current)
    )
    calculate_endpoint <- function(endpoint) {
      mapping <- mirrored_ate_gate_endpoint_map(gate, endpoint)
      signal_sum <- signal_square_sum <- numeric(current)
      replication_index <- seq_len(current)
      for (group in seq_len(gate$groups)) {
        component <- mapping$component[[group]]
        sign <- mapping$sign[group]
        local_count <- matrix(
          counts[, group, , ],
          nrow = category_count,
          ncol = folds * current
        )
        fold_loss <- array(
          t(local_count) %*% component$selection_loss,
          dim = c(folds, current, candidate_count)
        )
        fold_score <- array(
          t(local_count) %*% component$score,
          dim = c(folds, current, candidate_count)
        )
        fold_score_squared <- array(
          t(local_count) %*% component$score^2,
          dim = c(folds, current, candidate_count)
        )
        total_loss <- apply(fold_loss, c(2L, 3L), sum)
        for (fold in seq_len(folds)) {
          training_loss <- total_loss - matrix(
            fold_loss[fold, , ], nrow = current
          )
          selected <- max.col(-training_loss, ties.method = "first")
          index <- cbind(replication_index, selected)
          signal_sum <- signal_sum + sign * matrix(
            fold_score[fold, , ], nrow = current
          )[index]
          signal_square_sum <- signal_square_sum + matrix(
            fold_score_squared[fold, , ], nrow = current
          )[index]
        }
      }
      estimate <- signal_sum / sample_size
      standard_error <- sqrt(pmax(
        (signal_square_sum - sample_size * estimate^2) /
          (sample_size - 1),
        0
      ) / sample_size)
      list(estimate = estimate, standard_error = standard_error)
    }
    lower <- calculate_endpoint("lower")
    upper <- calculate_endpoint("upper")
    rows <- completed + seq_len(current)
    lower_estimate[rows] <- lower$estimate
    upper_estimate[rows] <- upper$estimate
    lower_standard_error[rows] <- lower$standard_error
    upper_standard_error[rows] <- upper$standard_error
    completed <- completed + current
  }
  mirrored_ate_gate_summarize(
    sample_size = sample_size,
    estimator = paste0("oracle_", folds, "_fold_crossfit"),
    lower_estimate = lower_estimate,
    upper_estimate = upper_estimate,
    lower_standard_error = lower_standard_error,
    upper_standard_error = upper_standard_error,
    gate = gate
  )
}

simulate_mirrored_ate_gate_online <- function(
    gate = build_mirrored_ate_nonunique_threshold_gate(),
    sample_size = 3000L,
    replications = 1000L,
    initial_exponent = 0.7,
    chunks = 10L,
    seed = 610901L
) {
  sample_size <- as.integer(sample_size)
  replications <- as.integer(replications)
  chunks <- as.integer(chunks)
  initial_size <- max(200L, ceiling(sample_size^initial_exponent))
  evaluation_count <- sample_size - initial_size
  if (evaluation_count < chunks || replications < 1L || chunks < 1L) {
    stop("Invalid mirrored ATE online simulation inputs.")
  }
  chunk_sizes <- tabulate(
    rep(seq_len(chunks), length.out = evaluation_count),
    nbins = chunks
  )
  category_count <- nrow(gate$upper$categories)
  lower_estimate <- upper_estimate <- numeric(replications)
  lower_standard_error <- upper_standard_error <- numeric(replications)
  set.seed(as.integer(seed))
  for (replication in seq_len(replications)) {
    training <- array(
      stats::rmultinom(1L, initial_size, gate$category_probability),
      dim = c(category_count, gate$groups)
    )
    evaluation <- lapply(chunk_sizes, function(size) {
      array(
        stats::rmultinom(1L, size, gate$category_probability),
        dim = c(category_count, gate$groups)
      )
    })
    calculate_endpoint <- function(endpoint) {
      mapping <- mirrored_ate_gate_endpoint_map(gate, endpoint)
      local_training <- training
      numerator <- denominator <- 0
      for (chunk in seq_len(chunks)) {
        chunk_signal_sum <- selected_second_moment <- 0
        for (group in seq_len(gate$groups)) {
          component <- mapping$component[[group]]
          sign <- mapping$sign[group]
          selected <- which.min(as.vector(
            local_training[, group] %*% component$selection_loss
          ))
          count <- evaluation[[chunk]][, group]
          chunk_signal_sum <- chunk_signal_sum + sign * sum(
            count * component$score[, selected]
          )
          selected_second_moment <- selected_second_moment + sum(
            component$categories$probability *
              component$score[, selected]^2
          )
        }
        selected_second_moment <- selected_second_moment / gate$groups
        selected_standard_deviation <- sqrt(pmax(
          selected_second_moment - mapping$target^2,
          .Machine$double.eps
        ))
        numerator <- numerator + chunk_signal_sum /
          selected_standard_deviation
        denominator <- denominator + chunk_sizes[chunk] /
          selected_standard_deviation
        local_training <- local_training + evaluation[[chunk]]
      }
      c(
        estimate = numerator / denominator,
        standard_error = (evaluation_count / denominator) /
          sqrt(evaluation_count)
      )
    }
    lower <- calculate_endpoint("lower")
    upper <- calculate_endpoint("upper")
    lower_estimate[replication] <- lower[["estimate"]]
    upper_estimate[replication] <- upper[["estimate"]]
    lower_standard_error[replication] <- lower[["standard_error"]]
    upper_standard_error[replication] <- upper[["standard_error"]]
  }
  mirrored_ate_gate_summarize(
    sample_size = sample_size,
    estimator = "oracle_online_predictable",
    lower_estimate = lower_estimate,
    upper_estimate = upper_estimate,
    lower_standard_error = lower_standard_error,
    upper_standard_error = upper_standard_error,
    gate = gate,
    initial_size = initial_size,
    chunks = chunks
  )
}
