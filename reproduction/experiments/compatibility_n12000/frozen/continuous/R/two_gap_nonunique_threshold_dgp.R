# Selected two-gap IV law for outward confidence-set stress testing.
#
# At Lambda = 4 the observed treated-arm law has six ordered outcome atoms.
# Its probabilities solve the dual lower- and upper-extremizer normalization
# equations and make both Z-specific outcome CDFs equal 0.20 after atom 2 and
# 0.80 after atom 4.  Consequently, the dual IV-MSM and both Tan comparators
# have exact lower and upper threshold ties.  Every envelope branch remains
# separated from its max/min kink, so this isolates threshold nonuniqueness
# rather than mixing it with an envelope tie.
#
# The `joint_adverse` profile is the compatible candidate selected as row 36
# of the broad joint dual/Tan search.  At 16 strata, its exact finite-sample
# oracle experiment gives approximately 89% five-fold Bonferroni outward-set
# coverage and 95% predictable-online coverage at n = 3000.  The profile is
# kept separate from the original law so the earlier audit remains exactly
# reproducible.
#
# The `joint_adverse_dense` profile is row 113 of the same search.  It retains
# the exact two-threshold tie and adverse joint five-fold/Tan geometry, but its
# observed outcome atoms are less sparse.  This makes it useful for checking
# whether nuisance-estimation bias is masking the nonregular cross-fit limit.
# The `joint_adverse_spread` profile keeps row 36's probabilities and density
# ratios but separates adjacent outcome atoms.  The population tie is
# unchanged, while the RKHS learner can distinguish the alternating ratios
# without resolving nearly coincident outcome values.

if (!exists("build_multibranch_nonunique_threshold_oracle", mode = "function")) {
  if (!file.exists("R/multibranch_nonunique_threshold_pilot.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/multibranch_nonunique_threshold_pilot.R")
}

make_selected_two_gap_nonunique_threshold_specs <- function(
    strata = 64L,
    atom_layout = c(
      "original", "learnable", "gamma_one", "joint_adverse",
      "joint_adverse_dense", "joint_adverse_spread"
    )
) {
  strata <- as.integer(strata)
  atom_layout <- match.arg(atom_layout)
  if (
    length(strata) != 1L || is.na(strata) || strata < 2L ||
      bitwAnd(strata, strata - 1L) != 0L
  ) {
    stop("strata must be a power of two of at least two.")
  }
  common <- list(
    sensitivity_lambda = 4,
    strata = strata,
    candidates = strata,
    candidate_design = "hadamard",
    treatment = 1L,
    probability_z0 = if (atom_layout == "joint_adverse_dense") 0.2 else 0.5,
    pi0 = if (atom_layout %in% c(
      "joint_adverse", "joint_adverse_spread"
    )) {
      0.15
    } else if (atom_layout == "joint_adverse_dense") {
      0.30
    } else {
      0.20
    },
    pi1 = if (atom_layout %in% c(
      "joint_adverse", "joint_adverse_spread"
    )) {
      0.45
    } else if (atom_layout == "joint_adverse_dense") {
      0.55
    } else {
      0.35
    },
    atoms = if (atom_layout == "joint_adverse_spread") {
      c(0, 1, 4, 5, 8, 9)
    } else if (atom_layout %in% c(
      "original", "joint_adverse", "joint_adverse_dense"
    )) {
      c(
        0,
        0.0367812743571667,
        1.56242006332992,
        1.63548558874866,
        4.42124442948423,
        4.56552575820334
      )
    } else {
      c(
        0,
        0.0658675774870749,
        0.466091290211156,
        0.532568350642681,
        0.928640765888784,
        1
      )
    },
    tail_probability = 0.05,
    threshold_gap_fractions = c(
      4.22684506259259e-06,
      0.999809663362851
    ),
    gamma_atoms = if (atom_layout == "gamma_one") {
      rep(1, 6L)
    } else if (atom_layout %in% c(
      "joint_adverse", "joint_adverse_spread"
    )) {
      c(
        1.50349184907417865,
        0.87210735344821044,
        0.19457616042003101,
        1.03320758352549857,
        0.33069659406632762,
        1.35667701324935974
      )
    } else if (atom_layout == "joint_adverse_dense") {
      c(
        0.291712545112464,
        1.55443890711659,
        1.38974639056042,
        0.333854854514729,
        1.98973589215759,
        0.657912202189576
      )
    } else {
      c(
        0.581441128873817,
        1.58587398962943,
        2.25833660522745,
        0.320828155669296,
        1.50118767162008,
        0.490475378749418
      )
    },
    atom_probability_override = if (atom_layout %in% c(
      "joint_adverse", "joint_adverse_spread"
    )) {
      c(
        0.040511810929089717,
        0.159488189070910308,
        0.023758411104509082,
        0.576241588895490930,
        0.069529009829130003,
        0.130470990170870021
      )
    } else if (atom_layout == "joint_adverse_dense") {
      c(
        0.0878161609355516,
        0.112183839064448,
        0.378530439582828,
        0.221469560417172,
        0.0513713339666814,
        0.148628666033319
      )
    } else {
      c(
        0.116657670715529,
        0.0833423292844705,
        0.210323266817935,
        0.389676733182065,
        0.100824858833648,
        0.0991751411663523
      )
    }
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

build_selected_two_gap_nonunique_threshold_oracle <- function(
    strata = 64L,
    atom_layout = c(
      "original", "learnable", "gamma_one", "joint_adverse",
      "joint_adverse_dense", "joint_adverse_spread"
    )
) {
  atom_layout <- match.arg(atom_layout)
  specifications <- make_selected_two_gap_nonunique_threshold_specs(
    strata, atom_layout = atom_layout
  )
  lower <- build_multibranch_nonunique_threshold_oracle(
    specifications$lower
  )
  upper <- build_multibranch_nonunique_threshold_oracle(
    specifications$upper
  )
  if (
    nrow(lower$categories) != nrow(upper$categories) ||
      max(abs(
        lower$categories$probability - upper$categories$probability
      )) > 1e-12
  ) {
    stop("The selected lower and upper endpoints do not share one law.")
  }
  list(
    name = paste0("selected_two_gap_nonunique_threshold_", atom_layout),
    atom_layout = atom_layout,
    specifications = specifications,
    categories = lower$categories,
    lower = lower,
    upper = upper
  )
}

selected_two_gap_x_columns <- function(oracle) {
  paste0("X", seq_len(oracle$lower$specification$strata - 1L))
}

# Finite-library threshold learner for the exact two-gap stress experiment.
#
# The library is the Hadamard collection used to define the population law.
# Every member is an exact population optimizer of the unrestricted sharp
# endpoint, so selecting any member has exactly zero population regret.  The
# selected member is nevertheless data adaptive: within each training fold we
# minimize the same weighted check loss used by fit_bound_arm_fold().  This
# deliberately isolates failure of stabilization from threshold approximation
# or optimization error while leaving all nuisance and EIF machinery in the
# official one-step estimator unchanged.
selected_two_gap_candidate_thresholds <- function(data, specification) {
  if (!is.data.frame(data) || !"X" %in% names(data)) {
    stop("Two-gap threshold predictions require the integer X stratum.")
  }
  x <- as.integer(data$X)
  if (
    anyNA(x) || any(!x %in% seq_len(specification$strata)) ||
      nrow(specification$candidate_right_endpoint) !=
        specification$candidates
  ) {
    stop("The two-gap threshold strata or candidate library are invalid.")
  }
  gap_left <- specification$atoms[specification$gap_index]
  gap_right <- specification$atoms[specification$gap_index + 1L]
  right <- specification$candidate_right_endpoint[, x, drop = FALSE]
  fraction <- ifelse(
    right,
    specification$threshold_gap_fractions[2L],
    specification$threshold_gap_fractions[1L]
  )
  t(gap_left + (gap_right - gap_left) * fraction)
}

select_two_gap_threshold_candidate <- function(
    train,
    nuisance_train,
    envelope_train,
    tau_train,
    treatment,
    specification
) {
  baseline <- train$Z == 0L & train$T_obs == treatment
  if (sum(baseline) < 2L) {
    stop("Too few baseline observations for finite-library selection.")
  }
  candidate <- selected_two_gap_candidate_thresholds(
    train, specification
  )
  target_odds <- compute_target_odds(
    nuisance = nuisance_train,
    t_value = treatment
  )
  weight <- target_odds[baseline] * envelope_train$D[baseline] /
    nuisance_train$e[baseline]
  level <- if (specification$endpoint_direction == "lower") {
    1 - tau_train[baseline]
  } else {
    tau_train[baseline]
  }
  risk <- vapply(seq_len(ncol(candidate)), function(index) {
    sum(
      weight * check_loss(
        train$Y[baseline],
        candidate[baseline, index],
        level
      )
    )
  }, numeric(1))
  which.min(risk)
}

make_selected_two_gap_finite_library_provider_factory <- function(
    oracle = build_selected_two_gap_nonunique_threshold_oracle(
      strata = 16L,
      atom_layout = "joint_adverse"
    )
) {
  if (
    !is.list(oracle) ||
      !all(c("lower", "upper") %in% names(oracle))
  ) {
    stop("A complete selected two-gap oracle is required.")
  }
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
          abs(Lambda0 - oracle$lower$specification$sensitivity_lambda) >
            1e-12 ||
            abs(Lambda1 - oracle$lower$specification$sensitivity_lambda) >
              1e-12
        ) {
          stop("The finite two-gap library was called at the wrong Lambda.")
        }

        # The control outcome is identically zero in the full observed-data
        # law.  Its contribution to both ATE endpoints is therefore zero for
        # the common threshold returned here.
        if (treatment == 0L) {
          lower_predict <- function(new_data) rep(0, nrow(new_data))
          upper_predict <- function(new_data) rep(0, nrow(new_data))
          return(list(
            q_lower_train = lower_predict(train),
            q_upper_train = upper_predict(train),
            q_lower_validation = lower_predict(validation),
            q_upper_validation = upper_predict(validation),
            lower_predict = lower_predict,
            upper_predict = upper_predict
          ))
        }

        lower_specification <- oracle$lower$specification
        upper_specification <- oracle$upper$specification
        selected_lower <- select_two_gap_threshold_candidate(
          train = train,
          nuisance_train = nuisance_train,
          envelope_train = envelope_train,
          tau_train = tau_train,
          treatment = treatment,
          specification = lower_specification
        )
        selected_upper <- select_two_gap_threshold_candidate(
          train = train,
          nuisance_train = nuisance_train,
          envelope_train = envelope_train,
          tau_train = tau_train,
          treatment = treatment,
          specification = upper_specification
        )
        lower_predict <- function(new_data) {
          selected_two_gap_candidate_thresholds(
            new_data, lower_specification
          )[, selected_lower]
        }
        upper_predict <- function(new_data) {
          selected_two_gap_candidate_thresholds(
            new_data, upper_specification
          )[, selected_upper]
        }
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

selected_two_gap_add_x_columns <- function(data, oracle) {
  strata <- oracle$lower$specification$strata
  for (stratum in seq_len(strata - 1L)) {
    data[[paste0("X", stratum)]] <- as.numeric(data$X == stratum)
  }
  data
}

# Full-estimation version of the selected law.  It preserves the treated-arm
# two-gap law exactly and fixes the control outcome at zero.  Consequently the
# treated lower/upper endpoints are also the lower/upper ATE endpoints, without
# regular control-arm noise diluting the nonregular threshold experiment.
generate_selected_two_gap_full_data <- function(
    n,
    oracle = build_selected_two_gap_nonunique_threshold_oracle(strata = 16L),
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
  treatment_probability <- ifelse(
    z == 0L, specification$pi0, specification$pi1
  )
  treatment <- stats::rbinom(n, 1L, treatment_probability)
  probability0 <- specification$atom_probability_override
  probability1 <- probability0 * specification$gamma_atoms
  atom <- rep(1L, n)
  treated <- treatment == 1L
  treated_groups <- split(which(treated), z[treated])
  for (group_name in names(treated_groups)) {
    rows <- treated_groups[[group_name]]
    atom[rows] <- sample.int(
      length(specification$atoms), length(rows), replace = TRUE,
      prob = if (identical(group_name, "1")) probability1 else probability0
    )
  }
  output <- data.frame(
    Y = ifelse(treated, specification$atoms[atom], 0),
    T_obs = treatment,
    Z = z,
    X = x,
    atom = atom,
    e_true = specification$probability_z0,
    pi0_true = specification$pi0,
    pi1_true = specification$pi1,
    gamma0_oracle = 1,
    gamma1_oracle = specification$gamma_atoms[atom],
    prob_Z = 1 - specification$probability_z0,
    prob_T = treatment_probability,
    stringsAsFactors = FALSE
  )
  output <- selected_two_gap_add_x_columns(output, oracle)
  attr(output, "selected_two_gap_oracle") <- oracle
  output
}

compute_selected_two_gap_full_dual_population <- function(
    oracle = build_selected_two_gap_nonunique_threshold_oracle(strata = 16L)
) {
  specification <- oracle$lower$specification
  endpoints <- data.frame(
    design = "selected_two_gap_full_ate_threshold_tie",
    Lambda = specification$sensitivity_lambda,
    mu1_lower = oracle$lower$target,
    mu1_upper = oracle$upper$target,
    mu0_lower = 0,
    mu0_upper = 0,
    stringsAsFactors = FALSE
  )
  endpoints$ate_lower <- endpoints$mu1_lower
  endpoints$ate_upper <- endpoints$mu1_upper
  list(
    endpoints = endpoints,
    diagnostics = audit_selected_two_gap_nonunique_threshold(oracle),
    compatible = TRUE
  )
}

compute_selected_two_gap_full_tan_population <- function(
    oracle = build_selected_two_gap_nonunique_threshold_oracle(strata = 16L)
) {
  arm1 <- compute_selected_two_gap_tan_population(oracle)
  estimator <- c(
    tan_cal_x_intercept = "tan_cal_x_saturated",
    tan_cal_xz = "tan_cal_xz_saturated"
  )
  output <- data.frame(
    estimator = unname(estimator[arm1$estimator]),
    Lambda = arm1$Lambda,
    mu1_lower = arm1$lower_target,
    mu1_upper = arm1$upper_target,
    mu0_lower = 0,
    mu0_upper = 0,
    stringsAsFactors = FALSE
  )
  output$ate_lower <- output$mu1_lower
  output$ate_upper <- output$mu1_upper
  output
}

audit_selected_two_gap_nonunique_threshold <- function(
    oracle = build_selected_two_gap_nonunique_threshold_oracle()
) {
  lower <- oracle$lower
  upper <- oracle$upper
  specification <- lower$specification
  atom_probability0 <- specification$atom_probability_override
  atom_probability1 <- atom_probability0 * specification$gamma_atoms
  atom_envelope <- compute_projected_envelopes(
    nuisance = list(
      e = rep(specification$probability_z0, length(atom_probability0)),
      pi0 = rep(specification$pi0, length(atom_probability0)),
      pi1 = rep(specification$pi1, length(atom_probability0))
    ),
    gamma = specification$gamma_atoms,
    t_value = specification$treatment,
    Lambda0 = specification$sensitivity_lambda,
    Lambda1 = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = 0,
    gamma_compatibility_action = "stop"
  )
  lower_tie_index <- specification$gap_index
  upper_tie_index <- upper$specification$gap_index
  data.frame(
    Lambda = lower$specification$sensitivity_lambda,
    strata = lower$specification$strata,
    candidates = lower$specification$candidates,
    treatment = lower$specification$treatment,
    probability_z0 = lower$specification$probability_z0,
    pi0 = lower$specification$pi0,
    pi1 = lower$specification$pi1,
    minimum_observed_tz_probability = min(c(
      lower$specification$probability_z0 * lower$specification$pi0,
      lower$specification$probability_z0 * (1 - lower$specification$pi0),
      (1 - lower$specification$probability_z0) * lower$specification$pi1,
      (1 - lower$specification$probability_z0) *
        (1 - lower$specification$pi1)
    )),
    lower_target = lower$target,
    upper_target = upper$target,
    lower_threshold_level = lower$threshold_level,
    upper_threshold_level = upper$threshold_level,
    lower_risk_range = diff(range(lower$population_risk)),
    upper_risk_range = diff(range(upper$population_risk)),
    lower_target_range = diff(range(lower$direct_target)),
    upper_target_range = diff(range(upper$direct_target)),
    lower_maximum_normalization_residual = max(
      abs(lower$direct_normalization_residual)
    ),
    upper_maximum_normalization_residual = max(
      abs(upper$direct_normalization_residual)
    ),
    tan_z0_lower_cdf_residual =
      sum(atom_probability0[seq_len(lower_tie_index)]) - 0.2,
    tan_z0_upper_cdf_residual =
      sum(atom_probability0[seq_len(upper_tie_index)]) - 0.8,
    tan_z1_lower_cdf_residual =
      sum(atom_probability1[seq_len(lower_tie_index)]) - 0.2,
    tan_z1_upper_cdf_residual =
      sum(atom_probability1[seq_len(upper_tie_index)]) - 0.8,
    minimum_envelope_width = min(atom_envelope$U - atom_envelope$L),
    minimum_envelope_branch_gap = min(c(
      abs(atom_envelope$candidate_lower -
        1 / specification$sensitivity_lambda),
      abs(atom_envelope$candidate_upper -
        specification$sensitivity_lambda)
    )),
    stringsAsFactors = FALSE
  )
}

two_gap_bind_rows <- function(rows) {
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

simulate_selected_two_gap_crossfit <- function(
    oracle = build_selected_two_gap_nonunique_threshold_oracle(),
    sample_sizes = c(1500L, 3000L, 6000L),
    replications = 1000L,
    folds = 5L,
    seed = 1400000L,
    batch_size = 25L
) {
  sample_sizes <- as.integer(sample_sizes)
  replications <- as.integer(replications)
  folds <- as.integer(folds)
  batch_size <- as.integer(batch_size)
  if (
    !length(sample_sizes) || any(sample_sizes < folds) ||
      replications < 1L || folds < 2L || batch_size < 1L ||
      any(sample_sizes %% folds != 0L)
  ) {
    stop("Invalid selected two-gap cross-fit simulation inputs.")
  }
  lower <- oracle$lower
  upper <- oracle$upper
  probability <- oracle$categories$probability
  candidate_count <- ncol(lower$score)
  simulate_n <- function(n, local_seed) {
    set.seed(as.integer(local_seed))
    z_lower <- z_upper <- numeric(replications)
    estimate_lower <- estimate_upper <- numeric(replications)
    standard_error_lower <- standard_error_upper <- numeric(replications)
    completed <- 0L
    while (completed < replications) {
      current <- min(batch_size, replications - completed)
      counts <- stats::rmultinom(
        current * folds, n / folds, probability
      )
      calculate_endpoint <- function(endpoint) {
        fold_loss <- array(
          t(counts) %*% endpoint$selection_loss,
          dim = c(folds, current, candidate_count)
        )
        fold_score <- array(
          t(counts) %*% endpoint$score,
          dim = c(folds, current, candidate_count)
        )
        fold_score_squared <- array(
          t(counts) %*% endpoint$score^2,
          dim = c(folds, current, candidate_count)
        )
        total_loss <- apply(fold_loss, c(2L, 3L), sum)
        selected <- matrix(NA_integer_, nrow = folds, ncol = current)
        for (fold in seq_len(folds)) {
          training_loss <- total_loss - matrix(
            fold_loss[fold, , ], nrow = current
          )
          selected[fold, ] <- max.col(
            -training_loss, ties.method = "first"
          )
        }
        signal_sum <- signal_square_sum <- numeric(current)
        replication_index <- seq_len(current)
        for (fold in seq_len(folds)) {
          local_score <- matrix(
            fold_score[fold, , ], nrow = current
          )
          local_square <- matrix(
            fold_score_squared[fold, , ], nrow = current
          )
          index <- cbind(replication_index, selected[fold, ])
          signal_sum <- signal_sum + local_score[index]
          signal_square_sum <- signal_square_sum + local_square[index]
        }
        estimate <- signal_sum / n
        standard_error <- sqrt(pmax(
          (signal_square_sum - n * estimate^2) / (n - 1), 0
        ) / n)
        cbind(
          estimate = estimate,
          standard_error = standard_error,
          z = (estimate - endpoint$target) / standard_error
        )
      }
      rows <- completed + seq_len(current)
      lower_result <- calculate_endpoint(lower)
      upper_result <- calculate_endpoint(upper)
      estimate_lower[rows] <- lower_result[, "estimate"]
      estimate_upper[rows] <- upper_result[, "estimate"]
      standard_error_lower[rows] <- lower_result[, "standard_error"]
      standard_error_upper[rows] <- upper_result[, "standard_error"]
      z_lower[rows] <- lower_result[, "z"]
      z_upper[rows] <- upper_result[, "z"]
      completed <- completed + current
    }
    one_sided <- stats::qnorm(0.95)
    bonferroni <- stats::qnorm(0.975)
    data.frame(
      sample_size = n,
      replications = replications,
      estimator = "oracle_five_fold_crossfit",
      lower_target = lower$target,
      upper_target = upper$target,
      lower_bias = mean(estimate_lower - lower$target),
      upper_bias = mean(estimate_upper - upper$target),
      lower_inward_bias = mean(estimate_lower - lower$target),
      upper_inward_bias = mean(upper$target - estimate_upper),
      lower_monte_carlo_sd = stats::sd(estimate_lower),
      upper_monte_carlo_sd = stats::sd(estimate_upper),
      lower_mean_standard_error = mean(standard_error_lower),
      upper_mean_standard_error = mean(standard_error_upper),
      lower_se_to_sd_ratio = mean(standard_error_lower) /
        stats::sd(estimate_lower),
      upper_se_to_sd_ratio = mean(standard_error_upper) /
        stats::sd(estimate_upper),
      lower_z_mean = mean(z_lower),
      lower_z_standard_deviation = stats::sd(z_lower),
      upper_z_mean = mean(z_upper),
      upper_z_standard_deviation = stats::sd(z_upper),
      lower_outward_coverage_95 = mean(z_lower <= one_sided),
      upper_outward_coverage_95 = mean(z_upper >= -one_sided),
      lower_one_sided_coverage_95 = mean(z_lower <= one_sided),
      upper_one_sided_coverage_95 = mean(z_upper >= -one_sided),
      lower_component_coverage_97_5 = mean(z_lower <= bonferroni),
      upper_component_coverage_97_5 = mean(z_upper >= -bonferroni),
      outward_set_coverage_pointwise_95 = mean(
        z_lower <= one_sided & z_upper >= -one_sided
      ),
      outward_set_coverage_bonferroni_95 = mean(
        z_lower <= bonferroni & z_upper >= -bonferroni
      ),
      lower_only_failure_bonferroni_95 = mean(
        z_lower > bonferroni & z_upper >= -bonferroni
      ),
      upper_only_failure_bonferroni_95 = mean(
        z_lower <= bonferroni & z_upper < -bonferroni
      ),
      both_fail_bonferroni_95 = mean(
        z_lower > bonferroni & z_upper < -bonferroni
      ),
      mean_point_set_width = mean(estimate_upper - estimate_lower),
      mean_outward_set_width = mean(
        estimate_upper + bonferroni * standard_error_upper -
          estimate_lower + bonferroni * standard_error_lower
      ),
      endpoint_z_correlation = stats::cor(z_lower, z_upper),
      stringsAsFactors = FALSE
    )
  }
  two_gap_bind_rows(lapply(seq_along(sample_sizes), function(index) {
    simulate_n(sample_sizes[index], seed + sample_sizes[index] + index)
  }))
}

simulate_selected_two_gap_online <- function(
    oracle = build_selected_two_gap_nonunique_threshold_oracle(),
    sample_sizes = c(1500L, 3000L, 6000L),
    replications = 1000L,
    initial_exponent = 0.7,
    chunks = 10L,
    seed = 1450000L
) {
  lower <- oracle$lower
  upper <- oracle$upper
  probability <- oracle$categories$probability
  sample_sizes <- as.integer(sample_sizes)
  replications <- as.integer(replications)
  chunks <- as.integer(chunks)
  simulate_n <- function(n, local_seed) {
    initial_size <- max(100L, ceiling(n^initial_exponent))
    evaluation_count <- n - initial_size
    chunk_sizes <- tabulate(
      sort(rep(seq_len(chunks), length.out = evaluation_count)),
      nbins = chunks
    )
    z_lower <- z_upper <- numeric(replications)
    estimate_lower <- estimate_upper <- numeric(replications)
    standard_error_lower <- standard_error_upper <- numeric(replications)
    set.seed(as.integer(local_seed))
    for (replication in seq_len(replications)) {
      initial_count <- as.vector(stats::rmultinom(
        1L, initial_size, probability
      ))
      chunk_count <- lapply(chunk_sizes, function(size) {
        as.vector(stats::rmultinom(1L, size, probability))
      })
      endpoint_result <- function(endpoint) {
        training_count <- initial_count
        numerator <- denominator <- 0
        for (chunk in seq_len(chunks)) {
          selected <- which.min(as.vector(
            training_count %*% endpoint$selection_loss
          ))
          standard_deviation <- endpoint$population_sd[selected]
          score_sum <- sum(
            chunk_count[[chunk]] * endpoint$score[, selected]
          )
          numerator <- numerator + score_sum / standard_deviation
          denominator <- denominator +
            chunk_sizes[chunk] / standard_deviation
          training_count <- training_count + chunk_count[[chunk]]
        }
        estimate <- numerator / denominator
        standard_error <- (evaluation_count / denominator) /
          sqrt(evaluation_count)
        c(
          estimate = estimate,
          standard_error = standard_error,
          z = (estimate - endpoint$target) / standard_error
        )
      }
      lower_result <- endpoint_result(lower)
      upper_result <- endpoint_result(upper)
      estimate_lower[replication] <- lower_result[["estimate"]]
      estimate_upper[replication] <- upper_result[["estimate"]]
      standard_error_lower[replication] <-
        lower_result[["standard_error"]]
      standard_error_upper[replication] <-
        upper_result[["standard_error"]]
      z_lower[replication] <- lower_result[["z"]]
      z_upper[replication] <- upper_result[["z"]]
    }
    one_sided <- stats::qnorm(0.95)
    bonferroni <- stats::qnorm(0.975)
    data.frame(
      sample_size = n,
      replications = replications,
      estimator = "oracle_online_predictable",
      initial_size = initial_size,
      chunks = chunks,
      lower_target = lower$target,
      upper_target = upper$target,
      lower_bias = mean(estimate_lower - lower$target),
      upper_bias = mean(estimate_upper - upper$target),
      lower_inward_bias = mean(estimate_lower - lower$target),
      upper_inward_bias = mean(upper$target - estimate_upper),
      lower_monte_carlo_sd = stats::sd(estimate_lower),
      upper_monte_carlo_sd = stats::sd(estimate_upper),
      lower_mean_standard_error = mean(standard_error_lower),
      upper_mean_standard_error = mean(standard_error_upper),
      lower_se_to_sd_ratio = mean(standard_error_lower) /
        stats::sd(estimate_lower),
      upper_se_to_sd_ratio = mean(standard_error_upper) /
        stats::sd(estimate_upper),
      lower_z_mean = mean(z_lower),
      lower_z_standard_deviation = stats::sd(z_lower),
      upper_z_mean = mean(z_upper),
      upper_z_standard_deviation = stats::sd(z_upper),
      lower_outward_coverage_95 = mean(z_lower <= one_sided),
      upper_outward_coverage_95 = mean(z_upper >= -one_sided),
      lower_one_sided_coverage_95 = mean(z_lower <= one_sided),
      upper_one_sided_coverage_95 = mean(z_upper >= -one_sided),
      lower_component_coverage_97_5 = mean(z_lower <= bonferroni),
      upper_component_coverage_97_5 = mean(z_upper >= -bonferroni),
      outward_set_coverage_pointwise_95 = mean(
        z_lower <= one_sided & z_upper >= -one_sided
      ),
      outward_set_coverage_bonferroni_95 = mean(
        z_lower <= bonferroni & z_upper >= -bonferroni
      ),
      lower_only_failure_bonferroni_95 = mean(
        z_lower > bonferroni & z_upper >= -bonferroni
      ),
      upper_only_failure_bonferroni_95 = mean(
        z_lower <= bonferroni & z_upper < -bonferroni
      ),
      both_fail_bonferroni_95 = mean(
        z_lower > bonferroni & z_upper < -bonferroni
      ),
      mean_point_set_width = mean(estimate_upper - estimate_lower),
      mean_outward_set_width = mean(
        estimate_upper + bonferroni * standard_error_upper -
          estimate_lower + bonferroni * standard_error_lower
      ),
      endpoint_z_correlation = stats::cor(z_lower, z_upper),
      stringsAsFactors = FALSE
    )
  }
  two_gap_bind_rows(lapply(seq_along(sample_sizes), function(index) {
    simulate_n(sample_sizes[index], seed + sample_sizes[index] + index)
  }))
}

generate_selected_two_gap_data <- function(
    n,
    oracle = build_selected_two_gap_nonunique_threshold_oracle(),
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

compute_selected_two_gap_tan_population <- function(
    oracle = build_selected_two_gap_nonunique_threshold_oracle()
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
  tan_x <- observed_dgp_standard_msm_arm(
    specification$atoms, marginal_probability, marginal_pi, 1L, Lambda
  )
  tan_xz <-
    e * observed_dgp_standard_msm_arm(
      specification$atoms, probability0, pi0, 1L, Lambda
    ) +
    (1 - e) * observed_dgp_standard_msm_arm(
      specification$atoms, probability1, pi1, 1L, Lambda
    )
  data.frame(
    estimator = c("tan_cal_x_intercept", "tan_cal_xz"),
    Lambda = Lambda,
    lower_target = c(tan_x[["lower"]], tan_xz[["lower"]]),
    upper_target = c(tan_x[["upper"]], tan_xz[["upper"]]),
    stringsAsFactors = FALSE
  )
}

simulate_selected_two_gap_tan <- function(
    oracle = build_selected_two_gap_nonunique_threshold_oracle(),
    sample_sizes = c(1500L, 3000L, 6000L),
    replications = 1000L,
    seed = 1500000L,
    probability_epsilon = 0.01
) {
  if (!exists("estimate_tan_cal_arm_endpoint", mode = "function")) {
    if (!file.exists("R/comparator_estimators.R")) {
      stop("Run this function from the repository root.")
    }
    source("R/comparator_estimators.R")
  }
  sample_sizes <- as.integer(sample_sizes)
  replications <- as.integer(replications)
  if (!length(sample_sizes) || any(sample_sizes < 10L) || replications < 1L) {
    stop("Invalid selected two-gap Tan simulation inputs.")
  }
  population <- compute_selected_two_gap_tan_population(oracle)
  Lambda <- oracle$lower$specification$sensitivity_lambda
  configurations <- list(
    tan_cal_x_intercept = character(),
    tan_cal_xz = "Z"
  )
  results <- list()
  result_index <- 0L
  for (sample_index in seq_along(sample_sizes)) {
    n <- sample_sizes[sample_index]
    z_values <- lapply(configurations, function(ignored) {
      list(
        lower = rep(NA_real_, replications),
        upper = rep(NA_real_, replications)
      )
    })
    successful <- setNames(integer(length(configurations)), names(configurations))
    for (replication in seq_len(replications)) {
      data <- generate_selected_two_gap_data(
        n, oracle, seed + n + sample_index + replication
      )
      arm <- as.numeric(data$T_obs == 1L)
      for (estimator in names(configurations)) {
        design <- comparator_design_matrix(
          data, configurations[[estimator]], degree = 1L
        )
        fits <- tryCatch(
          list(
            lower = estimate_tan_cal_arm_endpoint(
              data$Y, arm, design, Lambda, "lower", probability_epsilon
            ),
            upper = estimate_tan_cal_arm_endpoint(
              data$Y, arm, design, Lambda, "upper", probability_epsilon
            )
          ),
          error = identity
        )
        if (inherits(fits, "error")) next
        target <- population[population$estimator == estimator, , drop = FALSE]
        z_values[[estimator]]$lower[replication] <-
          (fits$lower$estimate - target$lower_target) /
          fits$lower$standard_error
        z_values[[estimator]]$upper[replication] <-
          (fits$upper$estimate - target$upper_target) /
          fits$upper$standard_error
        successful[[estimator]] <- successful[[estimator]] + 1L
      }
    }
    for (estimator in names(configurations)) {
      complete <- is.finite(z_values[[estimator]]$lower) &
        is.finite(z_values[[estimator]]$upper)
      success <- sum(complete)
      if (success < 2L) {
        stop("Fewer than two successful Tan replications for ", estimator, ".")
      }
      lower <- z_values[[estimator]]$lower[complete]
      upper <- z_values[[estimator]]$upper[complete]
      one_sided <- stats::qnorm(0.95)
      bonferroni <- stats::qnorm(0.975)
      result_index <- result_index + 1L
      results[[result_index]] <- data.frame(
        sample_size = n,
        requested_replications = replications,
        successful_replications = success,
        estimator = estimator,
        lower_target = population$lower_target[
          population$estimator == estimator
        ],
        upper_target = population$upper_target[
          population$estimator == estimator
        ],
        lower_z_mean = mean(lower),
        lower_z_standard_deviation = stats::sd(lower),
        upper_z_mean = mean(upper),
        upper_z_standard_deviation = stats::sd(upper),
        lower_outward_coverage_95 = mean(lower <= one_sided),
        upper_outward_coverage_95 = mean(upper >= -one_sided),
        outward_set_coverage_pointwise_95 = mean(
          lower <= one_sided & upper >= -one_sided
        ),
        outward_set_coverage_bonferroni_95 = mean(
          lower <= bonferroni & upper >= -bonferroni
        ),
        endpoint_z_correlation = stats::cor(lower, upper),
        stringsAsFactors = FALSE
      )
    }
  }
  two_gap_bind_rows(results)
}

simulate_selected_two_gap_tan_multinomial <- function(
    oracle = build_selected_two_gap_nonunique_threshold_oracle(),
    sample_sizes = c(1500L, 3000L, 6000L),
    replications = 10000L,
    seed = 1550000L
) {
  sample_sizes <- as.integer(sample_sizes)
  replications <- as.integer(replications)
  if (!length(sample_sizes) || any(sample_sizes < 10L) || replications < 1L) {
    stop("Invalid selected two-gap multinomial Tan simulation inputs.")
  }
  specification <- oracle$lower$specification
  Lambda <- specification$sensitivity_lambda
  multiplier <- Lambda - 1 / Lambda
  population <- compute_selected_two_gap_tan_population(oracle)
  collapsed <- stats::aggregate(
    probability ~ Z + T_obs + atom + Y,
    oracle$categories,
    sum
  )
  collapsed <- collapsed[order(
    collapsed$Z, collapsed$T_obs, collapsed$atom
  ), , drop = FALSE]

  endpoint_fit <- function(count, adjustment, direction) {
    sign <- if (direction == "lower") -1 else 1
    tau <- if (direction == "lower") {
      1 / (Lambda + 1)
    } else {
      Lambda / (Lambda + 1)
    }
    group <- if (adjustment == "tan_cal_x_intercept") {
      rep(1L, nrow(collapsed))
    } else {
      collapsed$Z + 1L
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
      treated_by_atom <- tapply(
        count[treated], collapsed$atom[treated], sum
      )
      treated_by_atom <- treated_by_atom[
        as.character(seq_along(specification$atoms))
      ]
      treated_by_atom[is.na(treated_by_atom)] <- 0
      quantile_index <- which(
        cumsum(treated_by_atom) / treated_count >= tau
      )[1L]
      q <- specification$atoms[quantile_index]
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
    z <- lapply(population$estimator, function(ignored) {
      list(
        lower = rep(NA_real_, replications),
        upper = rep(NA_real_, replications)
      )
    })
    names(z) <- population$estimator
    set.seed(seed + n + sample_index)
    counts <- stats::rmultinom(
      replications, n, collapsed$probability
    )
    for (replication in seq_len(replications)) {
      for (estimator in population$estimator) {
        target <- population[
          population$estimator == estimator, , drop = FALSE
        ]
        for (direction in c("lower", "upper")) {
          fit <- endpoint_fit(counts[, replication], estimator, direction)
          z[[estimator]][[direction]][replication] <-
            (fit[["estimate"]] - target[[paste0(direction, "_target")]]) /
            fit[["standard_error"]]
        }
      }
    }
    for (estimator in population$estimator) {
      lower <- z[[estimator]]$lower
      upper <- z[[estimator]]$upper
      complete <- is.finite(lower) & is.finite(upper)
      lower <- lower[complete]
      upper <- upper[complete]
      pointwise <- stats::qnorm(0.95)
      simultaneous <- stats::qnorm(0.975)
      output_index <- output_index + 1L
      output[[output_index]] <- data.frame(
        sample_size = n,
        requested_replications = replications,
        successful_replications = sum(complete),
        estimator = estimator,
        lower_z_mean = mean(lower),
        lower_z_standard_deviation = stats::sd(lower),
        upper_z_mean = mean(upper),
        upper_z_standard_deviation = stats::sd(upper),
        outward_set_coverage_pointwise_95 = mean(
          lower <= pointwise & upper >= -pointwise
        ),
        outward_set_coverage_bonferroni_95 = mean(
          lower <= simultaneous & upper >= -simultaneous
        ),
        endpoint_z_correlation = stats::cor(lower, upper),
        stringsAsFactors = FALSE
      )
    }
  }
  two_gap_bind_rows(output)
}
