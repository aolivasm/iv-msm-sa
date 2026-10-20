# Cross-fitted direct-RKHS density-ratio nuisances for binary outcomes.
#
# Required sources:
#   source("R/observed_binary_dgp.R")
#   source("R/sharp_bounds_plugin.R")
#
# The direct density-ratio fit targets
#
#   gamma_t(y,x) = f(y | T=t,X=x,Z=1) / f(y | T=t,X=x,Z=0)
#
# with the same reduced Riesz/RKHS system used by the continuous-outcome
# estimator.  For binary Y, a baseline probability a_t0 and the two fitted
# density-ratio values determine a normalized comparison probability:
#
#                   gamma_t(1,x) a_t0(x)
#   a_t1(x) = ------------------------------------------------.
#             gamma_t(1,x) a_t0(x) + gamma_t(0,x){1-a_t0(x)}
#
# This file deliberately does not use the discarded classification odds-ratio
# estimator.

binary_rkhs_gamma_requirements <- function(require_oracle = FALSE) {
  required <- c(
    "binary_fit_regression_with_learner",
    "binary_make_tz_folds",
    "binary_clip_probability",
    "default_gamma_lambda_grid",
    "fit_gamma_rkhs",
    "fit_gamma_rkhs_fixed"
  )
  if (isTRUE(require_oracle)) {
    required <- c(
      required,
      "observed_dgp_nuisances",
      "observed_binary_probabilities",
      "observed_binary_gamma"
    )
  }
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing)) {
    stop(
      "Source R/observed_binary_dgp.R and R/sharp_bounds_plugin.R first; ",
      "missing functions: ", paste(missing, collapse = ", "), "."
    )
  }
  invisible(TRUE)
}

binary_rkhs_gamma_oracle_predictions <- function(
    data,
    specification,
    x_cols
) {
  binary_rkhs_gamma_requirements(require_oracle = TRUE)
  x <- data[, x_cols, drop = FALSE]
  design <- observed_dgp_nuisances(x, specification)
  output <- data.frame(
    e = design$e,
    pi0 = design$pi0,
    pi1 = design$pi1,
    stringsAsFactors = FALSE
  )
  for (treatment in 0:1) {
    outcome <- observed_binary_probabilities(x, treatment, specification)
    output[[sprintf("a_t%d_z0", treatment)]] <- outcome$a0
    output[[sprintf("a_t%d_z1", treatment)]] <- outcome$a1
  }
  output
}

# Legacy fixed-tuning rule retained for reproducibility and fast debugging.
# The production default below fixes sigma^2 at one half of the median squared
# pairwise distance and tunes lambda over the generic C/sqrt(n_t) grid.
make_binary_rkhs_sqrt_n_lambda_rule <- function(
    reference_lambda_by_t = c(`0` = 0.03162278, `1` = 0.01),
    reference_training_size_by_t = c(`0` = 500, `1` = 700)
) {
  if (
    is.null(names(reference_lambda_by_t)) ||
      is.null(names(reference_training_size_by_t)) ||
      !all(c("0", "1") %in% names(reference_lambda_by_t)) ||
      !all(c("0", "1") %in% names(reference_training_size_by_t))
  ) {
    stop("Both RKHS reference vectors must have entries named '0' and '1'.")
  }
  if (
    any(!is.finite(reference_lambda_by_t[c("0", "1")])) ||
      any(reference_lambda_by_t[c("0", "1")] <= 0) ||
      any(!is.finite(reference_training_size_by_t[c("0", "1")])) ||
      any(reference_training_size_by_t[c("0", "1")] <= 0)
  ) {
    stop("RKHS reference penalties and training sizes must be positive.")
  }
  function(n_training, treatment) {
    key <- as.character(as.integer(treatment))
    if (!key %in% c("0", "1")) stop("treatment must be zero or one.")
    if (length(n_training) != 1L || !is.finite(n_training) || n_training < 1) {
      stop("n_training must be one positive number.")
    }
    unname(
      reference_lambda_by_t[key] *
        sqrt(reference_training_size_by_t[key] / n_training)
    )
  }
}

binary_resolve_rkhs_lambda <- function(
    rkhs_lambda_by_t,
    treatment,
    n_training
) {
  key <- as.character(as.integer(treatment))
  value <- if (is.function(rkhs_lambda_by_t)) {
    rkhs_lambda_by_t(
      n_training = as.integer(n_training),
      treatment = as.integer(treatment)
    )
  } else {
    if (
      is.null(names(rkhs_lambda_by_t)) ||
        !all(c("0", "1") %in% names(rkhs_lambda_by_t))
    ) {
      stop("rkhs_lambda_by_t must be a function or have entries '0' and '1'.")
    }
    unname(rkhs_lambda_by_t[key])
  }
  if (length(value) != 1L || !is.finite(value) || value <= 0) {
    stop("The resolved binary RKHS penalty must be positive and finite.")
  }
  as.numeric(value)
}

binary_rkhs_gamma_reconstruct <- function(
    a0,
    gamma_zero,
    gamma_one,
    probability_epsilon = 1e-4,
    gamma_floor = 1e-6
) {
  a0 <- as.numeric(a0)
  gamma_zero <- pmax(as.numeric(gamma_zero), gamma_floor)
  gamma_one <- pmax(as.numeric(gamma_one), gamma_floor)
  if (
    length(a0) != length(gamma_zero) ||
      length(a0) != length(gamma_one) ||
      any(!is.finite(c(a0, gamma_zero, gamma_one)))
  ) {
    stop("Invalid binary RKHS gamma reconstruction inputs.")
  }
  unnormalized_zero <- gamma_zero * (1 - a0)
  unnormalized_one <- gamma_one * a0
  normalizer <- unnormalized_zero + unnormalized_one
  if (any(!is.finite(normalizer)) || any(normalizer <= 0)) {
    stop("The binary RKHS gamma normalizer is nonpositive or non-finite.")
  }
  raw_a1 <- unnormalized_one / normalizer
  a1 <- binary_clip_probability(raw_a1, probability_epsilon)
  list(
    a1 = a1,
    raw_a1 = raw_a1,
    normalizer = normalizer,
    normalized_gamma_zero = gamma_zero / normalizer,
    normalized_gamma_one = gamma_one / normalizer,
    lower_probability_clipping_rate = mean(raw_a1 < probability_epsilon),
    upper_probability_clipping_rate = mean(
      raw_a1 > 1 - probability_epsilon
    )
  )
}

binary_rkhs_gamma_fit_component <- function(
    training,
    response,
    x_cols,
    seed,
    spline_df = 5L,
    gam_min_basis = 5L,
    gam_max_basis = 20L,
    gam_growth_exponent = 1 / 3,
    library_candidates = c("glm_spline", "gam", "glmnet", "ranger"),
    library_cv_folds = 3L,
    library_selection_rule = "one_se",
    library_ranger_num_trees = 500L,
    learner = "library_cv",
    super_learner_candidates = binary_default_super_learner_candidates(),
    super_learner_cv_folds = 3L,
    super_learner_ranger_num_trees = 500L,
    super_learner_hal_max_degree = 2L
) {
  binary_fit_regression_with_learner(
    data = training,
    response = response,
    x_cols = x_cols,
    spline_df = spline_df,
    learner = learner,
    gam_min_basis = gam_min_basis,
    gam_max_basis = gam_max_basis,
    gam_growth_exponent = gam_growth_exponent,
    library_candidates = library_candidates,
    library_cv_folds = library_cv_folds,
    library_selection_rule = library_selection_rule,
    ranger_num_trees = library_ranger_num_trees,
    super_learner_candidates = super_learner_candidates,
    super_learner_cv_folds = super_learner_cv_folds,
    super_learner_ranger_num_trees = super_learner_ranger_num_trees,
    super_learner_hal_max_degree = super_learner_hal_max_degree,
    seed = seed
  )
}

# Production nuisance fit for binary outcomes.  Unlike the paired simulation
# audit below, this function contains no oracle-DGP calls and is therefore safe
# to use in real-data analyses.  It fits only the primitive regressions needed
# by the direct-gamma parameterization: e, pi_0, pi_1, and a_{t0}.  The
# comparison regression a_{t1} is reconstructed from a_{t0} and the direct
# density-ratio fit rather than estimated as a separate probability model.
fit_binary_cross_fitted_direct_rkhs_nuisances <- function(
    data,
    x_cols,
    K = 5L,
    fold_id = NULL,
    seed = 1L,
    overlap_epsilon = 0.05,
    outcome_probability_epsilon = 1e-4,
    spline_df = 5L,
    gam_min_basis = 5L,
    gam_max_basis = 20L,
    gam_growth_exponent = 1 / 3,
    primitive_learner = c(
      "super_learner", "library_cv", "gam", "ranger", "glm_spline",
      "glmnet"
    ),
    library_candidates = c("glm_spline", "gam", "glmnet", "ranger"),
    library_cv_folds = 3L,
    library_selection_rule = c("one_se", "minimum"),
    library_ranger_num_trees = 500L,
    super_learner_candidates = binary_default_super_learner_candidates(),
    super_learner_cv_folds = 3L,
    super_learner_ranger_num_trees = 500L,
    super_learner_hal_max_degree = 2L,
    rkhs_lambda_by_t = NULL,
    rkhs_lambda_grid = default_gamma_lambda_grid,
    rkhs_bandwidth_multipliers = 1 / sqrt(2),
    rkhs_cv_folds = 3L,
    rkhs_cv_repeats = 2L,
    rkhs_selection_rule = c("minimum", "one_se"),
    rkhs_bandwidth_multiplier = 1 / sqrt(2),
    rkhs_gamma_floor = 1e-6,
    rkhs_gamma_cap = 100,
    rkhs_backend = c("nystrom", "exact"),
    rkhs_nystrom_rank = NULL,
    rkhs_kernel = c("gaussian", "sobolev_matern32", "sobolev_matern52"),
    rkhs_bandwidth_sample = 1000L
) {
  binary_rkhs_gamma_requirements()
  primitive_learner <- match.arg(primitive_learner)
  library_selection_rule <- match.arg(library_selection_rule)
  rkhs_selection_rule <- match.arg(rkhs_selection_rule)
  rkhs_backend <- match.arg(rkhs_backend)
  rkhs_kernel <- match.arg(rkhs_kernel)
  gamma_tuning_mode <- if (is.null(rkhs_lambda_by_t)) "joint_cv" else "fixed"
  if (
    gamma_tuning_mode == "joint_cv" &&
      (length(rkhs_cv_folds) != 1L || !is.finite(rkhs_cv_folds) ||
        rkhs_cv_folds < 2L || rkhs_cv_folds != as.integer(rkhs_cv_folds))
  ) {
    stop("rkhs_cv_folds must be one integer of at least two.")
  }
  if (
    gamma_tuning_mode == "joint_cv" &&
      (length(rkhs_cv_repeats) != 1L || !is.finite(rkhs_cv_repeats) ||
        rkhs_cv_repeats < 1L ||
        rkhs_cv_repeats != as.integer(rkhs_cv_repeats))
  ) {
    stop("rkhs_cv_repeats must be one positive integer.")
  }
  required <- c("Y", "T_obs", "Z", x_cols)
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  if (is.null(fold_id)) {
    fold_id <- binary_make_tz_folds(data, K = K, seed = seed)
  }
  if (length(fold_id) != nrow(data)) stop("fold_id has the wrong length.")
  folds <- sort(unique(fold_id))
  if (length(folds) < 2L) stop("At least two outer folds are required.")

  prediction_names <- c(
    "e", "pi0", "pi1",
    "a_t0_z0", "a_t0_z1", "a_t1_z0", "a_t1_z1"
  )
  predictions <- as.data.frame(setNames(
    replicate(
      length(prediction_names), rep(NA_real_, nrow(data)), simplify = FALSE
    ),
    prediction_names
  ))
  fit_diagnostics <- gamma_diagnostics <- reconstruction_diagnostics <- list()
  cv_tables <- gamma_cv_tables <- gamma_cv_resample_tables <- list()
  fit_counter <- gamma_counter <- reconstruction_counter <- 0L

  fit_component <- function(training, response, local_seed) {
    binary_rkhs_gamma_fit_component(
      training = training,
      response = response,
      x_cols = x_cols,
      seed = local_seed,
      spline_df = spline_df,
      gam_min_basis = gam_min_basis,
      gam_max_basis = gam_max_basis,
      gam_growth_exponent = gam_growth_exponent,
      library_candidates = library_candidates,
      library_cv_folds = library_cv_folds,
      library_selection_rule = library_selection_rule,
      library_ranger_num_trees = library_ranger_num_trees,
      learner = primitive_learner,
      super_learner_candidates = super_learner_candidates,
      super_learner_cv_folds = super_learner_cv_folds,
      super_learner_ranger_num_trees = super_learner_ranger_num_trees,
      super_learner_hal_max_degree = super_learner_hal_max_degree
    )
  }
  clip_design <- function(value) {
    binary_clip_probability(value, overlap_epsilon)
  }
  clip_outcome <- function(value) {
    binary_clip_probability(value, outcome_probability_epsilon)
  }

  for (fold in folds) {
    training <- data[fold_id != fold, , drop = FALSE]
    validation <- data[fold_id == fold, , drop = FALSE]
    validation_index <- which(fold_id == fold)

    fits <- list(
      e = fit_component(training, "Z", seed + 1000L * fold + 1L),
      pi0 = fit_component(
        training[training$Z == 0, , drop = FALSE],
        "T_obs", seed + 1000L * fold + 2L
      ),
      pi1 = fit_component(
        training[training$Z == 1, , drop = FALSE],
        "T_obs", seed + 1000L * fold + 3L
      )
    )
    for (treatment in 0:1) {
      baseline_name <- sprintf("a_t%d_z0", treatment)
      cell <- training$T_obs == treatment & training$Z == 0
      fits[[baseline_name]] <- fit_component(
        training[cell, , drop = FALSE],
        "Y",
        seed + 1000L * fold + 100L * treatment + 4L
      )
    }

    raw_training <- list(
      e = 1 - fits$e$predict(training),
      pi0 = fits$pi0$predict(training),
      pi1 = fits$pi1$predict(training)
    )
    raw_validation <- list(
      e = 1 - fits$e$predict(validation),
      pi0 = fits$pi0$predict(validation),
      pi1 = fits$pi1$predict(validation)
    )
    for (treatment in 0:1) {
      baseline_name <- sprintf("a_t%d_z0", treatment)
      raw_training[[baseline_name]] <- fits[[baseline_name]]$predict(training)
      raw_validation[[baseline_name]] <- fits[[baseline_name]]$predict(validation)
    }
    estimated_training <- data.frame(
      e = clip_design(raw_training$e),
      pi0 = clip_design(raw_training$pi0),
      pi1 = clip_design(raw_training$pi1),
      stringsAsFactors = FALSE
    )
    estimated_validation <- data.frame(
      e = clip_design(raw_validation$e),
      pi0 = clip_design(raw_validation$pi0),
      pi1 = clip_design(raw_validation$pi1),
      stringsAsFactors = FALSE
    )
    for (treatment in 0:1) {
      baseline_name <- sprintf("a_t%d_z0", treatment)
      estimated_training[[baseline_name]] <- clip_outcome(
        raw_training[[baseline_name]]
      )
      estimated_validation[[baseline_name]] <- clip_outcome(
        raw_validation[[baseline_name]]
      )
    }

    predictions[validation_index, c("e", "pi0", "pi1")] <-
      estimated_validation[c("e", "pi0", "pi1")]
    for (treatment in 0:1) {
      baseline_name <- sprintf("a_t%d_z0", treatment)
      comparison_name <- sprintf("a_t%d_z1", treatment)
      n_training_t <- sum(training$T_obs == treatment)
      gamma_seed <- seed + 100000L * fold + 1000L * treatment
      gamma_fit <- if (gamma_tuning_mode == "joint_cv") {
        fit_gamma_rkhs(
          train = training,
          nuisance_train = estimated_training,
          t_value = treatment,
          x_cols = x_cols,
          lambda_grid = rkhs_lambda_grid,
          bandwidth_multipliers = rkhs_bandwidth_multipliers,
          cv_folds = rkhs_cv_folds,
          cv_repeats = rkhs_cv_repeats,
          selection_rule = rkhs_selection_rule,
          gamma_floor = rkhs_gamma_floor,
          gamma_cap = rkhs_gamma_cap,
          bandwidth_sample = rkhs_bandwidth_sample,
          backend = rkhs_backend,
          nystrom_rank = rkhs_nystrom_rank,
          kernel_type = rkhs_kernel,
          seed = gamma_seed
        )
      } else {
        selected_lambda <- binary_resolve_rkhs_lambda(
          rkhs_lambda_by_t = rkhs_lambda_by_t,
          treatment = treatment,
          n_training = n_training_t
        )
        fit_gamma_rkhs_fixed(
          train = training,
          nuisance_train = estimated_training,
          t_value = treatment,
          x_cols = x_cols,
          lambda = selected_lambda,
          bandwidth_multiplier = rkhs_bandwidth_multiplier,
          gamma_floor = rkhs_gamma_floor,
          gamma_cap = rkhs_gamma_cap,
          bandwidth_sample = rkhs_bandwidth_sample,
          backend = rkhs_backend,
          nystrom_rank = rkhs_nystrom_rank,
          kernel_type = rkhs_kernel,
          seed = gamma_seed
        )
      }
      if (!is.null(gamma_fit$cv_table) && nrow(gamma_fit$cv_table)) {
        table <- gamma_fit$cv_table
        table$outer_fold <- fold
        table$treatment <- treatment
        table$n_training_t <- n_training_t
        table$lambda_sqrt_n_constant <- table$lambda * sqrt(n_training_t)
        table$selection_rule <- rkhs_selection_rule
        gamma_cv_tables[[length(gamma_cv_tables) + 1L]] <- table
      }
      if (
        !is.null(gamma_fit$cv_resample_table) &&
          nrow(gamma_fit$cv_resample_table)
      ) {
        resample_table <- gamma_fit$cv_resample_table
        resample_table$outer_fold <- fold
        resample_table$treatment <- treatment
        resample_table$n_training_t <- n_training_t
        resample_table$lambda_sqrt_n_constant <-
          resample_table$lambda * sqrt(n_training_t)
        resample_table$selection_rule <- rkhs_selection_rule
        gamma_cv_resample_tables[[
          length(gamma_cv_resample_tables) + 1L
        ]] <- resample_table
      }
      new_zero <- new_one <- validation
      new_zero$Y <- 0
      new_one$Y <- 1
      prediction_zero <- gamma_fit$predict(new_zero)
      prediction_one <- gamma_fit$predict(new_one)
      reconstructed <- binary_rkhs_gamma_reconstruct(
        a0 = estimated_validation[[baseline_name]],
        gamma_zero = prediction_zero$normalized,
        gamma_one = prediction_one$normalized,
        probability_epsilon = outcome_probability_epsilon,
        gamma_floor = rkhs_gamma_floor
      )
      predictions[[baseline_name]][validation_index] <-
        estimated_validation[[baseline_name]]
      predictions[[comparison_name]][validation_index] <- reconstructed$a1

      gamma_counter <- gamma_counter + 1L
      gamma_diagnostics[[gamma_counter]] <- data.frame(
        fold = fold,
        treatment = treatment,
        n_training_t = n_training_t,
        tuning_mode = gamma_tuning_mode,
        cv_folds = if (gamma_tuning_mode == "joint_cv") {
          as.integer(rkhs_cv_folds)
        } else NA_integer_,
        cv_repeats = if (gamma_tuning_mode == "joint_cv") {
          as.integer(rkhs_cv_repeats)
        } else NA_integer_,
        selected_lambda = gamma_fit$selected_lambda,
        selected_lambda_sqrt_n_constant =
          gamma_fit$selected_lambda * sqrt(n_training_t),
        selected_bandwidth_multiplier =
          gamma_fit$selected_bandwidth_multiplier,
        base_bandwidth = gamma_fit$base_bandwidth,
        selected_bandwidth = gamma_fit$selected_bandwidth,
        normalization_constant = gamma_fit$normalization_constant,
        training_negative_rate = gamma_fit$training_negative_rate,
        rkhs_backend = gamma_fit$rkhs_backend,
        rkhs_kernel = gamma_fit$rkhs_kernel,
        nystrom_rank_requested = gamma_fit$nystrom_rank_requested,
        nystrom_landmark_count = gamma_fit$nystrom_landmark_count,
        nystrom_effective_rank = gamma_fit$nystrom_effective_rank,
        stringsAsFactors = FALSE
      )
      reconstruction_counter <- reconstruction_counter + 1L
      reconstruction_diagnostics[[reconstruction_counter]] <- data.frame(
        fold = fold,
        treatment = treatment,
        mean_absolute_pre_normalization_error = mean(abs(
          reconstructed$normalizer - 1
        )),
        maximum_absolute_pre_normalization_error = max(abs(
          reconstructed$normalizer - 1
        )),
        lower_probability_clipping_rate =
          reconstructed$lower_probability_clipping_rate,
        upper_probability_clipping_rate =
          reconstructed$upper_probability_clipping_rate,
        stringsAsFactors = FALSE
      )
    }

    component_epsilon <- c(
      e = overlap_epsilon,
      pi0 = overlap_epsilon,
      pi1 = overlap_epsilon,
      a_t0_z0 = outcome_probability_epsilon,
      a_t1_z0 = outcome_probability_epsilon
    )
    for (component in names(fits)) {
      fit_counter <- fit_counter + 1L
      fit <- fits[[component]]
      fit_diagnostics[[fit_counter]] <- data.frame(
        fold = fold,
        component = component,
        learner = fit$learner,
        selected_learner = if (is.null(fit$selected_learner)) {
          fit$learner
        } else {
          fit$selected_learner
        },
        converged = isTRUE(fit$converged),
        rank = if (is.null(fit$rank)) NA_integer_ else as.integer(fit$rank),
        constant_fit = isTRUE(fit$constant),
        basis_dimension = if (is.null(fit$basis_dimension)) {
          NA_integer_
        } else {
          as.integer(fit$basis_dimension)
        },
        effective_df = if (is.null(fit$effective_df)) {
          NA_real_
        } else {
          as.numeric(fit$effective_df)
        },
        cv_mean_risk = if (is.null(fit$cv_mean_risk)) {
          NA_real_
        } else {
          as.numeric(fit$cv_mean_risk)
        },
        truncation_epsilon = unname(component_epsilon[component]),
        lower_truncation_rate = mean(
          raw_validation[[component]] < component_epsilon[[component]]
        ),
        upper_truncation_rate = mean(
          raw_validation[[component]] > 1 - component_epsilon[[component]]
        ),
        stringsAsFactors = FALSE
      )
      if (!is.null(fit$cv_table) && nrow(fit$cv_table)) {
        table <- fit$cv_table
        table$outer_fold <- fold
        table$component <- component
        cv_tables[[length(cv_tables) + 1L]] <- table
      }
    }
  }

  if (any(!is.finite(as.matrix(predictions)))) {
    stop("At least one direct-RKHS binary nuisance prediction is non-finite.")
  }
  list(
    predictions = predictions,
    fold_id = fold_id,
    learner = "direct_rkhs_gamma",
    diagnostics = do.call(rbind, fit_diagnostics),
    cv_table = if (length(cv_tables)) do.call(rbind, cv_tables) else data.frame(),
    gamma_cv_table = if (length(gamma_cv_tables)) {
      do.call(rbind, gamma_cv_tables)
    } else {
      data.frame()
    },
    gamma_cv_resample_table = if (length(gamma_cv_resample_tables)) {
      do.call(rbind, gamma_cv_resample_tables)
    } else {
      data.frame()
    },
    gamma_diagnostics = do.call(rbind, gamma_diagnostics),
    reconstruction_diagnostics = do.call(
      rbind, reconstruction_diagnostics
    ),
    configuration = list(
      K = length(folds),
      primitive_learner = primitive_learner,
      overlap_epsilon = overlap_epsilon,
      outcome_probability_epsilon = outcome_probability_epsilon,
      library_candidates = library_candidates,
      library_cv_folds = as.integer(library_cv_folds),
      library_selection_rule = library_selection_rule,
      super_learner_candidates = super_learner_candidates,
      super_learner_cv_folds = as.integer(super_learner_cv_folds),
      super_learner_ranger_num_trees = as.integer(
        super_learner_ranger_num_trees
      ),
      super_learner_hal_max_degree = as.integer(
        super_learner_hal_max_degree
      ),
      rkhs_tuning_mode = gamma_tuning_mode,
      rkhs_lambda_specification = if (gamma_tuning_mode == "joint_cv") {
        if (identical(rkhs_lambda_grid, default_gamma_lambda_grid)) {
          "C/sqrt(n_t)"
        } else if (is.function(rkhs_lambda_grid)) {
          "user_function_of_treatment_arm_training_size"
        } else {
          "user_candidate_grid"
        }
      } else if (is.function(rkhs_lambda_by_t)) {
        "fixed_function_of_treatment_arm_training_size"
      } else {
        rkhs_lambda_by_t
      },
      rkhs_lambda_grid_at_n1 = if (gamma_tuning_mode == "joint_cv") {
        resolve_gamma_lambda_grid(
          rkhs_lambda_grid, n_training = 1L, treatment = 0L
        )
      } else {
        NULL
      },
      rkhs_bandwidth_multipliers = if (gamma_tuning_mode == "joint_cv") {
        rkhs_bandwidth_multipliers
      } else {
        rkhs_bandwidth_multiplier
      },
      rkhs_cv_folds = if (gamma_tuning_mode == "joint_cv") {
        as.integer(rkhs_cv_folds)
      } else {
        NA_integer_
      },
      rkhs_cv_repeats = if (gamma_tuning_mode == "joint_cv") {
        as.integer(rkhs_cv_repeats)
      } else {
        NA_integer_
      },
      rkhs_selection_rule = if (gamma_tuning_mode == "joint_cv") {
        rkhs_selection_rule
      } else {
        NA_character_
      },
      rkhs_gamma_floor = rkhs_gamma_floor,
      rkhs_gamma_cap = rkhs_gamma_cap,
      rkhs_backend = rkhs_backend,
      rkhs_nystrom_rank = rkhs_nystrom_rank,
      rkhs_kernel = rkhs_kernel
    )
  )
}

# Public production wrapper.  The legacy estimate_binary_dual_bounds() in
# binary_outcome_bounds.R retains the separate-probability parameterization as
# an explicit comparator; new binary analyses should call this function.
estimate_binary_dual_bounds_rkhs <- function(
    data,
    x_cols = "X",
    Lambda0 = 2,
    Lambda1 = Lambda0,
    epsilon = 0.01,
    K = 5L,
    fold_id = NULL,
    confidence_level = 0.95,
    compatibility_action = c("project", "none"),
    compatibility_tolerance = 1e-12,
    seed = 1L,
    ...
) {
  compatibility_action <- match.arg(compatibility_action)
  nuisance_fit <- fit_binary_cross_fitted_direct_rkhs_nuisances(
    data = data,
    x_cols = x_cols,
    K = K,
    fold_id = fold_id,
    seed = seed,
    ...
  )
  result <- evaluate_binary_bounds_from_nuisances(
    data = data,
    nuisance_fit = nuisance_fit,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    epsilon = epsilon,
    confidence_level = confidence_level,
    compatibility_action = compatibility_action,
    compatibility_tolerance = compatibility_tolerance
  )
  result$nuisance_fit <- nuisance_fit
  result
}

# Once this production source is loaded, keep the former public implementation
# under an explicit comparator name and make the direct-RKHS wrapper the public
# default.  The guard makes repeated source() calls idempotent.
if (
  exists("estimate_binary_dual_bounds", mode = "function") &&
    !exists(
      "estimate_binary_dual_bounds_probability_regression",
      mode = "function"
    )
) {
  estimate_binary_dual_bounds_probability_regression <-
    estimate_binary_dual_bounds
}
estimate_binary_dual_bounds <- estimate_binary_dual_bounds_rkhs

fit_binary_cross_fitted_rkhs_gamma_nuisances <- function(
    data,
    specification,
    x_cols,
    K = 5L,
    fold_id = NULL,
    seed = 1L,
    overlap_epsilon = 0.05,
    outcome_probability_epsilon = 1e-4,
    spline_df = 5L,
    gam_min_basis = 5L,
    gam_max_basis = 20L,
    gam_growth_exponent = 1 / 3,
    library_candidates = c("glm_spline", "gam", "glmnet", "ranger"),
    library_cv_folds = 3L,
    library_selection_rule = "one_se",
    library_ranger_num_trees = 500L,
    rkhs_lambda_by_t = NULL,
    rkhs_lambda_grid = default_gamma_lambda_grid,
    rkhs_bandwidth_multipliers = 1 / sqrt(2),
    rkhs_cv_folds = 3L,
    rkhs_cv_repeats = 2L,
    rkhs_selection_rule = c("minimum", "one_se"),
    rkhs_bandwidth_multiplier = 1 / sqrt(2),
    rkhs_gamma_floor = 1e-6,
    rkhs_gamma_cap = 100,
    rkhs_backend = c("exact", "nystrom"),
    rkhs_nystrom_rank = NULL,
    rkhs_kernel = c("gaussian", "sobolev_matern32", "sobolev_matern52"),
    rkhs_bandwidth_sample = 1000L
) {
  binary_rkhs_gamma_requirements(require_oracle = TRUE)
  rkhs_selection_rule <- match.arg(rkhs_selection_rule)
  rkhs_backend <- match.arg(rkhs_backend)
  rkhs_kernel <- match.arg(rkhs_kernel)
  required <- c("Y", "T_obs", "Z", x_cols)
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  gamma_tuning_mode <- if (is.null(rkhs_lambda_by_t)) "joint_cv" else "fixed"
  if (
    gamma_tuning_mode == "joint_cv" &&
      (length(rkhs_cv_folds) != 1L || !is.finite(rkhs_cv_folds) ||
        rkhs_cv_folds < 2L || rkhs_cv_folds != as.integer(rkhs_cv_folds))
  ) {
    stop("rkhs_cv_folds must be one integer of at least two.")
  }
  if (
    gamma_tuning_mode == "joint_cv" &&
      (length(rkhs_cv_repeats) != 1L || !is.finite(rkhs_cv_repeats) ||
        rkhs_cv_repeats < 1L ||
        rkhs_cv_repeats != as.integer(rkhs_cv_repeats))
  ) {
    stop("rkhs_cv_repeats must be one positive integer.")
  }
  if (is.null(fold_id)) {
    fold_id <- binary_make_tz_folds(data, K = K, seed = seed)
  }
  if (length(fold_id) != nrow(data)) stop("fold_id has the wrong length.")
  folds <- sort(unique(fold_id))
  if (length(folds) < 2L) stop("At least two outer folds are required.")

  prediction_names <- c(
    "e", "pi0", "pi1",
    "a_t0_z0", "a_t0_z1", "a_t1_z0", "a_t1_z1"
  )
  make_prediction_frame <- function() {
    as.data.frame(setNames(
      replicate(
        length(prediction_names), rep(NA_real_, nrow(data)), simplify = FALSE
      ),
      prediction_names
    ))
  }
  variant_names <- c(
    "official_probability_ratio",
    "rkhs_gamma_estimated",
    "rkhs_gamma_oracle_weights",
    "rkhs_gamma_oracle_ep_a0",
    "oracle_gamma_estimated_a0",
    "oracle_all"
  )
  predictions <- setNames(
    replicate(length(variant_names), make_prediction_frame(), simplify = FALSE),
    variant_names
  )
  fit_diagnostics <- gamma_diagnostics <- reconstruction_diagnostics <- list()
  gamma_cv_tables <- gamma_cv_resample_tables <- list()
  fit_counter <- gamma_counter <- reconstruction_counter <- 0L

  clip_design <- function(value) {
    binary_clip_probability(value, overlap_epsilon)
  }
  clip_outcome <- function(value) {
    binary_clip_probability(value, outcome_probability_epsilon)
  }

  for (fold in folds) {
    training <- data[fold_id != fold, , drop = FALSE]
    validation <- data[fold_id == fold, , drop = FALSE]
    validation_index <- which(fold_id == fold)

    fits <- list()
    fits$e <- binary_rkhs_gamma_fit_component(
      training, "Z", x_cols, seed + 1000L * fold + 1L,
      spline_df, gam_min_basis, gam_max_basis, gam_growth_exponent,
      library_candidates, library_cv_folds, library_selection_rule,
      library_ranger_num_trees
    )
    fits$pi0 <- binary_rkhs_gamma_fit_component(
      training[training$Z == 0, , drop = FALSE],
      "T_obs", x_cols, seed + 1000L * fold + 2L,
      spline_df, gam_min_basis, gam_max_basis, gam_growth_exponent,
      library_candidates, library_cv_folds, library_selection_rule,
      library_ranger_num_trees
    )
    fits$pi1 <- binary_rkhs_gamma_fit_component(
      training[training$Z == 1, , drop = FALSE],
      "T_obs", x_cols, seed + 1000L * fold + 3L,
      spline_df, gam_min_basis, gam_max_basis, gam_growth_exponent,
      library_candidates, library_cv_folds, library_selection_rule,
      library_ranger_num_trees
    )
    for (treatment in 0:1) {
      for (z_value in 0:1) {
        key <- sprintf("a_t%d_z%d", treatment, z_value)
        cell <- training$T_obs == treatment & training$Z == z_value
        fits[[key]] <- binary_rkhs_gamma_fit_component(
          training[cell, , drop = FALSE],
          "Y", x_cols,
          seed + 1000L * fold + 100L * treatment + 10L * z_value + 4L,
          spline_df, gam_min_basis, gam_max_basis, gam_growth_exponent,
          library_candidates, library_cv_folds, library_selection_rule,
          library_ranger_num_trees
        )
      }
    }

    estimated_training <- data.frame(
      e = clip_design(1 - fits$e$predict(training)),
      pi0 = clip_design(fits$pi0$predict(training)),
      pi1 = clip_design(fits$pi1$predict(training)),
      stringsAsFactors = FALSE
    )
    estimated_validation <- data.frame(
      e = clip_design(1 - fits$e$predict(validation)),
      pi0 = clip_design(fits$pi0$predict(validation)),
      pi1 = clip_design(fits$pi1$predict(validation)),
      stringsAsFactors = FALSE
    )
    for (key in setdiff(prediction_names, c("e", "pi0", "pi1"))) {
      estimated_training[[key]] <- clip_outcome(fits[[key]]$predict(training))
      estimated_validation[[key]] <- clip_outcome(fits[[key]]$predict(validation))
    }
    oracle_training <- binary_rkhs_gamma_oracle_predictions(
      training, specification, x_cols
    )
    oracle_validation <- binary_rkhs_gamma_oracle_predictions(
      validation, specification, x_cols
    )

    predictions$official_probability_ratio[
      validation_index, prediction_names
    ] <- estimated_validation[prediction_names]
    predictions$oracle_all[validation_index, prediction_names] <-
      oracle_validation[prediction_names]

    rkhs_by_weight <- list()
    for (weight_source in c("estimated", "oracle")) {
      nuisance_training <- if (weight_source == "estimated") {
        estimated_training
      } else {
        oracle_training
      }
      rkhs_by_weight[[weight_source]] <- list()
      for (treatment in 0:1) {
        key <- as.character(treatment)
        n_training_t <- sum(training$T_obs == treatment)
        gamma_seed <- seed + 100000L * fold + 1000L * treatment +
          ifelse(weight_source == "oracle", 500000L, 0L)
        gamma_fit <- if (gamma_tuning_mode == "joint_cv") {
          fit_gamma_rkhs(
            train = training,
            nuisance_train = nuisance_training,
            t_value = treatment,
            x_cols = x_cols,
            lambda_grid = rkhs_lambda_grid,
            bandwidth_multipliers = rkhs_bandwidth_multipliers,
            cv_folds = rkhs_cv_folds,
            cv_repeats = rkhs_cv_repeats,
            selection_rule = rkhs_selection_rule,
            gamma_floor = rkhs_gamma_floor,
            gamma_cap = rkhs_gamma_cap,
            bandwidth_sample = rkhs_bandwidth_sample,
            backend = rkhs_backend,
            nystrom_rank = rkhs_nystrom_rank,
            kernel_type = rkhs_kernel,
            seed = gamma_seed
          )
        } else {
          selected_lambda <- binary_resolve_rkhs_lambda(
            rkhs_lambda_by_t = rkhs_lambda_by_t,
            treatment = treatment,
            n_training = n_training_t
          )
          fit_gamma_rkhs_fixed(
            train = training,
            nuisance_train = nuisance_training,
            t_value = treatment,
            x_cols = x_cols,
            lambda = selected_lambda,
            bandwidth_multiplier = rkhs_bandwidth_multiplier,
            gamma_floor = rkhs_gamma_floor,
            gamma_cap = rkhs_gamma_cap,
            bandwidth_sample = rkhs_bandwidth_sample,
            backend = rkhs_backend,
            nystrom_rank = rkhs_nystrom_rank,
            kernel_type = rkhs_kernel,
            seed = gamma_seed
          )
        }
        if (!is.null(gamma_fit$cv_table) && nrow(gamma_fit$cv_table)) {
          table <- gamma_fit$cv_table
          table$outer_fold <- fold
          table$treatment <- treatment
          table$weight_source <- weight_source
          table$n_training_t <- n_training_t
          table$lambda_sqrt_n_constant <- table$lambda * sqrt(n_training_t)
          table$selection_rule <- rkhs_selection_rule
          gamma_cv_tables[[length(gamma_cv_tables) + 1L]] <- table
        }
        if (
          !is.null(gamma_fit$cv_resample_table) &&
            nrow(gamma_fit$cv_resample_table)
        ) {
          resample_table <- gamma_fit$cv_resample_table
          resample_table$outer_fold <- fold
          resample_table$treatment <- treatment
          resample_table$weight_source <- weight_source
          resample_table$n_training_t <- n_training_t
          resample_table$lambda_sqrt_n_constant <-
            resample_table$lambda * sqrt(n_training_t)
          resample_table$selection_rule <- rkhs_selection_rule
          gamma_cv_resample_tables[[
            length(gamma_cv_resample_tables) + 1L
          ]] <- resample_table
        }
        new_zero <- new_one <- validation
        new_zero$Y <- 0
        new_one$Y <- 1
        prediction_zero <- gamma_fit$predict(new_zero)
        prediction_one <- gamma_fit$predict(new_one)
        true_zero <- observed_binary_gamma(
          0, validation[, x_cols, drop = FALSE], treatment, specification
        )
        true_one <- observed_binary_gamma(
          1, validation[, x_cols, drop = FALSE], treatment, specification
        )
        rkhs_by_weight[[weight_source]][[key]] <- list(
          zero = prediction_zero$normalized,
          one = prediction_one$normalized,
          fit = gamma_fit
        )
        gamma_counter <- gamma_counter + 1L
        gamma_diagnostics[[gamma_counter]] <- data.frame(
          fold = fold,
          treatment = treatment,
          weight_source = weight_source,
          tuning_mode = gamma_tuning_mode,
          cv_folds = if (gamma_tuning_mode == "joint_cv") {
            as.integer(rkhs_cv_folds)
          } else NA_integer_,
          cv_repeats = if (gamma_tuning_mode == "joint_cv") {
            as.integer(rkhs_cv_repeats)
          } else NA_integer_,
          n_training_t = n_training_t,
          selected_lambda = gamma_fit$selected_lambda,
          selected_lambda_sqrt_n_constant =
            gamma_fit$selected_lambda * sqrt(n_training_t),
          selected_bandwidth_multiplier =
            gamma_fit$selected_bandwidth_multiplier,
          base_bandwidth = gamma_fit$base_bandwidth,
          selected_bandwidth = gamma_fit$selected_bandwidth,
          normalization_constant = gamma_fit$normalization_constant,
          training_negative_rate = gamma_fit$training_negative_rate,
          nystrom_rank_requested = if (
            is.null(gamma_fit$nystrom_rank_requested)
          ) NA_integer_ else gamma_fit$nystrom_rank_requested,
          nystrom_landmark_count = if (
            is.null(gamma_fit$nystrom_landmark_count)
          ) NA_integer_ else gamma_fit$nystrom_landmark_count,
          gamma_zero_bias = mean(prediction_zero$normalized - true_zero),
          gamma_zero_rmse = sqrt(mean(
            (prediction_zero$normalized - true_zero)^2
          )),
          gamma_one_bias = mean(prediction_one$normalized - true_one),
          gamma_one_rmse = sqrt(mean(
            (prediction_one$normalized - true_one)^2
          )),
          gamma_zero_correlation = suppressWarnings(stats::cor(
            prediction_zero$normalized, true_zero
          )),
          gamma_one_correlation = suppressWarnings(stats::cor(
            prediction_one$normalized, true_one
          )),
          stringsAsFactors = FALSE
        )
      }
    }

    construct_variant <- function(
        variant,
        endpoint_design,
        baseline_source,
        gamma_source,
        gamma_weight_source = NULL
    ) {
      local <- if (endpoint_design == "estimated") {
        estimated_validation
      } else {
        oracle_validation
      }
      for (treatment in 0:1) {
        baseline_name <- sprintf("a_t%d_z0", treatment)
        comparison_name <- sprintf("a_t%d_z1", treatment)
        a0 <- if (baseline_source == "estimated") {
          estimated_validation[[baseline_name]]
        } else {
          oracle_validation[[baseline_name]]
        }
        if (gamma_source == "oracle") {
          gamma_zero <- observed_binary_gamma(
            0, validation[, x_cols, drop = FALSE], treatment, specification
          )
          gamma_one <- observed_binary_gamma(
            1, validation[, x_cols, drop = FALSE], treatment, specification
          )
        } else {
          gamma_zero <- rkhs_by_weight[[gamma_weight_source]][[
            as.character(treatment)
          ]]$zero
          gamma_one <- rkhs_by_weight[[gamma_weight_source]][[
            as.character(treatment)
          ]]$one
        }
        reconstructed <- binary_rkhs_gamma_reconstruct(
          a0 = a0,
          gamma_zero = gamma_zero,
          gamma_one = gamma_one,
          probability_epsilon = outcome_probability_epsilon,
          gamma_floor = rkhs_gamma_floor
        )
        local[[baseline_name]] <- a0
        local[[comparison_name]] <- reconstructed$a1
        reconstruction_counter <<- reconstruction_counter + 1L
        reconstruction_diagnostics[[reconstruction_counter]] <<- data.frame(
          fold = fold,
          variant = variant,
          treatment = treatment,
          mean_absolute_pre_normalization_error = mean(abs(
            reconstructed$normalizer - 1
          )),
          maximum_absolute_pre_normalization_error = max(abs(
            reconstructed$normalizer - 1
          )),
          reconstructed_a1_bias = mean(
            reconstructed$a1 - oracle_validation[[comparison_name]]
          ),
          reconstructed_a1_rmse = sqrt(mean(
            (reconstructed$a1 - oracle_validation[[comparison_name]])^2
          )),
          lower_probability_clipping_rate =
            reconstructed$lower_probability_clipping_rate,
          upper_probability_clipping_rate =
            reconstructed$upper_probability_clipping_rate,
          stringsAsFactors = FALSE
        )
      }
      predictions[[variant]][validation_index, prediction_names] <<-
        local[prediction_names]
      invisible(NULL)
    }

    construct_variant(
      "rkhs_gamma_estimated", "estimated", "estimated", "rkhs", "estimated"
    )
    construct_variant(
      "rkhs_gamma_oracle_weights", "estimated", "estimated", "rkhs", "oracle"
    )
    construct_variant(
      "rkhs_gamma_oracle_ep_a0", "oracle", "oracle", "rkhs", "oracle"
    )
    construct_variant(
      "oracle_gamma_estimated_a0", "estimated", "estimated", "oracle"
    )

    for (component in names(fits)) {
      fit_counter <- fit_counter + 1L
      fit <- fits[[component]]
      fit_diagnostics[[fit_counter]] <- data.frame(
        fold = fold,
        component = component,
        selected_learner = if (is.null(fit$selected_learner)) {
          fit$learner
        } else {
          fit$selected_learner
        },
        converged = isTRUE(fit$converged),
        cv_mean_risk = if (is.null(fit$cv_mean_risk)) {
          NA_real_
        } else {
          fit$cv_mean_risk
        },
        stringsAsFactors = FALSE
      )
    }
  }

  if (any(vapply(predictions, function(value) {
    any(!is.finite(as.matrix(value)))
  }, logical(1)))) {
    stop("At least one binary RKHS gamma nuisance prediction is non-finite.")
  }
  list(
    predictions = predictions,
    fold_id = fold_id,
    variants = data.frame(
      variant = variant_names,
      description = c(
        "Legacy comparator with separately estimated a_tz probabilities",
        "Official estimated e, pi, a_t0 and direct RKHS gamma",
        "Estimated e, pi, a_t0 and oracle-weight direct RKHS gamma",
        "Oracle e, pi, a_t0 and oracle-weight direct RKHS gamma",
        "Estimated e, pi, a_t0 and oracle gamma",
        "All oracle nuisances"
      ),
      stringsAsFactors = FALSE
    ),
    fit_diagnostics = do.call(rbind, fit_diagnostics),
    gamma_cv_table = if (length(gamma_cv_tables)) {
      do.call(rbind, gamma_cv_tables)
    } else {
      data.frame()
    },
    gamma_cv_resample_table = if (length(gamma_cv_resample_tables)) {
      do.call(rbind, gamma_cv_resample_tables)
    } else {
      data.frame()
    },
    gamma_diagnostics = do.call(rbind, gamma_diagnostics),
    reconstruction_diagnostics = do.call(rbind, reconstruction_diagnostics),
    configuration = list(
      K = length(folds),
      overlap_epsilon = overlap_epsilon,
      outcome_probability_epsilon = outcome_probability_epsilon,
      rkhs_tuning_mode = gamma_tuning_mode,
      rkhs_lambda_by_t = rkhs_lambda_by_t,
      rkhs_lambda_grid_at_n1 = if (gamma_tuning_mode == "joint_cv") {
        resolve_gamma_lambda_grid(
          rkhs_lambda_grid, n_training = 1L, treatment = 0L
        )
      } else {
        NULL
      },
      rkhs_bandwidth_multipliers = if (gamma_tuning_mode == "joint_cv") {
        rkhs_bandwidth_multipliers
      } else {
        rkhs_bandwidth_multiplier
      },
      rkhs_cv_folds = if (gamma_tuning_mode == "joint_cv") {
        as.integer(rkhs_cv_folds)
      } else {
        NA_integer_
      },
      rkhs_cv_repeats = if (gamma_tuning_mode == "joint_cv") {
        as.integer(rkhs_cv_repeats)
      } else {
        NA_integer_
      },
      rkhs_selection_rule = if (gamma_tuning_mode == "joint_cv") {
        rkhs_selection_rule
      } else {
        NA_character_
      },
      rkhs_gamma_floor = rkhs_gamma_floor,
      rkhs_gamma_cap = rkhs_gamma_cap,
      rkhs_backend = rkhs_backend,
      rkhs_nystrom_rank = rkhs_nystrom_rank,
      rkhs_kernel = rkhs_kernel
    )
  )
}
