# Honest learned-witness test for quantitative Dual-MSM compatibility.
#
# Required first:
#   source("R/sharp_bounds_plugin.R")
#   source("R/conservative_compatibility_test.R")
#
# For each cross-fitting fold and candidate Lambda, the training sample learns
# allocation functions r(Y,X) that concentrate on the four model-implied
# nonemptiness/normalization violations. For the two conditional-normalization
# restrictions it
# also learns a nonnegative localization h(X). The witnesses are frozen before
# their one-step moment signals are evaluated on the held-out fold.
#
# This construction has a deliberately limited inferential claim. Conditional
# on suitable nuisance-rate and stability conditions, every learned held-out
# moment remains a necessary implication of the candidate model. It is more
# targeted than a fixed omnibus dictionary, but it is not a uniform adaptive
# test over all measurable witnesses.

check_learned_witness_dependencies <- function() {
  required <- c(
    "fit_nuisance_models", "fit_gamma_rkhs_fixed",
    "make_tz_stratified_folds", "fit_series_conditional_allocations",
    "compute_dual_msm_moment_primitives", "make_arm_lambda_signal_blocks"
  )
  missing <- required[
    !vapply(required, exists, logical(1), mode = "function")
  ]
  if (length(missing) > 0L) {
    stop(
      "Source R/sharp_bounds_plugin.R and ",
      "R/conservative_compatibility_test.R first. Missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

make_tz_cluster_stratified_folds <- function(
    data,
    cluster_id,
    K = 5L,
    seed = 1L
) {
  n <- nrow(data)
  K <- as.integer(K)
  if (length(cluster_id) != n || anyNA(cluster_id)) {
    stop("cluster_id must be complete and have one value per observation.")
  }
  if (!is.finite(K) || K < 2L) stop("K must be at least two.")
  cluster_id <- as.character(cluster_id)
  cluster_rows <- split(seq_len(n), cluster_id)
  if (length(cluster_rows) < K) {
    stop("There must be at least K independent clusters.")
  }
  cluster_summary <- lapply(cluster_rows, function(index) {
    local_z <- unique(data$Z[index])
    if (length(local_z) != 1L) {
      stop("Every cluster must have one common instrument value.")
    }
    c(
      Z = local_z,
      size = length(index),
      treated = sum(data$T_obs[index] == 1)
    )
  })
  cluster_summary <- as.data.frame(do.call(rbind, cluster_summary))
  cluster_summary$cluster <- rownames(cluster_summary)
  # Exact strata preserve instrument status, cluster size, and treatment
  # composition. In the OHIE application these correspond to the randomized
  # household and its observed treatment pattern.
  stratum <- interaction(
    cluster_summary$Z,
    cluster_summary$size,
    cluster_summary$treated,
    drop = TRUE
  )
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  cluster_fold <- integer(nrow(cluster_summary))
  for (level in levels(stratum)) {
    index <- which(stratum == level)
    index <- sample(index, length(index), replace = FALSE)
    cluster_fold[index] <- rep(seq_len(K), length.out = length(index))
  }
  fold_id <- cluster_fold[match(cluster_id, cluster_summary$cluster)]
  if (anyNA(fold_id) || any(fold_id < 1L | fold_id > K)) {
    stop("Cluster-stratified fold construction failed.")
  }
  if (any(vapply(
      split(fold_id, cluster_id),
      function(value) length(unique(value)) != 1L,
      logical(1)
    ))) {
    stop("At least one cluster was split across folds.")
  }
  cell <- table(fold_id, data$T_obs, data$Z)
  if (any(cell == 0L)) {
    stop("A cluster-stratified validation fold has an empty T-by-Z cell.")
  }
  as.integer(fold_id)
}

learned_witness_probability_fit <- function(
    training,
    validation,
    x_cols,
    learner,
    positivity_epsilon,
    treatment_spline_df,
    gam_basis_dimension,
    integrated_confounder_strength,
    integrated_quadrature_nodes,
    super_learner_candidates = NULL,
    super_learner_cv_folds = 3L,
    super_learner_ranger_num_trees = 500L,
    super_learner_hal_max_degree = 2L,
    library_seed = 1L
) {
  if (learner == "oracle") {
    required <- c("e_true", "pi0_true", "pi1_true")
    if (!all(required %in% names(training)) ||
        !all(required %in% names(validation))) {
      stop("Oracle learned-witness testing requires true probability columns.")
    }
    return(list(
      training = list(
        e = training$e_true,
        pi0 = training$pi0_true,
        pi1 = training$pi1_true
      ),
      validation = list(
        e = validation$e_true,
        pi0 = validation$pi0_true,
        pi1 = validation$pi1_true
      ),
      fit = NULL
    ))
  }
  if (learner == "glm_main") {
    formula_for <- function(response) stats::reformulate(x_cols, response)
    instrument_fit <- stats::glm(
      formula_for("Z"), data = training, family = stats::binomial()
    )
    treatment_fit_z0 <- stats::glm(
      formula_for("T_obs"),
      data = training[training$Z == 0, , drop = FALSE],
      family = stats::binomial()
    )
    treatment_fit_z1 <- stats::glm(
      formula_for("T_obs"),
      data = training[training$Z == 1, , drop = FALSE],
      family = stats::binomial()
    )
    predict_probability <- function(newdata) list(
      e = clip_probability(
        1 - stats::predict(instrument_fit, newdata, type = "response"),
        positivity_epsilon
      ),
      pi0 = clip_probability(
        stats::predict(treatment_fit_z0, newdata, type = "response"),
        positivity_epsilon
      ),
      pi1 = clip_probability(
        stats::predict(treatment_fit_z1, newdata, type = "response"),
        positivity_epsilon
      )
    )
    return(list(
      training = predict_probability(training),
      validation = predict_probability(validation),
      fit = list(
        instrument = instrument_fit,
        treatment_z0 = treatment_fit_z0,
        treatment_z1 = treatment_fit_z1
      )
    ))
  }
  if (learner %in% c(
      "binary_glm_spline", "binary_gam", "binary_glmnet",
      "binary_library_cv", "binary_super_learner"
  )) {
    if (!exists("binary_fit_regression_with_learner", mode = "function")) {
      stop("Source R/binary_outcome_bounds.R for binary nuisance learners.")
    }
    binary_learner <- sub("^binary_", "", learner)
    fit_one <- function(local_data, response, seed_offset) {
      binary_fit_regression_with_learner(
        local_data,
        response,
        x_cols,
        if (binary_learner == "gam") gam_basis_dimension else treatment_spline_df,
        learner = binary_learner,
        super_learner_candidates = super_learner_candidates,
        super_learner_cv_folds = super_learner_cv_folds,
        super_learner_ranger_num_trees = super_learner_ranger_num_trees,
        super_learner_hal_max_degree = super_learner_hal_max_degree,
        seed = as.integer(library_seed) + seed_offset
      )
    }
    instrument_fit <- fit_one(training, "Z", 1L)
    treatment_fit_z0 <- fit_one(
      training[training$Z == 0, , drop = FALSE], "T_obs", 2L
    )
    treatment_fit_z1 <- fit_one(
      training[training$Z == 1, , drop = FALSE], "T_obs", 3L
    )
    predict_probability <- function(newdata) list(
      e = clip_probability(
        1 - instrument_fit$predict(newdata), positivity_epsilon
      ),
      pi0 = clip_probability(
        treatment_fit_z0$predict(newdata), positivity_epsilon
      ),
      pi1 = clip_probability(
        treatment_fit_z1$predict(newdata), positivity_epsilon
      )
    )
    return(list(
      training = predict_probability(training),
      validation = predict_probability(validation),
      fit = list(
        instrument = instrument_fit,
        treatment_z0 = treatment_fit_z0,
        treatment_z1 = treatment_fit_z1
      )
    ))
  }
  fit <- fit_nuisance_models(
    train = training,
    x_cols = x_cols,
    positivity_epsilon = positivity_epsilon,
    learner = learner,
    treatment_spline_df = treatment_spline_df,
    gam_basis_dimension = gam_basis_dimension,
    integrated_confounder_strength = integrated_confounder_strength,
    integrated_quadrature_nodes = integrated_quadrature_nodes,
    super_learner_candidates = super_learner_candidates,
    super_learner_cv_folds = super_learner_cv_folds,
    super_learner_ranger_num_trees = super_learner_ranger_num_trees,
    super_learner_hal_max_degree = super_learner_hal_max_degree,
    library_seed = as.integer(library_seed)
  )
  list(
    training = fit$predict(training),
    validation = fit$predict(validation),
    fit = fit
  )
}

fit_learned_witness_binary_gamma <- function(
    training,
    validation,
    treatment,
    x_cols,
    learner = c(
      "glm_spline", "gam", "glmnet", "library_cv", "super_learner"
    ),
    spline_df = 5L,
    probability_epsilon = 0.01,
    super_learner_candidates = NULL,
    super_learner_cv_folds = 3L,
    super_learner_ranger_num_trees = 500L,
    super_learner_hal_max_degree = 2L,
    seed = 1L
) {
  learner <- match.arg(learner)
  if (!exists("binary_fit_regression_with_learner", mode = "function")) {
    stop("Source R/binary_outcome_bounds.R for the binary gamma learner.")
  }
  if (!all(training$Y %in% c(0, 1)) ||
      !all(validation$Y %in% c(0, 1))) {
    stop("The binary gamma learner requires a binary outcome.")
  }
  fit <- lapply(0:1, function(z_value) {
    cell <- training$T_obs == treatment & training$Z == z_value
    binary_fit_regression_with_learner(
      training[cell, , drop = FALSE],
      "Y",
      x_cols,
      spline_df,
      learner = learner,
      super_learner_candidates = super_learner_candidates,
      super_learner_cv_folds = super_learner_cv_folds,
      super_learner_ranger_num_trees = super_learner_ranger_num_trees,
      super_learner_hal_max_degree = super_learner_hal_max_degree,
      seed = seed + z_value
    )
  })
  predict_gamma <- function(newdata) {
    a0 <- clip_probability(fit[[1L]]$predict(newdata), probability_epsilon)
    a1 <- clip_probability(fit[[2L]]$predict(newdata), probability_epsilon)
    ifelse(newdata$Y == 1, a1 / a0, (1 - a1) / (1 - a0))
  }
  list(
    training = predict_gamma(training),
    validation = predict_gamma(validation),
    fit = fit,
    learner = learner
  )
}

fit_learned_witness_discrete_gamma <- function(
    training,
    validation,
    treatment,
    x_cols,
    support = sort(unique(c(training$Y, validation$Y))),
    df_grid = 1:3,
    cv_folds = 3L,
    pmf_epsilon = 1e-5,
    multinomial_decay = 1e-4,
    seed = 1L
) {
  required <- c(
    "discrete_fit_multinomial_spline_cv", "discrete_normalize_pmf"
  )
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing)) {
    stop(
      "Source R/discrete_outcome_bounds.R for the finite-support gamma ",
      "learner. Missing: ", paste(missing, collapse = ", "), "."
    )
  }
  support <- sort(unique(as.numeric(support)))
  if (length(support) < 2L || any(!is.finite(support)) ||
      any(!training$Y %in% support) || any(!validation$Y %in% support)) {
    stop("The finite-support gamma learner received an invalid support.")
  }
  df_grid <- sort(unique(as.integer(df_grid)))
  if (!length(df_grid) || anyNA(df_grid) || any(df_grid < 1L)) {
    stop("df_grid must contain positive integers.")
  }
  if (length(cv_folds) != 1L || !is.finite(cv_folds) || cv_folds < 2L) {
    stop("cv_folds must be one integer of at least two.")
  }
  if (length(pmf_epsilon) != 1L || !is.finite(pmf_epsilon) ||
      pmf_epsilon <= 0 || pmf_epsilon >= 1 / length(support)) {
    stop("pmf_epsilon is invalid for the requested support.")
  }
  fit <- lapply(0:1, function(z_value) {
    cell <- training$T_obs == treatment & training$Z == z_value
    if (sum(cell) < max(50L, 5L * length(support))) {
      stop(
        "Too few observations in a treatment-instrument cell for the ",
        "finite-support gamma learner."
      )
    }
    discrete_fit_multinomial_spline_cv(
      data = training[cell, , drop = FALSE],
      response = "Y",
      x_cols = x_cols,
      support = support,
      weights = rep(1, sum(cell)),
      df_grid = df_grid,
      cv_folds = as.integer(cv_folds),
      decay = multinomial_decay,
      seed = as.integer(seed + z_value)
    )
  })
  predict_gamma <- function(newdata) {
    p0 <- discrete_normalize_pmf(fit[[1L]]$predict(newdata), pmf_epsilon)
    p1 <- discrete_normalize_pmf(fit[[2L]]$predict(newdata), pmf_epsilon)
    category <- match(newdata$Y, support)
    if (anyNA(category)) stop("A prediction outcome falls outside support.")
    index <- cbind(seq_len(nrow(newdata)), category)
    as.numeric(p1[index] / p0[index])
  }
  list(
    training = predict_gamma(training),
    validation = predict_gamma(validation),
    fit = fit,
    learner = "discrete_multinomial_spline",
    selected_df_z0 = fit[[1L]]$selected_df,
    selected_df_z1 = fit[[2L]]$selected_df,
    support = support,
    pmf_epsilon = pmf_epsilon
  )
}

compute_learned_witness_scores <- function(
    gamma,
    pi0,
    pi1,
    treatment,
    Lambda
) {
  p0 <- if (treatment == 1L) pi0 else 1 - pi0
  p1 <- if (treatment == 1L) pi1 else 1 - pi1
  primitive <- compute_dual_msm_moment_primitives(
    p0, p1, Lambda, Lambda
  )
  data.frame(
    pointwise_lower = primitive$gamma_lower - gamma,
    pointwise_upper = gamma - primitive$gamma_upper,
    normalization_lower = primitive$c_lower * gamma + primitive$b - 1 / Lambda,
    normalization_upper = Lambda - primitive$c_upper * gamma - primitive$b,
    stringsAsFactors = FALSE
  )
}

make_witness_cutoffs <- function(
    score,
    active,
    upper_tail_fractions = c(0.05, 0.10, 0.20, 0.40)
) {
  local <- score[active & is.finite(score)]
  if (length(local) < 20L) stop("Too few observations to learn a witness.")
  cutoffs <- c(
    zero = 0,
    stats::setNames(
      as.numeric(stats::quantile(
        local,
        probs = 1 - upper_tail_fractions,
        type = 8,
        names = FALSE
      )),
      paste0("top", format(100 * upper_tail_fractions, trim = TRUE))
    )
  )
  cutoffs[!duplicated(round(cutoffs, 12))]
}

fit_witness_allocation_means <- function(
    training,
    validation,
    response_training,
    treatment,
    x_cols,
    spline_df
) {
  output <- list()
  for (z_value in 0:1) {
    local <- training$T_obs == treatment & training$Z == z_value
    if (sum(local) <= spline_df * length(x_cols) + 3L) {
      stop("Too few observations for a learned-witness allocation regression.")
    }
    output[[paste0("train_z", z_value)]] <-
      fit_series_conditional_allocations(
        training = training[local, , drop = FALSE],
        validation = training,
        response = response_training[local, , drop = FALSE],
        x_cols = x_cols,
        spline_df = spline_df
      )
    output[[paste0("validation_z", z_value)]] <-
      fit_series_conditional_allocations(
        training = training[local, , drop = FALSE],
        validation = validation,
        response = response_training[local, , drop = FALSE],
        x_cols = x_cols,
        spline_df = spline_df
      )
  }
  output
}

compute_learned_family_base <- function(
    family,
    primitive,
    m0,
    m1,
    Lambda
) {
  switch(
    family,
    pointwise_lower = primitive$gamma_lower * m0 - m1,
    pointwise_upper = m1 - primitive$gamma_upper * m0,
    normalization_lower = primitive$c_lower * m1 +
      (primitive$b - 1 / Lambda) * m0 - (1 - 1 / Lambda),
    normalization_upper = -primitive$c_upper * m1 +
      (Lambda - primitive$b) * m0 - (Lambda - 1),
    stop("Unknown learned-witness family: ", family)
  )
}

make_localization_cutoffs <- function(
    base,
    family,
    upper_tail_fractions = c(0.10, 0.20, 0.40)
) {
  if (!grepl("^normalization_", family)) return(c(all = -Inf))
  local <- base[is.finite(base)]
  cutoffs <- c(
    positive = 0,
    stats::setNames(
      as.numeric(stats::quantile(
        local,
        probs = 1 - upper_tail_fractions,
        type = 8,
        names = FALSE
      )),
      paste0("top", format(100 * upper_tail_fractions, trim = TRUE))
    ),
    all = -Inf
  )
  cutoffs[!duplicated(round(cutoffs, 12))]
}

select_learned_witness <- function(
    family,
    response_training,
    response_validation,
    mean_fit,
    primitive_training,
    primitive_validation,
    Lambda,
    localization_tail_fractions
) {
  candidates <- list()
  counter <- 0L
  for (column in seq_len(ncol(response_training))) {
    base_training <- compute_learned_family_base(
      family,
      primitive_training,
      mean_fit$train_z0[, column],
      mean_fit$train_z1[, column],
      Lambda
    )
    base_validation <- compute_learned_family_base(
      family,
      primitive_validation,
      mean_fit$validation_z0[, column],
      mean_fit$validation_z1[, column],
      Lambda
    )
    h_cutoffs <- make_localization_cutoffs(
      base_training,
      family,
      upper_tail_fractions = localization_tail_fractions
    )
    for (h_name in names(h_cutoffs)) {
      h_training <- as.numeric(base_training > h_cutoffs[h_name])
      active_fraction <- mean(h_training)
      if (active_fraction < 0.025) next
      plug_in <- h_training * base_training
      standard_deviation <- stats::sd(plug_in)
      criterion <- if (is.finite(standard_deviation) && standard_deviation > 1e-10) {
        sqrt(length(plug_in)) * mean(plug_in) / standard_deviation
      } else {
        -Inf
      }
      counter <- counter + 1L
      candidates[[counter]] <- list(
        column = column,
        r_name = colnames(response_training)[column],
        h_name = h_name,
        h_cutoff = unname(h_cutoffs[h_name]),
        criterion = criterion,
        h_validation = as.numeric(base_validation > h_cutoffs[h_name]),
        base_training_mean = mean(plug_in),
        base_training_standard_deviation = standard_deviation,
        h_training_fraction = active_fraction
      )
    }
  }
  if (length(candidates) == 0L) stop("No learned witness candidate remained.")
  criterion <- vapply(candidates, `[[`, numeric(1), "criterion")
  selected <- candidates[[which.max(criterion)]]
  selected$R_validation <- response_validation[, selected$column, drop = FALSE]
  selected$m0_validation <-
    mean_fit$validation_z0[, selected$column, drop = FALSE]
  selected$m1_validation <-
    mean_fit$validation_z1[, selected$column, drop = FALSE]
  selected$H_validation <- matrix(
    selected$h_validation,
    ncol = 1L,
    dimnames = list(NULL, paste0("learned_", family))
  )
  selected
}

fit_fold_learned_witness_signals <- function(
    training,
    validation,
    probability_training,
    probability_validation,
    gamma_training,
    gamma_validation,
    treatment,
    Lambda,
    x_cols,
    allocation_spline_df,
    witness_tail_fractions,
    localization_tail_fractions,
    witness_shapes
) {
  score_training <- compute_learned_witness_scores(
    gamma_training,
    probability_training$pi0,
    probability_training$pi1,
    treatment,
    Lambda
  )
  score_validation <- compute_learned_witness_scores(
    gamma_validation,
    probability_validation$pi0,
    probability_validation$pi1,
    treatment,
    Lambda
  )
  families <- names(score_training)
  response_training <- response_validation <- list()
  response_family <- character()
  for (family in families) {
    family_response_count <- length(response_training)
    score_local <- score_training[[family]]
    active_local <- training$T_obs == treatment & is.finite(score_local)
    cutoffs <- make_witness_cutoffs(
      score_local,
      active = active_local,
      upper_tail_fractions = witness_tail_fractions
    )
    if ("indicator" %in% witness_shapes) {
      for (name in names(cutoffs)) {
        key <- paste(family, "indicator", name, sep = "__")
        response_training[[key]] <- as.numeric(
          score_training[[family]] > cutoffs[name]
        )
        response_validation[[key]] <- as.numeric(
          score_validation[[family]] > cutoffs[name]
        )
        response_family[key] <- family
      }
    }

    # Smooth nonnegative witnesses retain the necessary-moment implication
    # while concentrating more weight where the compatible-interval violation
    # is larger. The scale is learned only on the training fold and the held-
    # out response is clipped to [0,1], so the witness remains frozen and
    # bounded when evaluated in validation.
    positive_training <- pmax(score_training[[family]], 0)
    positive_validation <- pmax(score_validation[[family]], 0)
    # Scaling does not change the standardized moment, but clipping does. Use
    # the training maximum so the linear witness is proportional to the full
    # positive excess on the training fold; clip only genuinely new held-out
    # extremes.
    scale <- max(positive_training[active_local], na.rm = TRUE)
    if (is.finite(scale) && scale > 1e-10) {
      for (shape in intersect(witness_shapes, c("linear", "sqrt"))) {
        transform <- if (shape == "linear") identity else sqrt
        key <- paste(family, shape, "positive", sep = "__")
        response_training[[key]] <- pmin(
          transform(positive_training / scale), 1
        )
        response_validation[[key]] <- pmin(
          transform(positive_validation / scale), 1
        )
        response_family[key] <- family
      }
    }
    if (length(response_training) == family_response_count) {
      # A candidate Lambda may make this family nonviolating throughout the
      # training fold. Retain the zero-threshold indicator as a harmless
      # fallback so other informative arm/family moments can still be tested.
      key <- paste(family, "fallback", "zero", sep = "__")
      response_training[[key]] <- as.numeric(
        score_training[[family]] > 0
      )
      response_validation[[key]] <- as.numeric(
        score_validation[[family]] > 0
      )
      response_family[key] <- family
    }
  }
  response_training <- do.call(cbind, response_training)
  response_validation <- do.call(cbind, response_validation)
  storage.mode(response_training) <- "double"
  storage.mode(response_validation) <- "double"
  mean_fit <- fit_witness_allocation_means(
    training,
    validation,
    response_training,
    treatment,
    x_cols,
    allocation_spline_df
  )
  p0_training <- if (treatment == 1L) {
    probability_training$pi0
  } else 1 - probability_training$pi0
  p1_training <- if (treatment == 1L) {
    probability_training$pi1
  } else 1 - probability_training$pi1
  p0_validation <- if (treatment == 1L) {
    probability_validation$pi0
  } else 1 - probability_validation$pi0
  p1_validation <- if (treatment == 1L) {
    probability_validation$pi1
  } else 1 - probability_validation$pi1
  primitive_training <- compute_dual_msm_moment_primitives(
    p0_training, p1_training, Lambda, Lambda
  )
  primitive_validation <- compute_dual_msm_moment_primitives(
    p0_validation, p1_validation, Lambda, Lambda
  )

  output_signals <- output_metadata <- diagnostics <- list()
  for (family in families) {
    columns <- which(response_family == family)
    local_mean_fit <- lapply(mean_fit, function(value) {
      value[, columns, drop = FALSE]
    })
    selected <- select_learned_witness(
      family = family,
      response_training = response_training[, columns, drop = FALSE],
      response_validation = response_validation[, columns, drop = FALSE],
      mean_fit = local_mean_fit,
      primitive_training = primitive_training,
      primitive_validation = primitive_validation,
      Lambda = Lambda,
      localization_tail_fractions = localization_tail_fractions
    )
    block <- make_arm_lambda_signal_blocks(
      data = validation,
      H = selected$H_validation,
      R = selected$R_validation,
      m0 = selected$m0_validation,
      m1 = selected$m1_validation,
      e = probability_validation$e,
      pi0 = probability_validation$pi0,
      pi1 = probability_validation$pi1,
      t_value = treatment,
      Lambda0 = Lambda,
      Lambda1 = Lambda
    )
    keep <- block$metadata$family == family
    if (sum(keep) != 1L) stop("A learned witness did not map to one moment.")
    output_signals[[family]] <- block$signals[, keep, drop = FALSE]
    local_metadata <- block$metadata[keep, , drop = FALSE]
    local_metadata$r <- selected$r_name
    local_metadata$h <- paste0("learned_", selected$h_name)
    output_metadata[[family]] <- local_metadata
    diagnostics[[family]] <- data.frame(
      treatment = treatment,
      Lambda = Lambda,
      family = family,
      selected_r = selected$r_name,
      selected_h = selected$h_name,
      selection_criterion = selected$criterion,
      training_plugin_mean = selected$base_training_mean,
      training_plugin_standard_deviation =
        selected$base_training_standard_deviation,
      training_h_fraction = selected$h_training_fraction,
      validation_r_fraction = mean(selected$R_validation),
      validation_h_fraction = mean(selected$H_validation),
      stringsAsFactors = FALSE
    )
  }
  list(
    signals = do.call(cbind, output_signals),
    metadata = do.call(rbind, output_metadata),
    diagnostics = do.call(rbind, diagnostics)
  )
}

calibrate_learned_witness_maximum <- function(
    signals,
    metadata,
    bootstrap_replications = 999L,
    confidence_level = 0.95,
    seed = 1L,
    cluster_id = NULL
) {
  n <- nrow(signals)
  clustered <- !is.null(cluster_id)
  if (is.null(cluster_id)) cluster_id <- seq_len(n)
  if (length(cluster_id) != n || anyNA(cluster_id)) {
    stop("cluster_id must be complete and have one value per signal row.")
  }
  cluster_id <- as.character(cluster_id)
  G <- length(unique(cluster_id))
  if (G < 2L) stop("At least two independent clusters are required.")
  estimate <- colMeans(signals)
  centered_individual <- sweep(signals, 2L, estimate, "-")
  # For the individual-level mean, define the cluster influence contribution
  # psi_g=(G/n) sum_{i in g}{S_i-E(S)}. Then the estimator is the average of
  # G independent cluster contributions and its standard error is
  # sd(psi_g)/sqrt(G). With one row per cluster this is exactly the original
  # independent-row calculation.
  cluster_signal <- rowsum(
    centered_individual, cluster_id, reorder = FALSE
  ) * (G / n)
  standard_deviation <- apply(cluster_signal, 2L, stats::sd)
  keep <- is.finite(estimate) & is.finite(standard_deviation) &
    standard_deviation > 1e-8
  if (!any(keep)) stop("No nondegenerate learned-witness moments remain.")
  signals <- signals[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  estimate <- estimate[keep]
  standard_deviation <- standard_deviation[keep]
  cluster_signal <- cluster_signal[, keep, drop = FALSE]
  z_value <- sqrt(G) * estimate / standard_deviation
  lambda_key <- interaction(metadata$Lambda0, metadata$Lambda1, drop = TRUE)
  lambda_levels <- levels(lambda_key)
  observed <- vapply(
    lambda_levels,
    function(level) max(z_value[lambda_key == level]),
    numeric(1)
  )
  observed_cvm <- vapply(
    lambda_levels,
    function(level) mean(pmax(z_value[lambda_key == level], 0)^2),
    numeric(1)
  )
  centered <- sweep(cluster_signal, 2L, colMeans(cluster_signal), "-")
  centered <- sweep(centered, 2L, standard_deviation, "/")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  bootstrap <- matrix(
    NA_real_, bootstrap_replications, length(lambda_levels),
    dimnames = list(NULL, lambda_levels)
  )
  bootstrap_cvm <- bootstrap
  for (replication in seq_len(bootstrap_replications)) {
    draw <- as.vector(crossprod(stats::rnorm(G), centered)) / sqrt(G)
    bootstrap[replication, ] <- vapply(
      lambda_levels,
      function(level) max(draw[lambda_key == level]),
      numeric(1)
    )
    bootstrap_cvm[replication, ] <- vapply(
      lambda_levels,
      function(level) mean(pmax(draw[lambda_key == level], 0)^2),
      numeric(1)
    )
  }
  alpha <- 1 - confidence_level
  rows <- lapply(seq_along(lambda_levels), function(index) {
    local_metadata <- metadata[lambda_key == lambda_levels[index], ]
    local_bootstrap <- bootstrap[, index]
    local_cvm_bootstrap <- bootstrap_cvm[, index]
    data.frame(
      Lambda0 = local_metadata$Lambda0[1L],
      Lambda1 = local_metadata$Lambda1[1L],
      statistic = observed[index],
      critical_value = as.numeric(stats::quantile(
        local_bootstrap, confidence_level, type = 8, names = FALSE
      )),
      p_value = (1 + sum(local_bootstrap >= observed[index])) /
        (bootstrap_replications + 1),
      reject = FALSE,
      cvm_statistic = observed_cvm[index],
      cvm_critical_value = as.numeric(stats::quantile(
        local_cvm_bootstrap, confidence_level, type = 8, names = FALSE
      )),
      cvm_p_value =
        (1 + sum(local_cvm_bootstrap >= observed_cvm[index])) /
          (bootstrap_replications + 1),
      cvm_reject = FALSE,
      moments_tested = nrow(local_metadata),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  summary$reject <- summary$p_value <= alpha
  summary$cvm_reject <- summary$cvm_p_value <= alpha
  summary <- summary[order(summary$Lambda0, summary$Lambda1), ]

  # Fixed-sequence inversion is valid for the nested nulls M(Lambda): begin at
  # the smallest Lambda and stop at the first nonrejection. The first true null
  # is therefore tested at level alpha, without a blanket grid correction.
  summary$ordered_reject <- FALSE
  summary$cvm_ordered_reject <- FALSE
  continue <- TRUE
  for (index in seq_len(nrow(summary))) {
    if (continue && summary$reject[index]) {
      summary$ordered_reject[index] <- TRUE
    } else {
      continue <- FALSE
    }
  }
  continue <- TRUE
  for (index in seq_len(nrow(summary))) {
    if (continue && summary$cvm_reject[index]) {
      summary$cvm_ordered_reject[index] <- TRUE
    } else {
      continue <- FALSE
    }
  }
  simultaneous_bootstrap <- apply(bootstrap, 1L, max)
  simultaneous_cvm_bootstrap <- apply(bootstrap_cvm, 1L, max)
  summary$grid_simultaneous_p_value <- vapply(
    summary$statistic,
    function(value) (1 + sum(simultaneous_bootstrap >= value)) /
      (bootstrap_replications + 1),
    numeric(1)
  )
  summary$grid_simultaneous_reject <-
    summary$grid_simultaneous_p_value <= alpha
  summary$cvm_grid_simultaneous_p_value <- vapply(
    summary$cvm_statistic,
    function(value) {
      (1 + sum(simultaneous_cvm_bootstrap >= value)) /
        (bootstrap_replications + 1)
    },
    numeric(1)
  )
  summary$cvm_grid_simultaneous_reject <-
    summary$cvm_grid_simultaneous_p_value <= alpha
  moment_detail <- cbind(
    metadata,
    estimate = estimate,
    standard_error = standard_deviation / sqrt(G),
    z_value = z_value
  )
  moment_detail <- moment_detail[
    order(moment_detail$Lambda0, -moment_detail$z_value),
  ]
  list(
    summary = summary,
    moment_detail = moment_detail,
    bootstrap = bootstrap,
    cluster_count = G,
    clustered_inference = clustered
  )
}

test_learned_witness_dual_msm_compatibility <- function(
    data,
    lambda_grid,
    x_cols = "X",
    folds = 5L,
    primitive_learner = c(
      "integrated_parametric", "glm_spline", "glm_main", "gam",
      "ranger", "library_cv", "super_learner",
      "binary_glm_spline", "binary_gam", "binary_glmnet",
      "binary_library_cv", "binary_super_learner", "oracle"
    ),
    gamma_learner = c(
      "rkhs", "rkhs_cv", "binary_glm_spline", "binary_gam",
      "binary_glmnet", "binary_library_cv", "binary_super_learner",
      "discrete_multinomial_spline", "oracle"
    ),
    positivity_epsilon = 0.05,
    treatment_spline_df = 5L,
    gam_basis_dimension = 6L,
    allocation_spline_df = 5L,
    integrated_confounder_strength = 1,
    integrated_quadrature_nodes = 40L,
    super_learner_candidates = NULL,
    super_learner_cv_folds = 3L,
    super_learner_ranger_num_trees = 500L,
    super_learner_hal_max_degree = 2L,
    gamma_lambda_by_t = c("0" = 0.03162278, "1" = 0.01),
    gamma_bandwidth_multiplier = 4,
    gamma_cv_lambda_grid = function(n_training, treatment) {
      10^seq(-1, 0.5, by = 0.5) / sqrt(n_training)
    },
    gamma_cv_bandwidth_multipliers = c(1, 2, 4, 8),
    gamma_cv_folds = 3L,
    gamma_cv_repeats = 2L,
    gamma_cv_selection_rule = c("one_se", "minimum"),
    gamma_nystrom_rank = 80L,
    gamma_kernel = "gaussian",
    gamma_discrete_support = NULL,
    gamma_discrete_df_grid = 1:3,
    gamma_discrete_cv_folds = 3L,
    gamma_discrete_pmf_epsilon = 1e-5,
    gamma_discrete_decay = 1e-4,
    witness_tail_fractions = c(0.05, 0.10, 0.20, 0.40),
    localization_tail_fractions = c(0.10, 0.20, 0.40),
    witness_shapes = c("indicator", "linear", "sqrt"),
    bootstrap_replications = 999L,
    confidence_level = 0.95,
    seed = 1L,
    cluster_id = NULL,
    progress = FALSE
) {
  check_learned_witness_dependencies()
  validate_conservative_test_data(data, x_cols)
  primitive_learner <- match.arg(primitive_learner)
  gamma_learner <- match.arg(gamma_learner)
  gamma_cv_selection_rule <- match.arg(gamma_cv_selection_rule)
  witness_shapes <- unique(as.character(witness_shapes))
  if (!length(witness_shapes) ||
      !all(witness_shapes %in% c("indicator", "linear", "sqrt"))) {
    stop("witness_shapes must use indicator, linear, and/or sqrt.")
  }
  lambda_grid <- sort(unique(as.numeric(lambda_grid)))
  if (!length(lambda_grid) || any(!is.finite(lambda_grid)) ||
      any(lambda_grid < 1)) {
    stop("lambda_grid must contain finite values of at least one.")
  }
  if (gamma_learner == "oracle") {
    required <- c("gamma0_oracle", "gamma1_oracle")
    if (!all(required %in% names(data))) {
      stop("Oracle gamma witnesses require gamma0_oracle and gamma1_oracle.")
    }
  }
  clustered_inference <- !is.null(cluster_id)
  if (is.null(cluster_id)) {
    cluster_id <- seq_len(nrow(data))
  } else if (length(cluster_id) != nrow(data) || anyNA(cluster_id)) {
    stop("cluster_id must be complete and have one value per data row.")
  }
  fold_id <- if (clustered_inference) {
    make_tz_cluster_stratified_folds(
      data, cluster_id = cluster_id, K = folds, seed = seed
    )
  } else {
    make_tz_stratified_folds(data, K = folds, seed = seed)
  }
  signal_rows <- vector("list", folds)
  selected_signal_rows <- vector("list", folds)
  aggregate_signal_rows <- vector("list", folds)
  component_selected_signal_rows <- vector("list", folds)
  component_aggregate_signal_rows <- vector("list", folds)
  metadata_rows <- vector("list", folds)
  diagnostic_rows <- vector("list", folds)
  gamma_diagnostics <- list()
  gamma_cv_diagnostics <- list()
  gamma_cv_resample_diagnostics <- list()
  probability_diagnostics <- list()
  gamma_counter <- 0L
  gamma_cv_counter <- 0L
  for (fold in seq_len(folds)) {
    if (isTRUE(progress)) {
      cat("Starting learned-witness outer fold ", fold, "/", folds, ".\n", sep = "")
    }
    training_index <- which(fold_id != fold)
    validation_index <- which(fold_id == fold)
    training <- data[training_index, , drop = FALSE]
    validation <- data[validation_index, , drop = FALSE]
    probability <- learned_witness_probability_fit(
      training,
      validation,
      x_cols,
      primitive_learner,
      positivity_epsilon,
      treatment_spline_df,
      gam_basis_dimension,
      integrated_confounder_strength,
      integrated_quadrature_nodes,
      super_learner_candidates = super_learner_candidates,
      super_learner_cv_folds = super_learner_cv_folds,
      super_learner_ranger_num_trees =
        super_learner_ranger_num_trees,
      super_learner_hal_max_degree = super_learner_hal_max_degree,
      library_seed = seed + 10000L * fold
    )
    if (!is.null(probability$fit) &&
        is.data.frame(probability$fit$diagnostics)) {
      local_probability_diagnostics <- probability$fit$diagnostics
      local_probability_diagnostics$outer_fold <- fold
      probability_diagnostics[[length(probability_diagnostics) + 1L]] <-
        local_probability_diagnostics
    }
    fold_signals <- fold_metadata <- fold_diagnostics <- list()
    counter <- 0L
    for (treatment in 0:1) {
      n_treatment_training <- sum(training$T_obs == treatment)
      if (gamma_learner == "oracle") {
        gamma_training <- training[[paste0("gamma", treatment, "_oracle")]]
        gamma_validation <- validation[[paste0("gamma", treatment, "_oracle")]]
        gamma_fit <- NULL
      } else if (gamma_learner == "discrete_multinomial_spline") {
        gamma_fit <- fit_learned_witness_discrete_gamma(
          training = training,
          validation = validation,
          treatment = treatment,
          x_cols = x_cols,
          support = if (is.null(gamma_discrete_support)) {
            sort(unique(data$Y))
          } else {
            gamma_discrete_support
          },
          df_grid = gamma_discrete_df_grid,
          cv_folds = gamma_discrete_cv_folds,
          pmf_epsilon = gamma_discrete_pmf_epsilon,
          multinomial_decay = gamma_discrete_decay,
          seed = seed + 10000L * fold + 1000L * treatment
        )
        gamma_training <- gamma_fit$training
        gamma_validation <- gamma_fit$validation
      } else if (grepl("^binary_", gamma_learner)) {
        gamma_fit <- fit_learned_witness_binary_gamma(
          training,
          validation,
          treatment,
          x_cols,
          learner = sub("^binary_", "", gamma_learner),
          spline_df = if (gamma_learner == "binary_gam") {
            gam_basis_dimension
          } else treatment_spline_df,
          super_learner_candidates = super_learner_candidates,
          super_learner_cv_folds = super_learner_cv_folds,
          super_learner_ranger_num_trees =
            super_learner_ranger_num_trees,
          super_learner_hal_max_degree = super_learner_hal_max_degree,
          seed = seed + 10000L * fold + 1000L * treatment
        )
        gamma_training <- gamma_fit$training
        gamma_validation <- gamma_fit$validation
      } else {
        resolved_nystrom_rank <- if (is.function(gamma_nystrom_rank)) {
          gamma_nystrom_rank(
            n_training = n_treatment_training,
            treatment = treatment
          )
        } else if (is.null(gamma_nystrom_rank)) {
          NULL
        } else {
          gamma_nystrom_rank
        }
        if (!is.null(resolved_nystrom_rank)) {
          if (length(resolved_nystrom_rank) != 1L ||
              !is.finite(resolved_nystrom_rank) ||
              resolved_nystrom_rank < 1) {
            stop("The resolved gamma Nystrom rank must be positive.")
          }
          resolved_nystrom_rank <- min(
            as.integer(ceiling(resolved_nystrom_rank)),
            n_treatment_training
          )
        }
        if (gamma_learner == "rkhs_cv") {
          gamma_fit <- fit_gamma_rkhs(
            train = training,
            nuisance_train = probability$training,
            t_value = treatment,
            x_cols = x_cols,
            lambda_grid = gamma_cv_lambda_grid,
            bandwidth_multipliers = gamma_cv_bandwidth_multipliers,
            cv_folds = gamma_cv_folds,
            cv_repeats = gamma_cv_repeats,
            selection_rule = gamma_cv_selection_rule,
            backend = "nystrom",
            nystrom_rank = resolved_nystrom_rank,
            kernel_type = gamma_kernel,
            seed = seed + 10000L * fold + 1000L * treatment
          )
          gamma_cv_counter <- gamma_cv_counter + 1L
          gamma_cv_diagnostics[[gamma_cv_counter]] <- transform(
            gamma_fit$cv_table,
            outer_fold = fold,
            treatment = treatment,
            selection_rule = gamma_cv_selection_rule
          )
          gamma_cv_resample_diagnostics[[gamma_cv_counter]] <- transform(
            gamma_fit$cv_resample_table,
            outer_fold = fold,
            treatment = treatment,
            selection_rule = gamma_cv_selection_rule
          )
        } else {
          gamma_lambda <- if (is.function(gamma_lambda_by_t)) {
            gamma_lambda_by_t(
              n_training = n_treatment_training,
              treatment = treatment
            )
          } else {
            gamma_lambda_by_t[as.character(treatment)]
          }
          if (length(gamma_lambda) != 1L || !is.finite(gamma_lambda) ||
              gamma_lambda <= 0) {
            stop(
              "The resolved fixed gamma RKHS penalty must be one positive number."
            )
          }
          gamma_fit <- fit_gamma_rkhs_fixed(
            train = training,
            nuisance_train = probability$training,
            t_value = treatment,
            x_cols = x_cols,
            lambda = gamma_lambda,
            bandwidth_multiplier = gamma_bandwidth_multiplier,
            backend = "nystrom",
            nystrom_rank = resolved_nystrom_rank,
            kernel_type = gamma_kernel,
            seed = seed + 10000L * fold + 1000L * treatment
          )
        }
        gamma_training <- gamma_fit$predict(training)$normalized
        gamma_validation <- gamma_fit$predict(validation)$normalized
      }
      gamma_counter <- gamma_counter + 1L
      gamma_diagnostics[[gamma_counter]] <- data.frame(
        fold = fold,
        treatment = treatment,
        learner = gamma_learner,
        lambda = if (is.null(gamma_fit) ||
            is.null(gamma_fit$selected_lambda)) {
          NA_real_
        } else gamma_fit$selected_lambda,
        lambda_constant = if (is.null(gamma_fit) ||
            is.null(gamma_fit$selected_lambda)) {
          NA_real_
        } else gamma_fit$selected_lambda * sqrt(n_treatment_training),
        bandwidth_multiplier = if (is.null(gamma_fit) ||
            is.null(gamma_fit$selected_bandwidth_multiplier)) {
          NA_real_
        } else gamma_fit$selected_bandwidth_multiplier,
        selected_bandwidth = if (is.null(gamma_fit) ||
            is.null(gamma_fit$selected_bandwidth)) {
          NA_real_
        } else gamma_fit$selected_bandwidth,
        nystrom_rank_requested = if (is.null(gamma_fit) ||
            is.null(gamma_fit$nystrom_rank_requested)) {
          NA_integer_
        } else gamma_fit$nystrom_rank_requested,
        nystrom_effective_rank = if (is.null(gamma_fit) ||
            is.null(gamma_fit$nystrom_effective_rank)) {
          NA_integer_
        } else gamma_fit$nystrom_effective_rank,
        selected_df_z0 = if (is.null(gamma_fit) ||
            is.null(gamma_fit$selected_df_z0)) {
          NA_integer_
        } else gamma_fit$selected_df_z0,
        selected_df_z1 = if (is.null(gamma_fit) ||
            is.null(gamma_fit$selected_df_z1)) {
          NA_integer_
        } else gamma_fit$selected_df_z1,
        discrete_pmf_epsilon = if (is.null(gamma_fit) ||
            is.null(gamma_fit$pmf_epsilon)) {
          NA_real_
        } else gamma_fit$pmf_epsilon,
        training_mean = mean(gamma_training),
        validation_mean = mean(gamma_validation),
        training_oracle_bias = if (paste0(
            "gamma", treatment, "_oracle"
          ) %in% names(training)) {
          active <- training$T_obs == treatment
          mean(gamma_training[active] - training[[paste0(
            "gamma", treatment, "_oracle"
          )]][active])
        } else NA_real_,
        training_oracle_rmse = if (paste0(
            "gamma", treatment, "_oracle"
          ) %in% names(training)) {
          active <- training$T_obs == treatment
          sqrt(mean((gamma_training[active] - training[[paste0(
            "gamma", treatment, "_oracle"
          )]][active])^2))
        } else NA_real_,
        validation_oracle_bias = if (paste0(
            "gamma", treatment, "_oracle"
          ) %in% names(validation)) {
          active <- validation$T_obs == treatment
          mean(gamma_validation[active] - validation[[paste0(
            "gamma", treatment, "_oracle"
          )]][active])
        } else NA_real_,
        validation_oracle_rmse = if (paste0(
            "gamma", treatment, "_oracle"
          ) %in% names(validation)) {
          active <- validation$T_obs == treatment
          sqrt(mean((gamma_validation[active] - validation[[paste0(
            "gamma", treatment, "_oracle"
          )]][active])^2))
        } else NA_real_,
        stringsAsFactors = FALSE
      )
      if (isTRUE(progress)) {
        cat(
          "Completed gamma nuisance for fold ", fold,
          ", treatment ", treatment, ".\n", sep = ""
        )
      }
      for (Lambda in lambda_grid) {
        counter <- counter + 1L
        local <- fit_fold_learned_witness_signals(
          training,
          validation,
          probability$training,
          probability$validation,
          gamma_training,
          gamma_validation,
          treatment,
          Lambda,
          x_cols,
          allocation_spline_df,
          witness_tail_fractions,
          localization_tail_fractions,
          witness_shapes
        )
        fold_signals[[counter]] <- local$signals
        fold_metadata[[counter]] <- local$metadata
        local$diagnostics$fold <- fold
        fold_diagnostics[[counter]] <- local$diagnostics
      }
    }
    local_signals <- do.call(cbind, fold_signals)
    local_metadata <- do.call(rbind, fold_metadata)
    moment_key <- paste(
      local_metadata$Lambda0,
      local_metadata$treatment,
      local_metadata$family,
      sep = "::"
    )
    ordering <- order(moment_key)
    local_signals <- local_signals[, ordering, drop = FALSE]
    local_metadata <- local_metadata[ordering, , drop = FALSE]
    local_diagnostics <- do.call(rbind, fold_diagnostics)
    if (fold == 1L) {
      reference_key <- moment_key[ordering]
      metadata_reference <- local_metadata
      metadata_reference$r <- "cross_fitted_learned"
      metadata_reference$h <- "cross_fitted_learned"
    } else if (!identical(moment_key[ordering], reference_key)) {
      stop("Learned-witness moment columns differ across folds.")
    }
    # A training-directed nonnegative combination accumulates several modest
    # violations. Every component is a necessary inequality and every weight
    # is nonnegative and learned outside validation, so the combined held-out
    # signal is also a necessary inequality. No validation outcome is used to
    # choose or rescale these weights.
    diagnostic_key <- paste(
      local_diagnostics$Lambda,
      local_diagnostics$treatment,
      local_diagnostics$family,
      sep = "::"
    )
    metadata_key <- paste(
      local_metadata$Lambda0,
      local_metadata$treatment,
      local_metadata$family,
      sep = "::"
    )
    criterion <- local_diagnostics$selection_criterion[
      match(metadata_key, diagnostic_key)
    ]
    training_standard_deviation <-
      local_diagnostics$training_plugin_standard_deviation[
        match(metadata_key, diagnostic_key)
      ]
    if (anyNA(criterion)) stop("Could not match learned aggregation weights.")
    if (anyNA(training_standard_deviation)) {
      stop("Could not match learned aggregation standard deviations.")
    }

    # A more focused version lets the training sample choose both the witness
    # and one of the eight arm-by-family inequalities.  Only that frozen
    # choice is evaluated in validation.  Under the candidate model every
    # possible choice has nonpositive expectation, so honest selection does
    # not require an eight-moment multiplicity correction.  The selected arm
    # and family may differ across folds; the target is their fold-average
    # necessary inequality.
    selected_columns <- vapply(lambda_grid, function(Lambda) {
      index <- which(local_metadata$Lambda0 == Lambda)
      index[which.max(criterion[index])]
    }, integer(1))
    selected_matrix <- local_signals[, selected_columns, drop = FALSE]
    colnames(selected_matrix) <- paste0("Lambda_", lambda_grid)
    selected_key <- metadata_key[selected_columns]
    local_diagnostics$selected_across_families <-
      diagnostic_key %in% selected_key
    aggregate_columns <- lapply(lambda_grid, function(Lambda) {
      index <- which(local_metadata$Lambda0 == Lambda)
      # For signals on different scales, weights proportional only to a
      # training z-score are not scale invariant. The nonnegative
      # estimate-to-variance rule is proportional to z/sd and maximizes the
      # training signal-to-noise ratio under a diagonal covariance
      # approximation. All weights remain training-measurable and
      # nonnegative, preserving the necessary-inequality implication.
      weight <- pmax(criterion[index], 0) /
        pmax(training_standard_deviation[index], 1e-10)
      if (!any(weight > 0)) weight <- rep(1, length(index))
      weight <- weight / sum(weight)
      as.vector(local_signals[, index, drop = FALSE] %*% weight)
    })
    aggregate_matrix <- do.call(cbind, aggregate_columns)
    colnames(aggregate_matrix) <- paste0("Lambda_", lambda_grid)

    # Keep pointwise nonemptiness and conditional normalization in separate
    # components. Under compatibility every arm-by-family moment has
    # nonpositive expectation, so a nonnegative training-directed combination
    # within either component is itself a necessary inequality. Calibrating
    # the maximum of the two held-out component signals tests the single global
    # compatibility null while preventing slack normalization moments from
    # cancelling a pointwise violation (or conversely).
    component_names <- c("pointwise", "normalization")
    component_selected_columns <- list()
    component_aggregate_columns <- list()
    component_counter <- 0L
    for (Lambda in lambda_grid) {
      for (component in component_names) {
        component_counter <- component_counter + 1L
        index <- which(
          local_metadata$Lambda0 == Lambda &
            startsWith(local_metadata$family, paste0(component, "_"))
        )
        if (!length(index)) {
          stop("No learned-witness moments remained for component ", component)
        }

        selected_index <- index[which.max(criterion[index])]
        component_selected_columns[[component_counter]] <-
          local_signals[, selected_index]

        weight <- pmax(criterion[index], 0) /
          pmax(training_standard_deviation[index], 1e-10)
        if (!any(weight > 0)) {
          # A deterministic training-selected fallback avoids averaging several
          # apparently slack inequalities. The choice remains frozen before
          # held-out evaluation and therefore preserves validity.
          weight <- numeric(length(index))
          weight[which.max(criterion[index])] <- 1
        } else {
          weight <- weight / sum(weight)
        }
        component_aggregate_columns[[component_counter]] <- as.vector(
          local_signals[, index, drop = FALSE] %*% weight
        )
      }
    }
    component_column_names <- unlist(lapply(lambda_grid, function(Lambda) {
      paste0("Lambda_", Lambda, "__", component_names)
    }))
    component_selected_matrix <- do.call(cbind, component_selected_columns)
    component_aggregate_matrix <- do.call(cbind, component_aggregate_columns)
    colnames(component_selected_matrix) <- component_column_names
    colnames(component_aggregate_matrix) <- component_column_names
    signal_rows[[fold]] <- data.frame(
      row_id = validation_index,
      local_signals,
      check.names = FALSE
    )
    selected_signal_rows[[fold]] <- data.frame(
      row_id = validation_index,
      selected_matrix,
      check.names = FALSE
    )
    aggregate_signal_rows[[fold]] <- data.frame(
      row_id = validation_index,
      aggregate_matrix,
      check.names = FALSE
    )
    component_selected_signal_rows[[fold]] <- data.frame(
      row_id = validation_index,
      component_selected_matrix,
      check.names = FALSE
    )
    component_aggregate_signal_rows[[fold]] <- data.frame(
      row_id = validation_index,
      component_aggregate_matrix,
      check.names = FALSE
    )
    diagnostic_rows[[fold]] <- local_diagnostics
    if (isTRUE(progress)) {
      cat("Completed learned-witness outer fold ", fold, "/", folds, ".\n", sep = "")
    }
  }
  stacked <- do.call(rbind, signal_rows)
  stacked <- stacked[order(stacked$row_id), , drop = FALSE]
  signals <- as.matrix(stacked[, -1L, drop = FALSE])
  selected_stacked <- do.call(rbind, selected_signal_rows)
  selected_stacked <- selected_stacked[
    order(selected_stacked$row_id), , drop = FALSE
  ]
  selected_signals <- as.matrix(
    selected_stacked[, -1L, drop = FALSE]
  )
  aggregate_stacked <- do.call(rbind, aggregate_signal_rows)
  aggregate_stacked <- aggregate_stacked[
    order(aggregate_stacked$row_id), , drop = FALSE
  ]
  aggregate_signals <- as.matrix(
    aggregate_stacked[, -1L, drop = FALSE]
  )
  component_selected_stacked <- do.call(
    rbind, component_selected_signal_rows
  )
  component_selected_stacked <- component_selected_stacked[
    order(component_selected_stacked$row_id), , drop = FALSE
  ]
  component_selected_signals <- as.matrix(
    component_selected_stacked[, -1L, drop = FALSE]
  )
  component_aggregate_stacked <- do.call(
    rbind, component_aggregate_signal_rows
  )
  component_aggregate_stacked <- component_aggregate_stacked[
    order(component_aggregate_stacked$row_id), , drop = FALSE
  ]
  component_aggregate_signals <- as.matrix(
    component_aggregate_stacked[, -1L, drop = FALSE]
  )
  calibration <- calibrate_learned_witness_maximum(
    signals,
    metadata_reference,
    bootstrap_replications,
    confidence_level,
    seed + 900000L,
    cluster_id = cluster_id
  )
  selected_metadata <- data.frame(
    h = "cross_fitted_learned",
    r = "cross_fitted_learned",
    family = "training_selected_family",
    treatment = -1L,
    Lambda0 = lambda_grid,
    Lambda1 = lambda_grid,
    stringsAsFactors = FALSE
  )
  selected_calibration <- calibrate_learned_witness_maximum(
    selected_signals,
    selected_metadata,
    bootstrap_replications,
    confidence_level,
    seed + 905000L,
    cluster_id = cluster_id
  )
  selected_summary <- selected_calibration$summary
  names(selected_summary)[names(selected_summary) == "statistic"] <-
    "selected_statistic"
  names(selected_summary)[names(selected_summary) == "critical_value"] <-
    "selected_critical_value"
  names(selected_summary)[names(selected_summary) == "p_value"] <-
    "selected_p_value"
  names(selected_summary)[names(selected_summary) == "reject"] <-
    "selected_reject"
  names(selected_summary)[names(selected_summary) == "ordered_reject"] <-
    "selected_ordered_reject"
  names(selected_summary)[
    names(selected_summary) == "grid_simultaneous_p_value"
  ] <- "selected_grid_simultaneous_p_value"
  names(selected_summary)[
    names(selected_summary) == "grid_simultaneous_reject"
  ] <- "selected_grid_simultaneous_reject"
  selected_summary$moments_tested <- NULL
  selected_summary <- selected_summary[
    , !startsWith(names(selected_summary), "cvm_"), drop = FALSE
  ]
  aggregate_metadata <- data.frame(
    h = "cross_fitted_learned",
    r = "cross_fitted_learned",
    family = "directed_aggregate",
    treatment = -1L,
    Lambda0 = lambda_grid,
    Lambda1 = lambda_grid,
    stringsAsFactors = FALSE
  )
  aggregate_calibration <- calibrate_learned_witness_maximum(
    aggregate_signals,
    aggregate_metadata,
    bootstrap_replications,
    confidence_level,
    seed + 910000L,
    cluster_id = cluster_id
  )
  aggregate_summary <- aggregate_calibration$summary
  names(aggregate_summary)[names(aggregate_summary) == "statistic"] <-
    "aggregate_statistic"
  names(aggregate_summary)[names(aggregate_summary) == "critical_value"] <-
    "aggregate_critical_value"
  names(aggregate_summary)[names(aggregate_summary) == "p_value"] <-
    "aggregate_p_value"
  names(aggregate_summary)[names(aggregate_summary) == "reject"] <-
    "aggregate_reject"
  names(aggregate_summary)[names(aggregate_summary) == "ordered_reject"] <-
    "aggregate_ordered_reject"
  names(aggregate_summary)[
    names(aggregate_summary) == "grid_simultaneous_p_value"
  ] <- "aggregate_grid_simultaneous_p_value"
  names(aggregate_summary)[
    names(aggregate_summary) == "grid_simultaneous_reject"
  ] <- "aggregate_grid_simultaneous_reject"
  aggregate_summary$moments_tested <- NULL
  aggregate_summary <- aggregate_summary[
    , !startsWith(names(aggregate_summary), "cvm_"), drop = FALSE
  ]
  component_metadata <- do.call(rbind, lapply(lambda_grid, function(Lambda) {
    data.frame(
      h = "cross_fitted_learned",
      r = "cross_fitted_learned",
      family = c("pointwise_component", "normalization_component"),
      treatment = -1L,
      Lambda0 = Lambda,
      Lambda1 = Lambda,
      stringsAsFactors = FALSE
    )
  }))
  component_selected_calibration <- calibrate_learned_witness_maximum(
    component_selected_signals,
    component_metadata,
    bootstrap_replications,
    confidence_level,
    seed + 915000L,
    cluster_id = cluster_id
  )
  component_selected_summary <- component_selected_calibration$summary
  names(component_selected_summary)[
    names(component_selected_summary) == "statistic"
  ] <- "component_selected_statistic"
  names(component_selected_summary)[
    names(component_selected_summary) == "critical_value"
  ] <- "component_selected_critical_value"
  names(component_selected_summary)[
    names(component_selected_summary) == "p_value"
  ] <- "component_selected_p_value"
  names(component_selected_summary)[
    names(component_selected_summary) == "reject"
  ] <- "component_selected_reject"
  names(component_selected_summary)[
    names(component_selected_summary) == "ordered_reject"
  ] <- "component_selected_ordered_reject"
  names(component_selected_summary)[
    names(component_selected_summary) == "grid_simultaneous_p_value"
  ] <- "component_selected_grid_simultaneous_p_value"
  names(component_selected_summary)[
    names(component_selected_summary) == "grid_simultaneous_reject"
  ] <- "component_selected_grid_simultaneous_reject"
  component_selected_summary$moments_tested <- NULL
  component_selected_summary <- component_selected_summary[
    , !startsWith(names(component_selected_summary), "cvm_"), drop = FALSE
  ]

  component_aggregate_calibration <- calibrate_learned_witness_maximum(
    component_aggregate_signals,
    component_metadata,
    bootstrap_replications,
    confidence_level,
    seed + 920000L,
    cluster_id = cluster_id
  )
  component_aggregate_summary <- component_aggregate_calibration$summary
  names(component_aggregate_summary)[
    names(component_aggregate_summary) == "statistic"
  ] <- "component_aggregate_statistic"
  names(component_aggregate_summary)[
    names(component_aggregate_summary) == "critical_value"
  ] <- "component_aggregate_critical_value"
  names(component_aggregate_summary)[
    names(component_aggregate_summary) == "p_value"
  ] <- "component_aggregate_p_value"
  names(component_aggregate_summary)[
    names(component_aggregate_summary) == "reject"
  ] <- "component_aggregate_reject"
  names(component_aggregate_summary)[
    names(component_aggregate_summary) == "ordered_reject"
  ] <- "component_aggregate_ordered_reject"
  names(component_aggregate_summary)[
    names(component_aggregate_summary) == "grid_simultaneous_p_value"
  ] <- "component_aggregate_grid_simultaneous_p_value"
  names(component_aggregate_summary)[
    names(component_aggregate_summary) == "grid_simultaneous_reject"
  ] <- "component_aggregate_grid_simultaneous_reject"
  component_aggregate_summary$moments_tested <- NULL
  component_aggregate_summary <- component_aggregate_summary[
    , !startsWith(names(component_aggregate_summary), "cvm_"), drop = FALSE
  ]
  calibration$summary <- merge(
    calibration$summary,
    selected_summary,
    by = c("Lambda0", "Lambda1"),
    sort = TRUE
  )
  calibration$summary <- merge(
    calibration$summary,
    aggregate_summary,
    by = c("Lambda0", "Lambda1"),
    sort = TRUE
  )
  calibration$summary <- merge(
    calibration$summary,
    component_selected_summary,
    by = c("Lambda0", "Lambda1"),
    sort = TRUE
  )
  calibration$summary <- merge(
    calibration$summary,
    component_aggregate_summary,
    by = c("Lambda0", "Lambda1"),
    sort = TRUE
  )
  rejected <- calibration$summary$Lambda0[
    calibration$summary$ordered_reject
  ]
  cvm_rejected <- calibration$summary$Lambda0[
    calibration$summary$cvm_ordered_reject
  ]
  selected_rejected <- calibration$summary$Lambda0[
    calibration$summary$selected_ordered_reject
  ]
  component_selected_rejected <- calibration$summary$Lambda0[
    calibration$summary$component_selected_ordered_reject
  ]
  component_aggregate_rejected <- calibration$summary$Lambda0[
    calibration$summary$component_aggregate_ordered_reject
  ]
  list(
    summary = calibration$summary,
    moment_detail = calibration$moment_detail,
    selected_moment_detail = selected_calibration$moment_detail,
    component_selected_moment_detail =
      component_selected_calibration$moment_detail,
    component_aggregate_moment_detail =
      component_aggregate_calibration$moment_detail,
    largest_ordered_refuted_grid_value = if (length(rejected)) {
      max(rejected)
    } else NA_real_,
    largest_cvm_ordered_refuted_grid_value = if (length(cvm_rejected)) {
      max(cvm_rejected)
    } else NA_real_,
    largest_selected_ordered_refuted_grid_value =
      if (length(selected_rejected)) max(selected_rejected) else NA_real_,
    largest_component_selected_ordered_refuted_grid_value =
      if (length(component_selected_rejected)) {
        max(component_selected_rejected)
      } else NA_real_,
    largest_component_aggregate_ordered_refuted_grid_value =
      if (length(component_aggregate_rejected)) {
        max(component_aggregate_rejected)
      } else NA_real_,
    witness_diagnostics = do.call(rbind, diagnostic_rows),
    aggregate_moment_detail = aggregate_calibration$moment_detail,
    gamma_diagnostics = do.call(rbind, gamma_diagnostics),
    gamma_cv = if (length(gamma_cv_diagnostics)) {
      do.call(rbind, gamma_cv_diagnostics)
    } else data.frame(),
    gamma_cv_resamples = if (length(gamma_cv_resample_diagnostics)) {
      do.call(rbind, gamma_cv_resample_diagnostics)
    } else data.frame(),
    probability_diagnostics = if (length(probability_diagnostics)) {
      do.call(rbind, probability_diagnostics)
    } else data.frame(),
    fold_id = fold_id,
    settings = list(
      lambda_grid = lambda_grid,
      folds = folds,
      primitive_learner = primitive_learner,
      super_learner_candidates = super_learner_candidates,
      super_learner_cv_folds = super_learner_cv_folds,
      super_learner_ranger_num_trees = super_learner_ranger_num_trees,
      gamma_learner = gamma_learner,
      gamma_cv_selection_rule = gamma_cv_selection_rule,
      gamma_cv_folds = gamma_cv_folds,
      gamma_cv_repeats = gamma_cv_repeats,
      gamma_discrete_support = gamma_discrete_support,
      gamma_discrete_df_grid = gamma_discrete_df_grid,
      gamma_discrete_cv_folds = gamma_discrete_cv_folds,
      gamma_discrete_pmf_epsilon = gamma_discrete_pmf_epsilon,
      gamma_discrete_decay = gamma_discrete_decay,
      bootstrap_replications = bootstrap_replications,
      confidence_level = confidence_level,
      clustered_inference = clustered_inference,
      cluster_count = length(unique(cluster_id)),
      progress = isTRUE(progress),
      inferential_scope = paste(
        "Honest cross-fitted learned-witness test of necessary compatibility",
        "conditions; not a uniform adaptive test over all witnesses."
      )
    )
  )
}
