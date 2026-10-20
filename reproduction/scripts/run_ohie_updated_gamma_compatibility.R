#!/usr/bin/env Rscript

# Formal learned-witness compatibility test for the binary OHIE ED outcome
# under the finalized direct-RKHS gamma specification.  The analysis uses one
# adult per lottery household, no survey weights, and the initial 11-variable
# adjustment set.  Historical probability-nuisance predictions can be reused,
# but gamma is re-estimated with the current tuning contract.

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)

sources <- c(
  "R/sharp_bounds_plugin.R",
  "R/binary_outcome_bounds.R",
  "R/binary_rkhs_gamma_nuisances.R",
  "R/conservative_compatibility_test.R",
  "R/learned_witness_compatibility_test.R"
)
input_file <- "data/ohie_application.rds"
if (!all(file.exists(c(sources, input_file)))) {
  stop("Run this script from the repository root after preparing OHIE data.")
}
invisible(lapply(sources, source))

packages <- c(
  "SuperLearner", "glmnet", "ranger", "xgboost", "mgcv",
  "hal9001", "nloptr", "nnls"
)
available <- vapply(packages, requireNamespace, logical(1L), quietly = TRUE)
if (!all(available)) {
  stop("Missing packages: ", paste(packages[!available], collapse = ", "))
}

save_atomic <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(".pending_", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(value, temporary, compress = "gzip")
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Could not save ", path, ".")
}

write_csv <- function(value, filename) {
  path <- file.path(output_directory, filename)
  temporary <- tempfile(".pending_", tmpdir = output_directory)
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(value, temporary, row.names = FALSE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Could not save ", path, ".")
}

prepared <- readRDS(input_file)
study_data <- prepared$ed_household_index_overlap
x_cols <- prepared$ed_overlap_expanded_x_cols
study_seed <- 20260814L
fold_id <- make_tz_stratified_folds(study_data, K = 5L, seed = study_seed)
candidates <- c(
  "glm_main", "glm_interactions", "glmnet_interactions",
  "ranger", "xgboost", "gam", "hal"
)
lambda_grid <- c(1, 1.1, 1.2, 1.4, 1.6, 1.8, 2)
lambda_constants <- c(
  0.03162278, 0.1, 0.31622777, 1, 3.16227766, 10
)
bandwidth_multiplier <- 1 / sqrt(2)
version <- "ed_household_index_overlap_expanded_gam_hal_direct_rkhs_half_median_minimum_r2_b999_lambda1_20_v1"
output_directory <- file.path(
  "results", "ohie_learned_witness_compatibility", version
)
checkpoint_directory <- file.path(output_directory, "checkpoints")
probability_cache_source <- file.path(
  "results", "ohie_learned_witness_compatibility",
  "ed_household_index_overlap_expanded_gam_hal_direct_rkhs_b999_grid15_v1",
  "checkpoints"
)

stopifnot(
  is.data.frame(study_data),
  nrow(study_data) == 22170L,
  length(x_cols) == 11L,
  identical(x_cols, c(
    "age_centered", "female", "baseline_ed_use", "household_size_2",
    paste0("lottery_draw_", 2:8)
  )),
  !anyDuplicated(study_data$person_id),
  !anyDuplicated(study_data$household_id),
  all(study_data$household_size_3 == 0),
  identical(
    fold_id,
    binary_make_tz_folds(study_data, K = 5L, seed = study_seed)
  )
)

runtime_configuration <- list(
  lambda_grid = lambda_grid,
  x_cols = x_cols,
  folds = 5L,
  primitive_learner = "binary_super_learner",
  gamma_learner = "rkhs_cv",
  positivity_epsilon = 0.05,
  treatment_spline_df = 3L,
  gam_basis_dimension = 5L,
  allocation_spline_df = 3L,
  super_learner_candidates = candidates,
  super_learner_cv_folds = 3L,
  super_learner_ranger_num_trees = 500L,
  super_learner_hal_max_degree = 2L,
  gamma_cv_lambda_grid = function(n_training, treatment) {
    lambda_constants / sqrt(n_training)
  },
  gamma_cv_bandwidth_multipliers = bandwidth_multiplier,
  gamma_cv_folds = 3L,
  gamma_cv_repeats = 2L,
  gamma_cv_selection_rule = "minimum",
  gamma_nystrom_rank = NULL,
  gamma_kernel = "gaussian",
  witness_shapes = c("indicator", "linear", "sqrt"),
  bootstrap_replications = 999L,
  confidence_level = 0.95,
  seed = study_seed
)
contract_configuration <- runtime_configuration
contract_configuration$gamma_cv_lambda_grid <- "C/sqrt(n_t)"
contract <- list(
  schema_version = 1L,
  version = version,
  outcome = "any emergency-department use during the study period",
  sample = "one deterministic outcome-blind index adult per household",
  weighting = "survey weights ignored",
  covariate_specification = "initial 11-variable adjustment set",
  data_md5 = unname(tools::md5sum(input_file)),
  source_md5 = tools::md5sum(sources),
  person_id = study_data$person_id,
  fold_id = fold_id,
  configuration = contract_configuration,
  gamma_method = "normalized direct RKHS/Riesz density ratio",
  gamma_lambda_constants = lambda_constants,
  gamma_bandwidth_rule = paste(
    "sigma^2=M/2, with M the training-fold median pairwise squared",
    "distance in standardized (Y,X)"
  ),
  gamma_bandwidth_multiplier = bandwidth_multiplier,
  gamma_cv_folds = 3L,
  gamma_cv_repeats = 2L,
  gamma_selection_rule = "minimum mean Riesz validation risk",
  gamma_rank_rule = "min(n_t,ceiling(10*sqrt(n_t)))",
  package_versions = vapply(
    packages, function(package) as.character(packageVersion(package)),
    character(1L)
  )
)

if (identical(Sys.getenv("IVMSM_OHIE_UPDATED_COMPAT_PREFLIGHT"), "1")) {
  cat(
    "PASS: binary OHIE sample, initial covariates, GAM/HAL library, ",
    "fixed half-median-squared gamma bandwidth, two repeated three-fold ",
    "CV partitions, and seven-value compatibility grid.\n",
    sep = ""
  )
  quit(save = "no", status = 0L, runLast = FALSE)
}

dir.create(checkpoint_directory, recursive = TRUE, showWarnings = FALSE)
contract_path <- file.path(output_directory, "run_contract.rds")
if (file.exists(contract_path)) {
  if (!identical(readRDS(contract_path), contract)) {
    stop("The compatibility output directory has a different contract.")
  }
} else {
  save_atomic(contract, contract_path)
}

restore_prediction_fit <- function(checkpoint) {
  value <- checkpoint$metadata
  prediction <- checkpoint$prediction
  reference_id <- contract$person_id
  value$predict <- function(newdata) {
    index <- match(newdata$person_id, reference_id)
    if (anyNA(index)) {
      stop("A cached prediction received an unknown person identifier.")
    }
    if (is.list(prediction)) {
      lapply(prediction, function(component) component[index])
    } else {
      prediction[index]
    }
  }
  value
}

# Probability nuisances are unchanged from the earlier GAM/HAL run. Reuse
# them only if their complete call-level contract matches; otherwise refit and
# save a new checkpoint in the current output directory.
component_progress <- list()
base_binary_fit <- binary_fit_regression_with_learner
binary_fit_regression_with_learner <- function(
    data, response, x_cols, spline_df, ..., seed = 1L
) {
  arguments <- list(...)
  local_contract <- list(
    training_ids = data$person_id,
    response = response,
    x_cols = x_cols,
    spline_df = spline_df,
    arguments = arguments,
    seed = seed
  )
  filename <- sprintf("probability_seed%d.rds", seed)
  source_path <- file.path(probability_cache_source, filename)
  current_path <- file.path(checkpoint_directory, filename)
  checkpoint <- NULL
  cache_source <- "new"
  for (candidate_path in c(current_path, source_path)) {
    if (!file.exists(candidate_path)) next
    candidate_checkpoint <- readRDS(candidate_path)
    if (identical(candidate_checkpoint$contract, local_contract)) {
      checkpoint <- candidate_checkpoint
      cache_source <- if (identical(candidate_path, current_path)) {
        "current"
      } else {
        "historical_probability_only"
      }
      break
    }
  }
  if (is.null(checkpoint)) {
    started <- proc.time()[["elapsed"]]
    value <- do.call(base_binary_fit, c(list(
      data = data,
      response = response,
      x_cols = x_cols,
      spline_df = spline_df,
      seed = seed
    ), arguments))
    failed <- any(as.logical(unlist(
      value$super_learner_fit$errorsInLibrary
    ))) || any(as.logical(unlist(
      value$super_learner_fit$errorsInCVLibrary
    )))
    if (failed || !identical(value$cv_table$candidate, candidates) ||
        any(!is.finite(value$cv_table$mean_risk)) ||
        !isTRUE(value$converged)) {
      stop("A requested probability Super Learner candidate failed.")
    }
    prediction <- value$predict(study_data)
    stopifnot(
      length(prediction) == nrow(study_data),
      all(is.finite(prediction))
    )
    fields <- c(
      "converged", "rank", "constant", "learner", "selected_learner",
      "cv_table", "ensemble_weights"
    )
    checkpoint <- list(
      contract = local_contract,
      metadata = value[intersect(fields, names(value))],
      prediction = prediction,
      elapsed_seconds = proc.time()[["elapsed"]] - started
    )
    save_atomic(checkpoint, current_path)
  }
  progress <- checkpoint$metadata$cv_table
  progress$component_seed <- seed
  progress$outer_fold <- as.integer((seed - study_seed) %/% 10000L)
  progress$response <- response
  progress$training_n <- nrow(data)
  progress$elapsed_seconds <- checkpoint$elapsed_seconds
  progress$cache_source <- cache_source
  component_progress[[as.character(seed)]] <<- progress
  write_csv(do.call(rbind, component_progress), "component_fit_progress.csv")
  message(
    "Probability nuisance ", length(component_progress), "/15: ",
    response, "; cache=", cache_source, "."
  )
  restore_prediction_fit(checkpoint)
}

# Gamma requires a new fit because the earlier compatibility run used a fixed
# penalty and sigma^2=M. Each current fold/arm fit is independently resumable.
gamma_progress <- list()
base_gamma_fit <- fit_gamma_rkhs
fit_gamma_rkhs <- function(
    train, nuisance_train, t_value, x_cols = "X", lambda_grid,
    bandwidth_multipliers, cv_folds, cv_repeats = 2L,
    selection_rule, ..., seed = 1L
) {
  arguments <- list(...)
  resolved_lambda_grid <- resolve_gamma_lambda_grid(
    lambda_grid,
    n_training = sum(train$T_obs == t_value),
    treatment = t_value
  )
  local_contract <- list(
    training_ids = train$person_id,
    nuisance_train = nuisance_train,
    treatment = t_value,
    x_cols = x_cols,
    lambda_grid = resolved_lambda_grid,
    bandwidth_multipliers = bandwidth_multipliers,
    cv_folds = cv_folds,
    cv_repeats = cv_repeats,
    selection_rule = selection_rule,
    arguments = arguments,
    seed = seed
  )
  path <- file.path(
    checkpoint_directory,
    sprintf("gamma_cv_seed%d_t%d.rds", seed, t_value)
  )
  if (file.exists(path)) {
    checkpoint <- readRDS(path)
    if (!identical(checkpoint$contract, local_contract)) {
      stop("Gamma checkpoint contract mismatch: ", basename(path), ".")
    }
    cache_reused <- TRUE
  } else {
    started <- proc.time()[["elapsed"]]
    value <- do.call(base_gamma_fit, c(list(
      train = train,
      nuisance_train = nuisance_train,
      t_value = t_value,
      x_cols = x_cols,
      lambda_grid = resolved_lambda_grid,
      bandwidth_multipliers = bandwidth_multipliers,
      cv_folds = cv_folds,
      cv_repeats = cv_repeats,
      selection_rule = selection_rule,
      seed = seed
    ), arguments))
    prediction <- value$predict(study_data)
    stopifnot(all(vapply(
      prediction,
      function(component) {
        length(component) == nrow(study_data) && all(is.finite(component))
      },
      logical(1L)
    )))
    checkpoint <- list(
      contract = local_contract,
      metadata = value[setdiff(names(value), "predict")],
      prediction = prediction,
      elapsed_seconds = proc.time()[["elapsed"]] - started
    )
    save_atomic(checkpoint, path)
    cache_reused <- FALSE
  }
  row <- data.frame(
    outer_fold = as.integer((seed - study_seed) %/% 10000L),
    treatment = t_value,
    training_n = sum(train$T_obs == t_value),
    selected_lambda = checkpoint$metadata$selected_lambda,
    selected_C = checkpoint$metadata$selected_lambda *
      sqrt(sum(train$T_obs == t_value)),
    bandwidth_multiplier =
      checkpoint$metadata$selected_bandwidth_multiplier,
    base_bandwidth = checkpoint$metadata$base_bandwidth,
    selected_bandwidth = checkpoint$metadata$selected_bandwidth,
    nystrom_landmarks = checkpoint$metadata$nystrom_landmark_count,
    nystrom_effective_rank = checkpoint$metadata$nystrom_effective_rank,
    cv_repeats = cv_repeats,
    cache_reused = cache_reused,
    elapsed_seconds = checkpoint$elapsed_seconds,
    stringsAsFactors = FALSE
  )
  gamma_progress[[paste(seed, t_value, sep = "_")]] <<- row
  write_csv(do.call(rbind, gamma_progress), "gamma_fit_progress.csv")
  message(
    "Gamma fit ", length(gamma_progress), "/10: fold=", row$outer_fold,
    ", treatment=", t_value, ", C=", format(row$selected_C, digits = 5),
    ", cache=", cache_reused, "."
  )
  restore_prediction_fit(checkpoint)
}

completion_file <- file.path(output_directory, "ohie_compatibility_audit.rds")
if (file.exists(completion_file)) {
  result <- readRDS(completion_file)
  if (!identical(result$contract, contract)) {
    stop("The completed compatibility result has a different contract.")
  }
  message("Reusing completed updated-gamma compatibility analysis.")
} else {
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    do.call(
      test_learned_witness_dual_msm_compatibility,
      c(list(data = study_data), runtime_configuration)
    ),
    error = function(error) {
      save_atomic(
        list(
          success = FALSE,
          error = conditionMessage(error),
          completed_at = as.character(Sys.time())
        ),
        file.path(output_directory, "run_failure.rds")
      )
      stop(error)
    }
  )
  result <- list(
    success = TRUE,
    contract = contract,
    fit = fit,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    completed_at = as.character(Sys.time()),
    session_info = utils::sessionInfo()
  )
  save_atomic(result, completion_file)
}

fit <- result$fit
write_csv(fit$summary, "compatibility_test.csv")
write_csv(fit$moment_detail, "moment_detail.csv")
write_csv(fit$selected_moment_detail, "selected_moment_detail.csv")
write_csv(fit$component_selected_moment_detail,
          "component_selected_moment_detail.csv")
write_csv(fit$component_aggregate_moment_detail,
          "component_aggregate_moment_detail.csv")
write_csv(fit$witness_diagnostics, "witness_diagnostics.csv")
write_csv(fit$gamma_diagnostics, "gamma_nuisance_diagnostics.csv")
write_csv(fit$gamma_cv, "gamma_cv.csv")
if (is.data.frame(fit$gamma_cv_resamples) &&
    nrow(fit$gamma_cv_resamples)) {
  write_csv(fit$gamma_cv_resamples, "gamma_cv_resamples.csv")
}
write_csv(fit$probability_diagnostics,
          "probability_nuisance_diagnostics.csv")
write_csv(
  prepared$sample_summary[
    prepared$sample_summary$sample == "ed_household_index_overlap",
    , drop = FALSE
  ],
  "sample_summary.csv"
)
if (file.exists(file.path(output_directory, "component_fit_progress.csv"))) {
  progress <- utils::read.csv(file.path(
    output_directory, "component_fit_progress.csv"
  ))
  write_csv(
    stats::aggregate(coefficient ~ candidate, progress, mean),
    "mean_ensemble_weights.csv"
  )
}

columns <- c(
  "Lambda0", "Lambda1", "statistic", "p_value", "reject",
  "ordered_reject", "grid_simultaneous_p_value",
  "grid_simultaneous_reject"
)
print(fit$summary[, columns], row.names = FALSE)
cat(
  "Completed updated-gamma OHIE compatibility tests in ",
  round(result$elapsed_seconds / 60, 1), " minutes.\nResults: ",
  output_directory, "\n",
  sep = ""
)
