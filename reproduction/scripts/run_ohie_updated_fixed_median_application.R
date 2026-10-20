#!/usr/bin/env Rscript

# Oregon Health Insurance Experiment application with the generic direct-RKHS
# gamma estimator. The Gaussian-kernel bandwidth is a prespecified multiple of
# the fold-specific median pairwise distance of standardized (Y,X); the
# multiplier and RKHS penalty are tuned jointly by repeated inner CV.
# This is a complete, standalone analysis driver: it loads the prepared data,
# fits all nuisances once, evaluates exact and epsilon-smoothed Dual-MSM bounds,
# and writes auditable diagnostics.

source_files <- c(
  "R/binary_outcome_bounds.R",
  "R/sharp_bounds_plugin.R",
  "R/binary_rkhs_gamma_nuisances.R"
)
covariate_specification <- tolower(trimws(Sys.getenv(
  "IVMSM_OHIE_COVARIATE_SPEC", unset = "expanded"
)))
if (!covariate_specification %in% c("expanded", "augmented")) {
  stop("IVMSM_OHIE_COVARIATE_SPEC must be expanded or augmented.")
}
input_file <- trimws(Sys.getenv(
  "IVMSM_OHIE_INPUT_FILE",
  unset = if (covariate_specification == "augmented") {
    "data/ohie_application_augmented.rds"
  } else {
    "data/ohie_application.rds"
  }
))
if (!all(file.exists(c(input_file, source_files)))) {
  stop("Run this script from the repository root after preparing OHIE data.")
}
for (path in source_files) source(path)

bind_rows <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(row) {
    row <- as.data.frame(row, stringsAsFactors = FALSE)
    for (column in setdiff(columns, names(row))) row[[column]] <- NA
    row[, columns, drop = FALSE]
  })
  value <- do.call(rbind, rows)
  rownames(value) <- NULL
  value
}

save_atomic <- function(object, path) {
  temporary <- tempfile(pattern = ".pending_", tmpdir = dirname(path))
  saveRDS(object, temporary, compress = "gzip")
  if (!file.rename(temporary, path)) stop("Could not save ", path, ".")
}

write_table <- function(value, filename) {
  utils::write.csv(
    value,
    file.path(output_directory, filename),
    row.names = FALSE
  )
}

lambda_text <- trimws(Sys.getenv(
  "IVMSM_OHIE_LAMBDAS", unset = "1.2,1.4,1.6,1.8,2.0"
))
Lambda_grid <- suppressWarnings(as.numeric(strsplit(
  lambda_text, ",", fixed = TRUE
)[[1L]]))
if (!length(Lambda_grid) || any(!is.finite(Lambda_grid)) ||
    any(Lambda_grid < 1)) {
  stop("IVMSM_OHIE_LAMBDAS must contain finite values of at least one.")
}
Lambda_grid <- unique(Lambda_grid)
epsilon_grid <- c(0, 0.01)
confidence_level <- 0.95
folds <- 5L
seed <- 20260814L
exclude_hal <- tolower(trimws(Sys.getenv(
  "IVMSM_OHIE_EXCLUDE_HAL", unset = "false"
))) %in% c("1", "true", "yes", "y")
exclude_gam <- tolower(trimws(Sys.getenv(
  "IVMSM_OHIE_EXCLUDE_GAM", unset = "false"
))) %in% c("1", "true", "yes", "y")
rkhs_selection_rule <- tolower(trimws(Sys.getenv(
  "IVMSM_OHIE_RKHS_SELECTION_RULE", unset = "minimum"
)))
if (!rkhs_selection_rule %in% c("minimum", "one_se")) {
  stop("IVMSM_OHIE_RKHS_SELECTION_RULE must be minimum or one_se.")
}
rkhs_cv_repeats <- suppressWarnings(as.integer(trimws(Sys.getenv(
  "IVMSM_OHIE_RKHS_CV_REPEATS", unset = "2"
))))
if (
  length(rkhs_cv_repeats) != 1L || is.na(rkhs_cv_repeats) ||
    rkhs_cv_repeats < 1L
) {
  stop("IVMSM_OHIE_RKHS_CV_REPEATS must be one positive integer.")
}
nystrom_sqrt_n_multiplier <- suppressWarnings(as.numeric(trimws(Sys.getenv(
  "IVMSM_OHIE_NYSTROM_SQRT_N_MULTIPLIER", unset = "10"
))))
if (
  length(nystrom_sqrt_n_multiplier) != 1L ||
    !is.finite(nystrom_sqrt_n_multiplier) ||
    nystrom_sqrt_n_multiplier <= 0
) {
  stop("IVMSM_OHIE_NYSTROM_SQRT_N_MULTIPLIER must be positive and finite.")
}
library_label <- if (exclude_hal && exclude_gam) {
  "no_gam_no_hal"
} else if (exclude_hal) {
  "gam_no_hal"
} else if (exclude_gam) {
  "hal_no_gam"
} else {
  "gam_hal"
}
default_version <- if (
  rkhs_selection_rule == "minimum" && rkhs_cv_repeats == 1L &&
    nystrom_sqrt_n_multiplier == 10
) {
  paste0(
    "ohie_ed_", covariate_specification, "_", library_label,
    "_gamma_fixed_median_lambda12_20_v1"
  )
} else {
  paste0(
    "ohie_ed_", covariate_specification, "_", library_label,
    "_gamma_fixed_median_",
    rkhs_selection_rule, "_rank",
    gsub("\\.", "p", format(nystrom_sqrt_n_multiplier, trim = TRUE)),
    "_r", rkhs_cv_repeats, "_lambda12_20_v1"
  )
}
version <- trimws(Sys.getenv(
  "IVMSM_OHIE_UPDATED_VERSION",
  unset = default_version
))
if (!nzchar(version)) stop("IVMSM_OHIE_UPDATED_VERSION cannot be empty.")
reuse_nuisance_checkpoint <- trimws(Sys.getenv(
  "IVMSM_OHIE_REUSE_NUISANCE_CHECKPOINT", unset = ""
))
if (nzchar(reuse_nuisance_checkpoint) &&
    !file.exists(reuse_nuisance_checkpoint)) {
  stop(
    "IVMSM_OHIE_REUSE_NUISANCE_CHECKPOINT does not exist: ",
    reuse_nuisance_checkpoint, "."
  )
}

prepared <- readRDS(input_file)
if (covariate_specification == "augmented") {
  data <- prepared$ed_household_index_overlap_augmented
  x_cols <- prepared$ed_overlap_augmented_x_cols
  expected_covariates <- 20L
} else {
  data <- prepared$ed_household_index_overlap
  x_cols <- prepared$ed_overlap_expanded_x_cols
  expected_covariates <- 11L
}
if (
  !is.data.frame(data) || nrow(data) != 22170L ||
    length(x_cols) != expected_covariates ||
    anyDuplicated(data$person_id) || anyDuplicated(data$household_id) ||
    !all(c("Y", "T_obs", "Z", x_cols) %in% names(data)) ||
    anyNA(data[, c("Y", "T_obs", "Z", x_cols)]) ||
    any(data$household_size_3 != 0)
) {
  stop("The prepared OHIE analysis sample does not satisfy its frozen contract.")
}

candidates <- c(
  "glm_main",
  "glm_interactions",
  "glmnet_interactions",
  "ranger",
  "xgboost"
)
candidates <- c(candidates, if (!exclude_gam) "gam")
candidates <- c(candidates, if (!exclude_hal) "hal")
required_packages <- c(
  "glmnet", "ranger", "xgboost", "nnls",
  if (!exclude_gam) "mgcv",
  if (!exclude_hal) "hal9001",
  "SuperLearner", "nloptr"
)
available <- vapply(
  required_packages, requireNamespace, logical(1L), quietly = TRUE
)
if (!all(available)) {
  stop("Missing packages: ", paste(required_packages[!available], collapse = ", "))
}

gamma_lambda_constants <- 10^seq(-1.5, 1, by = 0.5)
gamma_lambda_grid <- function(n_training, treatment) {
  gamma_lambda_constants / sqrt(n_training)
}
bandwidth_text <- trimws(Sys.getenv(
  "IVMSM_OHIE_RKHS_BANDWIDTH_MULTIPLIERS",
  unset = as.character(1 / sqrt(2))
))
gamma_bandwidth_multipliers <- suppressWarnings(as.numeric(strsplit(
  bandwidth_text, ",", fixed = TRUE
)[[1L]]))
if (
  !length(gamma_bandwidth_multipliers) ||
    any(!is.finite(gamma_bandwidth_multipliers)) ||
    any(gamma_bandwidth_multipliers <= 0)
) {
  stop(
    "IVMSM_OHIE_RKHS_BANDWIDTH_MULTIPLIERS must contain positive finite values."
  )
}
gamma_bandwidth_multipliers <- sort(unique(gamma_bandwidth_multipliers))
# The base bandwidth is the median pairwise distance of standardized (Y,X),
# recomputed in every training fold.  The production multiplier 1/sqrt(2)
# fixes sigma^2 at one half of the median squared distance, for which the
# median-pair Gaussian-kernel similarity is exp(-1).

configuration <- list(
  K = folds,
  spline_df = 5L,
  primitive_learner = "super_learner",
  overlap_epsilon = 0.05,
  outcome_probability_epsilon = 1e-4,
  gam_min_basis = 5L,
  gam_max_basis = 20L,
  gam_growth_exponent = 1 / 3,
  library_candidates = c("glm_spline", "gam", "glmnet", "ranger"),
  library_cv_folds = 3L,
  library_selection_rule = "minimum",
  library_ranger_num_trees = 500L,
  super_learner_candidates = candidates,
  super_learner_cv_folds = 3L,
  super_learner_ranger_num_trees = 500L,
  super_learner_hal_max_degree = 2L,
  rkhs_lambda_by_t = NULL,
  rkhs_lambda_grid = gamma_lambda_grid,
  rkhs_bandwidth_multipliers = gamma_bandwidth_multipliers,
  rkhs_cv_folds = 3L,
  rkhs_cv_repeats = rkhs_cv_repeats,
  rkhs_selection_rule = rkhs_selection_rule,
  rkhs_bandwidth_multiplier = 1 / sqrt(2),
  rkhs_backend = "nystrom",
  rkhs_nystrom_rank = if (nystrom_sqrt_n_multiplier == 10) {
    NULL
  } else {
    function(n_training) {
      ceiling(nystrom_sqrt_n_multiplier * sqrt(n_training))
    }
  },
  rkhs_kernel = "gaussian",
  seed = seed
)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)

output_directory <- file.path("results", "ohie_application", version)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
fold_id <- binary_make_tz_folds(data, K = folds, seed = seed)
contract_configuration <- configuration
contract_configuration$rkhs_lambda_grid <- "C/sqrt(n_t)"
contract_configuration$rkhs_nystrom_rank <- sprintf(
  "min(n_t,ceiling(%g*sqrt(n_t)))", nystrom_sqrt_n_multiplier
)
contract <- list(
  version = version,
  data_md5 = unname(tools::md5sum(input_file)),
  source_md5 = tools::md5sum(source_files),
  observations = nrow(data),
  x_cols = x_cols,
  person_id = data$person_id,
  fold_id = fold_id,
  Lambda_grid = Lambda_grid,
  epsilon_grid = epsilon_grid,
  confidence_level = confidence_level,
  configuration = contract_configuration,
  gamma_lambda_constants = gamma_lambda_constants,
  gamma_bandwidth_multipliers = gamma_bandwidth_multipliers,
  package_versions = vapply(required_packages, function(package) {
    as.character(utils::packageVersion(package))
  }, character(1L))
)

primitive_cache_directory <- file.path(
  "results", "ohie_application",
  if (exclude_hal && exclude_gam) {
    paste0("ohie_ed_", covariate_specification,
           "_no_gam_no_hal_primitive_cache_v1")
  } else if (exclude_hal) {
    paste0("ohie_ed_", covariate_specification,
           "_gam_no_hal_primitive_cache_v1")
  } else if (exclude_gam) {
    paste0("ohie_ed_", covariate_specification,
           "_hal_no_gam_primitive_cache_v1")
  } else {
    paste0("ohie_ed_", covariate_specification,
           "_gam_hal_primitive_cache_v1")
  }
)
dir.create(primitive_cache_directory, recursive = TRUE, showWarnings = FALSE)
primitive_cache_contract <- list(
  schema_version = 1L,
  data_md5 = contract$data_md5,
  observations = nrow(data),
  person_id = data$person_id,
  fold_id = fold_id,
  x_cols = x_cols,
  seed = seed,
  candidates = candidates,
  spline_df = configuration$spline_df,
  gam_min_basis = configuration$gam_min_basis,
  gam_max_basis = configuration$gam_max_basis,
  gam_growth_exponent = configuration$gam_growth_exponent,
  super_learner_cv_folds = configuration$super_learner_cv_folds,
  super_learner_ranger_num_trees =
    configuration$super_learner_ranger_num_trees,
  super_learner_hal_max_degree =
    configuration$super_learner_hal_max_degree,
  package_versions = contract$package_versions
)
primitive_cache_manifest <- file.path(
  primitive_cache_directory, "contract.rds"
)
if (file.exists(primitive_cache_manifest)) {
  cached_contract <- readRDS(primitive_cache_manifest)
  # Historical manifests included a hash of the combined primitive-and-gamma
  # source files. Gamma-only implementation changes do not invalidate the
  # cached HAL/GAM primitive predictions.
  cached_contract$source_md5 <- NULL
  if (!identical(cached_contract, primitive_cache_contract)) {
    stop("The shared primitive-fit cache has a different analysis contract.")
  }
} else {
  save_atomic(primitive_cache_contract, primitive_cache_manifest)
}
contract$primitive_cache_contract <- primitive_cache_contract

hydrate_primitive_fit <- function(record) {
  output <- record$fit_summary
  reference_id <- data$person_id
  full_prediction <- record$full_prediction
  output$predict <- function(newdata) {
    index <- match(newdata$person_id, reference_id)
    if (anyNA(index)) stop("Cached primitive prediction received an unknown ID.")
    full_prediction[index]
  }
  output
}

component_progress <- list()
gamma_progress <- list()
component_call_index <- 0L
base_component_fitter <- binary_rkhs_gamma_fit_component
binary_rkhs_gamma_fit_component <- function(
    training, response, x_cols, seed, ...
) {
  local_started <- proc.time()[["elapsed"]]
  component_call_index <<- component_call_index + 1L
  component_cache_file <- file.path(
    primitive_cache_directory,
    sprintf("component_%02d.rds", component_call_index)
  )
  component_metadata <- list(
    index = component_call_index,
    response = response,
    x_cols = x_cols,
    seed = seed,
    training_person_id = training$person_id
  )
  cache_reused <- file.exists(component_cache_file)
  if (cache_reused) {
    record <- readRDS(component_cache_file)
    if (!identical(record$metadata, component_metadata)) {
      stop("Primitive-fit cache metadata mismatch at component ",
           component_call_index, ".")
    }
    candidate_failed <- isTRUE(record$candidate_failed)
    value <- hydrate_primitive_fit(record)
  } else {
    fitted <- base_component_fitter(
      training = training,
      response = response,
      x_cols = x_cols,
      seed = seed,
      ...
    )
    full_errors <- fitted$super_learner_fit$errorsInLibrary
    cv_errors <- fitted$super_learner_fit$errorsInCVLibrary
    candidate_failed <- any(as.logical(unlist(full_errors))) ||
      any(as.logical(unlist(cv_errors))) ||
      any(!is.finite(fitted$cv_table$mean_risk))
    if (candidate_failed) {
      stop("A requested Super Learner candidate failed.")
    }
    full_prediction <- fitted$predict(data)
    if (length(full_prediction) != nrow(data) ||
        any(!is.finite(full_prediction))) {
      stop("A primitive learner produced invalid reusable predictions.")
    }
    record <- list(
      metadata = component_metadata,
      candidate_failed = candidate_failed,
      full_prediction = as.numeric(full_prediction),
      fit_summary = fitted[setdiff(
        names(fitted), c("predict", "super_learner_fit")
      )]
    )
    save_atomic(record, component_cache_file)
    value <- hydrate_primitive_fit(record)
    rm(fitted)
  }
  table <- value$cv_table
  table$response <- response
  table$component_seed <- seed
  table$training_n <- nrow(training)
  table$elapsed_seconds <- proc.time()[["elapsed"]] - local_started
  table$any_candidate_failed <- candidate_failed
  table$cache_reused <- cache_reused
  component_progress[[length(component_progress) + 1L]] <<- table
  write_table(bind_rows(component_progress), "component_fit_progress.csv")
  if (candidate_failed) {
    stop("A requested Super Learner candidate failed; see component progress.")
  }
  message(sprintf(
    "Completed primitive nuisance %d/25 (%s; n=%d; cache=%s) in %.1f seconds.",
    length(component_progress), response, nrow(training),
    if (cache_reused) "reused" else "new",
    proc.time()[["elapsed"]] - local_started
  ))
  value
}

base_gamma_fitter <- fit_gamma_rkhs
fit_gamma_rkhs <- function(
    train, nuisance_train, t_value, seed, ...
) {
  local_started <- proc.time()[["elapsed"]]
  message(sprintf(
    "Starting joint gamma penalty/bandwidth CV for treatment %d (training n=%d).",
    t_value, sum(train$T_obs == t_value)
  ))
  value <- base_gamma_fitter(
    train = train,
    nuisance_train = nuisance_train,
    t_value = t_value,
    seed = seed,
    ...
  )
  row <- data.frame(
    treatment = t_value,
    gamma_seed = seed,
    training_n = sum(train$T_obs == t_value),
    selected_lambda = value$selected_lambda,
    selected_lambda_sqrt_n_constant =
      value$selected_lambda * sqrt(sum(train$T_obs == t_value)),
    selected_bandwidth_multiplier = value$selected_bandwidth_multiplier,
    base_bandwidth = value$base_bandwidth,
    selected_bandwidth = value$selected_bandwidth,
    nystrom_landmark_count = value$nystrom_landmark_count,
    nystrom_effective_rank = value$nystrom_effective_rank,
    elapsed_seconds = proc.time()[["elapsed"]] - local_started,
    stringsAsFactors = FALSE
  )
  gamma_progress[[length(gamma_progress) + 1L]] <<- row
  write_table(bind_rows(gamma_progress), "gamma_fit_progress.csv")
  message(sprintf(
    paste0(
      "Completed gamma fit %d/10 in %.1f seconds; ",
      "lambda=%.6g, bandwidth multiplier=%g."
    ),
    length(gamma_progress), row$elapsed_seconds,
    row$selected_lambda, row$selected_bandwidth_multiplier
  ))
  value
}

cache_file <- file.path(output_directory, "nuisance_checkpoint.rds")
if (file.exists(cache_file)) {
  cache <- readRDS(cache_file)
  if (!identical(cache$contract, contract)) {
    stop("The nuisance checkpoint has a different analysis contract.")
  }
  nuisance_fit <- cache$nuisance_fit
  nuisance_seconds <- cache$nuisance_seconds
  message("Reusing the validated updated-estimator nuisance checkpoint.")
} else if (nzchar(reuse_nuisance_checkpoint)) {
  source_cache <- readRDS(reuse_nuisance_checkpoint)
  source_contract <- source_cache$contract
  configuration_keys <- c(
    "K", "spline_df", "primitive_learner", "overlap_epsilon",
    "outcome_probability_epsilon", "gam_min_basis", "gam_max_basis",
    "gam_growth_exponent", "library_candidates", "library_cv_folds",
    "library_selection_rule", "library_ranger_num_trees",
    "super_learner_candidates", "super_learner_cv_folds",
    "super_learner_ranger_num_trees", "super_learner_hal_max_degree",
    "rkhs_lambda_grid", "rkhs_bandwidth_multipliers", "rkhs_cv_folds",
    "rkhs_cv_repeats", "rkhs_selection_rule", "rkhs_backend",
    "rkhs_nystrom_rank", "rkhs_kernel", "seed"
  )
  source_signature <- list(
    data_md5 = source_contract$data_md5,
    observations = source_contract$observations,
    x_cols = source_contract$x_cols,
    person_id = source_contract$person_id,
    fold_id = source_contract$fold_id,
    configuration = source_contract$configuration[configuration_keys],
    gamma_lambda_constants = source_contract$gamma_lambda_constants,
    gamma_bandwidth_multipliers =
      source_contract$gamma_bandwidth_multipliers,
    package_versions = source_contract$package_versions
  )
  current_signature <- list(
    data_md5 = contract$data_md5,
    observations = contract$observations,
    x_cols = contract$x_cols,
    person_id = contract$person_id,
    fold_id = contract$fold_id,
    configuration = contract$configuration[configuration_keys],
    gamma_lambda_constants = contract$gamma_lambda_constants,
    gamma_bandwidth_multipliers = contract$gamma_bandwidth_multipliers,
    package_versions = contract$package_versions
  )
  signatures_match <- isTRUE(all.equal(
    source_signature,
    current_signature,
    tolerance = 1e-12,
    check.attributes = TRUE
  ))
  if (!signatures_match) {
    mismatch <- names(source_signature)[!vapply(
      names(source_signature),
      function(name) isTRUE(all.equal(
        source_signature[[name]], current_signature[[name]],
        tolerance = 1e-12, check.attributes = TRUE
      )),
      logical(1L)
    )]
    stop(
      "The requested reusable nuisance checkpoint does not have the same ",
      "data, folds, learners, or gamma-tuning contract. Mismatched fields: ",
      paste(mismatch, collapse = ", "), "."
    )
  }
  nuisance_fit <- source_cache$nuisance_fit
  nuisance_seconds <- source_cache$nuisance_seconds
  if (
    !identical(nuisance_fit$fold_id, fold_id) ||
      nrow(nuisance_fit$predictions) != nrow(data) ||
      any(!is.finite(as.matrix(nuisance_fit$predictions))) ||
      any(abs(
        nuisance_fit$gamma_cv_table$bandwidth_multiplier - 1 / sqrt(2)
      ) > 1e-12) ||
      any(nuisance_fit$gamma_cv_table$cv_repeats != 2L)
  ) {
    stop("The reusable nuisance object failed the current gamma validation.")
  }
  save_atomic(
    list(
      contract = contract,
      nuisance_fit = nuisance_fit,
      nuisance_seconds = nuisance_seconds,
      reused_from = normalizePath(reuse_nuisance_checkpoint),
      completed_at = as.character(Sys.time())
    ),
    cache_file
  )
  message(
    "Reused and revalidated the nuisance fit from: ",
    reuse_nuisance_checkpoint
  )
} else {
  message(
    paste(
      sprintf("Fitting OHIE nuisances with %d-learner Super Learner,",
              length(candidates)),
      sprintf(
        "scaled-median bandwidth grid {%s}, %s selection, %d CV repeat(s),",
        paste(gamma_bandwidth_multipliers, collapse = ","),
        rkhs_selection_rule, rkhs_cv_repeats
      ),
      sprintf("and Nyström rank ceiling(%g*sqrt(n_t)).",
              nystrom_sqrt_n_multiplier)
    )
  )
  started <- proc.time()[["elapsed"]]
  nuisance_fit <- do.call(
    fit_binary_cross_fitted_direct_rkhs_nuisances,
    c(list(data = data, x_cols = x_cols, fold_id = fold_id), configuration)
  )
  nuisance_seconds <- proc.time()[["elapsed"]] - started
  save_atomic(
    list(
      contract = contract,
      nuisance_fit = nuisance_fit,
      nuisance_seconds = nuisance_seconds,
      completed_at = as.character(Sys.time())
    ),
    cache_file
  )
  message(sprintf("Saved nuisance checkpoint after %.1f seconds.", nuisance_seconds))
}

if (
  !identical(nuisance_fit$fold_id, fold_id) ||
    nuisance_fit$learner != "direct_rkhs_gamma" ||
    nrow(nuisance_fit$predictions) != nrow(data) ||
    any(!is.finite(as.matrix(nuisance_fit$predictions))) ||
    nrow(nuisance_fit$diagnostics) != 25L ||
    nrow(nuisance_fit$gamma_diagnostics) != 10L ||
    nrow(nuisance_fit$gamma_cv_table) !=
      10L * length(gamma_lambda_constants) *
        length(gamma_bandwidth_multipliers) ||
    nrow(nuisance_fit$gamma_cv_resample_table) !=
      10L * length(gamma_lambda_constants) *
        length(gamma_bandwidth_multipliers) * 3L * rkhs_cv_repeats
) {
  stop("The fitted nuisance object failed its application validation contract.")
}

endpoint_rows <- list()
confidence_rows <- list()
diagnostic_rows <- list()
summary_rows <- list()
for (Lambda in Lambda_grid) {
  for (epsilon in epsilon_grid) {
    value <- evaluate_binary_bounds_from_nuisances(
      data = data,
      nuisance_fit = nuisance_fit,
      Lambda0 = Lambda,
      Lambda1 = Lambda,
      epsilon = epsilon,
      confidence_level = confidence_level,
      compatibility_action = "project"
    )
    endpoint <- value$endpoint_intervals
    endpoint$Lambda <- Lambda
    endpoint$epsilon <- epsilon
    endpoint_rows[[length(endpoint_rows) + 1L]] <- endpoint

    confidence <- data.frame(
      Lambda = Lambda,
      epsilon = epsilon,
      target = if (epsilon == 0) "exact" else "smoothed",
      confidence_level = confidence_level,
      confidence_lower = value$ate_confidence_set$lower,
      confidence_upper = value$ate_confidence_set$upper,
      confidence_width =
        value$ate_confidence_set$upper - value$ate_confidence_set$lower,
      stringsAsFactors = FALSE
    )
    confidence_rows[[length(confidence_rows) + 1L]] <- confidence

    diagnostic <- value$diagnostics
    diagnostic$Lambda <- Lambda
    diagnostic$epsilon <- epsilon
    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- diagnostic

    estimates <- value$estimates[
      value$estimates$method == "dual_smooth_one_step", , drop = FALSE
    ]
    get_estimate <- function(endpoint_name) {
      estimates$estimate[estimates$endpoint == endpoint_name]
    }
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      Lambda = Lambda,
      epsilon = epsilon,
      target = if (epsilon == 0) "exact" else "smoothed",
      mu0_lower = get_estimate("mu0_lower"),
      mu0_upper = get_estimate("mu0_upper"),
      mu1_lower = get_estimate("mu1_lower"),
      mu1_upper = get_estimate("mu1_upper"),
      ate_lower = get_estimate("ate_lower"),
      ate_upper = get_estimate("ate_upper"),
      ate_width = get_estimate("ate_upper") - get_estimate("ate_lower"),
      outward_95_lower = confidence$confidence_lower,
      outward_95_upper = confidence$confidence_upper,
      outward_95_width = confidence$confidence_width,
      stringsAsFactors = FALSE
    )
    message(sprintf(
      paste0(
        "Lambda=%g, %s: ATE bounds [%.6f, %.6f], ",
        "95%% outward set [%.6f, %.6f]."
      ),
      Lambda, if (epsilon == 0) "exact" else "smoothed",
      get_estimate("ate_lower"), get_estimate("ate_upper"),
      confidence$confidence_lower, confidence$confidence_upper
    ))
  }
}

endpoint_detail <- bind_rows(endpoint_rows)
confidence_detail <- bind_rows(confidence_rows)
envelope_diagnostics <- bind_rows(diagnostic_rows)
bounds_summary <- bind_rows(summary_rows)
if (
  nrow(endpoint_detail) != 6L * length(Lambda_grid) * length(epsilon_grid) ||
    nrow(bounds_summary) != length(Lambda_grid) * length(epsilon_grid) ||
    any(bounds_summary$ate_lower > bounds_summary$ate_upper) ||
    any(bounds_summary$outward_95_lower > bounds_summary$outward_95_upper)
) {
  stop("The final application bound table failed validation.")
}

mean_weights <- stats::aggregate(
  coefficient ~ component + candidate,
  nuisance_fit$cv_table,
  mean
)
gamma_diagnostics <- nuisance_fit$gamma_diagnostics
gamma_diagnostics$selected_median_pair_kernel_similarity <- exp(
  -1 / (2 * gamma_diagnostics$selected_bandwidth_multiplier^2)
)
gamma_diagnostics$nystrom_effective_rank_fraction <-
  gamma_diagnostics$nystrom_effective_rank /
    gamma_diagnostics$nystrom_landmark_count
gamma_candidate_geometry <- unique(nuisance_fit$gamma_cv_table[, c(
  "outer_fold", "treatment", "bandwidth_multiplier",
  "nystrom_effective_rank_min", "nystrom_effective_rank_max"
)])
gamma_candidate_geometry$median_pair_kernel_similarity <- exp(
  -1 / (2 * gamma_candidate_geometry$bandwidth_multiplier^2)
)
gamma_candidate_geometry <- gamma_candidate_geometry[order(
  gamma_candidate_geometry$outer_fold,
  gamma_candidate_geometry$treatment,
  gamma_candidate_geometry$bandwidth_multiplier
), , drop = FALSE]
gamma_cv_repeat_selection <- stats::aggregate(
  risk ~ outer_fold + treatment + cv_repeat + lambda +
    bandwidth_multiplier + n_training_t + lambda_sqrt_n_constant,
  nuisance_fit$gamma_cv_resample_table,
  mean
)
gamma_cv_repeat_selection$selected_within_repeat <- FALSE
repeat_groups <- split(
  seq_len(nrow(gamma_cv_repeat_selection)),
  interaction(
    gamma_cv_repeat_selection$outer_fold,
    gamma_cv_repeat_selection$treatment,
    gamma_cv_repeat_selection$cv_repeat,
    drop = TRUE
  )
)
for (indices in repeat_groups) {
  gamma_cv_repeat_selection$selected_within_repeat[
    indices[which.min(gamma_cv_repeat_selection$risk[indices])]
  ] <- TRUE
}
gamma_cv_repeat_selection <- gamma_cv_repeat_selection[
  gamma_cv_repeat_selection$selected_within_repeat,
  , drop = FALSE
]
gamma_cv_repeat_selection <- gamma_cv_repeat_selection[order(
  gamma_cv_repeat_selection$outer_fold,
  gamma_cv_repeat_selection$treatment,
  gamma_cv_repeat_selection$cv_repeat
), , drop = FALSE]
rownames(gamma_cv_repeat_selection) <- NULL
if (
  any(gamma_diagnostics$selected_bandwidth_multiplier >
        max(gamma_bandwidth_multipliers) + 1e-12) ||
    any(gamma_diagnostics$selected_median_pair_kernel_similarity >
        exp(-1 / (2 * max(gamma_bandwidth_multipliers)^2)) + 1e-12)
) {
  stop("The selected gamma bandwidth failed its prespecified flatness cap.")
}
settings <- data.frame(
  version = version,
  specification = covariate_specification,
  observations = nrow(data),
  covariates = paste(x_cols, collapse = ","),
  Lambda = paste(Lambda_grid, collapse = ","),
  epsilon = paste(epsilon_grid, collapse = ","),
  confidence_level = confidence_level,
  folds = folds,
  seed = seed,
  gamma_method = "direct_rkhs_riesz",
  gamma_tuning = "penalty_cv_fixed_half_median_squared_bandwidth",
  gamma_lambda_rule = "C/sqrt(n_t)",
  gamma_lambda_constants = paste(gamma_lambda_constants, collapse = ","),
  gamma_bandwidth_rule = paste(
    "sigma^2 equals one half the median pairwise squared distance",
    "of standardized (Y,X) in each training fold"
  ),
  gamma_bandwidth_multipliers =
    paste(gamma_bandwidth_multipliers, collapse = ","),
  gamma_cv_folds = 3L,
  gamma_cv_repeats = rkhs_cv_repeats,
  gamma_backend = "nystrom",
  gamma_rank_rule = sprintf(
    "min(n_t,ceiling(%g*sqrt(n_t)))", nystrom_sqrt_n_multiplier
  ),
  nuisance_learner = "SuperLearner_convex_stack",
  nuisance_candidates = paste(candidates, collapse = ","),
  exclude_hal_sensitivity = exclude_hal,
  exclude_gam_sensitivity = exclude_gam,
  nuisance_cv_folds = 3L,
  compatibility_action = "project",
  nuisance_seconds = nuisance_seconds,
  stringsAsFactors = FALSE
)

write_table(settings, "settings.csv")
write_table(bounds_summary, "bounds_summary.csv")
write_table(endpoint_detail, "endpoint_detail.csv")
write_table(confidence_detail, "confidence_set_detail.csv")
write_table(envelope_diagnostics, "envelope_diagnostics.csv")
write_table(nuisance_fit$diagnostics, "nuisance_diagnostics.csv")
write_table(nuisance_fit$cv_table, "nuisance_cv.csv")
write_table(mean_weights, "mean_ensemble_weights.csv")
write_table(gamma_diagnostics, "gamma_diagnostics.csv")
write_table(gamma_candidate_geometry, "gamma_candidate_geometry.csv")
write_table(gamma_cv_repeat_selection, "gamma_cv_repeat_selection.csv")
write_table(nuisance_fit$gamma_cv_table, "gamma_cv.csv")
write_table(nuisance_fit$gamma_cv_resample_table, "gamma_cv_resamples.csv")
write_table(
  nuisance_fit$reconstruction_diagnostics,
  "reconstruction_diagnostics.csv"
)
save_atomic(
  list(
    settings = settings,
    contract = contract,
    session_info = sessionInfo(),
    bounds_summary = bounds_summary,
    endpoint_detail = endpoint_detail,
    confidence_set_detail = confidence_detail,
    envelope_diagnostics = envelope_diagnostics,
    nuisance_checkpoint = basename(cache_file)
  ),
  file.path(output_directory, "ohie_updated_fixed_median_results.rds")
)

print(bounds_summary, row.names = FALSE)
cat("\nResults saved in: ", output_directory, "\n", sep = "")
