# Runtime adapters for the two immutable compatibility-test parent bundles.
# Only n and the Lambda grid change. The parent new_session functions remain
# the estimator entry points; the wrappers below collect compact diagnostics.

n12_adapter_contract <- function(outcome) {
  outcome <- match.arg(outcome, c('continuous', 'binary'))
  list(
    outcome = outcome,
    frozen_parent = if (outcome == 'continuous') 'x5_compatibility_lpc_v1' else
      'binary_x5_compatibility_mac_v1',
    methods = 'rkhs',
    parent_method_names = if (outcome == 'continuous')
      c(rkhs = 'rkhs_tuned') else
      c(rkhs = 'rkhs'),
    folds = 5L, bootstrap = 399L, multiplier = 'Gaussian',
    allocation_df = 5L, positivity_epsilon = .05,
    primitive_learner = if (outcome == 'continuous') 'super_learner' else
      'binary_super_learner',
    primitives = c('glm_main', 'glm_interactions', 'glmnet_interactions',
                   'ranger', 'xgboost'),
    primitive_cv_folds = 3L, ranger_num_trees = 500L,
    rkhs = list(
      backend = 'nystrom', nystrom_rank = NULL,
      coordinate_standardization = 'Each training sample only',
      base_bandwidth = 'sqrt(median pairwise squared distance), frozen parent implementation',
      bandwidth_multipliers = if (outcome == 'continuous') c(.25, .5, 1, 2, 4) else
        1/sqrt(2),
      bandwidth_rule = if (outcome == 'continuous')
        'Joint bandwidth and penalty CV over the frozen parent grid' else
        'Fixed sigma^2 = M/2; no bandwidth CV',
      penalty = 'lambda = C/sqrt(n_t)', C = 10^seq(-1.5, 1, by = .5),
      cv_folds = 3L, cv_repeats = if (outcome == 'continuous') 1L else 2L,
      selection_rule = 'minimum'),
    normalization = 'Frozen positive prediction and baseline-stratum mean normalization',
    projection_applied = FALSE, q_fits = 0L,
    classification_gamma = FALSE, pmf_ratio_gamma = FALSE,
    seed_rule = if (outcome == 'continuous') '820000 + replication' else
      '261000 + source_setting + 100*(replication-1); null source_setting = 5',
    null_rule = if (outcome == 'continuous') 'a=1.1, b=0, stress=FALSE' else
      'Frozen generator with config$settings row 5 (a=1.1, b=0)')
}

n12_adapter_fingerprint <- function(x) {
  if (!requireNamespace('digest', quietly = TRUE)) stop('Package digest is required.')
  digest::digest(x, algo = 'sha256', serialize = TRUE)
}

n12_adapter_in_directory <- function(path, thunk) {
  previous <- getwd()
  on.exit(setwd(previous), add = TRUE)
  setwd(path)
  thunk()
}

# SuperLearner predicts via stats::predict from its namespace, which cannot
# discover unregistered methods in a private source environment. Register the
# exact frozen methods for this call and restore every original registry entry
# on exit. No learner definition is exported into the user's global workspace.
n12_adapter_with_methods <- function(parent, thunk) {
  methods <- ls(parent, pattern = '^predict\\.', all.names = TRUE)
  registry <- get('.__S3MethodsTable__.', envir = asNamespace('stats'), inherits = FALSE)
  present <- vapply(methods, function(name)
    exists(name, envir = registry, inherits = FALSE), logical(1))
  previous <- lapply(methods[present], get, envir = registry, inherits = FALSE)
  names(previous) <- methods[present]
  on.exit({
    for (name in methods) {
      if (name %in% names(previous)) assign(name, previous[[name]], envir = registry) else
        if (exists(name, envir = registry, inherits = FALSE)) rm(list = name, envir = registry)
    }
  }, add = TRUE)
  for (name in methods) base::registerS3method('predict',
    sub('^predict\\.', '', name), get(name, envir = parent, inherits = FALSE),
    envir = asNamespace('stats'))
  thunk()
}

# Parent source() calls default to .GlobalEnv, and their new_session() calls
# explicitly use globalenv(). Redirect both names in this private lexical
# environment. Its parent is the package search path, excluding .GlobalEnv,
# so a current-project function cannot silently satisfy a missing dependency.
# .GlobalEnv itself is not masked: frozen RNG preservation must continue to
# address the real R random-number state.
n12_adapter_load_parent <- function(parent_root) {
  parent_root <- normalizePath(parent_root, mustWork = TRUE)
  private <- new.env(parent = parent.env(globalenv()))
  private$globalenv <- function() private
  private$source <- function(file, local = FALSE, ...) {
    target <- if (is.environment(local)) local else if (isTRUE(local))
      parent.frame() else private
    base::source(file, local = target, ...)
  }
  # Several frozen dependency guards use vapply(names, exists, ...). R's
  # default exists() lookup then starts in vapply's base frame, losing this
  # private lexical scope. Keep those guards scoped to their actual caller;
  # explicit environment lookups (including RNG and cache access) are intact.
  private$exists <- function(x, where = -1,
      envir = if (missing(frame)) as.environment(where) else sys.frame(frame),
      frame, mode = 'any', inherits = TRUE) {
    if (missing(envir) && missing(where) && missing(frame)) {
      lookup <- if (identical(sys.function(sys.parent()), base::vapply))
        parent.frame(2L) else parent.frame()
      return(base::exists(x, envir = lookup, mode = mode, inherits = inherits))
    }
    base::exists(x, envir = envir, mode = mode, inherits = inherits)
  }
  n12_adapter_in_directory(parent_root, function() {
    base::sys.source('scripts/common.R', envir = private)
  })
  private
}

n12_adapter_new <- function(bundle_root, outcome, Lambda = c(1, 1.1), n = 12000L) {
  outcome <- match.arg(outcome, c('continuous', 'binary'))
  stopifnot(is.numeric(Lambda), length(Lambda) > 0L,
            all(is.finite(Lambda)), all(Lambda >= 1),
            !anyDuplicated(Lambda), identical(sort(Lambda), Lambda),
            length(n) == 1L, is.finite(n), n > 0, n == as.integer(n))
  contract <- n12_adapter_contract(outcome)
  parent_root <- normalizePath(file.path(bundle_root, 'frozen', outcome), mustWork = TRUE)
  parent <- n12_adapter_load_parent(parent_root)
  stopifnot(parent$config$folds == contract$folds,
            parent$config$bootstrap == contract$bootstrap,
            parent$config$allocation_df == contract$allocation_df,
            identical(parent$config$primitives, contract$primitives))
  # Preserve every other frozen parent option, including outcome-specific CV.
  parent$config$n <- as.integer(n)
  parent$config$Lambda <- as.numeric(Lambda)
  session <- n12_adapter_in_directory(parent_root, function() parent$new_session())
  engine <- session$engine
  state <- new.env(parent = emptyenv())
  state$probability_hashes <- character()
  state$fold_id <- NULL
  state$previous <- NULL
  state$data_hash <- NULL

  # The parent's probability cache and fitted values are untouched. Hash the
  # values returned by that cache, including both training and held-out folds.
  probability_original <- engine$learned_witness_probability_fit
  engine$learned_witness_probability_fit <- function(...) {
    args <- list(...)
    fit <- do.call(probability_original, args)
    key <- as.character(args$library_seed)
    hash <- n12_adapter_fingerprint(fit[c('training', 'validation')])
    if (key %in% names(state$probability_hashes))
      stopifnot(identical(unname(state$probability_hashes[[key]]), hash))
    state$probability_hashes[key] <- hash
    fit
  }
  test_original <- engine$test_learned_witness_dual_msm_compatibility
  engine$test_learned_witness_dual_msm_compatibility <- function(...) {
    fit <- test_original(...)
    state$fold_id <- fit$fold_id
    fit
  }
  forbidden <- function(...) stop('Forbidden compatibility projection, q fitting, or alternative gamma learner.')
  forbidden_names <- c('project_gamma_to_compatible_interval',
    'fit_learned_witness_binary_gamma', 'fit_learned_witness_discrete_gamma',
    'fit_threshold_by_check_loss', 'fit_threshold_spline_cv',
    'fit_threshold_rkhs_cv', 'fit_threshold_model')
  for (name in forbidden_names) {
    assign(name, forbidden, envir = engine)
    assign(name, forbidden, envir = parent)
  }

  generate <- function(task, n_used = n) {
    stopifnot(is.data.frame(task), nrow(task) == 1L,
      all(c('outcome', 'setting', 'type', 'a', 'b', 'replication', 'seed') %in% names(task)),
      identical(as.character(task$outcome), outcome),
      length(n_used) == 1L, is.finite(n_used), n_used > 0,
      n_used == as.integer(n_used), task$type %in% c('stress', 'null'))
    source_setting <- if (task$type == 'null') 5L else as.integer(task$setting)
    stopifnot(source_setting %in% 1:5,
      if (task$type == 'null') task$setting == 0L else task$setting %in% 1:4)
    expected_a <- c(.55, 1.1, 1.65, 1.1, 1.1)[source_setting]
    expected_b <- c(.075, .075, .075, .15, 0)[source_setting]
    expected_seed <- if (outcome == 'continuous') 820000L + task$replication else
      261000L + source_setting + 100L * (task$replication - 1L)
    stopifnot(task$a == expected_a, task$b == expected_b, task$seed == expected_seed)
    n12_adapter_in_directory(parent_root, function() {
      if (outcome == 'continuous') {
        parent$x5g_generate(as.integer(n_used), task$a, task$b, task$seed,
                            stress = task$type == 'stress')
      } else {
        original_task <- task
        original_task$setting <- source_setting
        parent$generate(original_task, n = as.integer(n_used))
      }
    })
  }

  run <- function(data, method, seed) {
    stopifnot(length(method) == 1L, method %in% contract$methods,
              is.data.frame(data), nrow(data) > 0L)
    data_hash <- n12_adapter_fingerprint(data)
    if (!is.null(state$data_hash) && !identical(state$data_hash, data_hash))
      stop('An adapter session is paired within one task; create a new adapter for another data set.')
    state$data_hash <- data_hash
    parent_method <- unname(contract$parent_method_names[[method]])
    fit <- n12_adapter_in_directory(parent_root, function()
      n12_adapter_with_methods(parent, function() {
        if (outcome == 'continuous') session$run(data, parent_method, seed) else
          session$run(data, parent_method, seed, progress = TRUE)
      }))
    hashes <- state$probability_hashes[sort(names(state$probability_hashes))]
    paired <- list(data = data_hash, folds = n12_adapter_fingerprint(state$fold_id),
                   primitives = hashes)
    stopifnot(length(state$fold_id) == nrow(data), length(hashes) == contract$folds,
      setequal(names(hashes), as.character(seed + 10000L * seq_len(contract$folds))))
    if (!is.null(state$previous)) stopifnot(identical(state$previous, paired))
    if (!is.null(fit$paired)) stopifnot(identical(fit$paired, paired))
    state$previous <- paired

    if (outcome == 'continuous') {
      # Parent cache keys encode gamma seed and treatment. Add labels and C
      # without changing or reordering any retained parent CV table.
      selections <- fit$gamma_selection
      for (key in names(selections)) {
        coordinates <- as.integer(strsplit(key, '_', fixed = TRUE)[[1L]])
        stopifnot(length(coordinates) == 2L)
        gamma_seed <- coordinates[1L]
        treatment <- coordinates[2L]
        fold <- (gamma_seed - seed - 1000L * treatment)/10000L
        stopifnot(treatment %in% 0:1, fold %in% seq_len(contract$folds))
        selections[[key]]$treatment <- treatment
        selections[[key]]$seed <- gamma_seed
        selections[[key]]$outer_fold <- as.integer(fold)
        selections[[key]]$n_arm <- sum(data$T_obs == treatment & state$fold_id != fold)
        selections[[key]]$selected_C <- selections[[key]]$selected_lambda *
          sqrt(selections[[key]]$n_arm)
      }
      fit$gamma_selection <- selections
      fit$witness_candidate_counts <- as.data.frame(table(fit$witness_diagnostics$family))
    }
    stopifnot(length(fit$gamma_selection) == 2L * contract$folds)
    fit$settings$gamma_learner <- method
    fit$settings$frozen_parent <- contract$frozen_parent
    fit$settings$gamma_tuning_contract <- contract[[method]]
    fit$gamma_diagnostics$learner <- method
    list(summary = fit$summary, gamma_selection = fit$gamma_selection,
      gamma_diagnostics = fit$gamma_diagnostics,
      probability_diagnostics = fit$probability_diagnostics, settings = fit$settings,
      paired = paired, fold_sizes = fit$fold_sizes,
      projection_applied = FALSE, q_fits = 0L,
      witness_candidate_counts = fit$witness_candidate_counts)
  }
  list(generate = generate, run = run, contract = contract)
}
