# Frozen binary X5 stress design; honest learned-witness compatibility tests.
source('R/binary_outcome_bounds.R')
source('R/sharp_bounds_plugin.R')
source('R/x5_adaptive_gamma_provider.R')
source('R/conservative_compatibility_test.R')
config <- readRDS('study_config.rds')
stopifnot(config$n == 6000L, config$reps == 500L,
          identical(config$methods, c('rkhs', 'growing_direct')),
          identical(config$Lambda, c(1, 1.1, 1.3, 1.5, 2)))
tasks <- data.frame(id = seq_len(nrow(config$settings) * config$reps),
  setting = rep(seq_len(nrow(config$settings)), each = config$reps),
  replication = rep(seq_len(config$reps), nrow(config$settings)))
tasks$seed <- 261000L + tasks$setting + 100L * (tasks$replication - 1L)
tasks$n <- config$n
gate_ids <- (seq_len(nrow(config$settings)) - 1L) * config$reps + 1L

# Evaluate only the archived generator's function definition. Do not execute
# its old runner or inherit its endpoint-estimation stages.
extract_generator <- function(path) {
  expressions <- as.list(parse(path))
  hit <- vapply(expressions, function(x) is.call(x) &&
    identical(x[[1L]], as.name('<-')) && identical(x[[2L]], as.name('generate')), logical(1))
  stopifnot(sum(hit) == 1L)
  expressions[[which(hit)]][[3L]]
}
generate <- eval(extract_generator('provenance/binary_stress_common.R'))
fingerprint <- function(x) digest::digest(x, algo = 'sha256', serialize = TRUE)
atomic <- function(object, path, overwrite = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!overwrite && file.exists(path)) stop('Refusing to overwrite: ', path)
  tmp <- tempfile('unpublished_', dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp))
  withCallingHandlers(saveRDS(object, tmp, compress = 'gzip'),
    warning = function(w) stop(conditionMessage(w)))
  stopifnot(file.info(tmp)$size > 0, identical(readRDS(tmp), object))
  if (!overwrite && file.exists(path)) stop('Destination appeared: ', path)
  if (!file.rename(tmp, path)) stop('Atomic publication failed: ', path)
  stopifnot(identical(readRDS(path), object))
  invisible(path)
}
binary_truth <- function(d, setting) {
  g <- config$settings[setting, ]
  x <- as.matrix(d[paste0('X', 1:5)])
  h <- drop(x %*% c(.14, .12, .10, .08, .06))
  gy <- drop(x %*% c(.20, -.25, .30, .15, -.20))
  shape <- 1 + .6*sin(x[,1]*x[,2]) + .4*sin(x[,3]*x[,4])
  list(e = rep(.5, nrow(d)), pi0 = plogis(h-g$a/2), pi1 = plogis(h+g$a/2),
    outcome = lapply(0:1, function(t) list(p0 = plogis(-.4+1.5*t+gy),
      p1 = plogis(-.4+1.5*t+gy+(2*t-1)*g$b*shape))))
}
population_labels <- function() {
  audit <- readRDS('provenance/population_frontier_audit.rds')
  g <- expand.grid(setting = seq_len(nrow(config$settings)), Lambda0 = config$Lambda)
  g$a <- config$settings$a[g$setting]; g$b <- config$settings$b[g$setting]
  g$numerical_frontier <- audit$summary$frontier[match(g$setting, audit$summary$setting)]
  g$label <- ifelse(g$Lambda0 < g$numerical_frontier,
                    'incompatible: numerical positive witness', 'numerically compatible, not globally certified')
  g$numerical_frontier[g$b == 0] <- 1
  g$label[g$b == 0] <- 'analytically compatible b=0 control'
  g
}

new_session <- function() {
  engine <- new.env(parent = globalenv())
  sys.source('R/learned_witness_compatibility_test.R', engine)
  probability_original <- engine$learned_witness_probability_fit
  cache <- new.env(parent = emptyenv())
  engine$learned_witness_probability_fit <- function(...) {
    args <- list(...)
    key <- as.character(args$library_seed)
    hash <- fingerprint(args)
    if (exists(key, cache, inherits = FALSE)) {
      z <- get(key, cache); stopifnot(identical(z$hash, hash)); return(z$value)
    }
    z <- do.call(probability_original, args)
    # Cache predictions, not fitted forests. Each result is paired across gamma
    # learners; predictions are never written into checkpoint files.
    diagnostics <- data.frame(component = c('e','pi0','pi1'),
      training_clip_rate = vapply(c('e','pi0','pi1'), function(k)
        mean(z$training[[k]] <= config$positivity_epsilon |
             z$training[[k]] >= 1-config$positivity_epsilon), numeric(1)))
    z$fit <- list(diagnostics = diagnostics)
    assign(key, list(hash = hash, value = z), cache)
    z
  }
  forbidden <- function(...) stop('Forbidden: compatibility projection or q fitting in formal testing.')
  engine$project_gamma_to_compatible_interval <- forbidden
  engine$fit_learned_witness_binary_gamma <- function(...) stop('Classification gamma is forbidden.')
  engine$fit_learned_witness_discrete_gamma <- function(...) stop('PMF-ratio gamma is not this experiment.')
  previous <- NULL
  run <- function(data, method, seed, progress = TRUE) {
    stopifnot(method %in% config$methods)
    selection <- list()
    original_gamma <- fit_gamma_rkhs
    engine$fit_gamma_rkhs <- function(train, nuisance_train, t_value, x_cols, seed, ...) {
      fit <- if (method == 'rkhs') {
        original_gamma(train = train, nuisance_train = nuisance_train,
          t_value = t_value, x_cols = x_cols, seed = seed,
          lambda_grid = default_gamma_lambda_grid,
          bandwidth_multipliers = 1/sqrt(2), cv_folds = 3L, cv_repeats = 2L,
          selection_rule = 'minimum', backend = 'nystrom', nystrom_rank = NULL)
      } else x5a_fit(train, nuisance_train, t_value, seed, log = FALSE)
      record <- fit[intersect(names(fit), c('selected_lambda','selected_level',
        'feature_count','base_bandwidth','selected_bandwidth','selected_bandwidth_multiplier',
        'cv_table','cv_resample_table','rkhs_backend','normalization_constant',
        'training_negative_rate','clip_rate','gradient'))]
      record$treatment <- t_value; record$seed <- seed
      record$n_arm <- sum(train$T_obs == t_value)
      record$selected_C <- fit$selected_lambda * sqrt(record$n_arm)
      selection[[length(selection)+1L]] <<- record
      fit
    }
    fit <- engine$test_learned_witness_dual_msm_compatibility(
      data = data, lambda_grid = config$Lambda, x_cols = paste0('X',1:5),
      folds = config$folds, primitive_learner = 'binary_super_learner',
      gamma_learner = 'rkhs_cv', positivity_epsilon = config$positivity_epsilon,
      super_learner_candidates = config$primitives,
      super_learner_cv_folds = 3L, super_learner_ranger_num_trees = 500L,
      gamma_cv_lambda_grid = default_gamma_lambda_grid,
      gamma_cv_bandwidth_multipliers = 1/sqrt(2), gamma_cv_folds = 3L,
      gamma_cv_repeats = 2L, gamma_cv_selection_rule = 'minimum',
      gamma_nystrom_rank = NULL, allocation_spline_df = config$allocation_df,
      bootstrap_replications = config$bootstrap, seed = seed, progress = progress)
    hashes <- vapply(sort(ls(cache)), function(key) {
      z <- get(key, cache)$value
      fingerprint(z[c('training','validation')])
    }, character(1))
    paired <- list(data = fingerprint(data), folds = fingerprint(fit$fold_id), primitives = hashes)
    if (!is.null(previous)) stopifnot(identical(previous, paired))
    previous <<- paired
    stopifnot(length(selection) == 2L*config$folds, length(hashes) == config$folds)
    fit$settings$gamma_learner <- method
    fit$gamma_diagnostics$learner <- method
    fit$settings$gamma_normalization <- 'Positive prediction, training baseline-stratum mean normalization only; no compatibility projection or binary PMF reconstruction'
    # Retain fold-level selections and small diagnostic tables. Never save
    # observation-level data, predictions, multiplier draws, or model objects.
    list(summary = fit$summary, gamma_diagnostics = fit$gamma_diagnostics,
      gamma_selection = selection, probability_diagnostics = fit$probability_diagnostics,
      settings = fit$settings, paired = paired,
      fold_sizes = tabulate(fit$fold_id, config$folds),
      witness_candidate_counts = as.data.frame(table(fit$witness_diagnostics$family)),
      projection_applied = FALSE, q_fits = 0L)
  }
  list(run = run, engine = engine)
}
check_record <- function(v, task, method, signature, n_used) {
  stopifnot(identical(v$task, task), identical(v$method, method),
    identical(v$signature_hash, fingerprint(signature)), v$n_used == n_used,
    is.logical(v$success), length(v$success) == 1L, !is.na(v$success))
  if (v$success) {
    s <- v$summary
    valid_hash <- function(x) is.character(x) && length(x) == 1L &&
      !is.na(x) && grepl('^[0-9a-f]{64}$', x)
    p <- v$paired
    stopifnot(is.list(p), valid_hash(p$data), valid_hash(p$folds),
      is.character(p$primitives), length(p$primitives) == config$folds,
      !is.null(names(p$primitives)), !anyDuplicated(names(p$primitives)),
      setequal(names(p$primitives), as.character(task$seed + 10000L*seq_len(config$folds))),
      all(vapply(p$primitives, valid_hash, logical(1))))
    stopifnot(nrow(s) == length(config$Lambda),
      identical(as.numeric(s$Lambda0), config$Lambda), all(s$Lambda1 == s$Lambda0),
      all(is.finite(s$p_value)), all(s$p_value >= 0 & s$p_value <= 1),
      is.logical(s$reject), !anyNA(s$reject), isFALSE(v$projection_applied),
      v$q_fits == 0L, length(v$gamma_selection) == 2*config$folds,
      sum(v$fold_sizes) == n_used)
  } else stopifnot(length(v$error) == 1L, nzchar(v$error))
}
load_task <- function(path, task, signature, n_used = config$n) {
  if (!file.exists(path)) return(list(task = task, signature_hash = fingerprint(signature),
                                     n_used = n_used, records = list()))
  z <- readRDS(path)
  stopifnot(identical(z$task, task), identical(z$signature_hash, fingerprint(signature)),
    z$n_used == n_used, is.list(z$records), !anyDuplicated(names(z$records)),
    all(names(z$records) %in% config$methods))
  for (m in names(z$records)) check_record(z$records[[m]], task, m, signature, n_used)
  ok <- vapply(z$records, function(v) isTRUE(v$success), logical(1))
  if (sum(ok) == 2L) stopifnot(identical(z$records[[1]]$paired, z$records[[2]]$paired))
  z
}
run_task <- function(id, signature, smoke = FALSE, session_factory = new_session) {
  stopifnot(length(id) == 1L, id %in% tasks$id)
  task <- tasks[id, , drop = FALSE]
  n_used <- if (smoke) 600L else config$n
  root <- file.path(config$result_root, if (smoke) 'smoke' else 'runs')
  path <- file.path(root, sprintf('task_%04d.rds', id))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  lock <- paste0(path, '.lock')
  if (!dir.create(lock, showWarnings = FALSE)) stop('Task lock exists; inspect before resuming: ', lock)
  on.exit(unlink(lock, recursive = TRUE))
  z <- load_task(path, task, signature, n_used)
  if (length(z$records) == length(config$methods)) {
    cat('ALREADY COMPLETE', id, '; checkpoint unchanged\n'); return(invisible(z))
  }
  data <- generate(task, n = n_used)
  session <- session_factory()
  for (method in config$methods) {
    if (method %in% names(z$records)) next
    warnings <- character(); started <- proc.time()[3]
    cat(format(Sys.time()), 'START task', id, 'setting', task$setting,
        'rep', task$replication, method, 'n', n_used, '\n'); flush.console()
    v <- tryCatch(withCallingHandlers({
      fit <- session$run(data, method, task$seed)
      c(list(success = TRUE, error = ''), fit)
    }, warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')}),
      error = function(e) list(success = FALSE, error = conditionMessage(e)))
    v$task <- task; v$method <- method; v$signature_hash <- fingerprint(signature)
    v$n_used <- n_used; v$seconds <- unname(proc.time()[3]-started); v$warnings <- unique(warnings)
    # Schema and pairing errors are retained like estimator errors. A resume
    # must never silently refit a learner because its output failed validation.
    v <- tryCatch({
      check_record(v, task, method, signature, n_used)
      if (isTRUE(v$success)) for (prior in z$records) {
        if (isTRUE(prior$success)) stopifnot(identical(prior$paired, v$paired))
      }
      v
    }, error = function(e) {
      v$success <- FALSE
      v$error <- paste('Output validation failed:', conditionMessage(e))
      v
    })
    check_record(v, task, method, signature, n_used)
    z$records[[method]] <- v
    atomic(z, path)
    cat(format(Sys.time()), 'task', id, 'setting', task$setting, 'rep', task$replication,
        method, 'success', v$success, 'seconds', v$seconds, v$error, '\n'); flush.console()
  }
  invisible(z)
}
