#!/usr/bin/env Rscript
# Run from a built bundle root. All estimator calls are mocked: these tests
# exercise generators, frozen calibration, schema checks, and persistence only.
source('scripts/common.R')
source('scripts/adapters.R')

n12_test_expect_error <- function(expr, pattern = NULL) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, 'error'))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(error)))
  invisible(error)
}
n12_test_extract <- function(path, name, environment) {
  expressions <- as.list(parse(path))
  keep <- vapply(expressions, function(x) is.call(x) &&
    identical(x[[1L]], as.name('<-')) && identical(x[[2L]], as.name(name)),
    logical(1))
  stopifnot(sum(keep) == 1L)
  value <- eval(expressions[[which(keep)]][[3L]], envir = environment)
  assign(name, value, envir = environment)
  invisible(value)
}

config <- n12_config()
tasks <- n12_tasks(config)
stopifnot(config$n == 12000L, config$reps == 500L,
  identical(config$Lambda, c(1, 1.1)),
  identical(config$methods, c('rkhs', 'growing_direct')),
  identical(tasks$id, seq_len(5000L)), nrow(tasks) == 5000L,
  all(tasks$n == 12000L),
  identical(as.character(tasks$outcome), rep(c('continuous', 'binary'), each = 2500L)),
  identical(as.integer(tasks$setting), rep(rep(c(1:4, 0L), each = 500L), 2L)),
  identical(as.integer(tasks$replication), rep(seq_len(500L), 10L)),
  all(table(tasks$outcome, tasks$setting) == 500L),
  !anyDuplicated(tasks[c('outcome', 'setting', 'replication')]))
for (outcome in c('continuous', 'binary')) {
  local <- tasks[tasks$outcome == outcome, , drop = FALSE]
  source_setting <- ifelse(local$setting == 0L, 5L, local$setting)
  stopifnot(identical(as.numeric(local$a), rep(c(.55, 1.1, 1.65, 1.1, 1.1), each = 500L)),
    identical(as.numeric(local$b), rep(c(.075, .075, .075, .15, 0), each = 500L)),
    all(local$type == ifelse(local$setting == 0L, 'null', 'stress')))
  if (outcome == 'continuous') stopifnot(
    all(local$seed == 820000L + local$replication),
    identical(as.integer(local$parent_task_id), seq.int(2001L, 4500L))) else stopifnot(
    all(local$seed == 261000L + source_setting + 100L * (local$replication - 1L)),
    !anyDuplicated(local$seed), identical(as.integer(local$parent_task_id), seq_len(2500L)))
  contract <- n12_adapter_contract(outcome)
  stopifnot(contract$folds == 5L, contract$bootstrap == 399L,
    contract$allocation_df == 5L, contract$positivity_epsilon == .05,
    identical(contract$primitives, c('glm_main', 'glm_interactions',
      'glmnet_interactions', 'ranger', 'xgboost')),
    identical(contract$methods, config$methods),
    isFALSE(contract$projection_applied), contract$q_fits == 0L,
    isFALSE(contract$classification_gamma), isFALSE(contract$pmf_ratio_gamma),
    contract$growing_direct$cv_repeats == 2L)
  if (outcome == 'continuous') stopifnot(
    identical(contract$rkhs$bandwidth_multipliers, c(.25, .5, 1, 2, 4)),
    contract$rkhs$cv_repeats == 1L) else stopifnot(
    identical(contract$rkhs$bandwidth_multipliers, 1/sqrt(2)),
    contract$rkhs$cv_repeats == 2L)
}

# Evaluate original parent function definitions independently of the adapter.
# This is a generator parity check on small samples, not a fitted smoke run.
parent_generators <- list()
continuous <- new.env(parent = globalenv())
continuous_root <- file.path('frozen', 'continuous', 'R')
for (name in c('siv_normalizer', 'siv_gamma'))
  n12_test_extract(file.path(continuous_root, 'scalar_iv_population_experiment.R'), name, continuous)
for (name in c('x5s_contract', 'x5s_oracle', 'x5s_generate'))
  n12_test_extract(file.path(continuous_root, 'x5_simple_study.R'), name, continuous)
n12_test_extract(file.path(continuous_root, 'x5_growing_gamma.R'), 'x5g_generate', continuous)
binary <- new.env(parent = globalenv())
binary$config <- readRDS('frozen/binary/study_config.rds')
n12_test_extract('frozen/binary/provenance/binary_stress_common.R', 'generate', binary)
for (outcome in c('continuous', 'binary')) {
  adapter <- n12_adapter_new(getwd(), outcome)
  ids <- tasks$id[tasks$outcome == outcome & tasks$replication == 1L]
  for (id in ids) {
    task <- tasks[id, , drop = FALSE]
    observed <- adapter$generate(task, n_used = 120L)
    if (outcome == 'continuous') expected <- continuous$x5g_generate(
      120L, task$a, task$b, task$seed, stress = task$type == 'stress') else {
      parent_task <- task
      parent_task$setting <- if (task$setting == 0L) 5L else task$setting
      expected <- binary$generate(parent_task, 120L)
    }
    stopifnot(identical(observed, expected), nrow(observed) == 120L,
      all(abs(as.matrix(observed[paste0('X', 1:5)])) <= 2),
      all(observed$T_obs %in% 0:1), all(observed$Z %in% 0:1))
    if (outcome == 'binary') stopifnot(all(observed$Y %in% 0:1)) else {
      residual <- observed$Y - drop(as.matrix(observed[paste0('X', 1:5)]) %*%
        c(.20, -.25, .30, .15, -.20)) - .75 * observed$T_obs
      stopifnot(all(abs(residual) <= 2 + 1e-12))
      if (task$setting == 0L) stopifnot(all(observed$gamma0_oracle == 1),
        all(observed$gamma1_oracle == 1))
    }
  }
}

# The parent pointwise multiplier law is invariant to discarding larger
# Lambdas: draws are G-dimensional, and calibration resets its own seed.
# Grid-simultaneous p-values deliberately refer to the active grid instead.
for (outcome in c('continuous', 'binary')) {
  environment <- new.env(parent = globalenv())
  calibrate <- n12_test_extract(file.path('frozen', outcome, 'R',
    'learned_witness_compatibility_test.R'), 'calibrate_learned_witness_maximum', environment)
  grid <- if (outcome == 'continuous') c(1, 1.1, 1.25, 1.5, 1.8, 2) else
    c(1, 1.1, 1.3, 1.5, 2)
  metadata <- data.frame(Lambda0 = rep(grid, each = 8L),
    Lambda1 = rep(grid, each = 8L), treatment = rep(rep(0:1, each = 4L), length(grid)),
    family = rep(c('pointwise_lower', 'pointwise_upper',
      'normalization_lower', 'normalization_upper'), 2L * length(grid)))
  signals <- outer(seq_len(80L), seq_len(nrow(metadata)),
    function(i, j) sin(i * (.03 + j/90)) + cos(i/7 + j/4) + j/200)
  keep <- metadata$Lambda0 %in% c(1, 1.1)
  set.seed(991L); rng <- .Random.seed
  full <- calibrate(signals, metadata, bootstrap_replications = 39L, seed = 9912L)
  stopifnot(identical(.Random.seed, rng))
  subset <- calibrate(signals[, keep, drop = FALSE], metadata[keep, , drop = FALSE],
    bootstrap_replications = 39L, seed = 9912L)
  columns <- setdiff(names(full$summary), grep('grid_simultaneous', names(full$summary), value = TRUE))
  expected <- full$summary[full$summary$Lambda0 %in% c(1, 1.1), columns, drop = FALSE]
  obtained <- subset$summary[, columns, drop = FALSE]
  rownames(expected) <- rownames(obtained) <- NULL
  stopifnot(isTRUE(all.equal(expected, obtained, tolerance = 1e-13)),
    isTRUE(all.equal(unname(full$bootstrap[, colnames(subset$bootstrap), drop = FALSE]),
      unname(subset$bootstrap), tolerance = 1e-13)), identical(.Random.seed, rng))
}

# Account-quota tests use the actual LPC Fileset layout and no external query.
quota_env <- new.env(parent=globalenv())
sys.source('scripts/storage.R',quota_env)
quota_env$file.exists <- function(...)TRUE
quota_env$file.access <- function(...)0L
quota_env$quota_row <- 'mmfs1 sas1-homeFS USR 67.15G 100G 104.6G 1.426G none 614966 0 0 335 none arcamgr01.pixstor'
quota_env$system2 <- function(...)quota_env$quota_row
invisible(capture.output(q <- quota_env$n12_check_storage()))
stopifnot(abs(q$available_bytes/1024^3-31.424)<1e-8)
legacy <- quota_env$n12_parse_quota('mmfs1 USR 67.15G 100G 104.6G 1.426G none')
stopifnot(identical(q$available_bytes,legacy$available_bytes))
quota_env$quota_row <- 'mmfs1 sas1-homeFS USR 97G 100G 104.6G 2G none 614966 0 0 335 none arcamgr01.pixstor'
invisible(capture.output(n12_test_expect_error(quota_env$n12_check_storage(),'Insufficient GPFS')))
n12_test_expect_error(quota_env$n12_parse_quota('mmfs1 a b USR 1G 100G 100G 0 none'))

test_directory <- tempfile('x5_n12000_unit_')
dir.create(test_directory)
signature <- list(config = config, manifest = c('unit-fixture' = 'test-only'),
  runtime = list(R = as.character(getRversion()), purpose = 'mock persistence tests'))
result_root <- file.path(test_directory, 'results')
checkpoint_path <- function(id) file.path(result_root, 'runs', sprintf('task_%04d.rds', id))
calls <- new.env(parent = emptyenv())
calls$fit <- character(); calls$generate <- 0L
mock_factory <- function(bundle_root, outcome, Lambda, n) list(
  generate = function(task, n_used = n) {
    calls$generate <- calls$generate + 1L
    data.frame(row = seq_len(n_used), seed = rep(task$seed, n_used))
  },
  run = function(data, method, seed) {
    calls$fit <- c(calls$fit, method)
    list(summary = data.frame(Lambda0 = Lambda, Lambda1 = Lambda,
      statistic = c(.1, .2), critical_value = c(2, 2),
      p_value = c(.5, .4), reject = c(FALSE, FALSE), moments_tested = c(8L, 8L)),
      gamma_selection = rep(list(list(selected_C = 1, selected_lambda = .02,
        cv_table = data.frame(C = 1, risk = 0),
        cv_resample_table = data.frame(repeat_id = 1L, fold = 1L, risk = 0))), 10L),
      paired = list(data = n12_fingerprint(data),
        folds = n12_fingerprint(rep(seq_len(5L), length.out = nrow(data))),
        primitives = setNames(rep(n12_fingerprint('mock primitive predictions'), 5L),
          as.character(seed + 10000L * seq_len(5L)))),
      fold_sizes = rep(nrow(data)/5L, 5L), projection_applied = FALSE, q_fits = 0L)
  })
quiet_run <- function(id, factory = mock_factory) {
  invisible(capture.output(value <- n12_run_task(id, signature,
    bundle_root = getwd(), result_root = result_root, adapter_factory = factory)))
  value
}
never_fit <- function(...) stop('Completed checkpoints must never construct an adapter.')

# Atomic publication verifies its temporary file before replacing a checkpoint.
atomic_path <- file.path(test_directory, 'atomic.rds')
atomic_object <- list(values = 1:5, diagnostic = data.frame(p = .5))
n12_atomic(atomic_object, atomic_path, overwrite = FALSE)
stopifnot(identical(readRDS(atomic_path), atomic_object))
atomic_hash <- tools::md5sum(atomic_path)
n12_test_expect_error(n12_atomic(list(replacement = TRUE), atomic_path, overwrite = FALSE))
stopifnot(identical(tools::md5sum(atomic_path), atomic_hash))
# Environments cannot pass an identical() serialization round trip; rejection
# must leave the previously published file byte-for-byte intact.
n12_test_expect_error(n12_atomic(new.env(parent = emptyenv()), atomic_path))
stopifnot(identical(tools::md5sum(atomic_path), atomic_hash))
n12_test_expect_error(n12_atomic(atomic_object, atomic_path,
  writer = function(object, path) {
    writeBin(charToRaw('interrupted temporary file'), path)
    stop('Intentional interrupted write')
  }))
stopifnot(identical(tools::md5sum(atomic_path), atomic_hash),
  !length(list.files(test_directory, pattern = '^unpublished_')))

first <- quiet_run(1L)
stopifnot(identical(calls$fit, config$methods), calls$generate == 1L,
  all(vapply(first$records, function(v) isTRUE(v$success), logical(1))),
  file.info(checkpoint_path(1L))$size < 50000L)
hash <- tools::md5sum(checkpoint_path(1L))
again <- quiet_run(1L, never_fit)
stopifnot(identical(first, again), identical(hash, tools::md5sum(checkpoint_path(1L))))
partial <- first
partial$records <- partial$records['rkhs']
n12_atomic(partial, checkpoint_path(1L))
calls$fit <- character()
resumed <- quiet_run(1L)
stopifnot(identical(calls$fit, 'growing_direct'),
  identical(resumed$records$rkhs, first$records$rkhs),
  identical(resumed$records$rkhs$paired, resumed$records$growing_direct$paired))
hash <- tools::md5sum(checkpoint_path(1L))
stopifnot(identical(quiet_run(1L, never_fit), resumed),
  identical(hash, tools::md5sum(checkpoint_path(1L))))

failed_factory <- function(...) {
  adapter <- mock_factory(...)
  adapter$run <- function(...) stop('Intentional retained terminal failure')
  adapter
}
failed <- quiet_run(2L, failed_factory)
hash <- tools::md5sum(checkpoint_path(2L))
stopifnot(length(failed$records) == 2L,
  all(!vapply(failed$records, function(v) isTRUE(v$success), logical(1))),
  identical(failed, quiet_run(2L, never_fit)),
  identical(hash, tools::md5sum(checkpoint_path(2L))))

# Malformed successes and pairing errors are recorded as terminal failures,
# including a resumed task whose new adapter disagrees with the prior learner.
malformed_factory <- function(...) {
  adapter <- mock_factory(...); original <- adapter$run
  adapter$run <- function(...) { fit <- original(...); fit$summary$p_value[1L] <- NA_real_; fit }
  adapter
}
malformed <- quiet_run(3L, malformed_factory)
stopifnot(all(!vapply(malformed$records, function(v) isTRUE(v$success), logical(1))),
  identical(malformed, quiet_run(3L, never_fit)))
for (index in seq_along(c('data', 'folds', 'primitives'))) {
  component <- c('data', 'folds', 'primitives')[index]
  mismatch_factory <- local({ field <- component; function(...) {
    adapter <- mock_factory(...); original <- adapter$run
    adapter$run <- function(data, method, seed) {
      fit <- original(data, method, seed)
      if (method == 'growing_direct') {
        if (field == 'primitives') fit$paired$primitives[1L] <- n12_fingerprint('different') else
          fit$paired[[field]] <- n12_fingerprint('different')
      }
      fit
    }
    adapter
  }})
  mismatch <- quiet_run(3L + index, mismatch_factory)
  stopifnot(isTRUE(mismatch$records$rkhs$success),
    isFALSE(mismatch$records$growing_direct$success),
    identical(mismatch, quiet_run(3L + index, never_fit)))
}
prior <- quiet_run(7L)
prior$records <- prior$records['rkhs']
n12_atomic(prior, checkpoint_path(7L))
resumed_mismatch <- quiet_run(7L, mismatch_factory)
stopifnot(identical(resumed_mismatch$records$rkhs, prior$records$rkhs),
  isFALSE(resumed_mismatch$records$growing_direct$success),
  identical(resumed_mismatch, quiet_run(7L, never_fit)))
warning_factory <- function(...) {
  adapter <- mock_factory(...); original <- adapter$run
  adapter$run <- function(...) { warning('Retained diagnostic warning'); original(...) }
  adapter
}
warned <- quiet_run(8L, warning_factory)
stopifnot(all(vapply(warned$records, function(v) isTRUE(v$success) &&
  identical(v$warnings, 'Retained diagnostic warning'), logical(1))))
prepare_failed <- quiet_run(9L, function(...) stop('Intentional adapter construction failure'))
stopifnot(all(!vapply(prepare_failed$records, function(v) isTRUE(v$success), logical(1))),
  identical(prepare_failed, quiet_run(9L, never_fit)))

# Lock rejection, malformed records, and foreign provenance never modify a
# published checkpoint. The deliberately bad files all live in this fixture.
lock <- paste0(checkpoint_path(1L), '.lock')
dir.create(lock)
n12_test_expect_error(quiet_run(1L, never_fit))
unlink(lock, recursive = TRUE)
bad_path <- file.path(test_directory, 'invalid.rds')
invalid <- resumed; invalid$signature_hash <- n12_fingerprint('foreign signature')
n12_atomic(invalid, bad_path)
n12_test_expect_error(n12_load_task(bad_path, resumed$task, signature))
invalid <- resumed; invalid$records[[2L]]$paired$data <- n12_fingerprint('unpaired')
n12_atomic(invalid, bad_path)
n12_test_expect_error(n12_load_task(bad_path, resumed$task, signature))
invalid <- resumed; names(invalid$records) <- rep('rkhs', 2L)
n12_atomic(invalid, bad_path)
n12_test_expect_error(n12_load_task(bad_path, resumed$task, signature))
invalid <- resumed; invalid$records$rkhs$summary$Lambda0[1L] <- 1.01
n12_test_expect_error(n12_check_record(invalid$records$rkhs, resumed$task, 'rkhs', signature))
for (field in c('statistic', 'critical_value', 'p_value', 'moments_tested')) {
  invalid <- resumed$records$rkhs; invalid$summary[[field]] <- NULL
  n12_test_expect_error(n12_check_record(invalid, resumed$task, 'rkhs', signature))
}
for (field in c('n_used', 'seconds', 'warnings')) {
  invalid <- resumed$records$rkhs; invalid[[field]] <- NULL
  n12_test_expect_error(n12_check_record(invalid, resumed$task, 'rkhs', signature))
}
corrupt <- file.path(test_directory, 'corrupt.rds')
writeBin(charToRaw('This is deliberately not RDS.'), corrupt)
n12_test_expect_error(n12_load_task(corrupt, resumed$task, signature))
zero <- file.path(test_directory, 'zero.rds')
stopifnot(file.create(zero))
n12_test_expect_error(n12_load_task(zero, resumed$task, signature))
stopifnot(file.copy(corrupt, checkpoint_path(10L)))
hash <- tools::md5sum(checkpoint_path(10L))
n12_test_expect_error(quiet_run(10L, never_fit))
stopifnot(identical(hash, tools::md5sum(checkpoint_path(10L))))

# Exercise the actual CLI in its own sealed fixture with retained complete
# records. No frozen estimator files are present, so any refit would fail.
cli_root <- file.path(test_directory, 'cli_bundle')
dir.create(file.path(cli_root, 'scripts'), recursive = TRUE)
cli_files <- file.path('scripts', c('common.R', 'adapters.R', 'run.R'))
stopifnot(all(file.copy(cli_files, file.path(cli_root, cli_files))))
stopifnot(file.copy('paper_config.R', file.path(cli_root, 'paper_config.R')))
n12_atomic(config, file.path(cli_root, 'study_config.rds'))
cli_manifest <- tools::md5sum(file.path(cli_root, c(cli_files, 'study_config.rds')))
names(cli_manifest) <- c(cli_files, 'study_config.rds')
n12_atomic(cli_manifest, file.path(cli_root, 'manifest.rds'))
cli_signature <- n12_signature(cli_root)
cli_result <- file.path(cli_root, config$result_root)
n12_atomic(cli_signature, file.path(cli_result, 'signature.rds'))
for (value in list(resumed, failed)) {
  value$signature_hash <- n12_fingerprint(cli_signature)
  for (method in names(value$records))
    value$records[[method]]$signature_hash <- value$signature_hash
  n12_atomic(value, file.path(cli_result, 'runs', sprintf('task_%04d.rds', value$task$id)))
}
for (id in 1:2) {
  path <- file.path(cli_result, 'runs', sprintf('task_%04d.rds', id))
  hash <- tools::md5sum(path)
  log <- file.path(cli_root, paste0('task_', id, '.log'))
  status <- suppressWarnings(n12_adapter_in_directory(cli_root, function()
    system2(file.path(R.home('bin'), 'Rscript'), c('scripts/run.R', 'task', id),
      stdout = log, stderr = log)))
  output <- paste(readLines(log, warn = FALSE), collapse = '\n')
  stopifnot(identical(hash, tools::md5sum(path)), grepl('ALREADY ATTEMPTED', output))
  if (id == 1L) stopifnot(status == 0L, grepl('TASK COMPLETE', output)) else
    stopifnot(status != 0L, grepl('Failed records retained', output))
}

# Only the freshly-created unit-test fixture directory is removed.
unlink(test_directory, recursive = TRUE)
cat('TESTS PASSED: 5000-task contract; parent generator parity; frozen pointwise subgrid calibration; compact atomic persistence; completed and partial resume; retained failures and nonzero CLI; paired fingerprints; corrupt/foreign/duplicate checkpoint rejection.\n')
