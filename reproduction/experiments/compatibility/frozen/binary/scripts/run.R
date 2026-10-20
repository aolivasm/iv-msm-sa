Sys.setenv(OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1,
           VECLIB_MAXIMUM_THREADS=1, RCPP_PARALLEL_NUM_THREADS=1)
source('scripts/common.R')
args <- commandArgs(TRUE)
stage <- if (length(args)) args[1L] else 'preflight'
stopifnot(stage %in% c('preflight','smoke','task','check_gate','collect','status'))
manifest <- readRDS('manifest.rds')
stopifnot(identical(unname(manifest), unname(tools::md5sum(names(manifest)))))
pk <- c('SuperLearner','glmnet','ranger','xgboost','nloptr','nnls','digest')
stopifnot(all(vapply(pk, requireNamespace, logical(1), quietly=TRUE)))
signature <- list(config=config, manifest=manifest, R=as.character(getRversion()),
                  packages=sapply(pk, function(p) as.character(packageVersion(p))))
meta <- file.path(config$result_root, 'signature.rds')
if (file.exists(meta)) stopifnot(identical(readRDS(meta), signature)) else {
  stopifnot(stage == 'preflight'); atomic(signature, meta, overwrite=FALSE)
}
if (stage == 'preflight') {
  cat('PREFLIGHT:', nrow(tasks), 'datasets; n6000;', nrow(config$settings),
      'designs x 500 repetitions; two gamma learners; five Lambdas.\n')
  cat('Primary: pointwise maximum, nominal 5%; five outer folds; 399 multipliers.\n')
  cat('RKHS: sigma^2=M/2, penalty-only two-repeat three-fold CV. Growing-direct unchanged.\n')
  cat('No compatibility projection, classifier gamma, binary PMF reconstruction, or q fitting.\n')
  quit(status=0)
}
if (stage %in% c('smoke','task')) {
  id <- if (length(args)>1L) as.integer(args[2L]) else 1501L
  z <- run_task(id, signature, smoke=stage=='smoke')
  if (!all(vapply(z$records, function(v) isTRUE(v$success), logical(1))))
    stop('Failed learner retained; inspect checkpoint. No automatic replacement.')
  quit(status=0)
}
ids <- if (stage=='check_gate') gate_ids else tasks$id
status <- detail <- diagnostics <- list()
problems <- character()
for (id in ids) {
  task <- tasks[id, , drop=FALSE]
  path <- file.path(config$result_root, 'runs', sprintf('task_%04d.rds',id))
  z <- tryCatch({
    if (!file.exists(path)) stop('missing checkpoint')
    v <- load_task(path, task, signature)
    if (!setequal(names(v$records), config$methods)) stop('incomplete checkpoint')
    v
  }, error=function(e) {problems <<- c(problems, paste(id, conditionMessage(e))); NULL})
  if (is.null(z)) next
  for (method in config$methods) {
    v <- z$records[[method]]
    status[[length(status)+1L]] <- cbind(task, method=method,
      success=v$success, error=v$error, seconds=v$seconds, warnings=length(v$warnings))
    if (v$success) {
      detail[[length(detail)+1L]] <- cbind(task, method=method, v$summary)
      diagnostics[[length(diagnostics)+1L]] <- list(id=id, method=method,
        gamma_selection=v$gamma_selection, gamma_diagnostics=v$gamma_diagnostics,
        probability_diagnostics=v$probability_diagnostics, warnings=v$warnings,
        paired=v$paired, settings=v$settings)
    }
  }
}
status <- if(length(status)) do.call(rbind,status) else data.frame()
successes <- if(nrow(status)) sum(status$success) else 0L
cat('Successful learners:', successes, 'of', length(ids)*length(config$methods), '\n')
if(length(problems)) {cat('Missing/unreadable/incomplete datasets:',length(problems),'\n'); print(head(problems,30))}
if(stage=='status') quit(status=0)
if(length(problems)) stop('Collection incomplete. Checkpoints preserved; inspect listed tasks.')
if(stage=='check_gate') {
  if(!all(status$success)) {print(status[!status$success,]); stop('Numerical gate failed; full study blocked.')}
  elapsed <- aggregate(seconds ~ id + setting, status, sum)
  sizes <- file.info(file.path(config$result_root,'runs',sprintf('task_%04d.rds',ids)))$size
  report <- list(signature_hash=fingerprint(signature), status=status,
    task_seconds=elapsed, checkpoint_bytes=sizes,
    projected_content_bytes=ceiling(max(sizes)*nrow(tasks)*1.5),
    note='Small timing gate; extrapolation is approximate and excludes queue/startup/I/O variability.')
  atomic(report,file.path(config$result_root,'gate_report.rds'))
  cat('GATE PASSED:',length(ids),'actual n6000 datasets;',nrow(status),'paired learner runs.\n')
  cat('Observed seconds per dataset:',paste(round(elapsed$seconds,1),collapse=', '),'\n')
  cat('Projected checkpoint content with 50% allowance:',round(report$projected_content_bytes/1024^2,1),'MiB.\n')
  quit(status=0)
}
detail <- if(length(detail)) do.call(rbind,detail) else data.frame()
labels <- population_labels()
summary <- if(nrow(detail)) do.call(rbind,lapply(split(detail,
  interaction(detail$setting,detail$method,detail$Lambda0,drop=TRUE)),function(g) {
  ci<-binom.test(sum(g$reject),nrow(g))$conf.int
  data.frame(g[1,c('setting','n','method','Lambda0')],attempted=config$reps,
    successful=nrow(g), rejection=mean(g$reject), mc_lower=ci[1],mc_upper=ci[2],
    ordered_rejection=mean(g$ordered_reject),
    grid_simultaneous_rejection=mean(g$grid_simultaneous_reject))
})) else data.frame()
if(nrow(summary)) summary <- merge(summary,labels,by=c('setting','Lambda0'))
audit <- list(expected_datasets=nrow(tasks), expected_learners=2*nrow(tasks),
  expected_tests=2*nrow(tasks)*length(config$Lambda), saved_learners=nrow(status),
  successful_learners=successes, failures=if(nrow(status)) sum(!status$success) else 0L,
  all_tasks_complete=!length(problems), all_learners_successful=all(status$success),
  source_manifest_verified=TRUE, paired_data_folds_primitives_verified=TRUE)
result <- list(signature=signature,status=status,detail=detail,summary=summary,
  diagnostics=diagnostics,audit=audit,population_labels=labels,
  note='Primary rejection rate is pointwise maximum; secondary variants are not selected post hoc. Nonrejection does not prove compatibility. Population compatible labels are numerical except the analytic b=0 control.')
atomic(result,file.path(config$result_root,'summary.rds'))
print(summary,row.names=FALSE)
cat(successes,'of',nrow(status),'learner runs successful. Summary retains failures.\n')
if(!all(status$success)) stop('Failed learners retained. Summary is partial scientifically; inspect before reporting complete.')
