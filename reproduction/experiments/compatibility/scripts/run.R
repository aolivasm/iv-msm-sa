#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,
           VECLIB_MAXIMUM_THREADS=1,RCPP_PARALLEL_NUM_THREADS=1)
source('scripts/common.R')
source('scripts/adapters.R')
args <- commandArgs(TRUE); stage <- if(length(args))args[1L] else 'preflight'
stopifnot(stage%in%c('preflight','task','smoke','gate_check','status','collect','pack'))
bundle_root <- normalizePath('.')
cfg <- n12_config(); tasks <- n12_tasks(cfg)
root <- cfg$result_root; signature <- n12_signature(bundle_root)
meta <- file.path(root,'signature.rds')
if(file.exists(meta))stopifnot(identical(readRDS(meta),signature)) else {
  if(stage!='preflight')stop('Run preflight once before any tasks.')
  n12_atomic(signature,meta,overwrite=FALSE)
}
if(stage=='preflight') {
  stopifnot(nrow(tasks)==paper_profile$tasks,all(tasks$n==cfg$n))
  cat('PREFLIGHT:',nrow(tasks),'datasets; n =',cfg$n,'; learners:',cfg$methods,
      '; sensitivity values:',cfg$Lambda,'\n')
  cat('Configuration and source fingerprints verified; no fits performed.\n')
  quit(status=0)
}
if(stage%in%c('task','smoke')) {
  if(length(args)!=2L)stop('Usage: Rscript scripts/run.R ',stage,' TASK_ID')
  id <- suppressWarnings(as.integer(args[2L]))
  z <- n12_run_task(id,signature,bundle_root,root,smoke=stage=='smoke')
  if(!setequal(names(z$records),cfg$methods)||
     !all(vapply(z$records,function(v)isTRUE(v$success),logical(1))))
    stop('Failed records retained. No automatic refit; inspect checkpoint before recovery.')
  cat('TASK COMPLETE',id,'; checkpoint verified.\n')
  quit(status=0)
}
gate_ids <- seq.int(1L,nrow(tasks),by=cfg$reps)
ids <- if(stage=='gate_check')gate_ids else tasks$id
status <- detail <- diagnostics <- list(); issues <- character()
for(id in ids) {
  task <- tasks[id,,drop=FALSE]; path <- file.path(root,'runs',sprintf('task_%04d.rds',id))
  z <- tryCatch({
    if(file.exists(paste0(path,'.lock')))stop('active or stale lock; inspect before collection')
    if(!file.exists(path))stop('missing checkpoint')
    v <- n12_load_task(path,task,signature)
    if(!setequal(names(v$records),cfg$methods))stop('incomplete checkpoint')
    v
  },error=function(e){issues<<-c(issues,paste(id,conditionMessage(e)));NULL})
  if(is.null(z))next
  for(method in cfg$methods) {
    v <- z$records[[method]]
    status[[length(status)+1L]] <- cbind(task,method=method,success=v$success,
      seconds=v$seconds,error=v$error,warnings=length(v$warnings))
    if(v$success) {
      detail[[length(detail)+1L]] <- cbind(task,method=method,v$summary)
      if(stage%in%c('collect','pack'))diagnostics[[length(diagnostics)+1L]] <-
        c(list(id=id,outcome=task$outcome,method=method),v[c('gamma_selection',
          'gamma_diagnostics','probability_diagnostics','paired','settings','warnings',
          'witness_candidate_counts','fold_sizes','projection_applied','q_fits')])
    }
  }
}
status <- if(length(status))do.call(rbind,status) else data.frame()
successes <- if(nrow(status))sum(status$success) else 0L
cat('Successful learner runs:',successes,'of',length(cfg$methods)*length(ids),'\n')
if(length(issues)){cat('Datasets requiring attention:',length(issues),'\n');print(issues)}
if(nrow(status)&&any(!status$success))print(status[!status$success,],row.names=FALSE)
if(stage=='status')quit(status=0)
if(length(issues)||nrow(status)!=length(cfg$methods)*length(ids)||!all(status$success))
  stop('Not complete and successful. No release/archive; existing checkpoints retained.')
detail <- do.call(rbind,detail)
if(stage=='gate_check') {
  cat('Configured gate checkpoints verified.\n')
  quit(status=0)
}
summary <- do.call(rbind,lapply(split(detail,interaction(detail$outcome,detail$setting,
                                                        detail$method,detail$Lambda0,drop=TRUE)),function(g) {
  stopifnot(nrow(g)==cfg$reps,!anyDuplicated(g$replication),setequal(g$replication,seq_len(cfg$reps)))
  ci <- binom.test(sum(g$reject),nrow(g))$conf.int
  data.frame(g[1,c('outcome','setting','type','a','b','n','method','Lambda0')],
    attempted=cfg$reps,successful=nrow(g),rejected=sum(g$reject),rejection=mean(g$reject),
    mc_low=ci[1],mc_high=ci[2])
}))
expected_learners <- nrow(tasks)*length(cfg$methods)
expected_tests <- expected_learners*length(cfg$Lambda)
stopifnot(nrow(detail)==expected_tests,nrow(status)==expected_learners)
if(paper_profile$schema_n24) {
  summary <- summary[c('n','setting','method','Lambda0','attempted','successful','rejected','rejection','mc_low','mc_high')]
  names(summary)[names(summary)=='Lambda0'] <- 'Lambda'
}
audit <- list(expected_datasets=nrow(tasks),expected_learners=expected_learners,expected_tests=expected_tests,
  successful_learners=successes,failures=0L,all_tasks_complete=TRUE,
  paired_data_folds_primitives_verified=TRUE,source_manifest_verified=TRUE)
result <- list(signature=signature,status=status,detail=detail,summary=summary,audit=audit,
  note=paste('Pointwise maximum is primary. The main-text RKHS learner is used.',
    'Complete retained diagnostics are separate from this compact summary.',
    'Source tuning differs by outcome and is preserved from n6000.',
    'Only stress and null designs; no post-hoc test selection.',
    'Other grid-simultaneous outputs apply to Lambda1,1.1 only. Nonrejection does not prove compatibility.'))
n12_atomic(result,file.path(root,'summary.rds'))
n12_atomic(diagnostics,file.path(root,'diagnostics.rds'))
print(summary,row.names=FALSE)
cat('COMPLETE:',nrow(tasks),'datasets /',expected_learners,'successful learners /',expected_tests,'primary tests.\n')
