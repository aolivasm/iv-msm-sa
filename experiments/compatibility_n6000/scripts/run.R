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
  stopifnot(nrow(tasks)==5000L,all(table(tasks$outcome,tasks$setting)==500L))
  cat('PREFLIGHT: 5000 datasets; n6000; 500 repetitions; 4 stress + null per outcome.\n')
  cat('Continuous and binary; both gamma learners; Lambda=1,1.1: 10000 learner runs / 20000 pointwise tests.\n')
  cat('Five outer folds; 399 Gaussian multipliers; no compatibility projection or q fitting.\n')
  cat('Frozen outcome-specific RKHS tuning, not silently harmonized; growing-sieve rule unchanged.\n')
  cat('Ten actual-size gates are included in the 5000; full arrays exclude those tasks.\n')
  quit(status=0)
}
if(stage%in%c('task','smoke')) {
  if(length(args)!=2L)stop('Usage: Rscript scripts/run.R ',stage,' TASK_ID')
  id <- suppressWarnings(as.integer(args[2L]))
  z <- n12_run_task(id,signature,bundle_root,root,smoke=stage=='smoke')
  if(!setequal(names(z$records),cfg$methods)||
     !all(vapply(z$records,function(v)isTRUE(v$success),logical(1))))
    stop('Failed records retained. No automatic refit; inspect checkpoint before recovery.')
  cat('TASK COMPLETE',id,'; 2 successful learners; 4 pointwise tests; checkpoint verified.\n')
  quit(status=0)
}
gate_ids <- seq.int(1L,4501L,by=500L)
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
cat('Successful learner runs:',successes,'of',2L*length(ids),'\n')
if(length(issues)){cat('Datasets requiring attention:',length(issues),'\n');print(issues)}
if(nrow(status)&&any(!status$success))print(status[!status$success,],row.names=FALSE)
if(stage=='status')quit(status=0)
if(length(issues)||nrow(status)!=2L*length(ids)||!all(status$success))
  stop('Not complete and successful. No release/archive; existing checkpoints retained.')
detail <- do.call(rbind,detail)
if(stage=='gate_check') {
  source('scripts/storage.R')
  elapsed <- aggregate(seconds~id+outcome+setting,status,sum)
  paths <- file.path(root,'runs',sprintf('task_%04d.rds',ids))
  sizes <- file.info(paths)$size
  allocated <- vapply(paths,function(path) {
    raw <- system2('du',c('-k',shQuote(path)),stdout=TRUE,stderr=TRUE)
    if(!is.null(attr(raw,'status'))||length(raw)!=1L)stop('Cannot measure allocated checkpoint space: ',path)
    kb <- suppressWarnings(as.numeric(strsplit(trimws(raw),'[[:space:]]+')[[1L]][1L]))
    if(!is.finite(kb)||kb<0)stop('Unrecognized allocated-size report: ',paste(raw,collapse=' '))
    1024*kb
  },numeric(1))
  # Actual allocated bytes matter on GPFS; content bytes alone understate use.
  projected <- ceiling(max(pmax(sizes,allocated,256*1024))*nrow(tasks)*2 + 1024^3)
  storage <- n12_check_storage(required_bytes=max(8*1024^3,projected))
  report <- list(signature_hash=n12_fingerprint(signature),status=status,
    task_seconds=elapsed,checkpoint_bytes=sizes,checkpoint_allocated_bytes=allocated,
    projected_allocated_bytes=projected,
    quota=storage,gate_ids=gate_ids,
    note='Timings include both fitted learners; forecast is not an ETA. LSF per-gate memory usage is in retained job logs.')
  n12_atomic(report,file.path(root,'gate_report.rds'))
  print(elapsed,row.names=FALSE)
  cat('Maximum checkpoint bytes:',max(sizes),'\n')
  cat('GATE PASSED: ten actual n6000 datasets / twenty learners / forty tests.\n')
  cat('Projected allocated checkpoint/log allowance:',round(projected/1024^3,2),'GiB.\n')
  cat('Remaining arrays may now start under the SAME existing 60-worker group.\n')
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
stopifnot(nrow(detail)==20000L,nrow(summary)==40L,nrow(status)==10000L)
audit <- list(expected_datasets=5000L,expected_learners=10000L,expected_tests=20000L,
  successful_learners=successes,failures=0L,all_tasks_complete=TRUE,
  paired_data_folds_primitives_verified=TRUE,source_manifest_verified=TRUE)
result <- list(signature=signature,status=status,detail=detail,summary=summary,audit=audit,
  note=paste('Pointwise maximum is primary. Both learners share each dataset/folds/primitives.',
    'Complete retained diagnostics are separate from this compact summary.',
    'Source tuning differs by outcome and is preserved from n6000.',
    'Only stress and null designs; no post-hoc test selection.',
    'Other grid-simultaneous outputs apply to Lambda1,1.1 only. Nonrejection does not prove compatibility.'))
n12_atomic(result,file.path(root,'summary.rds'))
n12_atomic(diagnostics,file.path(root,'diagnostics.rds'))
print(summary,row.names=FALSE)
cat('COMPLETE: 5000 datasets / 10000 successful learners / 20000 primary tests.\n')
