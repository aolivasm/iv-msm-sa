#!/usr/bin/env Rscript
source('scripts/common.R')
source('scripts/adapters.R')
args <- commandArgs(TRUE)
stage <- if(length(args)) args[1] else 'preflight'
stopifnot(stage %in% c('preflight','task','smoke','collect','gate_check'))
cfg <- n12_config(); signature <- n12_signature(); tasks <- n12_tasks()
root <- cfg$result_root
dir.create(root,recursive=TRUE,showWarnings=FALSE)
sf <- file.path(root,'signature.rds')
if(file.exists(sf)) stopifnot(identical(readRDS(sf),signature)) else n12_atomic(signature,sf,FALSE)
if(stage=='preflight') {
  stopifnot(nrow(tasks)==500L,all(tasks$setting==4L),all(tasks$n==24000L),
    identical(cfg$methods,'rkhs'),identical(cfg$Lambda,1))
  cat('PREFLIGHT: binary Setting 4 only; a=1.1 b=.15; n24000; RKHS only; Lambda1; 500 repetitions.\n')
  cat('Frozen binary tuning: M/2 kernel variance; penalty-only 2 x 3-fold CV; five outer folds; 399 multipliers.\n')
  quit(status=0)
}
if(stage %in% c('task','smoke')) {
  stopifnot(length(args)==2L)
  z <- n12_run_task(as.integer(args[2]),signature,result_root=root,smoke=stage=='smoke')
  if(!setequal(names(z$records),'rkhs') || !isTRUE(z$records$rkhs$success))
    stop('Failed record retained; inspect before retry. No automatic refit.')
  cat('TASK VERIFIED:',args[2],'\n');quit(status=0)
}
ids <- if(stage=='gate_check') 1L else seq_len(500L)
records <- lapply(ids,function(id) {
  f <- file.path(root,'runs',sprintf('task_%04d.rds',id))
  if(!file.exists(f)||file.exists(paste0(f,'.lock'))) stop('Missing or locked task: ',id)
  z <- n12_load_task(f,tasks[id,,drop=FALSE],signature)
  stopifnot(setequal(names(z$records),'rkhs'),isTRUE(z$records$rkhs$success))
  z$records$rkhs
})
if(stage=='gate_check') {
  cat('Replication 1 verified; seconds:',records[[1]]$seconds,'\n')
  cat('Checkpoint bytes:',file.info(file.path(root,'runs/task_0001.rds'))$size,'\n')
  cat('Review LSF peak memory and runtime before submitting replications 2-500.\n')
  quit(status=0)
}
detail <- do.call(rbind,lapply(records,function(v) {
  task <- v$task;rownames(task)<-NULL
  cbind(task,method=v$method,v$summary)
}))
stopifnot(nrow(detail)==500L,!anyDuplicated(detail$replication),!anyNA(detail$reject))
summary <- data.frame(n=24000,setting=4,method='rkhs',Lambda=1,
  attempted=500,successful=500,rejected=sum(detail$reject),rejection=mean(detail$reject))
ci <- binom.test(sum(detail$reject),500)$conf.int
summary$mc_low<-ci[1];summary$mc_high<-ci[2]
n12_atomic(list(signature=signature,summary=summary,detail=detail,records=records),
  file.path(root,'summary.rds'))
print(summary,row.names=FALSE)
cat('COMPLETE: 500 of 500 replications verified; no exclusions.\n')
