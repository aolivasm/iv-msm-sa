Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,RCPP_PARALLEL_NUM_THREADS=1)
source('scripts/common.R')
args<-commandArgs(TRUE);stage<-if(length(args))args[1] else 'preflight'
stopifnot(stage%in%c('preflight','smoke','task','gate_check','collect'))
m<-readRDS('manifest.rds');stopifnot(identical(unname(m),unname(tools::md5sum(names(m)))))
pk<-c('SuperLearner','glmnet','ranger','xgboost','nloptr','nnls')
stopifnot(all(vapply(pk,requireNamespace,logical(1),quietly=TRUE)))
signature<-list(config=config,manifest=m,R=as.character(getRversion()),packages=sapply(pk,function(p)as.character(packageVersion(p))))
root<-'results/x5_compatibility_v1';meta<-file.path(root,'signature.rds')
if(file.exists(meta))stopifnot(identical(readRDS(meta),signature)) else {stopifnot(stage=='preflight');atomic(signature,meta)}
if(stage=='preflight'){
 cat('PREFLIGHT: 4500 datasets; n6000; 500 reps; 4 baseline + 4 stress + b0 null; two gamma learners; six Lambdas.\n')
 cat('Pointwise maximum is primary. Five outer folds, 399 multiplier draws. No compatibility projection; no q fitting.\n')
 quit(status=0)
}
if(stage%in%c('task','smoke')){
 id<-if(length(args)>1)as.integer(args[2]) else 1501L
 z<-run_task(id,root,signature,smoke=stage=='smoke')
 if(stage=='smoke')stopifnot(all(vapply(z$records,`[[`,logical(1),'success')))
 quit(status=0)
}
ids<-if(stage=='gate_check')gate_ids else tasks$id
status<-detail<-list()
for(id in ids){
 t<-tasks[id,,drop=FALSE];z<-load_task(file.path(root,'runs',sprintf('task_%04d.rds',id)),t,signature)
 stopifnot(setequal(names(z$records),config$methods))
 for(method in config$methods){
  v<-z$records[[method]];j<-length(status)+1L
  status[[j]]<-cbind(t,method=method,success=v$success,error=v$error,seconds=v$seconds)
  if(v$success)detail[[length(detail)+1L]]<-cbind(t,method=method,v$summary)
 }
}
status<-do.call(rbind,status)
if(stage=='gate_check'){
 if(!all(status$success)){print(status[!status$success,]);stop('Gate failed: full experiment blocked. Failures retained.')}
 cat('GATE PASSED: 9 actual-size datasets, 18 learner runs, 108 tests.\n');quit(status=0)
}
detail<-if(length(detail))do.call(rbind,detail) else data.frame()
summary<-if(nrow(detail))do.call(rbind,lapply(split(detail,interaction(detail$design_id,detail$method,detail$Lambda0,drop=TRUE)),function(g){
 ci<-binom.test(sum(g$reject),nrow(g))$conf.int
 cbind(g[1,c('design_id','type','setting','a','b','n','method','Lambda0')],attempted=500L,successful=nrow(g),
  rejection=mean(g$reject),mc_low=ci[1],mc_high=ci[2],ordered_rejection=mean(g$ordered_reject),
  grid_simultaneous_rejection=mean(g$grid_simultaneous_reject))
})) else data.frame()
atomic(list(signature=signature,status=status,detail=detail,summary=summary,
 note='Rejection rates among successful tests, failures explicit. Population compatibility labels/frontiers require separate validation. b=0 null is compatible at all grid values; nonrejection is not proof of compatibility.'),file.path(root,'summary.rds'))
cat(nrow(status),'learner checkpoints;',sum(status$success),'successful;',sum(!status$success),'failed.\n')
