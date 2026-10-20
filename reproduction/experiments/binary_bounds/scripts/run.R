Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,RCPP_PARALLEL_NUM_THREADS=1)
source('scripts/common.R')
args<-commandArgs(TRUE);stage<-if(length(args))args[1] else 'preflight'
stopifnot(stage%in%c('preflight','smoke','task','check_gate','collect'))
m<-readRDS('manifest.rds');stopifnot(identical(unname(m),unname(tools::md5sum(names(m)))))
pk<-c('SuperLearner','glmnet','ranger','xgboost','nloptr','nnls');stopifnot(all(vapply(pk,requireNamespace,logical(1),quietly=TRUE)))
pop<-readRDS('targets/population.rds');stopifnot(identical(pop$config,config),pop$max_refinement<1e-4)
sig<-list(config=config,manifest=m,R=as.character(getRversion()),packages=sapply(pk,function(p)as.character(packageVersion(p))))
root<-'results/binary_x5_stress_v1';meta<-file.path(root,'signature.rds')
if(file.exists(meta))stopifnot(identical(readRDS(meta),sig)) else {stopifnot(stage=='preflight');atomic(sig,meta)}
if(stage=='preflight'){
 cat('PREFLIGHT: stress only, n6000, four settings x 500 reps = 2000 datasets.\n')
 cat('RKHS; Lambda 1.8,2,2.5,3; 8000 interval fits. Compact checkpoints.\n');quit(status=0)
}
if(stage%in%c('task','smoke')){
 id<-if(length(args)>1)as.integer(args[2]) else 1501L
 z<-run_task(id,file.path(root,if(stage=='smoke')'smoke' else 'runs'),sig,n=if(stage=='smoke')600L else config$n)
 # Failures remain in checkpoints, but process failure makes the gate block safely.
 if(!all(vapply(z$records,`[[`,logical(1),'success')))stop('Numerical failure retained; inspect checkpoint.')
 quit(status=0)
}
ids<-if(stage=='check_gate')c(1,501,1001,1501) else tasks$id
detail<-status<-list()
for(id in ids){
 p<-file.path(root,'runs',sprintf('task_%04d.rds',id));stopifnot(file.exists(p))
 z<-readRDS(p);t<-tasks[id,,drop=FALSE];check(z,t,sig);stopifnot(setequal(names(z$records),config$methods))
 for(method in config$methods){
  v<-z$records[[method]];status[[length(status)+1L]]<-cbind(t,method,success=v$success,error=v$error)
  if(v$success)for(f in v$fits){
   q<-subset(pop$targets,setting==t$setting & Lambda==f$Lambda0)
   e<-f$endpoint_intervals[match(c('ate_lower','ate_upper'),trimws(f$endpoint_intervals$endpoint)),]
   detail[[length(detail)+1L]]<-cbind(t,method,Lambda=f$Lambda0,
    estimate_lower=e$estimate[1],estimate_upper=e$estimate[2],se_lower=e$standard_error[1],se_upper=e$standard_error[2],
    truth_lower=q$lower,truth_upper=q$upper,lower_error=e$estimate[1]-q$lower,upper_error=e$estimate[2]-q$upper,
    point_width=e$estimate[2]-e$estimate[1],confidence_width=e$ci_upper[2]-e$ci_lower[1],
    outward=e$ci_lower[1]<=q$lower & e$ci_upper[2]>=q$upper,
    lower_centered=e$ci_lower[1]<=q$lower & e$ci_upper[1]>=q$lower,
    upper_centered=e$ci_lower[2]<=q$upper & e$ci_upper[2]>=q$upper,
    projection_max=max(f$diagnostics$compatibility_projection_rate),warnings=length(v$warnings))
  }
 }
}
status<-do.call(rbind,status)
if(stage=='check_gate'){stopifnot(all(status$success));cat('GATE PASSED: four n6000 datasets, four learner runs, 16 interval fits.\n');quit(status=0)}
detail<-do.call(rbind,detail)
summary<-do.call(rbind,lapply(split(detail,interaction(detail$setting,detail$method,detail$Lambda,drop=TRUE)),function(g){
 cbind(g[1,c('setting','method','Lambda')],attempted=500,successful=nrow(g),outward=mean(g$outward),
  lower_centered=mean(g$lower_centered),upper_centered=mean(g$upper_centered),width=mean(g$confidence_width),
  lower_bias=mean(g$lower_error),upper_bias=mean(g$upper_error),
  lower_rmse=sqrt(mean(g$lower_error^2)),upper_rmse=sqrt(mean(g$upper_error^2)))
}))
atomic(list(status=status,detail=detail,summary=summary,signature=sig,population=pop),file.path(root,'summary.rds'))
print(summary,row.names=FALSE);cat(nrow(status),'learner runs;',sum(status$success),'successful. Failures retained.\n')
