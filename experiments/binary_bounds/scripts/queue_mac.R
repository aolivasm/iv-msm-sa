workers<-as.integer(Sys.getenv('BINARY_X5_WORKERS','8'));stopifnot(workers>=1,workers<=8)
run_one<-function(id){
 log<-sprintf('logs/task_%04d.log',id)
 status<-system2(file.path(R.home('bin'),'Rscript'),c('scripts/run.R','task',id),stdout=log,stderr=log)
 cat(format(Sys.time()),'task',id,'exit',status,'\n');flush.console();status
}
gate<-c(1L,501L,1001L,1501L)
for(phase in c('gate','full')){
 ids<-if(phase=='gate')gate else setdiff(1:2000,gate)
 cores<-if(phase=='gate')min(2L,workers) else workers
 cat('Stage',phase,':',length(ids),'datasets;',cores,'workers.\n');flush.console()
 codes<-parallel::mclapply(ids,run_one,mc.cores=cores,mc.preschedule=FALSE,mc.set.seed=FALSE)
 if(any(!vapply(codes,function(x)is.numeric(x)&&length(x)==1&&x==0,logical(1))))stop('Process failures: inspect logs. Checkpoints retained.')
 if(phase=='gate')stopifnot(system2(file.path(R.home('bin'),'Rscript'),c('scripts/run.R','check_gate'))==0)
}
