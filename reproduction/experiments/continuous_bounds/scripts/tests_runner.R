# Fast control-flow integration test with an explicitly mocked estimator.
# This is not a numerical validation of the estimator; real fits run in the LPC gate.
source('scripts/common.R')
c<-stress_config();tasks<-stress_tasks();spec<-stress_spec()
t<-tasks[1501,];signature<-list(test=TRUE)
tmp<-tempfile('stress-runner-test-');dir.create(tmp)
e<-new.env(parent=globalenv())
e$c<-c;e$tasks<-tasks;e$spec<-spec;e$t<-t;e$id<-1501L;e$a<-1.1;e$b<-.15
e$signature<-signature;e$path<-file.path(tmp,'task.rds');e$lock<-file.path(tmp,'task.lock')
e$targets<-expand.grid(setting=1:4,Lambda=c(2,3),class=c$classes,stringsAsFactors=FALSE)
e$targets$lower<-.2;e$targets$upper<-.8
bad<-e$targets$setting==4 & e$targets$Lambda==2
e$targets$lower[bad]<-e$targets$upper[bad]<-NA_real_
e$x5g_generate<-function(n,...)data.frame(Y=rep(0,n))
e$make_tz_stratified_folds<-function(d,K,seed)rep(1:K,length.out=nrow(d))
e$x5_adaptive_provider<-function(...) {
 f<-function(...)NULL;attr(f,'gamma_cache')<-new.env();f
}
e$x5s_threshold_provider<-function(cl,store) {
 store$fit<-list(lower=list(convergence=0),upper=list(convergence=0));function(...)NULL
}
e$project_gamma_to_compatible_interval<-function(...)NULL
e$estimate_sharp_bounds_one_step<-function(d,Lambda0_T1,fold_id,bound_component_provider_factory,...) {
 force(bound_component_provider_factory)
 if(Lambda0_T1==2)stop('deliberate test numerical failure')
 n<-nrow(d);v<-rep(c(-.1,.1),length.out=n)
 signals<-data.frame(mode='rkhs',ate_lower=.2+v,ate_upper=.8+v)
 ep<-data.frame(mode='rkhs',endpoint=c('ate_lower','ate_upper'),one_step_estimate=c(.2,.8),
   standard_error=sd(v)/sqrt(n))
 ep$confidence_lower<-ep$one_step_estimate-qnorm(.975)*ep$standard_error
 ep$confidence_upper<-ep$one_step_estimate+qnorm(.975)*ep$standard_error
 list(fold_id=fold_id,one_step=ep,one_step_signals=signals,
   fold_diagnostics=data.frame(gamma_projection_rate_validation=0),eif_diagnostics=data.frame())
}
p<-parse('scripts/run.R')
r<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name('run')),as.list(p))
stopifnot(length(r)==1);eval(r[[1]],e);e$run()
z<-readRDS(e$path)
stopifnot(length(z$records)==8,sum(vapply(z$records,`[[`,logical(1),'success'))==4,
 all(vapply(z$records,function(v)if(!v$compatible)is.na(v$outward) else isTRUE(v$success),logical(1))),
 !dir.exists(e$lock))
hash<-tools::md5sum(e$path);e$run();stopifnot(identical(hash,tools::md5sum(e$path)))
cat('PASS: mocked runner retains four failed incompatible fits, completes four compatible fits, releases its lock and resumes without modifying the checkpoint.\n')
