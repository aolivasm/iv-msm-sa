source('R/x5_adaptive_gamma_provider.R')
source('R/conservative_compatibility_test.R')
settings<-data.frame(setting=1:4,a=c(.55,1.1,1.65,1.1),b=c(.075,.075,.075,.15))
design<-rbind(cbind(design_id=1:4,type='baseline',settings,stress=FALSE),
 cbind(design_id=5:8,type='stress',settings,stress=TRUE),
 data.frame(design_id=9L,type='null',setting=0L,a=1.1,b=0,stress=FALSE))
config<-list(n=6000L,repetitions=500L,design=design,Lambda=c(1,1.1,1.25,1.5,1.8,2),
 methods=c('rkhs_tuned','growing_direct'),folds=5L,bootstrap=399L,allocation_df=5L,
 primary='reject: pointwise maximum; ordered and grid-simultaneous rejection are secondary',
 seed_base=820000L,primitives=sfp_contract()$primitive_library)
tasks<-design[rep(1:9,each=500),];rownames(tasks)<-NULL
tasks$replication<-rep(1:500,9);tasks$id<-1:4500;tasks$n<-6000L;tasks$seed<-820000L+tasks$replication
gate_ids<-(0:8)*500+1L
atomic<-function(object,path){
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 tmp<-tempfile('unpublished_',dirname(path));on.exit(if(file.exists(tmp))unlink(tmp))
 withCallingHandlers(saveRDS(object,tmp,compress='gzip'),warning=function(w)stop(conditionMessage(w)))
 stopifnot(file.info(tmp)$size>0,identical(readRDS(tmp),object))
 if(!file.rename(tmp,path))stop('Cannot publish: ',path)
 stopifnot(identical(readRDS(path),object));invisible(path)
}
check_record<-function(v,t,method,signature){
 stopifnot(identical(v$task,t),identical(v$method,method),identical(v$signature,signature),
  is.logical(v$success),length(v$success)==1,!is.na(v$success))
 if(v$success){s<-v$summary;stopifnot(nrow(s)==length(config$Lambda),
   identical(as.numeric(s$Lambda0),config$Lambda),all(s$Lambda0==s$Lambda1),
   all(is.finite(s$p_value)),all(s$p_value>=0&s$p_value<=1),!anyNA(s$reject))}
 else stopifnot(length(v$error)==1,nzchar(v$error))
}
load_task<-function(path,t,signature){
 if(!file.exists(path))return(list(task=t,signature=signature,records=list()))
 z<-readRDS(path);stopifnot(identical(z$task,t),identical(z$signature,signature),is.list(z$records))
 if(length(z$records)){
  stopifnot(!anyDuplicated(names(z$records)),all(names(z$records)%in%config$methods))
  for(m in names(z$records))check_record(z$records[[m]],t,m,signature)
 }
 z
}
new_session<-function(){
 engine<-new.env(parent=globalenv());sys.source('R/learned_witness_compatibility_test.R',engine)
 original<-engine$learned_witness_probability_fit;cache<-new.env(parent=emptyenv())
 engine$learned_witness_probability_fit<-function(...){
  args<-list(...);key<-as.character(args$library_seed)
  if(exists(key,cache,inherits=FALSE)){z<-get(key,cache);stopifnot(identical(args,z$args));return(z$fit)}
  f<-do.call(original,args);assign(key,list(args=args,fit=f),cache);f
 }
 # Testing must never invoke compatibility projection, even inadvertently.
 engine$project_gamma_to_compatible_interval<-function(...)stop('Forbidden compatibility projection in testing.')
 last_folds<-NULL
 run<-function(d,method,seed,B=config$bootstrap){
  provider<-x5_adaptive_provider(method)
  engine$fit_gamma_rkhs<-function(train,nuisance_train,t_value,x_cols,seed,...)
   provider(fold=seed,treatment=t_value,train=train,nuisance_train=nuisance_train,x_cols=x_cols,seed=seed)
  f<-engine$test_learned_witness_dual_msm_compatibility(d,config$Lambda,
   x_cols=paste0('X',1:5),folds=config$folds,primitive_learner='super_learner',gamma_learner='rkhs_cv',
   super_learner_candidates=config$primitives,super_learner_cv_folds=3,super_learner_ranger_num_trees=500,
   gamma_cv_selection_rule='minimum',gamma_cv_folds=3,gamma_cv_repeats=1,
   gamma_cv_bandwidth_multipliers=c(.25,.5,1,2,4),gamma_cv_lambda_grid=default_gamma_lambda_grid,
   gamma_nystrom_rank=NULL,allocation_spline_df=config$allocation_df,
   bootstrap_replications=B,seed=seed)
  if(!is.null(last_folds))stopifnot(identical(last_folds,f$fold_id))
  last_folds<<-f$fold_id
  f$settings$gamma_learner<-method
  f$settings$gamma_provider<-'Frozen x5_adaptive_provider; actual CV selections retained separately'
  selections<-lapply(as.list(attr(provider,'gamma_cache')),function(z){
   f<-z$fit;f[intersect(names(f),c('selected_lambda','selected_level','selected_bandwidth_multiplier',
    'base_bandwidth','selected_bandwidth','feature_count','rkhs_backend','cv_table','cv_resample_table',
    'training_negative_rate','normalization_constant','clip_rate','gradient'))]
  })
  list(summary=f$summary,gamma_diagnostics=f$gamma_diagnostics,gamma_selection=selections,
   probability_diagnostics=f$probability_diagnostics,settings=f$settings,
   witness_diagnostics=f$witness_diagnostics,fold_sizes=tabulate(f$fold_id,config$folds))
 }
 list(run=run,engine=engine)
}
run_task<-function(id,root,signature,smoke=FALSE,session_factory=new_session){
 stopifnot(id%in%tasks$id);t<-tasks[id,,drop=FALSE]
 path<-file.path(root,if(smoke)'smoke' else 'runs',sprintf('task_%04d.rds',id))
 lock<-paste0(path,'.lock');dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 if(!dir.create(lock,showWarnings=FALSE))stop('Lock exists; inspect processes: ',lock)
 on.exit(unlink(lock,recursive=TRUE))
 z<-load_task(path,t,signature)
 if(length(z$records)==length(config$methods)){cat('ALREADY COMPLETE',id,'; unchanged\n');return(invisible(z))}
 d<-x5g_generate(if(smoke)600L else config$n,t$a,t$b,t$seed,t$stress)
 session<-session_factory()
 for(m in config$methods){
  if(m%in%names(z$records))next
  warnings<-character();start<-proc.time()[3]
  v<-tryCatch(withCallingHandlers({
   f<-session$run(d,m,t$seed,if(smoke)99L else config$bootstrap)
   c(list(success=TRUE,error=''),f)
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
   error=function(e)list(success=FALSE,error=conditionMessage(e)))
  v$task<-t;v$signature<-signature;v$method<-m;v$n_used<-nrow(d)
  v$seconds<-unname(proc.time()[3]-start);v$warnings<-unique(warnings)
  check_record(v,t,m,signature);z$records[[m]]<-v;atomic(z,path)
  cat(format(Sys.time()),'task',id,t$type,'setting',t$setting,'rep',t$replication,m,
   'success',v$success,'seconds',v$seconds,v$error,'\n');flush.console()
 }
 cat('TASK ATTEMPTS COMPLETE',id,';',length(z$records),'learners saved\n');invisible(z)
}
