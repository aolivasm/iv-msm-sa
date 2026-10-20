source('R/observed_binary_dgp.R')
source('R/sharp_bounds_plugin.R')
source('R/binary_rkhs_gamma_nuisances.R')
source('R/x5_adaptive_gamma_provider.R')
config<-list(n=6000L,reps=500L,Lambda=c(1.8,2,2.5,3),methods=c('rkhs','growing_direct'),
 settings=data.frame(a=c(.55,1.1,1.65,1.1),b=c(.075,.075,.075,.15)),K=5L,epsilon=0,
 primitives=c('glm_main','glm_interactions','glmnet_interactions','ranger','xgboost'))
tasks<-data.frame(id=1:2000,setting=rep(1:4,each=500),replication=rep(1:500,4))
tasks$seed<-261000L+tasks$setting+100L*(tasks$replication-1L)
generate<-function(t,n=config$n){
 set.seed(t$seed);x<-matrix(qnorm(pnorm(-2)+runif(n*5)*(pnorm(2)-pnorm(-2))),ncol=5)
 colnames(x)<-paste0('X',1:5);z<-rbinom(n,1,.5);g<-config$settings[t$setting,]
 h<-drop(x%*%c(.14,.12,.10,.08,.06));tr<-rbinom(n,1,plogis(h+g$a*(z-.5)))
 shape<-1+.6*sin(x[,1]*x[,2])+.4*sin(x[,3]*x[,4])
 p<-plogis(-.4+1.5*tr+drop(x%*%c(.20,-.25,.30,.15,-.20))+(2*tr-1)*g$b*shape*z)
 data.frame(Y=rbinom(n,1,p),T_obs=tr,Z=z,x)
}
atomic<-function(z,p){
 dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);tmp<-tempfile('unpublished_',dirname(p))
 on.exit(if(file.exists(tmp))unlink(tmp))
 withCallingHandlers(saveRDS(z,tmp),warning=function(w)stop(conditionMessage(w)))
 stopifnot(identical(z,readRDS(tmp)));stopifnot(file.rename(tmp,p),identical(z,readRDS(p)))
}
check<-function(z,t,sig){
 stopifnot(identical(z$task,t),identical(z$signature,sig),is.list(z$records),
 !anyDuplicated(names(z$records)),all(names(z$records)%in%config$methods))
 for(m in names(z$records)){
  v<-z$records[[m]];stopifnot(is.logical(v$success),length(v$success)==1,!is.na(v$success))
  if(v$success){stopifnot(length(v$fits)==4L)
   for(i in 1:4){e<-v$fits[[i]]$endpoint_intervals
    stopifnot(v$fits[[i]]$Lambda0==config$Lambda[i],all(is.finite(e$estimate)),
     all(is.finite(e$standard_error)),all(e$standard_error>=0),all(c('ate_lower','ate_upper')%in%trimws(e$endpoint)))}
  }else stopifnot(nzchar(v$error))
 }
 invisible(z)
}
run_task<-function(id,root,sig,n=config$n){
 stopifnot(id%in%tasks$id);t<-tasks[id,,drop=FALSE];p<-file.path(root,sprintf('task_%04d.rds',id))
 dir.create(root,recursive=TRUE,showWarnings=FALSE);lock<-paste0(p,'.lock')
 if(!dir.create(lock,showWarnings=FALSE))stop('Lock exists: ',lock)
 on.exit(unlink(lock,recursive=TRUE))
 z<-if(file.exists(p))readRDS(p) else list(task=t,signature=sig,records=list())
 check(z,t,sig);if(length(z$records)==2L){cat('ALREADY COMPLETE',id,'\n');return(invisible(z))}
 d<-generate(t,n);folds<-binary_make_tz_folds(d,5,t$seed)
 cache<-new.env(parent=emptyenv());original<-binary_rkhs_gamma_fit_component
 previous<-NULL
 for(m in config$methods){
  if(m%in%names(z$records))next
  env<-new.env(parent=globalenv());fun<-fit_binary_cross_fitted_direct_rkhs_nuisances;environment(fun)<-env
  env$binary_rkhs_gamma_fit_component<-function(...){
   a<-list(...);key<-as.character(a$seed)
   if(exists(key,cache,inherits=FALSE)){v<-get(key,cache);stopifnot(identical(a,v$args));return(v$fit)}
   f<-do.call(original,a);assign(key,list(args=a,fit=f),cache);f
  }
  selections<-list()
  env$fit_gamma_rkhs<-function(train,nuisance_train,t_value,x_cols,seed,...){
   f<-if(m=='rkhs')fit_gamma_rkhs(train,nuisance_train,t_value,x_cols,seed=seed,...) else {
    f<-x5a_fit(train,nuisance_train,t_value,seed,log=FALSE);f$rkhs_kernel<-NA_character_;f}
   selections[[length(selections)+1L]]<<-f[intersect(names(f),c('selected_lambda','selected_level','feature_count',
    'base_bandwidth','selected_bandwidth','selected_bandwidth_multiplier','cv_table','cv_resample_table','rkhs_backend'))]
   f
  }
  w<-character();start<-proc.time()[3]
  v<-tryCatch(withCallingHandlers({
   nf<-fun(d,paste0('X',1:5),K=5,fold_id=folds,seed=t$seed,primitive_learner='super_learner',
    super_learner_candidates=config$primitives,rkhs_backend='nystrom',rkhs_bandwidth_multipliers=1/sqrt(2),
    rkhs_cv_folds=3,rkhs_cv_repeats=2,rkhs_selection_rule='minimum')
   stopifnot(identical(nf$fold_id,folds))
   primitive_columns<-c('e','pi0','pi1','a_t0_z0','a_t1_z0')
   if(!is.null(previous))stopifnot(identical(previous,nf$predictions[,primitive_columns]))
   previous<-nf$predictions[,primitive_columns]
   nf$configuration$gamma_method<-m
   fits<-lapply(config$Lambda,function(L){
    f<-evaluate_binary_bounds_from_nuisances(d,nf,Lambda0=L,Lambda1=L,epsilon=0,compatibility_action='project')
    f[c('endpoint_intervals','ate_confidence_set','diagnostics','Lambda0','Lambda1','epsilon')]
   })
   list(success=TRUE,error='',fits=fits,nuisance=nf[c('diagnostics','cv_table','gamma_cv_table',
    'gamma_cv_resample_table','gamma_diagnostics','reconstruction_diagnostics','configuration')],gamma_selection=selections)
  },warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),
   error=function(e)list(success=FALSE,error=conditionMessage(e)))
  v$warnings<-unique(w);v$seconds<-unname(proc.time()[3]-start);v$n<-n
  z$records[[m]]<-v;check(z,t,sig);atomic(z,p)
  cat(format(Sys.time()),'task',id,'setting',t$setting,'rep',t$replication,m,'success',v$success,'seconds',v$seconds,v$error,'\n');flush.console()
 }
 invisible(z)
}
