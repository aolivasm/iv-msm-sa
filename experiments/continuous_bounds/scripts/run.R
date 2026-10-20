#!/usr/bin/env Rscript
source('scripts/common.R')
args<-commandArgs(TRUE);stage<-if(length(args))args[1] else 'preflight'
stopifnot(stage%in%c('preflight','task','check_gate','collect'))
manifest<-readRDS('manifest.rds');stopifnot(identical(manifest,tools::md5sum(names(manifest))))
source('R/x5_growing_gamma.R');source('R/x5_simple_thresholds.R');source('R/x5_adaptive_gamma_provider.R')
c<-stress_config();tasks<-stress_tasks();spec<-stress_spec();root<-'results/x5_stress_v1'
dir.create(root,recursive=TRUE,showWarnings=FALSE)
packages<-c('glmnet','ranger','xgboost','SuperLearner','nloptr')
stopifnot(all(vapply(packages,requireNamespace,logical(1),quietly=TRUE)))
target_object<-readRDS('results/population_targets.rds');targets<-target_object$targets
stopifnot(identical(target_object$config,c),nrow(targets)==16,
 all(targets$refinement[targets$compatible]<5e-4),all(is.na(targets$lower[!targets$compatible])))
signature<-list(config=c,manifest=manifest,targets_md5=unname(tools::md5sum('results/population_targets.rds')),
 R=as.character(getRversion()),packages=sapply(packages,function(p)as.character(packageVersion(p))))
meta<-file.path(root,'signature.rds')
if(stage=='preflight') {
 if(file.exists(meta))stopifnot(identical(readRDS(meta),signature)) else stress_atomic(signature,meta,FALSE)
 cat('PREFLIGHT PASSED: 2000 datasets, n=6000, 500 repetitions, 16000 attempted fits.\n',
 'Setting 4 / Lambda 2 is deliberately incompatible; coverage is NA.\n',
 'Frozen main-study learners, penalties, folds and seeds; one compact checkpoint per dataset.\n')
 quit(save='no')
}
stopifnot(identical(readRDS(meta),signature))
path_for<-function(id)file.path(root,sprintf('task_%04d.rds',id))
if(stage%in%c('check_gate','collect')) {
 ids<-if(stage=='check_gate')c(1L,501L,1001L,1501L) else 1:2000
 rows<-list();missing<-integer()
 for(id in ids) {
  t<-tasks[id,];path<-path_for(id)
  if(!file.exists(path)){missing<-c(missing,id);next}
  z<-stress_load(path,t,signature)
  if(!all(spec$key%in%names(z$records))){missing<-c(missing,id);next}
  for(k in spec$key) {
   v<-z$records[[k]]
   if(stage=='check_gate'&&v$compatible&&!v$success)stop('Compatible gate fit failed: task ',id,' ',k,' ',v$error)
   lo<-if(v$success)subset(v$ate,endpoint=='ate_lower') else NULL
   hi<-if(v$success)subset(v$ate,endpoint=='ate_upper') else NULL
   rows[[length(rows)+1]]<-data.frame(id=id,setting=t$setting,replication=t$replication,n=t$n,
    method=v$method,class=v$class,Lambda=v$Lambda,compatible=v$compatible,success=v$success,
    outward=v$outward,lower=if(v$success)lo$one_step_estimate else NA_real_,
    upper=if(v$success)hi$one_step_estimate else NA_real_,
    lower_centered=if(v$success)lo$centered_coverage else NA,
    upper_centered=if(v$success)hi$centered_coverage else NA,
    point_width=if(v$success)v$point_width else NA_real_,
    confidence_width=if(v$success)v$confidence_width else NA_real_,
    crossed_point=if(v$success)v$crossed_point else NA,
    gamma_projection_max=if(v$success)max(v$diagnostics$gamma_projection_rate_validation) else NA_real_,
    gamma_projection_magnitude=if(nrow(v$projection_events))max(v$projection_events$maximum_absolute_change) else NA_real_,
    warnings=length(v$warnings),error=v$error)
  }
 }
 cat('Complete datasets:',length(ids)-length(missing),'of',length(ids),'\n')
 if(length(missing)){print(missing);stop('Missing or incomplete tasks; not complete.')}
 d<-do.call(rbind,rows)
 print(aggregate(success~setting+method+class+Lambda,d,sum))
 if(stage=='collect') {
  groups<-split(d,interaction(d[c('setting','method','class','Lambda')],drop=TRUE))
  summaries<-do.call(rbind,lapply(groups,function(g) {
   stopifnot(nrow(g)==500);ok<-g$success
   cbind(g[1,c('setting','method','class','Lambda','compatible')],attempted=nrow(g),successful=sum(ok),
    failed=sum(!ok),coverage_among_successes=if(g$compatible[1]&&any(ok))mean(g$outward[ok]) else NA_real_,
    fraction_returning_covering_set=if(g$compatible[1])sum(g$outward[ok])/nrow(g) else NA_real_,
    mean_width_among_successes=if(any(ok))mean(g$confidence_width[ok]) else NA_real_)
  }))
  stress_atomic(list(detail=d,summary=summaries,signature=signature),file.path(root,'summary.rds'))
  print(summaries,row.names=FALSE)
 } else cat('GATE PASSED: all 32 attempts saved; compatible fits successful; incompatible coverage NA.\n')
 quit(save='no')
}
id<-as.integer(args[2]);stopifnot(length(id)==1L,!is.na(id),id%in%tasks$id)
t<-tasks[id,];a<-c$settings$a[t$setting];b<-c$settings$b[t$setting];path<-path_for(id)
lock<-sub('[.]rds$','.lock',path)
run<-function() {
 if(!dir.create(lock,showWarnings=FALSE))stop('Existing task lock: ',lock)
 on.exit(unlink(lock,recursive=TRUE))
 checkpoint<-stress_load(path,t,signature)
 if(all(spec$key%in%names(checkpoint$records))){cat('ALREADY COMPLETE:',id,'; no changes\n');return(invisible(NULL))}
 d<-x5g_generate(t$n,a,b,t$seed,TRUE);folds<-make_tz_stratified_folds(d,K=5,seed=t$seed)
 primitive_cache<-new.env(parent=emptyenv())
 primitive<-function(fold,train,x_cols,seed,...) {
  key<-as.character(fold)
  if(exists(key,primitive_cache,inherits=FALSE)) {
   v<-get(key,primitive_cache);stopifnot(identical(v$train,train),identical(v$seed,seed));return(v$fit)
  }
  f<-fit_nuisance_models(train,x_cols,positivity_epsilon=.05,learner='super_learner',
   treatment_spline_df=5,super_learner_candidates=sfp_contract()$primitive_library,
   super_learner_cv_folds=3,super_learner_ranger_num_trees=500,library_seed=seed)
  assign(key,list(train=train,seed=seed,fit=f),primitive_cache);f
 }
 # Record magnitudes without changing the projection function's returned value.
 original_projection<-project_gamma_to_compatible_interval;events<-list()
 monitored_projection<-function(gamma,nuisance,t_value,Lambda0,Lambda1,margin_fraction=1e-4) {
  v<-original_projection(gamma,nuisance,t_value,Lambda0,Lambda1,margin_fraction)
  change<-abs(v$gamma-v$original)
  events[[length(events)+1L]]<<-data.frame(treatment=t_value,n=length(change),
    projected=sum(v$projected),mean_absolute_change=mean(change),maximum_absolute_change=max(change))
  v
 }
 assign('project_gamma_to_compatible_interval',monitored_projection,envir=.GlobalEnv)
 on.exit(assign('project_gamma_to_compatible_interval',original_projection,envir=.GlobalEnv),add=TRUE)
 for(m in c$methods) {
  provider<-x5_adaptive_provider(m)
  for(L in c$Lambdas)for(cl in c$classes) {
   k<-paste0(m,'_',cl,'_L',L);if(k%in%names(checkpoint$records))next
   compatible<-stress_compatible(t$setting,L);events<-list();warnings<-character();qstore<-new.env(parent=emptyenv())
   start<-proc.time()[3]
   result<-tryCatch(withCallingHandlers({
    fit<-estimate_sharp_bounds_one_step(d,Lambda0_T1=L,Lambda1_T1=L,Lambda0_T0=L,Lambda1_T0=L,
      x_cols=paste0('X',1:5),K=5,fold_id=folds,nuisance_learner='super_learner',
      nuisance_fit_provider_factory=primitive,gamma_fit_provider_factory=provider,
      threshold_learner='external',bound_component_provider_factory=x5s_threshold_provider(cl,qstore),
      eif_conditional_learner='library_cv',
      eif_conditional_control=list(spline_df_grid=c(3,5,7,9),cv_folds=3,
        selection_rule='one_se',ranger_num_trees=500),
      estimated_compatibility_action='project',max_gamma_projection_rate=.10,
      compute_oracle_layers=FALSE,seed=t$seed)
    stopifnot(identical(folds,fit$fold_id),length(as.list(qstore))>0,
      all(vapply(as.list(qstore),function(v)v$lower$convergence==0&&v$upper$convergence==0,logical(1))))
    e<-fit$one_step
    for(i in seq_len(nrow(e))) {
     signal<-fit$one_step_signals[fit$one_step_signals$mode==e$mode[i],e$endpoint[i]]
     stopifnot(abs(mean(signal)-e$one_step_estimate[i])<1e-10,abs(sd(signal)/sqrt(t$n)-e$standard_error[i])<1e-10)
    }
    target<-targets[targets$setting==t$setting&targets$Lambda==L&targets$class==cl,]
    c(list(success=TRUE,certified=TRUE,error='',endpoints=e,diagnostics=fit$fold_diagnostics,
      eif_diagnostics=fit$eif_diagnostics),stress_measure(e,compatible,target))
   },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
    error=function(e)list(success=FALSE,certified=FALSE,error=conditionMessage(e),outward=NA))
   result$attempted<-TRUE;result$compatible<-compatible;result$task<-t;result$signature<-signature
   result$method<-m;result$class<-cl;result$Lambda<-L;result$a<-a;result$b<-b
   result$warnings<-unique(warnings);result$seconds<-unname(proc.time()[3]-start)
   result$projection_events<-if(length(events))do.call(rbind,events) else data.frame()
   result$projection_note<-'Per-call summaries include repeated training/validation evaluations; not unique-observation rates. Use fold diagnostics for rates.'
   stress_check_record(result,t,spec[spec$key==k,,drop=FALSE],signature)
   checkpoint$records[k]<-list(result)
   checkpoint$gamma_selection[[m]]<-lapply(as.list(attr(provider,'gamma_cache')),function(z){f<-z$fit
    list(fold=z$fold,treatment=z$treatment,lambda=f$selected_lambda,level=f$selected_level,
     dimension=f$feature_count,multiplier=f$selected_bandwidth_multiplier,bandwidth=f$selected_bandwidth,
     clip_rate=f$clip_rate,gradient=f$gradient,cv=f$cv_table,cv_resamples=f$cv_resample_table)})
   stress_atomic(checkpoint,path)
   cat(format(Sys.time()),'task',id,'setting',t$setting,'rep',t$replication,m,cl,'Lambda',L,
    'compatible',compatible,'success',result$success,'seconds',result$seconds,result$error,'\n');flush.console()
   if(exists('fit',inherits=FALSE))rm(fit);gc(FALSE)
  }
 }
 check<-stress_load(path,t,signature);stopifnot(all(spec$key%in%names(check$records)))
 cat('TASK ATTEMPTS COMPLETE',id,'; saved',length(check$records),'fits; successful',
  sum(vapply(check$records,function(z)isTRUE(z$success),logical(1))),'\n')
}
run()
