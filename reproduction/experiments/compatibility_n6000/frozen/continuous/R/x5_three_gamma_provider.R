source('R/x5_paired_gamma_pilot.R')
x5_three_provider <- function(method) {
 if(method=='spline3')return(x5s_gamma_provider('spline3'))
 stopifnot(method %in% c('rkhs','log_interaction'))
 cache<-new.env(parent=emptyenv())
 provider<-function(fold,treatment,train,nuisance_train,x_cols,seed,...){
  key<-paste(fold,treatment,sep='_')
  signature<-list(data=train[,c('Y',x_cols,'T_obs','Z')],nuisance=nuisance_train,seed=seed)
  if(exists(key,cache,inherits=FALSE)){
   old<-get(key,cache);stopifnot(identical(old$signature,signature));return(old$fit)
  }
  fit<-if(method=='rkhs')fit_gamma_rkhs(train,nuisance_train,treatment,x_cols,
    lambda_grid=default_gamma_lambda_grid,bandwidth_multipliers=1/sqrt(2),
    cv_folds=3,cv_repeats=2,selection_rule='minimum',backend='nystrom',
    nystrom_rank=ceiling(10*sqrt(sum(train$T_obs==treatment))),seed=seed)
   else x5p_log_fit(train,nuisance_train,treatment,seed,smoke=FALSE)
  if(method=='log_interaction') {
   # The one-step API expects kernel diagnostics even for non-kernel fits.
   # These are inapplicable metadata, not changes to the fit or its tuning.
   fit$rkhs_backend <- 'experimental_direct_log_interaction'
   for(k in c('selected_bandwidth_multiplier','base_bandwidth','selected_bandwidth',
     'nystrom_rank_requested','nystrom_landmark_count','nystrom_effective_rank',
     'nystrom_min_retained_eigenvalue','nystrom_max_eigenvalue'))fit[[k]]<-NA_real_
  }
  assign(key,list(signature=signature,fit=fit,fold=fold,treatment=treatment),cache)
  fit
 }
 attr(provider,'gamma_cache')<-cache;provider
}
