# Experimental refit only: CV is unchanged, final lambda is multiplied by .1.
source('R/scalar_spline_refinement.R')
sgu_provider<-function(factor=.1){
 stopifnot(is.finite(factor),factor>0,factor<=1)
 e<-ssr_engine('spline3');original<-e$sg_fit
 e$sg_fit<-function(train,nuisance_train,t,type,seed){
  base<-original(train,nuisance_train,t,type,seed)
  i<-which(train$T_obs==t)
  fit<-e$sg_fit_fixed(train[i,],lapply(nuisance_train,function(z)z[i]),t,type,
    factor*base$selected_lambda)
  fit$cv_selected_lambda<-base$selected_lambda
  fit$selected_lambda<-factor*base$selected_lambda
  fit$undersmoothing_factor<-factor
  fit$cv_table<-base$cv_table;fit$cv_resample_table<-base$cv_resample_table
  fit$rkhs_backend<-'experimental_spline3_undersmoothed'
  for(k in c('selected_bandwidth_multiplier','base_bandwidth','selected_bandwidth',
    'nystrom_rank_requested','nystrom_landmark_count','nystrom_effective_rank',
    'nystrom_min_retained_eigenvalue','nystrom_max_eigenvalue'))fit[[k]]<-NA_real_
  stopifnot(abs(fit$selected_lambda/base$selected_lambda-factor)<1e-12)
  fit
 }
 e$sg_provider('riesz_spline')
}
