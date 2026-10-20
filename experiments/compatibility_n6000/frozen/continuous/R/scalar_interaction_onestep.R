# Adapter only: preserve the existing interaction learner and one-step API.
source('R/scalar_interaction_penalty.R')
sio_provider <- function(method='log_interaction') {
 stopifnot(method=='log_interaction')
 e<-sls_engine('log_tensor')
 e$sg_fit<-function(train,nuisance_train,t,type,seed){
  f<-sip_fit(train,nuisance_train,t,seed)
  f$rkhs_backend<-'experimental_direct_log_interaction'
  for(k in c('selected_bandwidth_multiplier','base_bandwidth','selected_bandwidth',
   'nystrom_rank_requested','nystrom_landmark_count','nystrom_effective_rank',
   'nystrom_min_retained_eigenvalue','nystrom_max_eigenvalue')) f[[k]]<-NA_real_
  f
 }
 e$sg_provider('loglinear')
}
