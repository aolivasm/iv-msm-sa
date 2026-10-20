# Exact oracle ratio for the ORIGINAL scalar DGP. No empirical renormalization.
soc_provider <- function() {
 e<-sls_engine('loglinear')
 e$sg_fit<-function(train,nuisance_train,t,type,seed) {
  predict<-function(newdata){
   g<-siv_gamma(newdata$Y-.65*newdata$X-.75*t,(2*t-1)*.15)
   list(raw=g,positive=g,normalized=g)
  }
  f<-list(predict=predict,normalization_constant=1,training_negative_rate=0,
   training_clip_rate=0,gradient=0,selected_lambda=NA_real_,
   gamma_estimator='oracle_original_scalar',rkhs_backend='oracle_original_scalar',
   cv_table=data.frame(),cv_resample_table=data.frame())
  for(k in c('selected_bandwidth_multiplier','base_bandwidth','selected_bandwidth',
   'nystrom_rank_requested','nystrom_landmark_count','nystrom_effective_rank',
   'nystrom_min_retained_eigenvalue','nystrom_max_eigenvalue'))f[[k]]<-NA_real_
  f
 }
 e$sg_provider('loglinear')
}
