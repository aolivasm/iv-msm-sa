source('R/x5_growing_gamma.R')
x5g_provider<-function(){
 cache<-new.env(parent=emptyenv())
 provider<-function(fold,treatment,train,nuisance_train,x_cols,seed,...){
  stopifnot(identical(x_cols,paste0('X',1:5)))
  key<-paste(fold,treatment,sep='_');signature<-list(data=train[,c('Y',x_cols,'T_obs','Z')],nu=nuisance_train,seed=seed)
  if(exists(key,cache,inherits=FALSE)){old<-get(key,cache);stopifnot(identical(old$signature,signature));return(old$fit)}
  fit<-x5g_fit(train,nuisance_train,treatment,seed)
  fit$rkhs_backend<-'experimental_growing_riesz';fit$training_clip_rate<-fit$clip_rate
  fit$normalization_constant<-fit$normalizer
  fit$training_negative_rate<-mean(fit$predict(train[train$T_obs==treatment,])$raw<=0)
  for(k in c('selected_bandwidth_multiplier','base_bandwidth','selected_bandwidth','nystrom_rank_requested','nystrom_landmark_count','nystrom_effective_rank','nystrom_min_retained_eigenvalue','nystrom_max_eigenvalue'))fit[[k]]<-NA_real_
  assign(key,list(signature=signature,fit=fit,fold=fold,treatment=treatment),cache);fit
 }
 attr(provider,'gamma_cache')<-cache;provider
}
