source('R/x5_sparse_sieve_gamma.R')

# Exponential direct-Riesz loss, not a classification likelihood.
x5a_log_solve <- function(F,w0,w1,lambda){
 penalty<-c(0,rep(lambda,ncol(F)-1))
 objective<-function(beta){s<-drop(F%*%beta)
  mean(w0*exp(s)-w1*s)+sum(penalty*beta^2)/2}
 gradient<-function(beta)drop(crossprod(F,w0*exp(drop(F%*%beta))-w1))/nrow(F)+penalty*beta
 initial<-rep(0,ncol(F))
 fit<-optim(initial,objective,gradient,method='BFGS',control=list(maxit=1500,reltol=1e-11))
 g<-max(abs(gradient(fit$par)))
 if(fit$convergence!=0 || !is.finite(g) || g>1e-6){
  refined<-nlminb(fit$par,objective,gradient,control=list(iter.max=1500,eval.max=3000,rel.tol=1e-12))
  if(is.finite(refined$objective)&&refined$objective<=fit$value+1e-10){
   fit$par<-refined$par;fit$convergence<-refined$convergence
  }
  g<-max(abs(gradient(fit$par)))
 }
 # Deterministic Newton polishing with an Armijo line search.
 for(iter in seq_len(40)){
  grad<-gradient(fit$par);g<-max(abs(grad));if(g<1e-8)break
  H<-crossprod(F,F*(w0*exp(drop(F%*%fit$par))))/nrow(F)+diag(penalty)
  step<-drop(solve(H,grad));scale<-1;old<-objective(fit$par)
  while(scale>2^-30){
   candidate<-fit$par-scale*step;value<-objective(candidate)
   if(is.finite(value)&&value<=old-1e-4*scale*sum(grad*step))break
   scale<-scale/2
  }
  if(scale<=2^-30)break
  fit$par<-candidate
 }
 g<-max(abs(gradient(fit$par)))
 stopifnot(is.finite(g),g<1e-6,all(is.finite(fit$par)))
 raw<-exp(drop(F%*%fit$par));positive<-pmin(100,pmax(1e-6,raw))
 list(beta=fit$par,normalizer=mean(positive[w0==1]),gradient=g,
      convergence=fit$convergence,clip_rate=mean(raw<1e-6|raw>100))
}
x5a_log_predict<-function(F,fit){
 raw<-exp(drop(F%*%fit$beta));positive<-pmin(100,pmax(1e-6,raw))
 list(raw=raw,positive=positive,normalized=pmin(100,pmax(1e-6,positive/fit$normalizer)))
}
x5a_fit<-function(train,nu,t,seed,log=FALSE){
 if(!log)return(x5ss_fit(train,nu,t,seed))
 f<-x5ss_engine();e<-environment(f)
 e$x5g_solve<-x5a_log_solve;e$x5g_predict<-x5a_log_predict
 # Both candidates are selected using the same direct quadratic validation risk.
 fit<-f(train,nu,t,seed)
 fit$normalization_constant<-fit$normalizer
 fit$training_negative_rate<-0
 fit$rkhs_backend<-'growing_sparse_tensor_log_spline_riesz'
 for(k in c('selected_bandwidth_multiplier','base_bandwidth','selected_bandwidth',
 'nystrom_rank_requested','nystrom_landmark_count','nystrom_effective_rank',
 'nystrom_min_retained_eigenvalue','nystrom_max_eigenvalue'))fit[[k]]<-NA_real_
 fit
}
x5_adaptive_provider<-function(method){
 stopifnot(method %in% c('rkhs_tuned','growing_direct','growing_log'))
 cache<-new.env(parent=emptyenv())
 provider<-function(fold,treatment,train,nuisance_train,x_cols,seed,...){
  key<-paste(fold,treatment,sep='_')
  signature<-list(data=train[,c('Y',x_cols,'T_obs','Z')],nu=nuisance_train,seed=seed)
  if(exists(key,cache,inherits=FALSE)){
   z<-get(key,cache);stopifnot(identical(z$signature,signature));return(z$fit)
  }
  fit<-if(method=='rkhs_tuned')fit_gamma_rkhs(train,nuisance_train,treatment,x_cols,
   lambda_grid=default_gamma_lambda_grid,bandwidth_multipliers=c(.25,.5,1,2,4),
   cv_folds=3,cv_repeats=1,selection_rule='minimum',backend='nystrom',
   nystrom_rank=NULL,seed=seed) else
   x5a_fit(train,nuisance_train,treatment,seed,log=method=='growing_log')
  assign(key,list(signature=signature,fit=fit,fold=fold,treatment=treatment),cache)
  fit
 }
 attr(provider,'gamma_cache')<-cache;provider
}
