# Direct Riesz and exponential-tilt spline learners for (Y,X1,...,X5).
# Original scalar implementation remains unchanged.
source('R/scalar_factorial_pilot.R')
source('R/x5_simple_study.R')
x5s_gamma_engine<-function(method=c('spline3','log_interaction')){
 method<-match.arg(method);e<-new.env(parent=globalenv())
 sys.source('R/scalar_gamma_sieve_audit.R',e)
 e$sg_features<-function(d,type,spec=NULL){
  v<-as.matrix(d[,c('Y',paste0('X',1:5))])
  if(is.null(spec))spec<-list(center=colMeans(v),scale=apply(v,2,sd))
  stopifnot(all(is.finite(spec$scale)),all(spec$scale>0))
  v<-sweep(sweep(v,2,spec$center),2,spec$scale,'/')
  if(is.null(spec$basis))spec$basis<-lapply(1:6,function(j)splines::ns(v[,j],df=3))
  B<-lapply(1:6,function(j)predict(spec$basis[[j]],v[,j]))
  # 1 + 18 main spline terms + 45 outcome-by-covariate terms = 64.
  # This is not a full six-variable tensor or a universal fixed-size class.
  interactions<-lapply(2:6,function(j)do.call(cbind,lapply(1:3,function(k)B[[1]][,k]*B[[j]])))
  list(F=cbind(1,do.call(cbind,B),do.call(cbind,interactions)),spec=spec)
 }
 # Same fixed penalty grid and 2x3-fold CV as scalar direct learners.
 # Numerical recovery retains the original objective and certification rule.
 e$recovery<-list()
 e$optim<-function(par,fn,gr,...,method,control){
  z<-stats::optim(par,fn,gr,...,method=method,control=control)
  before<-max(abs(gr(z$par,...)))
  if(z$convergence==0 && is.finite(before) && before<=1e-4)return(z)
  r<-stats::nlminb(z$par,fn,gr,...,control=list(iter.max=3000,eval.max=6000,rel.tol=1e-14,x.tol=1e-12))
  mid<-list(code=r$convergence,gradient=max(abs(gr(r$par,...))))
  if(r$convergence!=0 || !is.finite(mid$gradient) || mid$gradient>1e-4){
   p<-stats::optim(r$par,fn,gr,...,method='BFGS',control=list(maxit=3000,reltol=1e-12))
   r<-list(par=p$par,objective=p$value,convergence=p$convergence)
  }
  after<-max(abs(gr(r$par,...)))
  e$recovery[[length(e$recovery)+1]]<-list(before=before,after=after,intermediate=mid,code=r$convergence)
  stopifnot(r$convergence==0,is.finite(after),after<=1e-4,fn(r$par,...)<=fn(z$par,...)+1e-12)
  z$par<-r$par;z$value<-r$objective;z$convergence<-r$convergence;z
 }
 e$type<-if(method=='spline3')'riesz_spline' else 'loglinear'
 e
}
x5s_gamma_provider<-function(method=c('spline3','log_interaction')){
 method<-match.arg(method);e<-x5s_gamma_engine(method);cache<-new.env(parent=emptyenv())
 provider<-function(fold,treatment,train,nuisance_train,x_cols,seed,...){
  stopifnot(identical(x_cols,paste0('X',1:5)))
  key<-paste(fold,treatment,sep='_')
  signature<-list(data=train[,c('Y',x_cols,'T_obs','Z')],nuisance=nuisance_train,seed=seed)
  if(exists(key,cache,inherits=FALSE)){
   old<-get(key,cache);stopifnot(identical(old$signature,signature));return(old$fit)
  }
  e$recovery<-list();f<-e$sg_fit(train,nuisance_train,treatment,e$type,seed)
  f$numerical_recovery<-e$recovery;f$learner<-method;f$feature_count<-64L
  assign(key,list(signature=signature,fit=f,fold=fold,treatment=treatment),cache);f
 }
 attr(provider,'gamma_cache')<-cache;provider
}
