# Binary/multivariable adaptation of sip_fit's direct interaction log-spline.
# Source binary_riesz_spline_application.R and sharp_bounds_plugin.R first.
bls_fit_fixed <- function(d,nu,t,x_cols,lambda,multiplier) {
  fs<-brs_features(d,x_cols);full<-fs$F
  nmain<-(ncol(full)+2)/2
  cols<-if(is.infinite(multiplier))seq_len(nmain) else seq_len(ncol(full))
  F<-full[,cols,drop=FALSE];n<-nrow(F)
  penalty<-c(0,rep(1,nmain-1),if(is.finite(multiplier))rep(multiplier,ncol(F)-nmain))
  w1<-as.numeric(d$Z==1)*compute_riesz_weight(nu,t)
  F0<-F[d$Z==0,,drop=FALSE];rhs<-drop(crossprod(F,w1))/n
  objective<-function(beta) {
    eta<-drop(F0%*%beta)
    if(any(eta>700))return(1e100)
    sum(exp(eta))/n-sum(rhs*beta)+lambda*sum(penalty*beta^2)/2
  }
  gradient<-function(beta)drop(crossprod(F0,exp(drop(F0%*%beta))))/n-rhs+lambda*penalty*beta
  opt<-optim(rep(0,ncol(F)),objective,gradient,method='BFGS',control=list(maxit=2000,reltol=1e-12))
  beta<-opt$par;gr<-max(abs(gradient(beta)));recovery<-'none'
  if(opt$convergence!=0 || !is.finite(gr) || gr>1e-4) {
    alt<-nlminb(beta,objective,gradient,control=list(iter.max=2000,eval.max=5000,rel.tol=1e-12))
    beta<-alt$par;gr<-max(abs(gradient(beta)));opt$convergence<-alt$convergence;recovery<-'nlminb'
  }
  if(opt$convergence!=0 || !is.finite(gr) || gr>1e-4)stop('Log-spline optimization failed certification.')
  spec<-fs$spec
  rawfun<-function(newdata)exp(drop(brs_features(newdata,x_cols,spec)$F[,cols,drop=FALSE]%*%beta))
  raw<-rawfun(d);stopifnot(all(is.finite(raw)))
  positive<-pmin(100,pmax(1e-6,raw));norm<-mean(positive[d$Z==0])
  # Keep only small objects in the prediction environment.
  predictor<-local({sp<-spec;cc<-cols;bb<-beta;xx<-x_cols;nn<-norm
    function(newdata) {
      r<-exp(drop(brs_features(newdata,xx,sp)$F[,cc,drop=FALSE]%*%bb))
      p<-pmin(100,pmax(1e-6,r))
      list(raw=r,positive=p,normalized=pmin(100,pmax(1e-6,p/nn)))
    }})
  list(predict=predictor,beta=beta,gradient=gr,convergence=opt$convergence,
    recovery=recovery,normalization_constant=norm,training_negative_rate=mean(raw<=0),
    training_clip_rate=mean(raw<1e-6 | raw>100),feature_spec=spec,
    selected_columns=cols,penalty_weights=penalty)
}
bls_fit <- function(train,nuisance_train,t,x_cols,seed) {
  i<-which(train$T_obs==t);d<-train[i,];nu<-lapply(nuisance_train,function(z)z[i])
  constants<-10^seq(-1.5,1,.5);multipliers<-c(1,10,100,Inf);rows<-list()
  for(rep in 1:2) {
    folds<-make_z_stratified_folds(d,K=3,seed=seed+100000L*(rep-1))
    for(k in 1:3)for(mult in multipliers)for(C in constants) {
      tr<-which(folds!=k);va<-which(folds==k);lam<-C/sqrt(nrow(d))
      f<-bls_fit_fixed(d[tr,],lapply(nu,function(z)z[tr]),t,x_cols,lam,mult)
      g<-f$predict(d[va,])$normalized
      w1<-as.numeric(d$Z[va]==1)*compute_riesz_weight(lapply(nu,function(z)z[va]),t)
      risk<-mean(as.numeric(d$Z[va]==0)*g-w1*log(g))
      rows[[length(rows)+1]]<-data.frame(repeat_id=rep,fold=k,C=C,multiplier=mult,
        lambda=lam,risk=risk,gradient=f$gradient,recovery=f$recovery)
    }
  }
  detail<-do.call(rbind,rows);cv<-aggregate(risk~C+multiplier+lambda,detail,mean)
  cv$selected<-seq_len(nrow(cv))==which.min(cv$risk);best<-cv[cv$selected,]
  f<-bls_fit_fixed(d,nu,t,x_cols,best$lambda,best$multiplier)
  f$selected_lambda<-best$lambda;f$selected_multiplier<-best$multiplier
  f$cv_table<-cv;f$cv_resample_table<-detail
  f$rkhs_backend<-'direct_binary_interaction_log_spline'
  for(k in c('selected_bandwidth_multiplier','base_bandwidth','selected_bandwidth',
    'nystrom_rank_requested','nystrom_landmark_count','nystrom_effective_rank',
    'nystrom_min_retained_eigenvalue','nystrom_max_eigenvalue'))f[[k]]<-NA_real_
  f
}
bls_engine <- function(x_cols) {
  e<-new.env(parent=globalenv())
  e$sg_fit<-function(train,nuisance_train,t,type,seed)bls_fit(train,nuisance_train,t,x_cols,seed)
  e
}
