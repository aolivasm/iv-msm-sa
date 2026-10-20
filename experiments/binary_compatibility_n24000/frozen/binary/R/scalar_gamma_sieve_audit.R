# Experimental direct density-ratio learners. Production defaults unchanged.
# Source after R/scalar_factorial_pilot.R. Scalar X only.
sg_features <- function(d, type, spec=NULL) {
  z <- as.matrix(d[,c('Y','X')])
  if (is.null(spec)) spec <- list(center=colMeans(z),scale=apply(z,2,sd))
  stopifnot(all(is.finite(spec$scale)),all(spec$scale>0))
  z <- sweep(sweep(z,2,spec$center),2,spec$scale,'/')
  if(type=='loglinear') return(list(F=cbind(1,z),spec=spec))
  if(is.null(spec$basis)) spec$basis <- lapply(1:2,function(j) splines::ns(z[,j],df=3))
  B <- lapply(1:2,function(j) predict(spec$basis[[j]],z[,j]))
  # Modest tensor spline; gamma class, NOT the fixed threshold class.
  F <- cbind(1,B[[1]],B[[2]],do.call(cbind,lapply(1:3,function(j) B[[1]][,j]*B[[2]])))
  list(F=F,spec=spec)
}
sg_objective <- function(beta,F,w0,w1,lambda,type) {
  eta <- drop(F%*%beta)
  if(type=='loglinear') mean(w0*exp(eta)-w1*eta)+lambda*sum(beta[-1]^2)/2
  else mean(w0*(1+eta)^2/2-w1*(1+eta))+lambda*sum(beta[-1]^2)/2
}
sg_gradient <- function(beta,F,w0,w1,lambda,type) {
  eta <- drop(F%*%beta)
  g <- if(type=='loglinear') exp(eta) else 1+eta
  drop(crossprod(F,w0*g-w1))/nrow(F)+lambda*c(0,beta[-1])
}
sg_fit_fixed <- function(d,nu,t,type,lambda) {
  features <- sg_features(d,type); F <- features$F
  w0 <- as.numeric(d$Z==0); w1 <- as.numeric(d$Z==1)*compute_riesz_weight(nu,t)
  if(type=='riesz_spline') {
    P <- diag(c(0,rep(lambda,ncol(F)-1)))
    beta <- solve(crossprod(F,F*w0)/nrow(F)+P, crossprod(F,w1-w0)/nrow(F))[,1]
    convergence <- 0L
  } else {
    opt <- optim(rep(0,ncol(F)),sg_objective,sg_gradient,F=F,w0=w0,w1=w1,
      lambda=lambda,type=type,method='BFGS',control=list(maxit=2000,reltol=1e-12))
    beta <- opt$par; convergence <- opt$convergence
  }
  gradient <- max(abs(sg_gradient(beta,F,w0,w1,lambda,type)))
  if(convergence!=0 || !is.finite(gradient) || gradient>1e-4) stop('Gamma optimization failed certification.')
  raw_predict <- function(newdata) {
    eta <- drop(sg_features(newdata,type,features$spec)$F%*%beta)
    if(type=='loglinear') exp(eta) else 1+eta
  }
  raw <- raw_predict(d)
  positive <- clip_positive(raw,floor=1e-6,cap=100)
  norm <- mean(positive[d$Z==0])
  predict <- function(newdata) {
    r <- raw_predict(newdata); p <- clip_positive(r,floor=1e-6,cap=100)
    list(raw=r,positive=p,normalized=clip_positive(p/norm,floor=1e-6,cap=100))
  }
  list(predict=predict,normalization_constant=norm,training_negative_rate=mean(raw<=0),
    training_clip_rate=mean(raw<1e-6 | raw>100),gradient=gradient,beta=beta)
}
sg_fit <- function(train,nuisance_train,t,type,seed) {
  i <- which(train$T_obs==t); d <- train[i,]; nu <- lapply(nuisance_train,function(z) z[i])
  constants <- 10^seq(-1.5,1,.5)
  rows <- list()
  # Same 2x3 stratified partitions as RKHS; transformations learned inside splits.
  for(rep in 1:2) {
    folds <- make_z_stratified_folds(d,K=3,seed=seed+100000L*(rep-1))
    for(k in 1:3) for(C in constants) {
      tr <- which(folds!=k); va <- which(folds==k)
      lambda <- C/sqrt(nrow(d)) # Match current RKHS outer-arm grid convention.
      f <- sg_fit_fixed(d[tr,],lapply(nu,function(z) z[tr]),t,type,lambda)
      g <- f$predict(d[va,])$normalized
      w0 <- as.numeric(d$Z[va]==0)
      w1 <- as.numeric(d$Z[va]==1)*compute_riesz_weight(lapply(nu,function(z) z[va]),t)
      risk <- if(type=='loglinear') mean(w0*g-w1*log(g)) else mean(w0*g^2/2-w1*g)
      rows[[length(rows)+1]] <- data.frame(repeat_id=rep,fold=k,C=C,lambda=lambda,risk=risk)
    }
  }
  detail <- do.call(rbind,rows)
  cv <- aggregate(risk~C+lambda,detail,mean)
  cv$selected <- seq_len(nrow(cv))==which.min(cv$risk)
  fit <- sg_fit_fixed(d,nu,t,type,cv$lambda[cv$selected])
  fit$selected_lambda <- cv$lambda[cv$selected]
  fit$cv_table <- cv; fit$cv_resample_table <- detail
  fit$rkhs_backend <- paste0('experimental_direct_',type)
  for(k in c('selected_bandwidth_multiplier','base_bandwidth','selected_bandwidth',
    'nystrom_rank_requested','nystrom_landmark_count','nystrom_effective_rank',
    'nystrom_min_retained_eigenvalue','nystrom_max_eigenvalue')) fit[[k]] <- NA_real_
  fit
}
sg_provider <- function(type) {
  cache <- new.env(parent=emptyenv())
  provider <- function(fold,treatment,train,nuisance_train,x_cols,seed,...) {
    stopifnot(identical(x_cols,'X'))
    key <- paste(fold,treatment,sep='_')
    signature <- list(data=train[,c('Y','X','T_obs','Z')],nuisance=nuisance_train)
    if(exists(key,cache,inherits=FALSE)) {
      entry <- get(key,cache);stopifnot(identical(signature,entry$signature));return(entry$fit)
    }
    fit <- sg_fit(train,nuisance_train,treatment,type,seed)
    assign(key,list(signature=signature,fit=fit,fold=fold,treatment=treatment),cache)
    fit
  }
  attr(provider,'gamma_cache') <- cache
  provider
}
