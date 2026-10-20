# Multivariable binary-outcome adaptation of the scalar direct Riesz spline.
# No classification model; no changes to production defaults.
# Source sharp_bounds_plugin.R and scalar_gamma_sieve_audit.R before use.
brs_features <- function(d, x_cols, spec=NULL) {
  stopifnot(all(d$Y %in% c(0,1)))
  vars <- c('Y',x_cols)
  if(is.null(spec)) {
    spec <- list(vars=vars,center=vapply(d[vars],mean,numeric(1)),
      scale=vapply(d[vars],sd,numeric(1)),basis=list(),active=vars)
    spec$active <- vars[is.finite(spec$scale) & spec$scale>1e-12]
    stopifnot('Y' %in% spec$active)
    for(v in spec$active) {
      z <- (d[[v]]-spec$center[v])/spec$scale[v]
      if(length(unique(z))>4L) spec$basis[[v]] <- splines::ns(z,df=3)
    }
  }
  blocks <- lapply(spec$active,function(v) {
    z <- (d[[v]]-spec$center[v])/spec$scale[v]
    if(is.null(spec$basis[[v]])) matrix(z,ncol=1) else predict(spec$basis[[v]],z)
  }); names(blocks)<-spec$active
  y<-blocks[['Y']]; bx<-do.call(cbind,blocks[setdiff(spec$active,'Y')])
  F<-cbind(1,y,bx,bx*drop(y))
  stopifnot(all(is.finite(F)))
  list(F=F,spec=spec)
}
brs_engine <- function(x_cols) {
  e<-new.env(parent=globalenv())
  sys.source('R/scalar_gamma_sieve_audit.R',envir=e)
  e$sg_features<-function(d,type,spec=NULL) {
    stopifnot(type=='riesz_spline')
    brs_features(d,x_cols,spec)
  }
  e
}

brs_calibrate <- function(scores, B=200000L, seed=20260919L) {
  # scores columns alternate lower, upper endpoints at each Lambda.
  sds<-apply(scores,2,sd);stopifnot(all(sds>0))
  S<-sweep(sweep(scores,2,colMeans(scores),'-'),2,rep(c(1,-1),ncol(scores)/2),'*')
  V<-cor(S)*(nrow(S)-1)/nrow(S)
  ee<-eigen(V,symmetric=TRUE);stopifnot(min(ee$values)>-1e-10)
  A<-ee$vectors%*%diag(sqrt(pmax(ee$values,0)))
  set.seed(seed)
  draws<-A%*%matrix(rnorm(B*ncol(S)),ncol=B)
  list(critical=unname(quantile(apply(draws,2,max),.95,type=8)),
       covariance=V,B=B,seed=seed)
}
