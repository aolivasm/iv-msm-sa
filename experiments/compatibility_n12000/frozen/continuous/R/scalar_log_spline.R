# Experimental direct exponential-tilt fits, not treatment classification.
source('R/scalar_spline_refinement.R')
sls_engine <- function(method) {
  stopifnot(method %in% c('spline3','loglinear','log_additive','log_tensor'))
  e <- ssr_engine(if(method=='spline3') 'spline3' else 'loglinear')
  if(method %in% c('spline3','loglinear')) return(e)
  tensor <- method=='log_tensor'
  e$sg_features <- function(d,type,spec=NULL) {
    z <- as.matrix(d[,c('Y','X')])
    if(is.null(spec)) spec <- list(center=colMeans(z),scale=apply(z,2,sd))
    stopifnot(all(is.finite(spec$scale)),all(spec$scale>0))
    z <- sweep(sweep(z,2,spec$center),2,spec$scale,'/')
    if(is.null(spec$basis)) spec$basis <- lapply(1:2,function(j)splines::ns(z[,j],df=3))
    B <- lapply(1:2,function(j)predict(spec$basis[[j]],z[,j]))
    F <- cbind(1,B[[1]],B[[2]])
    if(tensor) F <- cbind(F,do.call(cbind,lapply(1:3,function(j)B[[1]][,j]*B[[2]])))
    list(F=F,spec=spec)
  }
  e
}
sls_provider <- function(method) {
  e <- sls_engine(method)
  e$sg_provider(if(method=='spline3') 'riesz_spline' else 'loglinear')
}

# Nonlinear stress DGP used FIRST for gamma-only validation; no inferential
# targets from the linear-tilt DGP are recycled for this law.
# f1(r|x,t)/f0(r) = exp{s_t*b*h(r,x)}/M_t(x), r in [-2,2].
# h(r,x) = {r + .5*(r^2-E0[R^2]) + .5*x*r}/3.
# Division by three keeps the nonlinear score within [-2,2].
sls_h <- function(r,x) (r+.5*(r^2-(1-4*dnorm(2)/(pnorm(2)-pnorm(-2))))+.5*x*r)/3
sls_normalizer <- function(x,t,b=.15) {
  rule<-siv_x_rule(60)
  vapply(x,function(v)sum(rule$w*exp((2*t-1)*b*sls_h(rule$x,v))),numeric(1))
}
sls_stress_generate <- function(n,seed,b=.15,a=1.1) {
  set.seed(seed)
  x<-qnorm(runif(n,pnorm(-2),pnorm(2)));z<-rbinom(n,1,.5)
  p0<-plogis(.5*x-a/2);p1<-plogis(.5*x+a/2)
  t<-rbinom(n,1,ifelse(z==0,p0,p1))
  r<-qnorm(runif(n,pnorm(-2),pnorm(2)))
  # Conservative rejection envelope |h| <= 2.
  pending<-which(z==1)
  while(length(pending)) {
    proposal<-qnorm(runif(length(pending),pnorm(-2),pnorm(2)))
    accept<-log(runif(length(pending))) <= (2*t[pending]-1)*b*sls_h(proposal,x[pending])-2*b
    r[pending[accept]]<-proposal[accept];pending<-pending[!accept]
  }
  y<-.65*x+.75*t+r
  # Potential-arm residuals must use that arm's location, even off-arm.
  g<-lapply(0:1,function(arm)exp((2*arm-1)*b*sls_h(y-.65*x-.75*arm,x))/sls_normalizer(x,arm,b))
  data.frame(Y=y,X=x,Z=z,T_obs=t,e_true=.5,pi0_true=p0,pi1_true=p1,
    gamma0_oracle=g[[1]],gamma1_oracle=g[[2]])
}
