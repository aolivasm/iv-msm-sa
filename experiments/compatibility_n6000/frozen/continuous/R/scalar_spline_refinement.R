# Isolated gamma refinement: threshold class and production defaults unchanged.
source('R/scalar_gamma_sieve_audit.R')
ssr_spec <- function(method) switch(method,
  spline3=list(df=3L,constants=10^seq(-1.5,1,.5)),
  spline3_wide=list(df=3L,constants=10^seq(-4,1,.5)),
  spline5_wide=list(df=5L,constants=10^seq(-4,1,.5)),
  loglinear=list(df=3L,constants=10^seq(-1.5,1,.5)),stop('Unknown method'))
ssr_engine <- function(method) {
  config<-ssr_spec(method);e<-new.env(parent=globalenv())
  sys.source('R/scalar_gamma_sieve_audit.R',e)
  e$ssr_config<-config
  e$sg_features<-function(d,type,spec=NULL) {
    z<-as.matrix(d[,c('Y','X')])
    if(is.null(spec)) spec<-list(center=colMeans(z),scale=apply(z,2,sd))
    stopifnot(all(is.finite(spec$scale)),all(spec$scale>0))
    z<-sweep(sweep(z,2,spec$center),2,spec$scale,'/')
    if(type=='loglinear')return(list(F=cbind(1,z),spec=spec))
    if(is.null(spec$basis))spec$basis<-lapply(1:2,function(j)splines::ns(z[,j],df=config$df))
    B<-lapply(1:2,function(j)predict(spec$basis[[j]],z[,j]))
    F<-cbind(1,B[[1]],B[[2]],do.call(cbind,lapply(seq_len(ncol(B[[1]])),function(j)B[[1]][,j]*B[[2]])))
    list(F=F,spec=spec)
  }
  # Replace the one prespecified grid assignment; retain the existing CV rule.
  b<-body(e$sg_fit);hits<-which(vapply(as.list(b),function(x)
    is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name('constants')),logical(1)))
  stopifnot(length(hits)==1);b[[hits]]<-quote(constants<-ssr_config$constants);body(e$sg_fit)<-b
  e
}
ssr_provider <- function(method) {
  e<-ssr_engine(method);e$sg_provider(if(method=='loglinear')'loglinear' else 'riesz_spline')
}
