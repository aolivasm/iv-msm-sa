# Structural/numerical fixtures only. No estimator, simulation, or integration run.
settings <- dget('configurations/main_experiments.R')
for(p in settings$profiles) {
  if(p$engine=='continuous_bounds')
    stopifnot(identical(p$config$methods,'rkhs_tuned'),
              identical(p$config$classes,'fixed_spline5'))
  if(p$engine%in%c('binary_bounds','compatibility'))
    stopifnot(identical(p$config$methods,'rkhs'))
  if(grepl('comparators',p$engine))
    stopifnot(identical(p$config$adjustments,'x'))
}
load_modules <- function(path) {
  e<-new.env(parent=globalenv())
  for(p in list.files(path,pattern='[.]R$',full.names=TRUE))sys.source(p,e)
  e
}
e<-load_modules('experiments/continuous_bounds/R')
stopifnot(!exists('x5a_fit',e,inherits=FALSE),
          !exists('x5g_fit',e,inherits=FALSE),
          !exists('fit_threshold_growing_rbf_cv',e,inherits=FALSE))
# The exact fixed spline dictionary used in the paper.
d<-as.data.frame(matrix(seq(-1.9,1.9,length.out=100),ncol=5))
names(d)<-paste0('X',1:5)
expected<-do.call(cbind,lapply(d,function(x)
  splines::ns(x,knots=c(-1.2,-.4,.4,1.2),Boundary.knots=c(-2,2),intercept=FALSE)))
stopifnot(identical(e$restricted_x5_features(d,'fixed_spline5'),expected))
# Mock the gamma fit: verify the original tuning and caching without fitting.
calls<-0L
e$fit_gamma_rkhs<-function(train,nuisance_train,t_value,x_cols,...) {
  calls<<-calls+1L;list(options=list(...))
}
provider<-e$x5_adaptive_provider('rkhs_tuned')
train<-data.frame(Y=1:3,X1=c(-1,0,1),T_obs=c(0,1,0),Z=c(0,1,1))
first<-provider(1,0,train,list(e=rep(.5,3)),'X1',123L)
again<-provider(1,0,train,list(e=rep(.5,3)),'X1',123L)
stopifnot(calls==1L,identical(first,again),
          identical(first$options$bandwidth_multipliers,c(.25,.5,1,2,4)),
          first$options$cv_repeats==1,first$options$cv_folds==3,
          first$options$selection_rule=='minimum',first$options$backend=='nystrom')
# All binary point-bound primitives needed by the population-target script load.
b<-load_modules('experiments/binary_bounds/R')
stopifnot(is.function(b$binary_arm_endpoint),
          is.function(b$fit_binary_cross_fitted_direct_rkhs_nuisances),
          !exists('x5a_fit',b,inherits=FALSE))
cat('PASS: main-only scope, exact spline dictionary, RKHS tuning/cache, binary dependencies; no fits run.\n')
