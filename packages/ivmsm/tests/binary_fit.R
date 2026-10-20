library(ivmsm)
set.seed(8)
d <- data.frame(Y=rbinom(240,1,.5),T_obs=rep(0:1,120),Z=rep(c(0,0,1,1),60),X=rnorm(240))
f <- estimate_binary_dual_bounds(d,x_cols='X',K=2,
  primitive_learner='glm_spline',spline_df=2,
  rkhs_lambda_by_t=c('0'=.1,'1'=.1),rkhs_nystrom_rank=8)
stopifnot(is.list(f$nuisance_fit),is.data.frame(f$diagnostics),
          all(is.finite(as.matrix(f$ate_confidence_set))))
