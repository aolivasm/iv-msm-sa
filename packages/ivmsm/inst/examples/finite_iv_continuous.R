library(ivmsm)
trial <- mock_iv_data(outcome='continuous',iv_categories=3)
# Continuous Y; IV Z=0,1,2; treatment T_obs; covariate X.
Lambda <- c('0'=2,'1'=2,'2'=2)

# Prespecified spline class. Supply knots/boundaries appropriate to your scale;
# do not learn knots from the validation outcomes or call ns(df=...) here.
basis <- function(x) cbind(1,splines::ns(x$X,knots=0,Boundary.knots=c(-1,1)))
fit <- estimate_finite_iv_continuous(
  trial,outcome='Y',treatment='T_obs',instrument='Z',
  covariates='X',reference='0',Lambda=Lambda,folds=2,
  threshold_basis=basis,seed=4,
  # Fast demonstration only; omit these controls for repeated-CV gamma tuning.
  control=list(gamma_fixed_lambda=.1,gamma_rank=8,spline_df=2))
print(fit)
# This general finite-IV estimator uses unsmoothed envelopes.
# Warnings/projection diagnostics are not formal compatibility tests.
