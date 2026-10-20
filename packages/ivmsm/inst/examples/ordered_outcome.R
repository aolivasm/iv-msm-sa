library(ivmsm)
# A small fitted example. This is not a power/coverage study.
if (!requireNamespace('nnet',quietly=TRUE))
  stop("Install the example dependency first: install.packages('nnet')")
{
  # Ordered Y=0,1,2; binary IV Z=0,1; treatment T_obs=0,1; covariate X.
  d <- mock_iv_data(outcome='ordered',iv_categories=2)
  Lambda0 <- 4; Lambda1 <- 4 # sensitivity values for IV strata, not treatment arms
  nuisances <- fit_discrete_cross_fitted_nuisances(d,support=0:2,x_cols='X',
                   folds=2,df_grid=1,inner_folds=2)
  fit <- estimate_discrete_dual_bounds(d,nuisances,Lambda0=Lambda0,Lambda1=Lambda1,
    envelope_epsilon=.002,threshold_epsilon=.01)
  # Smoothing changes the evaluated bounds; retain these settings in reporting.
  print(fit)
}
