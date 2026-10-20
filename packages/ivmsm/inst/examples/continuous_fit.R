library(ivmsm)
# Continuous Y; binary IV Z; treatment T_obs; covariate X.
d <- mock_iv_data(outcome='continuous',iv_categories=2)
Lambda0 <- 2; Lambda1 <- 2 # sensitivity by IV stratum
f <- estimate_finite_iv_continuous(d,'Y','T_obs','Z','X',
  Lambda=c('0'=Lambda0,'1'=Lambda1),reference='0',folds=2,
  threshold_class='spline',
  spline_spec=list(X=list(knots=0,boundary=c(-1,1))),seed=4,
  control=list(gamma_fixed_lambda=.1,gamma_rank=8,spline_df=2))
print(f)
# Fixed spline class; unsmoothed envelopes. Small rank/fixed penalty are for
# demonstration. The same interface accepts more than two IV categories.
