library(ivmsm)
# A small, deterministic teaching dataset, not a simulation or study result.
# Each row represents an independent person in this artificial illustration.
trial <- mock_iv_data(outcome='binary',iv_categories=3)
# Binary Y; IV Z=0,1,2; treatment T_obs; covariate X.
Lambda <- c('0'=2,'1'=2,'2'=2) # sensitivity by IV category

fit <- estimate_finite_iv_binary(
  trial, outcome='Y', treatment='T_obs', instrument='Z',
  covariates='X', reference='0', Lambda=Lambda, folds=2, seed=4,
  # Fast demonstration only. Omit fixed lambda/rank for production CV defaults.
  control=list(gamma_fixed_lambda=.1,gamma_rank=8,spline_df=2))
print(fit)
# Smoothed envelopes are not implemented for this general finite-IV estimator.
# No threshold class is needed for a binary outcome.
# Projection rates are estimation diagnostics, not compatibility-test results.
