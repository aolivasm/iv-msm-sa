library(ivmsm)
# Binary outcome Y, binary IV Z=0/1, binary treatment T_obs, covariate X.
d <- mock_iv_data(outcome='binary',iv_categories=2)
head(d)
# Lambda0 and Lambda1 index IV strata Z=0 and Z=1, NOT treatment arms.
Lambda0 <- 2
Lambda1 <- 2
# epsilon=0 uses unsmoothed envelopes; a positive value widens the envelope
# relaxation. Smoothing changes the evaluated bound, not just computation.
epsilon <- 0
# Fixed penalties and rank eight are for a fast demonstration only.
# Production defaults retain repeated CV and should be used for real analyses.
fit <- estimate_binary_dual_bounds(d,x_cols='X',Lambda0=Lambda0,Lambda1=Lambda1,
  epsilon=epsilon,K=2,
  primitive_learner='glm_spline',spline_df=2,
  rkhs_lambda_by_t=c('0'=.1,'1'=.1),rkhs_nystrom_rank=8)
print(fit)
# summary(fit) gives a compact structured summary; fit$diagnostics keeps details.
