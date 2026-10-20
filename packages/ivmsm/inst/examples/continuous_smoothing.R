library(ivmsm)
d <- mock_iv_data('continuous',2) # Y continuous; Z binary; T_obs binary; X covariate
Lambda0 <- 2; Lambda1 <- 2
fit <- estimate_sharp_bounds_one_step(d,x_cols='X',K=2,
  Lambda0_T0=Lambda0,Lambda1_T0=Lambda1,
  Lambda0_T1=Lambda0,Lambda1_T1=Lambda1,
  envelope_smoothing_epsilon=.01, # set to 0 for unsmoothed envelopes
  compute_oracle_layers=FALSE,nuisance_learner='glm_spline',
  treatment_spline_df=2,conditional_mean_spline_df=2,threshold_spline_df=2,
  threshold_learner='spline_fixed',rkhs_lambda_by_t=c('0'=.1,'1'=.1),
  gamma_rkhs_backend='nystrom',gamma_nystrom_rank=8)
print(fit)
# Lambda0/Lambda1 index Z strata; T0/T1 select the potential-outcome arm.
# Smoothing relaxes the envelopes and changes the evaluated bounds.
# Fixed penalties and small rank are for teaching, not production fitting.
