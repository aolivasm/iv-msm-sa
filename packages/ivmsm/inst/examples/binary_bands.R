library(ivmsm)
# Binary Y and binary IV Z; treatment T_obs and covariate X.
# Simulated teaching data, not OHIE data or the paper's simulation design.
# Randomized encouragement Z affects treatment, but not the outcome directly.
# U affects both treatment and outcome and is withheld from the fitted analysis.
set.seed(72)
n <- 2000
X <- runif(n,-1,1)
U <- rnorm(n)
Z <- rbinom(n,1,.5)
T_obs <- rbinom(n,1,plogis(-.5+.8*Z+.5*X+.6*U))
Y <- rbinom(n,1,plogis(-1+.8*T_obs+.4*X+.5*U))
d <- data.frame(Y,T_obs,Z,X)
# Common sensitivity path Lambda0 = Lambda1 over the entire interval.
fit <- estimate_binary_sensitivity_band(d,x_cols='X',Lambda_range=c(1.2,2),
  grid_points=21,draws=1000,seed=7,K=2,
  primitive_learner='glm_spline',spline_df=2,
  rkhs_lambda_by_t=c('0'=.1,'1'=.1),rkhs_nystrom_rank=8)
print(fit)
plot(fit)
# Retrieve estimates and simultaneous confidence limits at calibration points:
grid_results <- fit$table[,c('Lambda','lower','upper',
  'simultaneous_lower','simultaneous_upper')]
# Example: the grid value nearest Lambda0 = Lambda1 = 1.6.
grid_results[which.min(abs(grid_results$Lambda-1.6)),]
# Values where the reported grid band excludes zero:
subset(grid_results,simultaneous_lower<=simultaneous_upper &
  (simultaneous_lower>0 | simultaneous_upper<0))
# First evaluated continuum location where the band includes zero:
head(subset(fit$continuum,lower<=0 & upper>=0),1)
# fit$continuum: continuum-valid limits evaluated at 201 locations, not linear
# interpolation. This uses the binary endpoint monotonicity/Lipschitz result.
# Nuisances are fitted once. This wrapper uses unsmoothed envelopes and a common
# Lambda, not arbitrary (Lambda0,Lambda1) pairs. Full checks: fit$diagnostics.
# Fixed penalties, a small rank/grid and 1,000 draws keep this example small.
# Use production tuning and a denser grid for an analysis. Valid bands require
# uniform nuisance-rate/stability conditions; projection does not certify them.
