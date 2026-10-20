# Extension code

Follow the [download and installation instructions](../README.md#getting-started),
then load `library(ivmsm)` in the R console.
The functions below are exported unless explicitly identified as historical
experiment scripts. Loading the package launches no simulations.

## Discrete ordered outcomes

```r
library(ivmsm)
```

`discrete_arm_endpoint()` evaluates finite-support endpoints from arm-specific
probabilities. `fit_discrete_cross_fitted_nuisances()` fits nuisance models;
`estimate_discrete_dual_bounds()` evaluates bounds and one-step inference.
Supply the ordered numeric outcome support explicitly. These functions use a
binary instrument and must not be mistaken for the multiple-instrument-category
extension. Smoothing parameters change the endpoints being evaluated.

Focused validation is included in the package tests run by `R CMD check`.

## Finite discrete instruments

For binary outcomes, the data-based cross-fitted estimator is:

```r
fit <- estimate_finite_iv_binary(
  data, outcome = 'Y', treatment = 'T', instrument = 'Z',
  covariates = c('X1', 'X2'), Lambda = 2, folds = 5)
fit$endpoints
fit$outward_sets
fit$diagnostics
```

It fits primitive nuisance functions and pairwise direct RKHS/Riesz ratios,
reconstructs binary probabilities, and evaluates selected-branch one-step
signals on held-out folds. It requires no threshold class. Sensitivity values
can be scalar or named by instrument category. The installed
`examples/finite_iv_binary.R` shows a complete dataset example. See
`?estimate_finite_iv_binary` for controls and inference conditions. Projection
and tie diagnostics must not be interpreted as formal compatibility tests.
The default primitive learner is spline logistic regression; Super Learner
can be selected explicitly. Inference assumes independent observations.

The package provides `finite_iv_bounds(pi,p0,p1,support,Lambda)` for arbitrary
instrument-category counts and externally supplied conditional probabilities,
and `fit_finite_iv_bounds()` for a simple empirical discrete-strata fit.
Both provide finite-support outcome plug-in bounds and explicit feasibility
diagnostics, not confidence intervals. See `?finite_iv_bounds` and the installed
its help page for shapes and category ordering. These are lower-level utilities,
not the fitted-estimator examples.

For continuous outcomes, `estimate_finite_iv_continuous()` provides an end-to-end
cross-fitted estimator with arbitrary IV labels and a chosen reference category:

```r
fit <- estimate_finite_iv_continuous(
  data, outcome = 'Y', treatment = 'T', instrument = 'Z',
  covariates = c('X1', 'X2'), Lambda = 2)
fit$endpoints
fit$outward_sets
fit$diagnostics
```

It fits direct RKHS/Riesz density ratios and fixed-class thresholds. The default
threshold class is linear with an intercept; a user-specified fixed spline basis
is illustrated in the installed `examples/finite_iv_continuous.R`. Knots and
boundaries must be prespecified rather than recomputed on evaluation samples.
Direct plug-in estimates and restricted-class one-step estimates are reported
separately, along with projection, normalization, and optimization diagnostics.
The confidence sets are pointwise, assume independent observations, and require
the documented nuisance-rate, margin, and influence-function stabilization
conditions. The function does not certify those conditions, perform formal
compatibility testing, or provide simultaneous bands. See
`?estimate_finite_iv_continuous` for controls and output fields.
`finite_iv_envelopes()` remains available as a lower-level envelope calculation.

## Binary outcome confidence bands

```r
library(ivmsm)
# lower_signal and upper_signal: n independent units by length(grid) matrices
# of cross-fitted one-step ATE signals, with identical row order.
band <- binary_grid_band(grid, lower_signal, upper_signal)
continuous_band <- binary_continuum_band(
  grid, band$table$simultaneous_lower, band$table$simultaneous_upper,
  at = seq(min(grid), max(grid), length.out = 101))
```

The first function exposes the Gaussian calibration used in the historical
OHIE workflow, without application paths or a LogSpline dependency. It does
not fit nuisances. It assumes independent analysis units and is not a
clustered or sampling-weighted calibration. Statistical validity requires the
paper's regularity conditions; this function cannot establish them.

The second function extends valid grid limits along the common sensitivity
path using monotonicity and the binary ATE endpoint Lipschitz bound. It cannot
turn invalid grid limits into valid inference. Neither function projects
incompatibilities or repairs crossed limits. Retain diagnostics from the
upstream estimator.
