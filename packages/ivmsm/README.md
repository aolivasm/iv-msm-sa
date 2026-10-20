# ivmsm

R package for instrumental-variable bounds under the marginal sensitivity
model. See [the paper](https://arxiv.org/abs/2610.22513) for assumptions and theory.

## Fit the estimator to your data

Follow the [download and installation instructions](../../README.md#getting-started).
Start with a continuous outcome and binary instrument. Supply a data frame
`my_data` with columns `Y` (continuous outcome), `T_obs` (0/1 treatment received),
`Z` (0/1 instrument), and covariate `X`. Prespecify a threshold class; for example,
if `X` lies in [-1, 1], select the spline option:

```r
library(ivmsm)
Lambda0 <- 2
Lambda1 <- 2
fit <- estimate_finite_iv_continuous(
  my_data, outcome = "Y", treatment = "T_obs", instrument = "Z",
  covariates = "X", Lambda = c('0' = Lambda0, '1' = Lambda1),
  reference = '0', folds = 5, threshold_class = "spline",
  spline_spec = list(X = list(knots = 0, boundary = c(-1, 1))))
print(fit)
summary(fit)
```

Use `threshold_class="linear"` for linear thresholds. For several covariates,
give one specification per covariate, or use `spline_reference` with independent
reference covariate data to set the spline basis automatically. See the function
help for options. This fits nuisance functions,
including the direct RKHS/Riesz density-ratio learner, and evaluates cross-fitted
bounds and inference. It is a fitted analysis, not a calculation from supplied
population probabilities. Retain and inspect the diagnostics.

The threshold is learned within this fixed class. The resulting restricted-class
bounds can be wider than the unrestricted sharp bounds. The outward confidence
sets require the paper's regularity conditions.

For a binary outcome and binary instrument, no threshold class is needed:

```r
fit <- estimate_binary_dual_bounds(
  my_data, x_cols = c("X1", "X2"),
  Lambda0 = 2, Lambda1 = 2, K = 5)
print(fit)
```

Here encode `Y`, `T_obs`, and `Z` as 0/1 and replace the covariate names as needed.

For a **binary outcome and a finite-category instrument**, use:

```r
fit <- estimate_finite_iv_binary(
  my_data, outcome = "Y", treatment = "T_obs", instrument = "Z",
  covariates = c("X1", "X2"), Lambda = 2, folds = 5)
print(fit)
```

This fits the primitive nuisance models and direct RKHS/Riesz ratios with
cross-fitting. No threshold class is needed for binary outcomes. The default
primitive learner is spline logistic regression; use
`control = list(primitive_learner = "super_learner")` to select the Super Learner
family used in the experiments. These are independent-observation, pointwise
confidence sets, not simultaneous bands or a formal compatibility test.

## Run a complete example

`Lambda0` and `Lambda1` specify sensitivity in IV strata `Z=0` and `Z=1`,
not treatment arms. For more IV categories use named values, for example
`Lambda = c('0'=2, '1'=2, '2'=3)`. Changing these changes the assumptions.

| Example | Outcome | Instrument | Smoothing |
|---|---|---|---|
| `continuous_fit.R` | Continuous | Binary | Unsmoothed |
| `continuous_smoothing.R` | Continuous | Binary | `envelope_smoothing_epsilon` |
| `binary_fit.R` | Binary | Binary | `epsilon=0` or positive |
| `ordered_outcome.R` | Ordered 0/1/2 | Binary | Envelope and threshold controls |
| `finite_iv_binary.R` | Binary | Three categories | Unsmoothed only |
| `finite_iv_continuous.R` | Continuous | Three categories | Unsmoothed only |
| `binary_bands.R` | Binary | Binary | Unsmoothed; simultaneous band |

Create and inspect the primary teaching data with `d <- mock_iv_data('continuous', 2)`
and `head(d)`. Columns are outcome `Y`, treatment `T_obs`, instrument `Z`,
and covariate `X`. The ordered-outcome fitted interface currently supports a
binary IV, not the ordered-outcome/multi-category-IV combination.

The examples create small artificial teaching datasets and call the actual
estimators. They are not paper simulation results. Install the example
dependencies once in the R console:

```r
install.packages(c("nnet", "nloptr"))
```

Choose **one** example below and run its line in the R console. Unlike installing
or loading the package, these commands run model fits. Small ranks, fixed
penalties and fewer folds keep the demonstrations manageable; use the documented
production settings for an analysis.

```r
source(system.file('examples/continuous_fit.R',package='ivmsm'))
source(system.file('examples/continuous_smoothing.R',package='ivmsm'))
source(system.file('examples/finite_iv_continuous.R',package='ivmsm'))
source(system.file('examples/finite_iv_binary.R',package='ivmsm'))
source(system.file('examples/ordered_outcome.R',package='ivmsm'))
source(system.file('examples/binary_fit.R',package='ivmsm'))
source(system.file('examples/binary_bands.R',package='ivmsm'))
```

Use `?estimate_binary_dual_bounds`, `?estimate_finite_iv_binary` and
`?estimate_finite_iv_continuous` for input contracts. In the data templates above,
`my_data` means your own data frame; it is not a dataset installed by the package.

The continuous smoothing example disables simulation-specific oracle calculations
with `compute_oracle_layers=FALSE`. Historical truth/benchmark fields in that
research interface must not be interpreted as causal truths for user data.
For multiple-category IVs with continuous outcomes,
`estimate_finite_iv_continuous()` fits the nuisances and fixed-class thresholds,
and returns cross-fitted one-step estimates, pointwise outward confidence sets,
and diagnostics. It accepts arbitrary column names and IV labels; see its help
page for assumptions and controls. `finite_iv_envelopes()` also exposes the
underlying envelope calculation. Frozen experiment pipelines remain unchanged.

Additional learner dependencies are listed in DESCRIPTION. The confidence-band
example fits the binary estimator and uses its one-step signals across a common
sensitivity grid; it does not fabricate signals. Valid inference requires the
paper's regularity conditions and independent analysis units.

Low-level probability and envelope utilities remain available for development
and checking identities. They are not substitutes for the fitted estimators above.

## Results and sensitivity intervals

`print(fit)` shows the analysis settings, endpoint estimates and outward confidence
sets. `summary(fit)` provides a compact structured result with the same clean
default display. Numerical checks are optional: use `print(fit, diagnostics=TRUE)`.
The original fields, including `fit$diagnostics`, remain available. “Adjusted fitted values” means
numerical safeguards were needed to enforce constraints; it is not a formal
compatibility test or evidence that confidence coverage is valid.

For binary outcomes and a binary IV, the band interface accepts an interval:

```r
band <- estimate_binary_sensitivity_band(
  my_data, x_cols=c('X1','X2'), Lambda_range=c(1.2,2))
print(band)
plot(band)
```

This follows `Lambda0 = Lambda1`, fits nuisances once, and constructs a simultaneous
outward ATE band using a calibration grid and a continuum extension. It is not a
two-dimensional region for independently varying sensitivities. The wrapper uses
unsmoothed envelopes; smoothing is not silently applied. Positive smoothing in
the pointwise examples changes the bounds being evaluated.

The band example calibrates jointly at 21 sensitivity values and plots continuum
limits evaluated at 201 locations. `fit$table` retains grid results and
`fit$continuum` retains the continuum limits. Crossed limits are flagged;
the internal `empty` indicator is not printed as a routine results column.
The binary-band example simulates 2,000 teaching records with randomized
encouragement, nonadherence, and an unmeasured common cause of treatment and
outcome. Its range [1.2, 2] illustrates widening limits and inclusion of zero
near Lambda=1.56 in this selected example. These are not OHIE data, the paper's
simulation design, or a coverage assessment.
Retrieve the grid values with:

```r
grid_results <- fit$table[, c("Lambda", "lower", "upper",
                             "simultaneous_lower", "simultaneous_upper")]
grid_results
grid_results[which.min(abs(grid_results$Lambda - 1.6)), ]
```

`lower` and `upper` are endpoint estimates; `simultaneous_lower` and
`simultaneous_upper` are the jointly calibrated confidence limits. The nearest
grid-point command returns the actual `Lambda` used, not an interpolated value.
