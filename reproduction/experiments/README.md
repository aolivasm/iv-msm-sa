# Main-text experiment code

Run these experiments through `../scripts/paper.R`, from `reproduction/`.
Choose settings in `../configurations/main_experiments.R`. Do not source a
development-era study or change the package defaults to reproduce the paper.

The continuous-outcome engine has six files:

| File | Purpose |
|---|---|
| `design.R` | Generate the paper's covariates, instrument, treatment, and outcome |
| `nuisance_models.R` | Estimate the treatment/instrument probabilities and conditional regressions |
| `rkhs_gamma.R` | Fit and tune the RKHS/Riesz density ratio γ |
| `threshold_splines.R` | Fit the paper's fixed spline threshold class |
| `estimation_inference.R` | Evaluate bounds, influence-function corrections, and confidence sets |
| `population_targets.R` | Calculate the numerical population targets used to assess coverage |

The binary engine does not need threshold fitting. Comparator engines use the
official Tan and Dorn implementations, with the paper's X-adjusted specification.
The compatibility engine keeps outcome-specific implementations because their
published nuisance-fitting recipes differ; all sample sizes use the same driver.

There are no growing-spline γ or RBF-threshold experiments in this workflow.
Shared numerical helpers remain where the main estimator calls them. Internal
function names are retained where practical to make comparison with the original
implementation straightforward.
