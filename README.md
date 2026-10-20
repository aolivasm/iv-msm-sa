# IV–MSM sensitivity analysis

R package accompanying [Sharp Instrumental Variable Bounds under the
Marginal Sensitivity Model](https://arxiv.org/abs/2610.22513).

The paper combines the marginal sensitivity model with instrumental-variable
restrictions to obtain sharp bounds on potential-outcome means and average
treatment effects. The package implements the proposed estimation and inference
methods. Examples and code for reproducing the main numerical experiments and
the Oregon Health Insurance Experiment application are also included.

## Getting started

You do not need Git or to clone the repository. With R 4.4 or later:

1. On this GitHub page, click **Code → Download ZIP**.
2. Unzip the downloaded file.
3. Open RStudio and choose **Session → Set Working Directory → Choose Directory**.
   Select the extracted folder containing this `README.md` and the `packages` folder.
4. Run the following in the **R console** (not the terminal):

```r
install.packages("packages/ivmsm", repos = NULL, type = "source")
library(ivmsm)
```

Installation is needed once; use `library(ivmsm)` in each new R session.

### Continuous outcomes

Start with a continuous outcome, binary instrument, and a prespecified spline
threshold class, as in the paper's restricted-class approach:

```r
source(system.file("examples/continuous_fit.R", package = "ivmsm"))
```

The example explicitly defines the spline class and sets `Lambda0 = 2` and
`Lambda1 = 2` for the two instrument strata. It prints
endpoint estimates and outward confidence sets. Numerical diagnostics are optional.

### Binary outcomes and confidence bands

Binary outcomes are also studied in the main numerical experiments and the
Oregon Health Insurance Experiment application. For endpoint estimates and a
pointwise outward confidence set with a binary instrument, run:

```r
source(system.file("examples/binary_fit.R", package = "ivmsm"))
```

This example sets `Lambda0 = Lambda1 = 2`. For a simultaneous outward ATE
confidence band over the common sensitivity path `Lambda0 = Lambda1` in
`[1.2, 2]`, run:

```r
source(system.file("examples/binary_bands.R", package = "ivmsm"))
```

These examples fit models to small teaching datasets; they do not reproduce the
paper's simulations or application. The [package guide](packages/ivmsm/README.md)
provides controls and additional examples for smoothing, finite instruments, and
ordered outcomes. See the [installation guide](docs/installation.md) for optional
fitting dependencies.

## Numerical experiments and application

The [main-text reproduction guide](docs/reproduction.md) describes the continuous-
and binary-outcome numerical experiments and the Oregon Health Insurance Experiment
application. Experiment settings and computational instructions are kept
outside this README. Application data must be obtained separately following
the [provider access instructions](data/README.md).

Results and fitted objects are not distributed. The reproduction guide explains
how to run the analyses, collect their outputs, and generate the figures locally.

Code for ordered outcomes, finite discrete instruments, and binary-outcome
confidence bands is described in the [extension guide](docs/extensions_api.md).

## Citation

Please cite the [paper, arXiv:2610.22513](https://arxiv.org/abs/2610.22513),
and identify the code commit or release used in your analysis.

## License

Original code is distributed under the [MIT license](LICENSE). Third-party code
and datasets retain their own terms; this license does not grant permission to
redistribute the Oregon participant data.
