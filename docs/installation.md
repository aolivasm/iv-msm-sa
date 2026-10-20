# Installation and dependencies

Use R 4.4 or later. The **repository root** is the downloaded folder containing
`README.md`, `packages/`, and `reproduction/`. In R or RStudio, set this folder as
your working directory, then run:

```r
install.packages("packages/ivmsm", repos = NULL, type = "source")
library(ivmsm)
```

The first command installs the local package into your R library. It is needed
once, or again after updating the package. The second loads the installed
package in the current R session; it does not run simulations or the application.

Alternatively, from a **terminal** whose current directory is the repository
root, install with `R CMD INSTALL packages/ivmsm`. This is a terminal command,
not an R-console command. Then open R and run `library(ivmsm)`.

The package contains no participant data. Its examples are available through
`system.file('examples', package='ivmsm')`.

For a source archive, run `R CMD build packages/ivmsm`, then install the generated
archive with `install.packages('ivmsm_0.1.1.tar.gz', repos=NULL, type='source')`.
Only base/recommended libraries are required to load the package. Install optional
learners explicitly before using them; ordered-outcome fitting also uses `nnet`.

## Data-free checks

Run these commands inside `reproduction/` (from the repository root,
first run `cd reproduction`):

```sh
Rscript --vanilla scripts/public/check_environment.R minimal
Rscript --vanilla scripts/public/run_smoke_tests.R
```

These checks use base R. They validate selected binary population identities
and deterministic continuum-band properties, not fitted nuisance models or
Monte Carlo coverage.

## Fitting and plotting dependencies

The main application and plotting workflows use the following packages.
Install them explicitly in an appropriate R library if needed:

```r
install.packages(c('mgcv', 'glmnet', 'ranger', 'xgboost', 'SuperLearner',
                  'nnls', 'nloptr', 'hal9001', 'haven', 'digest',
                  'ggplot2', 'patchwork'))
```

```sh
Rscript --vanilla scripts/public/check_environment.R fitting
```

Some packages require system compilers or platform-specific dependencies. The
check reports availability and versions; it does not install or update anything.
The authors' comparator software has separate installation instructions in
[comparators](comparators.md).

Installing current packages is not exact paper reproduction. Preserve the
package versions, seeds, configuration, and signatures recorded by each frozen
experiment. A validated repository-wide dependency lockfile is not yet supplied.
Do not generate one from an unrelated personal library and call it the paper's
environment.
