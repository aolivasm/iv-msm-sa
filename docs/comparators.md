# Comparator software

The main-text workflows call the authors' original implementations, not local
substitutes. Install these before running a comparator analysis. Their licenses
and access terms remain separate from this repository's MIT license.

## Tan

Install the dependencies in R:

```r
install.packages(c("RCAL", "lpSolveAPI"))
```

Obtain **RCALsa version 1.0** from the author's
[software distribution](https://sites.sas.rutgers.edu/ztan/Publication/RCALsa_1.0.zip).
Extract the archive and install its R source package with `R CMD INSTALL`
from a terminal. If installed in a separate library, add that library to
`R_LIBS_USER` before starting R. `RCALsa` is not replaced by `RCAL`.

## Dorn–Guo–Kallus

Install the R dependencies:

```r
install.packages(c("grf", "dplyr", "progress", "ggplot2"))
```

From a terminal, obtain the audited author implementation:

```sh
git clone https://github.com/CausalML/DoublyValidDoublySharpSensitivity.git external/DVDS
git -C external/DVDS checkout 26cbbf6e90dcd0648c837c0d58568375291dd522
export IVMSM_DVDS_OFFICIAL_SOURCE="$(pwd)/external/DVDS"
```

Run these commands from the repository root. The absolute environment-variable
path allows experiment scripts to change working directories safely. The adapter
checks the commit before fitting. Keep the Git metadata in this external checkout.
No external source is redistributed inside the release.

## Check without fitting

```sh
Rscript scripts/check_comparator_dependencies.R
```

Historical runs used R 4.4.3 and RCALsa 1.0. The runner records installed package
versions in each run signature. Different dependency versions may change random
fits; the supplied aggregate results preserve the reported numerical values.
