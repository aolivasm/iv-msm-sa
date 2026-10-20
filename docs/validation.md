# Release validation

## Code-only release, 2026-10-19

Saved results and all RDS files are excluded. Configuration metadata is generated
locally from readable settings. Collector tests use artificial deterministic
fixtures, not the paper's estimates. The updated collector test and compatibility
bookkeeping tests passed; preparation was checked for idempotence. All 403 R
sources parsed (404 after adding the small target-integration test). The small
deterministic target-integration fixture and public smoke tests also passed.
The full experiments and population-target integrations were not
rerun for this packaging change. The notes below describe earlier validation,
including checks performed when reference aggregates were available locally.

Validated locally on 2026-10-18–19 using R 4.4.3 on macOS arm64.

- Package 0.1.1 built and checked from an isolated copy: `R CMD check --no-manual`
  returned Status: OK. All seven package test scripts passed.
- Public population/band smoke tests passed without participant data.
- All 400 distributed R source files parsed successfully.
- All five experiment preflights passed from an isolated copy, including
  source manifests, population-target/configuration checks, and dependencies.
- The n=12000 compatibility mock suite passed: generator parity, task contract,
  calibration subgrid, atomic persistence/resume, failure retention, paired
  fingerprints, and corrupt/foreign/duplicate checkpoint rejection.
- RKHS-only OHIE band evaluation reproduced all five saved application endpoint
  estimates and standard errors within 1e-10. This check used existing authorized
  local inputs and 399 calibration draws, with no nuisance refit. The public
  band script defaults to 200,000 draws on the full 81-point grid.
- The n=6,000 compatibility mock suite passed the same bookkeeping checks.
- The main-results collector passed complete-input and failed-replication tests.
- Figures 3–6 regenerated from saved aggregates and were visually inspected.
- ED-only preparation preserved the original 22,170 analysis rows and covariates
  exactly, using a read-only comparison with the authorized local input.
- Official comparator dependencies were available and the DVDS source matched
  the pinned commit. This was a dependency check, not a new comparator fit.
- Continuous comparator fitting functions matched the original source exactly;
  binary comparator changes were limited to dependency discovery and source
  verification. The analytic introductory illustration also rendered successfully.
- All 37 distributed RDS objects passed recursive checks for participant-ID
  fields, executable functions and environments; nested source manifests matched.
- The binary finite-IV data interface installed successfully. Its selected-branch
  bounds and scores matched the original experimental implementation on deterministic
  inputs. Tests passed for derivatives, category-specific sensitivity values,
  projection, input validation, cross-fitting orchestration and standard errors.
  Nuisance learners were mocked; repeated-CV and fixed-bandwidth dispatch were
  checked without fitting. The installed dataset example was not executed.
  An updated `R CMD check --no-manual --no-tests --no-examples` returned Status: OK;
  the new deterministic/mock test script was run separately and passed.
- User examples now call fitted estimators; the confidence-band example reuses
  estimated nuisances and their one-step signals, rather than artificial signals.
  Installed examples parsed and the band-example data flow passed with mocked
  nuisance fitting and calibration. Population-identity checks are kept under
  tests, not presented as estimator examples. Example model fits were not rerun.
- Readable print/summary methods, sensitivity labels, the ATE confidence-set
  row name, and teaching-data constructors passed deterministic tests. The
  sensitivity-interval wrapper and installed examples passed mocked data-flow
  checks without fitting learners or running multiplier calibration. Smoothing
  examples explicitly distinguish supported binary-IV methods from the
  unsmoothed finite-IV interfaces.
- The named spline option passed deterministic checks for additive multidimensional
  bases, fixed rowwise evaluation across folds, independent-reference basis
  construction, and invalid inputs. No nuisance fits or simulations were run.
- Distributed RDS objects include simulation endpoint estimates, aggregate
  application results, population targets, configuration and provenance. No
  participant records or fitted application caches are included.

These checks do not establish statistical assumptions or rerun the full Monte
Carlo experiment. No new 500-replication study or full OHIE nuisance fit was run
for this release. The GitHub Actions workflow is provided but has not yet run
on GitHub. A network-restricted local check could not query current CRAN indexes;
it used the installed dependencies below.

## Tested dependency versions

The updated binary-band teaching example was subsequently run locally to check
its display and grid-value extraction. The revised example simulates 2,000
records with randomized encouragement, nonadherence, and latent confounding
(seed 72). Over [1.2, 2], the upper 95% grid limit rose from 0.364 to 0.484 and
the lower limit fell from 0.102 to -0.098; the continuum limits first included
zero at the evaluated location 1.564. This selected teaching example fit is
not a Monte Carlo coverage validation or an application result.

| Package | Version |
|---|---|
| nnet | 7.3.20 |
| nloptr | 2.2.1 |
| mgcv | 1.9.1 |
| glmnet | 4.1.8 |
| ranger | 0.17.0 |
| xgboost | 3.2.1.1 |
| SuperLearner | 2.0.29 |
| nnls | 1.6 |
| hal9001 | 0.4.6 |
| haven | 2.5.4 |
| digest | 0.6.37 |

This records the validation environment, not a claim that all historical runs
used these versions. Frozen runs record their own runtime signatures.

## Before uploading

From the repository root, run `python3 reproduction/scripts/public/audit_release.py` on the clean, unmodified release
folder to verify `SHA256SUMS`, unexpected files, and selected credential patterns.
The release deliberately has no Git history. Publish this curated source tree,
not the development repository's historical branches, tags, or output files.
Do not merely copy it onto the old history and change that repository's visibility.
Choose a new repository or explicitly review a history-replacement plan first.

The distribution supplies main-study estimation, comparator, compatibility,
collection and plotting workflows. The newly connected workflows were checked
without rerunning fitted experiments; this is not a claim of an independently
reproduced full Monte Carlo study. Authorized OHIE inputs and separately licensed
comparator software must be obtained by the user. See [reproduction](reproduction.md).
## Shared experiment configuration (October 20, 2026)

The nine paper configurations now use five shared source engines. All nine
configuration objects were compared with the previous publication drivers and
matched exactly. Statistical dependencies and compatibility adapters matched
byte-for-byte in 279 comparisons. Sample sizes, sensitivity grids, seeds,
learner/class choices, and outcome-specific tuning were retained.

Verification covered configuration dispatch and task contracts, repeated
preparation, changed-source rejection, compatibility checkpoint/resume and
failure handling, main-results collection fixtures, and existing smoke tests.
These are software checks, not new Monte Carlo studies. No simulations or
application fits were run for this consolidation.
