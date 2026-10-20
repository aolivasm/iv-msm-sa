\name{estimate_finite_iv_binary}
\alias{estimate_finite_iv_binary}
\title{Cross-fitted binary-outcome bounds with a finite-category instrument}
\usage{
estimate_finite_iv_binary(data, outcome, treatment, instrument,
  covariates = character(), Lambda = 2, reference = NULL,
  folds = 5L, fold_id = NULL, control = list(),
  confidence_level = .95, seed = 1L, keep_signals = TRUE)
}
\arguments{
\item{data}{Data frame with one independent observation per row. Missing values are rejected.}
\item{outcome, treatment, instrument}{Column names. Outcome and treatment must be numeric 0/1. Instrument labels may be arbitrary, with at least two observed categories.}
\item{covariates}{Numeric covariate column names. Encode categorical covariates explicitly.}
\item{Lambda}{Common sensitivity value at least one, or one value per instrument category. Named values must exactly match instrument labels. Unnamed vectors follow sorted factor levels before reference reordering.}
\item{reference}{Instrument label for baseline outcome fitting and pairwise density ratios. Default: first sorted factor level.}
\item{folds}{Number of outer folds.}
\item{fold_id}{Optional integer fold assignments in 1 through folds. Otherwise folds are stratified by treatment and instrument.}
\item{control}{Named controls: primitive_learner (default glm_spline), spline_df (3), super_learner_candidates, primitive_cv_folds (3), probability_floor (.05), instrument_floor (min(.005,.5/J)), gamma_C (10 to the powers seq(-1.5,1,.5)), gamma_cv_folds (3), gamma_cv_repeats (2), gamma_backend (nystrom), gamma_rank (NULL), gamma_fixed_lambda (NULL), gamma_bandwidth_sample (1000), compatibility_action (project, none, or stop), and projection_warning (.05). Unknown controls are rejected.}
\item{confidence_level}{Pointwise confidence level, between zero and one.}
\item{seed}{Nonnegative integer seed. The caller's random-number state is restored.}
\item{keep_signals}{Retain observation-level one-step signals for downstream inference. Keep these private when using participant data.}
}
\description{Fits nuisance functions from training folds and evaluates selected-branch one-step estimates on held-out observations. Supports binary outcomes and any finite number of observed instrument categories, without oracle inputs or a threshold class.}
\details{The instrument distribution uses sequential binary regressions and treatment probabilities are fitted within instrument categories. Baseline outcome probabilities are fitted in the reference category separately by treatment. Direct RKHS/Riesz ratios for each reference/nonreference pair reconstruct the other outcome probabilities. Kernel variance is fixed at half the median squared training distance; only the penalty is tuned by repeated inner CV. Conditional normalization and clipping are recorded.

For each treatment, conditional lower means maximize the 2J binary probability restrictions; upper means minimize them. The selected-branch influence correction uses instrument and treatment inverse probabilities. Setting compatibility_action to project mixes reconstructed probabilities toward the reference until feasible, as in the experimental implementation; it never replaces crossed bounds with their midpoint. Raw crossing, projection, shrinkage, normalization and tie diagnostics are returned. Setting none preserves crossings; stop raises an error. These are not formal compatibility tests.

Inference requires appropriate nuisance rates, branch-margin/no-tie and influence-function stabilization conditions. No numerical check certifies these assumptions. Stabilization alone does not make confidence sets valid under misspecification. Confidence sets are pointwise in Lambda and separate for each parameter, not simultaneous across parameters. Clustered and sampling-weighted inference are not implemented. The general interface defaults to spline-logistic primitive fits; select super_learner explicitly to use the experiment's primitive learner family.}
\value{A list with endpoints (mu0, mu1 and ATE endpoints, standard errors and two-sided intervals), outward_sets, plugin (uncorrected averaged conditional endpoints), diagnostics, tuning, fold_id, instrument_levels, Lambda, control, n, seed, signals, inference and warnings. Empty outward sets and reversed estimates are not silently repaired.}
\examples{
\dontrun{
source(system.file("examples/finite_iv_binary.R", package="ivmsm"))
}
}
\seealso{\code{\link{estimate_finite_iv_continuous}}, \code{\link{finite_iv_bounds}}}
