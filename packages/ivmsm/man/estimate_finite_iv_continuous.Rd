\name{estimate_finite_iv_continuous}
\alias{estimate_finite_iv_continuous}
\title{Cross-Fitted Continuous Outcome Bounds with a Finite Instrument}
\description{Fits primitive probabilities, direct RKHS density ratios, fixed-class
thresholds, and conditional means, then evaluates plug-in and one-step bounds
on held-out observations. Supports arbitrary finite instrument labels.}
\usage{
estimate_finite_iv_continuous(data, outcome, treatment, instrument,
  covariates = character(), Lambda = 2, reference = NULL, folds = 5L,
  fold_id = NULL, threshold_basis = NULL, control = list(),
  confidence_level = .95, seed = 1L, keep_signals = TRUE,
  threshold_class = c("linear", "spline"), spline_spec = NULL,
  spline_reference = NULL)
}
\arguments{
\item{data}{Complete data frame of independent observations.}
\item{outcome}{Name of a finite numeric outcome column.}
\item{treatment}{Name of a treatment column coded zero/one.}
\item{instrument}{Name of an IV column with at least two observed categories.}
\item{covariates}{Numeric covariate column names. Encode categorical covariates explicitly. Empty means no covariates.}
\item{Lambda}{Finite scalar or category-specific values at least one, common across treatment arms. Named values must match IV labels exactly. Unnamed vectors follow factor level order before moving the reference to the first position.}
\item{reference}{Reference IV label. Defaults to the first factor level. Returned instrument_levels gives internal order.}
\item{folds}{Number of stratified outer folds.}
\item{fold_id}{Optional numeric vector assigning every observation to an outer fold. Every training IV/treatment cell must have sufficient observations.}
\item{threshold_basis}{Prespecified rowwise function of the covariate data frame returning a finite numeric design matrix, with an intercept of ones as its first column. NULL uses a fixed linear class. A spline basis must supply fixed knots and boundaries; data-adaptive basis construction is rejected.}
\item{control}{Named list of computational options; see Details. Unknown names are rejected.}
\item{confidence_level}{Pointwise confidence level.}
\item{seed}{Nonnegative integer seed. Caller RNG state is restored.}
\item{keep_signals}{Retain observation-level signals. These may contain sensitive information for real applications.}
\item{threshold_class}{Linear (default) or additive natural cubic spline thresholds. A custom threshold_basis remains available with the default class argument.}
\item{spline_spec}{Named list, one entry per covariate, containing knots and boundary (two values), or linear=TRUE for an indicator/linear term. Used with threshold_class="spline".}
\item{spline_reference}{Alternative to spline_spec: independent reference covariate data. Freezes quartile interior knots and range boundary knots once before cross-fitting; binary indicators remain linear. Do not pass the analysis sample as reference data. Covariates must be numeric; encode categorical covariates as indicators.}
}
\details{
All nuisance and threshold fits use outer-training observations only. Each
nonreference IV category is paired with the reference to estimate the observed
arm density ratio using direct Riesz risk, not classification. Gamma bandwidth
is fixed at the training median-distance variance M/2. By default penalty
C/sqrt(n_arm) is selected by minimum mean risk across two independent three-fold
partitions. n_arm counts observations of the chosen arm in that pair.

Primitive defaults are logistic splines; Gaussian spline regressions estimate
conditional envelope and score means. These are configurable learners, not a
claim of correct specification. The function does not automatically reproduce
historical experiment configurations. No oracle functions or population targets
are used.

Control options are: primitive_learner (glm_spline by default, also
super_learner, library_cv, gam, ranger, glmnet); spline_df (3);
super_learner_candidates (the simulation library); primitive_cv_folds (3);
mean_learner (spline_fixed, also spline_cv, gam, ranger, library_cv);
mean_control (NULL, passed to the existing conditional-mean fitter);
probability_floor (.05); instrument_floor (min(.005,.5/J)); mean_margin (1e-4);
gamma_C (10 raised to seq(-1.5,1,.5)); gamma_cv_folds (3); gamma_cv_repeats (2);
gamma_backend (nystrom or exact); gamma_rank (NULL, existing adaptive rank rule);
gamma_fixed_lambda (NULL; a positive value bypasses penalty CV, for debugging);
gamma_bandwidth_sample (1000); threshold_maxit (60000);
compatibility_action (project or stop); projection_warning (.05).

Projection mixes nonreference density ratios toward one just enough to prevent
pointwise envelope crossing. Conditional envelope means are stabilized around
one for threshold fitting. Both operations are reported, not treated as formal
compatibility tests. The stop option rejects estimated pointwise crossings.
No crossed final interval is repaired. Threshold rank deficiency or optimizer
failure stops the fit rather than silently changing the threshold class.

The main one-step estimand is the restricted-class pseudo-bound endpoint.
Direct extremizer evaluation is also returned, along with pseudo plug-in
evaluation and their normalization discrepancy. Outward confidence sets require
the paper's nuisance-rate, margin, and influence-function stabilization
conditions. Two-sided endpoint intervals additionally require negligible
optimization error and regularity at the optimizer. Neither validity nor
efficiency is certified by successful computation. Confidence sets are pointwise
in Lambda and separate for each target, not jointly calibrated across targets
or a sensitivity grid. Clustered and sampling-weighted inference are not
implemented by this interface.
}
\value{A list containing endpoints, outward_sets for mu0/mu1/ATE, direct_plugin,
pseudo_plugin, diagnostics, tuning (including gamma repeat-level CV), fold_id,
instrument_levels, Lambda, control, threshold class and dimension, warnings,
inference qualifications, and optionally observation-level signals.}
\examples{
set.seed(9)
d <- expand.grid(t=0:1,z=c('low','mid','high'),rep=1:30)
d$x <- runif(nrow(d),-1,1)
d$y <- .4*d$t+.1*d$x+rnorm(nrow(d),sd=.3)
# Fast demonstration, not the production tuning configuration.
fit <- estimate_finite_iv_continuous(d,'y','t','z','x',folds=2,
  control=list(gamma_fixed_lambda=.1,gamma_rank=6,spline_df=2))
fit$outward_sets
}
