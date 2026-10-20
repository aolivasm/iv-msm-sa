\name{finite_iv_bounds}
\alias{finite_iv_bounds}
\alias{fit_finite_iv_bounds}
\alias{finite_iv_envelopes}
\title{General Finite Instrument Interfaces}
\description{Compute finite-support outcome bounds using any finite number of
instrument categories, or common-density envelopes for continuous outcomes.}
\usage{
finite_iv_bounds(pi, p0, p1, support, Lambda = 2, weights = NULL)
fit_finite_iv_bounds(data, outcome, treatment, instrument, support,
                     Lambda = 2, strata = character())
finite_iv_envelopes(r, gamma, Lambda = 2, reference = 1L)
}
\arguments{
\item{pi}{Numeric n by J matrix of treatment-one probabilities given covariates and IV category.}
\item{p0, p1}{Numeric n by J by K arrays of observed outcome probabilities conditional on treatment zero or one, covariates, and IV. Last dimension follows support.}
\item{support}{Strictly increasing numeric outcome support, of length K.}
\item{Lambda}{Finite scalar or length-J vector of sensitivity parameters, at least one. Category order must match the arrays; the same parameters are applied to both arms.}
\item{weights}{Nonnegative averaging weights across covariate rows. Not survey inference weights.}
\item{data}{Complete data frame.}
\item{outcome, treatment, instrument}{Column names. Treatment is coded zero/one.}
\item{strata}{Column names defining discrete covariate strata. Empty means no covariates. Continuous covariates should instead use externally fitted probabilities with finite_iv_bounds.}
\item{r}{Numeric n by J matrix of probabilities of the selected treatment arm.}
\item{gamma}{Numeric n by J density ratios of the observed arm outcome distribution in each category relative to the reference category, evaluated at common outcome/covariate pairs. Reference column equals one.}
\item{reference}{Reference IV column index.}
}
\details{The finite-support interface intersects category-specific probability
boxes and optimizes the common outcome mass over its simplex. It checks both
pointwise and normalization feasibility. Incompatible rows are retained, with
aggregate estimates missing; there is no projection or deletion of observations.
The data interface uses empirical cell probabilities and requires both treatments
in every stratum/instrument cell. It provides no cross-fitting or confidence
intervals. General machine-learning predictions can be supplied through the
probability interface. These bounds invoke the marginal IV model, not additional
joint potential-outcome restrictions.

For continuous outcomes, finite_iv_envelopes returns density envelopes relative
to the reference observed arm distribution, not latent-selection-ratio envelopes.
Pointwise noncrossing alone does not establish normalization feasibility. This
function does not fit thresholds or deliver a complete finite-IV estimator.}
\value{A list of estimates, conditional bounds, compatibility diagnostics, and
settings; or a list of density envelopes and crossing diagnostics. No formal
falsification test is performed.}
\seealso{For fitting the cross-fitted estimator to individual records, use
\code{\link{estimate_finite_iv_binary}} or
\code{\link{estimate_finite_iv_continuous}}. These utilities are not replacements
for the fitted estimators.}
