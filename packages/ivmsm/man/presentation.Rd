\name{ivmsm_fit}
\alias{ivmsm_fit}
\alias{print.ivmsm_fit}
\alias{summary.ivmsm_fit}
\alias{print.summary.ivmsm_fit}
\alias{print.ivmsm_band}
\alias{summary.ivmsm_band}
\alias{plot.ivmsm_band}
\alias{mock_iv_data}
\alias{estimate_binary_sensitivity_band}
\title{Readable fitted results, teaching data, and binary sensitivity bands}
\usage{
mock_iv_data(outcome = c("binary", "ordered", "continuous"), iv_categories = 2L)
estimate_binary_sensitivity_band(data, x_cols = "X", Lambda_range = c(1.2, 2),
  grid_points = 21L, confidence_level = .95, draws = 10000L, seed = 1L, ...)
\method{print}{ivmsm_fit}(x, ...)
\method{summary}{ivmsm_fit}(object, ...)
\method{print}{summary.ivmsm_fit}(x, digits = 4, diagnostics = FALSE, ...)
\method{print}{ivmsm_band}(x, digits = 4, ...)
\method{summary}{ivmsm_band}(object, ...)
\method{plot}{ivmsm_band}(x, xlab = expression(Lambda[0] == Lambda[1]),
  ylab = "ATE", main = "Simultaneous outward confidence band", ...)
}
\arguments{
\item{outcome}{Type of outcome in the deterministic teaching dataset.}
\item{iv_categories}{Number of IV categories, at least two.}
\item{data}{Individual-level data with Y and T_obs coded 0/1, binary Z coded 0/1, and numeric covariates.}
\item{x_cols}{Covariate column names.}
\item{Lambda_range}{Increasing pair of endpoints, at least one, for the common sensitivity path Lambda0 = Lambda1.}
\item{grid_points}{Number of equally spaced calibration grid points.}
\item{confidence_level}{Simultaneous confidence level.}
\item{draws}{Gaussian multiplier draws.}
\item{seed}{Random seed.}
\item{x, object}{A fitted result or summary.}
\item{digits}{Printed precision.}
\item{diagnostics}{Whether to display optional numerical checks. Defaults to FALSE. Also accepted by print.ivmsm_fit through ... .}
\item{xlab, ylab, main}{Plot labels. Additional plot arguments are passed through ... .}
\item{...}{For fitting, controls passed to estimate_binary_dual_bounds, such as K and primitive_learner. Do not override the data, sensitivity, smoothing, confidence, or seed arguments. Otherwise reserved display arguments.}
}
\description{Prints one-step endpoint estimates and outward sets. Numerical checks are opt-in via print(fit, diagnostics=TRUE); full diagnostic objects are always retained. Creates clearly artificial teaching records, or fits a binary estimator once and calibrates a sensitivity band.}
\details{Sensitivity parameters index instrument categories, not treatment arms. Printed numerical adjustments refer to safeguards imposed on fitted nuisance values, not formal compatibility testing. When many diagnostic groups are present, the summary displays the maximum percentage, not an average over unequal folds. No violations are removed from the full diagnostics.

Band printing reports the interval and calibration-grid size. The plot displays
continuum limits at the evaluation locations. Full grid results, including the
empty-set indicator, remain in table; continuum results remain in continuum.
Crossed confidence limits are flagged and are never repaired or shaded as valid intervals.

The band wrapper is only for binary outcomes and a binary IV along a common sensitivity path. It uses unsmoothed envelopes (epsilon=0), reuses fitted nuisances at every grid value, and applies Gaussian multiplier calibration. Continuum limits use monotonicity and the binary ATE endpoint Lipschitz bound, not interpolation. They require valid uniform inference assumptions and are evaluated at 201 locations across the requested interval. Arbitrary two-dimensional sensitivity regions, finite-category-IV bands and smoothed-envelope bands are not implemented by this wrapper.

Mock data are deterministic records for learning the interface, not an empirical dataset, a Monte Carlo experiment, or evidence for inferential validity. Outcome is binary 0/1, ordered 0/1/2, or continuous; treatment is binary, IV categories run from zero to iv_categories minus one, and X is a numeric covariate.}
\value{mock_iv_data returns a data frame. Fitting returns an ivmsm_band list containing the grid table, critical value, continuum limits and diagnostics. Summary methods return compact display objects; print methods return their input invisibly.}
\examples{
head(mock_iv_data("binary", 3))
\dontrun{
source(system.file("examples/binary_fit.R", package="ivmsm"))
source(system.file("examples/binary_bands.R", package="ivmsm"))
}
}
