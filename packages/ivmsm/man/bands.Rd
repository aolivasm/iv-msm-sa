\name{binary_grid_band}
\alias{binary_grid_band}
\alias{binary_continuum_band}
\alias{binary_flat_arm_bounds}
\title{Binary Population Bounds and Outward Band Calibration}
\description{Population binary bounds, simultaneous Gaussian grid calibration,
and a deterministic continuum extension along a common sensitivity path.}
\usage{
binary_grid_band(grid, lower_signal, upper_signal, level = .95,
                 draws = 10000L, seed = 1L)
binary_continuum_band(grid, lower, upper, at = grid, lipschitz = 2)
binary_flat_arm_bounds(r0, r1, a0, a1, Lambda)
}
\arguments{
\item{grid}{Strictly increasing sensitivity grid, at least one.}
\item{lower_signal, upper_signal}{n by grid-length matrices of cross-fitted ATE one-step signals. Rows are independent units aligned across columns.}
\item{level}{Confidence level.}
\item{draws}{Number of Gaussian draws, at least 100.}
\item{seed}{Random seed; caller RNG state is restored.}
\item{lower, upper}{Simultaneous grid limits.}
\item{at}{Evaluation points inside the grid range.}
\item{lipschitz}{Valid binary ATE endpoint Lipschitz constant, at least two.}
\item{r0, r1}{Selected treatment-arm probabilities in the two IV strata.}
\item{a0, a1}{Observed binary success probabilities in that arm and IV stratum.}
\item{Lambda}{Common finite sensitivity parameter, at least one.}
}
\details{Grid calibration assumes the paper's uniform expansion and stability
conditions. It does not support cluster or sampling-weighted inference. The
continuum extension uses monotonicity and a valid Lipschitz bound; it cannot
establish coverage of invalid input limits. No crossed intervals are repaired.}
\value{Grid calibration returns a list containing a table, critical value, and
covariance. Continuum extension returns a data frame. Population bounds return
a list with lower and upper endpoints.}
\examples{
\dontrun{
# Fit nuisances and construct bands from the estimator's one-step signals.
source(system.file("examples/binary_bands.R", package="ivmsm"))
}
}
