\name{ivmsm_estimators}
\alias{estimate_sharp_bounds_plugin}
\alias{estimate_sharp_bounds_one_step}
\alias{estimate_binary_dual_bounds}
\alias{fit_binary_cross_fitted_direct_rkhs_nuisances}
\alias{evaluate_binary_bounds_from_nuisances}
\alias{fit_discrete_cross_fitted_nuisances}
\alias{estimate_discrete_dual_bounds}
\alias{discrete_arm_endpoint}
\title{Cross-Fitted Binary-Instrument Estimators}
\description{Packaged versions of the existing continuous, binary, and ordered
outcome estimators. These entry points retain their research interfaces and
defaults. Inspect function formals before launching an expensive fit.}
\details{Inputs use columns Y, T_obs, Z, and named covariate columns. Z and
T_obs are coded zero/one. The continuous entry points are
estimate_sharp_bounds_plugin and estimate_sharp_bounds_one_step. Binary
estimate_binary_dual_bounds uses direct RKHS/Riesz nuisance estimation.
For ordered outcomes use fit_discrete_cross_fitted_nuisances followed by
estimate_discrete_dual_bounds, with numeric support explicitly supplied.
discrete_arm_endpoint also accepts specified arm probabilities directly.

The ordered estimator's envelope and threshold smoothing defaults change its
target; consult its returned exact-plugin and smooth-one-step rows separately.
Projection and feasibility diagnostics are not formal compatibility tests.
Oracle modes are for the built-in simulation designs, not arbitrary user data.
Optional fitting libraries must be installed separately. These functions do not
automatically reproduce the frozen paper experiment configurations.}
\value{Lists containing estimates, fitted nuisances or endpoint components,
and diagnostics, depending on the entry point.}
