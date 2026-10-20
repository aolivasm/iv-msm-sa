# Bounded continuous-outcome DGP with uniform structural variables.
#
# X, U, and the potential-outcome shocks are independent Uniform(-1, 1).
# The instrument satisfies Z | X ~ Bernoulli(expit(z_x_coefficient * X)),
# treatment follows a logistic model in (X, Z, U), and
#
#   Y^t = treatment_effect * t + outcome_x * X + outcome_u * U
#         + outcome_error_scale * epsilon_t.
#
# The file also supplies exact oracle e(X), pi_z(X), and gamma_t(Y, X). The
# oracle density ratio uses analytic integration over the interval of U values
# compatible with a given (Y, X); no numerical quadrature is required.

uniform_dgp_softplus <- function(value) {
  log1p(exp(-abs(value))) + pmax(value, 0)
}

uniform_dgp_mean_expit <- function(eta, confounder_strength) {
  if (abs(confounder_strength) <= 1e-10) return(stats::plogis(eta))
  (
    uniform_dgp_softplus(eta + confounder_strength) -
      uniform_dgp_softplus(eta - confounder_strength)
  ) / (2 * confounder_strength)
}

uniform_dgp_treatment_integral <- function(
    eta,
    confounder_strength,
    lower,
    upper,
    treatment
) {
  width <- pmax(upper - lower, 0)
  if (abs(confounder_strength) <= 1e-10) {
    treated <- width * stats::plogis(eta)
  } else {
    treated <- (
      uniform_dgp_softplus(eta + confounder_strength * upper) -
        uniform_dgp_softplus(eta + confounder_strength * lower)
    ) / confounder_strength
  }
  if (treatment == 1L) treated else width - treated
}

uniform_compatibility_designs <- function() {
  list(
    uniform_same_signal = list(
      outcome_x = 1.5,
      outcome_u = 1.5,
      outcome_error_scale = 1
    ),
    uniform_broad_signal = list(
      outcome_x = 1,
      outcome_u = 2,
      outcome_error_scale = 0.5
    )
  )
}

generate_uniform_iv_data <- function(
    n,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 0.75,
    treatment_x = 1,
    instrument_x = 0.5,
    outcome_x = 1,
    outcome_u = 2,
    outcome_error_scale = 0.5,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  arguments <- c(
    n = n,
    confounder_strength = confounder_strength,
    treatment_intercept = treatment_intercept,
    instrument_strength = instrument_strength,
    treatment_effect = treatment_effect,
    treatment_x = treatment_x,
    instrument_x = instrument_x,
    outcome_x = outcome_x,
    outcome_u = outcome_u,
    outcome_error_scale = outcome_error_scale
  )
  if (any(!is.finite(arguments)) || n < 1 || n != as.integer(n)) {
    stop("Uniform-DGP arguments must be finite and n a positive integer.")
  }
  if (outcome_u <= 0 || outcome_error_scale <= 0) {
    stop("outcome_u and outcome_error_scale must be positive.")
  }

  X <- stats::runif(n, -1, 1)
  U <- stats::runif(n, -1, 1)
  epsilon0 <- stats::runif(n, -1, 1)
  epsilon1 <- stats::runif(n, -1, 1)
  probability_Z1 <- stats::plogis(instrument_x * X)
  Z <- stats::rbinom(n, 1, probability_Z1)
  treatment_linear_predictor <-
    treatment_intercept + treatment_x * X +
    instrument_strength * Z + confounder_strength * U
  probability_T1 <- stats::plogis(treatment_linear_predictor)
  T_obs <- stats::rbinom(n, 1, probability_T1)
  Y0 <- outcome_x * X + outcome_u * U +
    outcome_error_scale * epsilon0
  Y1 <- treatment_effect + outcome_x * X + outcome_u * U +
    outcome_error_scale * epsilon1
  Y <- ifelse(T_obs == 1L, Y1, Y0)

  data.frame(
    Y = Y,
    T_obs = T_obs,
    Z = Z,
    X = X,
    U = U,
    Y0 = Y0,
    Y1 = Y1,
    probability_Z1 = probability_Z1,
    probability_T1 = probability_T1,
    stringsAsFactors = FALSE
  )
}

uniform_dgp_oracle_pi <- function(
    x,
    z,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_x = 1
) {
  eta <- treatment_intercept + treatment_x * x + instrument_strength * z
  uniform_dgp_mean_expit(eta, confounder_strength)
}

uniform_dgp_oracle_gamma <- function(
    y,
    x,
    treatment,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 0.75,
    treatment_x = 1,
    outcome_x = 1,
    outcome_u = 2,
    outcome_error_scale = 0.5
) {
  treatment <- as.integer(treatment)
  if (!treatment %in% 0:1) stop("treatment must be zero or one.")
  centered_y <- y - treatment_effect * treatment - outcome_x * x
  lower <- pmax(
    -1,
    (centered_y - outcome_error_scale) / outcome_u
  )
  upper <- pmin(
    1,
    (centered_y + outcome_error_scale) / outcome_u
  )
  valid <- is.finite(lower) & is.finite(upper) & upper > lower
  eta0 <- treatment_intercept + treatment_x * x
  eta1 <- eta0 + instrument_strength
  integral0 <- integral1 <- numeric(length(y))
  integral0[valid] <- uniform_dgp_treatment_integral(
    eta0[valid], confounder_strength, lower[valid], upper[valid], treatment
  )
  integral1[valid] <- uniform_dgp_treatment_integral(
    eta1[valid], confounder_strength, lower[valid], upper[valid], treatment
  )
  pi0 <- uniform_dgp_oracle_pi(
    x, 0, confounder_strength, treatment_intercept,
    instrument_strength, treatment_x
  )
  pi1 <- uniform_dgp_oracle_pi(
    x, 1, confounder_strength, treatment_intercept,
    instrument_strength, treatment_x
  )
  probability0 <- if (treatment == 1L) pi0 else 1 - pi0
  probability1 <- if (treatment == 1L) pi1 else 1 - pi1
  output <- rep(NA_real_, length(y))
  stable <- valid & integral0 > 1e-14 & integral1 > 1e-14 &
    probability0 > 0 & probability1 > 0
  output[stable] <-
    (integral1[stable] / probability1[stable]) /
    (integral0[stable] / probability0[stable])
  output
}

add_uniform_oracle_quantities <- function(
    data,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 0.75,
    treatment_x = 1,
    instrument_x = 0.5,
    outcome_x = 1,
    outcome_u = 2,
    outcome_error_scale = 0.5
) {
  data$e_true <- 1 - stats::plogis(instrument_x * data$X)
  data$pi0_true <- uniform_dgp_oracle_pi(
    data$X, 0, confounder_strength, treatment_intercept,
    instrument_strength, treatment_x
  )
  data$pi1_true <- uniform_dgp_oracle_pi(
    data$X, 1, confounder_strength, treatment_intercept,
    instrument_strength, treatment_x
  )
  # Gamma_t is only used on observations with T=t. Potential outcomes are
  # available in simulation, so evaluating at Y^t keeps the stored oracle
  # columns finite outside the observed arm without changing any test moment.
  data$gamma0_oracle <- uniform_dgp_oracle_gamma(
    data$Y0, data$X, 0, confounder_strength, treatment_intercept,
    instrument_strength, treatment_effect, treatment_x,
    outcome_x, outcome_u, outcome_error_scale
  )
  data$gamma1_oracle <- uniform_dgp_oracle_gamma(
    data$Y1, data$X, 1, confounder_strength, treatment_intercept,
    instrument_strength, treatment_effect, treatment_x,
    outcome_x, outcome_u, outcome_error_scale
  )
  required <- c(
    "e_true", "pi0_true", "pi1_true",
    "gamma0_oracle", "gamma1_oracle"
  )
  if (any(!is.finite(as.matrix(data[, required, drop = FALSE])))) {
    stop("Uniform-DGP oracle construction produced non-finite values.")
  }
  data
}
