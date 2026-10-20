# Comparator estimators for the IV marginal-sensitivity project
#
# This standalone file implements:
# 1. Tan's low-dimensional calibrated (CAL) marginal-sensitivity estimator.
# 2. Binary-IV partial-identification bounds with a known bounded outcome.
#
# The second comparator requires scientifically valid bounds on the potential
# outcomes. The current truncated-normal-noise DGP supplies treatment-specific
# supports Y^0 in [-8,8] and Y^1 in [-6,10].

comparator_check_binary <- function(x, name) {
  values <- sort(unique(stats::na.omit(x)))
  if (length(values) == 0L || !all(values %in% c(0, 1))) {
    stop(name, " must be coded 0/1.")
  }
}

comparator_clip_probability <- function(p, epsilon = 0.01) {
  pmin(1 - epsilon, pmax(epsilon, p))
}

comparator_check_loss <- function(y, q, tau) {
  residual <- y - q
  tau * pmax(residual, 0) + (1 - tau) * pmax(-residual, 0)
}

comparator_design_matrix <- function(data, x_cols, degree = 2L) {
  if (!all(x_cols %in% names(data))) {
    stop("Missing covariates: ", paste(setdiff(x_cols, names(data)), collapse = ", "))
  }
  if (length(x_cols) == 0L) {
    return(matrix(1, nrow = nrow(data), ncol = 1L,
                  dimnames = list(NULL, "(Intercept)")))
  }
  degree <- as.integer(degree)
  if (degree < 1L) stop("degree must be at least 1.")

  pieces <- list(`(Intercept)` = rep(1, nrow(data)))
  for (variable in x_cols) {
    x <- as.numeric(data[[variable]])
    if (any(!is.finite(x))) stop(variable, " contains non-finite values.")
    observed_values <- sort(unique(x))
    variable_degree <- if (
      length(observed_values) <= 2L &&
        all(observed_values %in% c(0, 1))
    ) {
      1L
    } else {
      degree
    }
    for (power in seq_len(variable_degree)) {
      label <- if (power == 1L) variable else paste0(variable, "^", power)
      pieces[[label]] <- x^power
    }
  }
  design <- do.call(cbind, pieces)
  storage.mode(design) <- "double"
  design
}

fit_calibrated_propensity <- function(
    arm, design, probability_epsilon = 0.01, maxit = 2000L
) {
  comparator_check_binary(arm, "arm")
  initial_fit <- suppressWarnings(stats::glm.fit(
    x = design, y = arm, family = stats::binomial()
  ))
  initial <- initial_fit$coefficients
  initial[!is.finite(initial)] <- 0

  objective <- function(coef) {
    eta <- pmax(-35, pmin(35, drop(design %*% coef)))
    mean(arm * exp(-eta) + (1 - arm) * eta)
  }
  fit <- stats::optim(
    initial, objective, method = "BFGS",
    control = list(maxit = maxit, reltol = 1e-10)
  )
  if (fit$convergence != 0L || any(!is.finite(fit$par))) {
    stop("Calibrated propensity optimization failed for one treatment arm.")
  }
  probability_raw <- stats::plogis(drop(design %*% fit$par))
  probability <- comparator_clip_probability(probability_raw, probability_epsilon)
  list(
    coefficients = fit$par,
    probability = probability,
    probability_raw = probability_raw,
    objective = fit$value,
    convergence = fit$convergence,
    truncation_rate = mean(probability != probability_raw)
  )
}

fit_weighted_quantile_nelder_mead <- function(
    y, design, weights, tau, maxit = 5000L
) {
  active <- is.finite(y) & is.finite(weights) & weights > 0
  if (sum(active) <= ncol(design)) {
    stop("Too few observations for weighted quantile regression.")
  }
  ya <- y[active]
  xa <- design[active, , drop = FALSE]
  wa <- weights[active]
  if (requireNamespace("quantreg", quietly = TRUE)) {
    quantile_fit <- tryCatch(
      quantreg::rq.wfit(
        x = xa,
        y = ya,
        tau = tau,
        weights = wa,
        method = "fn"
      ),
      error = identity
    )
    if (!inherits(quantile_fit, "error") &&
        all(is.finite(quantile_fit$coefficients))) {
      coefficients <- as.numeric(quantile_fit$coefficients)
      fitted <- drop(design %*% coefficients)
      return(list(
        coefficients = coefficients,
        fitted = fitted,
        objective = sum(
          wa * comparator_check_loss(ya, drop(xa %*% coefficients), tau)
        ) / sum(wa),
        convergence = 0L,
        counts = NA_integer_,
        solver = "quantreg_fn"
      ))
    }
  }
  initial <- stats::lm.wfit(xa, ya, wa)$coefficients
  initial[!is.finite(initial)] <- 0
  objective <- function(coef) {
    sum(wa * comparator_check_loss(ya, drop(xa %*% coef), tau)) / sum(wa)
  }
  fit <- stats::optim(
    initial, objective, method = "Nelder-Mead",
    control = list(maxit = maxit, reltol = 1e-10)
  )
  if (fit$convergence != 0L || any(!is.finite(fit$par))) {
    warning("Weighted quantile optimization returned its best finite solution without convergence.")
  }
  list(
    coefficients = fit$par,
    fitted = drop(design %*% fit$par),
    objective = fit$value,
    convergence = fit$convergence,
    counts = fit$counts,
    solver = "nelder_mead_fallback"
  )
}

fit_weighted_mean_closed_form <- function(y, design, weights) {
  active <- is.finite(y) & is.finite(weights) & weights > 0
  fit <- stats::lm.wfit(
    design[active, , drop = FALSE], y[active], weights[active]
  )
  coefficients <- fit$coefficients
  coefficients[!is.finite(coefficients)] <- 0
  list(
    coefficients = coefficients,
    fitted = drop(design %*% coefficients),
    rank = fit$rank
  )
}

estimate_tan_cal_arm_endpoint <- function(
    y, arm, design, Lambda, direction = c("lower", "upper"),
    probability_epsilon = 0.01
) {
  direction <- match.arg(direction)
  if (!is.finite(Lambda) || Lambda < 1) {
    stop("Lambda must be finite and at least 1.")
  }
  propensity <- fit_calibrated_propensity(
    arm, design, probability_epsilon = probability_epsilon
  )
  probability <- propensity$probability
  odds_weight <- (1 - probability) / probability
  tau_upper <- Lambda / (Lambda + 1)
  tau <- if (direction == "upper") tau_upper else 1 - tau_upper
  quantile <- fit_weighted_quantile_nelder_mead(
    y, design, arm * odds_weight, tau
  )
  rho <- comparator_check_loss(y, quantile$fitted, tau)
  multiplier <- Lambda - 1 / Lambda
  sign <- if (direction == "upper") 1 else -1
  transformed_y <- y + sign * multiplier * rho
  outcome_mean <- fit_weighted_mean_closed_form(
    transformed_y, design, arm * odds_weight
  )
  phi <- arm / probability * y +
    sign * multiplier * arm * odds_weight * rho -
    (arm / probability - 1) * outcome_mean$fitted
  list(
    estimate = mean(phi),
    standard_error = stats::sd(phi) / sqrt(length(phi)),
    influence_value = phi,
    tau = tau,
    propensity = propensity,
    quantile = quantile,
    outcome_mean = outcome_mean
  )
}

estimate_tan_cal_bounds <- function(
    data, Lambda = 5, x_cols = "X", treatment_col = "T_obs",
    outcome_col = "Y", degree = 2L, probability_epsilon = 0.01,
    confidence_level = 0.90
) {
  required <- c(outcome_col, treatment_col, x_cols)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }
  if (!is.finite(confidence_level) || confidence_level <= 0 ||
      confidence_level >= 1) {
    stop("confidence_level must lie strictly between 0 and 1.")
  }
  y <- as.numeric(data[[outcome_col]])
  treatment <- as.numeric(data[[treatment_col]])
  comparator_check_binary(treatment, treatment_col)
  if (any(!is.finite(y))) stop(outcome_col, " contains non-finite values.")
  design <- comparator_design_matrix(data, x_cols, degree)

  arm_results <- list()
  for (treatment_value in c(0, 1)) {
    arm <- as.numeric(treatment == treatment_value)
    for (direction in c("lower", "upper")) {
      key <- paste0("t", treatment_value, "_", direction)
      arm_results[[key]] <- estimate_tan_cal_arm_endpoint(
        y, arm, design, Lambda, direction, probability_epsilon
      )
    }
  }
  estimates <- c(
    mu1_lower = arm_results$t1_lower$estimate,
    mu1_upper = arm_results$t1_upper$estimate,
    mu0_lower = arm_results$t0_lower$estimate,
    mu0_upper = arm_results$t0_upper$estimate
  )
  estimates <- c(
    estimates,
    ate_lower = unname(estimates["mu1_lower"] - estimates["mu0_upper"]),
    ate_upper = unname(estimates["mu1_upper"] - estimates["mu0_lower"])
  )
  ate_lower_phi <- arm_results$t1_lower$influence_value -
    arm_results$t0_upper$influence_value
  ate_upper_phi <- arm_results$t1_upper$influence_value -
    arm_results$t0_lower$influence_value
  standard_errors <- c(
    mu1_lower = arm_results$t1_lower$standard_error,
    mu1_upper = arm_results$t1_upper$standard_error,
    mu0_lower = arm_results$t0_lower$standard_error,
    mu0_upper = arm_results$t0_upper$standard_error,
    ate_lower = stats::sd(ate_lower_phi) / sqrt(nrow(data)),
    ate_upper = stats::sd(ate_upper_phi) / sqrt(nrow(data))
  )
  alpha <- 1 - confidence_level
  critical <- stats::qnorm(1 - alpha / 2)

  list(
    point_bounds = data.frame(
      method = "Tan CAL", Lambda = Lambda,
      mu1_lower = unname(estimates["mu1_lower"]),
      mu1_upper = unname(estimates["mu1_upper"]),
      mu0_lower = unname(estimates["mu0_lower"]),
      mu0_upper = unname(estimates["mu0_upper"]),
      ate_lower = unname(estimates["ate_lower"]),
      ate_upper = unname(estimates["ate_upper"]),
      stringsAsFactors = FALSE
    ),
    endpoint_standard_errors = data.frame(
      endpoint = names(estimates), estimate = unname(estimates),
      standard_error = unname(standard_errors[names(estimates)])
    ),
    ate_sensitivity_interval = data.frame(
      confidence_level = confidence_level,
      lower = estimates["ate_lower"] - critical * standard_errors["ate_lower"],
      upper = estimates["ate_upper"] + critical * standard_errors["ate_upper"]
    ),
    diagnostics = data.frame(
      treatment = c(0, 1),
      propensity_truncation_rate = c(
        arm_results$t0_lower$propensity$truncation_rate,
        arm_results$t1_lower$propensity$truncation_rate
      ),
      lower_quantile_convergence = c(
        arm_results$t0_lower$quantile$convergence,
        arm_results$t1_lower$quantile$convergence
      ),
      upper_quantile_convergence = c(
        arm_results$t0_upper$quantile$convergence,
        arm_results$t1_upper$quantile$convergence
      )
    ),
    arm_results = arm_results,
    specification = list(
      estimator = "Tan low-dimensional calibrated estimator (CAL)",
      regularized = FALSE,
      Lambda = Lambda,
      x_cols = x_cols,
      polynomial_degree = degree,
      probability_epsilon = probability_epsilon,
      note = "This is CAL from Tan (2024), not high-dimensional Lasso RCAL."
    )
  )
}

fit_instrument_propensity <- function(
    data, instrument_col, x_cols, spline_df = 4L,
    probability_epsilon = 0.01
) {
  z <- as.numeric(data[[instrument_col]])
  comparator_check_binary(z, instrument_col)
  if (length(x_cols) == 0L) {
    probability <- rep(mean(z), nrow(data))
    return(list(probability = probability, fit = NULL, formula = NULL,
                truncation_rate = 0))
  }
  terms <- vapply(x_cols, function(variable) {
    values <- sort(unique(stats::na.omit(data[[variable]])))
    if (length(values) <= 2L) {
      sprintf("`%s`", variable)
    } else {
      sprintf("splines::ns(`%s`, df = %d)", variable, spline_df)
    }
  }, character(1))
  formula <- stats::as.formula(
    paste0("`", instrument_col, "` ~ ", paste(terms, collapse = " + "))
  )
  fit <- stats::glm(formula, data = data, family = stats::binomial())
  probability_raw <- stats::predict(fit, newdata = data, type = "response")
  probability <- comparator_clip_probability(probability_raw, probability_epsilon)
  list(
    probability = probability,
    probability_raw = probability_raw,
    fit = fit,
    formula = formula,
    truncation_rate = mean(probability != probability_raw)
  )
}

standardized_iv_moments <- function(y, treatment, z, probability_z1) {
  moments <- lapply(c(0, 1), function(z_value) {
    probability_z <- if (z_value == 1) probability_z1 else 1 - probability_z1
    weight <- as.numeric(z == z_value) / probability_z
    weight <- weight / mean(weight)
    c(
      treatment_probability = mean(weight * treatment),
      treated_outcome_total = mean(weight * treatment * y),
      untreated_outcome_total = mean(weight * (1 - treatment) * y),
      outcome_mean = mean(weight * y),
      effective_sample_size = sum(weight)^2 / sum(weight^2)
    )
  })
  result <- do.call(rbind, moments)
  rownames(result) <- c("Z=0", "Z=1")
  result
}

comparator_line_value <- function(line, s) line[1] + line[2] * s

resolve_iv_outcome_support <- function(
    outcome_lower = -Inf,
    outcome_upper = Inf,
    outcome_support_by_t = NULL
) {
  if (is.null(outcome_support_by_t)) {
    support <- list(
      `0` = c(outcome_lower, outcome_upper),
      `1` = c(outcome_lower, outcome_upper)
    )
  } else {
    if (!is.list(outcome_support_by_t) ||
        !all(c("0", "1") %in% names(outcome_support_by_t))) {
      stop(
        "outcome_support_by_t must be a named list with elements '0' and '1'."
      )
    }
    support <- outcome_support_by_t[c("0", "1")]
  }

  for (treatment_value in c("0", "1")) {
    interval <- as.numeric(support[[treatment_value]])
    if (length(interval) != 2L || any(is.na(interval))) {
      stop("Each treatment-specific outcome support must contain two numbers.")
    }
    if (interval[1] >= interval[2]) {
      stop("Each outcome-support lower endpoint must be below its upper endpoint.")
    }
    support[[treatment_value]] <- interval
  }
  support
}

binary_iv_response_type_bounds <- function(
    moments,
    outcome_lower = -Inf,
    outcome_upper = Inf,
    outcome_support_by_t = NULL
) {
  support <- resolve_iv_outcome_support(
    outcome_lower = outcome_lower,
    outcome_upper = outcome_upper,
    outcome_support_by_t = outcome_support_by_t
  )
  a0 <- support[["0"]][1]
  b0 <- support[["0"]][2]
  a1 <- support[["1"]][1]
  b1 <- support[["1"]][2]
  p0 <- moments["Z=0", "treatment_probability"]
  p1 <- moments["Z=1", "treatment_probability"]
  m10 <- moments["Z=0", "treated_outcome_total"]
  m11 <- moments["Z=1", "treated_outcome_total"]
  m00 <- moments["Z=0", "untreated_outcome_total"]
  m01 <- moments["Z=1", "untreated_outcome_total"]
  s_lower <- max(0, p0 + p1 - 1)
  s_upper <- min(p0, p1)

  r11_lower_lines <- rbind(
    c(0, a1), c(m10 - b1 * p0, b1), c(m11 - b1 * p1, b1)
  )
  r11_upper_lines <- rbind(
    c(0, b1), c(m10 - a1 * p0, a1), c(m11 - a1 * p1, a1)
  )
  r00_lower_lines <- rbind(
    c(a0 * (1 - p0 - p1), a0),
    c(m00 - b0 * p1, b0),
    c(m01 - b0 * p0, b0)
  )
  r00_upper_lines <- rbind(
    c(b0 * (1 - p0 - p1), b0),
    c(m00 - a0 * p1, a0),
    c(m01 - a0 * p0, a0)
  )
  all_lines <- rbind(
    r11_lower_lines, r11_upper_lines, r00_lower_lines, r00_upper_lines
  )
  candidates <- c(s_lower, s_upper)
  pairs <- utils::combn(seq_len(nrow(all_lines)), 2L)
  for (column in seq_len(ncol(pairs))) {
    first <- all_lines[pairs[1, column], ]
    second <- all_lines[pairs[2, column], ]
    denominator <- first[2] - second[2]
    if (abs(denominator) > 1e-12) {
      crossing <- (second[1] - first[1]) / denominator
      if (crossing >= s_lower - 1e-10 && crossing <= s_upper + 1e-10) {
        candidates <- c(candidates, crossing)
      }
    }
  }
  candidates <- sort(unique(pmin(s_upper, pmax(s_lower, candidates))))
  if (length(candidates) > 1L) {
    candidates <- sort(unique(c(
      candidates, (candidates[-1] + candidates[-length(candidates)]) / 2
    )))
  }

  evaluate <- function(s) {
    q11 <- s
    q10 <- p0 - s
    q01 <- p1 - s
    q00 <- 1 - p0 - p1 + s
    r11_lower <- max(apply(r11_lower_lines, 1L, comparator_line_value, s = s))
    r11_upper <- min(apply(r11_upper_lines, 1L, comparator_line_value, s = s))
    r00_lower <- max(apply(r00_lower_lines, 1L, comparator_line_value, s = s))
    r00_upper <- min(apply(r00_upper_lines, 1L, comparator_line_value, s = s))
    feasible <- r11_lower <= r11_upper + 1e-9 &&
      r00_lower <= r00_upper + 1e-9 && min(q00, q01, q10, q11) >= -1e-9
    if (!feasible) return(NULL)
    mu1_lower <- m10 + m11 - r11_upper + a1 * q00
    mu1_upper <- m10 + m11 - r11_lower + b1 * q00
    mu0_lower <- m00 + m01 - r00_upper + a0 * q11
    mu0_upper <- m00 + m01 - r00_lower + b0 * q11
    c(
      s = s, mu1_lower = mu1_lower, mu1_upper = mu1_upper,
      mu0_lower = mu0_lower, mu0_upper = mu0_upper,
      ate_lower = mu1_lower - mu0_upper,
      ate_upper = mu1_upper - mu0_lower
    )
  }
  values <- lapply(candidates, evaluate)
  values <- values[!vapply(values, is.null, logical(1))]
  if (length(values) == 0L) {
    stop("Observed moments are incompatible with the supplied support and binary-IV assumptions.")
  }
  values <- do.call(rbind, values)
  list(
    bounds = c(
      mu1_lower = min(values[, "mu1_lower"]),
      mu1_upper = max(values[, "mu1_upper"]),
      mu0_lower = min(values[, "mu0_lower"]),
      mu0_upper = max(values[, "mu0_upper"]),
      ate_lower = min(values[, "ate_lower"]),
      ate_upper = max(values[, "ate_upper"])
    ),
    feasible_response_type_grid = values,
    always_taker_probability_range = range(values[, "s"])
  )
}

estimate_iv_partial_identification <- function(
    data, x_cols = "X", treatment_col = "T_obs", instrument_col = "Z",
    outcome_col = "Y", outcome_lower = -Inf, outcome_upper = Inf,
    outcome_support_by_t = NULL,
    instrument_spline_df = 4L, probability_epsilon = 0.01
) {
  required <- c(outcome_col, treatment_col, instrument_col, x_cols)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }
  y <- as.numeric(data[[outcome_col]])
  treatment <- as.numeric(data[[treatment_col]])
  z <- as.numeric(data[[instrument_col]])
  comparator_check_binary(treatment, treatment_col)
  comparator_check_binary(z, instrument_col)
  if (any(!is.finite(y))) stop(outcome_col, " contains non-finite values.")

  instrument_propensity <- fit_instrument_propensity(
    data, instrument_col, x_cols, instrument_spline_df, probability_epsilon
  )
  moments <- standardized_iv_moments(
    y, treatment, z, instrument_propensity$probability
  )
  first_stage <- moments["Z=1", "treatment_probability"] -
    moments["Z=0", "treatment_probability"]

  support <- resolve_iv_outcome_support(
    outcome_lower = outcome_lower,
    outcome_upper = outcome_upper,
    outcome_support_by_t = outcome_support_by_t
  )
  finite_support <- all(is.finite(unlist(support, use.names = FALSE)))
  if (finite_support) {
    for (treatment_value in c(0, 1)) {
      observed <- y[treatment == treatment_value]
      interval <- support[[as.character(treatment_value)]]
      if (min(observed) < interval[1] - 1e-10 ||
          max(observed) > interval[2] + 1e-10) {
        stop(
          "Observed outcomes for treatment ", treatment_value,
          " fall outside the supplied treatment-specific support."
        )
      }
    }
  }
  if (!finite_support) {
    bound_vector <- c(
      mu1_lower = -Inf, mu1_upper = Inf, mu0_lower = -Inf,
      mu0_upper = Inf, ate_lower = -Inf, ate_upper = Inf
    )
    response_type_details <- NULL
    note <- paste(
      "With unbounded potential outcomes, binary-IV independence and exclusion",
      "alone do not give finite population-mean or ATE bounds."
    )
  } else {
    response_type_details <- binary_iv_response_type_bounds(
      moments = moments,
      outcome_support_by_t = support
    )
    bound_vector <- response_type_details$bounds
    note <- paste(
      "Finite response-type bounds require the supplied outcome support,",
      "conditional IV independence, exclusion, and consistency. Monotonicity",
      "is not imposed."
    )
  }
  list(
    point_bounds = data.frame(
      method = "IV partial identification",
      outcome_lower = min(vapply(support, `[`, numeric(1), 1)),
      outcome_upper = max(vapply(support, `[`, numeric(1), 2)),
      outcome_lower_t0 = support[["0"]][1],
      outcome_upper_t0 = support[["0"]][2],
      outcome_lower_t1 = support[["1"]][1],
      outcome_upper_t1 = support[["1"]][2],
      mu1_lower = unname(bound_vector["mu1_lower"]),
      mu1_upper = unname(bound_vector["mu1_upper"]),
      mu0_lower = unname(bound_vector["mu0_lower"]),
      mu0_upper = unname(bound_vector["mu0_upper"]),
      ate_lower = unname(bound_vector["ate_lower"]),
      ate_upper = unname(bound_vector["ate_upper"]),
      stringsAsFactors = FALSE
    ),
    standardized_moments = moments,
    instrument_diagnostics = data.frame(
      treatment_first_stage = first_stage,
      effective_sample_size_z0 = moments["Z=0", "effective_sample_size"],
      effective_sample_size_z1 = moments["Z=1", "effective_sample_size"]
    ),
    instrument_propensity = instrument_propensity,
    response_type_details = response_type_details,
    specification = list(
      finite_support = finite_support,
      outcome_support_by_t = support,
      note = note,
      observed_range_is_not_assumed_support = TRUE
    )
  )
}
