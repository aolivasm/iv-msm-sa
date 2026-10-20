# Sharp Dual-MSM bounds for a finite-support outcome.
#
# This standalone implementation covers Y in {y_1,...,y_K}. It provides:
#   * exact pointwise compatibility and water-filling sharp bounds;
#   * smooth outer endpoints with analytic gradients;
#   * cross-fitted, CV-selected additive spline nuisance regressions;
#   * one-step endpoint and ATE confidence sets, including known sampling
#     weights and household-clustered standard errors.
#
# The binary-outcome implementation is the K=2 special case. No conditional
# quantile or threshold learner is needed.

discrete_check_indicator <- function(x, name) {
  if (anyNA(x) || !all(x %in% c(0, 1))) {
    stop(name, " must contain only zero and one.")
  }
}

discrete_clip_probability <- function(x, epsilon) {
  pmin(pmax(as.numeric(x), epsilon), 1 - epsilon)
}

discrete_normalize_pmf <- function(probability, epsilon = 1e-6) {
  probability <- as.matrix(probability)
  if (any(!is.finite(probability))) stop("A fitted pmf is non-finite.")
  probability <- pmax(probability, epsilon)
  probability / rowSums(probability)
}

discrete_logmeanexp_rows <- function(values, direction, epsilon) {
  values <- as.matrix(values)
  direction <- match.arg(direction, c("max_below", "min_above"))
  n_candidates <- ncol(values)
  if (n_candidates == 1L) {
    return(list(value = values[, 1L], weights = matrix(1, nrow(values), 1L)))
  }
  if (epsilon <= 0) {
    target <- if (direction == "max_below") {
      apply(values, 1L, max)
    } else {
      apply(values, 1L, min)
    }
    tie <- abs(values - target) <= 100 * .Machine$double.eps
    weights <- tie / rowSums(tie)
    return(list(value = target, weights = weights))
  }
  bandwidth <- epsilon / log(n_candidates)
  if (direction == "max_below") {
    center <- apply(values, 1L, max)
    exponential <- exp((values - center) / bandwidth)
    denominator <- rowSums(exponential)
    value <- center + bandwidth * log(denominator / n_candidates)
  } else {
    center <- apply(values, 1L, min)
    exponential <- exp(-(values - center) / bandwidth)
    denominator <- rowSums(exponential)
    value <- center - bandwidth * log(denominator / n_candidates)
  }
  list(value = value, weights = exponential / denominator)
}

discrete_waterfill <- function(lower, upper, support, direction) {
  direction <- match.arg(direction, c("lower", "upper"))
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  support <- as.numeric(support)
  tolerance <- 1e-10
  if (length(lower) != length(upper) || length(lower) != length(support)) {
    stop("lower, upper, and support must have the same length.")
  }
  if (any(lower > upper + tolerance) || sum(lower) > 1 + tolerance ||
      sum(upper) < 1 - tolerance) {
    return(list(feasible = FALSE, value = NA_real_, mass = rep(NA_real_, length(lower))))
  }
  mass <- lower
  remaining <- max(0, 1 - sum(mass))
  order_index <- order(support, decreasing = direction == "upper")
  for (index in order_index) {
    addition <- min(remaining, max(0, upper[index] - mass[index]))
    mass[index] <- mass[index] + addition
    remaining <- remaining - addition
    if (remaining <= tolerance) break
  }
  if (remaining > 1e-7) {
    return(list(feasible = FALSE, value = NA_real_, mass = mass))
  }
  list(
    feasible = TRUE,
    value = sum(support * mass),
    mass = mass / sum(mass)
  )
}

discrete_threshold_coefficients <- function(support, direction) {
  direction <- match.arg(direction, c("lower", "upper"))
  support <- as.numeric(support)
  K <- length(support)
  coefficient_lower <- coefficient_upper <- matrix(0, K, K)
  # Rows index candidate thresholds; columns index outcome categories.
  for (j in seq_len(K)) {
    difference <- support - support[j]
    if (direction == "lower") {
      coefficient_lower[j, support > support[j]] <- difference[support > support[j]]
      coefficient_upper[j, support < support[j]] <- difference[support < support[j]]
    } else {
      coefficient_lower[j, support < support[j]] <- difference[support < support[j]]
      coefficient_upper[j, support > support[j]] <- difference[support > support[j]]
    }
  }
  list(lower = coefficient_lower, upper = coefficient_upper)
}

discrete_arm_endpoint <- function(
    nuisance,
    support,
    Lambda0,
    Lambda1 = Lambda0,
    envelope_epsilon = 0,
    threshold_epsilon = envelope_epsilon,
    compute_gradient = FALSE
) {
  support <- sort(unique(as.numeric(support)))
  K <- length(support)
  required <- c("r0", "r1", "p0", "p1")
  if (!all(required %in% names(nuisance))) stop("Incomplete arm nuisance.")
  r0 <- as.numeric(nuisance$r0)
  r1 <- as.numeric(nuisance$r1)
  p0 <- as.matrix(nuisance$p0)
  p1 <- as.matrix(nuisance$p1)
  n <- length(r0)
  if (length(r1) != n || nrow(p0) != n || nrow(p1) != n ||
      ncol(p0) != K || ncol(p1) != K) {
    stop("Discrete arm nuisance dimensions do not match the support.")
  }
  if (any(r0 <= 0 | r0 >= 1 | r1 <= 0 | r1 >= 1) ||
      any(p0 <= 0) || any(p1 <= 0)) {
    stop("Propensities and fitted category probabilities must be positive and interior.")
  }
  p0 <- p0 / rowSums(p0)
  p1 <- p1 / rowSums(p1)
  finite0 <- is.finite(Lambda0)
  finite1 <- is.finite(Lambda1)
  if ((!finite0 && Lambda0 != Inf) || (!finite1 && Lambda1 != Inf) ||
      Lambda0 < 1 || Lambda1 < 1) {
    stop("Sensitivity parameters must lie in [1,Inf].")
  }

  one_minus_r0 <- 1 - r0
  zero_matrix <- matrix(0, n, K)
  one_matrix <- matrix(1, n, K)
  baseline_lower <- if (finite0) p0 / Lambda0 else zero_matrix
  baseline_upper <- if (finite0) Lambda0 * p0 else one_matrix

  comparison_candidate <- function(multiplier) {
    if (!is.finite(multiplier)) return(matrix(Inf, n, K))
    h <- r1 + (1 - r1) * multiplier
    (p1 * h - r0 * p0) / one_minus_r0
  }
  comparison_derivative <- function(multiplier) {
    h <- r1 + (1 - r1) * multiplier
    list(
      p0 = -r0 / one_minus_r0,
      p1 = h / one_minus_r0,
      r0 = (p1 * h - p0) / one_minus_r0^2,
      r1 = p1 * (1 - multiplier) / one_minus_r0
    )
  }
  comparison_probability_upper <-
    (1 - r1 - r0 * p0 + r1 * p1) / one_minus_r0
  comparison_probability_upper_derivative <- list(
    p0 = -r0 / one_minus_r0,
    p1 = r1 / one_minus_r0,
    r0 = (1 - r1 + r1 * p1 - p0) / one_minus_r0^2,
    r1 = (p1 - 1) / one_minus_r0
  )

  lower_candidates <- list(
    baseline_probability = zero_matrix,
    comparison_probability = comparison_candidate(0)
  )
  lower_derivatives <- list(
    baseline_probability = list(
      p0 = zero_matrix, p1 = zero_matrix,
      r0 = zero_matrix, r1 = zero_matrix
    ),
    comparison_probability = comparison_derivative(0)
  )
  if (finite0) {
    lower_candidates$baseline_sensitivity <- baseline_lower
    lower_derivatives$baseline_sensitivity <- list(
      p0 = 1 / Lambda0, p1 = zero_matrix,
      r0 = zero_matrix, r1 = zero_matrix
    )
  }
  if (finite1) {
    lower_candidates$comparison_sensitivity <- comparison_candidate(1 / Lambda1)
    lower_derivatives$comparison_sensitivity <- comparison_derivative(1 / Lambda1)
  }

  upper_candidates <- list(
    baseline_probability = one_matrix,
    comparison_probability = comparison_probability_upper
  )
  upper_derivatives <- list(
    baseline_probability = list(
      p0 = zero_matrix, p1 = zero_matrix,
      r0 = zero_matrix, r1 = zero_matrix
    ),
    comparison_probability = comparison_probability_upper_derivative
  )
  if (finite0) {
    upper_candidates$baseline_sensitivity <- baseline_upper
    upper_derivatives$baseline_sensitivity <- list(
      p0 = Lambda0, p1 = zero_matrix,
      r0 = zero_matrix, r1 = zero_matrix
    )
  }
  if (finite1) {
    upper_candidates$comparison_sensitivity <- comparison_candidate(Lambda1)
    upper_derivatives$comparison_sensitivity <- comparison_derivative(Lambda1)
  }

  reduce_candidates <- function(candidates, operation) {
    output <- candidates[[1L]]
    if (length(candidates) > 1L) {
      for (index in 2:length(candidates)) {
        output <- operation(output, candidates[[index]])
      }
    }
    matrix(output, nrow = n, ncol = K)
  }
  hard_lower <- reduce_candidates(lower_candidates, pmax)
  hard_upper <- reduce_candidates(upper_candidates, pmin)
  hard_pointwise_failure <- hard_lower > hard_upper + 1e-10
  hard_lower_mass <- rowSums(hard_lower)
  hard_upper_mass <- rowSums(hard_upper)
  hard_feasible <- !apply(hard_pointwise_failure, 1L, any) &
    hard_lower_mass <= 1 + 1e-10 & hard_upper_mass >= 1 - 1e-10

  smooth_candidates <- function(candidates, direction) {
    aggregation <- discrete_logmeanexp_rows(
      do.call(cbind, lapply(candidates, as.vector)), direction,
      envelope_epsilon
    )
    weights <- array(
      aggregation$weights,
      dim = c(n, K, length(candidates)),
      dimnames = list(NULL, NULL, names(candidates))
    )
    list(value = matrix(aggregation$value, n, K), weights = weights)
  }
  lower_smooth <- smooth_candidates(lower_candidates, "max_below")
  upper_smooth <- smooth_candidates(upper_candidates, "min_above")
  lower <- lower_smooth$value
  upper <- upper_smooth$value

  smooth_feasible <- !apply(lower > upper + 1e-10, 1L, any) &
    rowSums(lower) <= 1 + 1e-10 & rowSums(upper) >= 1 - 1e-10

  evaluate_direction <- function(direction) {
    coefficient <- discrete_threshold_coefficients(support, direction)
    candidate <- matrix(NA_real_, n, K)
    for (j in seq_len(K)) {
      candidate[, j] <- support[j] +
        rowSums(sweep(lower, 2L, coefficient$lower[j, ], `*`)) +
        rowSums(sweep(upper, 2L, coefficient$upper[j, ], `*`))
    }
    aggregation <- discrete_logmeanexp_rows(
      candidate,
      if (direction == "lower") "max_below" else "min_above",
      threshold_epsilon
    )
    gradient_lower <- aggregation$weights %*% coefficient$lower
    gradient_upper <- aggregation$weights %*% coefficient$upper
    list(
      value = aggregation$value,
      threshold_weights = aggregation$weights,
      gradient_lower = gradient_lower,
      gradient_upper = gradient_upper,
      candidates = candidate
    )
  }
  lower_value <- evaluate_direction("lower")
  upper_value <- evaluate_direction("upper")

  exact_lower <- exact_upper <- rep(NA_real_, n)
  exact_lower_mass <- exact_upper_mass <- matrix(NA_real_, n, K)
  for (i in which(hard_feasible)) {
    lower_fit <- discrete_waterfill(hard_lower[i, ], hard_upper[i, ], support, "lower")
    upper_fit <- discrete_waterfill(hard_lower[i, ], hard_upper[i, ], support, "upper")
    exact_lower[i] <- lower_fit$value
    exact_upper[i] <- upper_fit$value
    exact_lower_mass[i, ] <- lower_fit$mass
    exact_upper_mass[i, ] <- upper_fit$mass
  }

  make_mean <- function(value) {
    r0 * drop(p0 %*% support) + (1 - r0) * value
  }
  output <- list(
    lower = make_mean(lower_value$value),
    upper = make_mean(upper_value$value),
    exact_lower = make_mean(exact_lower),
    exact_upper = make_mean(exact_upper),
    lower_missing_mean = lower_value$value,
    upper_missing_mean = upper_value$value,
    exact_lower_missing_mean = exact_lower,
    exact_upper_missing_mean = exact_upper,
    hard_lower = hard_lower,
    hard_upper = hard_upper,
    smooth_lower = lower,
    smooth_upper = upper,
    hard_feasible = hard_feasible,
    smooth_feasible = smooth_feasible,
    hard_pointwise_failure = apply(hard_pointwise_failure, 1L, any),
    hard_lower_mass_failure = hard_lower_mass > 1 + 1e-10,
    hard_upper_mass_failure = hard_upper_mass < 1 - 1e-10,
    exact_lower_mass = exact_lower_mass,
    exact_upper_mass = exact_upper_mass,
    envelope_epsilon = envelope_epsilon,
    threshold_epsilon = threshold_epsilon
  )

  if (compute_gradient) {
    combine_candidate_derivatives <- function(weights, derivatives, parameter) {
      output <- matrix(0, n, K)
      for (candidate in names(derivatives)) {
        derivative <- derivatives[[candidate]][[parameter]]
        if (length(derivative) == 1L) derivative <- matrix(derivative, n, K)
        output <- output + weights[, , candidate] * derivative
      }
      output
    }
    smooth_derivative <- list(
      lower = lapply(c("p0", "p1", "r0", "r1"), function(parameter) {
        combine_candidate_derivatives(
          lower_smooth$weights, lower_derivatives, parameter
        )
      }),
      upper = lapply(c("p0", "p1", "r0", "r1"), function(parameter) {
        combine_candidate_derivatives(
          upper_smooth$weights, upper_derivatives, parameter
        )
      })
    )
    names(smooth_derivative$lower) <- names(smooth_derivative$upper) <-
      c("p0", "p1", "r0", "r1")

    build_gradient <- function(endpoint) {
      gradient_L <- endpoint$gradient_lower
      gradient_U <- endpoint$gradient_upper
      dL_dp0 <- smooth_derivative$lower$p0
      dL_dp1 <- smooth_derivative$lower$p1
      dL_dr0 <- smooth_derivative$lower$r0
      dL_dr1 <- smooth_derivative$lower$r1
      dU_dp0 <- smooth_derivative$upper$p0
      dU_dp1 <- smooth_derivative$upper$p1
      dU_dr0 <- smooth_derivative$upper$r0
      dU_dr1 <- smooth_derivative$upper$r1

      dV_dp0 <- gradient_L * dL_dp0 + gradient_U * dU_dp0
      dV_dp1 <- gradient_L * dL_dp1 + gradient_U * dU_dp1
      dV_dr0 <- rowSums(gradient_L * dL_dr0 + gradient_U * dU_dr0)
      dV_dr1 <- rowSums(gradient_L * dL_dr1 + gradient_U * dU_dr1)
      missing_mean <- endpoint$value
      observed_mean <- drop(p0 %*% support)
      list(
        r0 = observed_mean - missing_mean + (1 - r0) * dV_dr0,
        r1 = (1 - r0) * dV_dr1,
        p0 = (1 - r0) * dV_dp0 + outer(r0, support),
        p1 = (1 - r0) * dV_dp1
      )
    }
    output$gradient_lower <- build_gradient(lower_value)
    output$gradient_upper <- build_gradient(upper_value)
  }
  output
}

discrete_make_spline_formula <- function(response, data, x_cols, spline_df) {
  terms <- vapply(x_cols, function(variable) {
    observed <- stats::na.omit(data[[variable]])
    values <- sort(unique(observed))
    if (length(values) <= 2L || spline_df <= 1L) {
      return(sprintf("`%s`", variable))
    }
    effective_df <- min(as.integer(spline_df), length(values) - 1L)
    while (effective_df > 1L) {
      valid <- tryCatch({
        basis <- splines::ns(observed, df = effective_df)
        all(is.finite(basis))
      }, error = function(error) FALSE)
      if (isTRUE(valid)) break
      effective_df <- effective_df - 1L
    }
    if (effective_df <= 1L) {
      sprintf("`%s`", variable)
    } else {
      sprintf("splines::ns(`%s`, df=%d)", variable, effective_df)
    }
  }, character(1))
  stats::as.formula(paste0("`", response, "` ~ ", paste(terms, collapse = " + ")))
}

discrete_cv_fold_id <- function(n, folds, seed) {
  set.seed(as.integer(seed))
  sample(rep(seq_len(folds), length.out = n))
}

discrete_fit_binary_spline_cv <- function(
    data, response, x_cols, weights, df_grid = 1:5,
    cv_folds = 3L, seed = 1L
) {
  outcome <- as.numeric(data[[response]])
  discrete_check_indicator(outcome, response)
  weights <- as.numeric(weights) / mean(weights)
  fold_id <- discrete_cv_fold_id(nrow(data), cv_folds, seed)
  risk <- matrix(NA_real_, length(df_grid), cv_folds)
  for (d in seq_along(df_grid)) for (fold in seq_len(cv_folds)) {
    training <- fold_id != fold
    training_data <- data[training, , drop = FALSE]
    training_data$.analysis_weight <- weights[training]
    fit <- suppressWarnings(stats::glm(
      discrete_make_spline_formula(response, training_data, x_cols, df_grid[d]),
      data = training_data, family = stats::binomial(),
      weights = .analysis_weight, control = stats::glm.control(maxit = 100L)
    ))
    probability <- discrete_clip_probability(
      stats::predict(fit, newdata = data[!training, , drop = FALSE], type = "response"),
      1e-8
    )
    risk[d, fold] <- stats::weighted.mean(
      -(outcome[!training] * log(probability) +
          (1 - outcome[!training]) * log(1 - probability)),
      weights[!training]
    )
  }
  mean_risk <- rowMeans(risk)
  se_risk <- apply(risk, 1L, stats::sd) / sqrt(cv_folds)
  minimum <- which.min(mean_risk)
  eligible <- which(mean_risk <= mean_risk[minimum] + se_risk[minimum])
  selected <- eligible[which.min(df_grid[eligible])]
  fit_data <- data
  fit_data$.analysis_weight <- weights
  fit <- suppressWarnings(stats::glm(
    discrete_make_spline_formula(response, fit_data, x_cols, df_grid[selected]),
    data = fit_data, family = stats::binomial(), weights = .analysis_weight,
    control = stats::glm.control(maxit = 100L)
  ))
  list(
    predict = function(newdata) as.numeric(stats::predict(
      fit, newdata = newdata, type = "response"
    )),
    selected_df = as.integer(df_grid[selected]),
    cv_table = data.frame(df = df_grid, mean_risk = mean_risk, standard_error = se_risk),
    converged = isTRUE(fit$converged)
  )
}

discrete_fit_multinomial_spline_cv <- function(
    data, response, x_cols, support, weights, df_grid = 1:5,
    cv_folds = 3L, decay = 1e-4, seed = 1L
) {
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("Finite-support nuisance estimation requires package 'nnet'.")
  }
  support <- sort(unique(as.numeric(support)))
  local <- data
  local[[response]] <- factor(local[[response]], levels = support)
  if (anyNA(local[[response]])) stop("Observed outcomes fall outside support.")
  weights <- as.numeric(weights) / mean(weights)
  fold_id <- discrete_cv_fold_id(nrow(local), cv_folds, seed)
  risk <- matrix(NA_real_, length(df_grid), cv_folds)
  max_weights <- max(2000L, 5L * length(support) * (length(x_cols) * max(df_grid) + 2L))
  for (d in seq_along(df_grid)) for (fold in seq_len(cv_folds)) {
    training <- fold_id != fold
    training_data <- local[training, , drop = FALSE]
    training_data$.analysis_weight <- weights[training]
    fit <- suppressWarnings(nnet::multinom(
      discrete_make_spline_formula(response, training_data, x_cols, df_grid[d]),
      data = training_data, weights = .analysis_weight,
      decay = decay, trace = FALSE, MaxNWts = max_weights, maxit = 300L
    ))
    probability <- predict(fit, newdata = local[!training, , drop = FALSE], type = "probs")
    if (is.null(dim(probability))) probability <- cbind(1 - probability, probability)
    aligned <- matrix(1e-12, sum(!training), length(support), dimnames = list(NULL, support))
    aligned[, colnames(probability)] <- probability
    aligned <- discrete_normalize_pmf(aligned, 1e-12)
    observed <- match(as.character(local[[response]][!training]), colnames(aligned))
    risk[d, fold] <- stats::weighted.mean(
      -log(aligned[cbind(seq_len(nrow(aligned)), observed)]), weights[!training]
    )
  }
  mean_risk <- rowMeans(risk)
  se_risk <- apply(risk, 1L, stats::sd) / sqrt(cv_folds)
  minimum <- which.min(mean_risk)
  eligible <- which(mean_risk <= mean_risk[minimum] + se_risk[minimum])
  selected <- eligible[which.min(df_grid[eligible])]
  fit_data <- local
  fit_data$.analysis_weight <- weights
  fit <- suppressWarnings(nnet::multinom(
    discrete_make_spline_formula(response, fit_data, x_cols, df_grid[selected]),
    data = fit_data, weights = .analysis_weight, decay = decay, trace = FALSE,
    MaxNWts = max_weights, maxit = 500L
  ))
  list(
    predict = function(newdata) {
      probability <- predict(fit, newdata = newdata, type = "probs")
      if (is.null(dim(probability))) probability <- cbind(1 - probability, probability)
      aligned <- matrix(1e-12, nrow(newdata), length(support), dimnames = list(NULL, support))
      aligned[, colnames(probability)] <- probability
      discrete_normalize_pmf(aligned, 1e-12)
    },
    selected_df = as.integer(df_grid[selected]),
    cv_table = data.frame(df = df_grid, mean_risk = mean_risk, standard_error = se_risk),
    converged = isTRUE(fit$convergence == 0)
  )
}

discrete_make_cluster_folds <- function(data, clusters, K, seed) {
  clusters <- as.character(clusters)
  cluster_data <- data[!duplicated(clusters), c("Z"), drop = FALSE]
  cluster_data$cluster <- clusters[!duplicated(clusters)]
  set.seed(as.integer(seed))
  mapping <- integer(nrow(cluster_data))
  for (z in c(0, 1)) {
    index <- which(cluster_data$Z == z)
    mapping[index] <- sample(rep(seq_len(K), length.out = length(index)))
  }
  mapping[match(clusters, cluster_data$cluster)]
}

fit_discrete_cross_fitted_nuisances <- function(
    data,
    support = sort(unique(data$Y)),
    x_cols,
    folds = 5L,
    fold_id = NULL,
    cluster_id = NULL,
    sampling_weights = NULL,
    df_grid = 1:5,
    inner_folds = 3L,
    overlap_epsilon = 0.05,
    pmf_epsilon = 1e-5,
    multinomial_decay = 1e-4,
    seed = 1L
) {
  required <- c("Y", "T_obs", "Z", x_cols)
  if (!all(required %in% names(data))) stop("Discrete data are missing required columns.")
  discrete_check_indicator(data$T_obs, "T_obs")
  discrete_check_indicator(data$Z, "Z")
  support <- sort(unique(as.numeric(support)))
  if (!all(data$Y %in% support)) stop("Y contains values outside support.")
  n <- nrow(data)
  if (is.null(sampling_weights)) sampling_weights <- rep(1, n)
  if (is.null(cluster_id)) cluster_id <- seq_len(n)
  if (is.null(fold_id)) {
    fold_id <- discrete_make_cluster_folds(data, cluster_id, folds, seed)
  }
  Kcat <- length(support)
  predictions <- list(
    e = rep(NA_real_, n), pi0 = rep(NA_real_, n), pi1 = rep(NA_real_, n),
    p_t0_z0 = matrix(NA_real_, n, Kcat), p_t0_z1 = matrix(NA_real_, n, Kcat),
    p_t1_z0 = matrix(NA_real_, n, Kcat), p_t1_z1 = matrix(NA_real_, n, Kcat)
  )
  diagnostics <- cv_tables <- list()
  counter <- 0L
  for (fold in sort(unique(fold_id))) {
    training_index <- fold_id != fold
    validation_index <- which(fold_id == fold)
    training <- data[training_index, , drop = FALSE]
    validation <- data[validation_index, , drop = FALSE]
    training_weights <- sampling_weights[training_index]
    fits <- list()
    fits$z <- discrete_fit_binary_spline_cv(
      training, "Z", x_cols, training_weights, df_grid, inner_folds,
      seed + 10000L * fold + 1L
    )
    for (z in c(0, 1)) {
      local <- training$Z == z
      fits[[paste0("pi", z)]] <- discrete_fit_binary_spline_cv(
        training[local, , drop = FALSE], "T_obs", x_cols,
        training_weights[local], df_grid, inner_folds,
        seed + 10000L * fold + 10L * z + 2L
      )
      for (t in c(0, 1)) {
        cell <- local & training$T_obs == t
        key <- paste0("p_t", t, "_z", z)
        fits[[key]] <- discrete_fit_multinomial_spline_cv(
          training[cell, , drop = FALSE], "Y", x_cols, support,
          training_weights[cell], df_grid, inner_folds, multinomial_decay,
          seed + 10000L * fold + 100L * t + 10L * z + 3L
        )
      }
    }
    predictions$e[validation_index] <- discrete_clip_probability(
      1 - fits$z$predict(validation), overlap_epsilon
    )
    predictions$pi0[validation_index] <- discrete_clip_probability(
      fits$pi0$predict(validation), overlap_epsilon
    )
    predictions$pi1[validation_index] <- discrete_clip_probability(
      fits$pi1$predict(validation), overlap_epsilon
    )
    for (key in grep("^p_t", names(fits), value = TRUE)) {
      predictions[[key]][validation_index, ] <- discrete_normalize_pmf(
        fits[[key]]$predict(validation), pmf_epsilon
      )
    }
    for (component in names(fits)) {
      counter <- counter + 1L
      diagnostics[[counter]] <- data.frame(
        fold = fold, component = component,
        selected_df = fits[[component]]$selected_df,
        converged = fits[[component]]$converged,
        stringsAsFactors = FALSE
      )
      cv_tables[[counter]] <- cbind(
        data.frame(fold = fold, component = component), fits[[component]]$cv_table
      )
    }
  }
  if (any(!is.finite(unlist(predictions)))) stop("Cross-fitted discrete nuisances are non-finite.")
  list(
    predictions = predictions,
    support = support,
    fold_id = fold_id,
    sampling_weights = sampling_weights,
    cluster_id = cluster_id,
    diagnostics = do.call(rbind, diagnostics),
    cv_tables = do.call(rbind, cv_tables),
    configuration = list(
      df_grid = df_grid, inner_folds = inner_folds,
      overlap_epsilon = overlap_epsilon, pmf_epsilon = pmf_epsilon,
      multinomial_decay = multinomial_decay
    )
  )
}

discrete_arm_nuisance <- function(predictions, treatment) {
  list(
    r0 = if (treatment == 1) predictions$pi0 else 1 - predictions$pi0,
    r1 = if (treatment == 1) predictions$pi1 else 1 - predictions$pi1,
    p0 = predictions[[paste0("p_t", treatment, "_z0")]],
    p1 = predictions[[paste0("p_t", treatment, "_z1")]]
  )
}

discrete_arm_one_step_signal <- function(
    data, predictions, arm, endpoint, treatment, direction, support
) {
  gradient <- if (direction == "lower") endpoint$gradient_lower else endpoint$gradient_upper
  value <- if (direction == "lower") endpoint$lower else endpoint$upper
  observed_arm <- as.numeric(data$T_obs == treatment)
  observed_category <- match(data$Y, support)
  signal <- value
  for (z in c(0, 1)) {
    z_probability <- if (z == 0) predictions$e else 1 - predictions$e
    z_indicator <- as.numeric(data$Z == z)
    r <- arm[[paste0("r", z)]]
    p <- arm[[paste0("p", z)]]
    gp <- gradient[[paste0("p", z)]]
    category_residual <- gp[cbind(seq_len(nrow(data)), observed_category)] -
      rowSums(gp * p)
    signal <- signal + gradient[[paste0("r", z)]] * z_indicator /
      z_probability * (observed_arm - r)
    signal <- signal + z_indicator * observed_arm /
      (z_probability * r) * category_residual
  }
  signal
}

discrete_weighted_cluster_summary <- function(signal, weights, clusters) {
  weights <- as.numeric(weights)
  estimate <- sum(weights * signal) / sum(weights)
  cluster_contribution <- tapply(
    weights * (signal - estimate) / sum(weights), clusters, sum
  )
  G <- length(cluster_contribution)
  standard_error <- sqrt(G / max(1, G - 1) * sum(cluster_contribution^2))
  c(estimate = estimate, standard_error = standard_error, clusters = G)
}

estimate_discrete_dual_bounds <- function(
    data,
    nuisance_fit,
    Lambda0,
    Lambda1 = Lambda0,
    envelope_epsilon = 0.002,
    threshold_epsilon = 0.01,
    confidence_level = 0.90
) {
  predictions <- nuisance_fit$predictions
  support <- nuisance_fit$support
  weights <- nuisance_fit$sampling_weights
  clusters <- nuisance_fit$cluster_id
  arm_fit <- signals <- list()
  for (treatment in c(0, 1)) {
    key <- as.character(treatment)
    arm <- discrete_arm_nuisance(predictions, treatment)
    arm_fit[[key]] <- discrete_arm_endpoint(
      arm, support, Lambda0, Lambda1,
      envelope_epsilon, threshold_epsilon, compute_gradient = TRUE
    )
    if (any(!arm_fit[[key]]$smooth_feasible)) {
      warning(sprintf(
        "Smoothed finite-support envelopes are infeasible for treatment %d at %.3f%% of rows.",
        treatment, 100 * mean(!arm_fit[[key]]$smooth_feasible)
      ))
    }
    signals[[key]] <- list(
      lower = discrete_arm_one_step_signal(
        data, predictions, arm, arm_fit[[key]], treatment, "lower", support
      ),
      upper = discrete_arm_one_step_signal(
        data, predictions, arm, arm_fit[[key]], treatment, "upper", support
      )
    )
  }
  endpoint_signal <- list(
    mu1_lower = signals[["1"]]$lower,
    mu1_upper = signals[["1"]]$upper,
    mu0_lower = signals[["0"]]$lower,
    mu0_upper = signals[["0"]]$upper,
    ate_lower = signals[["1"]]$lower - signals[["0"]]$upper,
    ate_upper = signals[["1"]]$upper - signals[["0"]]$lower
  )
  smooth_plugin_value <- list(
    mu1_lower = arm_fit[["1"]]$lower,
    mu1_upper = arm_fit[["1"]]$upper,
    mu0_lower = arm_fit[["0"]]$lower,
    mu0_upper = arm_fit[["0"]]$upper,
    ate_lower = arm_fit[["1"]]$lower - arm_fit[["0"]]$upper,
    ate_upper = arm_fit[["1"]]$upper - arm_fit[["0"]]$lower
  )
  exact_plugin_value <- list(
    mu1_lower = arm_fit[["1"]]$exact_lower,
    mu1_upper = arm_fit[["1"]]$exact_upper,
    mu0_lower = arm_fit[["0"]]$exact_lower,
    mu0_upper = arm_fit[["0"]]$exact_upper,
    ate_lower = arm_fit[["1"]]$exact_lower - arm_fit[["0"]]$exact_upper,
    ate_upper = arm_fit[["1"]]$exact_upper - arm_fit[["0"]]$exact_lower
  )
  summary <- do.call(rbind, lapply(names(endpoint_signal), function(name) {
    value <- discrete_weighted_cluster_summary(endpoint_signal[[name]], weights, clusters)
    data.frame(
      endpoint = name,
      estimate = value[["estimate"]],
      standard_error = value[["standard_error"]],
      clusters = value[["clusters"]],
      stringsAsFactors = FALSE
    )
  }))
  summary$method <- "smooth_one_step"
  critical <- stats::qnorm(1 - (1 - confidence_level) / 2)
  summary$confidence_level <- confidence_level
  summary$ci_lower <- summary$estimate - critical * summary$standard_error
  summary$ci_upper <- summary$estimate + critical * summary$standard_error

  plugin_summary <- function(values, method) {
    do.call(rbind, lapply(names(values), function(name) {
      value <- values[[name]]
      estimate <- if (any(!is.finite(value))) {
        NA_real_
      } else {
        sum(weights * value) / sum(weights)
      }
      data.frame(
        endpoint = name, estimate = estimate,
        standard_error = NA_real_, clusters = length(unique(clusters)),
        method = method, confidence_level = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }
  estimates <- rbind(
    plugin_summary(exact_plugin_value, "exact_plugin"),
    plugin_summary(smooth_plugin_value, "smooth_plugin"),
    summary
  )

  point <- stats::setNames(summary$estimate, summary$endpoint)
  se <- stats::setNames(summary$standard_error, summary$endpoint)
  confidence_set <- data.frame(
    confidence_level = confidence_level,
    lower = point[["ate_lower"]] - critical * se[["ate_lower"]],
    upper = point[["ate_upper"]] + critical * se[["ate_upper"]]
  )
  diagnostics <- do.call(rbind, lapply(c(0, 1), function(treatment) {
    fit <- arm_fit[[as.character(treatment)]]
    weighted_rate <- function(value) sum(weights * value) / sum(weights)
    data.frame(
      treatment = treatment,
      hard_pointwise_failure_rate = weighted_rate(fit$hard_pointwise_failure),
      hard_lower_mass_failure_rate = weighted_rate(fit$hard_lower_mass_failure),
      hard_upper_mass_failure_rate = weighted_rate(fit$hard_upper_mass_failure),
      hard_any_failure_rate = weighted_rate(!fit$hard_feasible),
      smooth_any_failure_rate = weighted_rate(!fit$smooth_feasible),
      maximum_pointwise_crossing = max(
        0, fit$hard_lower - fit$hard_upper, na.rm = TRUE
      ),
      maximum_lower_mass_excess = max(
        0, rowSums(fit$hard_lower) - 1, na.rm = TRUE
      ),
      maximum_upper_mass_deficit = max(
        0, 1 - rowSums(fit$hard_upper), na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }))
  list(
    estimates = estimates,
    endpoint_summary = summary,
    ate_confidence_set = confidence_set,
    diagnostics = diagnostics,
    arm_fits = arm_fit,
    signals = endpoint_signal,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    envelope_epsilon = envelope_epsilon,
    threshold_epsilon = threshold_epsilon
  )
}
