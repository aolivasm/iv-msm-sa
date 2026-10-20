# Population evaluation of learned threshold functions under an observed-data
# DGP with known nuisance functions.
#
# This module is intentionally diagnostic.  It evaluates the fixed-threshold
# pseudo-bound functional nu_P(q) that is centered by the cross-fitted
# uncentered-EIF score.  Comparing that value with the optimized population
# endpoint separates threshold optimization/target drift from estimation error
# around the data-adaptive fold-average target.

observed_target_decomposition_check_scalar <- function(
    value, name, lower = -Inf, strict_lower = FALSE
) {
  valid_lower <- if (strict_lower) value > lower else value >= lower
  if (length(value) != 1L || !is.finite(value) || !valid_lower) {
    comparator <- if (strict_lower) "greater than" else "at least"
    stop(name, " must be one finite number ", comparator, " ", lower, ".")
  }
  value
}

observed_target_decomposition_check_integer <- function(
    value, name, minimum = 1L
) {
  value <- suppressWarnings(as.integer(value))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(name, " must be one integer of at least ", minimum, ".")
  }
  value
}

observed_target_decomposition_bind_rows <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    row <- as.data.frame(row, stringsAsFactors = FALSE)
    for (column in setdiff(columns, names(row))) row[[column]] <- NA
    row[, columns, drop = FALSE]
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

observed_data_oracle_thresholds <- function(
    newdata,
    treatment,
    specification,
    Lambda,
    epsilon = 0,
    residual_nodes = 801L,
    chunk_size = 256L,
    compatibility_tolerance = 1e-7
) {
  treatment <- as.integer(treatment)
  if (length(treatment) != 1L || !treatment %in% 0:1) {
    stop("treatment must be zero or one.")
  }
  Lambda <- observed_target_decomposition_check_scalar(
    Lambda, "Lambda", lower = 1
  )
  epsilon <- observed_target_decomposition_check_scalar(
    epsilon, "epsilon", lower = 0
  )
  residual_nodes <- observed_target_decomposition_check_integer(
    residual_nodes, "residual_nodes", 101L
  )
  chunk_size <- observed_target_decomposition_check_integer(
    chunk_size, "chunk_size", 1L
  )
  newdata <- as.data.frame(newdata)
  if (!nrow(newdata)) {
    return(data.frame(
      q_lower = numeric(), q_upper = numeric(), tau = numeric()
    ))
  }

  residual_rule <- observed_dgp_trapezoid_rule(
    specification$residual_lower,
    specification$residual_upper,
    residual_nodes
  )
  residual <- as.numeric(residual_rule$nodes)
  base_mass <- observed_dgp_dtruncnorm(
    residual,
    mean = 0,
    lower = specification$residual_lower,
    upper = specification$residual_upper
  ) * residual_rule$weights
  base_mass <- base_mass / sum(base_mass)
  q_lower <- q_upper <- tau_output <- rep(NA_real_, nrow(newdata))

  interpolate_quantile <- function(mass, probability) {
    mass <- pmax(as.numeric(mass), 0)
    total <- sum(mass)
    if (!is.finite(total) || total <= 0) {
      stop("The oracle threshold distribution has zero weighted mass.")
    }
    mass <- mass / total
    cdf_at_node <- cumsum(mass) - mass / 2
    stats::approx(
      x = c(0, cdf_at_node, 1),
      y = c(
        specification$residual_lower,
        residual,
        specification$residual_upper
      ),
      xout = pmin(1, pmax(0, probability)),
      ties = "ordered",
      rule = 2
    )$y
  }

  starts <- seq.int(1L, nrow(newdata), by = chunk_size)
  for (start in starts) {
    index <- start:min(nrow(newdata), start + chunk_size - 1L)
    local <- newdata[index, , drop = FALSE]
    local_n <- nrow(local)
    nuisance <- observed_dgp_nuisances(local, specification)
    pi0 <- as.numeric(nuisance$pi0)
    pi1 <- as.numeric(nuisance$pi1)
    location <- observed_dgp_outcome_location(
      local, treatment, specification
    )
    y <- matrix(
      location + rep(residual, each = local_n),
      nrow = local_n,
      ncol = length(residual)
    )
    expanded_x <- local[
      rep(seq_len(local_n), times = length(residual)), , drop = FALSE
    ]
    gamma <- matrix(
      observed_dgp_gamma(
        as.vector(y), expanded_x, treatment, specification
      ),
      nrow = local_n,
      ncol = length(residual)
    )
    envelope <- observed_dgp_envelopes(
      gamma = as.vector(gamma),
      pi0 = rep(pi0, times = length(residual)),
      pi1 = rep(pi1, times = length(residual)),
      treatment = treatment,
      Lambda0 = Lambda,
      Lambda1 = Lambda,
      epsilon = epsilon
    )
    L <- matrix(envelope$L, nrow = local_n)
    U <- matrix(envelope$U, nrow = local_n)
    D <- U - L
    if (any(D < -compatibility_tolerance)) {
      stop("The oracle threshold envelopes cross.")
    }
    D <- pmax(D, 0)
    mass <- matrix(rep(base_mass, each = local_n), nrow = local_n)
    Lbar <- rowSums(mass * L)
    Ubar <- rowSums(mass * U)
    if (
      any(Lbar > 1 + compatibility_tolerance) ||
        any(Ubar < 1 - compatibility_tolerance)
    ) {
      stop("The oracle envelope means do not bracket one.")
    }
    tau <- (Ubar - 1) / (Ubar - Lbar)
    tau <- pmin(1, pmax(0, tau))
    weighted_mass <- mass * D
    lower_residual <- vapply(seq_len(local_n), function(row) {
      interpolate_quantile(weighted_mass[row, ], 1 - tau[row])
    }, numeric(1L))
    upper_residual <- vapply(seq_len(local_n), function(row) {
      interpolate_quantile(weighted_mass[row, ], tau[row])
    }, numeric(1L))
    q_lower[index] <- location + lower_residual
    q_upper[index] <- location + upper_residual
    tau_output[index] <- tau
  }

  data.frame(
    q_lower = q_lower,
    q_upper = q_upper,
    tau = tau_output,
    stringsAsFactors = FALSE
  )
}

make_observed_data_oracle_threshold_provider_factory <- function(
    specification,
    residual_nodes = 801L,
    chunk_size = 256L,
    compatibility_tolerance = 1e-7
) {
  cache <- new.env(parent = emptyenv())
  factory <- function(
      mode, fold, treatment, gamma_fit, nuisance_fit, x_cols
  ) {
    thresholds <- function(
        train, validation, nuisance_train, nuisance_validation,
        envelope_train, envelope_validation, treatment, Lambda0, Lambda1,
        tau_train, tau_validation, compatibility_action,
        gamma_projection_margin, envelope_smoothing_epsilon, tau_epsilon
    ) {
      if (abs(Lambda0 - Lambda1) > sqrt(.Machine$double.eps)) {
        stop("The observed-DGP oracle provider requires Lambda0 = Lambda1.")
      }
      key <- paste(
        fold, treatment,
        format(Lambda0, digits = 17),
        format(envelope_smoothing_epsilon, digits = 17),
        sep = "::"
      )
      if (!exists(key, envir = cache, inherits = FALSE)) {
        value <- list(
          train = observed_data_oracle_thresholds(
            train, treatment, specification, Lambda0,
            envelope_smoothing_epsilon, residual_nodes, chunk_size,
            compatibility_tolerance
          ),
          validation = observed_data_oracle_thresholds(
            validation, treatment, specification, Lambda0,
            envelope_smoothing_epsilon, residual_nodes, chunk_size,
            compatibility_tolerance
          )
        )
        assign(key, value, envir = cache)
      }
      value <- get(key, envir = cache, inherits = FALSE)
      predict_all <- function(newdata) {
        observed_data_oracle_thresholds(
          newdata, treatment, specification, Lambda0,
          envelope_smoothing_epsilon, residual_nodes, chunk_size,
          compatibility_tolerance
        )
      }
      list(
        q_lower_train = value$train$q_lower,
        q_upper_train = value$train$q_upper,
        q_lower_validation = value$validation$q_lower,
        q_upper_validation = value$validation$q_upper,
        lower_predict = function(newdata) predict_all(newdata)$q_lower,
        upper_predict = function(newdata) predict_all(newdata)$q_upper
      )
    }
    list(thresholds = thresholds)
  }
  attr(factory, "cache") <- cache
  factory
}

evaluate_observed_data_fixed_threshold_arm <- function(
    lower_predict,
    upper_predict,
    treatment,
    specification,
    Lambda,
    epsilon,
    x_rule,
    residual_rule,
    chunk_size = 256L,
    compatibility_tolerance = 1e-7
) {
  if (!is.function(lower_predict) || !is.function(upper_predict)) {
    stop("lower_predict and upper_predict must be functions.")
  }
  treatment <- as.integer(treatment)
  if (length(treatment) != 1L || !treatment %in% 0:1) {
    stop("treatment must be zero or one.")
  }
  chunk_size <- observed_target_decomposition_check_integer(
    chunk_size, "chunk_size"
  )

  x <- as.data.frame(x_rule$x)
  x_weight <- as.numeric(x_rule$weights)
  n_x <- nrow(x)
  residual <- as.numeric(residual_rule$nodes)
  residual_weight <- as.numeric(residual_rule$weights)
  n_y <- length(residual)
  if (length(x_weight) != n_x || length(residual_weight) != n_y) {
    stop("The quadrature rules have inconsistent dimensions.")
  }

  q_lower <- as.numeric(lower_predict(x))
  q_upper <- as.numeric(upper_predict(x))
  if (
    length(q_lower) != n_x || length(q_upper) != n_x ||
      any(!is.finite(q_lower)) || any(!is.finite(q_upper))
  ) {
    stop("A learned threshold returned invalid population-grid predictions.")
  }

  base_density_mass <- observed_dgp_dtruncnorm(
    residual,
    mean = 0,
    lower = specification$residual_lower,
    upper = specification$residual_upper
  ) * residual_weight
  base_density_mass <- base_density_mass / sum(base_density_mass)

  totals <- c(
    fixed_q_lower = 0,
    fixed_q_upper = 0,
    direct_lower = 0,
    direct_upper = 0
  )
  maximum_lower_normalization_error <- 0
  maximum_upper_normalization_error <- 0
  maximum_envelope_mean_violation <- 0
  minimum_tau <- Inf
  maximum_tau <- -Inf

  starts <- seq.int(1L, n_x, by = chunk_size)
  for (start in starts) {
    index <- start:min(n_x, start + chunk_size - 1L)
    local_x <- x[index, , drop = FALSE]
    local_weight <- x_weight[index]
    local_n <- length(index)
    nuisance <- observed_dgp_nuisances(local_x, specification)
    pi0 <- as.numeric(nuisance$pi0)
    pi1 <- as.numeric(nuisance$pi1)
    factual_probability <- if (treatment == 1L) pi0 else 1 - pi0
    missing_probability <- 1 - factual_probability
    location <- observed_dgp_outcome_location(
      local_x, treatment, specification
    )
    y <- matrix(
      location + rep(residual, each = local_n),
      nrow = local_n,
      ncol = n_y
    )
    density_mass <- matrix(
      rep(base_density_mass, each = local_n),
      nrow = local_n,
      ncol = n_y
    )

    repeated_x <- local_x[
      rep(seq_len(local_n), times = n_y), , drop = FALSE
    ]
    gamma <- matrix(
      observed_dgp_gamma(
        as.vector(y), repeated_x, treatment, specification
      ),
      nrow = local_n,
      ncol = n_y
    )
    envelope <- observed_dgp_envelopes(
      gamma = as.vector(gamma),
      pi0 = rep(pi0, times = n_y),
      pi1 = rep(pi1, times = n_y),
      treatment = treatment,
      Lambda0 = Lambda,
      Lambda1 = Lambda,
      epsilon = epsilon
    )
    L <- matrix(envelope$L, nrow = local_n, ncol = n_y)
    U <- matrix(envelope$U, nrow = local_n, ncol = n_y)
    D <- U - L
    if (any(D < -compatibility_tolerance)) {
      stop("The true population envelopes cross on the quadrature grid.")
    }
    D <- pmax(D, 0)

    Lbar <- rowSums(density_mass * L)
    Ubar <- rowSums(density_mass * U)
    maximum_envelope_mean_violation <- max(
      maximum_envelope_mean_violation,
      Lbar - 1,
      1 - Ubar,
      0
    )
    if (maximum_envelope_mean_violation > compatibility_tolerance) {
      stop(
        "The true population envelope means do not bracket one; maximum ",
        "violation = ", signif(maximum_envelope_mean_violation, 6), "."
      )
    }
    Lbar <- pmin(Lbar, 1)
    Ubar <- pmax(Ubar, 1)
    denominator <- Ubar - Lbar
    if (any(denominator <= 0)) {
      stop("A true population envelope has zero or negative width.")
    }
    tau <- (Ubar - 1) / denominator
    if (any(!is.finite(tau)) || any(tau < -compatibility_tolerance) ||
        any(tau > 1 + compatibility_tolerance)) {
      stop("The true population interpolation weight is invalid.")
    }
    tau <- pmin(1, pmax(0, tau))
    minimum_tau <- min(minimum_tau, tau)
    maximum_tau <- max(maximum_tau, tau)

    interpolated <- tau * L + (1 - tau) * U
    observed_mean <- rowSums(density_mass * y)
    interpolated_mean <- rowSums(density_mass * y * interpolated)

    lower_residual <- y - q_lower[index]
    upper_residual <- y - q_upper[index]
    lower_level <- 1 - tau
    upper_level <- tau
    lower_loss <- lower_residual * (
      lower_level - (lower_residual < 0)
    )
    upper_loss <- upper_residual * (
      upper_level - (upper_residual < 0)
    )
    lower_missing <- interpolated_mean -
      rowSums(density_mass * D * lower_loss)
    upper_missing <- interpolated_mean +
      rowSums(density_mass * D * upper_loss)

    lambda_lower <- L + D * (y <= q_lower[index])
    lambda_upper <- U - D * (y <= q_upper[index])
    direct_lower_missing <- rowSums(density_mass * y * lambda_lower)
    direct_upper_missing <- rowSums(density_mass * y * lambda_upper)
    lower_normalization <- rowSums(density_mass * lambda_lower)
    upper_normalization <- rowSums(density_mass * lambda_upper)
    maximum_lower_normalization_error <- max(
      maximum_lower_normalization_error, abs(lower_normalization - 1)
    )
    maximum_upper_normalization_error <- max(
      maximum_upper_normalization_error, abs(upper_normalization - 1)
    )

    totals["fixed_q_lower"] <- totals["fixed_q_lower"] + sum(
      local_weight * (
        factual_probability * observed_mean +
          missing_probability * lower_missing
      )
    )
    totals["fixed_q_upper"] <- totals["fixed_q_upper"] + sum(
      local_weight * (
        factual_probability * observed_mean +
          missing_probability * upper_missing
      )
    )
    totals["direct_lower"] <- totals["direct_lower"] + sum(
      local_weight * (
        factual_probability * observed_mean +
          missing_probability * direct_lower_missing
      )
    )
    totals["direct_upper"] <- totals["direct_upper"] + sum(
      local_weight * (
        factual_probability * observed_mean +
          missing_probability * direct_upper_missing
      )
    )
  }

  data.frame(
    fixed_q_lower = unname(totals["fixed_q_lower"]),
    fixed_q_upper = unname(totals["fixed_q_upper"]),
    direct_lower = unname(totals["direct_lower"]),
    direct_upper = unname(totals["direct_upper"]),
    lower_direct_minus_fixed_q = unname(
      totals["direct_lower"] - totals["fixed_q_lower"]
    ),
    upper_direct_minus_fixed_q = unname(
      totals["direct_upper"] - totals["fixed_q_upper"]
    ),
    maximum_lower_normalization_error =
      maximum_lower_normalization_error,
    maximum_upper_normalization_error =
      maximum_upper_normalization_error,
    maximum_envelope_mean_violation = maximum_envelope_mean_violation,
    minimum_tau = minimum_tau,
    maximum_tau = maximum_tau,
    minimum_q_lower = min(q_lower),
    maximum_q_lower = max(q_lower),
    minimum_q_upper = min(q_upper),
    maximum_q_upper = max(q_upper),
    stringsAsFactors = FALSE
  )
}

make_observed_data_fixed_threshold_evaluator <- function(
    specification,
    Lambda,
    epsilon = 0,
    x_nodes = 120L,
    x_nodes_per_dimension = 7L,
    residual_nodes = 401L,
    chunk_size = 256L,
    compatibility_tolerance = 1e-7
) {
  Lambda <- observed_target_decomposition_check_scalar(
    Lambda, "Lambda", lower = 1
  )
  epsilon <- observed_target_decomposition_check_scalar(
    epsilon, "epsilon", lower = 0
  )
  if (epsilon > 0 && epsilon >= 1 / Lambda) {
    stop("epsilon must be smaller than 1/Lambda.")
  }
  x_nodes <- observed_target_decomposition_check_integer(
    x_nodes, "x_nodes", 2L
  )
  x_nodes_per_dimension <- observed_target_decomposition_check_integer(
    x_nodes_per_dimension, "x_nodes_per_dimension", 2L
  )
  residual_nodes <- observed_target_decomposition_check_integer(
    residual_nodes, "residual_nodes", 3L
  )
  chunk_size <- observed_target_decomposition_check_integer(
    chunk_size, "chunk_size"
  )
  compatibility_tolerance <- observed_target_decomposition_check_scalar(
    compatibility_tolerance,
    "compatibility_tolerance",
    lower = 0,
    strict_lower = TRUE
  )

  x_rule <- observed_dgp_covariate_rule(
    specification,
    x_nodes = x_nodes,
    x_nodes_per_dimension = x_nodes_per_dimension
  )
  residual_rule <- observed_dgp_trapezoid_rule(
    specification$residual_lower,
    specification$residual_upper,
    residual_nodes
  )
  population_result <- compute_observed_data_population_bounds(
    specification = specification,
    Lambda = Lambda,
    epsilon = epsilon,
    x_nodes = x_nodes,
    x_nodes_per_dimension = x_nodes_per_dimension,
    residual_nodes = residual_nodes,
    tolerance = compatibility_tolerance
  )
  population_bounds <- population_result$bounds[
    population_result$bounds$method == "Dual IV--MSM", , drop = FALSE
  ]
  if (nrow(population_bounds) != 1L) {
    stop("The optimized Dual IV--MSM population target lookup failed.")
  }

  evaluate <- function(lower_predict, upper_predict, treatment, mode) {
    output <- evaluate_observed_data_fixed_threshold_arm(
      lower_predict = lower_predict,
      upper_predict = upper_predict,
      treatment = treatment,
      specification = specification,
      Lambda = Lambda,
      epsilon = epsilon,
      x_rule = x_rule,
      residual_rule = residual_rule,
      chunk_size = chunk_size,
      compatibility_tolerance = compatibility_tolerance
    )
    output$population_grid_size <- nrow(x_rule$x)
    output$residual_grid_size <- length(residual_rule$nodes)
    output
  }

  list(
    evaluate = evaluate,
    population_bounds = population_bounds,
    population_diagnostics = population_result$diagnostics,
    settings = data.frame(
      design = specification$name,
      Lambda = Lambda,
      epsilon = epsilon,
      x_nodes = x_nodes,
      x_nodes_per_dimension = x_nodes_per_dimension,
      population_grid_size = nrow(x_rule$x),
      residual_nodes = residual_nodes,
      chunk_size = chunk_size,
      compatibility_tolerance = compatibility_tolerance,
      stringsAsFactors = FALSE
    )
  )
}

combine_fold_fixed_threshold_targets <- function(
    evaluation,
    fold_sizes,
    mode = "rkhs"
) {
  required <- c(
    "mode", "fold", "treatment", "fixed_q_lower", "fixed_q_upper",
    "direct_lower", "direct_upper"
  )
  missing <- setdiff(required, names(evaluation))
  if (length(missing)) {
    stop(
      "The fixed-threshold evaluation is missing: ",
      paste(missing, collapse = ", "), "."
    )
  }
  evaluation <- evaluation[evaluation$mode == mode, , drop = FALSE]
  fold_sizes <- as.numeric(fold_sizes)
  if (!length(fold_sizes) || any(!is.finite(fold_sizes)) ||
      any(fold_sizes <= 0)) {
    stop("fold_sizes must contain positive finite validation-fold sizes.")
  }
  folds <- seq_along(fold_sizes)
  expected <- expand.grid(
    fold = folds,
    treatment = 0:1,
    KEEP.OUT.ATTRS = FALSE
  )
  key <- paste(evaluation$fold, evaluation$treatment, sep = ":")
  expected_key <- paste(expected$fold, expected$treatment, sep = ":")
  if (anyDuplicated(key) || !setequal(key, expected_key)) {
    stop("Expected exactly one fixed-threshold evaluation per fold and arm.")
  }
  evaluation <- evaluation[match(expected_key, key), , drop = FALSE]
  fold_weight <- fold_sizes / sum(fold_sizes)

  weighted_arm <- function(treatment, column) {
    local <- evaluation[evaluation$treatment == treatment, , drop = FALSE]
    local <- local[match(folds, local$fold), , drop = FALSE]
    sum(fold_weight * local[[column]])
  }
  make_endpoints <- function(prefix) {
    mu1_lower <- weighted_arm(1L, paste0(prefix, "lower"))
    mu1_upper <- weighted_arm(1L, paste0(prefix, "upper"))
    mu0_lower <- weighted_arm(0L, paste0(prefix, "lower"))
    mu0_upper <- weighted_arm(0L, paste0(prefix, "upper"))
    c(
      mu1_lower = mu1_lower,
      mu1_upper = mu1_upper,
      mu0_lower = mu0_lower,
      mu0_upper = mu0_upper,
      ate_lower = mu1_lower - mu0_upper,
      ate_upper = mu1_upper - mu0_lower
    )
  }
  fixed_q <- make_endpoints("fixed_q_")
  direct <- make_endpoints("direct_")
  data.frame(
    endpoint = names(fixed_q),
    adaptive_fixed_q_target = unname(fixed_q),
    adaptive_direct_value = unname(direct),
    direct_minus_fixed_q = unname(direct - fixed_q),
    stringsAsFactors = FALSE
  )
}
