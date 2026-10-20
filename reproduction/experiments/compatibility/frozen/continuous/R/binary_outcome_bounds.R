# Binary-outcome sharp bounds under the dual marginal sensitivity model
#
# This file is a standalone implementation for binary outcomes.  It contains:
#   * a binary potential-outcome DGP with the same IV/treatment mechanism used
#     by the continuous-outcome simulations;
#   * numerical population references;
#   * cross-fitted Bernoulli nuisance regressions;
#   * exact plug-in sharp bounds;
#   * conservative smoothed outer bounds and their one-step estimators; and
#   * binary IV-only and standard-MSM plug-in comparators.
#
# No RKHS density-ratio or conditional-quantile learner is needed: for binary
# outcomes the density ratio is determined by two Bernoulli regressions and the
# extremization is algebraic.

binary_clip_probability <- function(value, epsilon = 0.05) {
  pmin(1 - epsilon, pmax(epsilon, value))
}

binary_check_indicator <- function(value, name) {
  observed <- sort(unique(stats::na.omit(value)))
  if (length(observed) == 0L || !all(observed %in% c(0, 1))) {
    stop(name, " must be coded 0/1.")
  }
}

binary_rtnorm <- function(n, lower = -2, upper = 2) {
  lower_probability <- stats::pnorm(lower)
  upper_probability <- stats::pnorm(upper)
  stats::qnorm(stats::runif(n, lower_probability, upper_probability))
}

binary_dtnorm <- function(value, lower = -2, upper = 2) {
  normalizer <- stats::pnorm(upper) - stats::pnorm(lower)
  ifelse(
    value >= lower & value <= upper,
    stats::dnorm(value) / normalizer,
    0
  )
}

binary_gauss_legendre_rule <- function(n_nodes = 80L, lower = -2, upper = 2) {
  n_nodes <- as.integer(n_nodes)
  if (n_nodes < 2L) stop("n_nodes must be at least two.")
  index <- seq_len(n_nodes - 1L)
  beta <- index / sqrt(4 * index^2 - 1)
  jacobi <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  jacobi[cbind(index, index + 1L)] <- beta
  jacobi[cbind(index + 1L, index)] <- beta
  decomposition <- eigen(jacobi, symmetric = TRUE)
  ordering <- order(decomposition$values)
  nodes <- decomposition$values[ordering]
  weights <- 2 * decomposition$vectors[1L, ordering]^2
  list(
    nodes = (upper + lower) / 2 + (upper - lower) * nodes / 2,
    weights = (upper - lower) * weights / 2
  )
}

generate_binary_iv_data <- function(
    n = 1500L,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    outcome_intercept = -0.5,
    treatment_log_odds_effect = 0.75,
    outcome_x_strength = 0.75,
    outcome_u_strength = 0.70,
    seed = NULL
) {
  n <- as.integer(n)
  if (n < 1L) stop("n must be positive.")
  if (!is.null(seed)) set.seed(seed)

  X <- binary_rtnorm(n)
  U <- binary_rtnorm(n)
  probability_Z1 <- stats::plogis(0.5 * X)
  Z <- stats::rbinom(n, 1, probability_Z1)
  probability_T1 <- stats::plogis(
    treatment_intercept + X + instrument_strength * Z +
      confounder_strength * U
  )
  T_obs <- stats::rbinom(n, 1, probability_T1)

  probability_Y0 <- stats::plogis(
    outcome_intercept + outcome_x_strength * X + outcome_u_strength * U
  )
  probability_Y1 <- stats::plogis(
    outcome_intercept + treatment_log_odds_effect +
      outcome_x_strength * X + outcome_u_strength * U
  )
  Y0 <- stats::rbinom(n, 1, probability_Y0)
  Y1 <- stats::rbinom(n, 1, probability_Y1)
  Y <- ifelse(T_obs == 1, Y1, Y0)

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
    probability_Y0 = probability_Y0,
    probability_Y1 = probability_Y1
  )
}

binary_make_tz_folds <- function(data, K = 5L, seed = 1L) {
  binary_check_indicator(data$T_obs, "T_obs")
  binary_check_indicator(data$Z, "Z")
  K <- as.integer(K)
  if (K < 2L) stop("K must be at least two.")
  set.seed(seed)
  strata <- interaction(data$T_obs, data$Z, drop = TRUE)
  fold_id <- integer(nrow(data))
  for (stratum in levels(strata)) {
    index <- which(strata == stratum)
    if (length(index) < K) {
      stop("A (T,Z) stratum contains fewer observations than folds.")
    }
    fold_id[index] <- sample(rep(seq_len(K), length.out = length(index)))
  }
  fold_id
}

binary_smooth_max_below <- function(first, second, epsilon) {
  if (epsilon == 0) {
    tie <- abs(first - second) <= 10 * .Machine$double.eps
    return(list(
      value = pmax(first, second),
      weight_first = ifelse(tie, 0.5, as.numeric(first > second)),
      weight_second = ifelse(tie, 0.5, as.numeric(second > first))
    ))
  }
  bandwidth <- epsilon / log(2)
  maximum <- pmax(first, second)
  exponential_first <- exp((first - maximum) / bandwidth)
  exponential_second <- exp((second - maximum) / bandwidth)
  denominator <- exponential_first + exponential_second
  list(
    value = maximum + bandwidth * log(denominator / 2),
    weight_first = exponential_first / denominator,
    weight_second = exponential_second / denominator
  )
}

binary_smooth_min_above <- function(first, second, epsilon) {
  if (epsilon == 0) {
    tie <- abs(first - second) <= 10 * .Machine$double.eps
    return(list(
      value = pmin(first, second),
      weight_first = ifelse(tie, 0.5, as.numeric(first < second)),
      weight_second = ifelse(tie, 0.5, as.numeric(second < first))
    ))
  }
  bandwidth <- epsilon / log(2)
  minimum <- pmin(first, second)
  exponential_first <- exp(-(first - minimum) / bandwidth)
  exponential_second <- exp(-(second - minimum) / bandwidth)
  denominator <- exponential_first + exponential_second
  list(
    value = minimum - bandwidth * log(denominator / 2),
    weight_first = exponential_first / denominator,
    weight_second = exponential_second / denominator
  )
}

binary_recycle_arm_nuisance <- function(nuisance) {
  required <- c("r0", "r1", "a0", "a1")
  if (!all(required %in% names(nuisance))) {
    stop("Arm nuisance must contain r0, r1, a0, and a1.")
  }
  target_length <- max(vapply(nuisance[required], length, integer(1)))
  output <- lapply(nuisance[required], rep_len, length.out = target_length)
  if (any(!is.finite(unlist(output, use.names = FALSE)))) {
    stop("Arm nuisance contains non-finite values.")
  }
  output
}

binary_gamma_compatible_interval <- function(r0, r1, Lambda0, Lambda1) {
  list(
    lower = ((1 - r0) / Lambda0 + r0) /
      ((1 - r1) * Lambda1 + r1),
    upper = (Lambda0 * (1 - r0) + r0) /
      ((1 - r1) / Lambda1 + r1)
  )
}

binary_arm_endpoint <- function(
    nuisance,
    Lambda0,
    Lambda1,
    epsilon = 0,
    compute_gradient = FALSE
) {
  if (!is.finite(Lambda0) || !is.finite(Lambda1) ||
      Lambda0 < 1 || Lambda1 < 1) {
    stop("Lambda0 and Lambda1 must be finite and at least one.")
  }
  if (!is.finite(epsilon) || epsilon < 0 ||
      (epsilon > 0 && epsilon >= 1 / Lambda0)) {
    stop("epsilon must be nonnegative and smaller than 1/Lambda0.")
  }
  nuisance <- binary_recycle_arm_nuisance(nuisance)
  r0 <- nuisance$r0
  r1 <- nuisance$r1
  a0 <- nuisance$a0
  a1 <- nuisance$a1
  n <- length(r0)

  if (any(r0 <= 0 | r0 >= 1 | r1 <= 0 | r1 >= 1 |
          a0 <= 0 | a0 >= 1 | a1 <= 0 | a1 >= 1)) {
    stop("Arm propensities and outcome regressions must lie strictly in (0,1).")
  }

  gradient_names <- c("r0", "r1", "a0", "a1")
  zero_gradient <- matrix(
    0, nrow = n, ncol = length(gradient_names),
    dimnames = list(NULL, gradient_names)
  )

  gamma_one <- a1 / a0
  gamma_zero <- (1 - a1) / (1 - a0)
  gamma_one_gradient <- zero_gradient
  gamma_one_gradient[, "a0"] <- -gamma_one / a0
  gamma_one_gradient[, "a1"] <- 1 / a0
  gamma_zero_gradient <- zero_gradient
  gamma_zero_gradient[, "a0"] <- gamma_zero / (1 - a0)
  gamma_zero_gradient[, "a1"] <- -1 / (1 - a0)

  construct_envelope <- function(gamma, gamma_gradient) {
    candidate <- function(c_value) {
      h_value <- r1 + (1 - r1) * c_value
      value <- (gamma * h_value - r0) / (1 - r0)
      derivative <- gamma_gradient * (h_value / (1 - r0))
      derivative[, "r0"] <- derivative[, "r0"] +
        (gamma * h_value - 1) / (1 - r0)^2
      derivative[, "r1"] <- derivative[, "r1"] +
        gamma * (1 - c_value) / (1 - r0)
      list(value = value, derivative = derivative)
    }

    if (Lambda0 == 1 && Lambda1 == 1) {
      return(list(
        L = rep(1, n), U = rep(1, n),
        derivative_L = zero_gradient,
        derivative_U = zero_gradient,
        candidate_lower = rep(1, n),
        candidate_upper = rep(1, n)
      ))
    }

    candidate_lower <- candidate(1 / Lambda1)
    candidate_upper <- candidate(Lambda1)
    lower <- binary_smooth_max_below(
      rep(1 / Lambda0, n), candidate_lower$value, epsilon
    )
    upper <- binary_smooth_min_above(
      rep(Lambda0, n), candidate_upper$value, epsilon
    )
    list(
      L = lower$value,
      U = upper$value,
      derivative_L = candidate_lower$derivative * lower$weight_second,
      derivative_U = candidate_upper$derivative * upper$weight_second,
      candidate_lower = candidate_lower$value,
      candidate_upper = candidate_upper$value
    )
  }

  envelope_one <- construct_envelope(gamma_one, gamma_one_gradient)
  envelope_zero <- construct_envelope(gamma_zero, gamma_zero_gradient)

  lower_first <- a0 * envelope_one$L
  lower_second <- 1 - (1 - a0) * envelope_zero$U
  upper_first <- a0 * envelope_one$U
  upper_second <- 1 - (1 - a0) * envelope_zero$L

  lower_pair <- binary_smooth_max_below(
    lower_first, lower_second, epsilon
  )
  upper_pair <- binary_smooth_min_above(
    upper_first, upper_second, epsilon
  )

  derivative_lower_first <- envelope_one$derivative_L * a0
  derivative_lower_first[, "a0"] <-
    derivative_lower_first[, "a0"] + envelope_one$L
  derivative_lower_second <- -envelope_zero$derivative_U * (1 - a0)
  derivative_lower_second[, "a0"] <-
    derivative_lower_second[, "a0"] + envelope_zero$U
  derivative_upper_first <- envelope_one$derivative_U * a0
  derivative_upper_first[, "a0"] <-
    derivative_upper_first[, "a0"] + envelope_one$U
  derivative_upper_second <- -envelope_zero$derivative_L * (1 - a0)
  derivative_upper_second[, "a0"] <-
    derivative_upper_second[, "a0"] + envelope_zero$L

  derivative_theta_lower <-
    derivative_lower_first * lower_pair$weight_first +
    derivative_lower_second * lower_pair$weight_second
  derivative_theta_upper <-
    derivative_upper_first * upper_pair$weight_first +
    derivative_upper_second * upper_pair$weight_second

  make_mean <- function(theta, derivative_theta) {
    value <- r0 * a0 + (1 - r0) * theta
    derivative <- derivative_theta * (1 - r0)
    derivative[, "r0"] <- derivative[, "r0"] + a0 - theta
    derivative[, "a0"] <- derivative[, "a0"] + r0
    list(value = value, derivative = derivative)
  }
  lower_mean <- make_mean(lower_pair$value, derivative_theta_lower)
  upper_mean <- make_mean(upper_pair$value, derivative_theta_upper)

  interval <- binary_gamma_compatible_interval(r0, r1, Lambda0, Lambda1)
  comparison_tolerance <- 100 * .Machine$double.eps
  lower_tolerance <- comparison_tolerance * pmax(1, abs(interval$lower))
  upper_tolerance <- comparison_tolerance * pmax(1, abs(interval$upper))
  gamma_violation_one <- gamma_one < interval$lower - lower_tolerance |
    gamma_one > interval$upper + upper_tolerance
  gamma_violation_zero <- gamma_zero < interval$lower - lower_tolerance |
    gamma_zero > interval$upper + upper_tolerance

  output <- list(
    lower = lower_mean$value,
    upper = upper_mean$value,
    theta_lower = lower_pair$value,
    theta_upper = upper_pair$value,
    gamma_one = gamma_one,
    gamma_zero = gamma_zero,
    compatible_gamma_lower = interval$lower,
    compatible_gamma_upper = interval$upper,
    L_one = envelope_one$L,
    U_one = envelope_one$U,
    L_zero = envelope_zero$L,
    U_zero = envelope_zero$U,
    candidate_lower_one = envelope_one$candidate_lower,
    candidate_upper_one = envelope_one$candidate_upper,
    candidate_lower_zero = envelope_zero$candidate_lower,
    candidate_upper_zero = envelope_zero$candidate_upper,
    gamma_violation = gamma_violation_one | gamma_violation_zero,
    envelope_crossing = envelope_one$L >
      envelope_one$U + comparison_tolerance * pmax(1, abs(envelope_one$U)) |
      envelope_zero$L >
        envelope_zero$U + comparison_tolerance * pmax(1, abs(envelope_zero$U)),
    theta_crossing = lower_pair$value >
      upper_pair$value + comparison_tolerance * pmax(1, abs(upper_pair$value)),
    epsilon = epsilon
  )
  if (compute_gradient) {
    output$gradient_lower <- lower_mean$derivative
    output$gradient_upper <- upper_mean$derivative
  }
  output
}

# Robins--Manski natural IV bounds for a binary potential outcome.  These are
# the Lambda -> Inf limit of the binary Dual-MSM bounds.  Smoothing is shifted
# outward: the lower endpoint is no larger than the hard max and the upper
# endpoint is no smaller than the hard min.
binary_iv_arm_endpoint <- function(
    nuisance,
    epsilon = 0,
    compute_gradient = FALSE
) {
  if (!is.finite(epsilon) || epsilon < 0) {
    stop("epsilon must be one nonnegative finite number.")
  }
  nuisance <- binary_recycle_arm_nuisance(nuisance)
  r0 <- nuisance$r0
  r1 <- nuisance$r1
  a0 <- nuisance$a0
  a1 <- nuisance$a1
  n <- length(r0)
  if (any(r0 <= 0 | r0 >= 1 | r1 <= 0 | r1 >= 1 |
          a0 <= 0 | a0 >= 1 | a1 <= 0 | a1 >= 1)) {
    stop("Arm propensities and outcome regressions must lie strictly in (0,1).")
  }

  lower0 <- r0 * a0
  lower1 <- r1 * a1
  upper0 <- 1 - r0 * (1 - a0)
  upper1 <- 1 - r1 * (1 - a1)
  lower <- binary_smooth_max_below(lower0, lower1, epsilon)
  upper <- binary_smooth_min_above(upper0, upper1, epsilon)
  output <- list(
    lower = lower$value,
    upper = upper$value,
    hard_lower = pmax(lower0, lower1),
    hard_upper = pmin(upper0, upper1),
    epsilon = epsilon
  )
  if (compute_gradient) {
    gradient_names <- c("r0", "r1", "a0", "a1")
    gradient_lower <- gradient_upper <- matrix(
      0, nrow = n, ncol = length(gradient_names),
      dimnames = list(NULL, gradient_names)
    )
    gradient_lower[, "r0"] <- lower$weight_first * a0
    gradient_lower[, "a0"] <- lower$weight_first * r0
    gradient_lower[, "r1"] <- lower$weight_second * a1
    gradient_lower[, "a1"] <- lower$weight_second * r1
    gradient_upper[, "r0"] <- upper$weight_first * (a0 - 1)
    gradient_upper[, "a0"] <- upper$weight_first * r0
    gradient_upper[, "r1"] <- upper$weight_second * (a1 - 1)
    gradient_upper[, "a1"] <- upper$weight_second * r1
    output$gradient_lower <- gradient_lower
    output$gradient_upper <- gradient_upper
  }
  output
}

binary_make_spline_formula <- function(
    response, data, x_cols, spline_df
) {
  terms <- unlist(lapply(x_cols, function(variable) {
    values <- sort(unique(stats::na.omit(data[[variable]])))
    if (length(values) <= 1L) {
      character()
    } else if (length(values) <= 2L) {
      sprintf("`%s`", variable)
    } else {
      sprintf("splines::ns(`%s`, df=%d)", variable, as.integer(spline_df))
    }
  }), use.names = FALSE)
  stats::as.formula(
    paste0(
      "`", response, "` ~ ",
      if (length(terms)) paste(terms, collapse = " + ") else "1"
    )
  )
}

binary_fit_regression <- function(data, response, x_cols, spline_df) {
  binary_fit_regression_with_learner(
    data = data,
    response = response,
    x_cols = x_cols,
    spline_df = spline_df,
    learner = "glm_spline",
    seed = 1L
  )
}

binary_make_gam_formula <- function(response, data, x_cols, basis_dimension) {
  terms <- unlist(lapply(x_cols, function(variable) {
    values <- sort(unique(stats::na.omit(data[[variable]])))
    if (length(values) <= 1L) {
      character()
    } else if (length(values) <= 2L) {
      sprintf("`%s`", variable)
    } else {
      local_k <- min(as.integer(basis_dimension), length(values) - 1L)
      local_k <- max(3L, local_k)
      sprintf("s(`%s`, k=%d, bs='tp')", variable, local_k)
    }
  }), use.names = FALSE)
  stats::as.formula(
    paste0(
      "`", response, "` ~ ",
      if (length(terms)) paste(terms, collapse = " + ") else "1"
    )
  )
}

binary_make_interaction_gam_formula <- function(
    response,
    data,
    x_cols,
    basis_dimension,
    interaction_basis_dimension = 3L
) {
  unique_counts <- vapply(x_cols, function(variable) {
    length(unique(stats::na.omit(data[[variable]])))
  }, integer(1))
  active <- x_cols[unique_counts > 1L]
  if (!length(active)) {
    return(stats::as.formula(sprintf("`%s` ~ 1", response)))
  }
  main_terms <- vapply(active, function(variable) {
    count <- unique_counts[[variable]]
    if (count <= 2L) return(sprintf("`%s`", variable))
    local_k <- min(as.integer(basis_dimension), count - 1L)
    local_k <- max(3L, local_k)
    sprintf("s(`%s`, k=%d, bs='tp')", variable, local_k)
  }, character(1))
  interaction_terms <- character()
  if (length(active) > 1L) {
    pairs <- utils::combn(active, 2L, simplify = FALSE)
    interaction_terms <- vapply(pairs, function(pair) {
      first_count <- unique_counts[[pair[1L]]]
      second_count <- unique_counts[[pair[2L]]]
      if (first_count <= 2L && second_count <= 2L) {
        return(sprintf("`%s`:`%s`", pair[1L], pair[2L]))
      }
      if (first_count <= 2L || second_count <= 2L) {
        binary_variable <- pair[if (first_count <= 2L) 1L else 2L]
        continuous_variable <- pair[if (first_count <= 2L) 2L else 1L]
        continuous_count <- unique_counts[[continuous_variable]]
        local_k <- min(
          as.integer(interaction_basis_dimension), continuous_count - 1L
        )
        local_k <- max(3L, local_k)
        return(sprintf(
          "s(`%s`, by=`%s`, k=%d, bs='tp')",
          continuous_variable, binary_variable, local_k
        ))
      }
      first_k <- max(3L, min(
        as.integer(interaction_basis_dimension), first_count - 1L
      ))
      second_k <- max(3L, min(
        as.integer(interaction_basis_dimension), second_count - 1L
      ))
      sprintf(
        "ti(`%s`, `%s`, k=c(%d,%d), bs=c('tp','tp'))",
        pair[1L], pair[2L], first_k, second_k
      )
    }, character(1))
  }
  stats::as.formula(paste0(
    "`", response, "` ~ ",
    paste(c(main_terms, interaction_terms), collapse = " + ")
  ))
}

binary_pairwise_interaction_design <- function(data) {
  data <- as.data.frame(data)
  if (!ncol(data)) stop("At least one predictor is required.")
  design <- as.matrix(data)
  storage.mode(design) <- "double"
  if (any(!is.finite(design))) {
    stop("The pairwise-interaction design contains non-finite values.")
  }
  if (ncol(design) == 1L) {
    squared <- matrix(
      design[, 1L]^2,
      ncol = 1L,
      dimnames = list(NULL, paste0(colnames(design), "_squared"))
    )
    return(cbind(design, squared))
  }
  pairs <- utils::combn(seq_len(ncol(design)), 2L)
  interactions <- vapply(seq_len(ncol(pairs)), function(index) {
    design[, pairs[1L, index]] * design[, pairs[2L, index]]
  }, numeric(nrow(design)))
  if (is.null(dim(interactions))) interactions <- matrix(interactions, ncol = 1L)
  colnames(interactions) <- vapply(seq_len(ncol(pairs)), function(index) {
    paste0(
      colnames(design)[pairs[1L, index]], "_X_",
      colnames(design)[pairs[2L, index]]
    )
  }, character(1))
  cbind(design, interactions)
}

binary_super_learner_numeric_design <- function(data, x_cols) {
  predictor_data <- data[, x_cols, drop = FALSE]
  formula <- stats::reformulate(x_cols, response = NULL, intercept = FALSE)
  terms <- stats::terms(formula, data = predictor_data)
  frame <- stats::model.frame(
    terms, data = predictor_data, na.action = stats::na.fail,
    drop.unused.levels = TRUE
  )
  xlevels <- stats::.getXlevels(terms, frame)
  contrasts <- lapply(frame[vapply(frame, is.factor, logical(1))], contrasts)
  design <- stats::model.matrix(terms, data = frame, contrasts.arg = contrasts)
  storage.mode(design) <- "double"
  if (!ncol(design) || any(!is.finite(design))) {
    stop("The Super Learner predictor design is empty or non-finite.")
  }
  training_columns <- colnames(design)
  transform <- function(newdata) {
    new_frame <- stats::model.frame(
      terms, data = newdata[, x_cols, drop = FALSE],
      na.action = stats::na.fail, xlev = xlevels
    )
    new_design <- stats::model.matrix(
      terms, data = new_frame, contrasts.arg = contrasts
    )
    missing <- setdiff(training_columns, colnames(new_design))
    if (length(missing)) {
      new_design <- cbind(
        new_design,
        matrix(
          0, nrow = nrow(new_design), ncol = length(missing),
          dimnames = list(NULL, missing)
        )
      )
    }
    extra <- setdiff(colnames(new_design), training_columns)
    if (length(extra)) {
      stop(
        "New Super Learner data generated unknown design columns: ",
        paste(extra, collapse = ", "), "."
      )
    }
    new_design <- new_design[, training_columns, drop = FALSE]
    storage.mode(new_design) <- "double"
    if (any(!is.finite(new_design))) {
      stop("New Super Learner predictor values are non-finite.")
    }
    new_design
  }
  list(training = design, transform = transform, columns = training_columns)
}

# SuperLearner-compatible wrappers.  Their names are deliberately specific to
# this implementation so the library and fitted weights are fully auditable.
SL.binary_glm_spline <- function(
    Y, X, newX, family, obsWeights, spline_df = 3L, ...
) {
  training <- as.data.frame(X)
  training$.outcome <- as.numeric(Y)
  formula <- binary_make_spline_formula(
    ".outcome", training, setdiff(names(training), ".outcome"), spline_df
  )
  environment(formula) <- environment()
  fit <- suppressWarnings(stats::glm(
    formula,
    data = training,
    family = family,
    weights = obsWeights,
    control = stats::glm.control(maxit = 100L)
  ))
  prediction <- as.numeric(stats::predict(
    fit, newdata = as.data.frame(newX), type = "response"
  ))
  output <- list(pred = prediction, fit = list(object = fit))
  class(output$fit) <- "SL.binary_glm_spline"
  output
}

predict.SL.binary_glm_spline <- function(object, newdata, ...) {
  as.numeric(stats::predict(
    object$object, newdata = as.data.frame(newdata), type = "response"
  ))
}

SL.binary_glmnet_interactions <- function(
    Y, X, newX, family, obsWeights, id,
    nfolds = 5L, nlambda = 75L, ...
) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("SL.binary_glmnet_interactions requires package 'glmnet'.")
  }
  design <- binary_pairwise_interaction_design(X)
  new_design <- binary_pairwise_interaction_design(newX)
  minimum_class_size <- min(tabulate(as.integer(Y) + 1L, nbins = 2L))
  local_folds <- min(as.integer(nfolds), minimum_class_size)
  if (local_folds < 4L) {
    stop("At least four observations per outcome class are required by glmnet.")
  }
  fit <- suppressWarnings(glmnet::cv.glmnet(
    x = design,
    y = as.numeric(Y),
    weights = obsWeights,
    family = family$family,
    alpha = 1,
    nfolds = local_folds,
    nlambda = as.integer(nlambda),
    type.measure = "deviance",
    standardize = TRUE,
    intercept = TRUE,
    parallel = FALSE
  ))
  prediction <- as.numeric(stats::predict(
    fit, newx = new_design, type = "response", s = "lambda.1se"
  ))
  output <- list(
    pred = prediction,
    fit = list(object = fit, use_one_se = TRUE)
  )
  class(output$fit) <- "SL.binary_glmnet_interactions"
  output
}

predict.SL.binary_glmnet_interactions <- function(object, newdata, ...) {
  design <- binary_pairwise_interaction_design(newdata)
  as.numeric(stats::predict(
    object$object, newx = design, type = "response", s = "lambda.1se"
  ))
}

SL.binary_mgcv_additive <- function(
    Y, X, newX, family, obsWeights, basis_dimension = NULL, ...
) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("SL.binary_mgcv_additive requires package 'mgcv'.")
  }
  training <- as.data.frame(X)
  training$.outcome <- as.numeric(Y)
  predictors <- setdiff(names(training), ".outcome")
  local_basis <- binary_resolve_gam_basis_dimension(
    training, predictors, basis_dimension = basis_dimension,
    minimum = 5L, maximum = 20L, growth_exponent = 1 / 3
  )
  formula <- binary_make_gam_formula(
    ".outcome", training, predictors, local_basis
  )
  environment(formula) <- environment()
  fit <- suppressWarnings(mgcv::gam(
    formula,
    data = training,
    family = family,
    weights = obsWeights,
    method = "REML",
    select = TRUE
  ))
  prediction <- as.numeric(stats::predict(
    fit, newdata = as.data.frame(newX), type = "response"
  ))
  output <- list(pred = prediction, fit = list(object = fit))
  class(output$fit) <- "SL.binary_mgcv_additive"
  output
}

predict.SL.binary_mgcv_additive <- function(object, newdata, ...) {
  as.numeric(stats::predict(
    object$object, newdata = as.data.frame(newdata), type = "response"
  ))
}

SL.binary_mgcv_interactions <- function(
    Y, X, newX, family, obsWeights,
    basis_dimension = NULL, interaction_basis_dimension = 3L, ...
) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("SL.binary_mgcv_interactions requires package 'mgcv'.")
  }
  training <- as.data.frame(X)
  training$.outcome <- as.numeric(Y)
  predictors <- setdiff(names(training), ".outcome")
  local_basis <- binary_resolve_gam_basis_dimension(
    training, predictors, basis_dimension = basis_dimension,
    minimum = 5L, maximum = 20L, growth_exponent = 1 / 3
  )
  formula <- binary_make_interaction_gam_formula(
    ".outcome", training, predictors, local_basis,
    interaction_basis_dimension = interaction_basis_dimension
  )
  environment(formula) <- environment()
  fit <- suppressWarnings(mgcv::gam(
    formula,
    data = training,
    family = family,
    weights = obsWeights,
    method = "REML",
    select = TRUE
  ))
  prediction <- as.numeric(stats::predict(
    fit, newdata = as.data.frame(newX), type = "response"
  ))
  output <- list(pred = prediction, fit = list(object = fit))
  class(output$fit) <- "SL.binary_mgcv_interactions"
  output
}

predict.SL.binary_mgcv_interactions <- function(object, newdata, ...) {
  as.numeric(stats::predict(
    object$object, newdata = as.data.frame(newdata), type = "response"
  ))
}

SL.binary_ranger <- function(
    Y, X, newX, family, obsWeights, id,
    num_trees = 500L, min_node_size = 20L, ...
) {
  SuperLearner::SL.ranger(
    Y = Y,
    X = X,
    newX = newX,
    family = family,
    obsWeights = obsWeights,
    id = id,
    num.trees = as.integer(num_trees),
    min.node.size = as.integer(min_node_size),
    num.threads = 1L,
    verbose = FALSE
  )
}

SL.binary_xgboost <- function(
    Y, X, newX, family, obsWeights, id,
    ntrees = 250L, max_depth = 3L, shrinkage = 0.05,
    minobspernode = 10L, ...
) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("SL.binary_xgboost requires package 'xgboost'.")
  }
  if (!identical(family$family, "binomial")) {
    stop("SL.binary_xgboost currently supports binomial outcomes only.")
  }
  training_design <- as.matrix(X)
  validation_design <- as.matrix(newX)
  storage.mode(training_design) <- "double"
  storage.mode(validation_design) <- "double"
  training_matrix <- xgboost::xgb.DMatrix(
    data = training_design,
    label = as.numeric(Y),
    weight = as.numeric(obsWeights)
  )
  fit <- xgboost::xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "logloss",
      max_depth = as.integer(max_depth),
      eta = as.numeric(shrinkage),
      min_child_weight = as.integer(minobspernode),
      tree_method = "hist",
      nthread = 1L
    ),
    data = training_matrix,
    nrounds = as.integer(ntrees),
    verbose = 0
  )
  prediction <- as.numeric(stats::predict(
    fit, newdata = validation_design
  ))
  output <- list(
    pred = prediction,
    fit = list(object = fit)
  )
  class(output$fit) <- "SL.binary_xgboost"
  output
}

predict.SL.binary_xgboost <- function(object, newdata, ...) {
  design <- as.matrix(newdata)
  storage.mode(design) <- "double"
  as.numeric(stats::predict(object$object, newdata = design))
}

SL.binary_hal <- function(
    Y, X, newX, family, obsWeights, id,
    max_degree = 2L, num_knots = NULL, ...
) {
  if (!requireNamespace("hal9001", quietly = TRUE)) {
    stop("SL.binary_hal requires package 'hal9001'.")
  }
  if (is.null(num_knots)) {
    # The basis grows with the training sample instead of fixing the
    # approximation class across all sample sizes.
    num_knots <- max(5L, ceiling(nrow(X)^(1 / 3)))
  }
  hal9001::SL.hal9001(
    Y = Y,
    X = X,
    newX = newX,
    family = family,
    obsWeights = obsWeights,
    id = id,
    max_degree = as.integer(max_degree),
    smoothness_orders = 1,
    num_knots = as.integer(num_knots)
  )
}

predict.SL.binary_hal <- function(object, newdata, ...) {
  as.numeric(stats::predict(object$object, new_data = as.matrix(newdata)))
}

binary_super_learner_method <- function() {
  if (!requireNamespace("nloptr", quietly = TRUE)) {
    stop("The convex log-loss Super Learner requires package 'nloptr'.")
  }
  compute_coef <- function(
      Z, Y, libraryNames, verbose, obsWeights, control,
      errorsInLibrary = NULL, ...
  ) {
    if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1L)
    Z <- as.matrix(Z)
    probability <- Z
    probability[] <- pmin(1 - 1e-8, pmax(1e-8, as.numeric(Z)))
    weights <- as.numeric(obsWeights)
    weights <- weights / sum(weights)
    individual_risk <- vapply(seq_len(ncol(probability)), function(index) {
      binary_cross_entropy_loss(Y, probability[, index])
    }, numeric(1))
    names(individual_risk) <- libraryNames
    unavailable <- apply(!is.finite(probability), 2L, any)
    if (!is.null(errorsInLibrary)) unavailable <- unavailable | errorsInLibrary
    available <- which(!unavailable)
    if (!length(available)) stop("Every Super Learner candidate failed.")
    local_probability <- probability[, available, drop = FALSE]
    objective <- function(beta) {
      combined <- pmin(
        1 - 1e-8,
        pmax(1e-8, as.numeric(local_probability %*% beta))
      )
      loss <- -Y * log(combined) - (1 - Y) * log1p(-combined)
      gradient_factor <- weights * (
        (combined - Y) / (combined * (1 - combined))
      )
      list(
        objective = sum(weights * loss),
        gradient = as.numeric(crossprod(local_probability, gradient_factor))
      )
    }
    start <- rep(1 / length(available), length(available))
    optimizer <- nloptr::nloptr(
      x0 = start,
      eval_f = objective,
      lb = rep(0, length(available)),
      ub = rep(1, length(available)),
      eval_g_eq = function(beta) sum(beta) - 1,
      eval_jac_g_eq = function(beta) rep(1, length(beta)),
      opts = list(
        algorithm = "NLOPT_LD_SLSQP",
        xtol_rel = 1e-9,
        maxeval = 1000L,
        print_level = 0L
      )
    )
    coefficient <- rep(0, ncol(probability))
    coefficient[available] <- optimizer$solution
    coefficient[coefficient < 1e-8] <- 0
    coefficient <- coefficient / sum(coefficient)
    names(coefficient) <- libraryNames
    if (verbose) {
      message(
        "Convex log-loss optimizer status: ", optimizer$status,
        " (", optimizer$message, ")"
      )
    }
    list(
      cvRisk = individual_risk,
      coef = coefficient,
      optimizer = optimizer
    )
  }
  compute_prediction <- function(predY, coef, control, ...) {
    if (is.null(dim(predY))) predY <- matrix(predY, ncol = 1L)
    prediction <- as.numeric(as.matrix(predY) %*% as.numeric(coef))
    pmin(1 - 1e-8, pmax(1e-8, prediction))
  }
  list(
    require = "nloptr",
    computeCoef = compute_coef,
    computePred = compute_prediction
  )
}

binary_super_learner_candidate_map <- function() {
  c(
    glm_main = "SL.glm",
    glm_spline = "SL.binary_glm_spline",
    glm_interactions = "SL.glm.interaction",
    glmnet = "SL.glmnet",
    glmnet_interactions = "SL.binary_glmnet_interactions",
    gam = "SL.binary_mgcv_additive",
    gam_interactions = "SL.binary_mgcv_interactions",
    ranger = "SL.binary_ranger",
    xgboost = "SL.binary_xgboost",
    hal = "SL.binary_hal"
  )
}

binary_default_super_learner_candidates <- function() {
  # The production library avoids redundant fits while spanning parametric
  # main effects, explicit pairwise interactions, recursive interactions, and
  # boosted trees. The full candidate map, including HAL, remains available
  # for sensitivity analyses; the production library avoids HAL's growing
  # degree-two basis because every outer-fold nuisance fit is itself
  # cross-validated.
  c(
    "glm_main", "glm_interactions", "glmnet_interactions",
    "ranger", "xgboost"
  )
}

binary_fit_regression_super_learner <- function(
    data,
    response,
    x_cols,
    candidates = binary_default_super_learner_candidates(),
    cv_folds = 3L,
    ranger_num_trees = 500L,
    hal_max_degree = 2L,
    seed = 1L
) {
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("The binary primitive Super Learner requires 'SuperLearner'.")
  }
  candidate_map <- binary_super_learner_candidate_map()
  candidates <- unique(as.character(candidates))
  invalid <- setdiff(candidates, names(candidate_map))
  if (!length(candidates) || length(invalid)) {
    stop(
      "Unknown binary Super Learner candidate(s): ",
      paste(invalid, collapse = ", "), "."
    )
  }
  package_by_candidate <- c(
    glmnet = "glmnet", glmnet_interactions = "glmnet",
    gam = "mgcv", gam_interactions = "mgcv",
    ranger = "ranger", xgboost = "xgboost", hal = "hal9001"
  )
  missing_packages <- unique(vapply(candidates, function(candidate) {
    package <- unname(package_by_candidate[candidate])
    if (is.na(package) || requireNamespace(package, quietly = TRUE)) {
      ""
    } else {
      package
    }
  }, character(1)))
  missing_packages <- missing_packages[nzchar(missing_packages)]
  if (length(missing_packages)) {
    stop(
      "The requested binary Super Learner library is incomplete; install: ",
      paste(missing_packages, collapse = ", "), "."
    )
  }
  outcome <- as.numeric(data[[response]])
  binary_check_indicator(outcome, response)
  if (length(unique(outcome)) == 1L) {
    probability <- (sum(outcome) + 0.5) / (length(outcome) + 1)
    return(list(
      predict = function(newdata) rep(probability, nrow(newdata)),
      converged = TRUE,
      rank = 1L,
      constant = TRUE,
      learner = "super_learner",
      selected_learner = "constant",
      cv_mean_risk = NA_real_,
      cv_standard_error = NA_real_,
      cv_table = data.frame(),
      ensemble_weights = setNames(1, "constant")
    ))
  }
  cv_folds <- as.integer(cv_folds)
  minimum_class_size <- min(tabulate(outcome + 1L, nbins = 2L))
  cv_folds <- min(cv_folds, minimum_class_size)
  if (cv_folds < 2L) stop("Super Learner requires both outcome classes.")
  design <- binary_super_learner_numeric_design(data, x_cols)

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))

  # Wrapper-specific tuning constants are local to this fit and are therefore
  # preserved by the environment passed to SuperLearner.
  All <- SuperLearner::All
  SL.glm <- SuperLearner::SL.glm
  SL.glm.interaction <- SuperLearner::SL.glm.interaction
  SL.glmnet <- SuperLearner::SL.glmnet
  base_binary_ranger <- SL.binary_ranger
  base_binary_xgboost <- SL.binary_xgboost
  base_binary_hal <- SL.binary_hal
  SL.binary_ranger <- function(
      Y, X, newX, family, obsWeights, id, ...
  ) {
    base_binary_ranger(
      Y = Y, X = X, newX = newX, family = family,
      obsWeights = obsWeights, id = id,
      num_trees = ranger_num_trees
    )
  }
  SL.binary_hal <- function(
      Y, X, newX, family, obsWeights, id, ...
  ) {
    base_binary_hal(
      Y = Y, X = X, newX = newX, family = family,
      obsWeights = obsWeights, id = id,
      max_degree = hal_max_degree
    )
  }
  SL.binary_xgboost <- function(
      Y, X, newX, family, obsWeights, id, ...
  ) {
    base_binary_xgboost(
      Y = Y, X = X, newX = newX, family = family,
      obsWeights = obsWeights, id = id
    )
  }
  fit_environment <- environment()
  fitted <- suppressWarnings(SuperLearner::SuperLearner(
    Y = outcome,
    X = as.data.frame(design$training),
    family = stats::binomial(),
    SL.library = unname(candidate_map[candidates]),
    method = binary_super_learner_method(),
    verbose = FALSE,
    control = list(saveFitLibrary = TRUE, trimLogit = 1e-5),
    cvControl = list(
      V = cv_folds, stratifyCV = TRUE, shuffle = TRUE
    ),
    env = fit_environment
  ))
  coefficients <- as.numeric(fitted$coef)
  names(coefficients) <- candidates
  optimizer_status <- fitted$method$optimizer$status
  converged <- is.null(optimizer_status) || optimizer_status %in% 1:4
  cv_table <- data.frame(
    candidate = candidates,
    wrapper = unname(candidate_map[candidates]),
    mean_risk = as.numeric(fitted$cvRisk),
    coefficient = coefficients,
    selected = coefficients > 1e-8,
    cv_available = TRUE,
    stringsAsFactors = FALSE
  )
  predict_probability <- function(newdata) {
    new_design <- as.data.frame(design$transform(newdata))
    prediction <- stats::predict(fitted, newdata = new_design, onlySL = TRUE)
    as.numeric(prediction$pred)
  }
  list(
    predict = predict_probability,
    converged = converged,
    rank = as.integer(sum(coefficients > 1e-8)),
    constant = FALSE,
    learner = "super_learner",
    selected_learner = "convex_stack",
    cv_mean_risk = binary_cross_entropy_loss(
      outcome, as.numeric(fitted$Z %*% coefficients)
    ),
    cv_standard_error = NA_real_,
    cv_table = cv_table,
    fold_risk = data.frame(),
    ensemble_weights = coefficients,
    super_learner_fit = fitted,
    library_configuration = list(
      candidates_requested = candidates,
      candidates_available = candidates,
      cv_folds = cv_folds,
      meta_learner = "simplex_probability_mixture_bernoulli_log_loss",
      ranger_num_trees = as.integer(ranger_num_trees),
      xgboost_num_trees = 250L,
      xgboost_max_depth = 3L,
      xgboost_shrinkage = 0.05,
      xgboost_minimum_child_weight = 10L,
      xgboost_threads = 1L,
      hal_max_degree = as.integer(hal_max_degree),
      hal_num_knots_rule = "max(5, ceiling(n_training^(1/3)))"
    )
  )
}

binary_resolve_gam_basis_dimension <- function(
    data,
    x_cols,
    basis_dimension = NULL,
    minimum = 5L,
    maximum = 20L,
    growth_exponent = 1 / 3
) {
  if (!is.null(basis_dimension)) {
    basis_dimension <- as.integer(basis_dimension)
    if (length(basis_dimension) != 1L || is.na(basis_dimension) ||
        basis_dimension < 3L) {
      stop("gam_basis_dimension must be NULL or one integer of at least 3.")
    }
    return(basis_dimension)
  }
  minimum <- as.integer(minimum)
  maximum <- as.integer(maximum)
  if (length(minimum) != 1L || is.na(minimum) || minimum < 3L ||
      length(maximum) != 1L || is.na(maximum) || maximum < minimum ||
      length(growth_exponent) != 1L || !is.finite(growth_exponent) ||
      growth_exponent <= 0 || growth_exponent > 1) {
    stop("Invalid adaptive GAM basis controls.")
  }

  # The basis dimension is only an upper limit on wiggliness; REML and
  # select=TRUE determine the effective degrees of freedom.  Letting this
  # upper limit grow slowly with the observations available per covariate
  # avoids fixing approximation bias as the sample size increases.
  observations_per_covariate <- nrow(data) / max(1L, length(x_cols))
  proposed <- ceiling(observations_per_covariate^growth_exponent)
  as.integer(min(maximum, max(minimum, proposed)))
}

binary_cross_entropy_loss <- function(outcome, probability) {
  outcome <- as.numeric(outcome)
  probability <- pmin(1 - 1e-8, pmax(1e-8, as.numeric(probability)))
  if (length(outcome) != length(probability) ||
      any(!is.finite(outcome)) || any(!is.finite(probability))) {
    stop("Invalid inputs to binary_cross_entropy_loss().")
  }
  mean(-outcome * log(probability) -
    (1 - outcome) * log1p(-probability))
}

binary_make_response_folds <- function(outcome, K = 3L, seed = 1L) {
  outcome <- as.numeric(outcome)
  binary_check_indicator(outcome, "outcome")
  K <- as.integer(K)
  class_sizes <- tabulate(outcome + 1L, nbins = 2L)
  if (length(unique(outcome)) < 2L || min(class_sizes) < 2L) {
    return(NULL)
  }
  K <- min(K, min(class_sizes))
  if (K < 2L) return(NULL)

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  fold_id <- integer(length(outcome))
  for (outcome_value in c(0, 1)) {
    index <- which(outcome == outcome_value)
    fold_id[index] <- sample(rep(seq_len(K), length.out = length(index)))
  }
  fold_id
}

binary_available_library_candidates <- function(candidates) {
  candidates <- unique(as.character(candidates))
  allowed <- c("glm_spline", "gam", "glmnet", "ranger")
  if (length(candidates) < 1L || any(!candidates %in% allowed)) {
    stop(
      "Binary nuisance-library candidates must be drawn from: ",
      paste(allowed, collapse = ", "), "."
    )
  }
  package_by_candidate <- c(gam = "mgcv", glmnet = "glmnet", ranger = "ranger")
  available <- vapply(candidates, function(candidate) {
    if (candidate == "glm_spline") return(TRUE)
    requireNamespace(package_by_candidate[[candidate]], quietly = TRUE)
  }, logical(1))
  candidates[available]
}

binary_fit_regression_library_cv <- function(
    data,
    response,
    x_cols,
    spline_df,
    gam_basis_dimension = NULL,
    gam_min_basis = 5L,
    gam_max_basis = 20L,
    gam_growth_exponent = 1 / 3,
    seed = 1L,
    candidates = c("glm_spline", "gam", "glmnet", "ranger"),
    cv_folds = 3L,
    selection_rule = c("one_se", "minimum"),
    ranger_num_trees = 500L
) {
  selection_rule <- match.arg(selection_rule)
  candidates_requested <- unique(as.character(candidates))
  candidates <- binary_available_library_candidates(candidates_requested)
  if (length(candidates) < 1L) {
    stop("No requested binary nuisance-library candidate is available.")
  }
  outcome <- as.numeric(data[[response]])
  binary_check_indicator(outcome, response)

  fit_candidate <- function(training, candidate, local_seed) {
    binary_fit_regression_with_learner(
      data = training,
      response = response,
      x_cols = x_cols,
      spline_df = spline_df,
      learner = candidate,
      gam_basis_dimension = gam_basis_dimension,
      gam_min_basis = gam_min_basis,
      gam_max_basis = gam_max_basis,
      gam_growth_exponent = gam_growth_exponent,
      ranger_num_trees = ranger_num_trees,
      seed = local_seed
    )
  }

  fold_id <- binary_make_response_folds(outcome, K = cv_folds, seed = seed)
  if (is.null(fold_id)) {
    selected <- if ("glm_spline" %in% candidates) "glm_spline" else candidates[1L]
    final_fit <- fit_candidate(data, selected, seed + 900000L)
    final_fit$learner <- "library_cv"
    final_fit$selected_learner <- selected
    final_fit$cv_mean_risk <- NA_real_
    final_fit$cv_standard_error <- NA_real_
    final_fit$cv_table <- data.frame(
      candidate = selected,
      mean_risk = NA_real_,
      standard_error = NA_real_,
      selected = TRUE,
      cv_available = FALSE,
      stringsAsFactors = FALSE
    )
    final_fit$library_configuration <- list(
      candidates_requested = candidates_requested,
      candidates_available = candidates,
      cv_folds = 0L,
      selection_rule = selection_rule,
      selected_learner = selected
    )
    return(final_fit)
  }

  folds <- sort(unique(fold_id))
  risk_rows <- list()
  counter <- 0L
  for (candidate_index in seq_along(candidates)) {
    candidate <- candidates[candidate_index]
    for (fold in folds) {
      validation <- fold_id == fold
      fitted <- tryCatch(
        fit_candidate(
          data[!validation, , drop = FALSE],
          candidate,
          seed + 10000L * candidate_index + fold
        ),
        error = identity
      )
      counter <- counter + 1L
      if (inherits(fitted, "error")) {
        risk_rows[[counter]] <- data.frame(
          candidate = candidate,
          fold = fold,
          risk = Inf,
          error = conditionMessage(fitted),
          stringsAsFactors = FALSE
        )
      } else {
        prediction <- fitted$predict(data[validation, , drop = FALSE])
        risk_rows[[counter]] <- data.frame(
          candidate = candidate,
          fold = fold,
          risk = binary_cross_entropy_loss(outcome[validation], prediction),
          error = "",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  fold_risk <- do.call(rbind, risk_rows)
  cv_table <- do.call(rbind, lapply(candidates, function(candidate) {
    local <- fold_risk[fold_risk$candidate == candidate, , drop = FALSE]
    finite <- is.finite(local$risk)
    data.frame(
      candidate = candidate,
      mean_risk = if (all(finite)) mean(local$risk) else Inf,
      standard_error = if (all(finite) && nrow(local) > 1L) {
        stats::sd(local$risk) / sqrt(nrow(local))
      } else if (all(finite)) {
        0
      } else {
        NA_real_
      },
      failed_folds = sum(!finite),
      error = paste(unique(local$error[nzchar(local$error)]), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }))
  if (!any(is.finite(cv_table$mean_risk))) {
    stop("Every binary nuisance-library candidate failed cross-validation.")
  }
  minimum <- which.min(cv_table$mean_risk)
  choice <- minimum
  if (selection_rule == "one_se") {
    cutoff <- cv_table$mean_risk[minimum] +
      ifelse(is.finite(cv_table$standard_error[minimum]),
        cv_table$standard_error[minimum], 0)
    eligible <- which(cv_table$mean_risk <= cutoff)
    # Candidate order is the prespecified stability/complexity ordering.
    choice <- eligible[1L]
  }
  cv_table$selected <- seq_len(nrow(cv_table)) == choice
  cv_table$cv_available <- TRUE
  selected <- cv_table$candidate[choice]
  final_fit <- fit_candidate(data, selected, seed + 900000L + choice)
  final_fit$learner <- "library_cv"
  final_fit$selected_learner <- selected
  final_fit$cv_mean_risk <- cv_table$mean_risk[choice]
  final_fit$cv_standard_error <- cv_table$standard_error[choice]
  final_fit$cv_table <- cv_table
  final_fit$fold_risk <- fold_risk
  final_fit$library_configuration <- list(
    candidates_requested = candidates_requested,
    candidates_available = candidates,
    cv_folds = length(folds),
    selection_rule = selection_rule,
    selected_learner = selected
  )
  final_fit
}

binary_fit_regression_with_learner <- function(
    data,
    response,
    x_cols,
    spline_df,
    learner = c(
      "glm_spline", "gam", "ranger", "glmnet", "library_cv",
      "super_learner"
    ),
    gam_basis_dimension = NULL,
    gam_min_basis = 5L,
    gam_max_basis = 20L,
    gam_growth_exponent = 1 / 3,
    library_candidates = c("glm_spline", "gam", "glmnet", "ranger"),
    library_cv_folds = 3L,
    library_selection_rule = c("one_se", "minimum"),
    ranger_num_trees = 500L,
    super_learner_candidates = binary_default_super_learner_candidates(),
    super_learner_cv_folds = 3L,
    super_learner_ranger_num_trees = 500L,
    super_learner_hal_max_degree = 2L,
    seed = 1L
) {
  learner <- match.arg(learner)
  library_selection_rule <- match.arg(library_selection_rule)
  outcome <- as.numeric(data[[response]])
  binary_check_indicator(outcome, response)
  if (length(unique(outcome)) == 1L) {
    probability <- (sum(outcome) + 0.5) / (length(outcome) + 1)
    return(list(
      predict = function(newdata) rep(probability, nrow(newdata)),
      converged = TRUE,
      rank = 1L,
      constant = TRUE,
      learner = learner,
      selected_learner = if (learner == "library_cv") "constant" else learner,
      cv_mean_risk = NA_real_,
      cv_standard_error = NA_real_,
      cv_table = data.frame()
    ))
  }
  if (learner == "library_cv") {
    return(binary_fit_regression_library_cv(
      data = data,
      response = response,
      x_cols = x_cols,
      spline_df = spline_df,
      gam_basis_dimension = gam_basis_dimension,
      gam_min_basis = gam_min_basis,
      gam_max_basis = gam_max_basis,
      gam_growth_exponent = gam_growth_exponent,
      seed = seed,
      candidates = library_candidates,
      cv_folds = library_cv_folds,
      selection_rule = library_selection_rule,
      ranger_num_trees = ranger_num_trees
    ))
  }
  if (learner == "super_learner") {
    return(binary_fit_regression_super_learner(
      data = data,
      response = response,
      x_cols = x_cols,
      candidates = super_learner_candidates,
      cv_folds = super_learner_cv_folds,
      ranger_num_trees = super_learner_ranger_num_trees,
      hal_max_degree = super_learner_hal_max_degree,
      seed = seed
    ))
  }
  if (learner == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("The binary GAM nuisance learner requires package 'mgcv'.")
    }
    selected_basis_dimension <- binary_resolve_gam_basis_dimension(
      data = data,
      x_cols = x_cols,
      basis_dimension = gam_basis_dimension,
      minimum = gam_min_basis,
      maximum = gam_max_basis,
      growth_exponent = gam_growth_exponent
    )
    fit <- suppressWarnings(mgcv::gam(
      binary_make_gam_formula(
        response, data, x_cols, selected_basis_dimension
      ),
      data = data,
      family = stats::binomial(),
      method = "REML",
      select = TRUE
    ))
    return(list(
      predict = function(newdata) {
        as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
      },
      converged = isTRUE(fit$converged),
      rank = fit$rank,
      constant = FALSE,
      learner = learner,
      basis_dimension = selected_basis_dimension,
      effective_df = sum(fit$edf)
    ))
  }
  if (learner == "ranger") {
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("The binary ranger nuisance learner requires package 'ranger'.")
    }
    ranger_data <- data[, c(response, x_cols), drop = FALSE]
    ranger_data[[response]] <- factor(outcome, levels = c(0, 1))
    fit <- ranger::ranger(
      formula = stats::reformulate(x_cols, response = response),
      data = ranger_data,
      probability = TRUE,
      num.trees = as.integer(ranger_num_trees),
      mtry = max(1L, floor(sqrt(length(x_cols)))),
      min.node.size = 20L,
      seed = as.integer(seed),
      respect.unordered.factors = "order"
    )
    return(list(
      predict = function(newdata) {
        prediction <- predict(fit, data = newdata)$predictions
        if (is.null(dim(prediction))) return(as.numeric(prediction))
        if ("1" %in% colnames(prediction)) {
          return(as.numeric(prediction[, "1"]))
        }
        as.numeric(prediction[, ncol(prediction)])
      },
      converged = TRUE,
      rank = as.integer(length(x_cols)),
      constant = FALSE,
      learner = learner
    ))
  }
  if (learner == "glmnet") {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("The binary glmnet nuisance learner requires package 'glmnet'.")
    }
    design <- as.matrix(data[, x_cols, drop = FALSE])
    storage.mode(design) <- "double"
    if (any(!is.finite(design))) {
      stop("The glmnet nuisance design contains non-finite values.")
    }
    single_predictor_continuous <- ncol(design) == 1L &&
      length(unique(design[, 1L])) > 2L
    if (ncol(design) == 1L) {
      auxiliary <- if (single_predictor_continuous) {
        design[, 1L]^2
      } else {
        rep(0, nrow(design))
      }
      design <- cbind(design, `.glmnet_auxiliary` = auxiliary)
    }
    minimum_class_size <- min(tabulate(outcome + 1L, nbins = 2L))
    cv_folds <- min(5L, minimum_class_size)
    if (cv_folds < 2L) {
      stop("The glmnet nuisance fit requires both outcome classes.")
    }

    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
    fold_id <- integer(length(outcome))
    for (outcome_value in c(0, 1)) {
      index <- which(outcome == outcome_value)
      fold_id[index] <- sample(rep(seq_len(cv_folds), length.out = length(index)))
    }
    fit <- suppressWarnings(glmnet::cv.glmnet(
      x = design,
      y = outcome,
      family = "binomial",
      alpha = 1,
      foldid = fold_id,
      type.measure = "deviance",
      standardize = TRUE,
      intercept = TRUE,
      parallel = FALSE
    ))
    coefficients <- stats::coef(fit, s = "lambda.1se")
    active_count <- sum(as.numeric(coefficients) != 0) - 1L
    return(list(
      predict = function(newdata) {
        new_design <- as.matrix(newdata[, x_cols, drop = FALSE])
        storage.mode(new_design) <- "double"
        if (ncol(new_design) == 1L) {
          auxiliary <- if (single_predictor_continuous) {
            new_design[, 1L]^2
          } else {
            rep(0, nrow(new_design))
          }
          new_design <- cbind(
            new_design, `.glmnet_auxiliary` = auxiliary
          )
        }
        as.numeric(stats::predict(
          fit, newx = new_design, s = "lambda.1se", type = "response"
        ))
      },
      converged = TRUE,
      rank = as.integer(max(0L, active_count)),
      constant = FALSE,
      learner = learner,
      selected_lambda = unname(fit$lambda.1se),
      minimum_lambda = unname(fit$lambda.min)
    ))
  }
  fit <- suppressWarnings(stats::glm(
    binary_make_spline_formula(response, data, x_cols, spline_df),
    data = data,
    family = stats::binomial(),
    control = stats::glm.control(maxit = 100L)
  ))
  list(
    predict = function(newdata) {
      as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
    },
    converged = isTRUE(fit$converged),
    rank = fit$rank,
    constant = FALSE,
    learner = learner
  )
}

fit_binary_cross_fitted_nuisances <- function(
    data,
    x_cols = "X",
    K = 5L,
    fold_id = NULL,
    spline_df = 5L,
    learner = c(
      "library_cv", "super_learner", "gam", "glm_spline", "ranger",
      "glmnet"
    ),
    overlap_epsilon = 0.05,
    outcome_probability_epsilon = 1e-4,
    probability_epsilon = NULL,
    gam_basis_dimension = NULL,
    gam_min_basis = 5L,
    gam_max_basis = 20L,
    gam_growth_exponent = 1 / 3,
    library_candidates = c("glm_spline", "gam", "glmnet", "ranger"),
    library_cv_folds = 3L,
    library_selection_rule = c("one_se", "minimum"),
    library_ranger_num_trees = 500L,
    super_learner_candidates = binary_default_super_learner_candidates(),
    super_learner_cv_folds = 3L,
    super_learner_ranger_num_trees = 500L,
    super_learner_hal_max_degree = 2L,
    seed = 1L
) {
  learner <- match.arg(learner)
  library_selection_rule <- match.arg(library_selection_rule)
  if (!is.null(probability_epsilon)) {
    # Backward-compatible interpretation: the historical argument controls
    # only treatment/instrument overlap.  It must not impose 0.05 truncation
    # on the observed outcome regressions a_{tz}.
    overlap_epsilon <- probability_epsilon
  }
  validate_epsilon <- function(value, name) {
    if (length(value) != 1L || !is.finite(value) ||
        value < 0 || value >= 0.5) {
      stop(name, " must be one number in [0, 0.5).")
    }
  }
  validate_epsilon(overlap_epsilon, "overlap_epsilon")
  validate_epsilon(
    outcome_probability_epsilon, "outcome_probability_epsilon"
  )
  required <- c("Y", "T_obs", "Z", x_cols)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }
  binary_check_indicator(data$Y, "Y")
  binary_check_indicator(data$T_obs, "T_obs")
  binary_check_indicator(data$Z, "Z")
  if (is.null(fold_id)) {
    fold_id <- binary_make_tz_folds(data, K = K, seed = seed)
  }
  if (length(fold_id) != nrow(data)) stop("fold_id has the wrong length.")
  folds <- sort(unique(fold_id))
  if (length(folds) < 2L) stop("At least two folds are required.")

  prediction_names <- c(
    "e", "pi0", "pi1",
    "a_t0_z0", "a_t0_z1", "a_t1_z0", "a_t1_z1"
  )
  predictions <- as.data.frame(
    setNames(replicate(length(prediction_names),
      rep(NA_real_, nrow(data)), simplify = FALSE), prediction_names)
  )
  diagnostics <- list()
  cv_tables <- list()

  for (fold in folds) {
    training <- data[fold_id != fold, , drop = FALSE]
    validation <- data[fold_id == fold, , drop = FALSE]
    validation_index <- which(fold_id == fold)

    fit_z <- binary_fit_regression_with_learner(
      training, "Z", x_cols, spline_df, learner,
      gam_basis_dimension = gam_basis_dimension,
      gam_min_basis = gam_min_basis,
      gam_max_basis = gam_max_basis,
      gam_growth_exponent = gam_growth_exponent,
      library_candidates = library_candidates,
      library_cv_folds = library_cv_folds,
      library_selection_rule = library_selection_rule,
      ranger_num_trees = library_ranger_num_trees,
      super_learner_candidates = super_learner_candidates,
      super_learner_cv_folds = super_learner_cv_folds,
      super_learner_ranger_num_trees = super_learner_ranger_num_trees,
      super_learner_hal_max_degree = super_learner_hal_max_degree,
      seed = seed + 1000L * fold + 1L
    )
    fit_pi0 <- binary_fit_regression_with_learner(
      training[training$Z == 0, , drop = FALSE],
      "T_obs", x_cols, spline_df, learner,
      gam_basis_dimension = gam_basis_dimension,
      gam_min_basis = gam_min_basis,
      gam_max_basis = gam_max_basis,
      gam_growth_exponent = gam_growth_exponent,
      library_candidates = library_candidates,
      library_cv_folds = library_cv_folds,
      library_selection_rule = library_selection_rule,
      ranger_num_trees = library_ranger_num_trees,
      super_learner_candidates = super_learner_candidates,
      super_learner_cv_folds = super_learner_cv_folds,
      super_learner_ranger_num_trees = super_learner_ranger_num_trees,
      super_learner_hal_max_degree = super_learner_hal_max_degree,
      seed = seed + 1000L * fold + 2L
    )
    fit_pi1 <- binary_fit_regression_with_learner(
      training[training$Z == 1, , drop = FALSE],
      "T_obs", x_cols, spline_df, learner,
      gam_basis_dimension = gam_basis_dimension,
      gam_min_basis = gam_min_basis,
      gam_max_basis = gam_max_basis,
      gam_growth_exponent = gam_growth_exponent,
      library_candidates = library_candidates,
      library_cv_folds = library_cv_folds,
      library_selection_rule = library_selection_rule,
      ranger_num_trees = library_ranger_num_trees,
      super_learner_candidates = super_learner_candidates,
      super_learner_cv_folds = super_learner_cv_folds,
      super_learner_ranger_num_trees = super_learner_ranger_num_trees,
      super_learner_hal_max_degree = super_learner_hal_max_degree,
      seed = seed + 1000L * fold + 3L
    )
    cell_fits <- list()
    for (t_value in c(0, 1)) {
      for (z_value in c(0, 1)) {
        key <- sprintf("a_t%d_z%d", t_value, z_value)
        cell <- training$T_obs == t_value & training$Z == z_value
        cell_fits[[key]] <- binary_fit_regression_with_learner(
          training[cell, , drop = FALSE], "Y", x_cols, spline_df, learner,
          gam_basis_dimension = gam_basis_dimension,
          gam_min_basis = gam_min_basis,
          gam_max_basis = gam_max_basis,
          gam_growth_exponent = gam_growth_exponent,
          library_candidates = library_candidates,
          library_cv_folds = library_cv_folds,
          library_selection_rule = library_selection_rule,
          ranger_num_trees = library_ranger_num_trees,
          super_learner_candidates = super_learner_candidates,
          super_learner_cv_folds = super_learner_cv_folds,
          super_learner_ranger_num_trees = super_learner_ranger_num_trees,
          super_learner_hal_max_degree = super_learner_hal_max_degree,
          seed = seed + 1000L * fold + 100L * t_value + 10L * z_value + 4L
        )
      }
    }

    raw_predictions <- list(
      e = 1 - fit_z$predict(validation),
      pi0 = fit_pi0$predict(validation),
      pi1 = fit_pi1$predict(validation)
    )
    for (key in names(cell_fits)) {
      raw_predictions[[key]] <- cell_fits[[key]]$predict(validation)
    }
    predictions$e[validation_index] <- binary_clip_probability(
      raw_predictions$e, overlap_epsilon
    )
    predictions$pi0[validation_index] <- binary_clip_probability(
      raw_predictions$pi0, overlap_epsilon
    )
    predictions$pi1[validation_index] <- binary_clip_probability(
      raw_predictions$pi1, overlap_epsilon
    )
    for (key in names(cell_fits)) {
      predictions[[key]][validation_index] <- binary_clip_probability(
        raw_predictions[[key]], outcome_probability_epsilon
      )
    }

    all_fits <- c(list(e = fit_z, pi0 = fit_pi0, pi1 = fit_pi1), cell_fits)
    component_epsilon <- c(
      e = overlap_epsilon,
      pi0 = overlap_epsilon,
      pi1 = overlap_epsilon,
      stats::setNames(
        rep(outcome_probability_epsilon, length(cell_fits)),
        names(cell_fits)
      )
    )
    diagnostics[[as.character(fold)]] <- data.frame(
      fold = fold,
      component = names(all_fits),
      learner = vapply(all_fits, `[[`, character(1), "learner"),
      selected_learner = vapply(all_fits, function(fit) {
        if (is.null(fit$selected_learner)) fit$learner else fit$selected_learner
      }, character(1)),
      converged = vapply(all_fits, `[[`, logical(1), "converged"),
      rank = vapply(all_fits, `[[`, integer(1), "rank"),
      constant_fit = vapply(all_fits, `[[`, logical(1), "constant"),
      basis_dimension = vapply(all_fits, function(fit) {
        if (is.null(fit$basis_dimension)) NA_integer_ else
          as.integer(fit$basis_dimension)
      }, integer(1)),
      effective_df = vapply(all_fits, function(fit) {
        if (is.null(fit$effective_df)) NA_real_ else
          as.numeric(fit$effective_df)
      }, numeric(1)),
      cv_mean_risk = vapply(all_fits, function(fit) {
        if (is.null(fit$cv_mean_risk)) NA_real_ else fit$cv_mean_risk
      }, numeric(1)),
      cv_standard_error = vapply(all_fits, function(fit) {
        if (is.null(fit$cv_standard_error)) NA_real_ else
          fit$cv_standard_error
      }, numeric(1)),
      truncation_epsilon = unname(component_epsilon[names(all_fits)]),
      lower_truncation_rate = vapply(names(all_fits), function(component) {
        epsilon <- component_epsilon[[component]]
        mean(raw_predictions[[component]] < epsilon)
      }, numeric(1)),
      upper_truncation_rate = vapply(names(all_fits), function(component) {
        epsilon <- component_epsilon[[component]]
        mean(raw_predictions[[component]] > 1 - epsilon)
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
    for (component in names(all_fits)) {
      table <- all_fits[[component]]$cv_table
      if (!is.null(table) && nrow(table) > 0L) {
        table$outer_fold <- fold
        table$component <- component
        cv_tables[[length(cv_tables) + 1L]] <- table
      }
    }
  }

  if (any(!is.finite(as.matrix(predictions)))) {
    stop("At least one cross-fitted nuisance prediction is non-finite.")
  }
  list(
    predictions = predictions,
    fold_id = fold_id,
    learner = learner,
    diagnostics = do.call(rbind, diagnostics),
    cv_table = if (length(cv_tables) > 0L) {
      do.call(rbind, cv_tables)
    } else {
      data.frame()
    },
    configuration = list(
      overlap_epsilon = overlap_epsilon,
      outcome_probability_epsilon = outcome_probability_epsilon,
      gam_basis_dimension = gam_basis_dimension,
      gam_min_basis = as.integer(gam_min_basis),
      gam_max_basis = as.integer(gam_max_basis),
      gam_growth_exponent = gam_growth_exponent,
      library_candidates = library_candidates,
      library_cv_folds = as.integer(library_cv_folds),
      library_selection_rule = library_selection_rule,
      library_ranger_num_trees = as.integer(library_ranger_num_trees),
      super_learner_candidates = super_learner_candidates,
      super_learner_cv_folds = as.integer(super_learner_cv_folds),
      super_learner_ranger_num_trees = as.integer(
        super_learner_ranger_num_trees
      ),
      super_learner_hal_max_degree = as.integer(
        super_learner_hal_max_degree
      )
    )
  )
}

binary_arm_nuisance_from_predictions <- function(predictions, t_value) {
  if (t_value == 1) {
    list(
      r0 = predictions$pi0,
      r1 = predictions$pi1,
      a0 = predictions$a_t1_z0,
      a1 = predictions$a_t1_z1
    )
  } else {
    list(
      r0 = 1 - predictions$pi0,
      r1 = 1 - predictions$pi1,
      a0 = predictions$a_t0_z0,
      a1 = predictions$a_t0_z1
    )
  }
}

binary_project_arm_outcome_probability <- function(
    arm_nuisance,
    Lambda0,
    Lambda1 = Lambda0,
    action = c("project", "none"),
    tolerance = 1e-12
) {
  action <- match.arg(action)
  arm_nuisance <- binary_recycle_arm_nuisance(arm_nuisance)
  interval <- binary_gamma_compatible_interval(
    arm_nuisance$r0, arm_nuisance$r1, Lambda0, Lambda1
  )

  # For a binary outcome, gamma(1) = a1/a0 and
  # gamma(0) = (1-a1)/(1-a0).  Intersecting both gamma restrictions gives
  # the exact feasible interval for a1 while holding (r0,r1,a0) fixed.
  lower <- pmax(
    interval$lower * arm_nuisance$a0,
    1 - interval$upper * (1 - arm_nuisance$a0)
  )
  upper <- pmin(
    interval$upper * arm_nuisance$a0,
    1 - interval$lower * (1 - arm_nuisance$a0)
  )
  if (any(lower > upper + tolerance)) {
    stop(
      "The binary compatible interval for a1 is empty; no midpoint repair ",
      "was applied."
    )
  }
  # Only absorb machine-precision reversals.  Substantive crossings stop.
  upper <- pmax(upper, lower)
  projected_a1 <- pmin(upper, pmax(lower, arm_nuisance$a1))
  magnitude <- abs(projected_a1 - arm_nuisance$a1)
  required <- magnitude > tolerance

  stabilized <- arm_nuisance
  if (action == "project") stabilized$a1 <- projected_a1
  list(
    nuisance = stabilized,
    original_a1 = arm_nuisance$a1,
    projected_a1 = projected_a1,
    lower = lower,
    upper = upper,
    required = required,
    magnitude = magnitude,
    action = action,
    tolerance = tolerance
  )
}

binary_arm_one_step_signal <- function(
    data, predictions, arm_nuisance, endpoint, t_value, direction
) {
  direction <- match.arg(direction, c("lower", "upper"))
  gradient <- if (direction == "lower") {
    endpoint$gradient_lower
  } else {
    endpoint$gradient_upper
  }
  value <- if (direction == "lower") endpoint$lower else endpoint$upper
  observed_arm <- as.numeric(data$T_obs == t_value)
  signal <- value
  for (z_value in c(0, 1)) {
    z_probability <- if (z_value == 0) predictions$e else 1 - predictions$e
    z_indicator <- as.numeric(data$Z == z_value)
    r_name <- paste0("r", z_value)
    a_name <- paste0("a", z_value)
    r_value <- arm_nuisance[[r_name]]
    a_value <- arm_nuisance[[a_name]]
    signal <- signal + gradient[, r_name] * z_indicator / z_probability *
      (observed_arm - r_value)
    signal <- signal + gradient[, a_name] * z_indicator * observed_arm /
      (z_probability * r_value) * (data$Y - a_value)
  }
  signal
}

binary_standard_msm_arm <- function(r, a, Lambda) {
  theta_lower <- pmax(a / Lambda, 1 - (1 - a) * Lambda)
  theta_upper <- pmin(a * Lambda, 1 - (1 - a) / Lambda)
  list(
    lower = r * a + (1 - r) * theta_lower,
    upper = r * a + (1 - r) * theta_upper
  )
}

binary_comparator_plugin_values <- function(predictions, Lambda) {
  e <- predictions$e
  output <- list()
  for (t_value in c(0, 1)) {
    arm <- binary_arm_nuisance_from_predictions(predictions, t_value)

    marginal_r <- e * arm$r0 + (1 - e) * arm$r1
    marginal_a <- (
      e * arm$r0 * arm$a0 + (1 - e) * arm$r1 * arm$a1
    ) / marginal_r
    tan_x <- binary_standard_msm_arm(marginal_r, marginal_a, Lambda)

    tan_z0 <- binary_standard_msm_arm(arm$r0, arm$a0, Lambda)
    tan_z1 <- binary_standard_msm_arm(arm$r1, arm$a1, Lambda)
    tan_xz <- list(
      lower = e * tan_z0$lower + (1 - e) * tan_z1$lower,
      upper = e * tan_z0$upper + (1 - e) * tan_z1$upper
    )

    iv_only <- list(
      lower = pmax(arm$r0 * arm$a0, arm$r1 * arm$a1),
      upper = pmin(
        arm$r0 * arm$a0 + 1 - arm$r0,
        arm$r1 * arm$a1 + 1 - arm$r1
      )
    )
    output[[as.character(t_value)]] <- list(
      tan_x = tan_x, tan_xz = tan_xz, iv_only = iv_only
    )
  }

  make_row <- function(method, component) {
    mu1_lower <- mean(output[["1"]][[component]]$lower)
    mu1_upper <- mean(output[["1"]][[component]]$upper)
    mu0_lower <- mean(output[["0"]][[component]]$lower)
    mu0_upper <- mean(output[["0"]][[component]]$upper)
    data.frame(
      method = method,
      mu1_lower = mu1_lower,
      mu1_upper = mu1_upper,
      mu0_lower = mu0_lower,
      mu0_upper = mu0_upper,
      ate_lower = mu1_lower - mu0_upper,
      ate_upper = mu1_upper - mu0_lower,
      stringsAsFactors = FALSE
    )
  }
  rbind(
    make_row("tan_x_plugin", "tan_x"),
    make_row("tan_xz_plugin", "tan_xz"),
    make_row("iv_only_plugin", "iv_only")
  )
}

evaluate_binary_bounds_from_nuisances <- function(
    data,
    nuisance_fit,
    Lambda0 = 2,
    Lambda1 = Lambda0,
    epsilon = 0.01,
    confidence_level = 0.95,
    compatibility_action = c("project", "none"),
    compatibility_tolerance = 1e-12
) {
  compatibility_action <- match.arg(compatibility_action)
  predictions <- nuisance_fit$predictions
  if (!is.finite(confidence_level) || confidence_level <= 0 ||
      confidence_level >= 1) {
    stop("confidence_level must lie in (0,1).")
  }
  hard <- smooth <- signals <- unstabilized_hard <-
    unstabilized_smooth <- stabilization <- list()
  iv_hard <- iv_smooth <- iv_signals <- list()
  for (t_value in c(0, 1)) {
    key <- as.character(t_value)
    raw_arm <- binary_arm_nuisance_from_predictions(predictions, t_value)
    iv_hard[[key]] <- binary_iv_arm_endpoint(
      raw_arm, epsilon = 0, compute_gradient = FALSE
    )
    iv_smooth[[key]] <- binary_iv_arm_endpoint(
      raw_arm, epsilon = epsilon, compute_gradient = TRUE
    )
    iv_signals[[key]] <- list(
      lower = binary_arm_one_step_signal(
        data, predictions, raw_arm, iv_smooth[[key]], t_value, "lower"
      ),
      upper = binary_arm_one_step_signal(
        data, predictions, raw_arm, iv_smooth[[key]], t_value, "upper"
      )
    )
    unstabilized_hard[[key]] <- binary_arm_endpoint(
      raw_arm, Lambda0, Lambda1, epsilon = 0, compute_gradient = FALSE
    )
    unstabilized_smooth[[key]] <- binary_arm_endpoint(
      raw_arm, Lambda0, Lambda1, epsilon = epsilon,
      compute_gradient = FALSE
    )
    stabilization[[key]] <- binary_project_arm_outcome_probability(
      raw_arm,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      action = compatibility_action,
      tolerance = compatibility_tolerance
    )
    arm <- stabilization[[key]]$nuisance
    hard[[key]] <- binary_arm_endpoint(
      arm, Lambda0, Lambda1, epsilon = 0, compute_gradient = FALSE
    )
    smooth[[key]] <- binary_arm_endpoint(
      arm, Lambda0, Lambda1, epsilon = epsilon, compute_gradient = TRUE
    )
    signals[[key]] <- list(
      lower = binary_arm_one_step_signal(
        data, predictions, arm, smooth[[key]], t_value, "lower"
      ),
      upper = binary_arm_one_step_signal(
        data, predictions, arm, smooth[[key]], t_value, "upper"
      )
    )
  }

  endpoint_vector <- function(container, use_signal = FALSE) {
    get_value <- function(treatment, direction) {
      value <- container[[as.character(treatment)]][[direction]]
      if (use_signal) mean(value) else mean(value)
    }
    mu1_lower <- get_value(1, "lower")
    mu1_upper <- get_value(1, "upper")
    mu0_lower <- get_value(0, "lower")
    mu0_upper <- get_value(0, "upper")
    c(
      mu1_lower = mu1_lower,
      mu1_upper = mu1_upper,
      mu0_lower = mu0_lower,
      mu0_upper = mu0_upper,
      ate_lower = mu1_lower - mu0_upper,
      ate_upper = mu1_upper - mu0_lower
    )
  }

  hard_values <- endpoint_vector(hard)
  smooth_values <- endpoint_vector(smooth)
  one_step_values <- endpoint_vector(signals, use_signal = TRUE)
  iv_hard_values <- endpoint_vector(iv_hard)
  iv_smooth_values <- endpoint_vector(iv_smooth)
  iv_one_step_values <- endpoint_vector(iv_signals, use_signal = TRUE)
  ate_lower_signal <- signals[["1"]]$lower - signals[["0"]]$upper
  ate_upper_signal <- signals[["1"]]$upper - signals[["0"]]$lower
  standard_errors <- c(
    mu1_lower = stats::sd(signals[["1"]]$lower) / sqrt(nrow(data)),
    mu1_upper = stats::sd(signals[["1"]]$upper) / sqrt(nrow(data)),
    mu0_lower = stats::sd(signals[["0"]]$lower) / sqrt(nrow(data)),
    mu0_upper = stats::sd(signals[["0"]]$upper) / sqrt(nrow(data)),
    ate_lower = stats::sd(ate_lower_signal) / sqrt(nrow(data)),
    ate_upper = stats::sd(ate_upper_signal) / sqrt(nrow(data))
  )
  iv_ate_lower_signal <-
    iv_signals[["1"]]$lower - iv_signals[["0"]]$upper
  iv_ate_upper_signal <-
    iv_signals[["1"]]$upper - iv_signals[["0"]]$lower
  iv_standard_errors <- c(
    mu1_lower = stats::sd(iv_signals[["1"]]$lower) / sqrt(nrow(data)),
    mu1_upper = stats::sd(iv_signals[["1"]]$upper) / sqrt(nrow(data)),
    mu0_lower = stats::sd(iv_signals[["0"]]$lower) / sqrt(nrow(data)),
    mu0_upper = stats::sd(iv_signals[["0"]]$upper) / sqrt(nrow(data)),
    ate_lower = stats::sd(iv_ate_lower_signal) / sqrt(nrow(data)),
    ate_upper = stats::sd(iv_ate_upper_signal) / sqrt(nrow(data))
  )

  long_rows <- function(method, values, se = rep(NA_real_, length(values))) {
    data.frame(
      method = method,
      endpoint = names(values),
      estimate = unname(values),
      standard_error = unname(se[names(values)]),
      stringsAsFactors = FALSE
    )
  }
  estimates <- rbind(
    long_rows("dual_hard_plugin", hard_values),
    long_rows("dual_smooth_plugin", smooth_values),
    long_rows("dual_smooth_one_step", one_step_values, standard_errors),
    long_rows("iv_only_smooth_plugin", iv_smooth_values),
    long_rows(
      "iv_only_smooth_one_step", iv_one_step_values, iv_standard_errors
    )
  )

  comparator_wide <- binary_comparator_plugin_values(predictions, Lambda0)
  comparator_long <- do.call(rbind, lapply(seq_len(nrow(comparator_wide)), function(i) {
    values <- unlist(comparator_wide[i, setdiff(names(comparator_wide), "method")])
    names(values) <- setdiff(names(comparator_wide), "method")
    long_rows(comparator_wide$method[i], values)
  }))
  estimates <- rbind(estimates, comparator_long)

  critical_two_sided <- stats::qnorm(1 - (1 - confidence_level) / 2)
  one_step_rows <- estimates$method == "dual_smooth_one_step"
  intervals <- estimates[one_step_rows, , drop = FALSE]
  intervals$confidence_level <- confidence_level
  intervals$ci_lower <- intervals$estimate -
    critical_two_sided * intervals$standard_error
  intervals$ci_upper <- intervals$estimate +
    critical_two_sided * intervals$standard_error
  ate_confidence_set <- data.frame(
    confidence_level = confidence_level,
    lower = one_step_values["ate_lower"] -
      critical_two_sided * standard_errors["ate_lower"],
    upper = one_step_values["ate_upper"] +
      critical_two_sided * standard_errors["ate_upper"],
    stringsAsFactors = FALSE
  )
  iv_one_step_rows <- estimates$method == "iv_only_smooth_one_step"
  iv_endpoint_intervals <- estimates[iv_one_step_rows, , drop = FALSE]
  iv_endpoint_intervals$confidence_level <- confidence_level
  iv_endpoint_intervals$ci_lower <- iv_endpoint_intervals$estimate -
    critical_two_sided * iv_endpoint_intervals$standard_error
  iv_endpoint_intervals$ci_upper <- iv_endpoint_intervals$estimate +
    critical_two_sided * iv_endpoint_intervals$standard_error
  iv_ate_confidence_set <- data.frame(
    confidence_level = confidence_level,
    lower = iv_one_step_values["ate_lower"] -
      critical_two_sided * iv_standard_errors["ate_lower"],
    upper = iv_one_step_values["ate_upper"] +
      critical_two_sided * iv_standard_errors["ate_upper"],
    stringsAsFactors = FALSE
  )

  diagnostics <- do.call(rbind, lapply(c(0, 1), function(t_value) {
    key <- as.character(t_value)
    raw_hard <- unstabilized_hard[[key]]
    raw_smooth <- unstabilized_smooth[[key]]
    post_hard <- hard[[key]]
    post_smooth <- smooth[[key]]
    hard_gamma_magnitude <- pmax(
      raw_hard$compatible_gamma_lower - raw_hard$gamma_one,
      raw_hard$gamma_one - raw_hard$compatible_gamma_upper,
      raw_hard$compatible_gamma_lower - raw_hard$gamma_zero,
      raw_hard$gamma_zero - raw_hard$compatible_gamma_upper,
      0
    )
    hard_envelope_magnitude <- pmax(
      raw_hard$L_one - raw_hard$U_one,
      raw_hard$L_zero - raw_hard$U_zero,
      0
    )
    hard_theta_magnitude <- pmax(
      raw_hard$theta_lower - raw_hard$theta_upper, 0
    )
    smooth_theta_magnitude <- pmax(
      raw_smooth$theta_lower - raw_smooth$theta_upper, 0
    )
    diagnostic_tolerance <- 1e-8
    post_gamma_violation <- post_hard$gamma_violation
    post_envelope_crossing <- post_hard$envelope_crossing
    post_hard_theta_crossing <- post_hard$theta_crossing
    post_smooth_theta_crossing <- post_smooth$theta_crossing
    data.frame(
      treatment = t_value,
      compatibility_action = compatibility_action,
      compatibility_projection_required_rate = mean(
        stabilization[[key]]$required
      ),
      compatibility_projection_rate = if (compatibility_action == "project") {
        mean(stabilization[[key]]$required)
      } else 0,
      mean_compatibility_projection_magnitude = mean(
        stabilization[[key]]$magnitude
      ),
      maximum_compatibility_projection_magnitude = max(
        stabilization[[key]]$magnitude
      ),
      hard_gamma_violation_rate = mean(raw_hard$gamma_violation),
      hard_envelope_crossing_rate = mean(raw_hard$envelope_crossing),
      hard_theta_crossing_rate = mean(raw_hard$theta_crossing),
      smooth_theta_crossing_rate = mean(raw_smooth$theta_crossing),
      post_hard_gamma_violation_rate = mean(post_gamma_violation),
      post_hard_envelope_crossing_rate = mean(post_envelope_crossing),
      post_hard_theta_crossing_rate = mean(post_hard_theta_crossing),
      post_smooth_theta_crossing_rate = mean(post_smooth_theta_crossing),
      diagnostic_tolerance = diagnostic_tolerance,
      hard_gamma_violation_rate_above_tolerance = mean(
        hard_gamma_magnitude > diagnostic_tolerance
      ),
      hard_envelope_crossing_rate_above_tolerance = mean(
        hard_envelope_magnitude > diagnostic_tolerance
      ),
      hard_theta_crossing_rate_above_tolerance = mean(
        hard_theta_magnitude > diagnostic_tolerance
      ),
      smooth_theta_crossing_rate_above_tolerance = mean(
        smooth_theta_magnitude > diagnostic_tolerance
      ),
      mean_hard_gamma_violation_magnitude = mean(hard_gamma_magnitude),
      maximum_hard_gamma_violation_magnitude = max(hard_gamma_magnitude),
      mean_hard_envelope_crossing_magnitude = mean(
        hard_envelope_magnitude
      ),
      maximum_hard_envelope_crossing_magnitude = max(
        hard_envelope_magnitude
      ),
      mean_hard_theta_crossing_magnitude = mean(hard_theta_magnitude),
      maximum_hard_theta_crossing_magnitude = max(hard_theta_magnitude),
      mean_smooth_theta_crossing_magnitude = mean(smooth_theta_magnitude),
      maximum_smooth_theta_crossing_magnitude = max(smooth_theta_magnitude),
      stringsAsFactors = FALSE
    )
  }))

  list(
    estimates = estimates,
    endpoint_intervals = intervals,
    ate_confidence_set = ate_confidence_set,
    iv_endpoint_intervals = iv_endpoint_intervals,
    iv_ate_confidence_set = iv_ate_confidence_set,
    diagnostics = diagnostics,
    hard = hard,
    smooth = smooth,
    unstabilized_hard = unstabilized_hard,
    unstabilized_smooth = unstabilized_smooth,
    stabilization = stabilization,
    signals = signals,
    iv_hard = iv_hard,
    iv_smooth = iv_smooth,
    iv_signals = iv_signals,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    epsilon = epsilon,
    compatibility_action = compatibility_action,
    compatibility_tolerance = compatibility_tolerance
  )
}

estimate_binary_dual_bounds <- function(
    data,
    x_cols = "X",
    Lambda0 = 2,
    Lambda1 = Lambda0,
    epsilon = 0.01,
    K = 5L,
    fold_id = NULL,
    spline_df = 5L,
    nuisance_learner = c(
      "library_cv", "super_learner", "gam", "glm_spline", "ranger",
      "glmnet"
    ),
    overlap_epsilon = 0.05,
    outcome_probability_epsilon = 1e-4,
    probability_epsilon = NULL,
    gam_basis_dimension = NULL,
    gam_min_basis = 5L,
    gam_max_basis = 20L,
    gam_growth_exponent = 1 / 3,
    nuisance_library_candidates = c(
      "glm_spline", "gam", "glmnet", "ranger"
    ),
    nuisance_library_cv_folds = 3L,
    nuisance_library_selection_rule = c("one_se", "minimum"),
    nuisance_library_ranger_num_trees = 500L,
    nuisance_super_learner_candidates =
      binary_default_super_learner_candidates(),
    nuisance_super_learner_cv_folds = 3L,
    nuisance_super_learner_ranger_num_trees = 500L,
    nuisance_super_learner_hal_max_degree = 2L,
    confidence_level = 0.95,
    compatibility_action = c("project", "none"),
    compatibility_tolerance = 1e-12,
    seed = 1L
) {
  nuisance_learner <- match.arg(nuisance_learner)
  nuisance_library_selection_rule <- match.arg(
    nuisance_library_selection_rule
  )
  compatibility_action <- match.arg(compatibility_action)
  nuisance_fit <- fit_binary_cross_fitted_nuisances(
    data = data,
    x_cols = x_cols,
    K = K,
    fold_id = fold_id,
    spline_df = spline_df,
    learner = nuisance_learner,
    overlap_epsilon = overlap_epsilon,
    outcome_probability_epsilon = outcome_probability_epsilon,
    probability_epsilon = probability_epsilon,
    gam_basis_dimension = gam_basis_dimension,
    gam_min_basis = gam_min_basis,
    gam_max_basis = gam_max_basis,
    gam_growth_exponent = gam_growth_exponent,
    library_candidates = nuisance_library_candidates,
    library_cv_folds = nuisance_library_cv_folds,
    library_selection_rule = nuisance_library_selection_rule,
    library_ranger_num_trees = nuisance_library_ranger_num_trees,
    super_learner_candidates = nuisance_super_learner_candidates,
    super_learner_cv_folds = nuisance_super_learner_cv_folds,
    super_learner_ranger_num_trees =
      nuisance_super_learner_ranger_num_trees,
    super_learner_hal_max_degree = nuisance_super_learner_hal_max_degree,
    seed = seed
  )
  result <- evaluate_binary_bounds_from_nuisances(
    data = data,
    nuisance_fit = nuisance_fit,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    epsilon = epsilon,
    confidence_level = confidence_level,
    compatibility_action = compatibility_action,
    compatibility_tolerance = compatibility_tolerance
  )
  result$nuisance_fit <- nuisance_fit
  result
}

# Optional wrapper around the calibrated estimator in comparator_estimators.R.
# Tan's estimator is deliberately not cross-fitted, matching the methodology
# implemented for the continuous-outcome comparator.  Source
# R/comparator_estimators.R before calling this function.
estimate_binary_tan_cal_pair <- function(
    data,
    Lambda,
    x_cols = "X",
    degree = 2L,
    probability_epsilon = 0.01,
    confidence_level = 0.95
) {
  if (!exists("estimate_tan_cal_bounds", mode = "function")) {
    stop(
      "Source R/comparator_estimators.R before fitting Tan CAL comparators."
    )
  }
  specifications <- list(
    tan_cal_x = x_cols,
    tan_cal_xz = unique(c(x_cols, "Z"))
  )
  fits <- lapply(specifications, function(x_cols) {
    estimate_tan_cal_bounds(
      data = data,
      Lambda = Lambda,
      x_cols = x_cols,
      degree = degree,
      probability_epsilon = probability_epsilon,
      confidence_level = confidence_level
    )
  })

  estimates <- do.call(rbind, lapply(names(fits), function(method) {
    standard_errors <- fits[[method]]$endpoint_standard_errors
    data.frame(
      method = method,
      endpoint = standard_errors$endpoint,
      estimate = standard_errors$estimate,
      standard_error = standard_errors$standard_error,
      stringsAsFactors = FALSE
    )
  }))
  confidence_sets <- do.call(rbind, lapply(names(fits), function(method) {
    interval <- fits[[method]]$ate_sensitivity_interval
    data.frame(
      method = method,
      confidence_level = interval$confidence_level,
      lower = as.numeric(interval$lower),
      upper = as.numeric(interval$upper),
      stringsAsFactors = FALSE
    )
  }))
  diagnostics <- do.call(rbind, lapply(names(fits), function(method) {
    cbind(
      data.frame(method = method, stringsAsFactors = FALSE),
      fits[[method]]$diagnostics
    )
  }))
  rownames(estimates) <- rownames(confidence_sets) <- rownames(diagnostics) <- NULL
  list(
    estimates = estimates,
    ate_confidence_sets = confidence_sets,
    diagnostics = diagnostics,
    fits = fits
  )
}

binary_oracle_nuisances <- function(
    x,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    outcome_intercept = -0.5,
    treatment_log_odds_effect = 0.75,
    outcome_x_strength = 0.75,
    outcome_u_strength = 0.70,
    u_nodes = 120L
) {
  rule <- binary_gauss_legendre_rule(u_nodes)
  u <- rule$nodes
  u_weight <- rule$weights * binary_dtnorm(u)
  u_weight <- u_weight / sum(u_weight)
  x <- as.numeric(x)
  e <- 1 - stats::plogis(0.5 * x)
  arms <- list()
  for (t_value in c(0, 1)) {
    arm <- list()
    outcome_probability <- stats::plogis(outer(
      outcome_intercept + treatment_log_odds_effect * t_value +
        outcome_x_strength * x,
      outcome_u_strength * u,
      "+"
    ))
    for (z_value in c(0, 1)) {
      treatment_probability <- stats::plogis(outer(
        treatment_intercept + x + instrument_strength * z_value,
        confounder_strength * u,
        "+"
      ))
      arm_probability <- if (t_value == 1) {
        treatment_probability
      } else {
        1 - treatment_probability
      }
      r_value <- drop(arm_probability %*% u_weight)
      a_value <- drop((arm_probability * outcome_probability) %*% u_weight) /
        r_value
      arm[[paste0("r", z_value)]] <- r_value
      arm[[paste0("a", z_value)]] <- a_value
    }
    arms[[as.character(t_value)]] <- arm
  }
  list(e = e, arms = arms)
}

binary_population_truth <- function(
    outcome_intercept = -0.5,
    treatment_log_odds_effect = 0.75,
    outcome_x_strength = 0.75,
    outcome_u_strength = 0.70,
    n_nodes = 160L
) {
  rule <- binary_gauss_legendre_rule(n_nodes)
  nodes <- rule$nodes
  weights <- rule$weights * binary_dtnorm(nodes)
  weights <- weights / sum(weights)
  grid <- expand.grid(X = nodes, U = nodes)
  joint_weight <- as.vector(outer(weights, weights))
  probability_Y0 <- stats::plogis(
    outcome_intercept + outcome_x_strength * grid$X +
      outcome_u_strength * grid$U
  )
  probability_Y1 <- stats::plogis(
    outcome_intercept + treatment_log_odds_effect +
      outcome_x_strength * grid$X + outcome_u_strength * grid$U
  )
  mu0 <- sum(joint_weight * probability_Y0)
  mu1 <- sum(joint_weight * probability_Y1)
  data.frame(mu1 = mu1, mu0 = mu0, ate = mu1 - mu0)
}

compute_binary_population_bounds <- function(
    Lambda0 = 2,
    Lambda1 = Lambda0,
    epsilon = 0.01,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    outcome_intercept = -0.5,
    treatment_log_odds_effect = 0.75,
    outcome_x_strength = 0.75,
    outcome_u_strength = 0.70,
    x_nodes = 160L,
    u_nodes = 160L
) {
  x_rule <- binary_gauss_legendre_rule(x_nodes)
  x <- x_rule$nodes
  x_weight <- x_rule$weights * binary_dtnorm(x)
  x_weight <- x_weight / sum(x_weight)
  oracle <- binary_oracle_nuisances(
    x = x,
    confounder_strength = confounder_strength,
    treatment_intercept = treatment_intercept,
    instrument_strength = instrument_strength,
    outcome_intercept = outcome_intercept,
    treatment_log_odds_effect = treatment_log_odds_effect,
    outcome_x_strength = outcome_x_strength,
    outcome_u_strength = outcome_u_strength,
    u_nodes = u_nodes
  )
  hard <- smooth <- list()
  for (t_value in c(0, 1)) {
    key <- as.character(t_value)
    hard[[key]] <- binary_arm_endpoint(
      oracle$arms[[key]], Lambda0, Lambda1, 0, FALSE
    )
    smooth[[key]] <- binary_arm_endpoint(
      oracle$arms[[key]], Lambda0, Lambda1, epsilon, FALSE
    )
  }
  weighted_mean <- function(value) sum(x_weight * value)
  make_dual_row <- function(method, container) {
    mu1_lower <- weighted_mean(container[["1"]]$lower)
    mu1_upper <- weighted_mean(container[["1"]]$upper)
    mu0_lower <- weighted_mean(container[["0"]]$lower)
    mu0_upper <- weighted_mean(container[["0"]]$upper)
    data.frame(
      method = method,
      mu1_lower = mu1_lower, mu1_upper = mu1_upper,
      mu0_lower = mu0_lower, mu0_upper = mu0_upper,
      ate_lower = mu1_lower - mu0_upper,
      ate_upper = mu1_upper - mu0_lower,
      stringsAsFactors = FALSE
    )
  }

  pseudo_predictions <- data.frame(
    e = oracle$e,
    pi0 = oracle$arms[["1"]]$r0,
    pi1 = oracle$arms[["1"]]$r1,
    a_t0_z0 = oracle$arms[["0"]]$a0,
    a_t0_z1 = oracle$arms[["0"]]$a1,
    a_t1_z0 = oracle$arms[["1"]]$a0,
    a_t1_z1 = oracle$arms[["1"]]$a1
  )
  comparator_pointwise <- list()
  for (t_value in c(0, 1)) {
    arm <- binary_arm_nuisance_from_predictions(pseudo_predictions, t_value)
    marginal_r <- oracle$e * arm$r0 + (1 - oracle$e) * arm$r1
    marginal_a <- (
      oracle$e * arm$r0 * arm$a0 +
        (1 - oracle$e) * arm$r1 * arm$a1
    ) / marginal_r
    comparator_pointwise[[as.character(t_value)]] <- list(
      tan_x = binary_standard_msm_arm(marginal_r, marginal_a, Lambda0),
      tan_xz0 = binary_standard_msm_arm(arm$r0, arm$a0, Lambda0),
      tan_xz1 = binary_standard_msm_arm(arm$r1, arm$a1, Lambda0),
      iv = list(
        lower = pmax(arm$r0 * arm$a0, arm$r1 * arm$a1),
        upper = pmin(
          arm$r0 * arm$a0 + 1 - arm$r0,
          arm$r1 * arm$a1 + 1 - arm$r1
        )
      )
    )
  }
  comparator_row <- function(method, component) {
    arm_values <- lapply(c("0", "1"), function(key) {
      if (component == "tan_xz") {
        list(
          lower = oracle$e * comparator_pointwise[[key]]$tan_xz0$lower +
            (1 - oracle$e) * comparator_pointwise[[key]]$tan_xz1$lower,
          upper = oracle$e * comparator_pointwise[[key]]$tan_xz0$upper +
            (1 - oracle$e) * comparator_pointwise[[key]]$tan_xz1$upper
        )
      } else {
        comparator_pointwise[[key]][[component]]
      }
    })
    names(arm_values) <- c("0", "1")
    mu1_lower <- weighted_mean(arm_values[["1"]]$lower)
    mu1_upper <- weighted_mean(arm_values[["1"]]$upper)
    mu0_lower <- weighted_mean(arm_values[["0"]]$lower)
    mu0_upper <- weighted_mean(arm_values[["0"]]$upper)
    data.frame(
      method = method,
      mu1_lower = mu1_lower, mu1_upper = mu1_upper,
      mu0_lower = mu0_lower, mu0_upper = mu0_upper,
      ate_lower = mu1_lower - mu0_upper,
      ate_upper = mu1_upper - mu0_lower,
      stringsAsFactors = FALSE
    )
  }

  bounds <- rbind(
    make_dual_row("dual_hard", hard),
    make_dual_row("dual_smooth", smooth),
    comparator_row("tan_x", "tan_x"),
    comparator_row("tan_xz", "tan_xz"),
    comparator_row("iv_only", "iv")
  )
  bounds$Lambda0 <- Lambda0
  bounds$Lambda1 <- Lambda1
  bounds$epsilon <- ifelse(bounds$method == "dual_smooth", epsilon, 0)

  diagnostics <- data.frame(
    treatment = c(0, 1),
    gamma_violation_rate = c(
      weighted_mean(hard[["0"]]$gamma_violation),
      weighted_mean(hard[["1"]]$gamma_violation)
    ),
    envelope_crossing_rate = c(
      weighted_mean(hard[["0"]]$envelope_crossing),
      weighted_mean(hard[["1"]]$envelope_crossing)
    ),
    theta_crossing_rate = c(
      weighted_mean(hard[["0"]]$theta_crossing),
      weighted_mean(hard[["1"]]$theta_crossing)
    )
  )
  list(
    truth = binary_population_truth(
      outcome_intercept = outcome_intercept,
      treatment_log_odds_effect = treatment_log_odds_effect,
      outcome_x_strength = outcome_x_strength,
      outcome_u_strength = outcome_u_strength,
      n_nodes = max(x_nodes, u_nodes)
    ),
    bounds = bounds,
    diagnostics = diagnostics,
    oracle_nuisances = oracle
  )
}

compute_binary_dgp_sensitivity_frontier <- function(
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    outcome_intercept = -0.5,
    treatment_log_odds_effect = 0.75,
    outcome_x_strength = 0.75,
    outcome_u_strength = 0.70,
    x_grid = seq(-2, 2, length.out = 801L),
    u_nodes = 200L
) {
  rule <- binary_gauss_legendre_rule(u_nodes)
  u <- rule$nodes
  u_weight <- rule$weights * binary_dtnorm(u)
  u_weight <- u_weight / sum(u_weight)
  maximum_ratio <- 1
  maximizing <- NULL
  for (z_value in c(0, 1)) {
    treatment_probability <- stats::plogis(outer(
      treatment_intercept + x_grid + instrument_strength * z_value,
      confounder_strength * u,
      "+"
    ))
    for (t_value in c(0, 1)) {
      outcome_probability <- stats::plogis(outer(
        outcome_intercept + treatment_log_odds_effect * t_value +
          outcome_x_strength * x_grid,
        outcome_u_strength * u,
        "+"
      ))
      observed_weight <- if (t_value == 1) {
        treatment_probability
      } else {
        1 - treatment_probability
      }
      missing_weight <- 1 - observed_weight
      observed_probability <- drop(
        (observed_weight * outcome_probability) %*% u_weight
      ) / drop(observed_weight %*% u_weight)
      missing_probability <- drop(
        (missing_weight * outcome_probability) %*% u_weight
      ) / drop(missing_weight %*% u_weight)
      ratios <- cbind(
        y1 = missing_probability / observed_probability,
        y0 = (1 - missing_probability) / (1 - observed_probability)
      )
      symmetric <- pmax(ratios, 1 / ratios)
      local_index <- which(symmetric == max(symmetric), arr.ind = TRUE)[1L, ]
      local_maximum <- symmetric[local_index[1L], local_index[2L]]
      if (local_maximum > maximum_ratio) {
        maximum_ratio <- local_maximum
        maximizing <- data.frame(
          treatment = t_value,
          instrument = z_value,
          x = x_grid[local_index[1L]],
          outcome = colnames(symmetric)[local_index[2L]],
          ratio = ratios[local_index[1L], local_index[2L]],
          symmetric_ratio = local_maximum,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  list(
    minimum_symmetric_Lambda = unname(maximum_ratio),
    maximizing_configuration = maximizing
  )
}
