# Cross-fitted doubly-valid/doubly-sharp (DVDS) comparator
#
# This file implements the estimating equations of Dorn, Guo, and Kallus
# (JASA, 2025) for the conventional marginal sensitivity model.  It is an
# independent implementation: no source code from the authors' replication
# repository is copied or vendored here.
#
# The treatment is T_obs.  The adjustment columns are supplied explicitly so
# the same implementation can represent either the conventional MSM adjusted
# for X or the conventional MSM adjusted for (X,Z).

dvds_check_binary <- function(value, name) {
  observed <- sort(unique(stats::na.omit(as.numeric(value))))
  if (!length(observed) || !all(observed %in% c(0, 1))) {
    stop(name, " must be coded 0/1.")
  }
  invisible(TRUE)
}

dvds_clip_probability <- function(value, epsilon = 0.05) {
  if (length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon <= 0 || epsilon >= 0.5) {
    stop("epsilon must lie strictly between zero and one half.")
  }
  pmin(1 - epsilon, pmax(epsilon, as.numeric(value)))
}

dvds_bind_rows <- function(rows) {
  rows <- Filter(function(value) is.data.frame(value) && nrow(value), rows)
  if (!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(value) {
    for (column in setdiff(columns, names(value))) value[[column]] <- NA
    value[, columns, drop = FALSE]
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

dvds_make_treatment_stratified_folds <- function(
    treatment, K = 5L, seed = 1L
) {
  treatment <- as.integer(treatment)
  dvds_check_binary(treatment, "treatment")
  K <- as.integer(K)
  if (length(K) != 1L || is.na(K) || K < 2L) {
    stop("K must be an integer of at least two.")
  }
  minimum_arm <- min(tabulate(treatment + 1L, nbins = 2L))
  if (minimum_arm < K) stop("Every treatment arm must contain at least K rows.")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  fold <- integer(length(treatment))
  for (arm in 0:1) {
    set.seed(as.integer(seed + 1009L * arm))
    index <- which(treatment == arm)
    labels <- rep(seq_len(K), length.out = length(index))
    fold[index] <- sample(labels, length(labels), replace = FALSE)
  }
  fold
}

dvds_fit_probability_forest <- function(
    x, y, num_trees = 500L, min_node_size = 20L, seed = 1L
) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("The DVDS forest implementation requires package 'ranger'.")
  }
  x <- as.data.frame(x)
  y <- as.integer(y)
  dvds_check_binary(y, "binary forest outcome")
  if (length(unique(y)) == 1L) {
    probability <- (sum(y) + 0.5) / (length(y) + 1)
    return(list(
      predict = function(newdata) rep(probability, nrow(newdata)),
      constant = TRUE,
      probability = probability
    ))
  }
  fit <- ranger::ranger(
    x = x,
    y = factor(y, levels = c(0, 1)),
    probability = TRUE,
    num.trees = as.integer(num_trees),
    min.node.size = as.integer(min_node_size),
    num.threads = 1L,
    seed = as.integer(seed),
    write.forest = TRUE,
    verbose = FALSE
  )
  predict_probability <- function(newdata) {
    prediction <- predict(
      fit, data = as.data.frame(newdata), num.threads = 1L
    )$predictions
    if (is.null(dim(prediction))) return(as.numeric(prediction))
    if ("1" %in% colnames(prediction)) {
      return(as.numeric(prediction[, "1"]))
    }
    stop("The probability forest did not return the probability of class 1.")
  }
  list(predict = predict_probability, fit = fit, constant = FALSE)
}

dvds_fit_quantile_forest <- function(
    x, y, num_trees = 500L, min_node_size = 10L, seed = 1L
) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("The DVDS forest implementation requires package 'ranger'.")
  }
  num_trees <- as.integer(num_trees)
  if (num_trees < 100L) {
    stop("ranger quantile forests require at least 100 trees here.")
  }
  ranger::ranger(
    x = as.data.frame(x),
    y = as.numeric(y),
    quantreg = TRUE,
    keep.inbag = TRUE,
    num.trees = num_trees,
    min.node.size = as.integer(min_node_size),
    num.threads = 1L,
    seed = as.integer(seed),
    write.forest = TRUE,
    verbose = FALSE
  )
}

dvds_predict_quantiles <- function(fit, newdata, probabilities) {
  probabilities <- as.numeric(probabilities)
  prediction <- predict(
    fit,
    data = as.data.frame(newdata),
    type = "quantiles",
    quantiles = probabilities,
    num.threads = 1L
  )$predictions
  if (is.null(dim(prediction))) {
    prediction <- matrix(prediction, ncol = length(probabilities))
  }
  prediction <- as.matrix(prediction)
  if (nrow(prediction) != nrow(newdata) ||
      ncol(prediction) != length(probabilities) ||
      any(!is.finite(prediction))) {
    stop("The quantile forest returned invalid predictions.")
  }
  prediction
}

dvds_fit_regression_weight_forest <- function(
    x, reference_outcome, num_trees = 500L,
    min_node_size = 10L, seed = 1L
) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("The DVDS forest implementation requires package 'ranger'.")
  }
  fit <- ranger::ranger(
    x = as.data.frame(x),
    y = as.numeric(reference_outcome),
    num.trees = as.integer(num_trees),
    min.node.size = as.integer(min_node_size),
    keep.inbag = TRUE,
    num.threads = 1L,
    seed = as.integer(seed),
    write.forest = TRUE,
    verbose = FALSE
  )
  training_nodes <- predict(
    fit, data = as.data.frame(x), type = "terminalNodes", num.threads = 1L
  )$predictions
  list(fit = fit, training_nodes = as.matrix(training_nodes))
}

# Reuse one forest's terminal-node weights to estimate several conditional
# expectations.  This follows the computational strategy used in the authors'
# public simplified implementation, while retaining explicit outer-fold sample
# splitting in this implementation.
dvds_predict_multiple_regressions <- function(
    forest, training_responses, newdata
) {
  training_responses <- as.matrix(training_responses)
  training_nodes <- as.matrix(forest$training_nodes)
  if (nrow(training_nodes) != nrow(training_responses)) {
    stop("The regression-forest response matrix has the wrong number of rows.")
  }
  new_nodes <- predict(
    forest$fit,
    data = as.data.frame(newdata),
    type = "terminalNodes",
    num.threads = 1L
  )$predictions
  new_nodes <- as.matrix(new_nodes)
  if (ncol(new_nodes) != ncol(training_nodes)) {
    stop("Training and prediction terminal-node matrices are incompatible.")
  }
  output <- matrix(0, nrow(new_nodes), ncol(training_responses))
  for (tree in seq_len(ncol(training_nodes))) {
    node <- training_nodes[, tree]
    inbag <- as.numeric(forest$fit$inbag.counts[[tree]])
    if (length(inbag) != nrow(training_responses)) {
      stop("The ranger in-bag counts have the wrong length.")
    }
    weighted_sum <- rowsum(
      training_responses * inbag,
      group = node,
      reorder = FALSE
    )
    weight_sum <- rowsum(inbag, group = node, reorder = FALSE)
    node_mean <- weighted_sum / as.numeric(weight_sum)
    match_index <- match(as.character(new_nodes[, tree]), rownames(node_mean))
    if (anyNA(match_index)) {
      stop("A prediction terminal node contains no in-bag training rows.")
    }
    output <- output + node_mean[match_index, , drop = FALSE]
  }
  output <- output / ncol(training_nodes)
  if (any(!is.finite(output))) {
    stop("The reusable regression forest returned nonfinite predictions.")
  }
  output
}

dvds_inner_training_split <- function(index, seed) {
  index <- as.integer(index)
  if (length(index) < 40L) {
    stop("Too few arm-specific training rows for the DVDS inner split.")
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
  shuffled <- sample(index, length(index), replace = FALSE)
  first <- seq_len(floor(length(shuffled) / 2L))
  list(quantile = shuffled[first], regression = shuffled[-first])
}

dvds_hinge <- function(y, q, Lambda, direction = c("lower", "upper")) {
  direction <- match.arg(direction)
  exponent <- if (direction == "upper") 1 else -1
  q + Lambda^(exponent * sign(y - q)) * (y - q)
}

dvds_uncentered_endpoint <- function(
    y, treatment, arm_probability, q, kappa, Lambda,
    direction = c("lower", "upper"), arm
) {
  direction <- match.arg(direction)
  observed_arm <- as.numeric(treatment == arm)
  hinge <- dvds_hinge(y, q, Lambda, direction)
  observed_arm * y + (1 - observed_arm) * kappa +
    observed_arm * (1 - arm_probability) / arm_probability * (hinge - kappa)
}

dvds_binary_q_kappa <- function(
    outcome_probability, Lambda, direction = c("lower", "upper")
) {
  direction <- match.arg(direction)
  tau <- if (direction == "upper") {
    Lambda / (Lambda + 1)
  } else {
    1 / (Lambda + 1)
  }
  q <- as.numeric(outcome_probability > 1 - tau)
  exponent <- if (direction == "upper") 1 else -1
  kappa <- ifelse(
    q == 1,
    1 - (1 - outcome_probability) * Lambda^(-exponent),
    outcome_probability * Lambda^exponent
  )
  list(q = q, kappa = kappa, tau = tau)
}

dvds_endpoint_keys <- function(Lambdas) {
  rows <- expand.grid(
    Lambda = as.numeric(Lambdas),
    direction = c("lower", "upper"),
    stringsAsFactors = FALSE
  )
  rows <- rows[order(rows$Lambda, match(rows$direction, c("lower", "upper"))), ]
  rownames(rows) <- NULL
  rows$key <- paste(rows$direction, format(rows$Lambda, trim = TRUE), sep = "_")
  rows$tau <- ifelse(
    rows$direction == "upper",
    rows$Lambda / (rows$Lambda + 1),
    1 / (rows$Lambda + 1)
  )
  rows
}

dvds_summarize_signals <- function(
    signals, Lambdas, method = "dvds_crossfit_forest_x"
) {
  endpoint_rows <- list()
  signal_rows <- list()
  for (Lambda in Lambdas) {
    tag <- format(Lambda, trim = TRUE)
    mu1_lower <- signals[[paste(1L, "lower", tag, sep = "_")]]
    mu1_upper <- signals[[paste(1L, "upper", tag, sep = "_")]]
    mu0_lower <- signals[[paste(0L, "lower", tag, sep = "_")]]
    mu0_upper <- signals[[paste(0L, "upper", tag, sep = "_")]]
    local_signals <- list(
      mu1_lower = mu1_lower,
      mu1_upper = mu1_upper,
      mu0_lower = mu0_lower,
      mu0_upper = mu0_upper,
      ate_lower = mu1_lower - mu0_upper,
      ate_upper = mu1_upper - mu0_lower
    )
    signal_rows[[tag]] <- local_signals
    endpoint_rows[[tag]] <- dvds_bind_rows(lapply(
      names(local_signals),
      function(endpoint) {
        value <- local_signals[[endpoint]]
        data.frame(
          method = method,
          Lambda = Lambda,
          endpoint = endpoint,
          estimate = mean(value),
          standard_error = stats::sd(value) / sqrt(length(value)),
          stringsAsFactors = FALSE
        )
      }
    ))
  }
  list(endpoint_detail = dvds_bind_rows(endpoint_rows), signals = signal_rows)
}

estimate_dvds_crossfit_forest <- function(
    data,
    x_cols = paste0("X", 1:5),
    treatment_col = "T_obs",
    outcome_col = "Y",
    outcome_type = c("continuous", "binary"),
    Lambdas = c(1.5, 2, 2.5, 3, 4, 5),
    K = 5L,
    probability_epsilon = 0.05,
    num_trees = 500L,
    probability_min_node_size = 20L,
    quantile_min_node_size = 10L,
    regression_min_node_size = 10L,
    seed = 1L,
    keep_signals = FALSE,
    method = "dvds_crossfit_forest_x",
    adjustment = "X only",
    allow_zero_standard_error = FALSE
) {
  outcome_type <- match.arg(outcome_type)
  required <- c(x_cols, treatment_col, outcome_col)
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Missing DVDS columns: ", paste(missing, collapse = ", "), ".")
  }
  Lambdas <- sort(unique(as.numeric(Lambdas)))
  if (!length(Lambdas) || any(!is.finite(Lambdas)) || any(Lambdas < 1)) {
    stop("Lambdas must be finite and at least one.")
  }
  x <- as.data.frame(data[, x_cols, drop = FALSE])
  if (any(!is.finite(as.matrix(x)))) stop("DVDS covariates must be finite.")
  y <- as.numeric(data[[outcome_col]])
  treatment <- as.integer(data[[treatment_col]])
  dvds_check_binary(treatment, treatment_col)
  if (outcome_type == "binary") dvds_check_binary(y, outcome_col)
  if (any(!is.finite(y))) stop("The DVDS outcome must be finite.")
  fold_id <- dvds_make_treatment_stratified_folds(treatment, K, seed)
  key_grid <- dvds_endpoint_keys(Lambdas)
  signals <- list()
  for (arm in 0:1) {
    for (index in seq_len(nrow(key_grid))) {
      key <- paste(
        arm, key_grid$direction[index],
        format(key_grid$Lambda[index], trim = TRUE), sep = "_"
      )
      signals[[key]] <- rep(NA_real_, nrow(data))
    }
  }
  diagnostics <- list()

  for (fold in seq_len(K)) {
    validation <- fold_id == fold
    training <- !validation
    propensity_fit <- dvds_fit_probability_forest(
      x[training, , drop = FALSE], treatment[training],
      num_trees = num_trees,
      min_node_size = probability_min_node_size,
      seed = seed + 10000L + fold
    )
    propensity_raw <- propensity_fit$predict(x[validation, , drop = FALSE])
    propensity <- dvds_clip_probability(propensity_raw, probability_epsilon)
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      fold = fold,
      component = "treatment_propensity",
      arm = NA_integer_,
      n_training = sum(training),
      n_quantile_training = NA_integer_,
      n_regression_training = NA_integer_,
      minimum_prediction = min(propensity_raw),
      maximum_prediction = max(propensity_raw),
      truncation_rate = mean(propensity != propensity_raw),
      stringsAsFactors = FALSE
    )

    for (arm in 0:1) {
      arm_training <- which(training & treatment == arm)
      p_arm <- if (arm == 1L) propensity else 1 - propensity
      if (outcome_type == "binary") {
        mean_fit <- dvds_fit_probability_forest(
          x[arm_training, , drop = FALSE], y[arm_training],
          num_trees = num_trees,
          min_node_size = probability_min_node_size,
          seed = seed + 20000L + 100L * fold + arm
        )
        mean_raw <- mean_fit$predict(x[validation, , drop = FALSE])
        outcome_probability <- dvds_clip_probability(mean_raw, 1e-6)
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(
          fold = fold,
          component = "binary_outcome_mean",
          arm = arm,
          n_training = length(arm_training),
          n_quantile_training = NA_integer_,
          n_regression_training = NA_integer_,
          minimum_prediction = min(mean_raw),
          maximum_prediction = max(mean_raw),
          truncation_rate = mean(outcome_probability != mean_raw),
          stringsAsFactors = FALSE
        )
        for (index in seq_len(nrow(key_grid))) {
          nuisance <- dvds_binary_q_kappa(
            outcome_probability,
            key_grid$Lambda[index],
            key_grid$direction[index]
          )
          key <- paste(
            arm, key_grid$direction[index],
            format(key_grid$Lambda[index], trim = TRUE), sep = "_"
          )
          signals[[key]][validation] <- dvds_uncentered_endpoint(
            y = y[validation],
            treatment = treatment[validation],
            arm_probability = p_arm,
            q = nuisance$q,
            kappa = nuisance$kappa,
            Lambda = key_grid$Lambda[index],
            direction = key_grid$direction[index],
            arm = arm
          )
        }
      } else {
        inner <- dvds_inner_training_split(
          arm_training,
          seed = seed + 30000L + 100L * fold + arm
        )
        quantile_fit <- dvds_fit_quantile_forest(
          x[inner$quantile, , drop = FALSE], y[inner$quantile],
          num_trees = num_trees,
          min_node_size = quantile_min_node_size,
          seed = seed + 40000L + 100L * fold + arm
        )
        q_regression <- dvds_predict_quantiles(
          quantile_fit,
          x[inner$regression, , drop = FALSE],
          key_grid$tau
        )
        q_validation <- dvds_predict_quantiles(
          quantile_fit,
          x[validation, , drop = FALSE],
          key_grid$tau
        )
        transformed_training <- vapply(
          seq_len(nrow(key_grid)),
          function(index) dvds_hinge(
            y[inner$regression], q_regression[, index],
            key_grid$Lambda[index], key_grid$direction[index]
          ),
          numeric(length(inner$regression))
        )
        if (is.null(dim(transformed_training))) {
          transformed_training <- matrix(transformed_training, ncol = 1L)
        }
        # The public simplified DVDS implementation reuses forest weights
        # across Lambda values within a tail direction.  Keep distinct forest
        # partitions for the lower and upper transformed outcomes.
        kappa_validation <- matrix(
          NA_real_, nrow = sum(validation), ncol = nrow(key_grid)
        )
        for (direction in c("lower", "upper")) {
          direction_index <- which(key_grid$direction == direction)
          reference_index <- direction_index[
            which.min(abs(key_grid$tau[direction_index] - 0.5))[1L]
          ]
          regression_forest <- dvds_fit_regression_weight_forest(
            x[inner$regression, , drop = FALSE],
            transformed_training[, reference_index],
            num_trees = num_trees,
            min_node_size = regression_min_node_size,
            seed = seed + 50000L + 1000L * (direction == "upper") +
              100L * fold + arm
          )
          kappa_validation[, direction_index] <-
            dvds_predict_multiple_regressions(
              regression_forest,
              transformed_training[, direction_index, drop = FALSE],
              x[validation, , drop = FALSE]
            )
          diagnostics[[length(diagnostics) + 1L]] <- data.frame(
            fold = fold,
            component = "continuous_quantile_tail_mean",
            arm = arm,
            direction = direction,
            n_training = length(arm_training),
            n_quantile_training = length(inner$quantile),
            n_regression_training = length(inner$regression),
            minimum_prediction = min(q_validation[, direction_index]),
            maximum_prediction = max(q_validation[, direction_index]),
            truncation_rate = 0,
            stringsAsFactors = FALSE
          )
        }
        # At Lambda=1 both transformed outcomes equal Y.  Use the same fitted
        # conditional mean in both directions so the finite-sample estimator,
        # not only its population limit, collapses to AIPW exactly.
        lambda_one <- which(abs(key_grid$Lambda - 1) < 1e-12)
        if (length(lambda_one) > 1L) {
          kappa_validation[, lambda_one] <-
            kappa_validation[, lambda_one[1L]]
        }
        for (index in seq_len(nrow(key_grid))) {
          key <- paste(
            arm, key_grid$direction[index],
            format(key_grid$Lambda[index], trim = TRUE), sep = "_"
          )
          signals[[key]][validation] <- dvds_uncentered_endpoint(
            y = y[validation],
            treatment = treatment[validation],
            arm_probability = p_arm,
            q = q_validation[, index],
            kappa = kappa_validation[, index],
            Lambda = key_grid$Lambda[index],
            direction = key_grid$direction[index],
            arm = arm
          )
        }
      }
    }
  }

  if (any(vapply(signals, function(value) any(!is.finite(value)), logical(1)))) {
    stop("At least one cross-fitted DVDS contribution is missing or nonfinite.")
  }
  if (length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("method must be one nonempty character string.")
  }
  if (length(adjustment) != 1L || is.na(adjustment) || !nzchar(adjustment)) {
    stop("adjustment must be one nonempty character string.")
  }
  summary <- dvds_summarize_signals(signals, Lambdas, method = method)
  endpoint <- summary$endpoint_detail
  if (any(!is.finite(endpoint$estimate)) ||
      any(!is.finite(endpoint$standard_error)) ||
      any(endpoint$standard_error < 0) ||
      (!isTRUE(allow_zero_standard_error) &&
       any(endpoint$standard_error == 0))) {
    stop("The DVDS endpoint summary contains an invalid estimate or standard error.")
  }
  output <- list(
    endpoint_detail = endpoint,
    diagnostics = dvds_bind_rows(diagnostics),
    fold_id = fold_id,
    settings = list(
      outcome_type = outcome_type,
      x_cols = x_cols,
      treatment_col = treatment_col,
      outcome_col = outcome_col,
      Lambdas = Lambdas,
      folds = as.integer(K),
      probability_epsilon = probability_epsilon,
      num_trees = as.integer(num_trees),
      probability_min_node_size = as.integer(probability_min_node_size),
      quantile_min_node_size = as.integer(quantile_min_node_size),
      regression_min_node_size = as.integer(regression_min_node_size),
      inner_quantile_regression_split = TRUE,
      reusable_tail_mean_forest_weights = TRUE,
      method = method,
      adjustment = adjustment,
      allow_zero_standard_error = isTRUE(allow_zero_standard_error)
    )
  )
  if (isTRUE(keep_signals)) output$signals <- summary$signals
  output
}

# Binary-outcome DVDS implementation using the learner collection employed in
# the authors' empirical application: gradient boosting, radial-basis SVM,
# random forest, and logistic regression.  Each base learner is cross-fitted;
# nonnegative convex ensemble weights minimize Bernoulli log loss on the
# resulting out-of-fold predictions, as in the public application code.
dvds_author_application_learners <- function() {
  c("boosting", "svm", "random_forest", "logistic")
}

dvds_require_author_application_packages <- function() {
  packages <- c("gbm", "e1071", "randomForest")
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    stop(
      "The authors' application learner library requires: ",
      paste(missing, collapse = ", "), "."
    )
  }
  invisible(packages)
}

dvds_fit_author_binary_learner <- function(
    x, y, learner, num_trees = 1000L, gbm_cv_folds = 5L, seed = 1L
) {
  learner <- match.arg(learner, dvds_author_application_learners())
  x <- as.data.frame(x)
  y <- as.integer(y)
  dvds_check_binary(y, "binary learner outcome")
  if (length(unique(y)) == 1L) {
    probability <- (sum(y) + 0.5) / (length(y) + 1)
    return(list(
      learner = learner,
      predict = function(newdata) rep(probability, nrow(newdata)),
      fit = NULL,
      constant = TRUE
    ))
  }

  set.seed(as.integer(seed))
  if (learner == "logistic") {
    design <- cbind(`(Intercept)` = 1, as.matrix(x))
    fit <- suppressWarnings(stats::glm.fit(
      x = design, y = y, family = stats::binomial()
    ))
    coefficients <- fit$coefficients
    coefficients[!is.finite(coefficients)] <- 0
    predict_probability <- function(newdata) {
      new_design <- cbind(`(Intercept)` = 1, as.matrix(newdata))
      stats::plogis(drop(new_design %*% coefficients))
    }
  } else if (learner == "boosting") {
    training <- x
    training$.dvds_outcome <- y
    fit <- gbm::gbm(
      formula = .dvds_outcome ~ .,
      distribution = "bernoulli",
      data = training,
      n.trees = as.integer(num_trees),
      interaction.depth = 2L,
      shrinkage = 0.01,
      bag.fraction = 0.5,
      train.fraction = 1,
      cv.folds = as.integer(gbm_cv_folds),
      n.cores = 1L,
      verbose = FALSE
    )
    best_trees <- suppressWarnings(gbm::gbm.perf(
      fit, plot.it = FALSE, method = "cv"
    ))
    if (!is.finite(best_trees) || best_trees < 1L) {
      best_trees <- as.integer(num_trees)
    }
    predict_probability <- function(newdata) {
      as.numeric(stats::predict(
        fit, newdata = as.data.frame(newdata),
        n.trees = best_trees, type = "response"
      ))
    }
  } else if (learner == "svm") {
    fit <- e1071::svm(
      x = as.matrix(x),
      y = factor(y, levels = c(0, 1)),
      kernel = "radial",
      probability = TRUE,
      type = "C-classification"
    )
    predict_probability <- function(newdata) {
      prediction <- stats::predict(
        fit, as.matrix(newdata), probability = TRUE
      )
      probabilities <- attr(prediction, "probabilities")
      if (is.null(probabilities) || !"1" %in% colnames(probabilities)) {
        stop("The SVM learner did not return class-1 probabilities.")
      }
      as.numeric(probabilities[, "1"])
    }
  } else {
    fit <- randomForest::randomForest(
      x = x,
      y = factor(y, levels = c(0, 1)),
      ntree = as.integer(num_trees),
      nodesize = 1L,
      replace = TRUE,
      na.action = stats::na.omit
    )
    predict_probability <- function(newdata) {
      prediction <- stats::predict(
        fit, newdata = as.data.frame(newdata), type = "prob"
      )
      if (is.null(dim(prediction)) || !"1" %in% colnames(prediction)) {
        stop("The random-forest learner did not return class-1 probabilities.")
      }
      as.numeric(prediction[, "1"])
    }
  }

  list(
    learner = learner,
    predict = function(newdata) {
      value <- predict_probability(newdata)
      if (length(value) != nrow(newdata) || any(!is.finite(value))) {
        stop(learner, " returned invalid probability predictions.")
      }
      value
    },
    fit = fit,
    constant = FALSE
  )
}

dvds_crossfit_author_binary_library <- function(
    x, y, fold_id, training_eligible = rep(TRUE, length(y)),
    num_trees = 1000L, gbm_cv_folds = 5L, seed = 1L
) {
  x <- as.data.frame(x)
  y <- as.integer(y)
  training_eligible <- as.logical(training_eligible)
  if (nrow(x) != length(y) || length(fold_id) != length(y) ||
      length(training_eligible) != length(y)) {
    stop("The cross-fitted learner inputs have incompatible lengths.")
  }
  eligible_outcomes <- y[training_eligible]
  dvds_check_binary(eligible_outcomes, "eligible binary learner outcome")
  learners <- dvds_author_application_learners()
  predictions <- matrix(
    NA_real_, nrow(x), length(learners),
    dimnames = list(NULL, learners)
  )

  for (fold in sort(unique(fold_id))) {
    validation <- fold_id == fold
    training <- !validation & training_eligible
    if (sum(training) < 20L) {
      stop("Too few training observations for one application-library fold.")
    }
    for (learner_index in seq_along(learners)) {
      learner <- learners[learner_index]
      fit <- dvds_fit_author_binary_learner(
        x = x[training, , drop = FALSE],
        y = y[training],
        learner = learner,
        num_trees = num_trees,
        gbm_cv_folds = gbm_cv_folds,
        seed = seed + 1009L * fold + 100003L * learner_index
      )
      predictions[validation, learner] <- fit$predict(
        x[validation, , drop = FALSE]
      )
    }
  }
  if (any(!is.finite(predictions))) {
    stop("The application learner library returned missing predictions.")
  }

  meta_predictions <- cbind(predictions, constant_zero = 0, constant_one = 1)
  meta_x <- meta_predictions[training_eligible, , drop = FALSE]
  meta_y <- y[training_eligible]
  objective <- function(raw_weight) {
    total <- sum(raw_weight)
    if (!is.finite(total) || total <= 0) return(.Machine$double.xmax)
    weight <- raw_weight / total
    probability <- drop(meta_x %*% weight)
    probability <- pmin(1 - 1e-8, pmax(1e-8, probability))
    -mean(meta_y * log(probability) + (1 - meta_y) * log(1 - probability))
  }
  initial <- rep(1 / ncol(meta_predictions), ncol(meta_predictions))
  fit <- stats::optim(
    initial,
    objective,
    method = "L-BFGS-B",
    lower = rep(0, length(initial)),
    control = list(maxit = 2000L, factr = 1e7)
  )
  if (fit$convergence != 0L || any(!is.finite(fit$par)) ||
      sum(fit$par) <= 0) {
    stop("The application-library ensemble-weight optimization failed.")
  }
  weights <- fit$par / sum(fit$par)
  names(weights) <- colnames(meta_predictions)
  combined <- drop(meta_predictions %*% weights)
  if (any(!is.finite(combined))) {
    stop("The application-library ensemble returned invalid predictions.")
  }
  list(
    predictions = combined,
    component_predictions = predictions,
    weights = weights,
    cross_validated_log_loss = objective(weights),
    convergence = fit$convergence
  )
}

estimate_dvds_binary_author_application_library <- function(
    data,
    x_cols,
    treatment_col = "T_obs",
    outcome_col = "Y",
    Lambdas = c(1, 1.2, 1.25, 1.4, 1.5, 1.6, 1.75, 1.8,
                2, 2.5, 3, 3.5, 4, 4.5, 5),
    K = 5L,
    probability_epsilon = 0.05,
    num_trees = 1000L,
    gbm_cv_folds = 5L,
    seed = 1L,
    keep_signals = FALSE,
    method = "dvds_author_application_library_x",
    adjustment = "X only"
) {
  dvds_require_author_application_packages()
  required <- c(x_cols, treatment_col, outcome_col)
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Missing DVDS columns: ", paste(missing, collapse = ", "), ".")
  }
  Lambdas <- sort(unique(as.numeric(Lambdas)))
  if (!length(Lambdas) || any(!is.finite(Lambdas)) || any(Lambdas < 1)) {
    stop("Lambdas must be finite and at least one.")
  }
  x <- as.data.frame(data[, x_cols, drop = FALSE])
  y <- as.integer(data[[outcome_col]])
  treatment <- as.integer(data[[treatment_col]])
  dvds_check_binary(y, outcome_col)
  dvds_check_binary(treatment, treatment_col)
  if (any(!is.finite(as.matrix(x)))) stop("DVDS covariates must be finite.")
  fold_id <- dvds_make_treatment_stratified_folds(treatment, K, seed)

  propensity_fit <- dvds_crossfit_author_binary_library(
    x, treatment, fold_id,
    num_trees = num_trees,
    gbm_cv_folds = gbm_cv_folds,
    seed = seed + 1000000L
  )
  propensity_raw <- propensity_fit$predictions
  propensity <- dvds_clip_probability(propensity_raw, probability_epsilon)
  outcome_fits <- lapply(0:1, function(arm) {
    dvds_crossfit_author_binary_library(
      x, y, fold_id,
      training_eligible = treatment == arm,
      num_trees = num_trees,
      gbm_cv_folds = gbm_cv_folds,
      seed = seed + 2000000L + 100000L * arm
    )
  })
  names(outcome_fits) <- c("0", "1")
  outcome_probability <- lapply(
    outcome_fits,
    function(fit) dvds_clip_probability(fit$predictions, 1e-6)
  )

  signals <- list()
  key_grid <- dvds_endpoint_keys(Lambdas)
  for (arm in 0:1) {
    p_arm <- if (arm == 1L) propensity else 1 - propensity
    for (index in seq_len(nrow(key_grid))) {
      nuisance <- dvds_binary_q_kappa(
        outcome_probability[[as.character(arm)]],
        key_grid$Lambda[index],
        key_grid$direction[index]
      )
      key <- paste(
        arm, key_grid$direction[index],
        format(key_grid$Lambda[index], trim = TRUE), sep = "_"
      )
      signals[[key]] <- dvds_uncentered_endpoint(
        y = y,
        treatment = treatment,
        arm_probability = p_arm,
        q = nuisance$q,
        kappa = nuisance$kappa,
        Lambda = key_grid$Lambda[index],
        direction = key_grid$direction[index],
        arm = arm
      )
    }
  }
  if (any(vapply(signals, function(value) any(!is.finite(value)), logical(1)))) {
    stop("At least one application-library DVDS contribution is nonfinite.")
  }
  summary <- dvds_summarize_signals(signals, Lambdas, method = method)

  fit_list <- c(list(treatment_propensity = propensity_fit), outcome_fits)
  component_names <- c("treatment_propensity", "outcome_mean_0", "outcome_mean_1")
  diagnostics <- dvds_bind_rows(lapply(seq_along(fit_list), function(index) {
    fit <- fit_list[[index]]
    data.frame(
      component = component_names[index],
      learner = names(fit$weights),
      ensemble_weight = unname(fit$weights),
      cross_validated_log_loss = fit$cross_validated_log_loss,
      ensemble_convergence = fit$convergence,
      minimum_prediction = min(fit$predictions),
      maximum_prediction = max(fit$predictions),
      stringsAsFactors = FALSE
    )
  }))
  diagnostics$truncation_rate <- c(
    rep(mean(propensity != propensity_raw), length(propensity_fit$weights)),
    rep(
      mean(outcome_probability[["0"]] != outcome_fits[["0"]]$predictions),
      length(outcome_fits[["0"]]$weights)
    ),
    rep(
      mean(outcome_probability[["1"]] != outcome_fits[["1"]]$predictions),
      length(outcome_fits[["1"]]$weights)
    )
  )

  output <- list(
    endpoint_detail = summary$endpoint_detail,
    diagnostics = diagnostics,
    fold_id = fold_id,
    settings = list(
      outcome_type = "binary",
      x_cols = x_cols,
      treatment_col = treatment_col,
      outcome_col = outcome_col,
      Lambdas = Lambdas,
      folds = as.integer(K),
      probability_epsilon = probability_epsilon,
      learner_library = dvds_author_application_learners(),
      constant_meta_learners = c(0, 1),
      ensemble_loss = "Bernoulli negative log likelihood",
      num_trees = as.integer(num_trees),
      gbm_cv_folds = as.integer(gbm_cv_folds),
      method = method,
      adjustment = adjustment
    )
  )
  if (isTRUE(keep_signals)) output$signals <- summary$signals
  output
}
