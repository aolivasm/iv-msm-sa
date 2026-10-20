# Canonical exact-tie submodel for five-fold threshold selection.
#
# This file isolates the weighted check-loss component of a sharp-bound
# threshold learner.  It is not, by itself, the complete IV sharp-bound DGP.
# Its purpose is to give a bounded, analytically transparent example in which
#
#   * Lambda is moderate (the default is Lambda = 4);
#   * every candidate threshold function is a population minimizer;
#   * the selected threshold function has no unique probability limit;
#   * ordinary pooled five-fold Wald inference undercovers; and
#   * predictable online evaluation has a martingale-difference structure.
#
# Let alpha = 1 / (Lambda + 1), X be uniform on S strata, and
#
#   P(Y = y_low) = alpha,   P(Y = y_high) = 1 - alpha.
#
# Every threshold q strictly between y_low and y_high minimizes the population
# alpha-check loss.  Candidate threshold functions use the two values
#
#   q_minus = c - h,  q_plus = c + h,
#   c = (1 - alpha) y_low + alpha y_high,
#
# according to rows of a Hadamard matrix.  Their centered loss signals are
# orthogonal in population.  With 16 candidates, the five-fold Gaussian limit
# reproduces the approximately 82% coverage geometry of a max/min of 16 tied
# one-sided differentiable components.

hadamard_threshold_matrix <- function(order) {
  order <- as.integer(order)
  if (length(order) != 1L || is.na(order) || order < 2L ||
      bitwAnd(order, order - 1L) != 0L) {
    stop("order must be a power of two of at least two.")
  }
  output <- matrix(1, nrow = 1L, ncol = 1L)
  while (nrow(output) < order) {
    output <- rbind(
      cbind(output, output),
      cbind(output, -output)
    )
  }
  output
}

make_hadamard_nonunique_threshold_spec <- function(
    sensitivity_lambda = 4,
    strata = 16L,
    candidates = 16L,
    lower_outcome = 0,
    upper_outcome = 1,
    threshold_separation_fraction = 0.8,
    folds = 5L
) {
  sensitivity_lambda <- as.numeric(sensitivity_lambda)
  strata <- as.integer(strata)
  candidates <- as.integer(candidates)
  folds <- as.integer(folds)
  numeric_values <- c(
    sensitivity_lambda, lower_outcome, upper_outcome,
    threshold_separation_fraction
  )
  if (
    any(!is.finite(numeric_values)) || sensitivity_lambda <= 1 ||
      strata < 2L || bitwAnd(strata, strata - 1L) != 0L ||
      candidates < 2L || candidates > strata ||
      lower_outcome >= upper_outcome ||
      threshold_separation_fraction <= 0 ||
      threshold_separation_fraction >= 1 ||
      folds < 2L
  ) {
    stop("Invalid Hadamard nonunique-threshold specification.")
  }

  alpha <- 1 / (sensitivity_lambda + 1)
  gap <- upper_outcome - lower_outcome
  center <- (1 - alpha) * lower_outcome + alpha * upper_outcome
  half_separation <- threshold_separation_fraction *
    min(alpha, 1 - alpha) * gap
  thresholds <- c(
    minus = center - half_separation,
    plus = center + half_separation
  )
  if (thresholds[["minus"]] <= lower_outcome ||
      thresholds[["plus"]] >= upper_outcome) {
    stop("The two candidate thresholds must lie strictly inside the gap.")
  }

  code <- hadamard_threshold_matrix(strata)[
    seq_len(candidates), , drop = FALSE
  ]
  threshold_matrix <- ifelse(
    code > 0,
    thresholds[["plus"]],
    thresholds[["minus"]]
  )
  dim(threshold_matrix) <- dim(code)
  dimnames(threshold_matrix) <- list(
    candidate = seq_len(candidates),
    stratum = seq_len(strata)
  )

  list(
    name = "hadamard_exact_nonunique_threshold",
    sensitivity_lambda = sensitivity_lambda,
    threshold_level = alpha,
    strata = strata,
    candidates = candidates,
    folds = folds,
    outcomes = c(lower_outcome, upper_outcome),
    threshold_center = center,
    threshold_half_separation = half_separation,
    threshold_separation_fraction = threshold_separation_fraction,
    threshold_matrix = threshold_matrix,
    code = code,
    target = alpha * (1 - alpha) * gap,
    candidate_standard_deviation =
      half_separation * sqrt(alpha * (1 - alpha))
  )
}

hadamard_threshold_check_loss <- function(y, q, alpha) {
  residual <- y - q
  residual * (alpha - as.numeric(residual < 0))
}

hadamard_nonunique_threshold_population_audit <- function(
    specification = make_hadamard_nonunique_threshold_spec()
) {
  x <- rep(seq_len(specification$strata), each = 2L)
  y <- rep(specification$outcomes, times = specification$strata)
  probability <- rep(
    c(
      specification$threshold_level,
      1 - specification$threshold_level
    ) / specification$strata,
    times = specification$strata
  )
  loss <- matrix(
    NA_real_, nrow = length(x), ncol = specification$candidates
  )
  for (candidate in seq_len(specification$candidates)) {
    q <- specification$threshold_matrix[candidate, x]
    loss[, candidate] <- hadamard_threshold_check_loss(
      y, q, specification$threshold_level
    )
  }
  population_mean <- colSums(probability * loss)
  centered <- sweep(loss, 2L, population_mean, "-")
  covariance <- crossprod(sqrt(probability) * centered)
  standard_deviation <- sqrt(diag(covariance))
  correlation <- covariance / tcrossprod(
    standard_deviation, standard_deviation
  )
  list(
    support = data.frame(X = x, Y = y, probability = probability),
    loss = loss,
    population_mean = population_mean,
    covariance = covariance,
    correlation = correlation,
    standard_deviation = standard_deviation,
    maximum_target_error = max(abs(
      population_mean - specification$target
    )),
    maximum_variance_error = max(abs(
      standard_deviation - specification$candidate_standard_deviation
    )),
    maximum_off_diagonal_correlation = max(abs(
      correlation[row(correlation) != col(correlation)]
    ))
  )
}

generate_hadamard_nonunique_threshold_data <- function(
    n,
    specification = make_hadamard_nonunique_threshold_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  x <- sample.int(specification$strata, n, replace = TRUE)
  low <- stats::rbinom(
    n, 1L, specification$threshold_level
  ) == 1L
  y <- ifelse(
    low,
    specification$outcomes[1L],
    specification$outcomes[2L]
  )
  data.frame(X = x, Y = y, low = low)
}

evaluate_hadamard_nonunique_threshold_losses <- function(
    data,
    specification
) {
  if (!is.data.frame(data) || !all(c("X", "Y") %in% names(data))) {
    stop("data must contain X and Y.")
  }
  x <- as.integer(data$X)
  if (any(!x %in% seq_len(specification$strata))) {
    stop("X contains a value outside the specified strata.")
  }
  output <- matrix(
    NA_real_, nrow = nrow(data), ncol = specification$candidates
  )
  for (candidate in seq_len(specification$candidates)) {
    q <- specification$threshold_matrix[candidate, x]
    output[, candidate] <- hadamard_threshold_check_loss(
      data$Y, q, specification$threshold_level
    )
  }
  output
}

hadamard_threshold_select <- function(loss, rows) {
  if (length(rows) < 1L) stop("Threshold selection needs training rows.")
  which.min(colMeans(loss[rows, , drop = FALSE]))
}

estimate_hadamard_threshold_crossfit <- function(
    loss,
    specification,
    fold_id = NULL,
    confidence_level = 0.95
) {
  n <- nrow(loss)
  folds <- specification$folds
  if (is.null(fold_id)) {
    fold_id <- sample(rep(seq_len(folds), length.out = n))
  }
  fold_id <- as.integer(fold_id)
  if (length(fold_id) != n ||
      !identical(sort(unique(fold_id)), seq_len(folds))) {
    stop("fold_id must contain every requested fold.")
  }
  signal <- numeric(n)
  selected <- integer(folds)
  for (fold in seq_len(folds)) {
    validation <- fold_id == fold
    selected[fold] <- hadamard_threshold_select(
      loss, which(!validation)
    )
    rows <- which(validation)
    signal[rows] <- loss[cbind(rows, selected[fold])]
  }
  estimate <- mean(signal)
  standard_error <- stats::sd(signal) / sqrt(n)
  quantile <- stats::qnorm(1 - (1 - confidence_level) / 2)
  list(
    estimate = estimate,
    standard_error = standard_error,
    z = (estimate - specification$target) / standard_error,
    covered = abs(estimate - specification$target) <=
      quantile * standard_error,
    selected = selected,
    distinct_selected = length(unique(selected)),
    signal = signal
  )
}

estimate_hadamard_threshold_online <- function(
    loss,
    specification,
    initial_size = NULL,
    chunks = 10L,
    confidence_level = 0.95
) {
  n <- nrow(loss)
  if (is.null(initial_size)) initial_size <- ceiling(sqrt(n))
  initial_size <- as.integer(initial_size)
  chunks <- as.integer(chunks)
  if (initial_size < 2L || initial_size >= n || chunks < 1L ||
      chunks > n - initial_size) {
    stop("Invalid online initial size or chunk count.")
  }
  evaluation_rows <- seq.int(initial_size + 1L, n)
  chunk_id <- sort(rep(seq_len(chunks), length.out = length(evaluation_rows)))
  signal <- numeric(length(evaluation_rows))
  selected <- integer(chunks)
  cursor <- 1L
  for (chunk in seq_len(chunks)) {
    rows <- evaluation_rows[chunk_id == chunk]
    selected[chunk] <- hadamard_threshold_select(
      loss, seq_len(min(rows) - 1L)
    )
    count <- length(rows)
    indices <- cursor:(cursor + count - 1L)
    signal[indices] <- loss[cbind(rows, selected[chunk])]
    cursor <- cursor + count
  }
  estimate <- mean(signal)
  standard_error <- specification$candidate_standard_deviation /
    sqrt(length(signal))
  quantile <- stats::qnorm(1 - (1 - confidence_level) / 2)
  list(
    estimate = estimate,
    standard_error = standard_error,
    z = (estimate - specification$target) / standard_error,
    covered = abs(estimate - specification$target) <=
      quantile * standard_error,
    selected = selected,
    distinct_selected = length(unique(selected)),
    signal = signal,
    initial_size = initial_size
  )
}

run_hadamard_nonunique_threshold_replication <- function(
    n,
    specification = make_hadamard_nonunique_threshold_spec(),
    seed,
    online_chunks = 10L
) {
  set.seed(as.integer(seed))
  data <- generate_hadamard_nonunique_threshold_data(n, specification)
  loss <- evaluate_hadamard_nonunique_threshold_losses(
    data, specification
  )
  fold_id <- sample(rep(
    seq_len(specification$folds), length.out = as.integer(n)
  ))
  crossfit <- estimate_hadamard_threshold_crossfit(
    loss, specification, fold_id = fold_id
  )
  online <- estimate_hadamard_threshold_online(
    loss, specification, chunks = online_chunks
  )
  data.frame(
    sample_size = as.integer(n),
    seed = as.integer(seed),
    Lambda = specification$sensitivity_lambda,
    estimator = c("five_fold_crossfit", "online_predictable"),
    target = specification$target,
    estimate = c(crossfit$estimate, online$estimate),
    standard_error = c(
      crossfit$standard_error, online$standard_error
    ),
    z = c(crossfit$z, online$z),
    covered = c(crossfit$covered, online$covered),
    distinct_selected = c(
      crossfit$distinct_selected, online$distinct_selected
    ),
    stringsAsFactors = FALSE
  )
}

summarize_hadamard_nonunique_threshold <- function(detail) {
  groups <- split(
    detail,
    interaction(detail$sample_size, detail$estimator, drop = TRUE)
  )
  output <- do.call(rbind, lapply(groups, function(group) {
    data.frame(
      sample_size = group$sample_size[1L],
      Lambda = group$Lambda[1L],
      estimator = group$estimator[1L],
      replications = nrow(group),
      target = group$target[1L],
      bias = mean(group$estimate - group$target),
      empirical_standard_deviation = stats::sd(group$estimate),
      mean_standard_error = mean(group$standard_error),
      z_mean = mean(group$z),
      z_standard_deviation = stats::sd(group$z),
      two_sided_95_coverage = mean(group$covered),
      mean_distinct_selected = mean(group$distinct_selected),
      stringsAsFactors = FALSE
    )
  }))
  rownames(output) <- NULL
  output[order(output$sample_size, output$estimator), , drop = FALSE]
}

simulate_hadamard_threshold_gaussian_limit <- function(
    specification = make_hadamard_nonunique_threshold_spec(),
    replications = 200000L,
    seed = 1L,
    batch_size = 5000L,
    confidence_level = 0.95
) {
  replications <- as.integer(replications)
  batch_size <- as.integer(batch_size)
  if (replications < 1L || batch_size < 1L) {
    stop("replications and batch_size must be positive integers.")
  }
  folds <- specification$folds
  candidates <- specification$candidates
  z <- numeric(replications)
  set.seed(as.integer(seed))
  completed <- 0L
  while (completed < replications) {
    current <- min(batch_size, replications - completed)
    draw <- array(
      stats::rnorm(current * folds * candidates),
      dim = c(current, folds, candidates)
    )
    for (replication in seq_len(current)) {
      local <- draw[replication, , , drop = FALSE]
      dim(local) <- c(folds, candidates)
      total <- colSums(local)
      selected <- vapply(seq_len(folds), function(fold) {
        which.min(total - local[fold, ])
      }, integer(1L))
      z[completed + replication] <- sum(
        local[cbind(seq_len(folds), selected)]
      ) / sqrt(folds)
    }
    completed <- completed + current
  }
  quantile <- stats::qnorm(1 - (1 - confidence_level) / 2)
  data.frame(
    Lambda = specification$sensitivity_lambda,
    folds = folds,
    candidates = candidates,
    replications = replications,
    z_mean = mean(z),
    z_standard_deviation = stats::sd(z),
    two_sided_95_coverage = mean(abs(z) <= quantile),
    stringsAsFactors = FALSE
  )
}
