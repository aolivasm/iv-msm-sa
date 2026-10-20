# Canonical local-to-tie submodel for one dual-MSM envelope branch.
#
# The lower envelope contains the operation
#
#   max{1 / Lambda0, candidate_lower(P)}.
#
# This file studies the least-favorable scalar submodel
#
#   candidate_lower(P_n) = 1 / Lambda0 + mu_n,
#   mu_n = local_coefficient / sqrt(n),
#
# where a bounded observation V = mu_n + sigma R, R in {-1, 1}, identifies
# mu_n.  The hard target is nonsmooth at the limiting law mu = 0.  The shifted
# log-mean-exp target uses exactly the conservative smoothing implemented in
# compute_projected_envelopes().  Five-fold estimation learns the branch on
# the training observations and evaluates its one-step signal on held-out
# observations.

make_envelope_kink_submodel_spec <- function(
    sample_size,
    sensitivity_lambda = 3,
    envelope_smoothing_epsilon = 0.01,
    local_coefficient = 0.0774596669241483,
    noise_scale = 0.1,
    folds = 5L
) {
  sample_size <- as.integer(sample_size)
  folds <- as.integer(folds)
  values <- c(
    sensitivity_lambda, envelope_smoothing_epsilon,
    local_coefficient, noise_scale
  )
  if (length(sample_size) != 1L || is.na(sample_size) || sample_size < 10L ||
      length(folds) != 1L || is.na(folds) || folds < 2L ||
      folds >= sample_size || any(!is.finite(values)) ||
      sensitivity_lambda <= 1 || envelope_smoothing_epsilon <= 0 ||
      envelope_smoothing_epsilon >= 1 / sensitivity_lambda ||
      local_coefficient <= 0 || noise_scale <= 0) {
    stop("Invalid envelope-kink submodel specification.")
  }
  list(
    name = "local_rademacher_envelope_kink",
    sample_size = sample_size,
    sensitivity_lambda = sensitivity_lambda,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    smoothing_bandwidth = envelope_smoothing_epsilon / log(2),
    local_coefficient = local_coefficient,
    noise_scale = noise_scale,
    folds = folds,
    baseline_branch = 1 / sensitivity_lambda,
    local_mean = local_coefficient / sqrt(sample_size)
  )
}

envelope_kink_shifted_softplus <- function(value, bandwidth) {
  maximum <- pmax(0, value)
  maximum + bandwidth * log(
    (exp(-maximum / bandwidth) +
      exp((value - maximum) / bandwidth)) / 2
  )
}

envelope_kink_targets <- function(specification) {
  mu <- specification$local_mean
  baseline <- specification$baseline_branch
  hard_increment <- max(0, mu)
  smooth_increment <- envelope_kink_shifted_softplus(
    mu, specification$smoothing_bandwidth
  )
  data.frame(
    sample_size = specification$sample_size,
    Lambda = specification$sensitivity_lambda,
    epsilon = specification$envelope_smoothing_epsilon,
    local_coefficient = specification$local_coefficient,
    noise_scale = specification$noise_scale,
    local_mean = mu,
    hard_target = baseline + hard_increment,
    smooth_target = baseline + smooth_increment,
    hard_minus_smooth = hard_increment - smooth_increment,
    stringsAsFactors = FALSE
  )
}

generate_envelope_kink_submodel_data <- function(
    specification,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  specification$local_mean + specification$noise_scale * sample(
    c(-1, 1),
    size = specification$sample_size,
    replace = TRUE
  )
}

estimate_envelope_kink_crossfit <- function(
    value,
    specification,
    fold_id = NULL
) {
  value <- as.numeric(value)
  n <- length(value)
  if (n != specification$sample_size || any(!is.finite(value))) {
    stop("value must be a finite vector with the specified sample size.")
  }
  if (is.null(fold_id)) {
    fold_id <- rep(seq_len(specification$folds), length.out = n)
  }
  if (length(fold_id) != n ||
      !all(sort(unique(fold_id)) == seq_len(specification$folds))) {
    stop("fold_id must contain every requested fold label.")
  }

  baseline <- specification$baseline_branch
  bandwidth <- specification$smoothing_bandwidth
  hard_signal <- smooth_signal <- numeric(n)
  fold_diagnostics <- vector("list", specification$folds)

  for (fold in seq_len(specification$folds)) {
    validation <- fold_id == fold
    training_mean <- mean(value[!validation])

    hard_weight <- as.numeric(training_mean > 0)
    smooth_weight <- stats::plogis(training_mean / bandwidth)
    hard_plugin <- baseline + max(0, training_mean)
    smooth_plugin <- baseline + envelope_kink_shifted_softplus(
      training_mean, bandwidth
    )

    hard_signal[validation] <- hard_plugin +
      hard_weight * (value[validation] - training_mean)
    smooth_signal[validation] <- smooth_plugin +
      smooth_weight * (value[validation] - training_mean)
    fold_diagnostics[[fold]] <- data.frame(
      fold = fold,
      training_mean = training_mean,
      hard_weight = hard_weight,
      smooth_weight = smooth_weight,
      stringsAsFactors = FALSE
    )
  }

  targets <- envelope_kink_targets(specification)
  summarize <- function(signal, target, method, epsilon) {
    estimate <- mean(signal)
    standard_error <- stats::sd(signal) / sqrt(n)
    target_shift_to_sharp <- targets$hard_target - target
    data.frame(
      method = method,
      epsilon = epsilon,
      target = target,
      sharp_reference_target = targets$hard_target,
      target_shift_to_sharp = target_shift_to_sharp,
      estimate = estimate,
      standard_error = standard_error,
      confidence_lower = estimate - stats::qnorm(0.975) * standard_error,
      confidence_upper = estimate + stats::qnorm(0.975) * standard_error,
      covered = estimate - stats::qnorm(0.975) * standard_error <= target &
        estimate + stats::qnorm(0.975) * standard_error >= target,
      z_score = if (standard_error > 0) {
        (estimate - target) / standard_error
      } else NA_real_,
      zero_standard_error = standard_error == 0,
      translated_estimate = estimate + target_shift_to_sharp,
      translated_confidence_lower = estimate + target_shift_to_sharp -
        stats::qnorm(0.975) * standard_error,
      translated_confidence_upper = estimate + target_shift_to_sharp +
        stats::qnorm(0.975) * standard_error,
      translated_sharp_covered = estimate + target_shift_to_sharp -
        stats::qnorm(0.975) * standard_error <= targets$hard_target &
        estimate + target_shift_to_sharp +
          stats::qnorm(0.975) * standard_error >= targets$hard_target,
      stringsAsFactors = FALSE
    )
  }

  list(
    estimates = rbind(
      summarize(hard_signal, targets$hard_target, "hard", 0),
      summarize(
        smooth_signal,
        targets$smooth_target,
        "smooth",
        specification$envelope_smoothing_epsilon
      )
    ),
    fold_diagnostics = do.call(rbind, fold_diagnostics),
    targets = targets
  )
}

run_envelope_kink_submodel_simulation <- function(
    sample_sizes = c(750L, 1500L, 3000L),
    replications = 20000L,
    seed_start = 76001L,
    sensitivity_lambda = 3,
    envelope_smoothing_epsilon = 0.01,
    local_coefficient = 0.0774596669241483,
    noise_scale = 0.1
) {
  sample_sizes <- as.integer(sample_sizes)
  replications <- as.integer(replications)
  if (any(is.na(sample_sizes)) || any(sample_sizes < 10L) ||
      length(replications) != 1L || is.na(replications) ||
      replications < 1L) {
    stop("Invalid simulation sizes.")
  }

  detail_rows <- list()
  target_rows <- list()
  counter <- 1L
  set.seed(as.integer(seed_start))
  for (sample_size in sample_sizes) {
    specification <- make_envelope_kink_submodel_spec(
      sample_size = sample_size,
      sensitivity_lambda = sensitivity_lambda,
      envelope_smoothing_epsilon = envelope_smoothing_epsilon,
      local_coefficient = local_coefficient,
      noise_scale = noise_scale
    )
    target_rows[[length(target_rows) + 1L]] <-
      envelope_kink_targets(specification)
    fold_id <- rep(
      seq_len(specification$folds),
      length.out = specification$sample_size
    )
    for (replication in seq_len(replications)) {
      value <- generate_envelope_kink_submodel_data(
        specification
      )
      fit <- estimate_envelope_kink_crossfit(
        value, specification, fold_id = fold_id
      )
      local <- fit$estimates
      local$sample_size <- sample_size
      local$replication <- replication
      local$random_stream_index <- counter
      local$mean_hard_weight <- mean(
        fit$fold_diagnostics$hard_weight
      )
      local$mean_smooth_weight <- mean(
        fit$fold_diagnostics$smooth_weight
      )
      detail_rows[[counter]] <- local
      counter <- counter + 1L
    }
  }
  detail <- do.call(rbind, detail_rows)
  rownames(detail) <- NULL
  targets <- do.call(rbind, target_rows)
  rownames(targets) <- NULL

  groups <- split(
    detail,
    interaction(
      detail$sample_size, detail$method, detail$epsilon, drop = TRUE
    )
  )
  summary <- do.call(rbind, lapply(groups, function(group) {
    finite_z <- is.finite(group$z_score)
    data.frame(
      sample_size = group$sample_size[1L],
      method = group$method[1L],
      epsilon = group$epsilon[1L],
      replications = nrow(group),
      target = group$target[1L],
      bias = mean(group$estimate - group$target),
      empirical_standard_deviation = stats::sd(group$estimate),
      mean_standard_error = mean(group$standard_error),
      z_mean = mean(group$z_score[finite_z]),
      z_standard_deviation = stats::sd(group$z_score[finite_z]),
      two_sided_95_coverage = mean(group$covered),
      translated_sharp_95_coverage = mean(
        group$translated_sharp_covered
      ),
      zero_standard_error_rate = mean(group$zero_standard_error),
      mean_hard_weight = mean(group$mean_hard_weight),
      mean_smooth_weight = mean(group$mean_smooth_weight),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  summary <- summary[order(
    summary$sample_size, summary$epsilon, summary$method
  ), ]
  list(targets = targets, detail = detail, summary = summary)
}

run_envelope_kink_epsilon_grid_simulation <- function(
    epsilon_grid = c(0.005, 0.01, 0.025, 0.05, 0.10, 0.20),
    ...
) {
  epsilon_grid <- sort(unique(as.numeric(epsilon_grid)))
  if (length(epsilon_grid) < 1L || any(!is.finite(epsilon_grid)) ||
      any(epsilon_grid <= 0)) {
    stop("epsilon_grid must contain positive finite values.")
  }
  runs <- lapply(epsilon_grid, function(epsilon) {
    run_envelope_kink_submodel_simulation(
      envelope_smoothing_epsilon = epsilon,
      ...
    )
  })

  hard_detail <- runs[[1L]]$detail[
    runs[[1L]]$detail$method == "hard", , drop = FALSE
  ]
  smooth_detail <- do.call(rbind, lapply(runs, function(run) {
    run$detail[run$detail$method == "smooth", , drop = FALSE]
  }))
  detail <- rbind(hard_detail, smooth_detail)
  rownames(detail) <- NULL

  hard_summary <- runs[[1L]]$summary[
    runs[[1L]]$summary$method == "hard", , drop = FALSE
  ]
  smooth_summary <- do.call(rbind, lapply(runs, function(run) {
    run$summary[run$summary$method == "smooth", , drop = FALSE]
  }))
  summary <- rbind(hard_summary, smooth_summary)
  summary <- summary[order(
    summary$sample_size, summary$epsilon, summary$method
  ), ]
  rownames(summary) <- NULL

  targets <- do.call(rbind, lapply(runs, `[[`, "targets"))
  rownames(targets) <- NULL
  list(targets = targets, detail = detail, summary = summary)
}
