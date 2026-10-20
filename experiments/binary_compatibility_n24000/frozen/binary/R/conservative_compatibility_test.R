# Conservative finite-moment compatibility test for the Dual MSM.
#
# Required first:
#   source("R/sharp_bounds_plugin.R")
#
# The test uses only necessary implications of compatibility.  For a finite,
# honestly selected collection of nonnegative functions h(X) and
# r(Y, X) in [0, 1], it tests (i) nonemptiness of the latent lambda interval
# L <= U, algebraically represented by a distinct observable interval for
# gamma_t, and (ii) conditional normalization of the latent lambda envelopes.
# The exact max/min envelopes are represented by their positive-part
# variational identities, so gamma_t is not estimated in the held-out signal.
#
# Inference uses cross-fitted one-step signals for each fixed moment and a
# least-favorable Gaussian multiplier maximum.  The pointwise test is
# asymptotically conservative under the stated nuisance-rate conditions.  A
# second, grid-simultaneous critical value permits reporting the largest value
# on a prespecified Lambda grid that is refuted while controlling familywise
# error across the grid.

check_conservative_compatibility_dependencies <- function() {
  required <- c(
    "check_binary", "clip_probability", "fit_nuisance_models",
    "make_tz_stratified_folds", "make_spline_design"
  )
  missing <- required[
    !vapply(required, exists, logical(1), mode = "function")
  ]
  if (length(missing) > 0L) {
    stop(
      "Source R/sharp_bounds_plugin.R first. Missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_conservative_test_data <- function(data, x_cols) {
  required <- c("Y", "T_obs", "Z", x_cols)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Data omit: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  check_binary(data$T_obs, "T_obs")
  check_binary(data$Z, "Z")
  numeric_columns <- c("Y", x_cols)
  invalid <- vapply(numeric_columns, function(name) {
    !is.numeric(data[[name]]) || any(!is.finite(data[[name]]))
  }, logical(1))
  if (any(invalid)) {
    stop(
      "Y and the test covariates must be finite numeric variables. Invalid: ",
      paste(numeric_columns[invalid], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

make_honest_discovery_split <- function(
    data,
    fraction = 0.25,
    minimum_per_tz_cell = 15L,
    seed = 1L
) {
  if (
    length(fraction) != 1L || !is.finite(fraction) ||
      fraction <= 0 || fraction >= 0.5
  ) {
    stop("fraction must lie strictly between zero and one half.")
  }
  cell <- interaction(data$T_obs, data$Z, drop = TRUE)
  if (any(table(cell) < 2L * minimum_per_tz_cell)) {
    stop("Every treatment-by-instrument cell is too small for honest splitting.")
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
  set.seed(seed)
  discovery <- integer()
  for (level in levels(cell)) {
    index <- which(cell == level)
    count <- max(minimum_per_tz_cell, floor(fraction * length(index)))
    count <- min(count, length(index) - minimum_per_tz_cell)
    discovery <- c(discovery, sample(index, count, replace = FALSE))
  }
  discovery <- sort(discovery)
  list(
    discovery = discovery,
    inference = setdiff(seq_len(nrow(data)), discovery),
    fraction_realized = length(discovery) / nrow(data)
  )
}

make_quantile_cell_spec <- function(
    discovery,
    x_cols,
    resolutions = c(2L, 4L)
) {
  resolutions <- sort(unique(as.integer(resolutions)))
  if (length(resolutions) == 0L || any(resolutions < 2L)) {
    stop("resolutions must contain integers of at least two.")
  }
  cells <- list()
  metadata <- list(data.frame(
    instrument = "one", covariate = "all", resolution = 1L, cell = 1L,
    stringsAsFactors = FALSE
  ))
  for (name in x_cols) {
    for (resolution in resolutions) {
      breaks <- unique(as.numeric(stats::quantile(
        discovery[[name]],
        probs = seq(0, 1, length.out = resolution + 1L),
        type = 8, names = FALSE
      )))
      if (length(breaks) < 3L) next
      breaks[1L] <- -Inf
      breaks[length(breaks)] <- Inf
      key <- paste(name, resolution, sep = "::")
      cells[[key]] <- list(
        covariate = name,
        resolution = resolution,
        breaks = breaks
      )
      for (value in seq_len(length(breaks) - 1L)) {
        metadata[[length(metadata) + 1L]] <- data.frame(
          instrument = paste0(name, "_r", resolution, "_c", value),
          covariate = name,
          resolution = resolution,
          cell = value,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  list(
    cells = cells,
    metadata = do.call(rbind, metadata),
    x_cols = x_cols,
    resolutions = resolutions
  )
}

evaluate_quantile_cell_spec <- function(spec, data) {
  columns <- list(one = rep(1, nrow(data)))
  if (length(spec$cells) > 0L) {
    for (cell_spec in spec$cells) {
      value <- cut(
        data[[cell_spec$covariate]],
        breaks = cell_spec$breaks,
        include.lowest = TRUE,
        labels = FALSE
      )
      for (cell in seq_len(length(cell_spec$breaks) - 1L)) {
        name <- paste0(
          cell_spec$covariate, "_r", cell_spec$resolution, "_c", cell
        )
        columns[[name]] <- as.numeric(value == cell)
      }
    }
  }
  matrix <- do.call(cbind, columns)
  storage.mode(matrix) <- "double"
  matrix
}

make_allocation_spec <- function(
    discovery,
    t_value,
    probabilities = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
) {
  y <- discovery$Y[discovery$T_obs == t_value]
  if (length(y) < 20L) stop("Too few discovery observations in one arm.")
  thresholds <- unique(as.numeric(stats::quantile(
    y, probs = probabilities, type = 8, names = FALSE
  )))
  list(
    treatment = as.integer(t_value),
    thresholds = thresholds,
    probabilities = probabilities
  )
}

evaluate_allocation_spec <- function(spec, data) {
  columns <- list(one = rep(1, nrow(data)))
  for (index in seq_along(spec$thresholds)) {
    threshold <- spec$thresholds[index]
    tag <- gsub("[^A-Za-z0-9]", "_", format(threshold, digits = 8))
    columns[[paste0("lower_", index, "_", tag)]] <-
      as.numeric(data$Y <= threshold)
    columns[[paste0("upper_", index, "_", tag)]] <-
      as.numeric(data$Y > threshold)
  }
  if (!is.null(spec$fixed_functions) && length(spec$fixed_functions) > 0L) {
    for (name in names(spec$fixed_functions)) {
      value <- spec$fixed_functions[[name]](data)
      if (
        length(value) != nrow(data) || any(!is.finite(value)) ||
          any(value < 0 | value > 1)
      ) {
        stop("A fixed allocation function returned values outside [0,1].")
      }
      columns[[name]] <- as.numeric(value)
    }
  }
  matrix <- do.call(cbind, columns)
  storage.mode(matrix) <- "double"
  matrix
}

add_gamma_targeted_allocation_functions <- function(
    discovery,
    allocation_specs,
    probability,
    probability_predict,
    lambda_grid,
    x_cols = "X",
    gamma_lambda_by_t = c("0" = 0.03162278, "1" = 0.01),
    gamma_bandwidth_multiplier = 4,
    gamma_nystrom_rank = 80L,
    seed = 1L
) {
  nuisance <- list(
    e = probability$e,
    pi0 = probability$pi0,
    pi1 = probability$pi1
  )
  output <- allocation_specs
  diagnostics <- list()
  for (t_value in 0:1) {
    rank <- min(
      as.integer(gamma_nystrom_rank),
      sum(discovery$T_obs == t_value)
    )
    fit <- fit_gamma_rkhs_fixed(
      train = discovery,
      nuisance_train = nuisance,
      t_value = t_value,
      x_cols = x_cols,
      lambda = gamma_lambda_by_t[as.character(t_value)],
      bandwidth_multiplier = gamma_bandwidth_multiplier,
      backend = "nystrom",
      nystrom_rank = rank,
      seed = seed + 1000L * t_value
    )
    # Convert the fitted continuous score to fixed threshold rules. The cutoffs
    # are learned in discovery and then frozen; the inference sample never
    # affects their choice. Score quantiles are more robust to calibration bias
    # than using the model boundary itself (for example gamma_hat < 1).
    gamma_discovery <- fit$predict(discovery)$normalized
    score_cutoffs <- unique(as.numeric(stats::quantile(
      gamma_discovery,
      probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
      type = 8,
      names = FALSE
    )))
    for (cutoff_index in seq_along(score_cutoffs)) {
      cutoff <- score_cutoffs[cutoff_index]
      output[[t_value + 1L]]$fixed_functions[[paste0(
        "target_score_below_", cutoff_index
      )]] <- local({
        local_fit <- fit
        local_cutoff <- cutoff
        function(newdata) {
          as.numeric(local_fit$predict(newdata)$normalized <= local_cutoff)
        }
      })
      output[[t_value + 1L]]$fixed_functions[[paste0(
        "target_score_above_", cutoff_index
      )]] <- local({
        local_fit <- fit
        local_cutoff <- cutoff
        function(newdata) {
          as.numeric(local_fit$predict(newdata)$normalized > local_cutoff)
        }
      })
    }
    predict_gamma <- local({
      local_fit <- fit
      function(newdata) local_fit$predict(newdata)$normalized
    })
    fixed <- list()
    for (Lambda in lambda_grid) {
      tag <- gsub("\\.", "p", format(Lambda, scientific = FALSE, trim = TRUE))
      fixed[[paste0("target_gamma_below_L", tag)]] <- local({
        local_Lambda <- Lambda
        local_t <- t_value
        local_predict <- predict_gamma
        local_probability_predict <- probability_predict
        function(newdata) {
          gamma <- local_predict(newdata)
          predicted_probability <- local_probability_predict(newdata)
          interval <- compute_compatible_gamma_interval(
            nuisance = list(
              pi0 = predicted_probability$pi0,
              pi1 = predicted_probability$pi1
            ),
            t_value = local_t,
            Lambda0 = local_Lambda,
            Lambda1 = local_Lambda
          )
          as.numeric(gamma < interval$lower)
        }
      })
      fixed[[paste0("target_gamma_above_L", tag)]] <- local({
        local_Lambda <- Lambda
        local_t <- t_value
        local_predict <- predict_gamma
        local_probability_predict <- probability_predict
        function(newdata) {
          gamma <- local_predict(newdata)
          predicted_probability <- local_probability_predict(newdata)
          interval <- compute_compatible_gamma_interval(
            nuisance = list(
              pi0 = predicted_probability$pi0,
              pi1 = predicted_probability$pi1
            ),
            t_value = local_t,
            Lambda0 = local_Lambda,
            Lambda1 = local_Lambda
          )
          as.numeric(gamma > interval$upper)
        }
      })
      fixed[[paste0("target_lower_branch_L", tag)]] <- local({
        local_Lambda <- Lambda
        local_t <- t_value
        local_predict <- predict_gamma
        local_probability_predict <- probability_predict
        function(newdata) {
          gamma <- local_predict(newdata)
          predicted_probability <- local_probability_predict(newdata)
          p0 <- if (local_t == 1L) {
            predicted_probability$pi0
          } else 1 - predicted_probability$pi0
          p1 <- if (local_t == 1L) {
            predicted_probability$pi1
          } else 1 - predicted_probability$pi1
          primitive <- compute_dual_msm_moment_primitives(
            p0, p1, local_Lambda, local_Lambda
          )
          as.numeric(
            primitive$c_lower * gamma + primitive$b > 1 / local_Lambda
          )
        }
      })
      fixed[[paste0("target_upper_branch_L", tag)]] <- local({
        local_Lambda <- Lambda
        local_t <- t_value
        local_predict <- predict_gamma
        local_probability_predict <- probability_predict
        function(newdata) {
          gamma <- local_predict(newdata)
          predicted_probability <- local_probability_predict(newdata)
          p0 <- if (local_t == 1L) {
            predicted_probability$pi0
          } else 1 - predicted_probability$pi0
          p1 <- if (local_t == 1L) {
            predicted_probability$pi1
          } else 1 - predicted_probability$pi1
          primitive <- compute_dual_msm_moment_primitives(
            p0, p1, local_Lambda, local_Lambda
          )
          as.numeric(
            primitive$c_upper * gamma + primitive$b < local_Lambda
          )
        }
      })
    }
    output[[t_value + 1L]]$fixed_functions <- c(
      output[[t_value + 1L]]$fixed_functions,
      fixed
    )
    diagnostics[[t_value + 1L]] <- data.frame(
      treatment = t_value,
      selected_lambda = gamma_lambda_by_t[as.character(t_value)],
      bandwidth_multiplier = gamma_bandwidth_multiplier,
      nystrom_rank = rank,
      gamma_mean = mean(gamma_discovery),
      gamma_minimum = min(gamma_discovery),
      gamma_maximum = max(gamma_discovery),
      stringsAsFactors = FALSE
    )
  }
  list(specifications = output, diagnostics = do.call(rbind, diagnostics))
}

screen_honest_moment_pairs <- function(
    discovery,
    H,
    allocation_specs,
    minimum_discovery_count = 5L
) {
  if (
    length(minimum_discovery_count) != 1L ||
      !is.finite(minimum_discovery_count) || minimum_discovery_count < 1L
  ) {
    stop("minimum_discovery_count must be a positive integer.")
  }
  rows <- list()
  counter <- 0L
  for (t_value in 0:1) {
    R <- evaluate_allocation_spec(allocation_specs[[t_value + 1L]], discovery)
    for (h_index in seq_len(ncol(H))) {
      for (r_index in seq_len(ncol(R))) {
        counts <- lapply(0:1, function(z_value) {
          cell <- discovery$T_obs == t_value & discovery$Z == z_value
          active <- H[, h_index] > 0
          c(
            one = sum(cell & active & R[, r_index] > 0),
            zero = sum(cell & active & R[, r_index] <= 0),
            total = sum(cell & active)
          )
        })
        counts <- do.call(rbind, counts)
        is_constant <- identical(colnames(R)[r_index], "one")
        eligible <- if (is_constant) {
          all(counts[, "total"] >= 2L * minimum_discovery_count)
        } else {
          all(counts[, "one"] >= minimum_discovery_count) &&
            all(counts[, "zero"] >= minimum_discovery_count)
        }
        counter <- counter + 1L
        rows[[counter]] <- data.frame(
          treatment = t_value,
          h = colnames(H)[h_index],
          r = colnames(R)[r_index],
          eligible = eligible,
          minimum_active_count = min(counts[, "one"]),
          minimum_inactive_count = min(counts[, "zero"]),
          minimum_cell_count = min(counts[, "total"]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

fit_series_conditional_allocations <- function(
    training,
    validation,
    response,
    x_cols,
    spline_df = 5L,
    probability_epsilon = 1e-4
) {
  if (nrow(training) != nrow(response)) {
    stop("The allocation response must align with training rows.")
  }
  output <- matrix(NA_real_, nrow(validation), ncol(response))
  colnames(output) <- colnames(response)
  constant <- apply(response, 2L, function(value) {
    max(value) - min(value) <= 10 * .Machine$double.eps
  })
  if (any(constant)) {
    output[, constant] <- matrix(
      rep(response[1L, constant], each = nrow(validation)),
      nrow = nrow(validation)
    )
  }
  if (any(!constant)) {
    variable_columns <- x_cols[vapply(
      training[, x_cols, drop = FALSE],
      function(value) length(unique(value)) >= 4L,
      logical(1)
    )]
    binary_columns <- setdiff(x_cols, variable_columns)
    training_continuous <- if (length(variable_columns)) {
      make_spline_design(
        data = training,
        x_cols = variable_columns,
        spline_df = spline_df
      )
    } else NULL
    make_design <- function(newdata) {
      continuous <- if (is.null(training_continuous)) {
        matrix(1, nrow(newdata), 1L, dimnames = list(NULL, "intercept"))
      } else training_continuous$predict_matrix(newdata)
      if (!length(binary_columns)) return(continuous)
      cbind(
        continuous,
        as.matrix(newdata[, binary_columns, drop = FALSE])
      )
    }
    design <- list(
      matrix = make_design(training),
      predict_matrix = make_design
    )
    fit <- stats::lm.fit(
      x = design$matrix,
      y = response[, !constant, drop = FALSE]
    )
    coefficients <- fit$coefficients
    coefficients[!is.finite(coefficients)] <- 0
    prediction <- design$predict_matrix(validation) %*% coefficients
    output[, !constant] <- pmin(
      1 - probability_epsilon,
      pmax(probability_epsilon, prediction)
    )
  }
  output
}

crossfit_conservative_test_nuisances <- function(
    data,
    allocation_specs,
    x_cols = "X",
    folds = 5L,
    primitive_learner = c("glm_spline", "gam", "integrated_parametric", "oracle"),
    positivity_epsilon = 0.05,
    treatment_spline_df = 5L,
    gam_basis_dimension = 6L,
    allocation_spline_df = 5L,
    integrated_confounder_strength = 1,
    integrated_quadrature_nodes = 40L,
    seed = 1L
) {
  primitive_learner <- match.arg(primitive_learner)
  n <- nrow(data)
  if (
    primitive_learner == "oracle" &&
      !all(c("e_true", "pi0_true", "pi1_true") %in% names(data))
  ) {
    stop("Oracle primitive nuisances require e_true, pi0_true, and pi1_true.")
  }
  fold_id <- make_tz_stratified_folds(data, K = folds, seed = seed)
  e <- pi0 <- pi1 <- rep(NA_real_, n)
  allocation_observed <- lapply(allocation_specs, evaluate_allocation_spec, data = data)
  allocation_mean <- lapply(seq_along(allocation_specs), function(index) {
    list(
      z0 = matrix(NA_real_, n, ncol(allocation_observed[[index]]),
                  dimnames = list(NULL, colnames(allocation_observed[[index]]))),
      z1 = matrix(NA_real_, n, ncol(allocation_observed[[index]]),
                  dimnames = list(NULL, colnames(allocation_observed[[index]])))
    )
  })
  diagnostics <- vector("list", folds)
  for (fold in seq_len(folds)) {
    training_index <- which(fold_id != fold)
    validation_index <- which(fold_id == fold)
    training <- data[training_index, , drop = FALSE]
    validation <- data[validation_index, , drop = FALSE]
    if (primitive_learner == "oracle") {
      e[validation_index] <- validation$e_true
      pi0[validation_index] <- validation$pi0_true
      pi1[validation_index] <- validation$pi1_true
      diagnostics[[fold]] <- data.frame(
        learner = "oracle", fold = fold, stringsAsFactors = FALSE
      )
    } else {
      fit <- fit_nuisance_models(
        train = training,
        x_cols = x_cols,
        positivity_epsilon = positivity_epsilon,
        learner = primitive_learner,
        treatment_spline_df = treatment_spline_df,
        gam_basis_dimension = gam_basis_dimension,
        integrated_confounder_strength = integrated_confounder_strength,
        integrated_quadrature_nodes = integrated_quadrature_nodes
      )
      prediction <- fit$predict(validation)
      e[validation_index] <- prediction$e
      pi0[validation_index] <- prediction$pi0
      pi1[validation_index] <- prediction$pi1
      local <- fit$diagnostics
      local$fold <- fold
      diagnostics[[fold]] <- local
    }
    for (t_value in 0:1) {
      response_all <- allocation_observed[[t_value + 1L]]
      for (z_value in 0:1) {
        local <- training$T_obs == t_value & training$Z == z_value
        if (sum(local) <= allocation_spline_df * length(x_cols) + 3L) {
          stop("Too few observations for an allocation conditional-mean fit.")
        }
        prediction <- fit_series_conditional_allocations(
          training = training[local, , drop = FALSE],
          validation = validation,
          response = response_all[training_index[local], , drop = FALSE],
          x_cols = x_cols,
          spline_df = allocation_spline_df
        )
        allocation_mean[[t_value + 1L]][[paste0("z", z_value)]][
          validation_index,
        ] <- prediction
      }
    }
  }
  probabilities <- list(
    e = clip_probability(e, positivity_epsilon),
    pi0 = clip_probability(pi0, positivity_epsilon),
    pi1 = clip_probability(pi1, positivity_epsilon)
  )
  if (
    any(!is.finite(unlist(probabilities))) ||
      any(vapply(allocation_mean, function(arm) {
        any(!is.finite(arm$z0)) || any(!is.finite(arm$z1))
      }, logical(1)))
  ) {
    stop("Cross-fitted conservative-test nuisance predictions are incomplete.")
  }
  list(
    probabilities = probabilities,
    allocation_observed = allocation_observed,
    allocation_mean = allocation_mean,
    fold_id = fold_id,
    diagnostics = do.call(rbind, diagnostics),
    primitive_learner = primitive_learner
  )
}

compute_dual_msm_moment_primitives <- function(p0, p1, Lambda0, Lambda1) {
  n_lower <- p0 + (1 - p0) / Lambda0
  d_lower <- p1 + (1 - p1) * Lambda1
  # These are observable bounds on gamma_t implied by nonemptiness of the
  # latent interval [L,U]. They are not the latent envelopes L and U.
  gamma_lower <- n_lower / d_lower
  d_gamma_lower_p0 <- (1 - 1 / Lambda0) / d_lower
  d_gamma_lower_p1 <- n_lower * (Lambda1 - 1) / d_lower^2

  n_upper <- p0 + Lambda0 * (1 - p0)
  d_upper <- p1 + (1 - p1) / Lambda1
  gamma_upper <- n_upper / d_upper
  d_gamma_upper_p0 <- (1 - Lambda0) / d_upper
  d_gamma_upper_p1 <- -n_upper * (1 - 1 / Lambda1) / d_upper^2

  denominator <- 1 - p0
  c_lower <- d_upper / denominator
  c_upper <- d_lower / denominator
  b <- -p0 / denominator
  list(
    gamma_lower = gamma_lower,
    gamma_upper = gamma_upper,
    d_gamma_lower_p0 = d_gamma_lower_p0,
    d_gamma_lower_p1 = d_gamma_lower_p1,
    d_gamma_upper_p0 = d_gamma_upper_p0,
    d_gamma_upper_p1 = d_gamma_upper_p1,
    c_lower = c_lower,
    c_upper = c_upper,
    b = b,
    d_c_lower_p0 = d_upper / denominator^2,
    d_c_lower_p1 = (1 - 1 / Lambda1) / denominator,
    d_c_upper_p0 = d_lower / denominator^2,
    d_c_upper_p1 = (1 - Lambda1) / denominator,
    d_b_p0 = -1 / denominator^2
  )
}

make_arm_lambda_signal_blocks <- function(
    data,
    H,
    R,
    m0,
    m1,
    e,
    pi0,
    pi1,
    t_value,
    Lambda0,
    Lambda1
) {
  n <- nrow(data)
  p0 <- if (t_value == 1L) pi0 else 1 - pi0
  p1 <- if (t_value == 1L) pi1 else 1 - pi1
  primitive <- compute_dual_msm_moment_primitives(
    p0, p1, Lambda0, Lambda1
  )
  residual_p0 <- as.numeric(data$T_obs == t_value) - p0
  residual_p1 <- as.numeric(data$T_obs == t_value) - p1
  score_p0 <- as.numeric(data$Z == 0) / e * residual_p0
  score_p1 <- as.numeric(data$Z == 1) / (1 - e) * residual_p1
  score_m0 <- as.numeric(data$T_obs == t_value & data$Z == 0) /
    (e * p0) * (R - m0)
  score_m1 <- as.numeric(data$T_obs == t_value & data$Z == 1) /
    ((1 - e) * p1) * (R - m1)
  a <- 1 / Lambda0

  build_family <- function(base, dm0, dm1, dp0, dp1, family) {
    unweighted <- base + dm0 * score_m0 + dm1 * score_m1 +
      dp0 * score_p0 + dp1 * score_p1
    blocks <- lapply(seq_len(ncol(unweighted)), function(index) {
      H * unweighted[, index]
    })
    matrix <- do.call(cbind, blocks)
    metadata <- expand.grid(
      h = colnames(H),
      r = colnames(R),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    # lapply above groups all H columns within each R column.
    metadata <- metadata[order(match(metadata$r, colnames(R)),
                               match(metadata$h, colnames(H))), , drop = FALSE]
    metadata$family <- family
    list(signals = matrix, metadata = metadata)
  }

  families <- list(
    pointwise_lower = build_family(
      base = primitive$gamma_lower * m0 - m1,
      dm0 = primitive$gamma_lower,
      dm1 = -1,
      dp0 = primitive$d_gamma_lower_p0 * m0,
      dp1 = primitive$d_gamma_lower_p1 * m0,
      family = "pointwise_lower"
    ),
    pointwise_upper = build_family(
      base = m1 - primitive$gamma_upper * m0,
      dm0 = -primitive$gamma_upper,
      dm1 = 1,
      dp0 = -primitive$d_gamma_upper_p0 * m0,
      dp1 = -primitive$d_gamma_upper_p1 * m0,
      family = "pointwise_upper"
    ),
    normalization_lower = build_family(
      base = primitive$c_lower * m1 + (primitive$b - a) * m0 - (1 - a),
      dm0 = primitive$b - a,
      dm1 = primitive$c_lower,
      dp0 = primitive$d_c_lower_p0 * m1 + primitive$d_b_p0 * m0,
      dp1 = primitive$d_c_lower_p1 * m1,
      family = "normalization_lower"
    ),
    normalization_upper = build_family(
      base = -primitive$c_upper * m1 + (Lambda0 - primitive$b) * m0 -
        (Lambda0 - 1),
      dm0 = Lambda0 - primitive$b,
      dm1 = -primitive$c_upper,
      dp0 = -primitive$d_c_upper_p0 * m1 - primitive$d_b_p0 * m0,
      dp1 = -primitive$d_c_upper_p1 * m1,
      family = "normalization_upper"
    )
  )
  signals <- do.call(cbind, lapply(families, `[[`, "signals"))
  metadata <- do.call(rbind, lapply(families, `[[`, "metadata"))
  metadata$treatment <- t_value
  metadata$Lambda0 <- Lambda0
  metadata$Lambda1 <- Lambda1
  rownames(metadata) <- NULL
  list(signals = signals, metadata = metadata)
}

calibrate_conservative_moment_maximum <- function(
    signals,
    metadata,
    bootstrap_replications = 999L,
    confidence_level = 0.95,
    minimum_standard_deviation = 1e-8,
    seed = 1L
) {
  n <- nrow(signals)
  estimate <- colMeans(signals)
  standard_deviation <- apply(signals, 2L, stats::sd)
  keep <- is.finite(estimate) & is.finite(standard_deviation) &
    standard_deviation > minimum_standard_deviation
  if (!any(keep)) stop("No nondegenerate compatibility moments remain.")
  signals <- signals[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  estimate <- estimate[keep]
  standard_deviation <- standard_deviation[keep]
  z_value <- sqrt(n) * estimate / standard_deviation
  observed_maximum_by_lambda <- tapply(
    z_value,
    interaction(metadata$Lambda0, metadata$Lambda1, drop = TRUE),
    max
  )
  lambda_key <- interaction(metadata$Lambda0, metadata$Lambda1, drop = TRUE)
  lambda_levels <- levels(lambda_key)
  centered_scaled <- sweep(signals, 2L, estimate, "-")
  centered_scaled <- sweep(centered_scaled, 2L, standard_deviation, "/")

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  bootstrap_maximum_by_lambda <- matrix(
    NA_real_, nrow = bootstrap_replications, ncol = length(lambda_levels),
    dimnames = list(NULL, lambda_levels)
  )
  bootstrap_cvm_by_lambda <- bootstrap_maximum_by_lambda
  observed_cvm_by_lambda <- vapply(
    lambda_levels,
    function(level) mean(pmax(z_value[lambda_key == level], 0)^2),
    numeric(1)
  )
  for (replication in seq_len(bootstrap_replications)) {
    multiplier <- stats::rnorm(n)
    draw <- as.vector(crossprod(multiplier, centered_scaled)) / sqrt(n)
    bootstrap_maximum_by_lambda[replication, ] <- vapply(
      lambda_levels,
      function(level) max(draw[lambda_key == level]),
      numeric(1)
    )
    bootstrap_cvm_by_lambda[replication, ] <- vapply(
      lambda_levels,
      function(level) mean(pmax(draw[lambda_key == level], 0)^2),
      numeric(1)
    )
  }
  alpha <- 1 - confidence_level
  component_confidence <- 1 - alpha / 2
  bootstrap_maximum_simultaneous <- apply(
    bootstrap_maximum_by_lambda, 1L, max
  )
  bootstrap_cvm_simultaneous <- apply(bootstrap_cvm_by_lambda, 1L, max)
  maximum_simultaneous_critical <- as.numeric(stats::quantile(
    bootstrap_maximum_simultaneous,
    component_confidence, type = 8, names = FALSE
  ))
  cvm_simultaneous_critical <- as.numeric(stats::quantile(
    bootstrap_cvm_simultaneous,
    component_confidence, type = 8, names = FALSE
  ))
  lambda_rows <- lapply(seq_along(lambda_levels), function(index) {
    local_maximum_bootstrap <- bootstrap_maximum_by_lambda[, index]
    local_cvm_bootstrap <- bootstrap_cvm_by_lambda[, index]
    observed_maximum <- observed_maximum_by_lambda[index]
    observed_cvm <- observed_cvm_by_lambda[index]
    local_maximum_critical <- as.numeric(stats::quantile(
      local_maximum_bootstrap,
      component_confidence, type = 8, names = FALSE
    ))
    local_cvm_critical <- as.numeric(stats::quantile(
      local_cvm_bootstrap,
      component_confidence, type = 8, names = FALSE
    ))
    local_metadata <- metadata[lambda_key == lambda_levels[index], , drop = FALSE]
    data.frame(
      Lambda0 = local_metadata$Lambda0[1L],
      Lambda1 = local_metadata$Lambda1[1L],
      statistic = observed_maximum,
      maximum_statistic = observed_maximum,
      cvm_statistic = observed_cvm,
      maximum_pointwise_critical_value = local_maximum_critical,
      cvm_pointwise_critical_value = local_cvm_critical,
      maximum_pointwise_p_value =
        (1 + sum(local_maximum_bootstrap >= observed_maximum)) /
          (bootstrap_replications + 1),
      cvm_pointwise_p_value =
        (1 + sum(local_cvm_bootstrap >= observed_cvm)) /
          (bootstrap_replications + 1),
      maximum_simultaneous_critical_value =
        maximum_simultaneous_critical,
      cvm_simultaneous_critical_value = cvm_simultaneous_critical,
      maximum_simultaneous_p_value =
        (1 + sum(bootstrap_maximum_simultaneous >= observed_maximum)) /
          (bootstrap_replications + 1),
      cvm_simultaneous_p_value =
        (1 + sum(bootstrap_cvm_simultaneous >= observed_cvm)) /
          (bootstrap_replications + 1),
      moments_tested = nrow(local_metadata),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, lambda_rows)
  summary$pointwise_p_value <- pmin(
    1,
    2 * pmin(
      summary$maximum_pointwise_p_value,
      summary$cvm_pointwise_p_value
    )
  )
  summary$pointwise_reject <- summary$pointwise_p_value <= alpha
  summary$simultaneous_p_value <- pmin(
    1,
    2 * pmin(
      summary$maximum_simultaneous_p_value,
      summary$cvm_simultaneous_p_value
    )
  )
  summary$simultaneous_reject <- summary$simultaneous_p_value <= alpha
  # Backward-compatible aliases refer to the maximum component. The actual
  # decision combines maximum and CvM components by Bonferroni.
  summary$pointwise_critical_value <-
    summary$maximum_pointwise_critical_value
  summary$simultaneous_critical_value <-
    summary$maximum_simultaneous_critical_value
  summary <- summary[order(summary$Lambda0, summary$Lambda1), , drop = FALSE]
  moment_detail <- cbind(
    metadata,
    estimate = estimate,
    standard_error = standard_deviation / sqrt(n),
    z_value = z_value
  )
  moment_detail <- moment_detail[order(moment_detail$z_value, decreasing = TRUE), ]
  list(
    summary = summary,
    moment_detail = moment_detail,
    bootstrap_by_lambda = bootstrap_maximum_by_lambda,
    bootstrap_cvm_by_lambda = bootstrap_cvm_by_lambda,
    bootstrap_simultaneous = bootstrap_maximum_simultaneous,
    bootstrap_cvm_simultaneous = bootstrap_cvm_simultaneous,
    simultaneous_critical_value = maximum_simultaneous_critical,
    cvm_simultaneous_critical_value = cvm_simultaneous_critical
  )
}

select_honest_moment_keys <- function(
    signals,
    metadata,
    eligible_pairs,
    maximum_per_lambda = 16L
) {
  if (
    length(maximum_per_lambda) != 1L || !is.finite(maximum_per_lambda) ||
      maximum_per_lambda < 1L
  ) {
    stop("maximum_per_lambda must be a positive integer.")
  }
  eligibility_key <- paste(
    eligible_pairs$treatment, eligible_pairs$h, eligible_pairs$r, sep = "::"
  )
  moment_pair_key <- paste(metadata$treatment, metadata$h, metadata$r, sep = "::")
  keep <- eligible_pairs$eligible[match(moment_pair_key, eligibility_key)]
  if (anyNA(keep) || !any(keep)) {
    stop("Honest screening did not retain any candidate moments.")
  }
  signals <- signals[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  estimate <- colMeans(signals)
  standard_deviation <- apply(signals, 2L, stats::sd)
  score <- sqrt(nrow(signals)) * estimate / pmax(standard_deviation, 1e-8)
  score[!is.finite(score)] <- -Inf
  lambda_key <- interaction(metadata$Lambda0, metadata$Lambda1, drop = TRUE)
  selected <- unlist(lapply(split(seq_len(nrow(metadata)), lambda_key), function(index) {
    index[order(score[index], decreasing = TRUE)][
      seq_len(min(length(index), as.integer(maximum_per_lambda)))
    ]
  }), use.names = FALSE)
  selected_metadata <- metadata[selected, , drop = FALSE]
  selected_metadata$discovery_estimate <- estimate[selected]
  selected_metadata$discovery_standard_error <-
    standard_deviation[selected] / sqrt(nrow(signals))
  selected_metadata$discovery_z_value <- score[selected]
  selected_metadata$moment_key <- paste(
    selected_metadata$Lambda0,
    selected_metadata$Lambda1,
    selected_metadata$treatment,
    selected_metadata$family,
    selected_metadata$h,
    selected_metadata$r,
    sep = "::"
  )
  selected_metadata
}

test_conservative_dual_msm_compatibility <- function(
    data,
    lambda_grid,
    x_cols = "X",
    discovery_fraction = 0.25,
    x_resolutions = c(2L, 4L),
    y_probabilities = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95),
    minimum_discovery_count = 5L,
    maximum_selected_moments_per_lambda = 16L,
    use_gamma_targeted_allocations = TRUE,
    discovery_gamma_lambda_by_t = c("0" = 0.03162278, "1" = 0.01),
    discovery_gamma_bandwidth_multiplier = 4,
    discovery_gamma_nystrom_rank = 80L,
    folds = 5L,
    primitive_learner = c("glm_spline", "gam", "integrated_parametric", "oracle"),
    positivity_epsilon = 0.05,
    treatment_spline_df = 5L,
    gam_basis_dimension = 6L,
    allocation_spline_df = 5L,
    integrated_confounder_strength = 1,
    integrated_quadrature_nodes = 40L,
    bootstrap_replications = 999L,
    confidence_level = 0.95,
    seed = 1L
) {
  check_conservative_compatibility_dependencies()
  validate_conservative_test_data(data, x_cols)
  primitive_learner <- match.arg(primitive_learner)
  lambda_grid <- sort(unique(as.numeric(lambda_grid)))
  if (length(lambda_grid) == 0L || any(!is.finite(lambda_grid)) ||
      any(lambda_grid < 1)) {
    stop("lambda_grid must contain finite values of at least one.")
  }
  if (bootstrap_replications < 19L) {
    stop("bootstrap_replications must be at least 19.")
  }
  split <- make_honest_discovery_split(
    data, fraction = discovery_fraction, seed = seed
  )
  discovery <- data[split$discovery, , drop = FALSE]
  inference <- data[split$inference, , drop = FALSE]
  h_spec <- make_quantile_cell_spec(
    discovery, x_cols = x_cols, resolutions = x_resolutions
  )
  H <- evaluate_quantile_cell_spec(h_spec, inference)
  H_discovery <- evaluate_quantile_cell_spec(h_spec, discovery)
  allocation_specs <- lapply(0:1, function(t_value) {
    make_allocation_spec(
      discovery, t_value = t_value, probabilities = y_probabilities
    )
  })
  gamma_targeting_diagnostics <- data.frame()
  if (isTRUE(use_gamma_targeted_allocations)) {
    discovery_probability_fit <- fit_nuisance_models(
      train = discovery,
      x_cols = x_cols,
      positivity_epsilon = positivity_epsilon,
      learner = if (primitive_learner == "oracle") {
        "integrated_parametric"
      } else primitive_learner,
      treatment_spline_df = treatment_spline_df,
      gam_basis_dimension = gam_basis_dimension,
      integrated_confounder_strength = integrated_confounder_strength,
      integrated_quadrature_nodes = integrated_quadrature_nodes
    )
    discovery_probability_predict <- function(newdata) {
      discovery_probability_fit$predict(newdata)
    }
    targeting <- add_gamma_targeted_allocation_functions(
      discovery = discovery,
      allocation_specs = allocation_specs,
      probability = discovery_probability_predict(discovery),
      probability_predict = discovery_probability_predict,
      lambda_grid = lambda_grid,
      x_cols = x_cols,
      gamma_lambda_by_t = discovery_gamma_lambda_by_t,
      gamma_bandwidth_multiplier = discovery_gamma_bandwidth_multiplier,
      gamma_nystrom_rank = discovery_gamma_nystrom_rank,
      seed = seed + 2500L
    )
    allocation_specs <- targeting$specifications
    gamma_targeting_diagnostics <- targeting$diagnostics
  }
  eligible_pairs <- screen_honest_moment_pairs(
    discovery = discovery,
    H = H_discovery,
    allocation_specs = allocation_specs,
    minimum_discovery_count = minimum_discovery_count
  )
  discovery_nuisance <- crossfit_conservative_test_nuisances(
    data = discovery,
    allocation_specs = allocation_specs,
    x_cols = x_cols,
    folds = min(3L, folds),
    primitive_learner = primitive_learner,
    positivity_epsilon = positivity_epsilon,
    treatment_spline_df = treatment_spline_df,
    gam_basis_dimension = gam_basis_dimension,
    allocation_spline_df = allocation_spline_df,
    integrated_confounder_strength = integrated_confounder_strength,
    integrated_quadrature_nodes = integrated_quadrature_nodes,
    seed = seed + 5000L
  )
  discovery_blocks <- list()
  discovery_counter <- 0L
  for (Lambda in lambda_grid) {
    for (t_value in 0:1) {
      discovery_counter <- discovery_counter + 1L
      discovery_blocks[[discovery_counter]] <- make_arm_lambda_signal_blocks(
        data = discovery,
        H = H_discovery,
        R = discovery_nuisance$allocation_observed[[t_value + 1L]],
        m0 = discovery_nuisance$allocation_mean[[t_value + 1L]]$z0,
        m1 = discovery_nuisance$allocation_mean[[t_value + 1L]]$z1,
        e = discovery_nuisance$probabilities$e,
        pi0 = discovery_nuisance$probabilities$pi0,
        pi1 = discovery_nuisance$probabilities$pi1,
        t_value = t_value,
        Lambda0 = Lambda,
        Lambda1 = Lambda
      )
    }
  }
  selected_moments <- select_honest_moment_keys(
    signals = do.call(cbind, lapply(discovery_blocks, `[[`, "signals")),
    metadata = do.call(rbind, lapply(discovery_blocks, `[[`, "metadata")),
    eligible_pairs = eligible_pairs,
    maximum_per_lambda = maximum_selected_moments_per_lambda
  )
  nuisance <- crossfit_conservative_test_nuisances(
    data = inference,
    allocation_specs = allocation_specs,
    x_cols = x_cols,
    folds = folds,
    primitive_learner = primitive_learner,
    positivity_epsilon = positivity_epsilon,
    treatment_spline_df = treatment_spline_df,
    gam_basis_dimension = gam_basis_dimension,
    allocation_spline_df = allocation_spline_df,
    integrated_confounder_strength = integrated_confounder_strength,
    integrated_quadrature_nodes = integrated_quadrature_nodes,
    seed = seed + 10000L
  )
  probability <- nuisance$probabilities
  blocks <- list()
  counter <- 0L
  for (Lambda in lambda_grid) {
    for (t_value in 0:1) {
      counter <- counter + 1L
      blocks[[counter]] <- make_arm_lambda_signal_blocks(
        data = inference,
        H = H,
        R = nuisance$allocation_observed[[t_value + 1L]],
        m0 = nuisance$allocation_mean[[t_value + 1L]]$z0,
        m1 = nuisance$allocation_mean[[t_value + 1L]]$z1,
        e = probability$e,
        pi0 = probability$pi0,
        pi1 = probability$pi1,
        t_value = t_value,
        Lambda0 = Lambda,
        Lambda1 = Lambda
      )
    }
  }
  signals <- do.call(cbind, lapply(blocks, `[[`, "signals"))
  metadata <- do.call(rbind, lapply(blocks, `[[`, "metadata"))
  eligibility_key <- paste(
    eligible_pairs$treatment, eligible_pairs$h, eligible_pairs$r, sep = "::"
  )
  moment_key <- paste(metadata$treatment, metadata$h, metadata$r, sep = "::")
  keep <- eligible_pairs$eligible[match(moment_key, eligibility_key)]
  if (anyNA(keep) || !any(keep)) {
    stop("Honest screening did not retain any conservative moments.")
  }
  signals <- signals[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  moment_key <- paste(
    metadata$Lambda0, metadata$Lambda1, metadata$treatment, metadata$family,
    metadata$h, metadata$r, sep = "::"
  )
  selected_match <- match(moment_key, selected_moments$moment_key)
  core_moment <- metadata$h == "one" |
    (metadata$r == "one" & grepl("^normalization_", metadata$family))
  keep_selected <- !is.na(selected_match) | core_moment
  if (!any(keep_selected)) {
    stop("No honestly selected moment could be matched in the inference sample.")
  }
  signals <- signals[, keep_selected, drop = FALSE]
  metadata <- metadata[keep_selected, , drop = FALSE]
  metadata$discovery_z_value <- NA_real_
  matched <- !is.na(selected_match[keep_selected])
  metadata$discovery_z_value[matched] <- selected_moments$discovery_z_value[
    selected_match[keep_selected][matched]
  ]
  calibration <- calibrate_conservative_moment_maximum(
    signals = signals,
    metadata = metadata,
    bootstrap_replications = bootstrap_replications,
    confidence_level = confidence_level,
    seed = seed + 900000L
  )
  rejected <- calibration$summary$Lambda0[
    calibration$summary$simultaneous_reject
  ]
  list(
    summary = calibration$summary,
    moment_detail = calibration$moment_detail,
    largest_simultaneously_refuted_grid_value = if (length(rejected)) {
      max(rejected)
    } else NA_real_,
    split = split,
    instrument_specification = list(
      h = h_spec,
      allocation = allocation_specs,
      honest_pair_screen = eligible_pairs,
      selected_moments = selected_moments
    ),
    gamma_targeting_diagnostics = gamma_targeting_diagnostics,
    nuisance = nuisance,
    calibration = calibration,
    settings = list(
      lambda_grid = lambda_grid,
      x_cols = x_cols,
      discovery_fraction = discovery_fraction,
      minimum_discovery_count = minimum_discovery_count,
      maximum_selected_moments_per_lambda =
        maximum_selected_moments_per_lambda,
      use_gamma_targeted_allocations = use_gamma_targeted_allocations,
      inference_sample_size = nrow(inference),
      primitive_learner = primitive_learner,
      bootstrap_replications = bootstrap_replications,
      confidence_level = confidence_level,
      inferential_scope = paste(
        "Conservative finite-moment test; rejection refutes compatibility,",
        "while nonrejection does not establish it."
      )
    )
  )
}
