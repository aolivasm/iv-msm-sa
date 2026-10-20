# Simulation-only numerical oracle for auditing the dual-MSM EIF.
#
# This file assumes that R/sharp_bounds_plugin.R has already been sourced. It
# uses the known scalar-X DGP to integrate the threshold and every conditional
# EIF expectation. These helpers are diagnostics, not public-data estimators.

audit_trapezoid_rule <- function(lower, upper, n_nodes) {
  if (n_nodes < 3L) stop("n_nodes must be at least three.")
  nodes <- seq(lower, upper, length.out = n_nodes)
  step <- nodes[2L] - nodes[1L]
  weights <- rep(step, n_nodes)
  weights[c(1L, n_nodes)] <- step / 2
  list(nodes = nodes, weights = weights)
}

audit_build_conditional_outcome_grid <- function(
    validation,
    t_value,
    x_col = "X",
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    y_nodes = 401L,
    full_interval_nodes = NULL,
    q_lower = NULL,
    q_upper = NULL,
    interval_nodes = 40L,
    u_nodes = 40L,
    batch_size = 25L
) {
  if (length(x_col) != 1L || !x_col %in% names(validation)) {
    stop("The DGP-integrated EIF audit requires one scalar X column.")
  }
  y_nodes <- as.integer(y_nodes)
  if (!is.null(full_interval_nodes)) {
    full_interval_nodes <- as.integer(full_interval_nodes)
  }
  interval_nodes <- as.integer(interval_nodes)
  u_nodes <- as.integer(u_nodes)
  batch_size <- as.integer(batch_size)
  if (
    is.null(full_interval_nodes) &&
      (y_nodes < 101L || y_nodes %% 2L != 1L)
  ) {
    stop("y_nodes must be an odd integer of at least 101.")
  }
  if (
    !is.null(full_interval_nodes) &&
      (length(full_interval_nodes) != 1L ||
        is.na(full_interval_nodes) ||
        full_interval_nodes < 8L)
  ) {
    stop("full_interval_nodes must be NULL or an integer of at least eight.")
  }
  if (u_nodes < 10L || batch_size < 1L) {
    stop("u_nodes must be at least 10 and batch_size must be positive.")
  }

  x <- as.numeric(validation[[x_col]])
  n <- length(x)
  piecewise <- !is.null(q_lower) || !is.null(q_upper)
  if (xor(is.null(q_lower), is.null(q_upper))) {
    stop("q_lower and q_upper must either both be supplied or both be NULL.")
  }
  if (piecewise) {
    if (
      length(q_lower) != n ||
        length(q_upper) != n ||
        interval_nodes < 4L
    ) {
      stop("Piecewise outcome quadrature inputs are invalid.")
    }
    center <- treatment_effect * t_value + 1.5 * x
    lower_relative <- pmin(5, pmax(-5, q_lower - center))
    upper_relative <- pmin(5, pmax(-5, q_upper - center))
    first_threshold <- pmin(lower_relative, upper_relative)
    second_threshold <- pmax(lower_relative, upper_relative)
    boundaries <- cbind(
      -5,
      first_threshold,
      second_threshold,
      5
    )
    standard_rule <- gauss_legendre_rule(interval_nodes, -1, 1)
    relative_y <- matrix(
      NA_real_, nrow = n, ncol = 3L * interval_nodes
    )
    y_weights <- matrix(
      NA_real_, nrow = n, ncol = 3L * interval_nodes
    )
    for (segment in seq_len(3L)) {
      columns <- (segment - 1L) * interval_nodes +
        seq_len(interval_nodes)
      lower <- boundaries[, segment]
      upper <- boundaries[, segment + 1L]
      half_width <- (upper - lower) / 2
      midpoint <- (upper + lower) / 2
      relative_y[, columns] <- matrix(
        midpoint,
        nrow = n,
        ncol = interval_nodes
      ) + outer(half_width, standard_rule$nodes)
      y_weights[, columns] <- outer(
        half_width,
        standard_rule$weights
      )
    }
    y_nodes <- ncol(relative_y)
    quadrature_type <- "piecewise_gauss_legendre"
  } else if (!is.null(full_interval_nodes)) {
    relative_rule <- gauss_legendre_rule(full_interval_nodes, -5, 5)
    y_nodes <- full_interval_nodes
    relative_y <- matrix(
      relative_rule$nodes,
      nrow = n,
      ncol = y_nodes,
      byrow = TRUE
    )
    y_weights <- matrix(
      relative_rule$weights,
      nrow = n,
      ncol = y_nodes,
      byrow = TRUE
    )
    quadrature_type <- "full_gauss_legendre"
  } else {
    relative_rule <- audit_trapezoid_rule(-5, 5, y_nodes)
    relative_y <- matrix(
      relative_rule$nodes,
      nrow = n,
      ncol = y_nodes,
      byrow = TRUE
    )
    y_weights <- matrix(
      relative_rule$weights,
      nrow = n,
      ncol = y_nodes,
      byrow = TRUE
    )
    quadrature_type <- "uniform_trapezoid"
  }
  probability0 <- matrix(NA_real_, nrow = n, ncol = y_nodes)
  probability1 <- matrix(NA_real_, nrow = n, ncol = y_nodes)
  raw_mass0 <- raw_mass1 <- rep(NA_real_, n)

  u_rule <- gauss_legendre_rule(u_nodes, -2, 2)
  u_weights <- u_rule$weights * d_u(u_rule$nodes)
  density_rule <- gauss_legendre_rule(u_nodes, -1, 1)

  for (start in seq.int(1L, n, by = batch_size)) {
    finish <- min(n, start + batch_size - 1L)
    index <- start:finish
    local_x <- x[index]
    local_n <- length(index)
    x_u <- matrix(local_x, nrow = local_n, ncol = u_nodes)
    u_full <- matrix(
      u_rule$nodes,
      nrow = local_n,
      ncol = u_nodes,
      byrow = TRUE
    )
    eta0 <- treatment_intercept + x_u + confounder_strength * u_full
    eta1 <- eta0 + instrument_strength
    p_t1_z0 <- stats::plogis(eta0)
    p_t1_z1 <- stats::plogis(eta1)
    pi0 <- as.vector(p_t1_z0 %*% u_weights)
    pi1 <- as.vector(p_t1_z1 %*% u_weights)
    r0 <- if (t_value == 1) pi0 else 1 - pi0
    r1 <- if (t_value == 1) pi1 else 1 - pi1

    outcome_center <- treatment_effect * t_value + 1.5 * local_x
    local_relative_y <- relative_y[index, , drop = FALSE]
    grid_x <- rep(local_x, each = y_nodes)
    grid_center <- rep(outcome_center, each = y_nodes)
    grid_y <- grid_center + as.vector(t(local_relative_y))
    grid_n <- length(grid_y)

    u_lower <- pmax(-2, (grid_y - grid_center - 2) / 1.5)
    u_upper <- pmin(2, (grid_y - grid_center + 2) / 1.5)
    half_width <- pmax(0, (u_upper - u_lower) / 2)
    midpoint <- (u_upper + u_lower) / 2
    u_adaptive <- matrix(midpoint, nrow = grid_n, ncol = u_nodes) +
      outer(half_width, density_rule$nodes)
    adaptive_weights <- outer(
      half_width,
      density_rule$weights
    ) * d_u(u_adaptive)
    x_adaptive <- matrix(grid_x, nrow = grid_n, ncol = u_nodes)
    y_adaptive <- matrix(grid_y, nrow = grid_n, ncol = u_nodes)
    outcome_mean <- treatment_effect * t_value +
      1.5 * x_adaptive + 1.5 * u_adaptive
    density_y_given_xu <- dtnormtc(
      y_adaptive,
      mu = outcome_mean,
      s = 1
    )
    density_eta0 <- treatment_intercept +
      x_adaptive + confounder_strength * u_adaptive
    density_eta1 <- density_eta0 + instrument_strength
    density_p1_z0 <- stats::plogis(density_eta0)
    density_p1_z1 <- stats::plogis(density_eta1)
    treatment_weight0 <- if (t_value == 1) {
      density_p1_z0
    } else {
      1 - density_p1_z0
    }
    treatment_weight1 <- if (t_value == 1) {
      density_p1_z1
    } else {
      1 - density_p1_z1
    }
    joint0 <- rowSums(
      density_y_given_xu * treatment_weight0 * adaptive_weights
    )
    joint1 <- rowSums(
      density_y_given_xu * treatment_weight1 * adaptive_weights
    )
    density0 <- joint0 / rep(pmax(r0, 1e-12), each = y_nodes)
    density1 <- joint1 / rep(pmax(r1, 1e-12), each = y_nodes)
    density0 <- matrix(density0, nrow = local_n, byrow = TRUE)
    density1 <- matrix(density1, nrow = local_n, byrow = TRUE)
    local_y_weights <- y_weights[index, , drop = FALSE]
    local_probability0 <- density0 * local_y_weights
    local_probability1 <- density1 * local_y_weights
    local_mass0 <- rowSums(local_probability0)
    local_mass1 <- rowSums(local_probability1)
    if (
      any(!is.finite(local_mass0)) ||
        any(!is.finite(local_mass1)) ||
        any(local_mass0 <= 0) ||
        any(local_mass1 <= 0)
    ) {
      stop("The DGP conditional outcome quadrature has invalid mass.")
    }
    probability0[index, ] <- local_probability0 / local_mass0
    probability1[index, ] <- local_probability1 / local_mass1
    raw_mass0[index] <- local_mass0
    raw_mass1[index] <- local_mass1
  }

  list(
    x = x,
    relative_y = relative_y,
    probability0 = probability0,
    probability1 = probability1,
    raw_mass0 = raw_mass0,
    raw_mass1 = raw_mass1,
    settings = list(
      treatment = t_value,
      x_col = x_col,
      treatment_effect = treatment_effect,
      y_nodes = y_nodes,
      full_interval_nodes = if (
        identical(quadrature_type, "full_gauss_legendre")
      ) {
        y_nodes
      } else {
        NA_integer_
      },
      interval_nodes = if (piecewise) interval_nodes else NA_integer_,
      quadrature_type = quadrature_type,
      u_nodes = u_nodes,
      batch_size = batch_size
    )
  )
}

make_dgp_oracle_gamma_normalizer_factory <- function(
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    outcome_nodes = 80L,
    u_nodes = 40L,
    x_batch_size = 25L
) {
  outcome_nodes <- as.integer(outcome_nodes)
  u_nodes <- as.integer(u_nodes)
  x_batch_size <- as.integer(x_batch_size)
  if (outcome_nodes < 8L || u_nodes < 10L || x_batch_size < 1L) {
    stop("Invalid DGP oracle gamma-normalizer integration settings.")
  }

  function(fold, treatment, gamma_fit, x_cols) {
    if (length(x_cols) != 1L) {
      stop("The DGP oracle gamma normalizer currently supports scalar X.")
    }
    local_treatment <- treatment
    local_x_col <- x_cols[1L]
    local_gamma_fit <- gamma_fit
    normalizer_cache <- new.env(hash = TRUE, parent = emptyenv())
    density_error_cache <- new.env(hash = TRUE, parent = emptyenv())

    key_x <- function(x) sprintf("%.17g", as.numeric(x))

    add_missing_normalizers <- function(x) {
      unique_x <- unique(as.numeric(x))
      keys <- key_x(unique_x)
      missing <- !vapply(
        keys,
        exists,
        logical(1),
        envir = normalizer_cache,
        inherits = FALSE
      )
      unique_x <- unique_x[missing]
      keys <- keys[missing]
      if (length(unique_x) == 0L) return(invisible(NULL))

      for (start in seq.int(1L, length(unique_x), by = x_batch_size)) {
        finish <- min(length(unique_x), start + x_batch_size - 1L)
        index <- start:finish
        local_x <- unique_x[index]
        local_keys <- keys[index]
        local_validation <- data.frame(X = local_x)
        names(local_validation) <- local_x_col
        conditional_grid <- audit_build_conditional_outcome_grid(
          validation = local_validation,
          t_value = local_treatment,
          x_col = local_x_col,
          confounder_strength = confounder_strength,
          treatment_intercept = treatment_intercept,
          instrument_strength = instrument_strength,
          treatment_effect = treatment_effect,
          full_interval_nodes = outcome_nodes,
          u_nodes = u_nodes,
          batch_size = x_batch_size
        )
        center <- treatment_effect * local_treatment + 1.5 * local_x
        grid_data <- data.frame(
          Y = rep(center, each = outcome_nodes) +
            as.vector(t(conditional_grid$relative_y))
        )
        grid_data[[local_x_col]] <- rep(local_x, each = outcome_nodes)
        raw_gamma <- local_gamma_fit$predict(grid_data)$normalized
        raw_gamma <- matrix(raw_gamma, nrow = length(local_x), byrow = TRUE)
        normalizer <- rowSums(conditional_grid$probability0 * raw_gamma)
        if (any(!is.finite(normalizer)) || any(normalizer <= 0)) {
          stop("The DGP oracle conditional gamma normalizer is invalid.")
        }
        density_error <- rowSums(abs(
          conditional_grid$probability0 * raw_gamma -
            conditional_grid$probability1
        ))
        for (local_index in seq_along(local_keys)) {
          assign(
            local_keys[local_index],
            normalizer[local_index],
            envir = normalizer_cache
          )
          assign(
            local_keys[local_index],
            density_error[local_index],
            envir = density_error_cache
          )
        }
      }
      invisible(NULL)
    }

    transform <- function(newdata, normalized_predictions) {
      if (!local_x_col %in% names(newdata)) {
        stop("Gamma transformation data are missing the X column.")
      }
      x <- as.numeric(newdata[[local_x_col]])
      add_missing_normalizers(x)
      normalizer <- vapply(
        key_x(x),
        get,
        numeric(1),
        envir = normalizer_cache,
        inherits = FALSE
      )
      as.numeric(normalized_predictions) / normalizer
    }

    diagnostics <- function() {
      keys <- ls(normalizer_cache, all.names = TRUE)
      if (length(keys) == 0L) return(data.frame())
      normalizer <- vapply(
        keys, get, numeric(1), envir = normalizer_cache, inherits = FALSE
      )
      density_error <- vapply(
        keys, get, numeric(1), envir = density_error_cache, inherits = FALSE
      )
      data.frame(
        gamma_transform = "dgp_oracle_conditional_normalization",
        outcome_nodes = outcome_nodes,
        u_nodes = u_nodes,
        unique_x_count = length(keys),
        mean_raw_conditional_normalization_error = mean(normalizer - 1),
        mean_abs_raw_conditional_normalization_error =
          mean(abs(normalizer - 1)),
        max_abs_raw_conditional_normalization_error =
          max(abs(normalizer - 1)),
        min_raw_conditional_normalizer = min(normalizer),
        max_raw_conditional_normalizer = max(normalizer),
        mean_raw_gamma_density_l1_error = mean(density_error),
        max_raw_gamma_density_l1_error = max(density_error),
        stringsAsFactors = FALSE
      )
    }

    list(transform = transform, diagnostics = diagnostics)
  }
}

audit_integrate_eif_means_with_gamma <- function(
    conditional_grid,
    validation,
    nuisance_validation,
    t_value,
    Lambda0,
    Lambda1,
    q_lower_validation,
    q_upper_validation,
    gamma_source = c("oracle", "fitted"),
    gamma_predictor = NULL,
    compatibility_action = c("stop", "project"),
    gamma_projection_margin = 1e-4,
    batch_size = 25L
) {
  gamma_source <- match.arg(gamma_source)
  compatibility_action <- match.arg(compatibility_action)
  if (gamma_source == "fitted" && !is.function(gamma_predictor)) {
    stop("A gamma_predictor is required for fitted-gamma integration.")
  }
  n <- nrow(validation)
  y_nodes <- ncol(conditional_grid$relative_y)
  if (
    length(conditional_grid$x) != n ||
      length(q_lower_validation) != n ||
      length(q_upper_validation) != n
  ) {
    stop("Conditional-grid and EIF validation inputs have inconsistent size.")
  }

  baseline_names <- c(
    ".m0",
    unlist(lapply(c("lower", "upper"), function(endpoint) {
      paste0(
        c(".theta_", ".D_r0_", ".D_r1_", ".H0_"),
        endpoint
      )
    }))
  )
  comparison_names <- paste0(".H1_", c("lower", "upper"))
  baseline_prediction <- matrix(
    NA_real_, nrow = n, ncol = length(baseline_names),
    dimnames = list(NULL, baseline_names)
  )
  comparison_prediction <- matrix(
    NA_real_, nrow = n, ncol = length(comparison_names),
    dimnames = list(NULL, comparison_names)
  )
  gamma_normalization_error <- gamma_density_l1_error <-
    projected_gamma_normalization_error <- rep(NA_real_, n)
  lower_lambda_normalization_error <-
    upper_lambda_normalization_error <- rep(NA_real_, n)
  gamma_projection_mass0 <- rep(NA_real_, n)

  x_col <- conditional_grid$settings$x_col
  treatment_effect <- conditional_grid$settings$treatment_effect
  for (start in seq.int(1L, n, by = as.integer(batch_size))) {
    finish <- min(n, start + as.integer(batch_size) - 1L)
    index <- start:finish
    local_n <- length(index)
    probability0 <- conditional_grid$probability0[index, , drop = FALSE]
    probability1 <- conditional_grid$probability1[index, , drop = FALSE]
    center <- treatment_effect * t_value +
      1.5 * conditional_grid$x[index]
    y <- rep(center, each = y_nodes) +
      as.vector(t(
        conditional_grid$relative_y[index, , drop = FALSE]
      ))
    x <- rep(conditional_grid$x[index], each = y_nodes)
    probability0_vector <- as.vector(t(probability0))
    probability1_vector <- as.vector(t(probability1))
    relevant <- probability0_vector > 1e-14

    gamma <- rep(1, length(y))
    if (gamma_source == "oracle") {
      gamma[relevant] <-
        probability1_vector[relevant] / probability0_vector[relevant]
    } else {
      prediction_data <- data.frame(Y = y)
      prediction_data[[x_col]] <- x
      gamma <- as.numeric(gamma_predictor(prediction_data))
    }
    if (any(!is.finite(gamma)) || any(gamma <= 0)) {
      stop("Integrated gamma predictions must be finite and positive.")
    }

    nuisance_grid <- lapply(nuisance_validation, function(value) {
      rep(value[index], each = y_nodes)
    })
    envelope <- compute_projected_envelopes(
      nuisance = nuisance_grid,
      gamma = gamma,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      gamma_compatibility_action = compatibility_action,
      gamma_projection_margin = gamma_projection_margin
    )
    primitive <- compute_eif_row_primitives(
      data = data.frame(Y = y),
      nuisance = nuisance_grid,
      envelope = envelope,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      q_lower = rep(q_lower_validation[index], each = y_nodes),
      q_upper = rep(q_upper_validation[index], each = y_nodes)
    )
    matrix_by_x <- function(value) {
      matrix(value, nrow = local_n, byrow = TRUE)
    }
    mean0 <- function(value) {
      rowSums(probability0 * matrix_by_x(value))
    }
    mean1 <- function(value) {
      rowSums(probability1 * matrix_by_x(value))
    }

    baseline_prediction[index, ".m0"] <-
      rowSums(probability0 * matrix(y, nrow = local_n, byrow = TRUE))
    for (endpoint in c("lower", "upper")) {
      local <- primitive[[endpoint]]
      baseline_prediction[index, paste0(".theta_", endpoint)] <-
        mean0(local$lambda * y)
      baseline_prediction[index, paste0(".D_r0_", endpoint)] <-
        mean0(local$D_r0)
      baseline_prediction[index, paste0(".D_r1_", endpoint)] <-
        mean0(local$D_r1)
      baseline_prediction[index, paste0(".H0_", endpoint)] <-
        mean0(local$H0)
      comparison_prediction[index, paste0(".H1_", endpoint)] <-
        mean1(local$H1)
    }

    gamma_matrix <- matrix_by_x(gamma)
    projected_gamma_matrix <- matrix_by_x(envelope$gamma_used)
    gamma_normalization_error[index] <-
      rowSums(probability0 * gamma_matrix) - 1
    gamma_density_l1_error[index] <- rowSums(abs(
      probability0 * gamma_matrix - probability1
    ))
    projected_gamma_normalization_error[index] <-
      rowSums(probability0 * projected_gamma_matrix) - 1
    gamma_projection_mass0[index] <- rowSums(
      probability0 * matrix_by_x(as.numeric(envelope$gamma_projected))
    )
    lower_lambda_normalization_error[index] <-
      mean0(primitive$lower$lambda) - 1
    upper_lambda_normalization_error[index] <-
      mean0(primitive$upper$lambda) - 1
  }

  summarize_error <- function(value) {
    c(mean_absolute = mean(abs(value)), maximum_absolute = max(abs(value)))
  }
  gamma_normalization <- summarize_error(gamma_normalization_error)
  projected_gamma_normalization <- summarize_error(
    projected_gamma_normalization_error
  )
  lower_lambda_normalization <- summarize_error(
    lower_lambda_normalization_error
  )
  upper_lambda_normalization <- summarize_error(
    upper_lambda_normalization_error
  )

  list(
    baseline_prediction = baseline_prediction,
    comparison_prediction = comparison_prediction,
    learner = paste0("dgp_integrated_", gamma_source, "_gamma"),
    diagnostics = data.frame(
      oracle_integration_y_nodes = y_nodes,
      oracle_integration_interval_nodes =
        conditional_grid$settings$interval_nodes,
      oracle_integration_quadrature_type =
        conditional_grid$settings$quadrature_type,
      oracle_integration_u_nodes = conditional_grid$settings$u_nodes,
      oracle_integration_max_density_mass_error_z0 = max(abs(
        conditional_grid$raw_mass0 - 1
      )),
      oracle_integration_max_density_mass_error_z1 = max(abs(
        conditional_grid$raw_mass1 - 1
      )),
      oracle_integration_mean_abs_gamma_normalization_error =
        unname(gamma_normalization["mean_absolute"]),
      oracle_integration_max_abs_gamma_normalization_error =
        unname(gamma_normalization["maximum_absolute"]),
      oracle_integration_mean_gamma_density_l1_error =
        mean(gamma_density_l1_error),
      oracle_integration_max_gamma_density_l1_error =
        max(gamma_density_l1_error),
      oracle_integration_mean_abs_projected_gamma_normalization_error =
        unname(projected_gamma_normalization["mean_absolute"]),
      oracle_integration_max_abs_projected_gamma_normalization_error =
        unname(projected_gamma_normalization["maximum_absolute"]),
      oracle_integration_mean_gamma_projection_mass0 =
        mean(gamma_projection_mass0),
      oracle_integration_max_gamma_projection_mass0 =
        max(gamma_projection_mass0),
      oracle_integration_mean_abs_lower_lambda_normalization_error =
        unname(lower_lambda_normalization["mean_absolute"]),
      oracle_integration_max_abs_lower_lambda_normalization_error =
        unname(lower_lambda_normalization["maximum_absolute"]),
      oracle_integration_mean_abs_upper_lambda_normalization_error =
        unname(upper_lambda_normalization["mean_absolute"]),
      oracle_integration_max_abs_upper_lambda_normalization_error =
        unname(upper_lambda_normalization["maximum_absolute"]),
      stringsAsFactors = FALSE
    )
  )
}

audit_gamma_endpoint_moment_template <- function() {
  output <- list(
    gamma_moment_mean_normalization_error = NA_real_,
    gamma_moment_mean_abs_error = NA_real_,
    gamma_moment_mean_signed_outcome_error = NA_real_,
    gamma_moment_mean_abs_outcome_error = NA_real_,
    gamma_moment_active_lower_disagreement = NA_real_,
    gamma_moment_active_upper_disagreement = NA_real_,
    gamma_moment_fitted_projection_mass = NA_real_,
    gamma_moment_oracle_projection_mass = NA_real_
  )
  for (direction in c("lower", "upper")) {
    for (name in c(
      "theta_exact_gap",
      "theta_linearized_gap",
      "theta_nonlinear_remainder",
      "mean_abs_linearized_theta_contribution",
      "lambda_normalization_gap",
      "D_r0_gap",
      "D_r1_gap",
      "H0_gap",
      "H1_gap",
      "linked_M_gap",
      "target_M_exact_gap",
      "target_M_linearized_gap",
      "linked_minus_target_M_gap",
      "C0_gap",
      "C1_gap"
    )) {
      output[[paste0("gamma_moment_", direction, "_", name)]] <- NA_real_
    }
  }
  as.data.frame(output, stringsAsFactors = FALSE)
}

audit_compare_fitted_oracle_gamma_moments <- function(
    conditional_grid,
    validation,
    nuisance_validation,
    t_value,
    Lambda0,
    Lambda1,
    q_lower_validation,
    q_upper_validation,
    gamma_predictor,
    gamma_projection_margin = 1e-4,
    batch_size = 25L
) {
  if (!is.function(gamma_predictor)) {
    stop("gamma_predictor must be a function for the gamma-moment audit.")
  }
  n <- nrow(validation)
  y_nodes <- ncol(conditional_grid$relative_y)
  if (
    length(conditional_grid$x) != n ||
      length(q_lower_validation) != n ||
      length(q_upper_validation) != n
  ) {
    stop("Gamma-moment audit inputs have inconsistent dimensions.")
  }

  storage <- lapply(
    names(audit_gamma_endpoint_moment_template()),
    function(name) rep(NA_real_, n)
  )
  names(storage) <- names(audit_gamma_endpoint_moment_template())
  x_col <- conditional_grid$settings$x_col
  treatment_effect <- conditional_grid$settings$treatment_effect

  for (start in seq.int(1L, n, by = as.integer(batch_size))) {
    finish <- min(n, start + as.integer(batch_size) - 1L)
    index <- start:finish
    local_n <- length(index)
    probability0 <- conditional_grid$probability0[index, , drop = FALSE]
    probability1 <- conditional_grid$probability1[index, , drop = FALSE]
    center <- treatment_effect * t_value +
      1.5 * conditional_grid$x[index]
    y <- rep(center, each = y_nodes) + as.vector(t(
      conditional_grid$relative_y[index, , drop = FALSE]
    ))
    x <- rep(conditional_grid$x[index], each = y_nodes)
    prediction_data <- data.frame(Y = y)
    prediction_data[[x_col]] <- x
    gamma_fitted <- as.numeric(gamma_predictor(prediction_data))
    gamma_oracle <- as.vector(t(probability1)) /
      pmax(as.vector(t(probability0)), 1e-300)
    if (
      any(!is.finite(gamma_fitted)) || any(gamma_fitted <= 0) ||
        any(!is.finite(gamma_oracle)) || any(gamma_oracle <= 0)
    ) {
      stop("Gamma-moment audit encountered invalid gamma values.")
    }

    nuisance_grid <- lapply(nuisance_validation, function(value) {
      rep(value[index], each = y_nodes)
    })
    fitted_envelope <- compute_projected_envelopes(
      nuisance = nuisance_grid,
      gamma = gamma_fitted,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      gamma_compatibility_action = "project",
      gamma_projection_margin = gamma_projection_margin
    )
    oracle_envelope <- compute_projected_envelopes(
      nuisance = nuisance_grid,
      gamma = gamma_oracle,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      gamma_compatibility_action = "project",
      gamma_projection_margin = gamma_projection_margin
    )
    fitted_primitive <- compute_eif_row_primitives(
      data = data.frame(Y = y),
      nuisance = nuisance_grid,
      envelope = fitted_envelope,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      q_lower = rep(q_lower_validation[index], each = y_nodes),
      q_upper = rep(q_upper_validation[index], each = y_nodes)
    )
    oracle_primitive <- compute_eif_row_primitives(
      data = data.frame(Y = y),
      nuisance = nuisance_grid,
      envelope = oracle_envelope,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      q_lower = rep(q_lower_validation[index], each = y_nodes),
      q_upper = rep(q_upper_validation[index], each = y_nodes)
    )
    fitted_derivative <- compute_eif_envelope_derivatives(
      nuisance = nuisance_grid,
      gamma = fitted_envelope$gamma_used,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1
    )
    oracle_derivative <- compute_eif_envelope_derivatives(
      nuisance = nuisance_grid,
      gamma = oracle_envelope$gamma_used,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1
    )

    matrix_by_x <- function(value) matrix(value, nrow = local_n, byrow = TRUE)
    mean0 <- function(value) rowSums(probability0 * matrix_by_x(value))
    mean1 <- function(value) rowSums(probability1 * matrix_by_x(value))
    delta_raw <- gamma_fitted - gamma_oracle
    delta_used <- fitted_envelope$gamma_used - oracle_envelope$gamma_used

    storage$gamma_moment_mean_normalization_error[index] <- mean0(delta_raw)
    storage$gamma_moment_mean_abs_error[index] <- mean0(abs(delta_raw))
    storage$gamma_moment_mean_signed_outcome_error[index] <-
      mean0(y * delta_raw)
    storage$gamma_moment_mean_abs_outcome_error[index] <-
      mean0(abs(y * delta_raw))
    storage$gamma_moment_active_lower_disagreement[index] <- mean0(
      as.numeric(
        fitted_derivative$active_lower != oracle_derivative$active_lower
      )
    )
    storage$gamma_moment_active_upper_disagreement[index] <- mean0(
      as.numeric(
        fitted_derivative$active_upper != oracle_derivative$active_upper
      )
    )
    storage$gamma_moment_fitted_projection_mass[index] <- mean0(
      as.numeric(fitted_envelope$gamma_projected)
    )
    storage$gamma_moment_oracle_projection_mass[index] <- mean0(
      as.numeric(oracle_envelope$gamma_projected)
    )

    for (direction in c("lower", "upper")) {
      fitted <- fitted_primitive[[direction]]
      oracle <- oracle_primitive[[direction]]
      q <- if (direction == "lower") {
        rep(q_lower_validation[index], each = y_nodes)
      } else {
        rep(q_upper_validation[index], each = y_nodes)
      }
      allocate_upper <- if (direction == "lower") y <= q else y > q
      allocate_upper <- as.numeric(allocate_upper)
      allocate_lower <- 1 - allocate_upper
      derivative_lambda <-
        allocate_upper * oracle_derivative$upper$gamma +
        allocate_lower * oracle_derivative$lower$gamma
      exact_theta_gap <- mean0(y * (fitted$lambda - oracle$lambda))
      linearized_theta_gap <- mean0(y * derivative_lambda * delta_used)
      prefix <- paste0("gamma_moment_", direction, "_")
      storage[[paste0(prefix, "theta_exact_gap")]][index] <- exact_theta_gap
      storage[[paste0(prefix, "theta_linearized_gap")]][index] <-
        linearized_theta_gap
      storage[[paste0(prefix, "theta_nonlinear_remainder")]][index] <-
        exact_theta_gap - linearized_theta_gap
      storage[[paste0(
        prefix, "mean_abs_linearized_theta_contribution"
      )]][index] <- mean0(abs(y * derivative_lambda * delta_used))
      storage[[paste0(prefix, "lambda_normalization_gap")]][index] <-
        mean0(fitted$lambda - oracle$lambda)
      storage[[paste0(prefix, "D_r0_gap")]][index] <-
        mean0(fitted$D_r0 - oracle$D_r0)
      storage[[paste0(prefix, "D_r1_gap")]][index] <-
        mean0(fitted$D_r1 - oracle$D_r1)
      H0_gap <- mean0(fitted$H0 - oracle$H0)
      H1_gap <- mean1(fitted$H1 - oracle$H1)
      storage[[paste0(prefix, "H0_gap")]][index] <- H0_gap
      storage[[paste0(prefix, "H1_gap")]][index] <- H1_gap
      r0 <- oracle_primitive$r0
      linked_M_gap <- H0_gap + H1_gap
      target_M_exact_gap <- mean0(
        (1 - r0) * y * (fitted$lambda - oracle$lambda)
      )
      target_M_linearized_gap <- mean0(
        (1 - r0) * y * derivative_lambda * delta_used
      )
      storage[[paste0(prefix, "linked_M_gap")]][index] <- linked_M_gap
      storage[[paste0(prefix, "target_M_exact_gap")]][index] <-
        target_M_exact_gap
      storage[[paste0(prefix, "target_M_linearized_gap")]][index] <-
        target_M_linearized_gap
      storage[[paste0(prefix, "linked_minus_target_M_gap")]][index] <-
        linked_M_gap - target_M_exact_gap
      storage[[paste0(prefix, "C0_gap")]][index] <-
        -exact_theta_gap +
        mean0((1 - r0) * (fitted$D_r0 - oracle$D_r0))
      storage[[paste0(prefix, "C1_gap")]][index] <-
        mean0((1 - r0) * (fitted$D_r1 - oracle$D_r1))
    }
  }

  output <- lapply(storage, mean)
  as.data.frame(output, stringsAsFactors = FALSE)
}

make_dgp_integrated_eif_provider_factory <- function(
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    y_nodes = 401L,
    interval_nodes = 40L,
    u_nodes = 40L,
    batch_size = 25L,
    gamma_projection_margin = 1e-4,
    gamma_moment_audit = FALSE
) {
  function(mode, fold, treatment, gamma_fit, x_cols) {
    if (!mode %in% c("oracle_all", "oracle_gamma", "rkhs")) {
      return(NULL)
    }
    if (length(x_cols) != 1L) {
      stop("The DGP-integrated EIF provider currently supports scalar X.")
    }
    local_mode <- mode
    local_fold <- fold
    local_treatment <- treatment
    local_gamma_fit <- gamma_fit
    local_x_col <- x_cols[1L]
    function(
        validation,
        nuisance_validation,
        envelope_validation,
        t_value,
        Lambda0,
        Lambda1,
        q_lower_validation,
        q_upper_validation,
        baseline_outcomes,
        comparison_outcomes
    ) {
      if (t_value != local_treatment) {
        stop("The DGP-integrated provider received the wrong treatment arm.")
      }
      conditional_grid <- audit_build_conditional_outcome_grid(
        validation = validation,
        t_value = local_treatment,
        x_col = local_x_col,
        confounder_strength = confounder_strength,
        treatment_intercept = treatment_intercept,
        instrument_strength = instrument_strength,
        treatment_effect = treatment_effect,
        y_nodes = y_nodes,
        q_lower = q_lower_validation,
        q_upper = q_upper_validation,
        interval_nodes = interval_nodes,
        u_nodes = u_nodes,
        batch_size = batch_size
      )
      gamma_source <- if (local_mode == "rkhs") "fitted" else "oracle"
      gamma_predictor <- if (gamma_source == "fitted") {
        function(newdata) {
          local_gamma_fit$predict(newdata)$normalized
        }
      } else {
        NULL
      }
      output <- audit_integrate_eif_means_with_gamma(
        conditional_grid = conditional_grid,
        validation = validation,
        nuisance_validation = nuisance_validation,
        t_value = local_treatment,
        Lambda0 = Lambda0,
        Lambda1 = Lambda1,
        q_lower_validation = q_lower_validation,
        q_upper_validation = q_upper_validation,
        gamma_source = gamma_source,
        gamma_predictor = gamma_predictor,
        compatibility_action = if (local_mode == "oracle_all") {
          "stop"
        } else {
          "project"
        },
        gamma_projection_margin = gamma_projection_margin,
        batch_size = batch_size
      )
      if (isTRUE(gamma_moment_audit)) {
        moment_diagnostics <- audit_gamma_endpoint_moment_template()
        if (local_mode == "rkhs") {
          moment_diagnostics <- audit_compare_fitted_oracle_gamma_moments(
            conditional_grid = conditional_grid,
            validation = validation,
            nuisance_validation = nuisance_validation,
            t_value = local_treatment,
            Lambda0 = Lambda0,
            Lambda1 = Lambda1,
            q_lower_validation = q_lower_validation,
            q_upper_validation = q_upper_validation,
            gamma_predictor = gamma_predictor,
            gamma_projection_margin = gamma_projection_margin,
            batch_size = batch_size
          )
        }
        output$diagnostics <- cbind(output$diagnostics, moment_diagnostics)
      }
      if (
        !identical(colnames(output$baseline_prediction), baseline_outcomes) ||
          !identical(
            colnames(output$comparison_prediction), comparison_outcomes
          )
      ) {
        stop("The DGP-integrated provider returned unexpected columns.")
      }
      output
    }
  }
}

audit_lower_tail_allocation <- function(y, mass, target_mass) {
  target_mass <- min(sum(mass), max(0, target_mass))
  allocation <- numeric(length(y))

  if (target_mass <= 0) {
    return(list(
      allocation = allocation,
      threshold = min(y),
      achieved_mass = 0
    ))
  }

  cumulative <- cumsum(mass)
  boundary <- which(cumulative >= target_mass)[1L]
  before <- if (boundary == 1L) 0 else cumulative[boundary - 1L]
  remaining <- target_mass - before
  fraction <- if (mass[boundary] > 0) {
    remaining / mass[boundary]
  } else {
    0
  }

  if (boundary > 1L) {
    allocation[seq_len(boundary - 1L)] <- 1
  }
  allocation[boundary] <- fraction

  list(
    allocation = allocation,
    threshold = y[boundary],
    achieved_mass = sum(allocation * mass)
  )
}

compute_dgp_oracle_eif_at_x <- function(
    x,
    t_value,
    u,
    u_weights,
    Lambda0 = 5,
    Lambda1 = 5,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    envelope_smoothing_epsilon = 0,
    y_nodes = 4001L,
    return_threshold_distribution = FALSE
) {
  eta0 <- treatment_intercept + x + confounder_strength * u
  eta1 <- eta0 + instrument_strength
  p_t1_z0 <- stats::plogis(eta0)
  p_t1_z1 <- stats::plogis(eta1)
  pi0 <- sum(u_weights * p_t1_z0)
  pi1 <- sum(u_weights * p_t1_z1)

  treatment_weight0 <- if (t_value == 1) {
    p_t1_z0
  } else {
    1 - p_t1_z0
  }
  treatment_weight1 <- if (t_value == 1) {
    p_t1_z1
  } else {
    1 - p_t1_z1
  }
  r0 <- if (t_value == 1) pi0 else 1 - pi0
  r1 <- if (t_value == 1) pi1 else 1 - pi1

  outcome_center <- treatment_effect * t_value + 1.5 * x
  y_rule <- audit_trapezoid_rule(
    outcome_center - 5,
    outcome_center + 5,
    y_nodes
  )
  y <- y_rule$nodes
  y_weights <- y_rule$weights
  outcome_mean <- outcome_center + 1.5 * u

  density_matrix <- outer(
    y,
    outcome_mean,
    function(y_value, mean_value) {
      dtnormtc(y_value, mu = mean_value, s = 1)
    }
  )

  joint0 <- as.vector(
    density_matrix %*% (u_weights * treatment_weight0)
  )
  joint1 <- as.vector(
    density_matrix %*% (u_weights * treatment_weight1)
  )
  f0 <- joint0 / r0
  f1 <- joint1 / r1
  f0 <- f0 / sum(y_weights * f0)
  f1 <- f1 / sum(y_weights * f1)

  relevant <- y_weights * f0 > 1e-13
  gamma <- rep(1, length(y))
  gamma[relevant] <- f1[relevant] / f0[relevant]

  nuisance <- list(
    e = rep(stats::plogis(-0.5 * x), length(y)),
    pi0 = rep(pi0, length(y)),
    pi1 = rep(pi1, length(y))
  )
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = gamma,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    gamma_compatibility_action = "stop"
  )

  if (any(envelope$L[relevant] > envelope$U[relevant] + 1e-8)) {
    stop(sprintf(
      "Oracle envelopes cross at x=%.6f, treatment=%d.",
      x,
      t_value
    ))
  }

  probability_mass0 <- y_weights * f0
  probability_mass1 <- y_weights * f1
  Lbar <- sum(probability_mass0 * envelope$L)
  Ubar <- sum(probability_mass0 * envelope$U)
  gap_mass <- probability_mass0 * pmax(envelope$D, 0)
  derivative <- compute_eif_envelope_derivatives(
    nuisance = nuisance,
    gamma = gamma,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    envelope_smoothing_epsilon = envelope_smoothing_epsilon
  )

  make_endpoint <- function(direction = c("lower", "upper")) {
    direction <- match.arg(direction)
    target_mass <- if (direction == "lower") {
      1 - Lbar
    } else {
      Ubar - 1
    }
    tail <- audit_lower_tail_allocation(y, gap_mass, target_mass)
    allocate_upper <- if (direction == "lower") {
      tail$allocation
    } else {
      1 - tail$allocation
    }
    allocate_lower <- 1 - allocate_upper
    lambda <- allocate_upper * envelope$U +
      allocate_lower * envelope$L

    envelope_derivative <- function(name) {
      (y - tail$threshold) * (
        allocate_upper * derivative$upper[[name]] +
          allocate_lower * derivative$lower[[name]]
      )
    }

    D_gamma <- envelope_derivative("gamma")
    D_r0 <- envelope_derivative("r0")
    D_r1 <- envelope_derivative("r1")
    H0 <- r0 * y + (1 - r0) * (
      lambda * (y - tail$threshold) - gamma * D_gamma
    )
    H1 <- (1 - r0) * D_gamma

    m0 <- sum(probability_mass0 * y)
    theta <- sum(probability_mass0 * lambda * y)
    mean_D_r0 <- sum(probability_mass0 * D_r0)
    mean_D_r1 <- sum(probability_mass0 * D_r1)

    c(
      q = tail$threshold,
      M = r0 * m0 + (1 - r0) * theta,
      C0 = m0 - theta + (1 - r0) * mean_D_r0,
      C1 = (1 - r0) * mean_D_r1,
      Hbar0 = sum(probability_mass0 * H0),
      Hbar1 = sum(probability_mass1 * H1),
      m0 = m0,
      theta = theta,
      lambda_normalization = sum(probability_mass0 * lambda),
      tail_mass_error = tail$achieved_mass - target_mass
    )
  }

  lower <- make_endpoint("lower")
  upper <- make_endpoint("upper")

  output <- data.frame(
    x = x,
    treatment = t_value,
    pi0 = pi0,
    pi1 = pi1,
    Lbar = Lbar,
    Ubar = Ubar,
    stringsAsFactors = FALSE
  )

  for (name in names(lower)) {
    output[[paste0("lower_", name)]] <- lower[[name]]
    output[[paste0("upper_", name)]] <- upper[[name]]
  }

  if (return_threshold_distribution) {
    return(list(
      summary = output,
      distribution = data.frame(
        y = y,
        probability_mass0 = probability_mass0,
        L = envelope$L,
        U = envelope$U,
        D = pmax(envelope$D, 0),
        stringsAsFactors = FALSE
      ),
      r0 = r0,
      r1 = r1
    ))
  }

  output
}

build_dgp_oracle_eif_lookup <- function(
    x_nodes = 401L,
    u_nodes = 80L,
    y_nodes = 4001L,
    Lambda0 = 5,
    Lambda1 = 5,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    envelope_smoothing_epsilon = 0,
    verbose = TRUE
) {
  x_grid <- seq(-2, 2, length.out = x_nodes)
  u_rule <- gauss_legendre_rule(u_nodes, -2, 2)
  u_weights <- u_rule$weights * d_u(u_rule$nodes)
  u_weights <- u_weights / sum(u_weights)
  rows <- vector("list", 2L * x_nodes)
  counter <- 1L

  for (t_value in c(0, 1)) {
    if (verbose) {
      message(sprintf(
        "Building DGP-oracle EIF lookup for treatment %d (%d X nodes).",
        t_value,
        x_nodes
      ))
    }

    for (x in x_grid) {
      rows[[counter]] <- compute_dgp_oracle_eif_at_x(
        x = x,
        t_value = t_value,
        u = u_rule$nodes,
        u_weights = u_weights,
        Lambda0 = Lambda0,
        Lambda1 = Lambda1,
        confounder_strength = confounder_strength,
        treatment_intercept = treatment_intercept,
        instrument_strength = instrument_strength,
        treatment_effect = treatment_effect,
        envelope_smoothing_epsilon = envelope_smoothing_epsilon,
        y_nodes = y_nodes
      )
      counter <- counter + 1L
    }
  }

  list(
    grid = do.call(rbind, rows),
    settings = list(
      x_nodes = x_nodes,
      u_nodes = u_nodes,
      y_nodes = y_nodes,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      confounder_strength = confounder_strength,
      treatment_intercept = treatment_intercept,
      instrument_strength = instrument_strength,
      treatment_effect = treatment_effect,
      envelope_smoothing_epsilon = envelope_smoothing_epsilon
    )
  )
}

audit_population_spline_spec <- function(
    t_value,
    u,
    u_weights,
    spline_df = 5L,
    confounder_strength = 1,
    treatment_intercept = 0,
    dense_x_nodes = 20001L,
    x_distribution_rule = NULL
) {
  spline_df <- as.integer(spline_df)
  dense_x_nodes <- as.integer(dense_x_nodes)
  if (spline_df < 1L) stop("spline_df must be positive.")
  if (dense_x_nodes < 101L) {
    stop("dense_x_nodes must be at least 101.")
  }

  if (is.null(x_distribution_rule)) {
    x <- seq(-2, 2, length.out = dense_x_nodes)
    marginal_weights <- NULL
  } else {
    if (!is.list(x_distribution_rule) ||
        !all(c("nodes", "weights") %in% names(x_distribution_rule))) {
      stop("x_distribution_rule must contain nodes and weights.")
    }
    x <- as.numeric(x_distribution_rule$nodes)
    marginal_weights <- as.numeric(x_distribution_rule$weights)
    if (length(x) < 2L || length(marginal_weights) != length(x) ||
        any(!is.finite(x)) || any(!is.finite(marginal_weights)) ||
        any(marginal_weights < 0) || sum(marginal_weights) <= 0) {
      stop("The supplied X distribution rule is invalid.")
    }
    ordering <- order(x)
    x <- x[ordering]
    marginal_weights <- marginal_weights[ordering] / sum(marginal_weights)
  }
  eta0 <- outer(
    x,
    u,
    function(x_value, u_value) {
      treatment_intercept + x_value +
        confounder_strength * u_value
    }
  )
  pi0 <- as.vector(stats::plogis(eta0) %*% u_weights)
  r0 <- if (t_value == 1) pi0 else 1 - pi0
  instrument0_probability <- stats::plogis(-0.5 * x)
  if (is.null(marginal_weights)) {
    baseline_density <- dtnormtc(x) * instrument0_probability * r0
    step <- x[2L] - x[1L]
    increments <- step * (
      baseline_density[-length(x)] + baseline_density[-1L]
    ) / 2
    cumulative <- c(0, cumsum(increments))
    baseline_probability <- cumulative[length(cumulative)]
    quantile_cumulative <- cumulative / baseline_probability
    quantile_x <- x
  } else {
    baseline_mass <- marginal_weights * instrument0_probability * r0
    baseline_probability <- sum(baseline_mass)
    quantile_cumulative <- c(0, cumsum(baseline_mass) / baseline_probability)
    quantile_x <- c(x[1L], x)
  }
  if (!is.finite(baseline_probability) || baseline_probability <= 0) {
    stop("The population baseline-stratum probability is invalid.")
  }

  knot_count <- spline_df - 1L
  knot_probabilities <- if (knot_count > 0L) {
    seq_len(knot_count) / (knot_count + 1L)
  } else {
    numeric()
  }
  knots <- if (knot_count > 0L) {
    stats::approx(
      x = quantile_cumulative,
      y = quantile_x,
      xout = knot_probabilities,
      method = "linear",
      ties = "ordered"
    )$y
  } else {
    numeric()
  }

  list(
    knots = knots,
    knot_probabilities = knot_probabilities,
    boundary_knots = range(x),
    spline_df = spline_df,
    baseline_probability = baseline_probability
  )
}

audit_population_spline_matrix <- function(x, spline_spec) {
  basis <- splines::ns(
    x,
    knots = spline_spec$knots,
    Boundary.knots = spline_spec$boundary_knots,
    intercept = FALSE
  )
  colnames(basis) <- paste0("X_s", seq_len(ncol(basis)))
  output <- cbind(intercept = 1, basis)
  storage.mode(output) <- "double"
  output
}

audit_invert_linear_cumulative <- function(
    cumulative,
    density,
    target,
    grid,
    outcome_center
) {
  row_count <- nrow(cumulative)
  output <- numeric(row_count)
  for (row in seq_len(row_count)) {
    local_target <- min(
      cumulative[row, length(grid)],
      max(0, target[row])
    )
    if (local_target <= 0) {
      output[row] <- outcome_center[row] + grid[1L]
      next
    }
    boundary <- which(cumulative[row, ] >= local_target)[1L]
    if (boundary <= 1L) {
      output[row] <- outcome_center[row] + grid[1L]
      next
    }
    left <- boundary - 1L
    interval_width <- grid[boundary] - grid[left]
    remaining <- local_target - cumulative[row, left]
    g0 <- density[row, left]
    slope <- (density[row, boundary] - g0) / interval_width
    root <- stats::uniroot(
      function(distance) {
        g0 * distance + slope * distance^2 / 2 - remaining
      },
      interval = c(0, interval_width),
      tol = 1e-12
    )$root
    output[row] <- outcome_center[row] + grid[left] + root
  }
  output
}

build_population_restricted_threshold_problem <- function(
    t_value,
    x_nodes = 80L,
    u_nodes = 80L,
    y_nodes = 4001L,
    spline_df = 5L,
    Lambda0 = 5,
    Lambda1 = 5,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    envelope_smoothing_epsilon = 0,
    x_distribution_rule = NULL,
    verbose = TRUE
) {
  if (!t_value %in% c(0, 1)) stop("t_value must be zero or one.")
  x_rule <- if (is.null(x_distribution_rule)) {
    local <- gauss_legendre_rule(x_nodes, -2, 2)
    local$weights <- local$weights * dtnormtc(local$nodes)
    local
  } else {
    if (!is.list(x_distribution_rule) ||
        !all(c("nodes", "weights") %in% names(x_distribution_rule))) {
      stop("x_distribution_rule must contain nodes and weights.")
    }
    list(
      nodes = as.numeric(x_distribution_rule$nodes),
      weights = as.numeric(x_distribution_rule$weights)
    )
  }
  if (length(x_rule$nodes) < 2L ||
      length(x_rule$weights) != length(x_rule$nodes) ||
      any(!is.finite(x_rule$nodes)) || any(!is.finite(x_rule$weights)) ||
      any(x_rule$weights < 0) || sum(x_rule$weights) <= 0) {
    stop("The population X integration rule is invalid.")
  }
  ordering <- order(x_rule$nodes)
  x_rule$nodes <- x_rule$nodes[ordering]
  x_rule$weights <- x_rule$weights[ordering]
  x_nodes <- length(x_rule$nodes)
  u_rule <- gauss_legendre_rule(u_nodes, -2, 2)
  x_weights <- x_rule$weights
  u_weights <- u_rule$weights * dtnormtc(u_rule$nodes)
  x_weights <- x_weights / sum(x_weights)
  u_weights <- u_weights / sum(u_weights)

  spline_spec <- audit_population_spline_spec(
    t_value = t_value,
    u = u_rule$nodes,
    u_weights = u_weights,
    spline_df = spline_df,
    confounder_strength = confounder_strength,
    treatment_intercept = treatment_intercept,
    x_distribution_rule = list(
      nodes = x_rule$nodes,
      weights = x_weights
    )
  )
  design <- audit_population_spline_matrix(
    x_rule$nodes,
    spline_spec
  )
  relative_y <- seq(-5, 5, length.out = y_nodes)
  y_step <- relative_y[2L] - relative_y[1L]
  y_weights <- rep(y_step, y_nodes)
  y_weights[c(1L, y_nodes)] <- y_step / 2
  gap_density_matrix <- matrix(0, x_nodes, y_nodes)
  cumulative_gap_mass <- matrix(0, x_nodes, y_nodes)
  cumulative_gap_y_moment <- matrix(0, x_nodes, y_nodes)
  r0 <- tau <- m0 <- Lbar <- Ubar <- total_gap_mass <-
    total_gap_y_moment <- mean_y_L <- mean_y_U <-
    exact_q_lower <- exact_q_upper <- sharp_lower_x <-
    sharp_upper_x <- base_x <- numeric(x_nodes)

  if (verbose) {
    message(sprintf(
      paste0(
        "Building population restricted-threshold problem for ",
        "treatment %d (X=%d, U=%d, Y=%d)."
      ),
      t_value,
      x_nodes,
      u_nodes,
      y_nodes
    ))
  }

  for (index in seq_len(x_nodes)) {
    x <- x_rule$nodes[index]
    local <- compute_dgp_oracle_eif_at_x(
      x = x,
      t_value = t_value,
      u = u_rule$nodes,
      u_weights = u_weights,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      confounder_strength = confounder_strength,
      treatment_intercept = treatment_intercept,
      instrument_strength = instrument_strength,
      treatment_effect = treatment_effect,
      envelope_smoothing_epsilon = envelope_smoothing_epsilon,
      y_nodes = y_nodes,
      return_threshold_distribution = TRUE
    )
    summary <- local$summary[1L, , drop = FALSE]
    distribution <- local$distribution
    outcome_center <- treatment_effect * t_value + 1.5 * x
    if (max(abs(distribution$y - (outcome_center + relative_y))) > 1e-10) {
      stop("The conditional outcome grids are not aligned.")
    }

    gap_density <- distribution$probability_mass0 /
      y_weights * distribution$D
    gap_density_matrix[index, ] <- gap_density
    density_left <- gap_density[-y_nodes]
    density_slope <- (
      gap_density[-1L] - density_left
    ) / y_step
    y_left <- distribution$y[-y_nodes]
    cumulative_gap_mass[index, ] <- c(
      0,
      cumsum(y_step * (
        density_left + gap_density[-1L]
      ) / 2)
    )
    cumulative_gap_y_moment[index, ] <- c(
      0,
      cumsum(
        y_left * density_left * y_step +
          (y_left * density_slope + density_left) *
          y_step^2 / 2 +
          density_slope * y_step^3 / 3
      )
    )
    r0[index] <- local$r0
    Lbar[index] <- summary$Lbar
    Ubar[index] <- summary$Ubar
    tau[index] <- (Ubar[index] - 1) / (Ubar[index] - Lbar[index])
    m0[index] <- sum(
      distribution$probability_mass0 * distribution$y
    )
    total_gap_mass[index] <-
      cumulative_gap_mass[index, y_nodes]
    total_gap_y_moment[index] <-
      cumulative_gap_y_moment[index, y_nodes]
    mean_y_L[index] <- sum(
      distribution$probability_mass0 * distribution$y * distribution$L
    )
    mean_y_U[index] <-
      mean_y_L[index] + total_gap_y_moment[index]
    exact_q_lower[index] <- summary$lower_q
    exact_q_upper[index] <- summary$upper_q
    sharp_lower_x[index] <- summary$lower_M
    sharp_upper_x[index] <- summary$upper_M
    base_x[index] <- r0[index] * m0[index] +
      (1 - r0[index]) * (
        tau[index] * mean_y_L[index] +
          (1 - tau[index]) * mean_y_U[index]
    )
  }

  outcome_center <- treatment_effect * t_value + 1.5 * x_rule$nodes
  exact_q_lower <- audit_invert_linear_cumulative(
    cumulative = cumulative_gap_mass,
    density = gap_density_matrix,
    target = 1 - Lbar,
    grid = relative_y,
    outcome_center = outcome_center
  )
  exact_q_upper <- audit_invert_linear_cumulative(
    cumulative = cumulative_gap_mass,
    density = gap_density_matrix,
    target = Ubar - 1,
    grid = relative_y,
    outcome_center = outcome_center
  )

  problem <- list(
    t_value = t_value,
    x = x_rule$nodes,
    x_weights = x_weights,
    r0 = r0,
    tau = tau,
    m0 = m0,
    Lbar = Lbar,
    Ubar = Ubar,
    total_gap_mass = total_gap_mass,
    total_gap_y_moment = total_gap_y_moment,
    mean_y_L = mean_y_L,
    mean_y_U = mean_y_U,
    gap_density = gap_density_matrix,
    cumulative_gap_mass = cumulative_gap_mass,
    cumulative_gap_y_moment = cumulative_gap_y_moment,
    relative_y = relative_y,
    outcome_center = outcome_center,
    base_x = base_x,
    exact_q_lower = exact_q_lower,
    exact_q_upper = exact_q_upper,
    sharp_lower_x = sharp_lower_x,
    sharp_upper_x = sharp_upper_x,
    design = design,
    spline_spec = spline_spec,
    settings = list(
      x_nodes = x_nodes,
      u_nodes = u_nodes,
      y_nodes = y_nodes,
      spline_df = spline_df,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      confounder_strength = confounder_strength,
      treatment_intercept = treatment_intercept,
      instrument_strength = instrument_strength,
      treatment_effect = treatment_effect,
      envelope_smoothing_epsilon = envelope_smoothing_epsilon
    )
  )
  problem$sharp_direct_lower_x <- problem$sharp_lower_x
  problem$sharp_direct_upper_x <- problem$sharp_upper_x
  exact_lower <- evaluate_population_restricted_threshold(
    problem,
    problem$exact_q_lower,
    "lower"
  )
  exact_upper <- evaluate_population_restricted_threshold(
    problem,
    problem$exact_q_upper,
    "upper"
  )
  problem$sharp_lower_x <- exact_lower$pseudo_x
  problem$sharp_upper_x <- exact_upper$pseudo_x
  problem$sharp_lower_direct_quadrature_error <-
    sum(problem$x_weights * (
      problem$sharp_direct_lower_x - problem$sharp_lower_x
    ))
  problem$sharp_upper_direct_quadrature_error <-
    sum(problem$x_weights * (
      problem$sharp_direct_upper_x - problem$sharp_upper_x
    ))
  problem
}

audit_integrate_gap_to_threshold <- function(
    cumulative_mass,
    cumulative_y_moment,
    density,
    value,
    grid,
    outcome_center
) {
  row_count <- nrow(cumulative_mass)
  if (length(value) != row_count) {
    stop("value must have one entry per cumulative-integral row.")
  }
  mass <- y_moment <- density_at_value <- numeric(row_count)
  above <- value >= grid[length(grid)]
  mass[above] <- cumulative_mass[above, length(grid)]
  y_moment[above] <- cumulative_y_moment[above, length(grid)]
  interior <- value > grid[1L] & value < grid[length(grid)]
  if (any(interior)) {
    rows <- which(interior)
    left <- findInterval(value[interior], grid)
    interval_width <- grid[left + 1L] - grid[left]
    distance <- value[interior] - grid[left]
    g0 <- density[cbind(rows, left)]
    slope <- (
      density[cbind(rows, left + 1L)] - g0
    ) / interval_width
    y0 <- outcome_center[rows] + grid[left]
    mass[interior] <- cumulative_mass[cbind(rows, left)] +
      g0 * distance + slope * distance^2 / 2
    density_at_value[interior] <- g0 + slope * distance
    y_moment[interior] <-
      cumulative_y_moment[cbind(rows, left)] +
      y0 * g0 * distance +
      (y0 * slope + g0) * distance^2 / 2 +
      slope * distance^3 / 3
  }
  list(
    mass = mass,
    y_moment = y_moment,
    density = density_at_value
  )
}

evaluate_population_restricted_threshold <- function(
    problem,
    q,
    direction = c("lower", "upper")
) {
  direction <- match.arg(direction)
  if (length(q) != length(problem$x) || any(!is.finite(q))) {
    stop("q must be finite with one value per population X node.")
  }

  relative_q <- q - problem$outcome_center
  lower_gap <- audit_integrate_gap_to_threshold(
    cumulative_mass = problem$cumulative_gap_mass,
    cumulative_y_moment = problem$cumulative_gap_y_moment,
    density = problem$gap_density,
    value = relative_q,
    grid = problem$relative_y,
    outcome_center = problem$outcome_center
  )
  lower_gap_mass <- lower_gap$mass
  lower_gap_y_moment <- lower_gap$y_moment
  level <- if (direction == "lower") {
    1 - problem$tau
  } else {
    problem$tau
  }
  upper_gap_mass <- problem$total_gap_mass - lower_gap_mass
  upper_gap_y_moment <-
    problem$total_gap_y_moment - lower_gap_y_moment
  check_risk_x <-
    (1 - level) * (q * lower_gap_mass - lower_gap_y_moment) +
    level * (upper_gap_y_moment - q * upper_gap_mass)
  risk_weight <- problem$x_weights * (1 - problem$r0)
  risk <- sum(risk_weight * check_risk_x)
  risk_normalizer <- sum(risk_weight * problem$total_gap_mass)

  if (direction == "lower") {
    pseudo_x <- problem$base_x -
      (1 - problem$r0) * check_risk_x
    theta_x <- problem$mean_y_L + lower_gap_y_moment
    lambda_normalization_x <- problem$Lbar + lower_gap_mass
  } else {
    pseudo_x <- problem$base_x +
      (1 - problem$r0) * check_risk_x
    theta_x <- problem$mean_y_U - lower_gap_y_moment
    lambda_normalization_x <- problem$Ubar - lower_gap_mass
  }
  direct_x <- problem$r0 * problem$m0 +
    (1 - problem$r0) * theta_x
  score_residual <- lower_gap_mass -
    level * problem$total_gap_mass
  score <- colSums(
    problem$design *
      (risk_weight * score_residual)
  )
  score_jacobian <- crossprod(
    problem$design,
    problem$design * (risk_weight * lower_gap$density)
  )

  list(
    q = q,
    pseudo_bound = sum(problem$x_weights * pseudo_x),
    direct_bound = sum(problem$x_weights * direct_x),
    risk = risk,
    risk_normalizer = risk_normalizer,
    normalized_risk = risk / risk_normalizer,
    score = score,
    score_jacobian = score_jacobian,
    max_abs_score = max(abs(score)),
    lambda_normalization_x = lambda_normalization_x,
    weighted_normalization_error = sum(
      risk_weight * (lambda_normalization_x - 1)
    ),
    q_weighted_normalization_error = sum(
      risk_weight * q * (lambda_normalization_x - 1)
    ),
    direct_minus_pseudo =
      sum(problem$x_weights * direct_x) -
      sum(problem$x_weights * pseudo_x),
    check_risk_x = check_risk_x,
    pseudo_x = pseudo_x,
    direct_x = direct_x
  )
}

fit_population_restricted_threshold <- function(
    problem,
    direction = c("lower", "upper"),
    maxit = 12000L,
    reltol = 1e-12
) {
  direction <- match.arg(direction)
  exact_q <- if (direction == "lower") {
    problem$exact_q_lower
  } else {
    problem$exact_q_upper
  }
  initial_weight <- problem$x_weights *
    (1 - problem$r0) * problem$total_gap_mass
  initial_fit <- stats::lm.wfit(
    x = problem$design,
    y = exact_q,
    w = pmax(initial_weight, 1e-12)
  )
  initial <- initial_fit$coefficients
  initial[!is.finite(initial)] <- 0

  objective <- function(beta) {
    q <- as.vector(problem$design %*% beta)
    evaluate_population_restricted_threshold(
      problem,
      q,
      direction
    )$normalized_risk
  }
  gradient <- function(beta) {
    q <- as.vector(problem$design %*% beta)
    evaluation <- evaluate_population_restricted_threshold(
      problem,
      q,
      direction
    )
    evaluation$score / evaluation$risk_normalizer
  }
  nelder_mead <- stats::optim(
    par = initial,
    fn = objective,
    method = "Nelder-Mead",
    control = list(maxit = maxit, reltol = reltol)
  )
  bfgs <- stats::optim(
    par = nelder_mead$par,
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(maxit = maxit, reltol = reltol)
  )
  newton_coefficients <- bfgs$par
  newton_convergence <- 1L
  newton_iterations <- 0L
  for (iteration in seq_len(100L)) {
    newton_iterations <- iteration
    q <- as.vector(problem$design %*% newton_coefficients)
    evaluation <- evaluate_population_restricted_threshold(
      problem,
      q,
      direction
    )
    if (evaluation$max_abs_score < 1e-11) {
      newton_convergence <- 0L
      break
    }
    step <- tryCatch(
      solve(
        evaluation$score_jacobian +
          diag(1e-12, ncol(problem$design)),
        evaluation$score
      ),
      error = function(condition) NULL
    )
    if (is.null(step) || any(!is.finite(step))) break
    current_objective <- evaluation$normalized_risk
    step_fraction <- 1
    accepted <- FALSE
    for (line_search in seq_len(30L)) {
      candidate <- newton_coefficients - step_fraction * step
      candidate_objective <- objective(candidate)
      if (
        is.finite(candidate_objective) &&
          candidate_objective <= current_objective
      ) {
        newton_coefficients <- candidate
        accepted <- TRUE
        break
      }
      step_fraction <- step_fraction / 2
    }
    if (!accepted) break
  }
  candidate_coefficients <- list(
    initial = initial,
    nelder_mead = nelder_mead$par,
    bfgs = bfgs$par,
    newton = newton_coefficients
  )
  candidate_objectives <- vapply(
    candidate_coefficients,
    objective,
    numeric(1)
  )
  selected_method <- names(which.min(candidate_objectives))[1L]
  coefficients <- candidate_coefficients[[selected_method]]
  q <- as.vector(problem$design %*% coefficients)
  evaluation <- evaluate_population_restricted_threshold(
    problem,
    q,
    direction
  )

  perturbation <- 1e-5
  local_objective_differences <- unlist(lapply(
    seq_along(coefficients),
    function(index) {
      step <- rep(0, length(coefficients))
      step[index] <- perturbation
      c(
        objective(coefficients + step) - evaluation$normalized_risk,
        objective(coefficients - step) - evaluation$normalized_risk
      )
    }
  ))

  list(
    direction = direction,
    coefficients = coefficients,
    q = q,
    evaluation = evaluation,
    selected_method = selected_method,
    candidate_objectives = candidate_objectives,
    nelder_mead_convergence = nelder_mead$convergence,
    bfgs_convergence = bfgs$convergence,
    newton_convergence = newton_convergence,
    newton_iterations = newton_iterations,
    minimum_local_objective_difference =
      min(local_objective_differences),
    exact_q_rmse = sqrt(sum(
      problem$x_weights * (q - exact_q)^2
    )),
    exact_q_risk_weighted_rmse = sqrt(
      sum(initial_weight * (q - exact_q)^2) /
        sum(initial_weight)
    )
  )
}

compute_population_restricted_spline_bounds <- function(
    x_nodes = 80L,
    u_nodes = 80L,
    y_nodes = 4001L,
    spline_df = 5L,
    Lambda0 = 5,
    Lambda1 = 5,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    envelope_smoothing_epsilon = 0,
    x_distribution_rule = NULL,
    verbose = TRUE
) {
  arms <- list()
  coefficient_rows <- list()
  diagnostic_rows <- list()

  for (t_value in c(0, 1)) {
    problem <- build_population_restricted_threshold_problem(
      t_value = t_value,
      x_nodes = x_nodes,
      u_nodes = u_nodes,
      y_nodes = y_nodes,
      spline_df = spline_df,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      confounder_strength = confounder_strength,
      treatment_intercept = treatment_intercept,
      instrument_strength = instrument_strength,
      treatment_effect = treatment_effect,
      envelope_smoothing_epsilon = envelope_smoothing_epsilon,
      x_distribution_rule = x_distribution_rule,
      verbose = verbose
    )
    lower <- fit_population_restricted_threshold(problem, "lower")
    upper <- fit_population_restricted_threshold(problem, "upper")
    sharp_lower <- sum(problem$x_weights * problem$sharp_lower_x)
    sharp_upper <- sum(problem$x_weights * problem$sharp_upper_x)

    arms[[t_value + 1L]] <- data.frame(
      treatment = t_value,
      restricted_lower = lower$evaluation$pseudo_bound,
      restricted_upper = upper$evaluation$pseudo_bound,
      restricted_lower_direct = lower$evaluation$direct_bound,
      restricted_upper_direct = upper$evaluation$direct_bound,
      sharp_lower = sharp_lower,
      sharp_upper = sharp_upper,
      lower_restriction_gap =
        lower$evaluation$pseudo_bound - sharp_lower,
      upper_restriction_gap =
        upper$evaluation$pseudo_bound - sharp_upper,
      lower_q_rmse = lower$exact_q_rmse,
      upper_q_rmse = upper$exact_q_rmse,
      lower_q_risk_weighted_rmse =
        lower$exact_q_risk_weighted_rmse,
      upper_q_risk_weighted_rmse =
        upper$exact_q_risk_weighted_rmse,
      lower_direct_minus_pseudo =
        lower$evaluation$direct_minus_pseudo,
      upper_direct_minus_pseudo =
        upper$evaluation$direct_minus_pseudo,
      lower_max_abs_score = lower$evaluation$max_abs_score,
      upper_max_abs_score = upper$evaluation$max_abs_score,
      lower_minimum_local_objective_difference =
        lower$minimum_local_objective_difference,
      upper_minimum_local_objective_difference =
        upper$minimum_local_objective_difference,
      lower_selected_method = lower$selected_method,
      upper_selected_method = upper$selected_method,
      lower_nelder_mead_convergence =
        lower$nelder_mead_convergence,
      upper_nelder_mead_convergence =
        upper$nelder_mead_convergence,
      lower_bfgs_convergence = lower$bfgs_convergence,
      upper_bfgs_convergence = upper$bfgs_convergence,
      lower_newton_convergence = lower$newton_convergence,
      upper_newton_convergence = upper$newton_convergence,
      lower_newton_iterations = lower$newton_iterations,
      upper_newton_iterations = upper$newton_iterations,
      stringsAsFactors = FALSE
    )

    coefficient_rows[[t_value + 1L]] <- rbind(
      data.frame(
        treatment = t_value,
        direction = "lower",
        coefficient = names(lower$coefficients),
        value = unname(lower$coefficients),
        stringsAsFactors = FALSE
      ),
      data.frame(
        treatment = t_value,
        direction = "upper",
        coefficient = names(upper$coefficients),
        value = unname(upper$coefficients),
        stringsAsFactors = FALSE
      )
    )
    diagnostic_rows[[t_value + 1L]] <- rbind(
      data.frame(
        treatment = t_value,
        direction = "lower",
        x = problem$x,
        x_weight = problem$x_weights,
        tau = problem$tau,
        q_exact = problem$exact_q_lower,
        q_restricted = lower$q,
        lambda_normalization =
          lower$evaluation$lambda_normalization_x,
        pseudo_x = lower$evaluation$pseudo_x,
        direct_x = lower$evaluation$direct_x,
        stringsAsFactors = FALSE
      ),
      data.frame(
        treatment = t_value,
        direction = "upper",
        x = problem$x,
        x_weight = problem$x_weights,
        tau = problem$tau,
        q_exact = problem$exact_q_upper,
        q_restricted = upper$q,
        lambda_normalization =
          upper$evaluation$lambda_normalization_x,
        pseudo_x = upper$evaluation$pseudo_x,
        direct_x = upper$evaluation$direct_x,
        stringsAsFactors = FALSE
      )
    )
    attr(coefficient_rows[[t_value + 1L]], "spline_spec") <-
      problem$spline_spec
  }

  arms <- do.call(rbind, arms)
  mu0 <- arms[arms$treatment == 0, , drop = FALSE]
  mu1 <- arms[arms$treatment == 1, , drop = FALSE]
  bounds <- data.frame(
    mu1_lower = mu1$restricted_lower,
    mu1_upper = mu1$restricted_upper,
    mu0_lower = mu0$restricted_lower,
    mu0_upper = mu0$restricted_upper,
    ate_lower = mu1$restricted_lower - mu0$restricted_upper,
    ate_upper = mu1$restricted_upper - mu0$restricted_lower,
    sharp_mu1_lower = mu1$sharp_lower,
    sharp_mu1_upper = mu1$sharp_upper,
    sharp_mu0_lower = mu0$sharp_lower,
    sharp_mu0_upper = mu0$sharp_upper,
    sharp_ate_lower = mu1$sharp_lower - mu0$sharp_upper,
    sharp_ate_upper = mu1$sharp_upper - mu0$sharp_lower,
    x_nodes = x_nodes,
    u_nodes = u_nodes,
    y_nodes = y_nodes,
    spline_df = spline_df,
    stringsAsFactors = FALSE
  )

  list(
    bounds = bounds,
    arms = arms,
    coefficients = do.call(rbind, coefficient_rows),
    diagnostics = do.call(rbind, diagnostic_rows),
    spline_specs = lapply(coefficient_rows, attr, "spline_spec")
  )
}

predict_dgp_oracle_lookup_column <- function(
    lookup,
    x,
    t_value,
    column
) {
  grid <- lookup$grid[
    lookup$grid$treatment == t_value,
    ,
    drop = FALSE
  ]

  if (!column %in% names(grid)) {
    stop("Unknown DGP-oracle lookup column: ", column)
  }

  stats::approx(
    x = grid$x,
    y = grid[[column]],
    xout = x,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )$y
}

predict_dgp_oracle_eif_lookup <- function(
    lookup,
    x,
    t_value,
    endpoint = c("lower", "upper"),
    component = c("q", "M", "C0", "C1", "Hbar0", "Hbar1")
) {
  endpoint <- match.arg(endpoint)
  component <- match.arg(component)
  column <- paste0(endpoint, "_", component)
  predict_dgp_oracle_lookup_column(
    lookup = lookup,
    x = x,
    t_value = t_value,
    column = column
  )
}

compute_dgp_oracle_eif_signals <- function(data, lookup) {
  required <- c(
    "Y", "T_obs", "Z", "X", "e_true", "pi0_true", "pi1_true",
    "gamma0_oracle", "gamma1_oracle"
  )
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Missing oracle columns: ", paste(missing, collapse = ", "))
  }

  Lambda0 <- lookup$settings$Lambda0
  Lambda1 <- lookup$settings$Lambda1
  arm_signals <- list()
  component_rows <- list()
  component_counter <- 1L

  for (t_value in c(0, 1)) {
    nuisance <- list(
      e = data$e_true,
      pi0 = data$pi0_true,
      pi1 = data$pi1_true
    )
    gamma <- if (t_value == 1) {
      data$gamma1_oracle
    } else {
      data$gamma0_oracle
    }
    envelope <- compute_projected_envelopes(
      nuisance = nuisance,
      gamma = gamma,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      gamma_compatibility_action = "stop"
    )
    q_lower <- predict_dgp_oracle_eif_lookup(
      lookup, data$X, t_value, "lower", "q"
    )
    q_upper <- predict_dgp_oracle_eif_lookup(
      lookup, data$X, t_value, "upper", "q"
    )
    primitive <- compute_eif_row_primitives(
      data = data,
      nuisance = nuisance,
      envelope = envelope,
      t_value = t_value,
      Lambda0 = Lambda0,
      Lambda1 = Lambda1,
      q_lower = q_lower,
      q_upper = q_upper
    )
    observed_arm <- as.numeric(data$T_obs == t_value)
    instrument0 <- as.numeric(data$Z == 0)
    instrument1 <- 1 - instrument0
    r0 <- primitive$r0
    r1 <- primitive$r1

    for (endpoint in c("lower", "upper")) {
      local <- primitive[[endpoint]]
      component <- function(name) {
        predict_dgp_oracle_eif_lookup(
          lookup, data$X, t_value, endpoint, name
        )
      }
      M <- component("M")
      C0 <- component("C0")
      C1 <- component("C1")
      Hbar0 <- component("Hbar0")
      Hbar1 <- component("Hbar1")
      treatment_correction <-
        instrument0 / nuisance$e * C0 * (observed_arm - r0) +
        instrument1 / (1 - nuisance$e) * C1 * (observed_arm - r1)
      outcome_correction <-
        instrument0 / nuisance$e * observed_arm / r0 *
        (local$H0 - Hbar0) +
        instrument1 / (1 - nuisance$e) * observed_arm / r1 *
        (local$H1 - Hbar1)
      signal <- M + treatment_correction + outcome_correction
      arm_signals[[paste0("mu", t_value, "_", endpoint)]] <- signal
      component_rows[[component_counter]] <- data.frame(
        row_id = seq_len(nrow(data)),
        treatment = t_value,
        endpoint = endpoint,
        M = M,
        treatment_correction = treatment_correction,
        outcome_correction = outcome_correction,
        signal = signal,
        stringsAsFactors = FALSE
      )
      component_counter <- component_counter + 1L
    }
  }

  list(
    arm_signals = as.data.frame(arm_signals),
    components = do.call(rbind, component_rows)
  )
}

compute_crossfit_oracle_eif_signals <- function(
    data,
    lookup,
    K = 5L,
    seed = 1L,
    conditional_spline_df = 5L,
    threshold_spline_df = 5L
) {
  n <- nrow(data)
  data$.row_id <- seq_len(n)
  fold_id <- make_tz_stratified_folds(data, K = K, seed = seed)
  versions <- c(
    "oracle_q_spline",
    "oracle_tau_crossfit",
    "current_crossfit"
  )
  storage <- lapply(versions, function(version) {
    list(
      mu0_lower = rep(NA_real_, n),
      mu0_upper = rep(NA_real_, n),
      mu1_lower = rep(NA_real_, n),
      mu1_upper = rep(NA_real_, n)
    )
  })
  names(storage) <- versions
  component_storage <- lapply(versions, function(version) {
    output <- list()
    for (t_value in c(0, 1)) {
      for (endpoint in c("lower", "upper")) {
        output[[paste0("mu", t_value, "_", endpoint)]] <- list(
          M = rep(NA_real_, n),
          treatment_correction = rep(NA_real_, n),
          outcome_correction = rep(NA_real_, n)
        )
      }
    }
    output
  })
  names(component_storage) <- versions
  diagnostic_rows <- list()
  diagnostic_counter <- 1L

  for (fold in seq_len(K)) {
    validation_index <- which(fold_id == fold)
    training_index <- which(fold_id != fold)
    train <- data[training_index, , drop = FALSE]
    validation <- data[validation_index, , drop = FALSE]

    for (t_value in c(0, 1)) {
      nuisance_train <- list(
        e = train$e_true,
        pi0 = train$pi0_true,
        pi1 = train$pi1_true
      )
      nuisance_validation <- list(
        e = validation$e_true,
        pi0 = validation$pi0_true,
        pi1 = validation$pi1_true
      )
      gamma_train <- if (t_value == 1) {
        train$gamma1_oracle
      } else {
        train$gamma0_oracle
      }
      gamma_validation <- if (t_value == 1) {
        validation$gamma1_oracle
      } else {
        validation$gamma0_oracle
      }
      envelope_train <- compute_projected_envelopes(
        nuisance = nuisance_train,
        gamma = gamma_train,
        t_value = t_value,
        Lambda0 = lookup$settings$Lambda0,
        Lambda1 = lookup$settings$Lambda1,
        gamma_compatibility_action = "stop"
      )
      envelope_validation <- compute_projected_envelopes(
        nuisance = nuisance_validation,
        gamma = gamma_validation,
        t_value = t_value,
        Lambda0 = lookup$settings$Lambda0,
        Lambda1 = lookup$settings$Lambda1,
        gamma_compatibility_action = "stop"
      )
      q_lower_train <- predict_dgp_oracle_eif_lookup(
        lookup, train$X, t_value, "lower", "q"
      )
      q_upper_train <- predict_dgp_oracle_eif_lookup(
        lookup, train$X, t_value, "upper", "q"
      )
      q_lower_validation <- predict_dgp_oracle_eif_lookup(
        lookup, validation$X, t_value, "lower", "q"
      )
      q_upper_validation <- predict_dgp_oracle_eif_lookup(
        lookup, validation$X, t_value, "upper", "q"
      )

      oracle_q <- fit_eif_arm_fold(
        train = train,
        validation = validation,
        row_id_validation = validation$.row_id,
        nuisance_train = nuisance_train,
        nuisance_validation = nuisance_validation,
        envelope_train = envelope_train,
        envelope_validation = envelope_validation,
        t_value = t_value,
        Lambda0 = lookup$settings$Lambda0,
        Lambda1 = lookup$settings$Lambda1,
        q_lower_train = q_lower_train,
        q_upper_train = q_upper_train,
        q_lower_validation = q_lower_validation,
        q_upper_validation = q_upper_validation,
        x_cols = "X",
        conditional_spline_df = conditional_spline_df
      )

      Lbar_oracle_train <- predict_dgp_oracle_lookup_column(
        lookup,
        train$X,
        t_value,
        "Lbar"
      )
      Ubar_oracle_train <- predict_dgp_oracle_lookup_column(
        lookup,
        train$X,
        t_value,
        "Ubar"
      )
      tau_oracle_train <- (
        Ubar_oracle_train - 1
      ) / (
        Ubar_oracle_train - Lbar_oracle_train
      )
      baseline_train <- train$T_obs == t_value & train$Z == 0
      target_odds_train <- compute_target_odds(
        nuisance_train,
        t_value
      )
      threshold_weight <-
        1 / nuisance_train$e[baseline_train] *
        target_odds_train[baseline_train] *
        envelope_train$D[baseline_train]

      oracle_tau_lower_fit <- fit_threshold_by_check_loss(
        data = train[baseline_train, , drop = FALSE],
        y = train$Y[baseline_train],
        level = 1 - tau_oracle_train[baseline_train],
        weight = threshold_weight,
        x_cols = "X",
        spline_df = threshold_spline_df
      )
      oracle_tau_upper_fit <- fit_threshold_by_check_loss(
        data = train[baseline_train, , drop = FALSE],
        y = train$Y[baseline_train],
        level = tau_oracle_train[baseline_train],
        weight = threshold_weight,
        x_cols = "X",
        spline_df = threshold_spline_df
      )
      oracle_tau_q_lower_train <- oracle_tau_lower_fit$predict(train)
      oracle_tau_q_upper_train <- oracle_tau_upper_fit$predict(train)
      oracle_tau_q_lower_validation <-
        oracle_tau_lower_fit$predict(validation)
      oracle_tau_q_upper_validation <-
        oracle_tau_upper_fit$predict(validation)

      oracle_tau <- fit_eif_arm_fold(
        train = train,
        validation = validation,
        row_id_validation = validation$.row_id,
        nuisance_train = nuisance_train,
        nuisance_validation = nuisance_validation,
        envelope_train = envelope_train,
        envelope_validation = envelope_validation,
        t_value = t_value,
        Lambda0 = lookup$settings$Lambda0,
        Lambda1 = lookup$settings$Lambda1,
        q_lower_train = oracle_tau_q_lower_train,
        q_upper_train = oracle_tau_q_upper_train,
        q_lower_validation = oracle_tau_q_lower_validation,
        q_upper_validation = oracle_tau_q_upper_validation,
        x_cols = "X",
        conditional_spline_df = conditional_spline_df
      )

      current <- fit_bound_arm_fold(
        train = train,
        validation = validation,
        row_id_validation = validation$.row_id,
        nuisance_train = nuisance_train,
        nuisance_validation = nuisance_validation,
        gamma_train = gamma_train,
        gamma_validation = gamma_validation,
        t_value = t_value,
        Lambda0 = lookup$settings$Lambda0,
        Lambda1 = lookup$settings$Lambda1,
        mode = "oracle_all",
        x_cols = "X",
        conditional_mean_spline_df = conditional_spline_df,
        threshold_spline_df = threshold_spline_df,
        compatibility_action = "stop",
        compute_eif = TRUE,
        eif_conditional_spline_df = conditional_spline_df
      )

      for (endpoint in c("lower", "upper")) {
        key <- paste0("mu", t_value, "_", endpoint)
        for (version in versions) {
          signal_data <- switch(
            version,
            oracle_q_spline = oracle_q$signals,
            oracle_tau_crossfit = oracle_tau$signals,
            current_crossfit = current$eif$signals
          )
          storage[[version]][[key]][validation_index] <-
            signal_data[[paste0(endpoint, "_signal")]]
          component_storage[[version]][[key]]$M[validation_index] <-
            signal_data[[paste0(endpoint, "_M")]]
          component_storage[[version]][[key]]$
            treatment_correction[validation_index] <-
            signal_data[[paste0(
              endpoint,
              "_treatment_correction"
            )]]
          component_storage[[version]][[key]]$
            outcome_correction[validation_index] <-
            signal_data[[paste0(
              endpoint,
              "_outcome_correction"
            )]]
        }
      }

      prediction_index <- match(
        current$predictions$row_id,
        validation$.row_id
      )
      estimated_q <- data.frame(
        q_lower = current$predictions$q_lower,
        q_upper = current$predictions$q_upper,
        q_lower_true = q_lower_validation[prediction_index],
        q_upper_true = q_upper_validation[prediction_index],
        q_lower_oracle_tau =
          oracle_tau_q_lower_validation[prediction_index],
        q_upper_oracle_tau =
          oracle_tau_q_upper_validation[prediction_index]
      )
      Lbar_oracle_validation <- predict_dgp_oracle_lookup_column(
        lookup,
        validation$X,
        t_value,
        "Lbar"
      )
      Ubar_oracle_validation <- predict_dgp_oracle_lookup_column(
        lookup,
        validation$X,
        t_value,
        "Ubar"
      )
      tau_oracle_validation <- (
        Ubar_oracle_validation - 1
      ) / (
        Ubar_oracle_validation - Lbar_oracle_validation
      )
      diagnostic_rows[[diagnostic_counter]] <- data.frame(
        fold = fold,
        treatment = t_value,
        tau_bias = mean(
          current$predictions$tau -
            tau_oracle_validation[prediction_index]
        ),
        tau_rmse = sqrt(mean(
          (current$predictions$tau -
             tau_oracle_validation[prediction_index])^2
        )),
        lower_q_bias = mean(
          estimated_q$q_lower - estimated_q$q_lower_true
        ),
        lower_q_rmse = sqrt(mean(
          (estimated_q$q_lower - estimated_q$q_lower_true)^2
        )),
        upper_q_bias = mean(
          estimated_q$q_upper - estimated_q$q_upper_true
        ),
        upper_q_rmse = sqrt(mean(
          (estimated_q$q_upper - estimated_q$q_upper_true)^2
        )),
        oracle_tau_lower_q_bias = mean(
          estimated_q$q_lower_oracle_tau -
            estimated_q$q_lower_true
        ),
        oracle_tau_lower_q_rmse = sqrt(mean(
          (estimated_q$q_lower_oracle_tau -
             estimated_q$q_lower_true)^2
        )),
        oracle_tau_upper_q_bias = mean(
          estimated_q$q_upper_oracle_tau -
            estimated_q$q_upper_true
        ),
        oracle_tau_upper_q_rmse = sqrt(mean(
          (estimated_q$q_upper_oracle_tau -
             estimated_q$q_upper_true)^2
        )),
        oracle_tau_lower_convergence =
          oracle_tau_lower_fit$convergence,
        oracle_tau_upper_convergence =
          oracle_tau_upper_fit$convergence,
        stringsAsFactors = FALSE
      )
      diagnostic_counter <- diagnostic_counter + 1L
    }
  }

  components <- lapply(versions, function(version) {
    rows <- list()
    counter <- 1L
    for (t_value in c(0, 1)) {
      for (endpoint in c("lower", "upper")) {
        key <- paste0("mu", t_value, "_", endpoint)
        local <- component_storage[[version]][[key]]
        rows[[counter]] <- data.frame(
          row_id = seq_len(n),
          treatment = t_value,
          endpoint = endpoint,
          M = local$M,
          treatment_correction = local$treatment_correction,
          outcome_correction = local$outcome_correction,
          signal = storage[[version]][[key]],
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }
    do.call(rbind, rows)
  })
  names(components) <- versions

  list(
    signals = lapply(storage, as.data.frame),
    components = components,
    threshold_diagnostics = do.call(rbind, diagnostic_rows),
    fold_id = fold_id
  )
}

compute_oracle_tau_threshold_learner_audit <- function(
    data,
    lookup,
    learners = c("spline_fixed", "spline_cv", "rkhs"),
    K = 5L,
    seed = 1L,
    conditional_spline_df = 5L,
    threshold_spline_df = 5L,
    threshold_control = NULL
) {
  learners <- unique(vapply(
    learners,
    match.arg,
    character(1),
    choices = c("spline_fixed", "spline_cv", "rkhs")
  ))
  required_oracle <- c(
    "e_true", "pi0_true", "pi1_true",
    "gamma0_oracle", "gamma1_oracle"
  )
  if (!all(required_oracle %in% names(data))) {
    stop(
      "The threshold learner audit requires oracle nuisance columns."
    )
  }
  n <- nrow(data)
  data$.row_id <- seq_len(n)
  fold_id <- make_tz_stratified_folds(data, K = K, seed = seed)
  endpoint_keys <- c(
    "mu0_lower", "mu0_upper", "mu1_lower", "mu1_upper"
  )
  initialize_storage <- function() {
    output <- lapply(endpoint_keys, function(name) rep(NA_real_, n))
    names(output) <- endpoint_keys
    output
  }
  signal_storage <- lapply(learners, function(x) initialize_storage())
  direct_storage <- lapply(learners, function(x) initialize_storage())
  pseudo_storage <- lapply(learners, function(x) initialize_storage())
  names(signal_storage) <- names(direct_storage) <-
    names(pseudo_storage) <- learners
  diagnostic_rows <- list()
  diagnostic_counter <- 1L
  cv_rows <- list()
  cv_counter <- 1L

  fit_iterations <- function(fit) {
    if ("iterations" %in% names(fit$counts)) {
      return(unname(fit$counts[["iterations"]]))
    }
    if ("function" %in% names(fit$counts)) {
      return(unname(fit$counts[["function"]]))
    }
    NA_real_
  }

  for (fold in seq_len(K)) {
    validation_index <- which(fold_id == fold)
    training_index <- which(fold_id != fold)
    train <- data[training_index, , drop = FALSE]
    validation <- data[validation_index, , drop = FALSE]

    for (t_value in c(0, 1)) {
      nuisance_train <- list(
        e = train$e_true,
        pi0 = train$pi0_true,
        pi1 = train$pi1_true
      )
      nuisance_validation <- list(
        e = validation$e_true,
        pi0 = validation$pi0_true,
        pi1 = validation$pi1_true
      )
      gamma_train <- if (t_value == 1) {
        train$gamma1_oracle
      } else {
        train$gamma0_oracle
      }
      gamma_validation <- if (t_value == 1) {
        validation$gamma1_oracle
      } else {
        validation$gamma0_oracle
      }
      envelope_train <- compute_projected_envelopes(
        nuisance_train,
        gamma_train,
        t_value,
        lookup$settings$Lambda0,
        lookup$settings$Lambda1,
        gamma_compatibility_action = "stop"
      )
      envelope_validation <- compute_projected_envelopes(
        nuisance_validation,
        gamma_validation,
        t_value,
        lookup$settings$Lambda0,
        lookup$settings$Lambda1,
        gamma_compatibility_action = "stop"
      )
      Lbar_train <- predict_dgp_oracle_lookup_column(
        lookup, train$X, t_value, "Lbar"
      )
      Ubar_train <- predict_dgp_oracle_lookup_column(
        lookup, train$X, t_value, "Ubar"
      )
      Lbar_validation <- predict_dgp_oracle_lookup_column(
        lookup, validation$X, t_value, "Lbar"
      )
      Ubar_validation <- predict_dgp_oracle_lookup_column(
        lookup, validation$X, t_value, "Ubar"
      )
      tau_train <- (Ubar_train - 1) / (Ubar_train - Lbar_train)
      tau_validation <-
        (Ubar_validation - 1) / (Ubar_validation - Lbar_validation)
      exact_q_lower_validation <- predict_dgp_oracle_eif_lookup(
        lookup, validation$X, t_value, "lower", "q"
      )
      exact_q_upper_validation <- predict_dgp_oracle_eif_lookup(
        lookup, validation$X, t_value, "upper", "q"
      )
      baseline_train <- train$T_obs == t_value & train$Z == 0
      baseline_validation <-
        validation$T_obs == t_value & validation$Z == 0
      target_odds_train <- compute_target_odds(nuisance_train, t_value)
      threshold_weight <- 1 / nuisance_train$e[baseline_train] *
        target_odds_train[baseline_train] *
        envelope_train$D[baseline_train]

      for (learner in learners) {
        started <- proc.time()[["elapsed"]]
        learner_offset <- match(
          learner,
          c("spline_fixed", "spline_cv", "rkhs")
        )
        learner_seed <- seed + 100000L * fold +
          1000L * t_value + learner_offset * 10L
        local_control <- resolve_threshold_control(threshold_control)
        local_rkhs_plan <- if (learner == "rkhs") {
          make_threshold_rkhs_plan(
            data = train[baseline_train, , drop = FALSE],
            x_cols = "X",
            backend = local_control$rkhs_backend,
            bandwidth_multipliers =
              local_control$rkhs_bandwidth_multipliers,
            rank = local_control$rkhs_rank,
            cv_folds = local_control$cv_folds,
            bandwidth_sample = local_control$rkhs_bandwidth_sample,
            seed = learner_seed
          )
        } else {
          NULL
        }
        lower_fit <- fit_threshold_learner(
          data = train[baseline_train, , drop = FALSE],
          y = train$Y[baseline_train],
          level = 1 - tau_train[baseline_train],
          weight = threshold_weight,
          x_cols = "X",
          learner = learner,
          spline_df = threshold_spline_df,
          control = local_control,
          rkhs_plan = local_rkhs_plan,
          seed = learner_seed
        )
        upper_fit <- fit_threshold_learner(
          data = train[baseline_train, , drop = FALSE],
          y = train$Y[baseline_train],
          level = tau_train[baseline_train],
          weight = threshold_weight,
          x_cols = "X",
          learner = learner,
          spline_df = threshold_spline_df,
          control = local_control,
          rkhs_plan = local_rkhs_plan,
          seed = learner_seed
        )
        q_lower_train <- lower_fit$predict(train)
        q_upper_train <- upper_fit$predict(train)
        q_lower_validation <- lower_fit$predict(validation)
        q_upper_validation <- upper_fit$predict(validation)
        eif <- fit_eif_arm_fold(
          train = train,
          validation = validation,
          row_id_validation = validation$.row_id,
          nuisance_train = nuisance_train,
          nuisance_validation = nuisance_validation,
          envelope_train = envelope_train,
          envelope_validation = envelope_validation,
          t_value = t_value,
          Lambda0 = lookup$settings$Lambda0,
          Lambda1 = lookup$settings$Lambda1,
          q_lower_train = q_lower_train,
          q_upper_train = q_upper_train,
          q_lower_validation = q_lower_validation,
          q_upper_validation = q_upper_validation,
          x_cols = "X",
          conditional_spline_df = conditional_spline_df
        )
        train_contributions <- evaluate_threshold_bound_contributions(
          train,
          nuisance_train,
          envelope_train,
          t_value,
          tau_train,
          q_lower_train,
          q_upper_train
        )
        validation_contributions <-
          evaluate_threshold_bound_contributions(
            validation,
            nuisance_validation,
            envelope_validation,
            t_value,
            tau_validation,
            q_lower_validation,
            q_upper_validation
          )
        lower_key <- paste0("mu", t_value, "_lower")
        upper_key <- paste0("mu", t_value, "_upper")
        signal_storage[[learner]][[lower_key]][validation_index] <-
          eif$signals$lower_signal
        signal_storage[[learner]][[upper_key]][validation_index] <-
          eif$signals$upper_signal
        direct_storage[[learner]][[lower_key]][validation_index] <-
          validation_contributions$direct_lower
        direct_storage[[learner]][[upper_key]][validation_index] <-
          validation_contributions$direct_upper
        pseudo_storage[[learner]][[lower_key]][validation_index] <-
          validation_contributions$pseudo_lower
        pseudo_storage[[learner]][[upper_key]][validation_index] <-
          validation_contributions$pseudo_upper

        diagnostic_rows[[diagnostic_counter]] <- data.frame(
          learner = learner,
          rkhs_backend = lower_fit$rkhs_backend,
          fold = fold,
          treatment = t_value,
          n_train_baseline = sum(baseline_train),
          n_validation_baseline = sum(baseline_validation),
          elapsed_seconds = proc.time()[["elapsed"]] - started,
          lower_q_bias = mean(
            q_lower_validation[baseline_validation] -
              exact_q_lower_validation[baseline_validation]
          ),
          lower_q_rmse = sqrt(mean(
            (q_lower_validation[baseline_validation] -
               exact_q_lower_validation[baseline_validation])^2
          )),
          upper_q_bias = mean(
            q_upper_validation[baseline_validation] -
              exact_q_upper_validation[baseline_validation]
          ),
          upper_q_rmse = sqrt(mean(
            (q_upper_validation[baseline_validation] -
               exact_q_upper_validation[baseline_validation])^2
          )),
          lower_convergence = lower_fit$convergence,
          upper_convergence = upper_fit$convergence,
          lower_iterations = fit_iterations(lower_fit),
          upper_iterations = fit_iterations(upper_fit),
          lower_selected_spline_df = lower_fit$selected_spline_df,
          upper_selected_spline_df = upper_fit$selected_spline_df,
          lower_selected_lambda = lower_fit$selected_lambda,
          upper_selected_lambda = upper_fit$selected_lambda,
          lower_selected_bandwidth_multiplier =
            lower_fit$selected_bandwidth_multiplier,
          upper_selected_bandwidth_multiplier =
            upper_fit$selected_bandwidth_multiplier,
          lower_selected_bandwidth = lower_fit$selected_bandwidth,
          upper_selected_bandwidth = upper_fit$selected_bandwidth,
          lower_base_bandwidth = lower_fit$base_bandwidth,
          upper_base_bandwidth = upper_fit$base_bandwidth,
          lower_selected_rank = lower_fit$selected_rank,
          upper_selected_rank = upper_fit$selected_rank,
          lower_primal_residual = lower_fit$primal_residual,
          upper_primal_residual = upper_fit$primal_residual,
          lower_dual_residual = lower_fit$dual_residual,
          upper_dual_residual = upper_fit$dual_residual,
          lower_reciprocal_condition =
            lower_fit$minimum_reciprocal_condition,
          upper_reciprocal_condition =
            upper_fit$minimum_reciprocal_condition,
          lower_numerical_jitter =
            lower_fit$maximum_numerical_jitter,
          upper_numerical_jitter =
            upper_fit$maximum_numerical_jitter,
          lower_final_rho = lower_fit$final_rho,
          upper_final_rho = upper_fit$final_rho,
          lower_direct_minus_pseudo_train = mean(
            train_contributions$direct_lower -
              train_contributions$pseudo_lower
          ),
          upper_direct_minus_pseudo_train = mean(
            train_contributions$direct_upper -
              train_contributions$pseudo_upper
          ),
          lower_direct_minus_pseudo_validation = mean(
            validation_contributions$direct_lower -
              validation_contributions$pseudo_lower
          ),
          upper_direct_minus_pseudo_validation = mean(
            validation_contributions$direct_upper -
              validation_contributions$pseudo_upper
          ),
          stringsAsFactors = FALSE
        )
        diagnostic_counter <- diagnostic_counter + 1L

        for (endpoint in c("lower", "upper")) {
          fit <- if (endpoint == "lower") lower_fit else upper_fit
          if (nrow(fit$cv_table) > 0L) {
            local <- fit$cv_table
            local$learner <- learner
            local$endpoint <- endpoint
            local$fold <- fold
            local$treatment <- t_value
            cv_rows[[cv_counter]] <- local
            cv_counter <- cv_counter + 1L
          }
        }
      }
    }
  }

  list(
    signals = lapply(signal_storage, as.data.frame),
    direct = lapply(direct_storage, as.data.frame),
    pseudo = lapply(pseudo_storage, as.data.frame),
    threshold_diagnostics = do.call(rbind, diagnostic_rows),
    threshold_cv = if (length(cv_rows) > 0L) {
      all_names <- unique(unlist(lapply(cv_rows, names)))
      completed <- lapply(cv_rows, function(local) {
        missing <- setdiff(all_names, names(local))
        for (name in missing) local[[name]] <- NA
        local[, all_names, drop = FALSE]
      })
      do.call(rbind, completed)
    } else {
      data.frame()
    },
    fold_id = fold_id
  )
}

complete_ate_signals <- function(arm_signals) {
  arm_signals$ate_lower <-
    arm_signals$mu1_lower - arm_signals$mu0_upper
  arm_signals$ate_upper <-
    arm_signals$mu1_upper - arm_signals$mu0_lower
  arm_signals
}

summarize_eif_signal_version <- function(
    arm_signals,
    estimator,
    population_targets,
    confidence_level = 0.95
) {
  signals <- complete_ate_signals(arm_signals)
  endpoints <- names(population_targets)
  critical_value <- stats::qnorm(1 - (1 - confidence_level) / 2)

  endpoint_rows <- lapply(endpoints, function(endpoint) {
    signal <- signals[[endpoint]]
    estimate <- mean(signal)
    standard_error <- stats::sd(signal) / sqrt(length(signal))
    target <- population_targets[[endpoint]]
    data.frame(
      estimator = estimator,
      endpoint = endpoint,
      estimate = estimate,
      population_target = target,
      standard_error = standard_error,
      confidence_lower = estimate - critical_value * standard_error,
      confidence_upper = estimate + critical_value * standard_error,
      contains_population_endpoint =
        estimate - critical_value * standard_error <= target &&
        target <= estimate + critical_value * standard_error,
      stringsAsFactors = FALSE
    )
  })

  decomposition_rows <- lapply(
    c("ate_lower", "ate_upper"),
    function(endpoint) {
      arm1 <- if (endpoint == "ate_lower") {
        signals$mu1_lower
      } else {
        signals$mu1_upper
      }
      arm0 <- if (endpoint == "ate_lower") {
        signals$mu0_upper
      } else {
        signals$mu0_lower
      }
      n <- length(arm1)
      covariance <- stats::cov(arm1, arm0)
      data.frame(
        estimator = estimator,
        endpoint = endpoint,
        arm1_variance_over_n = stats::var(arm1) / n,
        arm0_variance_over_n = stats::var(arm0) / n,
        covariance_over_n = covariance / n,
        negative_two_covariance_over_n = -2 * covariance / n,
        ate_variance_over_n = stats::var(arm1 - arm0) / n,
        stringsAsFactors = FALSE
      )
    }
  )

  list(
    endpoints = do.call(rbind, endpoint_rows),
    variance_decomposition = do.call(rbind, decomposition_rows)
  )
}

summarize_eif_component_variance <- function(components, estimator) {
  groups <- split(
    components,
    interaction(
      components$treatment,
      components$endpoint,
      drop = TRUE
    )
  )

  rows <- lapply(groups, function(group) {
    n <- nrow(group)
    covariance <- stats::cov(group[, c(
      "M",
      "treatment_correction",
      "outcome_correction"
    )])

    data.frame(
      estimator = estimator,
      treatment = group$treatment[1],
      endpoint = group$endpoint[1],
      M_variance_over_n = covariance["M", "M"] / n,
      treatment_variance_over_n =
        covariance["treatment_correction", "treatment_correction"] / n,
      outcome_variance_over_n =
        covariance["outcome_correction", "outcome_correction"] / n,
      two_M_treatment_covariance_over_n =
        2 * covariance["M", "treatment_correction"] / n,
      two_M_outcome_covariance_over_n =
        2 * covariance["M", "outcome_correction"] / n,
      two_treatment_outcome_covariance_over_n =
        2 * covariance[
          "treatment_correction",
          "outcome_correction"
        ] / n,
      signal_variance_over_n = stats::var(group$signal) / n,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
