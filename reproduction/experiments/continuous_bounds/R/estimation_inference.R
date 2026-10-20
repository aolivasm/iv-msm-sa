binary_check_indicator <- function (value, name)
{
    observed <- sort(unique(stats::na.omit(value)))
    if (length(observed) == 0L || !all(observed %in% c(0x0p+0, 0x1p+0))) {
        stop(name, " must be coded 0/1.")
    }
}

binary_cross_entropy <- function (outcome, probability)
{
    outcome <- as.numeric(outcome)
    probability <- pmin(0x1p+0 - 0x1.5798ee2308c3ap-27, pmax(0x1.5798ee2308c3ap-27, as.numeric(probability)))
    if (length(outcome) != length(probability) || any(!is.finite(outcome)) || any(!is.finite(probability))) {
        stop("Invalid inputs to binary_cross_entropy().")
    }
    mean(-outcome * log(probability) - (0x1p+0 - outcome) * log1p(-probability))
}

binary_cross_entropy_loss <- function (outcome, probability)
{
    outcome <- as.numeric(outcome)
    probability <- pmin(0x1p+0 - 0x1.5798ee2308c3ap-27, pmax(0x1.5798ee2308c3ap-27, as.numeric(probability)))
    if (length(outcome) != length(probability) || any(!is.finite(outcome)) || any(!is.finite(probability))) {
        stop("Invalid inputs to binary_cross_entropy_loss().")
    }
    mean(-outcome * log(probability) - (0x1p+0 - outcome) * log1p(-probability))
}

binary_make_gam_formula <- function (response, data, x_cols, basis_dimension)
{
    terms <- unlist(lapply(x_cols, function(variable) {
        values <- sort(unique(stats::na.omit(data[[variable]])))
        if (length(values) <= 1L) {
            character()
        }
        else if (length(values) <= 2L) {
            sprintf("`%s`", variable)
        }
        else {
            local_k <- min(as.integer(basis_dimension), length(values) - 1L)
            local_k <- max(3L, local_k)
            sprintf("s(`%s`, k=%d, bs='tp')", variable, local_k)
        }
    }), use.names = FALSE)
    stats::as.formula(paste0("`", response, "` ~ ", if (length(terms))
        paste(terms, collapse = " + ")
    else "1"))
}

binary_make_interaction_gam_formula <- function (response, data, x_cols, basis_dimension, interaction_basis_dimension = 3L)
{
    unique_counts <- vapply(x_cols, function(variable) {
        length(unique(stats::na.omit(data[[variable]])))
    }, integer(0x1p+0))
    active <- x_cols[unique_counts > 1L]
    if (!length(active)) {
        return(stats::as.formula(sprintf("`%s` ~ 1", response)))
    }
    main_terms <- vapply(active, function(variable) {
        count <- unique_counts[[variable]]
        if (count <= 2L)
            return(sprintf("`%s`", variable))
        local_k <- min(as.integer(basis_dimension), count - 1L)
        local_k <- max(3L, local_k)
        sprintf("s(`%s`, k=%d, bs='tp')", variable, local_k)
    }, character(0x1p+0))
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
                binary_variable <- pair[if (first_count <= 2L)
                  1L
                else 2L]
                continuous_variable <- pair[if (first_count <= 2L)
                  2L
                else 1L]
                continuous_count <- unique_counts[[continuous_variable]]
                local_k <- min(as.integer(interaction_basis_dimension), continuous_count - 1L)
                local_k <- max(3L, local_k)
                return(sprintf("s(`%s`, by=`%s`, k=%d, bs='tp')", continuous_variable, binary_variable,
                  local_k))
            }
            first_k <- max(3L, min(as.integer(interaction_basis_dimension), first_count - 1L))
            second_k <- max(3L, min(as.integer(interaction_basis_dimension), second_count - 1L))
            sprintf("ti(`%s`, `%s`, k=c(%d,%d), bs=c('tp','tp'))", pair[1L], pair[2L], first_k, second_k)
        }, character(0x1p+0))
    }
    stats::as.formula(paste0("`", response, "` ~ ", paste(c(main_terms, interaction_terms), collapse = " + ")))
}

binary_make_spline_formula <- function (response, data, x_cols, spline_df)
{
    terms <- unlist(lapply(x_cols, function(variable) {
        values <- sort(unique(stats::na.omit(data[[variable]])))
        if (length(values) <= 1L) {
            character()
        }
        else if (length(values) <= 2L) {
            sprintf("`%s`", variable)
        }
        else {
            sprintf("splines::ns(`%s`, df=%d)", variable, as.integer(spline_df))
        }
    }), use.names = FALSE)
    stats::as.formula(paste0("`", response, "` ~ ", if (length(terms))
        paste(terms, collapse = " + ")
    else "1"))
}

binary_pairwise_interaction_design <- function (data)
{
    data <- as.data.frame(data)
    if (!ncol(data))
        stop("At least one predictor is required.")
    design <- as.matrix(data)
    storage.mode(design) <- "double"
    if (any(!is.finite(design))) {
        stop("The pairwise-interaction design contains non-finite values.")
    }
    if (ncol(design) == 1L) {
        squared <- matrix(design[, 1L]^0x1p+1, ncol = 1L, dimnames = list(NULL, paste0(colnames(design),
            "_squared")))
        return(cbind(design, squared))
    }
    pairs <- utils::combn(seq_len(ncol(design)), 2L)
    interactions <- vapply(seq_len(ncol(pairs)), function(index) {
        design[, pairs[1L, index]] * design[, pairs[2L, index]]
    }, numeric(nrow(design)))
    if (is.null(dim(interactions)))
        interactions <- matrix(interactions, ncol = 1L)
    colnames(interactions) <- vapply(seq_len(ncol(pairs)), function(index) {
        paste0(colnames(design)[pairs[1L, index]], "_X_", colnames(design)[pairs[2L, index]])
    }, character(0x1p+0))
    cbind(design, interactions)
}

binary_resolve_gam_basis_dimension <- function (data, x_cols, basis_dimension = NULL, minimum = 5L, maximum = 20L, growth_exponent = 0x1p+0/0x1.8p+1)
{
    if (!is.null(basis_dimension)) {
        basis_dimension <- as.integer(basis_dimension)
        if (length(basis_dimension) != 1L || is.na(basis_dimension) || basis_dimension < 3L) {
            stop("gam_basis_dimension must be NULL or one integer of at least 3.")
        }
        return(basis_dimension)
    }
    minimum <- as.integer(minimum)
    maximum <- as.integer(maximum)
    if (length(minimum) != 1L || is.na(minimum) || minimum < 3L || length(maximum) != 1L || is.na(maximum) ||
        maximum < minimum || length(growth_exponent) != 1L || !is.finite(growth_exponent) || growth_exponent <=
        0x0p+0 || growth_exponent > 0x1p+0) {
        stop("Invalid adaptive GAM basis controls.")
    }
    observations_per_covariate <- nrow(data)/max(1L, length(x_cols))
    proposed <- ceiling(observations_per_covariate^growth_exponent)
    as.integer(min(maximum, max(minimum, proposed)))
}

bind_conditional_cv_tables <- function (baseline, comparison)
{
    rows <- list()
    if (nrow(baseline) > 0L) {
        baseline$stratum <- "baseline_z0"
        rows[[length(rows) + 1L]] <- baseline
    }
    if (nrow(comparison) > 0L) {
        comparison$stratum <- "comparison_z1"
        rows[[length(rows) + 1L]] <- comparison
    }
    if (length(rows) == 0L)
        return(data.frame())
    output <- do.call(rbind, rows)
    rownames(output) <- NULL
    output
}

build_normalization_by_x_bin <- function (predictions, x_bins = 0x1.4p+3)
{
    output <- list()
    counter <- 0x1p+0
    groups <- split(predictions, interaction(predictions$mode, predictions$treatment, drop = TRUE))
    for (group in groups) {
        ranks <- rank(group$X, ties.method = "average")
        bin <- pmin(x_bins, ceiling(x_bins * ranks/length(ranks)))
        for (b in seq_len(x_bins)) {
            local <- group[bin == b, , drop = FALSE]
            if (nrow(local) == 0x0p+0) {
                next
            }
            output[[counter]] <- data.frame(mode = local$mode[0x1p+0], treatment = local$treatment[0x1p+0],
                x_bin = b, n = nrow(local), x_min = min(local$X), x_max = max(local$X), x_mean = mean(local$X),
                mean_Lbar = mean(local$Lbar), mean_Ubar = mean(local$Ubar), mean_tau = mean(local$tau),
                mean_lambda_lower = mean(local$lambda_lower), mean_lambda_upper = mean(local$lambda_upper),
                mean_gamma = mean(local$gamma_used))
            counter <- counter + 0x1p+0
        }
    }
    do.call(rbind, output)
}

check_binary <- function (x, name)
{
    values <- sort(unique(stats::na.omit(x)))
    if (!all(values %in% c(0x0p+0, 0x1p+0))) {
        stop(sprintf("%s must be coded 0/1.", name))
    }
}

check_loss <- function (y, q, level)
{
    residual <- y - q
    residual * (level - as.numeric(residual < 0x0p+0))
}

clip_positive <- function (x, floor = 0x1.0c6f7a0b5ed8dp-20, cap = Inf)
{
    x[!is.finite(x)] <- floor
    x <- pmax(floor, x)
    if (is.finite(cap)) {
        x <- pmin(cap, x)
    }
    x
}

clip_probability <- function (p, epsilon = 0x1.999999999999ap-5)
{
    pmin(0x1p+0 - epsilon, pmax(epsilon, p))
}

compute_eif_envelope_derivatives <- function (nuisance, gamma, t_value, Lambda0, Lambda1, envelope_smoothing_epsilon = 0x0p+0)
{
    if (Lambda0 <= 0x1p+0 || Lambda1 <= 0x1p+0) {
        stop(paste0("The generic EIF envelope derivative requires ", "Lambda0 > 1 and Lambda1 > 1."))
    }
    r0 <- if (t_value == 0x1p+0) {
        nuisance$pi0
    }
    else {
        0x1p+0 - nuisance$pi0
    }
    r1 <- if (t_value == 0x1p+0) {
        nuisance$pi1
    }
    else {
        0x1p+0 - nuisance$pi1
    }
    c_lower <- 0x1p+0/Lambda1
    c_upper <- Lambda1
    h_lower <- r1 + (0x1p+0 - r1) * c_lower
    h_upper <- r1 + (0x1p+0 - r1) * c_upper
    candidate_lower <- (gamma * h_lower - r0)/(0x1p+0 - r0)
    candidate_upper <- (gamma * h_upper - r0)/(0x1p+0 - r0)
    active_lower <- candidate_lower > 0x1p+0/Lambda0
    active_upper <- candidate_upper < Lambda0
    if (length(envelope_smoothing_epsilon) != 1L || !is.finite(envelope_smoothing_epsilon) || envelope_smoothing_epsilon <
        0x0p+0) {
        stop("envelope_smoothing_epsilon must be one nonnegative finite number.")
    }
    if (envelope_smoothing_epsilon > 0x0p+0) {
        smoothing_bandwidth <- envelope_smoothing_epsilon/log(0x1p+1)
        lower_multiplier <- stats::plogis((candidate_lower - 0x1p+0/Lambda0)/smoothing_bandwidth)
        upper_multiplier <- stats::plogis((Lambda0 - candidate_upper)/smoothing_bandwidth)
    }
    else {
        smoothing_bandwidth <- 0x0p+0
        lower_multiplier <- as.numeric(active_lower)
        upper_multiplier <- as.numeric(active_upper)
    }
    derivative_candidate <- function(c_value, h_value, multiplier) {
        list(gamma = multiplier * h_value/(0x1p+0 - r0), r0 = multiplier * (gamma * h_value - 0x1p+0)/(0x1p+0 -
            r0)^0x1p+1, r1 = multiplier * gamma * (0x1p+0 - c_value)/(0x1p+0 - r0))
    }
    list(r0 = r0, r1 = r1, candidate_lower = candidate_lower, candidate_upper = candidate_upper, active_lower = active_lower,
        active_upper = active_upper, lower_derivative_weight = lower_multiplier, upper_derivative_weight = upper_multiplier,
        envelope_smoothing_epsilon = envelope_smoothing_epsilon, envelope_smoothing_bandwidth = smoothing_bandwidth,
        lower = derivative_candidate(c_value = c_lower, h_value = h_lower, multiplier = lower_multiplier),
        upper = derivative_candidate(c_value = c_upper, h_value = h_upper, multiplier = upper_multiplier))
}

compute_eif_row_primitives <- function (data, nuisance, envelope, t_value, Lambda0, Lambda1, q_lower, q_upper)
{
    required_lengths <- c(nrow(data), length(nuisance$pi0), length(nuisance$pi1), length(envelope$L),
        length(envelope$U), length(envelope$gamma_used), length(q_lower), length(q_upper))
    if (length(unique(required_lengths)) != 1L) {
        stop("EIF row inputs must have a common length.")
    }
    if (any(!is.finite(q_lower)) || any(!is.finite(q_upper))) {
        stop("EIF threshold predictions must be finite.")
    }
    derivative <- compute_eif_envelope_derivatives(nuisance = nuisance, gamma = envelope$gamma_used,
        t_value = t_value, Lambda0 = Lambda0, Lambda1 = Lambda1, envelope_smoothing_epsilon = if (is.null(envelope$envelope_smoothing_epsilon)) {
            0x0p+0
        }
        else {
            envelope$envelope_smoothing_epsilon
        })
    make_endpoint <- function(q, direction = c("lower", "upper")) {
        direction <- match.arg(direction)
        allocate_upper <- if (direction == "lower") {
            data$Y <= q
        }
        else {
            data$Y > q
        }
        allocate_upper <- as.numeric(allocate_upper)
        allocate_lower <- 0x1p+0 - allocate_upper
        lambda <- allocate_upper * envelope$U + allocate_lower * envelope$L
        envelope_derivative <- function(name) {
            (data$Y - q) * (allocate_upper * derivative$upper[[name]] + allocate_lower * derivative$lower[[name]])
        }
        D_gamma <- envelope_derivative("gamma")
        D_r0 <- envelope_derivative("r0")
        D_r1 <- envelope_derivative("r1")
        H0 <- derivative$r0 * data$Y + (0x1p+0 - derivative$r0) * (lambda * (data$Y - q) - envelope$gamma_used *
            D_gamma)
        H1 <- (0x1p+0 - derivative$r0) * D_gamma
        data.frame(q = q, lambda = lambda, D_gamma = D_gamma, D_r0 = D_r0, D_r1 = D_r1, H0 = H0, H1 = H1)
    }
    list(r0 = derivative$r0, r1 = derivative$r1, active_lower = derivative$active_lower, active_upper = derivative$active_upper,
        lower_derivative_weight = derivative$lower_derivative_weight, upper_derivative_weight = derivative$upper_derivative_weight,
        lower = make_endpoint(q_lower, "lower"), upper = make_endpoint(q_upper, "upper"))
}

compute_fixed_q_linked_components <- function (m0, Hbar0, Hbar1, q, r0, mean_D_r0)
{
    input_lengths <- lengths(list(m0 = m0, Hbar0 = Hbar0, Hbar1 = Hbar1, q = q, r0 = r0, mean_D_r0 = mean_D_r0))
    if (length(unique(input_lengths)) != 1L || input_lengths[1L] == 0L) {
        stop("Fixed-q linked components must have a common positive length.")
    }
    inputs <- c(m0, Hbar0, Hbar1, q, r0, mean_D_r0)
    if (any(!is.finite(inputs)) || any(r0 <= 0x0p+0) || any(r0 >= 0x1p+0)) {
        stop("Fixed-q linked component inputs must be finite with 0 < r0 < 1.")
    }
    M <- Hbar0 + Hbar1 + (0x1p+0 - r0) * q
    theta_fixed_q <- (M - r0 * m0)/(0x1p+0 - r0)
    C0 <- m0 - theta_fixed_q + (0x1p+0 - r0) * mean_D_r0
    list(M = M, theta_fixed_q = theta_fixed_q, C0 = C0)
}

compute_ipw_identities_simple <- function (data)
{
    c(t1_z0_estimated = mean((data$Z == 0x0p+0) * data$T_obs/(data$e_hat * data$pi0_hat)), t1_z0_oracle = mean((data$Z ==
        0x0p+0) * data$T_obs/(data$e_true * data$pi0_true)), t0_z0_estimated = mean((data$Z == 0x0p+0) *
        (0x1p+0 - data$T_obs)/(data$e_hat * (0x1p+0 - data$pi0_hat))), t0_z0_oracle = mean((data$Z ==
        0x0p+0) * (0x1p+0 - data$T_obs)/(data$e_true * (0x1p+0 - data$pi0_true))), t1_z1_estimated = mean(data$Z *
        data$T_obs/((0x1p+0 - data$e_hat) * data$pi1_hat)), t1_z1_oracle = mean(data$Z * data$T_obs/((0x1p+0 -
        data$e_true) * data$pi1_true)), t0_z1_estimated = mean(data$Z * (0x1p+0 - data$T_obs)/((0x1p+0 -
        data$e_hat) * (0x1p+0 - data$pi1_hat))), t0_z1_oracle = mean(data$Z * (0x1p+0 - data$T_obs)/((0x1p+0 -
        data$e_true) * (0x1p+0 - data$pi1_true))))
}

compute_oracle_quantities <- function (data, confounder_strength = 0x1p+0, treatment_intercept = 0x0p+0, instrument_strength = 0x1p+0,
    treatment_effect = 0x1p+1, n_nodes = 0x1.4p+6, batch_size = 0x1.f4p+9)
{
    required <- c("Y", "X")
    missing_columns <- setdiff(required, names(data))
    if (length(missing_columns) > 0x0p+0) {
        stop("Missing columns: ", paste(missing_columns, collapse = ", "))
    }
    rule <- gauss_legendre_rule(n_nodes = n_nodes, lower = -0x1p+1, upper = 0x1p+1)
    density_rule <- gauss_legendre_rule(n_nodes = n_nodes, lower = -0x1p+0, upper = 0x1p+0)
    u <- rule$nodes
    integration_weights <- rule$weights * d_u(u)
    n <- nrow(data)
    pi0_true <- rep(NA_real_, n)
    pi1_true <- rep(NA_real_, n)
    e_true <- rep(NA_real_, n)
    gamma0_oracle <- rep(NA_real_, n)
    gamma1_oracle <- rep(NA_real_, n)
    starts <- seq.int(0x1p+0, n, by = batch_size)
    for (start in starts) {
        finish <- min(start + batch_size - 0x1p+0, n)
        index <- start:finish
        x <- data$X[index]
        y <- data$Y[index]
        n_batch <- length(index)
        x_matrix <- matrix(x, nrow = n_batch, ncol = n_nodes)
        u_matrix <- matrix(u, nrow = n_batch, ncol = n_nodes, byrow = TRUE)
        eta_z0 <- treatment_intercept + x_matrix + confounder_strength * u_matrix
        eta_z1 <- eta_z0 + instrument_strength
        p_t1_z0 <- stats::plogis(eta_z0)
        p_t1_z1 <- stats::plogis(eta_z1)
        weighted_integral <- function(matrix_value) {
            as.vector(matrix_value %*% integration_weights)
        }
        pi0_batch <- weighted_integral(p_t1_z0)
        pi1_batch <- weighted_integral(p_t1_z1)
        pi0_true[index] <- pi0_batch
        pi1_true[index] <- pi1_batch
        e_true[index] <- stats::plogis(-0x1p-1 * x)
        for (t_value in c(0x0p+0, 0x1p+0)) {
            outcome_location <- treatment_effect * t_value + 0x1.8p+0 * x
            density_u_lower <- pmax(-0x1p+1, (y - outcome_location - 0x1p+1)/0x1.8p+0)
            density_u_upper <- pmin(0x1p+1, (y - outcome_location + 0x1p+1)/0x1.8p+0)
            density_half_width <- pmax(0x0p+0, (density_u_upper - density_u_lower)/0x1p+1)
            density_midpoint <- (density_u_upper + density_u_lower)/0x1p+1
            density_u_matrix <- matrix(density_midpoint, nrow = n_batch, ncol = n_nodes) + outer(density_half_width,
                density_rule$nodes)
            density_quadrature_weights <- outer(density_half_width, density_rule$weights) * d_u(density_u_matrix)
            density_x_matrix <- matrix(x, nrow = n_batch, ncol = n_nodes)
            density_y_matrix <- matrix(y, nrow = n_batch, ncol = n_nodes)
            mean_y <- treatment_effect * t_value + 0x1.8p+0 * density_x_matrix + 0x1.8p+0 * density_u_matrix
            density_y_given_xu <- dtnormtc(density_y_matrix, mu = mean_y, s = 0x1p+0)
            density_eta_z0 <- treatment_intercept + density_x_matrix + confounder_strength * density_u_matrix
            density_eta_z1 <- density_eta_z0 + instrument_strength
            density_p_t1_z0 <- stats::plogis(density_eta_z0)
            density_p_t1_z1 <- stats::plogis(density_eta_z1)
            treatment_weight_z0 <- if (t_value == 0x1p+0) {
                density_p_t1_z0
            }
            else {
                0x1p+0 - density_p_t1_z0
            }
            treatment_weight_z1 <- if (t_value == 0x1p+0) {
                density_p_t1_z1
            }
            else {
                0x1p+0 - density_p_t1_z1
            }
            joint_y_t_z0 <- rowSums(density_y_given_xu * treatment_weight_z0 * density_quadrature_weights)
            joint_y_t_z1 <- rowSums(density_y_given_xu * treatment_weight_z1 * density_quadrature_weights)
            probability_t_z0 <- if (t_value == 0x1p+0) {
                pi0_batch
            }
            else {
                0x1p+0 - pi0_batch
            }
            probability_t_z1 <- if (t_value == 0x1p+0) {
                pi1_batch
            }
            else {
                0x1p+0 - pi1_batch
            }
            density_y_t_z0 <- joint_y_t_z0/pmax(probability_t_z0, 0x1.19799812dea11p-40)
            density_y_t_z1 <- joint_y_t_z1/pmax(probability_t_z1, 0x1.19799812dea11p-40)
            gamma_value <- ifelse(density_y_t_z0 > 0x1.56e1fc2f8f355p-997, density_y_t_z1/density_y_t_z0,
                0x1p+0)
            if (t_value == 0x1p+0) {
                gamma1_oracle[index] <- gamma_value
            }
            else {
                gamma0_oracle[index] <- gamma_value
            }
        }
    }
    data.frame(pi0_true = pi0_true, pi1_true = pi1_true, e_true = e_true, gamma0_oracle = gamma0_oracle,
        gamma1_oracle = gamma1_oracle)
}

compute_projected_envelopes <- function (nuisance, gamma, t_value, Lambda0, Lambda1, envelope_smoothing_epsilon = 0x0p+0, gamma_compatibility_action = c("stop",
    "project"), gamma_projection_margin = 0x1.a36e2eb1c432dp-14)
{
    gamma_compatibility_action <- match.arg(gamma_compatibility_action)
    if (Lambda0 < 0x1p+0 || Lambda1 < 0x1p+0) {
        stop("Lambda0 and Lambda1 must both be at least one.")
    }
    if (length(envelope_smoothing_epsilon) != 1L || !is.finite(envelope_smoothing_epsilon) || envelope_smoothing_epsilon <
        0x0p+0) {
        stop("envelope_smoothing_epsilon must be one nonnegative finite number.")
    }
    if (envelope_smoothing_epsilon > 0x0p+0 && envelope_smoothing_epsilon >= 0x1p+0/Lambda0) {
        stop(paste0("envelope_smoothing_epsilon must be smaller than 1/Lambda0 so ", "the conservative smoothed lower envelope remains positive."))
    }
    if (Lambda0 == 0x1p+0 && Lambda1 == 0x1p+0) {
        interval <- compute_compatible_gamma_interval(nuisance = nuisance, t_value = t_value, Lambda0 = Lambda0,
            Lambda1 = Lambda1)
        incompatible_original <- !is.finite(gamma) | abs(gamma - 0x1p+0) > 0x1.4p+3 * .Machine$double.eps
        return(list(A = rep(NA_real_, length(gamma)), B = rep(NA_real_, length(gamma)), L = rep(0x1p+0,
            length(gamma)), U = rep(0x1p+0, length(gamma)), D = rep(0x0p+0, length(gamma)), hard_L = rep(0x1p+0,
            length(gamma)), hard_U = rep(0x1p+0, length(gamma)), candidate_lower = rep(0x1p+0, length(gamma)),
            candidate_upper = rep(0x1p+0, length(gamma)), lower_smoothing_gap = rep(0x0p+0, length(gamma)),
            upper_smoothing_gap = rep(0x0p+0, length(gamma)), envelope_smoothing_epsilon = 0x0p+0, envelope_smoothing_bandwidth = 0x0p+0,
            gamma_original = gamma, gamma_used = gamma, gamma_compatible_lower = interval$lower, gamma_compatible_upper = interval$upper,
            gamma_incompatible_original = incompatible_original, gamma_projected = rep(FALSE, length(gamma))))
    }
    interval <- compute_compatible_gamma_interval(nuisance = nuisance, t_value = t_value, Lambda0 = Lambda0,
        Lambda1 = Lambda1)
    incompatible_original <- (gamma < interval$lower | gamma > interval$upper | !is.finite(gamma))
    if (gamma_compatibility_action == "project") {
        gamma_projection <- project_gamma_to_compatible_interval(gamma = gamma, nuisance = nuisance,
            t_value = t_value, Lambda0 = Lambda0, Lambda1 = Lambda1, margin_fraction = gamma_projection_margin)
        gamma_used <- gamma_projection$gamma
        gamma_projected <- gamma_projection$projected
    }
    else {
        gamma_used <- gamma
        gamma_projected <- rep(FALSE, length(gamma))
    }
    pi0 <- nuisance$pi0
    pi1 <- nuisance$pi1
    if (t_value == 0x1p+0) {
        A <- ((0x1p+0 - pi1)/(0x1p+0 - pi0)) * gamma_used
        B <- (pi1 * gamma_used - pi0)/(0x1p+0 - pi0)
    }
    else {
        A <- (pi1/pi0) * gamma_used
        B <- ((0x1p+0 - pi1) * gamma_used - (0x1p+0 - pi0))/pi0
    }
    candidate_lower <- A/Lambda1 + B
    candidate_upper <- A * Lambda1 + B
    hard_L <- pmax(0x1p+0/Lambda0, candidate_lower)
    hard_U <- pmin(Lambda0, candidate_upper)
    smoothing_bandwidth <- if (envelope_smoothing_epsilon > 0x0p+0) {
        envelope_smoothing_epsilon/log(0x1p+1)
    }
    else {
        0x0p+0
    }
    if (envelope_smoothing_epsilon > 0x0p+0) {
        smooth_lower <- function(first, second, bandwidth) {
            maximum <- pmax(first, second)
            value <- maximum + bandwidth * log((exp((first - maximum)/bandwidth) + exp((second - maximum)/bandwidth))/0x1p+1)
            value[is.infinite(maximum) & maximum > 0x0p+0] <- Inf
            value
        }
        smooth_upper <- function(first, second, bandwidth) {
            minimum <- pmin(first, second)
            value <- minimum - bandwidth * log((exp(-(first - minimum)/bandwidth) + exp(-(second - minimum)/bandwidth))/0x1p+1)
            value[is.infinite(minimum) & minimum < 0x0p+0] <- -Inf
            value
        }
        L <- smooth_lower(0x1p+0/Lambda0, candidate_lower, smoothing_bandwidth)
        U <- smooth_upper(Lambda0, candidate_upper, smoothing_bandwidth)
    }
    else {
        L <- hard_L
        U <- hard_U
    }
    list(A = A, B = B, L = L, U = U, D = U - L, hard_L = hard_L, hard_U = hard_U, candidate_lower = candidate_lower,
        candidate_upper = candidate_upper, lower_smoothing_gap = hard_L - L, upper_smoothing_gap = U -
            hard_U, envelope_smoothing_epsilon = envelope_smoothing_epsilon, envelope_smoothing_bandwidth = smoothing_bandwidth,
        gamma_original = gamma, gamma_used = gamma_used, gamma_compatible_lower = interval$lower, gamma_compatible_upper = interval$upper,
        gamma_incompatible_original = incompatible_original, gamma_projected = gamma_projected)
}

compute_target_odds <- function (nuisance, t_value)
{
    if (t_value == 0x1p+0) {
        (0x1p+0 - nuisance$pi0)/nuisance$pi0
    }
    else {
        nuisance$pi0/(0x1p+0 - nuisance$pi0)
    }
}

d_u <- function (u)
{
    dtnormtc(u, mu = 0x0p+0, s = 0x1p+0)
}

default_eif_conditional_control <- function ()
{
    list(spline_df_grid = c(3L, 5L, 7L, 9L), cv_folds = 3L, selection_rule = "one_se", gam_k = 6L, gam_method = "REML",
        ranger_num_trees = 500L, ranger_min_node_size = NULL, library_candidates = c("linear", "gam",
            "ranger"), library_selection_rule = "one_se")
}

dtnormtc <- function (x, mu = 0x0p+0, s = 0x1p+0)
{
    standardized <- (x - mu)/s
    normalizing_constant <- stats::pnorm(0x1p+1) - stats::pnorm(-0x1p+1)
    ifelse(standardized >= -0x1p+1 & standardized <= 0x1p+1, stats::dnorm(standardized)/(s * normalizing_constant),
        0x0p+0)
}

estimate_sharp_bounds_one_step <- function (data, ..., compute_oracle_layers = FALSE)
{
    stopifnot(isFALSE(compute_oracle_layers))
    if (length(compute_oracle_layers) != 1L || is.na(compute_oracle_layers)) {
        stop("compute_oracle_layers must be TRUE or FALSE.")
    }
    arguments <- list(...)
    threshold_learner <- "external"
    if ("threshold_learner" %in% names(arguments)) {
        threshold_learner <- match.arg(arguments$threshold_learner, "external")
        arguments$threshold_learner <- NULL
    }
    reserved <- intersect(names(arguments), c("compute_oracle_eif", "compute_estimated_eif", "one_step_modes",
        "estimation_modes"))
    if (length(reserved) > 0L) {
        stop("Use compute_oracle_layers rather than supplying ", paste(reserved, collapse = ", "), ".")
    }
    result <- do.call(estimate_sharp_bounds_plugin, c(list(data = data, threshold_learner = threshold_learner,
        estimation_modes = if (isTRUE(compute_oracle_layers)) {
            c("oracle_all", "oracle_gamma", "rkhs")
        } else {
            "rkhs"
        }, one_step_modes = if (isTRUE(compute_oracle_layers)) {
            c("oracle_all", "oracle_gamma", "rkhs")
        } else {
            "rkhs"
        }), arguments))
    result$one_step_configuration <- list(threshold_learner = threshold_learner, fully_estimated_mode = "rkhs",
        nuisance_learner = if ("nuisance_learner" %in% names(arguments)) {
            match.arg(arguments$nuisance_learner, c("gam", "ranger", "glm_spline", "library_cv", "super_learner",
                "integrated_parametric"))
        } else {
            "gam"
        }, nuisance_gam_basis_selection = if ("nuisance_learner" %in% names(arguments) && !identical(arguments$nuisance_learner,
            "gam")) {
            NA_character_
        } else if ("nuisance_gam_k" %in% names(arguments) && !is.null(arguments$nuisance_gam_k)) {
            "fixed"
        } else {
            "adaptive_training_stratum_size"
        }, nuisance_gam_min_k = if ("nuisance_gam_min_k" %in% names(arguments)) arguments$nuisance_gam_min_k else 5L,
        nuisance_gam_max_k = if ("nuisance_gam_max_k" %in% names(arguments)) arguments$nuisance_gam_max_k else 20L,
        nuisance_gam_growth_exponent = if ("nuisance_gam_growth_exponent" %in% names(arguments)) arguments$nuisance_gam_growth_exponent else 0x1p+0/0x1.8p+1,
        gamma_estimator = if ("gamma_fit_provider_factory" %in% names(arguments) && !is.null(arguments$gamma_fit_provider_factory)) {
            unique(result$selected_lambdas$gamma_estimator)[1L]
        } else {
            "direct_rkhs_riesz"
        }, gamma_rkhs_backend = if ("gamma_rkhs_backend" %in% names(arguments)) {
            match.arg(arguments$gamma_rkhs_backend, c("exact", "nystrom"))
        } else {
            "exact"
        }, gamma_rkhs_kernel = if ("gamma_rkhs_kernel" %in% names(arguments)) {
            match.arg(arguments$gamma_rkhs_kernel, c("gaussian", "sobolev_matern32", "sobolev_matern52"))
        } else {
            "gaussian"
        }, eif_conditional_learner = if ("eif_conditional_learner" %in% names(arguments)) {
            match.arg(arguments$eif_conditional_learner, c("spline_fixed", "spline_cv", "gam", "ranger",
                "library_cv"))
        } else {
            "spline_fixed"
        }, eif_conditional_spline_df = if ("eif_conditional_spline_df" %in% names(arguments)) {
            arguments$eif_conditional_spline_df
        } else if ("conditional_mean_spline_df" %in% names(arguments)) {
            arguments$conditional_mean_spline_df
        } else {
            5L
        }, custom_eif_instrument_probability = "eif_instrument_probability_transform" %in% names(arguments) &&
            is.function(arguments$eif_instrument_probability_transform), envelope_smoothing_epsilon = if ("envelope_smoothing_epsilon" %in%
            names(arguments)) {
            arguments$envelope_smoothing_epsilon
        } else {
            0x0p+0
        }, oracle_layers_computed = isTRUE(compute_oracle_layers))
    result
}

estimate_sharp_bounds_plugin <- function (data, Lambda0_T1 = 0x1.4p+2, Lambda1_T1 = 0x1.4p+2, Lambda0_T0 = 0x1.4p+2, Lambda1_T0 = 0x1.4p+2,
    x_cols = "X", confounder_strength = 0x1p+0, treatment_intercept = 0x0p+0, instrument_strength = 0x1p+0,
    treatment_effect = 0x1p+1, K = 0x1.4p+2, fold_id = NULL, positivity_epsilon = 0x1.999999999999ap-5,
    nuisance_learner = c("gam", "ranger", "glm_spline", "library_cv", "super_learner", "integrated_parametric"),
    treatment_spline_df = 0x1.4p+2, nuisance_gam_k = NULL, nuisance_gam_min_k = 5L, nuisance_gam_max_k = 20L,
    nuisance_gam_growth_exponent = 0x1p+0/0x1.8p+1, nuisance_gam_method = "REML", nuisance_library_candidates = c("glm_spline",
        "gam", "ranger"), nuisance_library_cv_folds = 3L, nuisance_library_selection_rule = c("one_se",
        "minimum"), nuisance_library_ranger_num_trees = 500L, nuisance_super_learner_candidates = NULL,
    nuisance_super_learner_cv_folds = 3L, nuisance_super_learner_ranger_num_trees = 500L, nuisance_super_learner_hal_max_degree = 2L,
    nuisance_integrated_nodes = 40L, nuisance_integrated_maxit = 500L, gamma_floor = 0x1.0c6f7a0b5ed8ep-20,
    gamma_cap = 0x1.9p+6, oracle_nodes = 0x1.4p+6, oracle_batch_size = 0x1.f4p+9, rkhs_lambda_by_t = NULL,
    rkhs_lambda_grid = default_gamma_lambda_grid, rkhs_bandwidth_multipliers = 0x1p+0/sqrt(0x1p+1), rkhs_cv_folds = 0x1.8p+1,
    rkhs_cv_repeats = 2L, rkhs_selection_rule = c("minimum", "one_se"), gamma_rkhs_backend = c("exact",
        "nystrom"), gamma_rkhs_kernel = c("gaussian", "sobolev_matern32", "sobolev_matern52"), gamma_nystrom_rank = NULL,
    gamma_nystrom_eigen_tolerance = 0x1.b7cdfd9d7bdbbp-34, gamma_nystrom_jitter = 0x1.b7cdfd9d7bdbbp-34,
    bandwidth_sample = 0x1.f4p+9, conditional_mean_spline_df = 0x1.4p+2, threshold_spline_df = 0x1.4p+2,
    threshold_learner = "external", threshold_control = NULL, envelope_smoothing_epsilon = 0x0p+0, estimated_compatibility_action = c("project",
        "stop"), compatibility_tolerance = 0x1.0c6f7a0b5ed8ep-20, gamma_projection_margin = 0x1.a36e2eb1c432dp-14,
    max_gamma_projection_rate = 0x1.47ae147ae147bp-6, tau_epsilon = 0x1.a36e2eb1c432dp-14, threshold_maxit = 0x1.f4p+11,
    compute_oracle_eif = TRUE, compute_estimated_eif = FALSE, one_step_modes = NULL, estimation_modes = c("oracle_all",
        "oracle_gamma", "rkhs"), eif_conditional_spline_df = conditional_mean_spline_df, eif_conditional_learner = c("spline_fixed",
        "spline_cv", "gam", "ranger", "library_cv"), eif_conditional_control = NULL, eif_conditional_provider_factory = NULL,
    eif_instrument_probability_transform = NULL, nuisance_fit_provider_factory = NULL, gamma_fit_provider_factory = NULL,
    gamma_fit_transform_factory = NULL, bound_component_provider_factory = NULL, fixed_threshold_target_evaluator = NULL,
    x_bins = 0x1.4p+3, seed = 0x1p+0)
{
    threshold_learner <- match.arg(threshold_learner)
    nuisance_learner <- match.arg(nuisance_learner)
    nuisance_library_selection_rule <- match.arg(nuisance_library_selection_rule)
    eif_conditional_learner <- match.arg(eif_conditional_learner)
    if (!is.null(eif_conditional_provider_factory) && !is.function(eif_conditional_provider_factory)) {
        stop("eif_conditional_provider_factory must be NULL or a function.")
    }
    if (!is.null(eif_instrument_probability_transform) && !is.function(eif_instrument_probability_transform)) {
        stop("eif_instrument_probability_transform must be NULL or a function.")
    }
    if (!is.null(nuisance_fit_provider_factory) && !is.function(nuisance_fit_provider_factory)) {
        stop("nuisance_fit_provider_factory must be NULL or a function.")
    }
    if (!is.null(gamma_fit_provider_factory) && !is.function(gamma_fit_provider_factory)) {
        stop("gamma_fit_provider_factory must be NULL or a function.")
    }
    if (!is.null(gamma_fit_transform_factory) && !is.function(gamma_fit_transform_factory)) {
        stop("gamma_fit_transform_factory must be NULL or a function.")
    }
    if (!is.null(bound_component_provider_factory) && !is.function(bound_component_provider_factory)) {
        stop("bound_component_provider_factory must be NULL or a function.")
    }
    if (!is.null(fixed_threshold_target_evaluator) && !is.function(fixed_threshold_target_evaluator)) {
        stop("fixed_threshold_target_evaluator must be NULL or a function.")
    }
    rkhs_selection_rule <- match.arg(rkhs_selection_rule)
    gamma_rkhs_backend <- match.arg(gamma_rkhs_backend)
    gamma_rkhs_kernel <- match.arg(gamma_rkhs_kernel)
    estimated_compatibility_action <- match.arg(estimated_compatibility_action)
    if (length(envelope_smoothing_epsilon) != 1L || !is.finite(envelope_smoothing_epsilon) || envelope_smoothing_epsilon <
        0x0p+0) {
        stop("envelope_smoothing_epsilon must be one nonnegative finite number.")
    }
    required <- unique(c("Y", "T_obs", "Z", x_cols))
    missing_columns <- setdiff(required, names(data))
    if (length(missing_columns) > 0x0p+0) {
        stop("Missing columns: ", paste(missing_columns, collapse = ", "))
    }
    analysis_data <- data[stats::complete.cases(data[, required, drop = FALSE]), , drop = FALSE]
    analysis_data <- as.data.frame(analysis_data)
    analysis_data$.row_id <- seq_len(nrow(analysis_data))
    oracle_columns <- c("e_true", "pi0_true", "pi1_true", "gamma0_oracle", "gamma1_oracle")
    if (any(c("oracle_all", "oracle_gamma") %in% estimation_modes) && !all(oracle_columns %in% names(analysis_data))) {
        oracle <- compute_oracle_quantities(data = analysis_data, confounder_strength = confounder_strength,
            treatment_intercept = treatment_intercept, instrument_strength = instrument_strength, treatment_effect = treatment_effect,
            n_nodes = oracle_nodes, batch_size = oracle_batch_size)
        for (name in names(oracle)) {
            analysis_data[[name]] <- oracle[[name]]
        }
    }
    n <- nrow(analysis_data)
    if (is.null(fold_id)) {
        fold_id <- make_tz_stratified_folds(data = analysis_data, K = K, seed = seed)
    }
    else {
        if (length(fold_id) != n) {
            stop("fold_id must have one value per analysis row.")
        }
        K <- max(fold_id)
    }
    analysis_data$fold_id <- fold_id
    available_modes <- c("oracle_all", "oracle_gamma", "rkhs")
    modes <- match.arg(estimation_modes, choices = available_modes, several.ok = TRUE)
    modes <- unique(as.character(modes))
    eif_modes <- if (is.null(one_step_modes)) {
        intersect(c(if (isTRUE(compute_oracle_eif)) "oracle_all", if (isTRUE(compute_estimated_eif)) c("oracle_gamma",
            "rkhs")), modes)
    }
    else {
        invalid_eif_modes <- setdiff(one_step_modes, available_modes)
        if (length(invalid_eif_modes) > 0L) {
            stop("Unknown one_step_modes: ", paste(invalid_eif_modes, collapse = ", "))
        }
        missing_estimation_modes <- setdiff(one_step_modes, modes)
        if (length(missing_estimation_modes) > 0L) {
            stop("Every one_step_modes entry must also be included in estimation_modes: ", paste(missing_estimation_modes,
                collapse = ", "))
        }
        unique(as.character(one_step_modes))
    }
    contribution_storage <- list()
    for (mode in modes) {
        contribution_storage[[mode]] <- list()
        for (t_value in c(0x0p+0, 0x1p+0)) {
            contribution_storage[[mode]][[as.character(t_value)]] <- list(lower = rep(NA_real_, n), upper = rep(NA_real_,
                n), lower_pseudo = rep(NA_real_, n), upper_pseudo = rep(NA_real_, n), lambda_one = rep(NA_real_,
                n))
        }
    }
    e_hat <- pi0_hat <- pi1_hat <- rep(NA_real_, n)
    gamma0_rkhs <- gamma1_rkhs <- rep(NA_real_, n)
    gamma_accuracy_rows <- list()
    gamma_accuracy_counter <- 0x1p+0
    fold_diagnostics <- list()
    diagnostic_counter <- 0x1p+0
    prediction_tables <- list()
    prediction_counter <- 0x1p+0
    selected_lambdas <- list()
    selected_lambda_counter <- 0x1p+0
    gamma_fit_diagnostics <- list()
    gamma_fit_diagnostic_counter <- 1L
    gamma_transform_diagnostics <- list()
    gamma_transform_diagnostic_counter <- 1L
    cv_tables <- list()
    cv_counter <- 0x1p+0
    threshold_cv_tables <- list()
    threshold_cv_counter <- 1L
    fixed_threshold_target_evaluations <- list()
    fixed_threshold_target_evaluation_counter <- 1L
    nuisance_fit_diagnostics <- list()
    nuisance_fit_diagnostic_counter <- 1L
    eif_storage <- setNames(lapply(eif_modes, function(mode) {
        list(`0` = list(lower = rep(NA_real_, n), upper = rep(NA_real_, n)), `1` = list(lower = rep(NA_real_,
            n), upper = rep(NA_real_, n)))
    }), eif_modes)
    eif_predictions <- list()
    eif_prediction_counter <- 1L
    eif_diagnostics <- list()
    eif_diagnostic_counter <- 1L
    eif_conditional_cv_tables <- list()
    eif_conditional_cv_counter <- 1L
    for (fold in seq_len(K)) {
        validation_index <- which(fold_id == fold)
        training_index <- which(fold_id != fold)
        train <- analysis_data[training_index, , drop = FALSE]
        validation <- analysis_data[validation_index, , drop = FALSE]
        nuisance_fit <- if (is.null(nuisance_fit_provider_factory)) {
            fit_nuisance_models(train = train, x_cols = x_cols, positivity_epsilon = positivity_epsilon,
                learner = nuisance_learner, treatment_spline_df = treatment_spline_df, gam_basis_dimension = nuisance_gam_k,
                gam_minimum_basis_dimension = nuisance_gam_min_k, gam_maximum_basis_dimension = nuisance_gam_max_k,
                gam_basis_growth_exponent = nuisance_gam_growth_exponent, gam_method = nuisance_gam_method,
                library_candidates = nuisance_library_candidates, library_cv_folds = nuisance_library_cv_folds,
                library_selection_rule = nuisance_library_selection_rule, library_ranger_num_trees = nuisance_library_ranger_num_trees,
                library_seed = seed + 90000L * fold, super_learner_candidates = nuisance_super_learner_candidates,
                super_learner_cv_folds = nuisance_super_learner_cv_folds, super_learner_ranger_num_trees = nuisance_super_learner_ranger_num_trees,
                super_learner_hal_max_degree = nuisance_super_learner_hal_max_degree, integrated_confounder_strength = confounder_strength,
                integrated_quadrature_nodes = nuisance_integrated_nodes, integrated_maxit = nuisance_integrated_maxit)
        }
        else {
            nuisance_provider_arguments <- list(fold = fold, train = train, x_cols = x_cols, seed = seed +
                90000L * fold)
            nuisance_provider_formals <- names(formals(nuisance_fit_provider_factory))
            if ("..." %in% nuisance_provider_formals) {
                nuisance_provider_arguments$positivity_epsilon <- positivity_epsilon
                nuisance_provider_arguments$learner <- nuisance_learner
            }
            do.call(nuisance_fit_provider_factory, nuisance_provider_arguments)
        }
        if (!is.list(nuisance_fit) || !is.function(nuisance_fit$predict) || !is.data.frame(nuisance_fit$diagnostics)) {
            stop(paste0("The nuisance fit provider must return a list containing a ", "predict function and a diagnostics data frame."))
        }
        local_nuisance_diagnostics <- nuisance_fit$diagnostics
        local_nuisance_diagnostics$fold <- fold
        if (!"n_training" %in% names(local_nuisance_diagnostics)) {
            local_nuisance_diagnostics$n_training <- c(nrow(train), sum(train$Z == 0x0p+0), sum(train$Z ==
                0x1p+0))
        }
        nuisance_fit_diagnostics[[nuisance_fit_diagnostic_counter]] <- local_nuisance_diagnostics
        nuisance_fit_diagnostic_counter <- nuisance_fit_diagnostic_counter + 1L
        nuisance_train_estimated <- nuisance_fit$predict(train)
        nuisance_validation_estimated <- nuisance_fit$predict(validation)
        nuisance_train_oracle <- list(e = train$e_true, pi0 = train$pi0_true, pi1 = train$pi1_true)
        nuisance_validation_oracle <- list(e = validation$e_true, pi0 = validation$pi0_true, pi1 = validation$pi1_true)
        e_hat[validation_index] <- nuisance_validation_estimated$e
        pi0_hat[validation_index] <- nuisance_validation_estimated$pi0
        pi1_hat[validation_index] <- nuisance_validation_estimated$pi1
        for (t_value in c(0x0p+0, 0x1p+0)) {
            lambda_key <- as.character(t_value)
            if (!is.null(gamma_fit_provider_factory)) {
                gamma_provider_arguments <- list(fold = fold, treatment = t_value, train = train, nuisance_train = nuisance_train_estimated,
                  x_cols = x_cols, seed = seed + 0x1.86ap+16 * fold + 0x1.9p+6 * t_value)
                gamma_provider_formals <- names(formals(gamma_fit_provider_factory))
                if ("nuisance_predict" %in% gamma_provider_formals || "..." %in% gamma_provider_formals) {
                  gamma_provider_arguments$nuisance_predict <- nuisance_fit$predict
                }
                rkhs_fit <- do.call(gamma_fit_provider_factory, gamma_provider_arguments)
                required_gamma_fit_fields <- c("predict", "selected_lambda", "selected_bandwidth_multiplier",
                  "base_bandwidth", "selected_bandwidth", "normalization_constant", "training_negative_rate",
                  "rkhs_backend", "nystrom_rank_requested", "nystrom_landmark_count", "nystrom_effective_rank",
                  "nystrom_min_retained_eigenvalue", "nystrom_max_eigenvalue")
                if (!is.list(rkhs_fit) || !all(required_gamma_fit_fields %in% names(rkhs_fit)) || !is.function(rkhs_fit$predict)) {
                  stop(paste0("gamma_fit_provider_factory returned an incompatible gamma ", "fit object."))
                }
                if (!is.null(rkhs_fit$cv_table) && nrow(rkhs_fit$cv_table) > 0L) {
                  cv_table <- as.data.frame(rkhs_fit$cv_table, stringsAsFactors = FALSE)
                  cv_table$fold <- fold
                  cv_table$treatment <- t_value
                  cv_tables[[cv_counter]] <- cv_table
                  cv_counter <- cv_counter + 1L
                }
            }
            else if (is.null(rkhs_lambda_by_t)) {
                rkhs_fit <- fit_gamma_rkhs(train = train, nuisance_train = nuisance_train_estimated,
                  t_value = t_value, x_cols = x_cols, lambda_grid = rkhs_lambda_grid, bandwidth_multipliers = rkhs_bandwidth_multipliers,
                  cv_folds = rkhs_cv_folds, cv_repeats = rkhs_cv_repeats, selection_rule = rkhs_selection_rule,
                  gamma_floor = gamma_floor, gamma_cap = gamma_cap, bandwidth_sample = bandwidth_sample,
                  backend = gamma_rkhs_backend, nystrom_rank = gamma_nystrom_rank, nystrom_eigen_tolerance = gamma_nystrom_eigen_tolerance,
                  nystrom_jitter = gamma_nystrom_jitter, kernel_type = gamma_rkhs_kernel, seed = seed +
                    0x1.86ap+16 * fold + 0x1.9p+6 * t_value)
                cv_table <- rkhs_fit$cv_table
                cv_table$fold <- fold
                cv_table$treatment <- t_value
                cv_tables[[cv_counter]] <- cv_table
                cv_counter <- cv_counter + 0x1p+0
            }
            else {
                if (!lambda_key %in% names(rkhs_lambda_by_t)) {
                  stop("rkhs_lambda_by_t must be named with entries '0' and '1'.")
                }
                rkhs_fit <- fit_gamma_rkhs_fixed(train = train, nuisance_train = nuisance_train_estimated,
                  t_value = t_value, x_cols = x_cols, lambda = rkhs_lambda_by_t[lambda_key], gamma_floor = gamma_floor,
                  gamma_cap = gamma_cap, bandwidth_sample = bandwidth_sample, backend = gamma_rkhs_backend,
                  nystrom_rank = gamma_nystrom_rank, nystrom_eigen_tolerance = gamma_nystrom_eigen_tolerance,
                  nystrom_jitter = gamma_nystrom_jitter, kernel_type = gamma_rkhs_kernel, seed = seed +
                    0x1.86ap+16 * fold + 0x1.9p+6 * t_value)
            }
            gamma_fit_pipeline <- rkhs_fit
            gamma_transform <- NULL
            if (!is.null(gamma_fit_transform_factory)) {
                gamma_transform <- gamma_fit_transform_factory(fold = fold, treatment = t_value, gamma_fit = rkhs_fit,
                  x_cols = x_cols)
                if (!is.list(gamma_transform) || !is.function(gamma_transform$transform)) {
                  stop(paste0("gamma_fit_transform_factory must return a list containing ", "a transform function."))
                }
                original_predict <- rkhs_fit$predict
                transform_function <- gamma_transform$transform
                gamma_fit_pipeline$predict <- local({
                  local_predict <- original_predict
                  local_transform <- transform_function
                  function(newdata) {
                    prediction <- local_predict(newdata)
                    transformed <- as.numeric(local_transform(newdata = newdata, normalized_predictions = prediction$normalized))
                    if (length(transformed) != nrow(newdata) || any(!is.finite(transformed)) || any(transformed <=
                      0x0p+0)) {
                      stop("The fitted-gamma transformation returned invalid predictions.")
                    }
                    prediction$normalized <- transformed
                    prediction
                  }
                })
            }
            gamma_train_rkhs <- gamma_fit_pipeline$predict(train)$normalized
            gamma_validation_rkhs <- gamma_fit_pipeline$predict(validation)$normalized
            if (t_value == 0x1p+0) {
                gamma1_rkhs[validation_index] <- gamma_validation_rkhs
                gamma_oracle_train <- train$gamma1_oracle
                gamma_oracle_validation <- validation$gamma1_oracle
            }
            else {
                gamma0_rkhs[validation_index] <- gamma_validation_rkhs
                gamma_oracle_train <- train$gamma0_oracle
                gamma_oracle_validation <- validation$gamma0_oracle
            }
            baseline_validation <- (validation$T_obs == t_value & validation$Z == 0x0p+0)
            accuracy <- gamma_accuracy_metrics_simple(estimated = gamma_validation_rkhs[baseline_validation],
                oracle = gamma_oracle_validation[baseline_validation])
            accuracy$fold <- fold
            accuracy$treatment <- t_value
            gamma_accuracy_rows[[gamma_accuracy_counter]] <- accuracy
            gamma_accuracy_counter <- gamma_accuracy_counter + 0x1p+0
            selected_lambdas[[selected_lambda_counter]] <- data.frame(fold = fold, treatment = t_value,
                gamma_estimator = if (!is.null(rkhs_fit$gamma_estimator)) {
                  rkhs_fit$gamma_estimator
                }
                else {
                  "direct_rkhs_riesz"
                }, lambda = rkhs_fit$selected_lambda, bandwidth_multiplier = rkhs_fit$selected_bandwidth_multiplier,
                base_bandwidth = rkhs_fit$base_bandwidth, selected_bandwidth = rkhs_fit$selected_bandwidth,
                normalization_constant = rkhs_fit$normalization_constant, training_negative_rate = rkhs_fit$training_negative_rate,
                rkhs_backend = rkhs_fit$rkhs_backend, rkhs_kernel = if (!is.null(rkhs_fit$rkhs_kernel)) {
                  rkhs_fit$rkhs_kernel
                }
                else {
                  gamma_rkhs_kernel
                }, nystrom_rank_requested = rkhs_fit$nystrom_rank_requested, nystrom_landmark_count = rkhs_fit$nystrom_landmark_count,
                nystrom_effective_rank = rkhs_fit$nystrom_effective_rank, nystrom_min_retained_eigenvalue = rkhs_fit$nystrom_min_retained_eigenvalue,
                nystrom_max_eigenvalue = rkhs_fit$nystrom_max_eigenvalue)
            selected_lambda_counter <- selected_lambda_counter + 0x1p+0
            if (!is.null(rkhs_fit$fit_diagnostics)) {
                local_gamma_fit_diagnostics <- as.data.frame(rkhs_fit$fit_diagnostics, stringsAsFactors = FALSE)
                local_gamma_fit_diagnostics$fold <- fold
                local_gamma_fit_diagnostics$treatment <- t_value
                gamma_fit_diagnostics[[gamma_fit_diagnostic_counter]] <- local_gamma_fit_diagnostics
                gamma_fit_diagnostic_counter <- gamma_fit_diagnostic_counter + 1L
            }
            mode_inputs <- list(oracle_all = list(nuisance_train = nuisance_train_oracle, nuisance_validation = nuisance_validation_oracle,
                gamma_train = gamma_oracle_train, gamma_validation = gamma_oracle_validation), oracle_gamma = list(nuisance_train = nuisance_train_estimated,
                nuisance_validation = nuisance_validation_estimated, gamma_train = gamma_oracle_train,
                gamma_validation = gamma_oracle_validation), rkhs = list(nuisance_train = nuisance_train_estimated,
                nuisance_validation = nuisance_validation_estimated, gamma_train = gamma_train_rkhs,
                gamma_validation = gamma_validation_rkhs))
            Lambda0 <- if (t_value == 0x1p+0) {
                Lambda0_T1
            }
            else {
                Lambda0_T0
            }
            Lambda1 <- if (t_value == 0x1p+0) {
                Lambda1_T1
            }
            else {
                Lambda1_T0
            }
            for (mode in modes) {
                bound_component_provider <- NULL
                if (!is.null(bound_component_provider_factory)) {
                  bound_component_provider <- bound_component_provider_factory(mode = mode, fold = fold,
                    treatment = t_value, gamma_fit = gamma_fit_pipeline, nuisance_fit = nuisance_fit,
                    x_cols = x_cols)
                }
                eif_conditional_prediction_provider <- NULL
                if (mode %in% eif_modes && !is.null(eif_conditional_provider_factory)) {
                  eif_conditional_prediction_provider <- eif_conditional_provider_factory(mode = mode,
                    fold = fold, treatment = t_value, gamma_fit = gamma_fit_pipeline, x_cols = x_cols)
                  if (!is.null(eif_conditional_prediction_provider) && !is.function(eif_conditional_prediction_provider)) {
                    stop(paste0("eif_conditional_provider_factory must return NULL or a ", "prediction-provider function."))
                  }
                }
                eif_instrument_probability_validation <- NULL
                if (mode %in% eif_modes) {
                  eif_instrument_probability_validation <- mode_inputs[[mode]]$nuisance_validation$e
                  if (!is.null(eif_instrument_probability_transform)) {
                    eif_instrument_probability_validation <- as.numeric(eif_instrument_probability_transform(probability = eif_instrument_probability_validation,
                      data = validation, mode = mode, fold = fold, treatment = t_value))
                  }
                }
                fit <- fit_bound_arm_fold(train = train, validation = validation, row_id_validation = validation$.row_id,
                  nuisance_train = mode_inputs[[mode]]$nuisance_train, nuisance_validation = mode_inputs[[mode]]$nuisance_validation,
                  gamma_train = mode_inputs[[mode]]$gamma_train, gamma_validation = mode_inputs[[mode]]$gamma_validation,
                  t_value = t_value, Lambda0 = Lambda0, Lambda1 = Lambda1, mode = mode, x_cols = x_cols,
                  conditional_mean_spline_df = conditional_mean_spline_df, threshold_spline_df = threshold_spline_df,
                  threshold_learner = threshold_learner, threshold_control = threshold_control, threshold_seed = seed +
                    500000L * fold + 1000L * t_value + match(mode, modes), envelope_smoothing_epsilon = envelope_smoothing_epsilon,
                  compatibility_tolerance = compatibility_tolerance, compatibility_action = if (mode ==
                    "oracle_all") {
                    "stop"
                  }
                  else {
                    estimated_compatibility_action
                  }, gamma_projection_margin = gamma_projection_margin, max_gamma_projection_rate = max_gamma_projection_rate,
                  tau_epsilon = tau_epsilon, threshold_maxit = threshold_maxit, compute_eif = mode %in%
                    eif_modes && Lambda0 > 0x1p+0 && Lambda1 > 0x1p+0, allow_projected_eif = mode !=
                    "oracle_all", eif_conditional_spline_df = eif_conditional_spline_df, eif_conditional_learner = eif_conditional_learner,
                  eif_conditional_control = eif_conditional_control, eif_conditional_seed = seed + 700000L *
                    fold + 10000L * t_value + 100L * match(mode, modes), eif_conditional_prediction_provider = eif_conditional_prediction_provider,
                  eif_instrument_probability_validation = eif_instrument_probability_validation, bound_component_provider = bound_component_provider,
                  fixed_threshold_target_evaluator = fixed_threshold_target_evaluator)
                contribution_storage[[mode]][[lambda_key]]$lower[validation_index] <- fit$lower_contribution
                contribution_storage[[mode]][[lambda_key]]$upper[validation_index] <- fit$upper_contribution
                contribution_storage[[mode]][[lambda_key]]$lambda_one[validation_index] <- fit$lambda_one_contribution
                contribution_storage[[mode]][[lambda_key]]$lower_pseudo[validation_index] <- fit$lower_pseudo_contribution
                contribution_storage[[mode]][[lambda_key]]$upper_pseudo[validation_index] <- fit$upper_pseudo_contribution
                diagnostic <- fit$diagnostics
                diagnostic$mode <- mode
                diagnostic$fold <- fold
                fold_diagnostics[[diagnostic_counter]] <- diagnostic
                diagnostic_counter <- diagnostic_counter + 0x1p+0
                prediction <- fit$predictions
                prediction$mode <- mode
                prediction$fold <- fold
                prediction_tables[[prediction_counter]] <- prediction
                prediction_counter <- prediction_counter + 0x1p+0
                if (nrow(fit$threshold_cv) > 0L) {
                  local_threshold_cv <- fit$threshold_cv
                  local_threshold_cv$mode <- mode
                  local_threshold_cv$fold <- fold
                  local_threshold_cv$treatment <- t_value
                  threshold_cv_tables[[threshold_cv_counter]] <- local_threshold_cv
                  threshold_cv_counter <- threshold_cv_counter + 1L
                }
                if (nrow(fit$fixed_threshold_target_evaluation) > 0L) {
                  local_fixed_target <- fit$fixed_threshold_target_evaluation
                  local_fixed_target$mode <- mode
                  local_fixed_target$fold <- fold
                  local_fixed_target$treatment <- t_value
                  fixed_threshold_target_evaluations[[fixed_threshold_target_evaluation_counter]] <- local_fixed_target
                  fixed_threshold_target_evaluation_counter <- fixed_threshold_target_evaluation_counter +
                    1L
                }
                if (!is.null(fit$eif)) {
                  eif_storage[[mode]][[lambda_key]]$lower[validation_index] <- fit$eif$signals$lower_signal
                  eif_storage[[mode]][[lambda_key]]$upper[validation_index] <- fit$eif$signals$upper_signal
                  eif_prediction <- fit$eif$signals
                  eif_prediction$mode <- mode
                  eif_prediction$fold <- fold
                  eif_predictions[[eif_prediction_counter]] <- eif_prediction
                  eif_prediction_counter <- eif_prediction_counter + 1L
                  eif_diagnostic <- fit$eif$diagnostics
                  eif_diagnostic$mode <- mode
                  eif_diagnostic$fold <- fold
                  eif_diagnostics[[eif_diagnostic_counter]] <- eif_diagnostic
                  eif_diagnostic_counter <- eif_diagnostic_counter + 1L
                  if (nrow(fit$eif$conditional_cv) > 0L) {
                    local_conditional_cv <- fit$eif$conditional_cv
                    local_conditional_cv$mode <- mode
                    local_conditional_cv$fold <- fold
                    local_conditional_cv$treatment <- t_value
                    eif_conditional_cv_tables[[eif_conditional_cv_counter]] <- local_conditional_cv
                    eif_conditional_cv_counter <- eif_conditional_cv_counter + 1L
                  }
                }
            }
            if (!is.null(gamma_transform) && is.function(gamma_transform$diagnostics)) {
                local_gamma_transform_diagnostics <- gamma_transform$diagnostics()
                if (!is.data.frame(local_gamma_transform_diagnostics)) {
                  stop("gamma transformation diagnostics must be a data frame.")
                }
                local_gamma_transform_diagnostics$fold <- fold
                local_gamma_transform_diagnostics$treatment <- t_value
                gamma_transform_diagnostics[[gamma_transform_diagnostic_counter]] <- local_gamma_transform_diagnostics
                gamma_transform_diagnostic_counter <- gamma_transform_diagnostic_counter + 1L
            }
        }
    }
    analysis_data$e_hat <- e_hat
    analysis_data$pi0_hat <- pi0_hat
    analysis_data$pi1_hat <- pi1_hat
    analysis_data$gamma0_rkhs <- gamma0_rkhs
    analysis_data$gamma1_rkhs <- gamma1_rkhs
    summary_rows <- list()
    summary_counter <- 0x1p+0
    pseudo_summary_rows <- list()
    pseudo_summary_counter <- 0x1p+0
    unconfounded_rows <- list()
    unconfounded_counter <- 0x1p+0
    for (mode in modes) {
        mu1_lower <- mean(contribution_storage[[mode]][["1"]]$lower)
        mu1_upper <- mean(contribution_storage[[mode]][["1"]]$upper)
        mu0_lower <- mean(contribution_storage[[mode]][["0"]]$lower)
        mu0_upper <- mean(contribution_storage[[mode]][["0"]]$upper)
        summary_rows[[summary_counter]] <- data.frame(mode = mode, mu1_lower = mu1_lower, mu1_upper = mu1_upper,
            mu0_lower = mu0_lower, mu0_upper = mu0_upper, ate_lower = mu1_lower - mu0_upper, ate_upper = mu1_upper -
                mu0_lower)
        summary_counter <- summary_counter + 0x1p+0
        pseudo_mu1_lower <- mean(contribution_storage[[mode]][["1"]]$lower_pseudo)
        pseudo_mu1_upper <- mean(contribution_storage[[mode]][["1"]]$upper_pseudo)
        pseudo_mu0_lower <- mean(contribution_storage[[mode]][["0"]]$lower_pseudo)
        pseudo_mu0_upper <- mean(contribution_storage[[mode]][["0"]]$upper_pseudo)
        pseudo_summary_rows[[pseudo_summary_counter]] <- data.frame(mode = mode, mu1_lower = pseudo_mu1_lower,
            mu1_upper = pseudo_mu1_upper, mu0_lower = pseudo_mu0_lower, mu0_upper = pseudo_mu0_upper,
            ate_lower = pseudo_mu1_lower - pseudo_mu0_upper, ate_upper = pseudo_mu1_upper - pseudo_mu0_lower)
        pseudo_summary_counter <- pseudo_summary_counter + 1L
        unconfounded_rows[[unconfounded_counter]] <- data.frame(mode = mode, mu1_no_unmeasured_confounding_benchmark = mean(contribution_storage[[mode]][["1"]]$lambda_one),
            mu0_no_unmeasured_confounding_benchmark = mean(contribution_storage[[mode]][["0"]]$lambda_one))
        unconfounded_counter <- unconfounded_counter + 0x1p+0
    }
    prediction_table <- do.call(rbind, prediction_tables)
    fold_diagnostic_table <- do.call(rbind, fold_diagnostics)
    gamma_accuracy_table <- do.call(rbind, gamma_accuracy_rows)
    selected_lambda_table <- do.call(rbind, selected_lambdas)
    summary_table <- do.call(rbind, summary_rows)
    pseudo_summary_table <- do.call(rbind, pseudo_summary_rows)
    endpoint_names <- c("mu1_lower", "mu1_upper", "mu0_lower", "mu0_upper", "ate_lower", "ate_upper")
    eif_complete_modes <- eif_modes[vapply(eif_modes, function(mode) {
        all(vapply(eif_storage[[mode]], function(arm) {
            all(is.finite(arm$lower)) && all(is.finite(arm$upper))
        }, logical(0x1p+0)))
    }, logical(0x1p+0))]
    one_step_signal_tables <- lapply(eif_complete_modes, function(mode) {
        local <- data.frame(row_id = seq_len(n), mode = mode, mu1_lower = eif_storage[[mode]][["1"]]$lower,
            mu1_upper = eif_storage[[mode]][["1"]]$upper, mu0_lower = eif_storage[[mode]][["0"]]$lower,
            mu0_upper = eif_storage[[mode]][["0"]]$upper, stringsAsFactors = FALSE)
        local$ate_lower <- local$mu1_lower - local$mu0_upper
        local$ate_upper <- local$mu1_upper - local$mu0_lower
        local
    })
    one_step_signal_table <- if (length(one_step_signal_tables) > 0L) {
        do.call(rbind, one_step_signal_tables)
    }
    else {
        data.frame()
    }
    eif_prediction_table <- if (length(eif_predictions) > 0L) {
        do.call(rbind, eif_predictions)
    }
    else {
        data.frame()
    }
    eif_diagnostic_table <- if (length(eif_diagnostics) > 0L) {
        do.call(rbind, eif_diagnostics)
    }
    else {
        data.frame()
    }
    one_step_summary <- if (length(eif_complete_modes) > 0L) {
        do.call(rbind, lapply(eif_complete_modes, function(mode) {
            signals <- one_step_signal_table[one_step_signal_table$mode == mode, , drop = FALSE]
            plugin <- summary_table[summary_table$mode == mode, , drop = FALSE]
            local_diagnostics <- fold_diagnostic_table[fold_diagnostic_table$mode == mode, , drop = FALSE]
            local_eif_diagnostics <- eif_diagnostic_table[eif_diagnostic_table$mode == mode, , drop = FALSE]
            maximum_gamma_projection_rate <- max(local_eif_diagnostics$gamma_projection_rate_train, local_eif_diagnostics$gamma_projection_rate_validation)
            maximum_conditional_projection_rate <- max(local_diagnostics$conditional_projection_rate_train,
                local_diagnostics$conditional_projection_rate_validation)
            projection_active <- any(local_eif_diagnostics$gamma_projection_active)
            do.call(rbind, lapply(endpoint_names, function(endpoint) {
                signal <- signals[[endpoint]]
                estimate <- mean(signal)
                standard_error <- stats::sd(signal)/sqrt(length(signal))
                data.frame(mode = mode, gamma_nuisance = if (mode == "rkhs") {
                  paste0("estimated_", unique(selected_lambda_table$gamma_estimator)[1L])
                }
                else {
                  "oracle"
                }, threshold_learner = threshold_learner, endpoint = endpoint, plugin_estimate = plugin[[endpoint]][1L],
                  one_step_estimate = estimate, one_step_correction = estimate - plugin[[endpoint]][1L],
                  standard_error = standard_error, confidence_level = 0x1.e666666666666p-1, confidence_lower = estimate -
                    stats::qnorm(0x1.f333333333333p-1) * standard_error, confidence_upper = estimate +
                    stats::qnorm(0x1.f333333333333p-1) * standard_error, maximum_gamma_projection_rate = maximum_gamma_projection_rate,
                  maximum_conditional_projection_rate = maximum_conditional_projection_rate, projected_nuisance_inference_flag = projection_active,
                  inference_status = if (projection_active) {
                    "if_with_projected_gamma_diagnostic"
                  }
                  else {
                    "if_under_nuisance_rate_conditions"
                  }, stringsAsFactors = FALSE)
            }))
        }))
    }
    else {
        data.frame()
    }
    oracle_eif_signal_table <- if ("oracle_all" %in% eif_complete_modes) {
        local <- one_step_signal_table[one_step_signal_table$mode == "oracle_all", , drop = FALSE]
        local$mode <- NULL
        local
    }
    else {
        data.frame()
    }
    oracle_one_step_summary <- if (nrow(one_step_summary) > 0L) {
        one_step_summary[one_step_summary$mode == "oracle_all", , drop = FALSE]
    }
    else {
        data.frame()
    }
    list(summary = summary_table, pseudo_summary = pseudo_summary_table, truth = data.frame(sample_mu1 = if ("Y1" %in%
        names(analysis_data)) {
        mean(analysis_data$Y1)
    } else {
        NA_real_
    }, sample_mu0 = if ("Y0" %in% names(analysis_data)) {
        mean(analysis_data$Y0)
    } else {
        NA_real_
    }, sample_ate = if (all(c("Y1", "Y0") %in% names(analysis_data))) {
        mean(analysis_data$Y1 - analysis_data$Y0)
    } else {
        NA_real_
    }, population_mu1 = treatment_effect, population_mu0 = 0x0p+0, population_ate = treatment_effect),
        unconfounded_checks = do.call(rbind, unconfounded_rows), gamma_accuracy_by_fold = gamma_accuracy_table,
        gamma_accuracy = stats::aggregate(cbind(bias, rmse, mae, correlation, mean_estimated, mean_oracle,
            rmse_constant_one, explained_fraction) ~ treatment, data = gamma_accuracy_table, FUN = mean,
            na.action = stats::na.pass), nuisance = compute_nuisance_diagnostics_simple(analysis_data),
        nuisance_learner = nuisance_learner, nuisance_fit_diagnostics = do.call(rbind, nuisance_fit_diagnostics),
        ipw_identities = compute_ipw_identities_simple(analysis_data), selected_lambdas = selected_lambda_table,
        gamma_fit_provider_active = !is.null(gamma_fit_provider_factory), gamma_fit_diagnostics = if (length(gamma_fit_diagnostics) >
            0L) {
            do.call(rbind, gamma_fit_diagnostics)
        } else {
            data.frame()
        }, rkhs_cv = if (length(cv_tables) > 0x0p+0) {
            do.call(rbind, cv_tables)
        } else {
            data.frame()
        }, threshold_cv = if (length(threshold_cv_tables) > 0L) {
            do.call(rbind, threshold_cv_tables)
        } else {
            data.frame()
        }, fixed_threshold_target_evaluation = if (length(fixed_threshold_target_evaluations) > 0L) {
            do.call(rbind, fixed_threshold_target_evaluations)
        } else {
            data.frame()
        }, fold_diagnostics = fold_diagnostic_table, one_step = one_step_summary, fully_estimated_one_step = if (nrow(one_step_summary) >
            0L) {
            one_step_summary[one_step_summary$mode == "rkhs", , drop = FALSE]
        } else {
            data.frame()
        }, one_step_signals = one_step_signal_table, eif_predictions = eif_prediction_table, eif_diagnostics = eif_diagnostic_table,
        eif_conditional_cv = if (length(eif_conditional_cv_tables) > 0L) {
            do.call(rbind, eif_conditional_cv_tables)
        } else {
            data.frame()
        }, eif_conditional_provider_active = !is.null(eif_conditional_provider_factory), gamma_fit_transform_active = !is.null(gamma_fit_transform_factory),
        gamma_transform_diagnostics = if (length(gamma_transform_diagnostics) > 0L) {
            do.call(rbind, gamma_transform_diagnostics)
        } else {
            data.frame()
        }, oracle_one_step = oracle_one_step_summary, oracle_eif_signals = oracle_eif_signal_table, oracle_eif_predictions = if (nrow(eif_prediction_table) >
            0L) {
            eif_prediction_table[eif_prediction_table$mode == "oracle_all", , drop = FALSE]
        } else {
            data.frame()
        }, oracle_eif_diagnostics = if (nrow(eif_diagnostic_table) > 0L) {
            eif_diagnostic_table[eif_diagnostic_table$mode == "oracle_all", , drop = FALSE]
        } else {
            data.frame()
        }, bound_predictions = prediction_table, normalization_by_x_bin = build_normalization_by_x_bin(predictions = prediction_table,
            x_bins = x_bins), data = analysis_data, fold_id = fold_id, envelope_smoothing = data.frame(method = if (envelope_smoothing_epsilon >
            0x0p+0) {
            "conservative_log_mean_exp"
        } else {
            "none"
        }, epsilon = envelope_smoothing_epsilon, bandwidth = if (envelope_smoothing_epsilon > 0x0p+0) {
            envelope_smoothing_epsilon/log(0x1p+1)
        } else {
            0x0p+0
        }, stringsAsFactors = FALSE))
}

evaluate_threshold_bound_contributions <- function (data, nuisance, envelope, t_value, tau, q_lower, q_upper)
{
    n <- nrow(data)
    input_lengths <- c(n, length(nuisance$e), length(tau), length(q_lower), length(q_upper), length(envelope$L),
        length(envelope$U), length(envelope$D))
    if (length(unique(input_lengths)) != 1L) {
        stop("Threshold-bound contribution inputs must have a common length.")
    }
    baseline <- data$T_obs == t_value & data$Z == 0x0p+0
    odds <- compute_target_odds(nuisance, t_value)
    inverse_instrument <- 0x1p+0/nuisance$e
    lower_below <- as.numeric(data$Y <= q_lower)
    upper_below <- as.numeric(data$Y <= q_upper)
    lambda_lower <- envelope$L + envelope$D * lower_below
    lambda_upper <- envelope$U - envelope$D * upper_below
    lower_level <- 0x1p+0 - tau
    upper_level <- tau
    interpolated_lambda <- tau * envelope$L + (0x1p+0 - tau) * envelope$U
    direct_lower <- direct_upper <- pseudo_lower <- pseudo_upper <- lambda_one <- numeric(n)
    lower_score_residual <- upper_score_residual <- numeric(n)
    direct_lower[baseline] <- inverse_instrument[baseline] * data$Y[baseline] * (0x1p+0 + odds[baseline] *
        lambda_lower[baseline])
    direct_upper[baseline] <- inverse_instrument[baseline] * data$Y[baseline] * (0x1p+0 + odds[baseline] *
        lambda_upper[baseline])
    lambda_one[baseline] <- inverse_instrument[baseline] * data$Y[baseline] * (0x1p+0 + odds[baseline])
    pseudo_lower[baseline] <- inverse_instrument[baseline] * (data$Y[baseline] * (0x1p+0 + odds[baseline] *
        interpolated_lambda[baseline]) - odds[baseline] * envelope$D[baseline] * check_loss(data$Y[baseline],
        q_lower[baseline], lower_level[baseline]))
    pseudo_upper[baseline] <- inverse_instrument[baseline] * (data$Y[baseline] * (0x1p+0 + odds[baseline] *
        interpolated_lambda[baseline]) + odds[baseline] * envelope$D[baseline] * check_loss(data$Y[baseline],
        q_upper[baseline], upper_level[baseline]))
    lower_score_residual[baseline] <- inverse_instrument[baseline] * odds[baseline] * envelope$D[baseline] *
        q_lower[baseline] * (lower_below[baseline] - lower_level[baseline])
    upper_score_residual[baseline] <- -inverse_instrument[baseline] * odds[baseline] * envelope$D[baseline] *
        q_upper[baseline] * (upper_below[baseline] - upper_level[baseline])
    list(baseline = baseline, lambda_lower = lambda_lower, lambda_upper = lambda_upper, direct_lower = direct_lower,
        direct_upper = direct_upper, pseudo_lower = pseudo_lower, pseudo_upper = pseudo_upper, lambda_one = lambda_one,
        lower_score_residual = lower_score_residual, upper_score_residual = upper_score_residual)
}

first_primes <- function (count)
{
    count <- as.integer(count)
    if (length(count) != 1L || is.na(count) || count < 1L) {
        stop("count must be one positive integer.")
    }
    primes <- integer()
    candidate <- 2L
    while (length(primes) < count) {
        divisor <- primes[primes <= floor(sqrt(candidate))]
        is_prime <- !length(divisor) || all(candidate%%divisor != 0L)
        if (is_prime)
            primes <- c(primes, candidate)
        candidate <- candidate + 1L
    }
    primes
}

fit_bound_arm_fold <- function (train, validation, row_id_validation, nuisance_train, nuisance_validation, gamma_train, gamma_validation,
    t_value, Lambda0, Lambda1, mode, x_cols = "X", conditional_mean_spline_df = 0x1.4p+2, threshold_spline_df = 0x1.4p+2,
    threshold_learner = "external", threshold_control = NULL, threshold_seed = 1L, envelope_smoothing_epsilon = 0x0p+0,
    compatibility_tolerance = 0x1.0c6f7a0b5ed8ep-20, compatibility_action = c("stop", "project"), gamma_projection_margin = 0x1.a36e2eb1c432dp-14,
    max_gamma_projection_rate = 0x1.47ae147ae147bp-6, tau_epsilon = 0x1.a36e2eb1c432dp-14, threshold_maxit = 0x1.f4p+11,
    compute_eif = FALSE, allow_projected_eif = FALSE, eif_conditional_spline_df = conditional_mean_spline_df,
    eif_conditional_learner = c("spline_fixed", "spline_cv", "gam", "ranger", "library_cv"), eif_conditional_control = NULL,
    eif_conditional_seed = threshold_seed + 20000L, eif_conditional_prediction_provider = NULL, eif_instrument_probability_validation = NULL,
    bound_component_provider = NULL, fixed_threshold_target_evaluator = NULL)
{
    stopifnot(identical(threshold_learner, "external"))
    threshold_learner <- match.arg(threshold_learner)
    eif_conditional_learner <- match.arg(eif_conditional_learner)
    compatibility_action <- match.arg(compatibility_action)
    if (!is.null(bound_component_provider) && !is.list(bound_component_provider)) {
        stop("bound_component_provider must be NULL or a list.")
    }
    if (!is.null(bound_component_provider$envelope_means) && !is.function(bound_component_provider$envelope_means)) {
        stop("bound_component_provider$envelope_means must be a function.")
    }
    if (!is.null(bound_component_provider$thresholds) && !is.function(bound_component_provider$thresholds)) {
        stop("bound_component_provider$thresholds must be a function.")
    }
    if (threshold_learner == "external" && (is.null(bound_component_provider) || !is.function(bound_component_provider$thresholds))) {
        stop(paste0("threshold_learner='external' requires ", "bound_component_provider$thresholds."))
    }
    if (!is.null(fixed_threshold_target_evaluator) && !is.function(fixed_threshold_target_evaluator)) {
        stop("fixed_threshold_target_evaluator must be NULL or a function.")
    }
    envelope_train <- compute_projected_envelopes(nuisance = nuisance_train, gamma = gamma_train, t_value = t_value,
        Lambda0 = Lambda0, Lambda1 = Lambda1, envelope_smoothing_epsilon = envelope_smoothing_epsilon,
        gamma_compatibility_action = compatibility_action, gamma_projection_margin = gamma_projection_margin)
    envelope_validation <- compute_projected_envelopes(nuisance = nuisance_validation, gamma = gamma_validation,
        t_value = t_value, Lambda0 = Lambda0, Lambda1 = Lambda1, envelope_smoothing_epsilon = envelope_smoothing_epsilon,
        gamma_compatibility_action = compatibility_action, gamma_projection_margin = gamma_projection_margin)
    baseline_train <- (train$T_obs == t_value & train$Z == 0x0p+0)
    baseline_validation <- (validation$T_obs == t_value & validation$Z == 0x0p+0)
    pointwise_crossing_train <- (envelope_train$L > envelope_train$U + compatibility_tolerance)
    pointwise_crossing_validation <- (envelope_validation$L > envelope_validation$U + compatibility_tolerance)
    gamma_projection_rate_train <- mean(envelope_train$gamma_projected[baseline_train])
    gamma_projection_rate_validation <- mean(envelope_validation$gamma_projected[baseline_validation])
    gamma_incompatibility_rate_train <- mean(envelope_train$gamma_incompatible_original[baseline_train])
    gamma_incompatibility_rate_validation <- mean(envelope_validation$gamma_incompatible_original[baseline_validation])
    gamma_projection_threshold_exceeded <- (compatibility_action == "project" && (gamma_projection_rate_train >
        max_gamma_projection_rate || gamma_projection_rate_validation > max_gamma_projection_rate))
    if (any(pointwise_crossing_train[baseline_train]) || any(pointwise_crossing_validation[baseline_validation])) {
        stop(sprintf(paste0("Pointwise projected envelopes cross for mode %s, ", "treatment %d after the requested compatibility action. ",
            "Training rate %.4f; validation rate %.4f."), mode, t_value, mean(pointwise_crossing_train[baseline_train]),
            mean(pointwise_crossing_validation[baseline_validation])))
    }
    if (Lambda0 == 0x1p+0 && Lambda1 == 0x1p+0) {
        target_odds_validation <- compute_target_odds(nuisance = nuisance_validation, t_value = t_value)
        contribution <- rep(0x0p+0, nrow(validation))
        contribution[baseline_validation] <- validation$Y[baseline_validation]/nuisance_validation$e[baseline_validation] *
            (0x1p+0 + target_odds_validation[baseline_validation])
        prediction_table <- data.frame(row_id = row_id_validation[baseline_validation], treatment = t_value,
            Y = validation$Y[baseline_validation], X = validation[baseline_validation, x_cols[0x1p+0],
                drop = TRUE], L = 0x1p+0, U = 0x1p+0, D = 0x0p+0, hard_L = 0x1p+0, hard_U = 0x1p+0, lower_smoothing_gap = 0x0p+0,
            upper_smoothing_gap = 0x0p+0, Lbar = 0x1p+0, Ubar = 0x1p+0, tau_raw = NA_real_, tau = NA_real_,
            q_lower = NA_real_, q_upper = NA_real_, lambda_lower = 0x1p+0, lambda_upper = 0x1p+0, direct_lower_contribution = contribution[baseline_validation],
            direct_upper_contribution = contribution[baseline_validation], pseudo_lower_contribution = contribution[baseline_validation],
            pseudo_upper_contribution = contribution[baseline_validation], lower_q_score_residual = 0x0p+0,
            upper_q_score_residual = 0x0p+0, threshold_learner = threshold_learner, threshold_rkhs_backend = NA_character_,
            gamma_original = envelope_validation$gamma_original[baseline_validation], gamma_used = envelope_validation$gamma_used[baseline_validation],
            gamma_compatible_lower = envelope_validation$gamma_compatible_lower[baseline_validation],
            gamma_compatible_upper = envelope_validation$gamma_compatible_upper[baseline_validation],
            gamma_projected = envelope_validation$gamma_projected[baseline_validation])
        diagnostics <- data.frame(treatment = t_value, n_train_baseline = sum(baseline_train), n_validation_baseline = sum(baseline_validation),
            pointwise_crossing_rate_train = 0x0p+0, pointwise_crossing_rate_validation = 0x0p+0, gamma_incompatibility_rate_train = gamma_incompatibility_rate_train,
            gamma_incompatibility_rate_validation = gamma_incompatibility_rate_validation, gamma_projection_rate_train = 0x0p+0,
            gamma_projection_rate_validation = 0x0p+0, gamma_projection_threshold_exceeded = FALSE, envelope_smoothing_epsilon = 0x0p+0,
            mean_lower_smoothing_gap_validation = 0x0p+0, mean_upper_smoothing_gap_validation = 0x0p+0,
            maximum_lower_smoothing_gap_validation = 0x0p+0, maximum_upper_smoothing_gap_validation = 0x0p+0,
            conditional_projection_rate_train = 0x0p+0, conditional_projection_rate_validation = 0x0p+0,
            tau_invalid_rate_train = NA_real_, tau_invalid_rate_validation = NA_real_, lower_threshold_convergence = NA_integer_,
            upper_threshold_convergence = NA_integer_, lower_threshold_objective = NA_real_, upper_threshold_objective = NA_real_,
            mean_L_validation = 0x1p+0, mean_U_validation = 0x1p+0, mean_tau_validation = NA_real_, mean_lambda_lower_validation = 0x1p+0,
            mean_lambda_upper_validation = 0x1p+0, threshold_learner = threshold_learner, threshold_rkhs_backend = NA_character_,
            lower_selected_spline_df = NA_integer_, upper_selected_spline_df = NA_integer_, lower_selected_lambda = NA_real_,
            upper_selected_lambda = NA_real_, lower_selected_bandwidth_multiplier = NA_real_, upper_selected_bandwidth_multiplier = NA_real_,
            lower_selected_bandwidth = NA_real_, upper_selected_bandwidth = NA_real_, lower_base_bandwidth = NA_real_,
            upper_base_bandwidth = NA_real_, lower_selected_rank = NA_integer_, upper_selected_rank = NA_integer_,
            lower_threshold_iterations = NA_real_, upper_threshold_iterations = NA_real_, lower_threshold_initial_convergence = NA_integer_,
            upper_threshold_initial_convergence = NA_integer_, lower_threshold_retry_used = NA, upper_threshold_retry_used = NA,
            lower_threshold_retry_accepted = NA, upper_threshold_retry_accepted = NA, lower_threshold_retry_convergence = NA_integer_,
            upper_threshold_retry_convergence = NA_integer_, lower_threshold_initial_penalized_objective = NA_real_,
            upper_threshold_initial_penalized_objective = NA_real_, lower_threshold_penalized_objective = NA_real_,
            upper_threshold_penalized_objective = NA_real_, lower_threshold_tie_break_ridge = NA_real_,
            upper_threshold_tie_break_ridge = NA_real_, lower_threshold_tie_break_ridge_numerator = NA_real_,
            upper_threshold_tie_break_ridge_numerator = NA_real_, lower_threshold_tie_break_ridge_training_size = NA_integer_,
            upper_threshold_tie_break_ridge_training_size = NA_integer_, lower_threshold_primal_residual = NA_real_,
            upper_threshold_primal_residual = NA_real_, lower_threshold_dual_residual = NA_real_, upper_threshold_dual_residual = NA_real_,
            lower_threshold_final_rho = NA_real_, upper_threshold_final_rho = NA_real_, lower_threshold_reciprocal_condition = NA_real_,
            upper_threshold_reciprocal_condition = NA_real_, lower_threshold_numerical_jitter = NA_real_,
            upper_threshold_numerical_jitter = NA_real_, lower_kernel_minimum_raw_eigenvalue = NA_real_,
            upper_kernel_minimum_raw_eigenvalue = NA_real_, fold_lower_mean = mean(contribution), fold_upper_mean = mean(contribution),
            fold_lambda_one_mean = mean(contribution), fold_lower_pseudo_mean = mean(contribution), fold_upper_pseudo_mean = mean(contribution),
            lower_direct_minus_pseudo_train = 0x0p+0, upper_direct_minus_pseudo_train = 0x0p+0, lower_direct_minus_pseudo_validation = 0x0p+0,
            upper_direct_minus_pseudo_validation = 0x0p+0, lower_q_score_residual_train = 0x0p+0, upper_q_score_residual_train = 0x0p+0,
            lower_q_score_residual_validation = 0x0p+0, upper_q_score_residual_validation = 0x0p+0, lower_direct_pseudo_identity_error_validation = 0x0p+0,
            upper_direct_pseudo_identity_error_validation = 0x0p+0)
        return(list(lower_contribution = contribution, upper_contribution = contribution, lower_pseudo_contribution = contribution,
            upper_pseudo_contribution = contribution, lambda_one_contribution = contribution, predictions = prediction_table,
            diagnostics = diagnostics, threshold_cv = data.frame(), fixed_threshold_target_evaluation = data.frame(),
            eif = NULL))
    }
    if (sum(baseline_train) <= conditional_mean_spline_df + 0x1p+1) {
        stop(sprintf("Too few T=%d,Z=0 training observations.", t_value))
    }
    regression_data <- train[baseline_train, , drop = FALSE]
    regression_data$.L <- envelope_train$L[baseline_train]
    regression_data$.U <- envelope_train$U[baseline_train]
    Lbar_fit <- fit_spline_mean(data = regression_data, outcome = ".L", x_cols = x_cols, spline_df = conditional_mean_spline_df)
    Ubar_fit <- fit_spline_mean(data = regression_data, outcome = ".U", x_cols = x_cols, spline_df = conditional_mean_spline_df)
    Lbar_train <- Lbar_fit$predict(train)
    Ubar_train <- Ubar_fit$predict(train)
    Lbar_validation <- Lbar_fit$predict(validation)
    Ubar_validation <- Ubar_fit$predict(validation)
    if (!is.null(bound_component_provider) && is.function(bound_component_provider$envelope_means)) {
        provider_envelope_means <- bound_component_provider$envelope_means(train = train, validation = validation,
            nuisance_train = nuisance_train, nuisance_validation = nuisance_validation, envelope_train = envelope_train,
            envelope_validation = envelope_validation, treatment = t_value, Lambda0 = Lambda0, Lambda1 = Lambda1,
            compatibility_action = compatibility_action, gamma_projection_margin = gamma_projection_margin,
            envelope_smoothing_epsilon = envelope_smoothing_epsilon)
        required_envelope_means <- c("Lbar_train", "Ubar_train", "Lbar_validation", "Ubar_validation")
        if (!is.list(provider_envelope_means) || !all(required_envelope_means %in% names(provider_envelope_means))) {
            stop(paste0("bound_component_provider$envelope_means must return ", paste(required_envelope_means,
                collapse = ", "), "."))
        }
        provider_lengths <- vapply(provider_envelope_means[required_envelope_means], length, integer(0x1p+0))
        if (!identical(unname(provider_lengths), c(nrow(train), nrow(train), nrow(validation), nrow(validation))) ||
            any(!is.finite(unlist(provider_envelope_means[required_envelope_means])))) {
            stop("The external envelope-mean provider returned invalid predictions.")
        }
        Lbar_train <- as.numeric(provider_envelope_means$Lbar_train)
        Ubar_train <- as.numeric(provider_envelope_means$Ubar_train)
        Lbar_validation <- as.numeric(provider_envelope_means$Lbar_validation)
        Ubar_validation <- as.numeric(provider_envelope_means$Ubar_validation)
    }
    invalid_train <- (Lbar_train > 0x1p+0 + compatibility_tolerance | Ubar_train < 0x1p+0 - compatibility_tolerance |
        Ubar_train - Lbar_train <= compatibility_tolerance)
    invalid_validation <- (Lbar_validation > 0x1p+0 + compatibility_tolerance | Ubar_validation < 0x1p+0 -
        compatibility_tolerance | Ubar_validation - Lbar_validation <= compatibility_tolerance)
    projection_rate_train <- mean(invalid_train[baseline_train])
    projection_rate_validation <- mean(invalid_validation[baseline_validation])
    Lbar_train <- pmin(Lbar_train, 0x1p+0 - tau_epsilon)
    Ubar_train <- pmax(Ubar_train, 0x1p+0 + tau_epsilon)
    Lbar_validation <- pmin(Lbar_validation, 0x1p+0 - tau_epsilon)
    Ubar_validation <- pmax(Ubar_validation, 0x1p+0 + tau_epsilon)
    tau_train_raw <- (Ubar_train - 0x1p+0)/(Ubar_train - Lbar_train)
    tau_validation_raw <- (Ubar_validation - 0x1p+0)/(Ubar_validation - Lbar_validation)
    tau_train <- tau_train_raw
    tau_validation <- tau_validation_raw
    tau_invalid_train <- (!is.finite(tau_train) | tau_train <= 0x0p+0 | tau_train >= 0x1p+0)
    tau_invalid_validation <- (!is.finite(tau_validation) | tau_validation <= 0x0p+0 | tau_validation >=
        0x1p+0)
    tau_train[!is.finite(tau_train)] <- 0x1p-1
    tau_validation[!is.finite(tau_validation)] <- 0x1p-1
    tau_train <- pmin(0x1p+0 - tau_epsilon, pmax(tau_epsilon, tau_train))
    tau_validation <- pmin(0x1p+0 - tau_epsilon, pmax(tau_epsilon, tau_validation))
    target_odds_train <- compute_target_odds(nuisance = nuisance_train, t_value = t_value)
    threshold_weight <- (0x1p+0/nuisance_train$e[baseline_train]) * target_odds_train[baseline_train] *
        envelope_train$D[baseline_train]
    if (any(threshold_weight < 0x0p+0)) {
        stop(sprintf("Negative threshold weights for treatment %d.", t_value))
    }
    threshold_rkhs_plan <- {
        NULL
    }
    make_external_threshold_fit <- function(predict_function, q_train, level) {
        list(predict = predict_function, convergence = 0L, objective = sum(threshold_weight * check_loss(train$Y[baseline_train],
            q_train[baseline_train], level[baseline_train])), penalized_objective = NA_real_, counts = c(iterations = 0L),
            selected_spline_df = NA_integer_, selected_lambda = NA_real_, selected_bandwidth_multiplier = NA_real_,
            selected_bandwidth = NA_real_, base_bandwidth = NA_real_, selected_rank = NA_integer_, cv_table = data.frame(),
            primal_residual = NA_real_, dual_residual = NA_real_, final_rho = NA_real_, rkhs_backend = "external_parametric",
            minimum_reciprocal_condition = NA_real_, maximum_numerical_jitter = NA_real_, kernel_minimum_raw_eigenvalue = NA_real_,
            tie_break_ridge = NA_real_, tie_break_ridge_numerator = NA_real_, tie_break_ridge_training_size = sum(baseline_train))
    }
    {
        provider_thresholds <- bound_component_provider$thresholds(train = train, validation = validation,
            nuisance_train = nuisance_train, nuisance_validation = nuisance_validation, envelope_train = envelope_train,
            envelope_validation = envelope_validation, treatment = t_value, Lambda0 = Lambda0, Lambda1 = Lambda1,
            tau_train = tau_train, tau_validation = tau_validation, compatibility_action = compatibility_action,
            gamma_projection_margin = gamma_projection_margin, envelope_smoothing_epsilon = envelope_smoothing_epsilon,
            tau_epsilon = tau_epsilon)
        required_thresholds <- c("q_lower_train", "q_upper_train", "q_lower_validation", "q_upper_validation",
            "lower_predict", "upper_predict")
        if (!is.list(provider_thresholds) || !all(required_thresholds %in% names(provider_thresholds)) ||
            !is.function(provider_thresholds$lower_predict) || !is.function(provider_thresholds$upper_predict)) {
            stop(paste0("bound_component_provider$thresholds returned an incompatible ", "threshold object."))
        }
        q_lower_train <- as.numeric(provider_thresholds$q_lower_train)
        q_upper_train <- as.numeric(provider_thresholds$q_upper_train)
        q_lower_validation <- as.numeric(provider_thresholds$q_lower_validation)
        q_upper_validation <- as.numeric(provider_thresholds$q_upper_validation)
        if (length(q_lower_train) != nrow(train) || length(q_upper_train) != nrow(train) || length(q_lower_validation) !=
            nrow(validation) || length(q_upper_validation) != nrow(validation) || any(!is.finite(c(q_lower_train,
            q_upper_train, q_lower_validation, q_upper_validation)))) {
            stop("The external threshold provider returned invalid predictions.")
        }
        lower_threshold_fit <- make_external_threshold_fit(provider_thresholds$lower_predict, q_lower_train,
            0x1p+0 - tau_train)
        upper_threshold_fit <- make_external_threshold_fit(provider_thresholds$upper_predict, q_upper_train,
            tau_train)
    }
    fixed_threshold_target_evaluation <- if (is.null(fixed_threshold_target_evaluator)) {
        data.frame()
    }
    else {
        value <- fixed_threshold_target_evaluator(lower_predict = lower_threshold_fit$predict, upper_predict = upper_threshold_fit$predict,
            treatment = as.integer(t_value), mode = mode)
        if (!is.data.frame(value)) {
            stop("fixed_threshold_target_evaluator must return a data frame.")
        }
        reserved_columns <- intersect(names(value), c("fold", "treatment", "mode"))
        if (length(reserved_columns) > 0L) {
            stop("fixed_threshold_target_evaluator returned reserved columns: ", paste(reserved_columns,
                collapse = ", "), ".")
        }
        value
    }
    train_contributions <- evaluate_threshold_bound_contributions(data = train, nuisance = nuisance_train,
        envelope = envelope_train, t_value = t_value, tau = tau_train, q_lower = q_lower_train, q_upper = q_upper_train)
    validation_contributions <- evaluate_threshold_bound_contributions(data = validation, nuisance = nuisance_validation,
        envelope = envelope_validation, t_value = t_value, tau = tau_validation, q_lower = q_lower_validation,
        q_upper = q_upper_validation)
    lambda_lower <- validation_contributions$lambda_lower
    lambda_upper <- validation_contributions$lambda_upper
    contribution_lower <- validation_contributions$direct_lower
    contribution_upper <- validation_contributions$direct_upper
    contribution_lower_pseudo <- validation_contributions$pseudo_lower
    contribution_upper_pseudo <- validation_contributions$pseudo_upper
    contribution_lambda_one <- validation_contributions$lambda_one
    prediction_table <- data.frame(row_id = row_id_validation[baseline_validation], treatment = t_value,
        Y = validation$Y[baseline_validation], X = validation[baseline_validation, x_cols[0x1p+0], drop = TRUE],
        L = envelope_validation$L[baseline_validation], U = envelope_validation$U[baseline_validation],
        D = envelope_validation$D[baseline_validation], hard_L = envelope_validation$hard_L[baseline_validation],
        hard_U = envelope_validation$hard_U[baseline_validation], lower_smoothing_gap = envelope_validation$lower_smoothing_gap[baseline_validation],
        upper_smoothing_gap = envelope_validation$upper_smoothing_gap[baseline_validation], Lbar = Lbar_validation[baseline_validation],
        Ubar = Ubar_validation[baseline_validation], tau_raw = tau_validation_raw[baseline_validation],
        tau = tau_validation[baseline_validation], q_lower = q_lower_validation[baseline_validation],
        q_upper = q_upper_validation[baseline_validation], lambda_lower = lambda_lower[baseline_validation],
        lambda_upper = lambda_upper[baseline_validation], direct_lower_contribution = contribution_lower[baseline_validation],
        direct_upper_contribution = contribution_upper[baseline_validation], pseudo_lower_contribution = contribution_lower_pseudo[baseline_validation],
        pseudo_upper_contribution = contribution_upper_pseudo[baseline_validation], lower_q_score_residual = validation_contributions$lower_score_residual[baseline_validation],
        upper_q_score_residual = validation_contributions$upper_score_residual[baseline_validation],
        threshold_learner = threshold_learner, threshold_rkhs_backend = lower_threshold_fit$rkhs_backend,
        gamma_original = envelope_validation$gamma_original[baseline_validation], gamma_used = envelope_validation$gamma_used[baseline_validation],
        gamma_compatible_lower = envelope_validation$gamma_compatible_lower[baseline_validation], gamma_compatible_upper = envelope_validation$gamma_compatible_upper[baseline_validation],
        gamma_projected = envelope_validation$gamma_projected[baseline_validation])
    eif_fit <- if (compute_eif) {
        eif_arm_train <- train$T_obs == t_value
        eif_arm_validation <- validation$T_obs == t_value
        eif_gamma_projection_rate_train <- mean(envelope_train$gamma_projected[eif_arm_train])
        eif_gamma_projection_rate_validation <- mean(envelope_validation$gamma_projected[eif_arm_validation])
        eif_projection_active <- eif_gamma_projection_rate_train > 0x0p+0 || eif_gamma_projection_rate_validation >
            0x0p+0
        if (eif_projection_active && !allow_projected_eif) {
            stop(paste0("The EIF fit encountered projected gamma predictions. ", "Set allow_projected_eif=TRUE only when retaining the ",
                "projection rate as an inference diagnostic."))
        }
        result <- fit_eif_arm_fold(train = train, validation = validation, row_id_validation = row_id_validation,
            nuisance_train = nuisance_train, nuisance_validation = nuisance_validation, envelope_train = envelope_train,
            envelope_validation = envelope_validation, t_value = t_value, Lambda0 = Lambda0, Lambda1 = Lambda1,
            q_lower_train = q_lower_train, q_upper_train = q_upper_train, q_lower_validation = q_lower_validation,
            q_upper_validation = q_upper_validation, x_cols = x_cols, conditional_spline_df = eif_conditional_spline_df,
            conditional_learner = eif_conditional_learner, conditional_control = eif_conditional_control,
            conditional_seed = eif_conditional_seed, conditional_prediction_provider = eif_conditional_prediction_provider,
            instrument_probability_validation = eif_instrument_probability_validation)
        result$diagnostics$gamma_projection_active <- eif_projection_active
        result$diagnostics$gamma_projection_rate_train <- eif_gamma_projection_rate_train
        result$diagnostics$gamma_projection_rate_validation <- eif_gamma_projection_rate_validation
        result$diagnostics$conditional_projection_rate_train <- projection_rate_train
        result$diagnostics$conditional_projection_rate_validation <- projection_rate_validation
        result$diagnostics$projected_nuisance_inference_flag <- eif_projection_active
        result
    }
    else {
        NULL
    }
    threshold_cv <- list()
    if (nrow(lower_threshold_fit$cv_table) > 0L) {
        local <- lower_threshold_fit$cv_table
        local$endpoint <- "lower"
        local$learner <- threshold_learner
        threshold_cv[[length(threshold_cv) + 1L]] <- local
    }
    if (nrow(upper_threshold_fit$cv_table) > 0L) {
        local <- upper_threshold_fit$cv_table
        local$endpoint <- "upper"
        local$learner <- threshold_learner
        threshold_cv[[length(threshold_cv) + 1L]] <- local
    }
    threshold_cv <- if (length(threshold_cv) > 0L) {
        do.call(rbind, threshold_cv)
    }
    else {
        data.frame()
    }
    threshold_iteration_count <- function(fit) {
        if ("iterations" %in% names(fit$counts)) {
            return(unname(fit$counts[["iterations"]]))
        }
        if ("function" %in% names(fit$counts)) {
            return(unname(fit$counts[["function"]]))
        }
        NA_real_
    }
    threshold_fit_value <- function(fit, name, default) {
        value <- fit[[name]]
        if (is.null(value) || length(value) != 1L)
            default
        else value
    }
    lower_identity_error <- max(abs(contribution_lower - contribution_lower_pseudo - validation_contributions$lower_score_residual))
    upper_identity_error <- max(abs(contribution_upper - contribution_upper_pseudo - validation_contributions$upper_score_residual))
    list(lower_contribution = contribution_lower, upper_contribution = contribution_upper, lower_pseudo_contribution = contribution_lower_pseudo,
        upper_pseudo_contribution = contribution_upper_pseudo, lambda_one_contribution = contribution_lambda_one,
        predictions = prediction_table, diagnostics = data.frame(treatment = t_value, n_train_baseline = sum(baseline_train),
            n_validation_baseline = sum(baseline_validation), pointwise_crossing_rate_train = mean(pointwise_crossing_train[baseline_train]),
            pointwise_crossing_rate_validation = mean(pointwise_crossing_validation[baseline_validation]),
            gamma_incompatibility_rate_train = gamma_incompatibility_rate_train, gamma_incompatibility_rate_validation = gamma_incompatibility_rate_validation,
            gamma_projection_rate_train = gamma_projection_rate_train, gamma_projection_rate_validation = gamma_projection_rate_validation,
            gamma_projection_threshold_exceeded = gamma_projection_threshold_exceeded, envelope_smoothing_epsilon = envelope_smoothing_epsilon,
            mean_lower_smoothing_gap_validation = mean(envelope_validation$lower_smoothing_gap[baseline_validation]),
            mean_upper_smoothing_gap_validation = mean(envelope_validation$upper_smoothing_gap[baseline_validation]),
            maximum_lower_smoothing_gap_validation = max(envelope_validation$lower_smoothing_gap[baseline_validation]),
            maximum_upper_smoothing_gap_validation = max(envelope_validation$upper_smoothing_gap[baseline_validation]),
            conditional_projection_rate_train = projection_rate_train, conditional_projection_rate_validation = projection_rate_validation,
            tau_invalid_rate_train = mean(tau_invalid_train[baseline_train]), tau_invalid_rate_validation = mean(tau_invalid_validation[baseline_validation]),
            lower_threshold_convergence = lower_threshold_fit$convergence, upper_threshold_convergence = upper_threshold_fit$convergence,
            lower_threshold_objective = lower_threshold_fit$objective, upper_threshold_objective = upper_threshold_fit$objective,
            mean_L_validation = mean(envelope_validation$L[baseline_validation]), mean_U_validation = mean(envelope_validation$U[baseline_validation]),
            mean_tau_validation = mean(tau_validation[baseline_validation]), mean_lambda_lower_validation = mean(lambda_lower[baseline_validation]),
            mean_lambda_upper_validation = mean(lambda_upper[baseline_validation]), threshold_learner = threshold_learner,
            threshold_rkhs_backend = lower_threshold_fit$rkhs_backend, lower_selected_spline_df = lower_threshold_fit$selected_spline_df,
            upper_selected_spline_df = upper_threshold_fit$selected_spline_df, lower_selected_lambda = lower_threshold_fit$selected_lambda,
            upper_selected_lambda = upper_threshold_fit$selected_lambda, lower_selected_bandwidth_multiplier = lower_threshold_fit$selected_bandwidth_multiplier,
            upper_selected_bandwidth_multiplier = upper_threshold_fit$selected_bandwidth_multiplier,
            lower_selected_bandwidth = lower_threshold_fit$selected_bandwidth, upper_selected_bandwidth = upper_threshold_fit$selected_bandwidth,
            lower_base_bandwidth = lower_threshold_fit$base_bandwidth, upper_base_bandwidth = upper_threshold_fit$base_bandwidth,
            lower_selected_rank = lower_threshold_fit$selected_rank, upper_selected_rank = upper_threshold_fit$selected_rank,
            lower_threshold_iterations = threshold_iteration_count(lower_threshold_fit), upper_threshold_iterations = threshold_iteration_count(upper_threshold_fit),
            lower_threshold_initial_convergence = threshold_fit_value(lower_threshold_fit, "initial_convergence",
                NA_integer_), upper_threshold_initial_convergence = threshold_fit_value(upper_threshold_fit,
                "initial_convergence", NA_integer_), lower_threshold_retry_used = threshold_fit_value(lower_threshold_fit,
                "retry_used", NA), upper_threshold_retry_used = threshold_fit_value(upper_threshold_fit,
                "retry_used", NA), lower_threshold_retry_accepted = threshold_fit_value(lower_threshold_fit,
                "retry_accepted", NA), upper_threshold_retry_accepted = threshold_fit_value(upper_threshold_fit,
                "retry_accepted", NA), lower_threshold_retry_convergence = threshold_fit_value(lower_threshold_fit,
                "retry_convergence", NA_integer_), upper_threshold_retry_convergence = threshold_fit_value(upper_threshold_fit,
                "retry_convergence", NA_integer_), lower_threshold_initial_penalized_objective = threshold_fit_value(lower_threshold_fit,
                "initial_penalized_objective", NA_real_), upper_threshold_initial_penalized_objective = threshold_fit_value(upper_threshold_fit,
                "initial_penalized_objective", NA_real_), lower_threshold_penalized_objective = lower_threshold_fit$penalized_objective,
            upper_threshold_penalized_objective = upper_threshold_fit$penalized_objective, lower_threshold_tie_break_ridge = lower_threshold_fit$tie_break_ridge,
            upper_threshold_tie_break_ridge = upper_threshold_fit$tie_break_ridge, lower_threshold_tie_break_ridge_numerator = lower_threshold_fit$tie_break_ridge_numerator,
            upper_threshold_tie_break_ridge_numerator = upper_threshold_fit$tie_break_ridge_numerator,
            lower_threshold_tie_break_ridge_training_size = lower_threshold_fit$tie_break_ridge_training_size,
            upper_threshold_tie_break_ridge_training_size = upper_threshold_fit$tie_break_ridge_training_size,
            lower_threshold_primal_residual = lower_threshold_fit$primal_residual, upper_threshold_primal_residual = upper_threshold_fit$primal_residual,
            lower_threshold_dual_residual = lower_threshold_fit$dual_residual, upper_threshold_dual_residual = upper_threshold_fit$dual_residual,
            lower_threshold_final_rho = lower_threshold_fit$final_rho, upper_threshold_final_rho = upper_threshold_fit$final_rho,
            lower_threshold_reciprocal_condition = lower_threshold_fit$minimum_reciprocal_condition,
            upper_threshold_reciprocal_condition = upper_threshold_fit$minimum_reciprocal_condition,
            lower_threshold_numerical_jitter = lower_threshold_fit$maximum_numerical_jitter, upper_threshold_numerical_jitter = upper_threshold_fit$maximum_numerical_jitter,
            lower_kernel_minimum_raw_eigenvalue = lower_threshold_fit$kernel_minimum_raw_eigenvalue,
            upper_kernel_minimum_raw_eigenvalue = upper_threshold_fit$kernel_minimum_raw_eigenvalue,
            fold_lower_mean = mean(contribution_lower), fold_upper_mean = mean(contribution_upper), fold_lambda_one_mean = mean(contribution_lambda_one),
            fold_lower_pseudo_mean = mean(contribution_lower_pseudo), fold_upper_pseudo_mean = mean(contribution_upper_pseudo),
            lower_direct_minus_pseudo_train = mean(train_contributions$direct_lower - train_contributions$pseudo_lower),
            upper_direct_minus_pseudo_train = mean(train_contributions$direct_upper - train_contributions$pseudo_upper),
            lower_direct_minus_pseudo_validation = mean(contribution_lower - contribution_lower_pseudo),
            upper_direct_minus_pseudo_validation = mean(contribution_upper - contribution_upper_pseudo),
            lower_q_score_residual_train = mean(train_contributions$lower_score_residual), upper_q_score_residual_train = mean(train_contributions$upper_score_residual),
            lower_q_score_residual_validation = mean(validation_contributions$lower_score_residual),
            upper_q_score_residual_validation = mean(validation_contributions$upper_score_residual),
            lower_direct_pseudo_identity_error_validation = lower_identity_error, upper_direct_pseudo_identity_error_validation = upper_identity_error),
        threshold_cv = threshold_cv, fixed_threshold_target_evaluation = fixed_threshold_target_evaluation,
        eif = eif_fit)
}

fit_eif_arm_fold <- function (train, validation, row_id_validation, nuisance_train, nuisance_validation, envelope_train,
    envelope_validation, t_value, Lambda0, Lambda1, q_lower_train, q_upper_train, q_lower_validation,
    q_upper_validation, x_cols = "X", conditional_spline_df = 0x1.4p+2, conditional_learner = c("spline_fixed",
        "spline_cv", "gam", "ranger", "library_cv"), conditional_control = NULL, conditional_seed = 1L,
    conditional_prediction_provider = NULL, instrument_probability_validation = NULL)
{
    conditional_learner <- match.arg(conditional_learner)
    if (!is.null(conditional_prediction_provider) && !is.function(conditional_prediction_provider)) {
        stop("conditional_prediction_provider must be NULL or a function.")
    }
    resolved_conditional_control <- resolve_eif_conditional_control(conditional_control)
    if (is.null(instrument_probability_validation)) {
        instrument_probability_validation <- nuisance_validation$e
    }
    if (length(instrument_probability_validation) != nrow(validation) || any(!is.finite(instrument_probability_validation)) ||
        any(instrument_probability_validation <= 0x0p+0) || any(instrument_probability_validation >=
        0x1p+0)) {
        stop(paste0("instrument_probability_validation must contain one finite value ", "strictly between zero and one per validation row."))
    }
    primitive_train <- compute_eif_row_primitives(data = train, nuisance = nuisance_train, envelope = envelope_train,
        t_value = t_value, Lambda0 = Lambda0, Lambda1 = Lambda1, q_lower = q_lower_train, q_upper = q_upper_train)
    primitive_validation <- compute_eif_row_primitives(data = validation, nuisance = nuisance_validation,
        envelope = envelope_validation, t_value = t_value, Lambda0 = Lambda0, Lambda1 = Lambda1, q_lower = q_lower_validation,
        q_upper = q_upper_validation)
    baseline_train <- train$T_obs == t_value & train$Z == 0x0p+0
    comparison_train <- train$T_obs == t_value & train$Z == 0x1p+0
    learner_complexity <- switch(conditional_learner, spline_fixed = conditional_spline_df, spline_cv = max(resolved_conditional_control$spline_df_grid),
        gam = resolved_conditional_control$gam_k, ranger = 2L, library_cv = max(vapply(resolved_conditional_control$library_candidates,
            function(candidate) switch(candidate, linear = 1L, ranger = 2L, gam = resolved_conditional_control$gam_k),
            integer(0x1p+0))))
    minimum_size <- learner_complexity * length(x_cols) + 3L
    if (sum(baseline_train) <= minimum_size || sum(comparison_train) <= minimum_size) {
        stop(sprintf(paste0("Too few observations to fit EIF conditional means for ", "treatment %d. Baseline=%d; comparison=%d."),
            t_value, sum(baseline_train), sum(comparison_train)))
    }
    baseline_data <- train[baseline_train, , drop = FALSE]
    baseline_data$.m0 <- baseline_data$Y
    for (endpoint in c("lower", "upper")) {
        primitive <- primitive_train[[endpoint]][baseline_train, , drop = FALSE]
        baseline_data[[paste0(".theta_", endpoint)]] <- primitive$lambda * baseline_data$Y
        baseline_data[[paste0(".D_r0_", endpoint)]] <- primitive$D_r0
        baseline_data[[paste0(".D_r1_", endpoint)]] <- primitive$D_r1
        baseline_data[[paste0(".H0_", endpoint)]] <- primitive$H0
    }
    baseline_outcomes <- c(".m0", unlist(lapply(c("lower", "upper"), function(endpoint) {
        paste0(c(".theta_", ".D_r0_", ".D_r1_", ".H0_"), endpoint)
    })))
    baseline_fit <- if (is.null(conditional_prediction_provider)) {
        fit_eif_conditional_means(data = baseline_data, outcomes = baseline_outcomes, x_cols = x_cols,
            learner = conditional_learner, spline_df = conditional_spline_df, control = resolved_conditional_control,
            seed = conditional_seed)
    }
    else {
        NULL
    }
    comparison_data <- train[comparison_train, , drop = FALSE]
    for (endpoint in c("lower", "upper")) {
        comparison_data[[paste0(".H1_", endpoint)]] <- primitive_train[[endpoint]]$H1[comparison_train]
    }
    comparison_outcomes <- paste0(".H1_", c("lower", "upper"))
    comparison_fit <- if (is.null(conditional_prediction_provider)) {
        fit_eif_conditional_means(data = comparison_data, outcomes = comparison_outcomes, x_cols = x_cols,
            learner = conditional_learner, spline_df = conditional_spline_df, control = resolved_conditional_control,
            seed = conditional_seed + 1L)
    }
    else {
        NULL
    }
    baseline_prediction <- if (is.null(baseline_fit)) {
        NULL
    }
    else {
        baseline_fit$predict(validation)
    }
    comparison_prediction <- if (is.null(comparison_fit)) {
        NULL
    }
    else {
        comparison_fit$predict(validation)
    }
    provider_output <- NULL
    if (!is.null(conditional_prediction_provider)) {
        provider_output <- conditional_prediction_provider(validation = validation, nuisance_validation = nuisance_validation,
            envelope_validation = envelope_validation, t_value = t_value, Lambda0 = Lambda0, Lambda1 = Lambda1,
            q_lower_validation = q_lower_validation, q_upper_validation = q_upper_validation, baseline_outcomes = baseline_outcomes,
            comparison_outcomes = comparison_outcomes)
        if (!is.list(provider_output) || !all(c("baseline_prediction", "comparison_prediction") %in%
            names(provider_output))) {
            stop(paste0("conditional_prediction_provider must return baseline_prediction ", "and comparison_prediction matrices."))
        }
        baseline_prediction <- as.matrix(provider_output$baseline_prediction)
        comparison_prediction <- as.matrix(provider_output$comparison_prediction)
        if (nrow(baseline_prediction) != nrow(validation) || nrow(comparison_prediction) != nrow(validation) ||
            !all(baseline_outcomes %in% colnames(baseline_prediction)) || !all(comparison_outcomes %in%
            colnames(comparison_prediction)) || any(!is.finite(baseline_prediction[, baseline_outcomes])) ||
            any(!is.finite(comparison_prediction[, comparison_outcomes]))) {
            stop(paste0("The conditional prediction provider returned invalid dimensions, ", "columns, or non-finite values."))
        }
    }
    observed_arm <- as.numeric(validation$T_obs == t_value)
    instrument0 <- as.numeric(validation$Z == 0x0p+0)
    instrument1 <- 0x1p+0 - instrument0
    make_signal <- function(endpoint) {
        primitive <- primitive_validation[[endpoint]]
        m0 <- baseline_prediction[, ".m0"]
        theta_direct <- baseline_prediction[, paste0(".theta_", endpoint)]
        Hbar0 <- baseline_prediction[, paste0(".H0_", endpoint)]
        Hbar1 <- comparison_prediction[, paste0(".H1_", endpoint)]
        fixed_q <- compute_fixed_q_linked_components(m0 = m0, Hbar0 = Hbar0, Hbar1 = Hbar1, q = primitive$q,
            r0 = primitive_validation$r0, mean_D_r0 = baseline_prediction[, paste0(".D_r0_", endpoint)])
        M <- fixed_q$M
        C0 <- fixed_q$C0
        C1 <- (0x1p+0 - primitive_validation$r0) * baseline_prediction[, paste0(".D_r1_", endpoint)]
        treatment_correction <- instrument0/instrument_probability_validation * C0 * (observed_arm -
            primitive_validation$r0) + instrument1/(0x1p+0 - instrument_probability_validation) * C1 *
            (observed_arm - primitive_validation$r1)
        outcome_correction <- instrument0/instrument_probability_validation * observed_arm/primitive_validation$r0 *
            (primitive$H0 - Hbar0) + instrument1/(0x1p+0 - instrument_probability_validation) * observed_arm/primitive_validation$r1 *
            (primitive$H1 - Hbar1)
        data.frame(q = primitive$q, M = M, theta_direct = theta_direct, theta_fixed_q = fixed_q$theta_fixed_q,
            C0 = C0, C1 = C1, Hbar0 = Hbar0, Hbar1 = Hbar1, treatment_correction = treatment_correction,
            outcome_correction = outcome_correction, signal = M + treatment_correction + outcome_correction)
    }
    lower_signal <- make_signal("lower")
    upper_signal <- make_signal("upper")
    diagnostic <- data.frame(treatment = t_value, active_lower_rate_train = mean(primitive_train$active_lower),
        active_upper_rate_train = mean(primitive_train$active_upper), active_lower_rate_validation = mean(primitive_validation$active_lower),
        active_upper_rate_validation = mean(primitive_validation$active_upper), lower_mean_treatment_correction = mean(lower_signal$treatment_correction),
        lower_mean_outcome_correction = mean(lower_signal$outcome_correction), upper_mean_treatment_correction = mean(upper_signal$treatment_correction),
        upper_mean_outcome_correction = mean(upper_signal$outcome_correction), lower_mean_absolute_fixed_q_theta_adjustment = mean(abs(lower_signal$theta_fixed_q -
            lower_signal$theta_direct)), upper_mean_absolute_fixed_q_theta_adjustment = mean(abs(upper_signal$theta_fixed_q -
            upper_signal$theta_direct)), baseline_regression_rank = if (is.null(provider_output)) {
            baseline_fit$rank
        }
        else {
            NA_real_
        }, comparison_regression_rank = if (is.null(provider_output)) {
            comparison_fit$rank
        }
        else {
            NA_real_
        }, conditional_learner = if (!is.null(provider_output) && !is.null(provider_output$learner)) {
            as.character(provider_output$learner)[1L]
        }
        else {
            conditional_learner
        }, baseline_selected_spline_df = if (is.null(provider_output)) {
            baseline_fit$selected_spline_df
        }
        else {
            NA_integer_
        }, comparison_selected_spline_df = if (is.null(provider_output)) {
            comparison_fit$selected_spline_df
        }
        else {
            NA_integer_
        }, baseline_effective_df = if (is.null(provider_output)) {
            baseline_fit$effective_df
        }
        else {
            NA_real_
        }, comparison_effective_df = if (is.null(provider_output)) {
            comparison_fit$effective_df
        }
        else {
            NA_real_
        }, mean_instrument_probability_validation = mean(instrument_probability_validation), mean_absolute_instrument_probability_perturbation = mean(abs(instrument_probability_validation -
            nuisance_validation$e)), maximum_absolute_instrument_probability_perturbation = max(abs(instrument_probability_validation -
            nuisance_validation$e)))
    if (!is.null(provider_output) && !is.null(provider_output$diagnostics)) {
        provider_diagnostics <- as.data.frame(provider_output$diagnostics)
        if (nrow(provider_diagnostics) != 1L) {
            stop("Provider diagnostics must contain exactly one row.")
        }
        duplicate_names <- intersect(names(diagnostic), names(provider_diagnostics))
        if (length(duplicate_names) > 0L) {
            stop("Provider diagnostics duplicate core columns: ", paste(duplicate_names, collapse = ", "))
        }
        diagnostic <- cbind(diagnostic, provider_diagnostics)
    }
    list(signals = data.frame(row_id = row_id_validation, treatment = t_value, r0 = primitive_validation$r0,
        lower_q = lower_signal$q, lower_M = lower_signal$M, lower_theta_direct = lower_signal$theta_direct,
        lower_theta_fixed_q = lower_signal$theta_fixed_q, lower_C0 = lower_signal$C0, lower_Hbar0 = lower_signal$Hbar0,
        lower_Hbar1 = lower_signal$Hbar1, lower_treatment_correction = lower_signal$treatment_correction,
        lower_outcome_correction = lower_signal$outcome_correction, lower_signal = lower_signal$signal,
        upper_q = upper_signal$q, upper_M = upper_signal$M, upper_theta_direct = upper_signal$theta_direct,
        upper_theta_fixed_q = upper_signal$theta_fixed_q, upper_C0 = upper_signal$C0, upper_Hbar0 = upper_signal$Hbar0,
        upper_Hbar1 = upper_signal$Hbar1, upper_treatment_correction = upper_signal$treatment_correction,
        upper_outcome_correction = upper_signal$outcome_correction, upper_signal = upper_signal$signal),
        diagnostics = diagnostic, conditional_cv = if (is.null(provider_output)) {
            bind_conditional_cv_tables(baseline_fit$cv_table, comparison_fit$cv_table)
        } else {
            data.frame()
        })
}

fit_eif_conditional_library_cv <- function (data, outcomes, x_cols = "X", candidates = c("linear", "gam", "ranger"), cv_folds = 3L, selection_rule = c("one_se",
    "minimum"), gam_k = 6L, gam_method = "REML", ranger_num_trees = 500L, ranger_min_node_size = NULL,
    seed = 1L)
{
    selection_rule <- match.arg(selection_rule)
    candidates <- unique(as.character(candidates))
    allowed <- c("linear", "gam", "ranger")
    if (length(candidates) < 1L || any(!candidates %in% allowed)) {
        stop("Conditional EIF library candidates must be drawn from: ", paste(allowed, collapse = ", "),
            ".")
    }
    cv_folds <- as.integer(cv_folds)
    if (length(cv_folds) != 1L || is.na(cv_folds) || cv_folds < 2L || nrow(data) < cv_folds) {
        stop("Conditional EIF library CV requires at least two nonempty folds.")
    }
    missing_outcomes <- setdiff(outcomes, names(data))
    if (length(missing_outcomes) > 0L) {
        stop("Missing conditional EIF outcomes: ", paste(missing_outcomes, collapse = ", "))
    }
    fit_candidate <- function(training, candidate, local_seed) {
        switch(candidate, linear = {
            fit <- fit_spline_means(data = training, outcomes = outcomes, x_cols = x_cols, spline_df = 1L)
            fit$learner <- "linear"
            fit$effective_df <- fit$rank
            fit
        }, gam = fit_gam_means(data = training, outcomes = outcomes, x_cols = x_cols, basis_dimension = gam_k,
            method = gam_method), ranger = fit_ranger_means(data = training, outcomes = outcomes, x_cols = x_cols,
            num_trees = ranger_num_trees, min_node_size = ranger_min_node_size, seed = local_seed))
    }
    set.seed(as.integer(seed))
    fold_id <- sample(rep(seq_len(cv_folds), length.out = nrow(data)))
    outcome_matrix <- as.matrix(data[, outcomes, drop = FALSE])
    storage.mode(outcome_matrix) <- "double"
    outcome_scale <- apply(outcome_matrix, 2L, stats::sd)
    outcome_scale[!is.finite(outcome_scale) | outcome_scale < 0x1.5798ee2308c3ap-27] <- 0x1p+0
    risk_rows <- list()
    counter <- 1L
    for (candidate_index in seq_along(candidates)) {
        candidate <- candidates[candidate_index]
        for (fold in seq_len(cv_folds)) {
            validation <- fold_id == fold
            fit <- fit_candidate(data[!validation, , drop = FALSE], candidate, as.integer(seed + 10000L *
                candidate_index + fold))
            prediction <- fit$predict(data[validation, , drop = FALSE])
            residual <- sweep(outcome_matrix[validation, , drop = FALSE] - prediction, 2L, outcome_scale,
                "/")
            for (outcome in outcomes) {
                risk_rows[[counter]] <- data.frame(outcome = outcome, candidate = candidate, fold = fold,
                  risk = mean(residual[, outcome]^0x1p+1), observations = sum(validation), stringsAsFactors = FALSE)
                counter <- counter + 1L
            }
        }
    }
    fold_risk <- do.call(rbind, risk_rows)
    groups <- split(fold_risk, interaction(fold_risk$outcome, fold_risk$candidate, drop = TRUE, lex.order = TRUE))
    cv_table <- do.call(rbind, lapply(groups, function(group) {
        data.frame(outcome = group$outcome[1L], candidate = group$candidate[1L], mean_risk = mean(group$risk),
            standard_error = if (nrow(group) > 1L) {
                stats::sd(group$risk)/sqrt(nrow(group))
            }
            else {
                0x0p+0
            }, stringsAsFactors = FALSE)
    }))
    rownames(cv_table) <- NULL
    complexity <- c(linear = 1L, gam = 2L, ranger = 3L)
    selected <- stats::setNames(character(length(outcomes)), outcomes)
    cv_table$selected <- FALSE
    for (outcome in outcomes) {
        local <- which(cv_table$outcome == outcome)
        minimum <- local[which.min(cv_table$mean_risk[local])]
        choice <- minimum
        if (selection_rule == "one_se") {
            cutoff <- cv_table$mean_risk[minimum] + cv_table$standard_error[minimum]
            eligible <- local[cv_table$mean_risk[local] <= cutoff]
            choice <- eligible[which.min(complexity[cv_table$candidate[eligible]])]
        }
        cv_table$selected[choice] <- TRUE
        selected[outcome] <- cv_table$candidate[choice]
    }
    selected_candidates <- unique(unname(selected))
    final_fits <- stats::setNames(lapply(seq_along(selected_candidates), function(index) {
        fit_candidate(data, selected_candidates[index], as.integer(seed + 900000L + index))
    }), selected_candidates)
    list(predict = function(newdata) {
        candidate_predictions <- lapply(final_fits, function(fit) {
            fit$predict(newdata)
        })
        output <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(outcomes), dimnames = list(NULL,
            outcomes))
        for (outcome in outcomes) {
            output[, outcome] <- candidate_predictions[[selected[[outcome]]]][, outcome]
        }
        output
    }, fits = final_fits, rank = NA_real_, outcomes = outcomes, learner = "library_cv", selected_learners = selected,
        selected_spline_df = NA_integer_, effective_df = NA_real_, cv_table = cv_table, configuration = list(candidates = candidates,
            cv_folds = cv_folds, selection_rule = selection_rule, selected_learners = selected))
}

fit_eif_conditional_means <- function (data, outcomes, x_cols = "X", learner = c("spline_fixed", "spline_cv", "gam", "ranger", "library_cv"),
    spline_df = 5L, control = NULL, seed = 1L)
{
    learner <- match.arg(learner)
    control <- resolve_eif_conditional_control(control)
    if (learner == "spline_fixed") {
        fit <- fit_spline_means(data = data, outcomes = outcomes, x_cols = x_cols, spline_df = spline_df)
        fit$learner <- learner
        fit$selected_spline_df <- as.integer(spline_df)
        fit$effective_df <- fit$rank
        fit$cv_table <- data.frame()
        return(fit)
    }
    if (learner == "spline_cv") {
        return(fit_spline_means_cv(data = data, outcomes = outcomes, x_cols = x_cols, spline_df_grid = control$spline_df_grid,
            cv_folds = control$cv_folds, selection_rule = control$selection_rule, seed = seed))
    }
    if (learner == "gam") {
        return(fit_gam_means(data = data, outcomes = outcomes, x_cols = x_cols, basis_dimension = control$gam_k,
            method = control$gam_method))
    }
    if (learner == "library_cv") {
        return(fit_eif_conditional_library_cv(data = data, outcomes = outcomes, x_cols = x_cols, candidates = control$library_candidates,
            cv_folds = control$cv_folds, selection_rule = control$library_selection_rule, gam_k = control$gam_k,
            gam_method = control$gam_method, ranger_num_trees = control$ranger_num_trees, ranger_min_node_size = control$ranger_min_node_size,
            seed = seed))
    }
    fit_ranger_means(data = data, outcomes = outcomes, x_cols = x_cols, num_trees = control$ranger_num_trees,
        min_node_size = control$ranger_min_node_size, seed = seed)
}

fit_gam_means <- function (data, outcomes, x_cols = "X", basis_dimension = 6L, method = "REML")
{
    if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("The optional conditional EIF learner 'gam' requires mgcv.")
    }
    non_numeric <- x_cols[!vapply(data[, x_cols, drop = FALSE], is.numeric, logical(0x1p+0))]
    if (length(non_numeric) > 0L) {
        stop("Conditional EIF GAM predictors must be numeric: ", paste(non_numeric, collapse = ", "))
    }
    unique_counts <- vapply(data[, x_cols, drop = FALSE], function(value) length(unique(value[is.finite(value)])),
        integer(0x1p+0))
    smooth_cols <- x_cols[unique_counts >= as.integer(basis_dimension)]
    linear_cols <- setdiff(x_cols, smooth_cols)
    smooth_terms <- sprintf("s(`%s`, bs='cr', k=%d)", gsub("`", "", smooth_cols, fixed = TRUE), as.integer(basis_dimension))
    linear_terms <- sprintf("`%s`", gsub("`", "", linear_cols, fixed = TRUE))
    model_terms <- c(linear_terms, smooth_terms)
    if (!length(model_terms))
        model_terms <- "1"
    fits <- stats::setNames(lapply(outcomes, function(outcome) {
        response <- data[[outcome]]
        center <- mean(response)
        scale <- stats::sd(response)
        if (!is.finite(scale) || scale < 0x1.b7cdfd9d7bdbbp-34) {
            return(list(fit = NULL, center = center, scale = 0x0p+0, rank = 0x1p+0, effective_df = 0x1p+0))
        }
        local_data <- data
        local_data$.eif_gam_response <- (response - center)/scale
        formula <- stats::as.formula(paste0(".eif_gam_response ~ ", paste(model_terms, collapse = " + ")))
        environment(formula) <- asNamespace("mgcv")
        fit <- tryCatch(mgcv::gam(formula = formula, data = local_data, family = stats::gaussian(), method = method,
            select = TRUE), error = function(error) {
            stop("Conditional EIF GAM failed for outcome ", outcome, ": ", conditionMessage(error), call. = FALSE)
        })
        list(fit = fit, center = center, scale = scale, rank = fit$rank, effective_df = sum(fit$edf))
    }), outcomes)
    list(predict = function(newdata) {
        output <- vapply(fits, function(specification) {
            if (is.null(specification$fit)) {
                return(rep(specification$center, nrow(newdata)))
            }
            specification$center + specification$scale * as.numeric(stats::predict(specification$fit,
                newdata = newdata, type = "response"))
        }, numeric(nrow(newdata)))
        output <- as.matrix(output)
        colnames(output) <- outcomes
        output
    }, fits = fits, rank = sum(vapply(fits, function(specification) {
        specification$rank
    }, numeric(0x1p+0))), outcomes = outcomes, learner = "gam", selected_spline_df = NA_integer_, smooth_predictors = smooth_cols,
        linear_predictors = linear_cols, effective_df = sum(vapply(fits, function(specification) specification$effective_df,
            numeric(0x1p+0))), cv_table = data.frame())
}

fit_ranger_means <- function (data, outcomes, x_cols = "X", num_trees = 500L, min_node_size = NULL, seed = 1L)
{
    if (!requireNamespace("ranger", quietly = TRUE)) {
        stop("The optional conditional EIF learner 'ranger' requires ranger.")
    }
    non_numeric <- x_cols[!vapply(data[, x_cols, drop = FALSE], is.numeric, logical(0x1p+0))]
    if (length(non_numeric) > 0L) {
        stop("Conditional EIF ranger predictors must be numeric: ", paste(non_numeric, collapse = ", "))
    }
    num_trees <- as.integer(num_trees)
    local_min_node_size <- if (is.null(min_node_size)) {
        max(5L, as.integer(floor(nrow(data)^(0x1p+0/0x1p+2))))
    }
    else {
        as.integer(min_node_size)
    }
    if (length(num_trees) != 1L || is.na(num_trees) || num_trees < 1L || length(local_min_node_size) !=
        1L || is.na(local_min_node_size) || local_min_node_size < 1L) {
        stop("Invalid conditional EIF ranger controls.")
    }
    fits <- stats::setNames(lapply(seq_along(outcomes), function(index) {
        outcome <- outcomes[index]
        response <- as.numeric(data[[outcome]])
        if (all(abs(response - response[1L]) < 0x1.19799812dea11p-40)) {
            return(list(fit = NULL, constant = response[1L]))
        }
        ranger_data <- data[, x_cols, drop = FALSE]
        ranger_data$.eif_ranger_response <- response
        fit <- ranger::ranger(formula = stats::reformulate(x_cols, response = ".eif_ranger_response"),
            data = ranger_data, num.trees = num_trees, mtry = max(1L, floor(sqrt(length(x_cols)))), min.node.size = local_min_node_size,
            seed = as.integer(seed + index - 1L), num.threads = 1L)
        list(fit = fit, constant = NA_real_)
    }), outcomes)
    list(predict = function(newdata) {
        output <- vapply(fits, function(specification) {
            if (is.null(specification$fit)) {
                return(rep(specification$constant, nrow(newdata)))
            }
            as.numeric(predict(specification$fit, data = newdata)$predictions)
        }, numeric(nrow(newdata)))
        output <- as.matrix(output)
        colnames(output) <- outcomes
        output
    }, fits = fits, rank = length(x_cols) * length(outcomes), outcomes = outcomes, learner = "ranger",
        selected_spline_df = NA_integer_, effective_df = NA_real_, cv_table = data.frame(), num_trees = num_trees,
        min_node_size = local_min_node_size)
}

fit_spline_mean <- function (data, outcome, x_cols = "X", spline_df = 0x1.4p+2)
{
    design <- make_spline_design(data = data, x_cols = x_cols, spline_df = spline_df)
    outcome_value <- data[[outcome]]
    fit <- stats::lm.fit(x = design$matrix, y = outcome_value)
    coefficients <- fit$coefficients
    coefficients[!is.finite(coefficients)] <- 0x0p+0
    list(predict = function(newdata) {
        as.vector(design$predict_matrix(newdata) %*% coefficients)
    }, coefficients = coefficients, rank = fit$rank)
}

fit_spline_means <- function (data, outcomes, x_cols = "X", spline_df = 0x1.4p+2)
{
    if (length(outcomes) == 0L) {
        stop("outcomes must contain at least one column name.")
    }
    missing_outcomes <- setdiff(outcomes, names(data))
    if (length(missing_outcomes) > 0L) {
        stop("Missing spline outcomes: ", paste(missing_outcomes, collapse = ", "))
    }
    design <- make_spline_design(data = data, x_cols = x_cols, spline_df = spline_df)
    outcome_matrix <- as.matrix(data[, outcomes, drop = FALSE])
    storage.mode(outcome_matrix) <- "double"
    fit <- stats::lm.fit(x = design$matrix, y = outcome_matrix)
    coefficients <- fit$coefficients
    if (is.null(dim(coefficients))) {
        coefficients <- matrix(coefficients, ncol = 1L, dimnames = list(NULL, outcomes))
    }
    coefficients[!is.finite(coefficients)] <- 0x0p+0
    list(predict = function(newdata) {
        output <- design$predict_matrix(newdata) %*% coefficients
        output <- as.matrix(output)
        colnames(output) <- outcomes
        output
    }, coefficients = coefficients, rank = fit$rank, outcomes = outcomes)
}

fit_spline_means_cv <- function (data, outcomes, x_cols = "X", spline_df_grid = c(3L, 5L, 7L, 9L), cv_folds = 3L, selection_rule = c("one_se",
    "minimum"), seed = 1L)
{
    selection_rule <- match.arg(selection_rule)
    spline_df_grid <- sort(unique(as.integer(spline_df_grid)))
    if (length(spline_df_grid) == 0L || any(spline_df_grid < 1L)) {
        stop("spline_df_grid must contain positive integers.")
    }
    cv_folds <- as.integer(cv_folds)
    if (cv_folds < 2L || nrow(data) < cv_folds) {
        stop("Conditional-mean CV requires at least two nonempty folds.")
    }
    set.seed(seed)
    fold_id <- sample(rep(seq_len(cv_folds), length.out = nrow(data)))
    outcome_matrix <- as.matrix(data[, outcomes, drop = FALSE])
    outcome_scale <- apply(outcome_matrix, 2L, stats::sd)
    outcome_scale[!is.finite(outcome_scale) | outcome_scale < 0x1.5798ee2308c3ap-27] <- 0x1p+0
    risk_rows <- list()
    counter <- 1L
    for (spline_df in spline_df_grid) {
        for (fold in seq_len(cv_folds)) {
            validation <- fold_id == fold
            fit <- fit_spline_means(data = data[!validation, , drop = FALSE], outcomes = outcomes, x_cols = x_cols,
                spline_df = spline_df)
            prediction <- fit$predict(data[validation, , drop = FALSE])
            residual <- sweep(outcome_matrix[validation, , drop = FALSE] - prediction, 2L, outcome_scale,
                "/")
            risk_rows[[counter]] <- data.frame(spline_df = spline_df, fold = fold, risk = mean(residual^0x1p+1),
                stringsAsFactors = FALSE)
            counter <- counter + 1L
        }
    }
    fold_risk <- do.call(rbind, risk_rows)
    cv_table <- do.call(rbind, lapply(split(fold_risk, fold_risk$spline_df), function(group) data.frame(spline_df = group$spline_df[1L],
        mean_risk = mean(group$risk), standard_error = stats::sd(group$risk)/sqrt(nrow(group)), stringsAsFactors = FALSE)))
    rownames(cv_table) <- NULL
    minimum <- which.min(cv_table$mean_risk)
    selected <- if (selection_rule == "minimum") {
        minimum
    }
    else {
        cutoff <- cv_table$mean_risk[minimum] + cv_table$standard_error[minimum]
        eligible <- which(cv_table$mean_risk <= cutoff)
        eligible[which.min(cv_table$spline_df[eligible])]
    }
    cv_table$selected <- seq_len(nrow(cv_table)) == selected
    final <- fit_spline_means(data = data, outcomes = outcomes, x_cols = x_cols, spline_df = cv_table$spline_df[selected])
    final$learner <- "spline_cv"
    final$selected_spline_df <- cv_table$spline_df[selected]
    final$effective_df <- final$rank
    final$cv_table <- cv_table
    final
}

fit_weighted_quantile_ridge_admm <- function (features, y, level, weight, lambda, rho = 0x1p+0, maxit = 1000L, abs_tol = 0x1.4f8b588e368f1p-17,
    rel_tol = 0x1.a36e2eb1c432dp-14)
{
    features <- as.matrix(features)
    design <- cbind(intercept = 0x1p+0, features)
    n <- nrow(design)
    p <- ncol(design)
    normalized_weight <- weight/sum(weight)
    penalty <- diag(c(0x0p+0, rep(lambda, p - 1L)), p)
    current_rho <- rho
    update_factor <- function() {
        system <- current_rho * crossprod(design) + penalty + diag(0x1.b7cdfd9d7bdbbp-34, p)
        chol(system)
    }
    factor <- update_factor()
    solve_system <- function(rhs) {
        backsolve(factor, forwardsolve(t(factor), rhs))
    }
    beta <- stats::lm.wfit(design, y, pmax(weight, 0x1.5798ee2308c3ap-27))$coefficients
    beta[!is.finite(beta)] <- 0x0p+0
    residual <- y - as.vector(design %*% beta)
    dual <- numeric(n)
    convergence <- 1L
    primal_norm <- dual_norm <- Inf
    for (iteration in seq_len(as.integer(maxit))) {
        residual_previous <- residual
        rhs <- current_rho * crossprod(design, y - residual - dual)
        beta <- as.vector(solve_system(rhs))
        proximal_center <- y - as.vector(design %*% beta) - dual
        threshold_positive <- normalized_weight * level/current_rho
        threshold_negative <- normalized_weight * (0x1p+0 - level)/current_rho
        residual <- ifelse(proximal_center > threshold_positive, proximal_center - threshold_positive,
            ifelse(proximal_center < -threshold_negative, proximal_center + threshold_negative, 0x0p+0))
        primal <- as.vector(design %*% beta) + residual - y
        dual <- dual + primal
        primal_norm <- sqrt(sum(primal^0x1p+1))
        dual_vector <- current_rho * crossprod(design, residual - residual_previous)
        dual_norm <- sqrt(sum(dual_vector^0x1p+1))
        primal_tolerance <- sqrt(n) * abs_tol + rel_tol * max(sqrt(sum((design %*% beta)^0x1p+1)), sqrt(sum(residual^0x1p+1)),
            sqrt(sum(y^0x1p+1)))
        dual_tolerance <- sqrt(p) * abs_tol + rel_tol * sqrt(sum((current_rho * crossprod(design, dual))^0x1p+1))
        if (primal_norm <= primal_tolerance && dual_norm <= dual_tolerance) {
            convergence <- 0L
            break
        }
        if (iteration%%25L == 0L) {
            new_rho <- current_rho
            if (primal_norm > 0x1.4p+3 * dual_norm) {
                new_rho <- min(current_rho * 0x1p+1, 0x1.9p+6)
            }
            else if (dual_norm > 0x1.4p+3 * primal_norm) {
                new_rho <- max(current_rho/0x1p+1, 0x1.4f8b588e368f1p-17)
            }
            if (new_rho != current_rho) {
                dual <- dual * current_rho/new_rho
                current_rho <- new_rho
                factor <- update_factor()
            }
        }
    }
    fitted <- as.vector(design %*% beta)
    unpenalized_risk <- sum(weight * check_loss(y, fitted, level))/sum(weight)
    penalized_objective <- unpenalized_risk + lambda * sum(beta[-1L]^0x1p+1)/0x1p+1
    list(coefficients = beta, fitted = fitted, convergence = convergence, iterations = iteration, primal_residual = primal_norm,
        dual_residual = dual_norm, final_rho = current_rho, objective = unpenalized_risk, penalized_objective = penalized_objective,
        predict_features = function(new_features) {
            as.vector(cbind(0x1p+0, as.matrix(new_features)) %*% beta)
        })
}

gauss_legendre_rule <- function (n_nodes = 0x1.4p+6, lower = -0x1p+1, upper = 0x1p+1)
{
    if (n_nodes < 0x1p+1) {
        stop("n_nodes must be at least 2.")
    }
    index <- seq_len(n_nodes - 0x1p+0)
    beta <- index/sqrt(0x1p+2 * index^0x1p+1 - 0x1p+0)
    jacobi <- matrix(0x0p+0, nrow = n_nodes, ncol = n_nodes)
    jacobi[cbind(index, index + 0x1p+0)] <- beta
    jacobi[cbind(index + 0x1p+0, index)] <- beta
    eig <- eigen(jacobi, symmetric = TRUE)
    ordering <- order(eig$values)
    nodes_standard <- eig$values[ordering]
    weights_standard <- 0x1p+1 * eig$vectors[0x1p+0, ordering]^0x1p+1
    list(nodes = (upper - lower) * nodes_standard/0x1p+1 + (upper + lower)/0x1p+1, weights = (upper -
        lower) * weights_standard/0x1p+1)
}

make_additive_spline_formula <- function (response, predictors, spline_df = 0x1.4p+2, data = NULL)
{
    terms <- vapply(predictors, function(name) {
        use_linear <- !is.null(data) && length(unique(data[[name]][is.finite(data[[name]])])) <= spline_df +
            1L
        if (use_linear)
            name
        else sprintf("splines::ns(%s, df=%d)", name, spline_df)
    }, character(0x1p+0))
    stats::as.formula(paste(response, "~", paste(terms, collapse = " + ")))
}

make_spline_design <- function (data, x_cols = "X", spline_df = 0x1.4p+2)
{
    non_numeric <- x_cols[!vapply(data[, x_cols, drop = FALSE], is.numeric, logical(0x1p+0))]
    if (length(non_numeric) > 0x0p+0) {
        stop("All x_cols must be numeric for the spline learners: ", paste(non_numeric, collapse = ", "))
    }
    basis_specs <- vector("list", length(x_cols))
    names(basis_specs) <- x_cols
    training_parts <- list(intercept = rep(0x1p+0, nrow(data)))
    for (name in x_cols) {
        unique_values <- sort(unique(data[[name]][is.finite(data[[name]])]))
        if (length(unique_values) < 2L) {
            warning("Dropping constant spline predictor in this training subset: ", name, call. = FALSE)
            basis_specs[[name]] <- list(type = "constant", training_value = if (length(unique_values)) unique_values[1L] else NA_real_)
            next
        }
        if (length(unique_values) <= spline_df + 1L) {
            basis <- matrix(data[[name]], ncol = 1L)
            basis_specs[[name]] <- list(type = "linear")
        }
        else {
            basis <- splines::ns(data[[name]], df = spline_df, intercept = FALSE)
            basis_specs[[name]] <- list(type = "spline", knots = attr(basis, "knots"), boundary_knots = attr(basis,
                "Boundary.knots"))
        }
        colnames(basis) <- paste0(name, "_s", seq_len(ncol(basis)))
        training_parts[[name]] <- basis
    }
    training_matrix <- do.call(cbind, training_parts)
    storage.mode(training_matrix) <- "double"
    predict_matrix <- function(newdata) {
        new_parts <- list(intercept = rep(0x1p+0, nrow(newdata)))
        for (name in x_cols) {
            spec <- basis_specs[[name]]
            if (identical(spec$type, "constant"))
                next
            basis <- if (identical(spec$type, "linear")) {
                matrix(newdata[[name]], ncol = 1L)
            }
            else {
                splines::ns(newdata[[name]], knots = spec$knots, Boundary.knots = spec$boundary_knots,
                  intercept = FALSE)
            }
            colnames(basis) <- paste0(name, "_s", seq_len(ncol(basis)))
            new_parts[[name]] <- basis
        }
        output <- do.call(cbind, new_parts)
        storage.mode(output) <- "double"
        output
    }
    list(matrix = training_matrix, predict_matrix = predict_matrix, basis_specs = basis_specs)
}

make_tz_stratified_folds <- function (data, K = 0x1.4p+2, seed = 0x1p+0)
{
    check_binary(data$T_obs, "T_obs")
    check_binary(data$Z, "Z")
    if (K < 0x1p+1) {
        stop("K must be at least 2.")
    }
    set.seed(seed)
    strata <- interaction(data$T_obs, data$Z, drop = TRUE)
    fold_id <- integer(nrow(data))
    for (stratum in levels(strata)) {
        index <- which(strata == stratum)
        if (length(index) < K) {
            stop(sprintf("Stratum %s has only %d observations, fewer than K=%d.", stratum, length(index),
                K))
        }
        fold_id[index] <- sample(rep(seq_len(K), length.out = length(index)))
    }
    fold_id
}

make_z_stratified_folds <- function (data, K = 0x1.8p+1, seed = 0x1p+0)
{
    if (K < 0x1p+1) {
        stop("K must be at least 2.")
    }
    set.seed(seed)
    fold_id <- integer(nrow(data))
    for (z_value in c(0x0p+0, 0x1p+0)) {
        index <- which(data$Z == z_value)
        if (length(index) < K) {
            stop(sprintf("Only %d observations with Z=%d, fewer than K=%d.", length(index), z_value,
                K))
        }
        fold_id[index] <- sample(rep(seq_len(K), length.out = length(index)))
    }
    fold_id
}

radical_inverse <- function (index, base)
{
    index <- as.integer(index)
    base <- as.integer(base)
    output <- numeric(length(index))
    fraction <- 0x1p+0/base
    remaining <- index
    while (any(remaining > 0L)) {
        output <- output + fraction * (remaining%%base)
        remaining <- remaining%/%base
        fraction <- fraction/base
    }
    output
}

resolve_eif_conditional_control <- function (control = NULL)
{
    output <- default_eif_conditional_control()
    if (!is.null(control)) {
        unknown <- setdiff(names(control), names(output))
        if (length(unknown) > 0L) {
            stop("Unknown eif_conditional_control entries: ", paste(unknown, collapse = ", "))
        }
        output[names(control)] <- control
    }
    output$spline_df_grid <- sort(unique(as.integer(output$spline_df_grid)))
    output$cv_folds <- as.integer(output$cv_folds)
    output$selection_rule <- match.arg(output$selection_rule, c("one_se", "minimum"))
    output$gam_k <- as.integer(output$gam_k)
    output$ranger_num_trees <- as.integer(output$ranger_num_trees)
    output$library_candidates <- unique(as.character(output$library_candidates))
    output$library_selection_rule <- match.arg(output$library_selection_rule, c("one_se", "minimum"))
    if (!is.null(output$ranger_min_node_size)) {
        output$ranger_min_node_size <- as.integer(output$ranger_min_node_size)
    }
    if (length(output$spline_df_grid) == 0L || anyNA(output$spline_df_grid) || any(output$spline_df_grid <
        1L) || length(output$cv_folds) != 1L || is.na(output$cv_folds) || output$cv_folds < 2L || length(output$gam_k) !=
        1L || is.na(output$gam_k) || output$gam_k < 3L || length(output$gam_method) != 1L || !nzchar(output$gam_method) ||
        length(output$ranger_num_trees) != 1L || is.na(output$ranger_num_trees) || output$ranger_num_trees <
        1L || length(output$library_candidates) < 1L || any(!output$library_candidates %in% c("linear",
        "gam", "ranger")) || (!is.null(output$ranger_min_node_size) && (length(output$ranger_min_node_size) !=
        1L || is.na(output$ranger_min_node_size) || output$ranger_min_node_size < 1L))) {
        stop("Invalid eif_conditional_control values.")
    }
    output
}

x5s_oracle <- function (data, a, b, contract = x5s_contract())
{
    x <- as.matrix(data[, contract$x_cols])
    h <- drop(x %*% contract$beta)
    location <- drop(x %*% contract$theta)
    list(e = rep(0x1p-1, nrow(data)), pi0 = plogis(h - a/0x1p+1), pi1 = plogis(h + a/0x1p+1), gamma0 = siv_gamma(data$Y -
        location, -b), gamma1 = siv_gamma(data$Y - location - contract$shift, b))
}

