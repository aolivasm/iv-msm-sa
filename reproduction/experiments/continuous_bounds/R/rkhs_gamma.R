# Main-text reproduction dependencies; published formulas retained.
compute_compatible_gamma_interval <- function(nuisance, t_value, Lambda0, Lambda1) {
    pi0 <- nuisance$pi0
    pi1 <- nuisance$pi1
    if (t_value == 1) {
        gamma_lower <- ((1 - pi0)/Lambda0 + pi0)/((1 - pi1) * Lambda1 + pi1)
        gamma_upper <- (Lambda0 * (1 - pi0) + pi0)/((1 - pi1)/Lambda1 + pi1)
    }
    else {
        gamma_lower <- (pi0/Lambda0 + 1 - pi0)/(pi1 * Lambda1 + 1 - pi1)
        gamma_upper <- (Lambda0 * pi0 + 1 - pi0)/(pi1/Lambda1 + 1 - pi1)
    }
    list(lower = gamma_lower, upper = gamma_upper)
}

compute_riesz_weight <- function(nuisance, t_value) {
    if (t_value == 1) {
        nuisance$e * nuisance$pi0/((1 - nuisance$e) * nuisance$pi1)
    }
    else {
        nuisance$e * (1 - nuisance$pi0)/((1 - nuisance$e) * (1 - nuisance$pi1))
    }
}

default_gamma_lambda_grid <- function(n_training, treatment = NULL, constants = 10^seq(-1.5, 1, by = 0.5)) {
    if (length(n_training) != 1L || !is.finite(n_training) || n_training < 1L) {
        stop("n_training must be a positive integer.")
    }
    constants <- sort(unique(as.numeric(constants)))
    if (length(constants) == 0L || any(!is.finite(constants)) || any(constants <= 0)) {
        stop("constants must contain positive finite values.")
    }
    constants/sqrt(n_training)
}

fit_gamma_rkhs <- function(train, nuisance_train, t_value, x_cols = "X", lambda_grid = default_gamma_lambda_grid,
    bandwidth_multipliers = 1/sqrt(2), cv_folds = 3, cv_repeats = 2L, selection_rule = c("minimum", "one_se"),
    gamma_floor = 1e-06, gamma_cap = 100, bandwidth_sample = 1000, cv_plan = NULL, backend = c("exact",
        "nystrom"), nystrom_rank = NULL, nystrom_eigen_tolerance = 1e-10, nystrom_jitter = 1e-10, kernel_type = c("gaussian",
        "sobolev_matern32", "sobolev_matern52"), seed = 1) {
    backend <- match.arg(backend)
    kernel_type <- match.arg(kernel_type)
    selection_rule <- match.arg(selection_rule)
    treatment_index <- which(train$T_obs == t_value)
    train_t <- train[treatment_index, , drop = FALSE]
    nuisance_t <- subset_nuisance(nuisance_train, treatment_index)
    lambda_grid <- resolve_gamma_lambda_grid(lambda_grid = lambda_grid, n_training = nrow(train_t), treatment = t_value)
    if (is.null(cv_plan)) {
        cv_plan <- make_rkhs_cv_plan(train_t = train_t, x_cols = x_cols, cv_folds = cv_folds, cv_repeats = cv_repeats,
            bandwidth_sample = bandwidth_sample, backend = backend, nystrom_rank = nystrom_rank, seed = seed)
    }
    else if (!identical(cv_plan$backend, backend)) {
        stop("cv_plan backend does not match the requested gamma backend.")
    }
    tuning <- select_rkhs_lambda(train_t = train_t, nuisance_t = nuisance_t, t_value = t_value, x_cols = x_cols,
        lambda_grid = lambda_grid, bandwidth_multipliers = bandwidth_multipliers, cv_plan = cv_plan,
        selection_rule = selection_rule, bandwidth_sample = bandwidth_sample, backend = backend, nystrom_rank = nystrom_rank,
        nystrom_eigen_tolerance = nystrom_eigen_tolerance, nystrom_jitter = nystrom_jitter, kernel_type = kernel_type)
    final_nystrom_feature_spec <- if (backend == "nystrom") {
        make_gamma_nystrom_features(coordinate_spec = cv_plan$final_coordinate_spec, landmark_spec = cv_plan$final_nystrom_landmark_spec,
            bandwidth_multiplier = tuning$selected_bandwidth_multiplier, kernel_type = kernel_type, eigen_tolerance = nystrom_eigen_tolerance,
            numerical_jitter = nystrom_jitter)
    }
    else {
        NULL
    }
    final_fit <- fit_gamma_rkhs_reduced(train_t = train_t, nuisance_t = nuisance_t, t_value = t_value,
        x_cols = x_cols, lambda = tuning$selected_lambda, bandwidth_multiplier = tuning$selected_bandwidth_multiplier,
        coordinate_spec = cv_plan$final_coordinate_spec, bandwidth_sample = bandwidth_sample, backend = backend,
        nystrom_rank = nystrom_rank, nystrom_landmark_spec = cv_plan$final_nystrom_landmark_spec, nystrom_feature_spec = final_nystrom_feature_spec,
        nystrom_eigen_tolerance = nystrom_eigen_tolerance, nystrom_jitter = nystrom_jitter, kernel_type = kernel_type)
    raw_train <- final_fit$predict_raw(train_t)
    positive_train <- clip_positive(raw_train, floor = gamma_floor, cap = gamma_cap)
    normalization_constant <- mean(positive_train[train_t$Z == 0])
    normalization_constant <- max(normalization_constant, gamma_floor)
    list(predict = function(newdata) {
        raw <- final_fit$predict_raw(newdata)
        positive <- clip_positive(raw, floor = gamma_floor, cap = gamma_cap)
        normalized <- clip_positive(positive/normalization_constant, floor = gamma_floor, cap = gamma_cap)
        list(raw = raw, positive = positive, normalized = normalized)
    }, selected_lambda = tuning$selected_lambda, selected_bandwidth_multiplier = tuning$selected_bandwidth_multiplier,
        base_bandwidth = final_fit$base_sigma, selected_bandwidth = final_fit$sigma, cv_table = tuning$cv_table,
        cv_resample_table = tuning$cv_resample_table, normalization_constant = normalization_constant,
        training_negative_rate = mean(raw_train <= 0), rkhs_backend = final_fit$backend, rkhs_kernel = final_fit$kernel_type,
        nystrom_rank_requested = final_fit$nystrom_rank_requested, nystrom_landmark_count = final_fit$nystrom_landmark_count,
        nystrom_effective_rank = final_fit$nystrom_effective_rank, nystrom_min_retained_eigenvalue = final_fit$nystrom_min_retained_eigenvalue,
        nystrom_max_eigenvalue = final_fit$nystrom_max_eigenvalue, nystrom_landmark_indices = final_fit$nystrom_landmark_indices,
        exact_solver_method = final_fit$exact_solver_method, exact_solver_minimum_eigenvalue = final_fit$exact_solver_minimum_eigenvalue,
        exact_solver_eigenvalue_cutoff = final_fit$exact_solver_eigenvalue_cutoff, cv_plan = cv_plan)
}

fit_gamma_rkhs_fixed <- function(train, nuisance_train, t_value, x_cols = "X", lambda, bandwidth_multiplier = 1/sqrt(2),
    gamma_floor = 1e-06, gamma_cap = 100, bandwidth_sample = 1000, backend = c("exact", "nystrom"), nystrom_rank = NULL,
    nystrom_eigen_tolerance = 1e-10, nystrom_jitter = 1e-10, kernel_type = c("gaussian", "sobolev_matern32",
        "sobolev_matern52"), seed = 1) {
    backend <- match.arg(backend)
    kernel_type <- match.arg(kernel_type)
    treatment_index <- which(train$T_obs == t_value)
    train_t <- train[treatment_index, , drop = FALSE]
    nuisance_t <- subset_nuisance(nuisance_train, treatment_index)
    coordinate_spec <- prepare_rkhs_coordinates(train_t = train_t, x_cols = x_cols, bandwidth_sample = bandwidth_sample,
        seed = seed)
    nystrom_landmark_spec <- if (backend == "nystrom") {
        prepare_gamma_nystrom_landmarks(coordinate_spec = coordinate_spec, rank = nystrom_rank, seed = seed +
            10000)
    }
    else {
        NULL
    }
    fit <- fit_gamma_rkhs_reduced(train_t = train_t, nuisance_t = nuisance_t, t_value = t_value, x_cols = x_cols,
        lambda = lambda, bandwidth_multiplier = bandwidth_multiplier, coordinate_spec = coordinate_spec,
        bandwidth_sample = bandwidth_sample, backend = backend, nystrom_rank = nystrom_rank, nystrom_landmark_spec = nystrom_landmark_spec,
        nystrom_eigen_tolerance = nystrom_eigen_tolerance, nystrom_jitter = nystrom_jitter, kernel_type = kernel_type,
        seed = seed)
    raw_train <- fit$predict_raw(train_t)
    positive_train <- clip_positive(raw_train, floor = gamma_floor, cap = gamma_cap)
    normalization_constant <- mean(positive_train[train_t$Z == 0])
    normalization_constant <- max(normalization_constant, gamma_floor)
    list(predict = function(newdata) {
        raw <- fit$predict_raw(newdata)
        positive <- clip_positive(raw, floor = gamma_floor, cap = gamma_cap)
        normalized <- clip_positive(positive/normalization_constant, floor = gamma_floor, cap = gamma_cap)
        list(raw = raw, positive = positive, normalized = normalized)
    }, selected_lambda = lambda, selected_bandwidth_multiplier = bandwidth_multiplier, base_bandwidth = fit$base_sigma,
        selected_bandwidth = fit$sigma, cv_table = NULL, normalization_constant = normalization_constant,
        training_negative_rate = mean(raw_train <= 0), sigma = fit$sigma, rkhs_backend = fit$backend,
        rkhs_kernel = fit$kernel_type, nystrom_rank_requested = fit$nystrom_rank_requested, nystrom_landmark_count = fit$nystrom_landmark_count,
        nystrom_effective_rank = fit$nystrom_effective_rank, nystrom_min_retained_eigenvalue = fit$nystrom_min_retained_eigenvalue,
        nystrom_max_eigenvalue = fit$nystrom_max_eigenvalue, nystrom_landmark_indices = fit$nystrom_landmark_indices,
        exact_solver_method = fit$exact_solver_method, exact_solver_minimum_eigenvalue = fit$exact_solver_minimum_eigenvalue,
        exact_solver_eigenvalue_cutoff = fit$exact_solver_eigenvalue_cutoff)
}

fit_gamma_rkhs_reduced <- function(train_t, nuisance_t, t_value, x_cols = "X", lambda, bandwidth_multiplier = 1/sqrt(2),
    coordinate_spec = NULL, bandwidth_sample = 1000, numerical_jitter = 1e-10, backend = c("exact", "nystrom"),
    nystrom_rank = NULL, nystrom_landmark_spec = NULL, nystrom_feature_spec = NULL, nystrom_eigen_tolerance = 1e-10,
    nystrom_jitter = 1e-10, kernel_type = c("gaussian", "sobolev_matern32", "sobolev_matern52"), seed = 1) {
    backend <- match.arg(backend)
    kernel_type <- match.arg(kernel_type)
    if (length(bandwidth_multiplier) != 1 || !is.finite(bandwidth_multiplier) || bandwidth_multiplier <=
        0) {
        stop("bandwidth_multiplier must be one positive finite value.")
    }
    if (is.null(coordinate_spec)) {
        coordinate_spec <- prepare_rkhs_coordinates(train_t = train_t, x_cols = x_cols, bandwidth_sample = bandwidth_sample,
            seed = seed)
    }
    selected_sigma <- coordinate_spec$sigma * bandwidth_multiplier
    W0 <- as.numeric(train_t$Z == 0)
    w1 <- as.numeric(train_t$Z == 1) * compute_riesz_weight(nuisance = nuisance_t, t_value = t_value)
    n_t <- nrow(train_t)
    if (backend == "nystrom") {
        if (is.null(nystrom_feature_spec)) {
            if (is.null(nystrom_landmark_spec)) {
                nystrom_landmark_spec <- prepare_gamma_nystrom_landmarks(coordinate_spec = coordinate_spec,
                  rank = nystrom_rank, seed = seed)
            }
            nystrom_feature_spec <- make_gamma_nystrom_features(coordinate_spec = coordinate_spec, landmark_spec = nystrom_landmark_spec,
                bandwidth_multiplier = bandwidth_multiplier, kernel_type = kernel_type, eigen_tolerance = nystrom_eigen_tolerance,
                numerical_jitter = nystrom_jitter)
        }
        features <- nystrom_feature_spec$features_train
        weighted_features <- sweep(features, 1L, W0, "*")
        system_matrix <- crossprod(features, weighted_features)
        diag(system_matrix) <- diag(system_matrix) + n_t * lambda + numerical_jitter
        right_hand_side <- as.vector(crossprod(features, w1))
        theta <- tryCatch(solve(system_matrix, right_hand_side), error = function(e) {
            qr.solve(system_matrix, right_hand_side, tol = 1e-10)
        })
        return(list(predict_raw = function(newdata) {
            as.vector(nystrom_feature_spec$predict_features(newdata) %*% theta)
        }, theta = theta, sigma = selected_sigma, base_sigma = coordinate_spec$sigma, bandwidth_multiplier = bandwidth_multiplier,
            coordinate_spec = coordinate_spec, backend = backend, kernel_type = kernel_type, nystrom_rank_requested = nystrom_feature_spec$rank_requested,
            nystrom_landmark_count = nystrom_feature_spec$landmark_count, nystrom_effective_rank = nystrom_feature_spec$effective_rank,
            nystrom_min_retained_eigenvalue = nystrom_feature_spec$min_retained_eigenvalue, nystrom_max_eigenvalue = nystrom_feature_spec$max_eigenvalue,
            nystrom_landmark_indices = nystrom_feature_spec$landmark_indices, nystrom_feature_spec = nystrom_feature_spec))
    }
    K <- rkhs_kernel(coordinate_spec$transformed_train, coordinate_spec$transformed_train, selected_sigma,
        kernel_type = kernel_type)
    exact_solution <- solve_gamma_exact_reduced_system(kernel = K, baseline_indicator = W0, right_hand_side = w1,
        ridge = n_t * lambda + numerical_jitter)
    alpha <- exact_solution$solution
    list(predict_raw = function(newdata) {
        K_new <- rkhs_kernel(coordinate_spec$transform(newdata), coordinate_spec$transformed_train, selected_sigma,
            kernel_type = kernel_type)
        as.vector(K_new %*% alpha)
    }, alpha = alpha, sigma = selected_sigma, base_sigma = coordinate_spec$sigma, bandwidth_multiplier = bandwidth_multiplier,
        coordinate_spec = coordinate_spec, backend = backend, kernel_type = kernel_type, nystrom_rank_requested = NA_integer_,
        nystrom_landmark_count = NA_integer_, nystrom_effective_rank = NA_integer_, nystrom_min_retained_eigenvalue = NA_real_,
        nystrom_max_eigenvalue = NA_real_, nystrom_landmark_indices = integer(), exact_solver_method = exact_solution$method,
        exact_solver_minimum_eigenvalue = exact_solution$minimum_eigenvalue, exact_solver_eigenvalue_cutoff = exact_solution$eigenvalue_cutoff)
}

gamma_accuracy_metrics_simple <- function(estimated, oracle) {
    keep <- is.finite(estimated) & is.finite(oracle)
    estimated <- estimated[keep]
    oracle <- oracle[keep]
    if (length(estimated) < 1L) {
        return(data.frame(n = 0L, bias = NA_real_, rmse = NA_real_, mae = NA_real_, correlation = NA_real_,
            mean_estimated = NA_real_, mean_oracle = NA_real_, rmse_constant_one = NA_real_, explained_fraction = NA_real_))
    }
    data.frame(n = length(estimated), bias = mean(estimated - oracle), rmse = sqrt(mean((estimated -
        oracle)^2)), mae = mean(abs(estimated - oracle)), correlation = if (length(estimated) > 1L &&
        is.finite(stats::sd(estimated)) && is.finite(stats::sd(oracle)) && stats::sd(estimated) > 0 &&
        stats::sd(oracle) > 0) {
        stats::cor(estimated, oracle)
    }
    else {
        NA_real_
    }, mean_estimated = mean(estimated), mean_oracle = mean(oracle), rmse_constant_one = sqrt(mean((1 -
        oracle)^2)), explained_fraction = 1 - mean((estimated - oracle)^2)/mean((1 - oracle)^2))
}

make_gamma_nystrom_features <- function(coordinate_spec, landmark_spec, bandwidth_multiplier = 1, kernel_type = c("gaussian",
    "sobolev_matern32", "sobolev_matern52"), eigen_tolerance = 1e-10, numerical_jitter = 1e-10) {
    kernel_type <- match.arg(kernel_type)
    if (length(eigen_tolerance) != 1L || !is.finite(eigen_tolerance) || eigen_tolerance < 0) {
        stop("gamma_nystrom_eigen_tolerance must be nonnegative and finite.")
    }
    if (length(numerical_jitter) != 1L || !is.finite(numerical_jitter) || numerical_jitter < 0) {
        stop("gamma_nystrom_jitter must be nonnegative and finite.")
    }
    sigma <- coordinate_spec$sigma * bandwidth_multiplier
    landmarks <- landmark_spec$transformed_landmarks
    landmark_kernel <- rkhs_kernel(landmarks, landmarks, sigma, kernel_type = kernel_type)
    landmark_kernel <- (landmark_kernel + t(landmark_kernel))/2
    if (numerical_jitter > 0) {
        diag(landmark_kernel) <- diag(landmark_kernel) + numerical_jitter
    }
    decomposition <- eigen(landmark_kernel, symmetric = TRUE)
    maximum_eigenvalue <- max(decomposition$values)
    cutoff <- max(eigen_tolerance * max(maximum_eigenvalue, 0), 0)
    retained <- which(decomposition$values > cutoff)
    if (length(retained) == 0L) {
        retained <- which.max(decomposition$values)
    }
    vectors <- decomposition$vectors[, retained, drop = FALSE]
    values <- decomposition$values[retained]
    if (any(!is.finite(values)) || any(values <= 0)) {
        stop("Nyström landmark kernel has no stable positive eigenspace.")
    }
    feature_map <- function(newdata) {
        kernel_to_landmarks <- rkhs_kernel(coordinate_spec$transform(newdata), landmarks, sigma, kernel_type = kernel_type)
        sweep(kernel_to_landmarks %*% vectors, 2L, sqrt(values), "/")
    }
    list(features_train = sweep(rkhs_kernel(coordinate_spec$transformed_train, landmarks, sigma, kernel_type = kernel_type) %*%
        vectors, 2L, sqrt(values), "/"), predict_features = feature_map, rank_requested = landmark_spec$rank_requested,
        landmark_count = nrow(landmarks), effective_rank = length(values), min_retained_eigenvalue = min(values),
        max_eigenvalue = maximum_eigenvalue, landmark_indices = landmark_spec$landmark_indices, sigma = sigma,
        kernel_type = kernel_type)
}

make_rkhs_cv_plan <- function(train_t, x_cols = "X", cv_folds = 3, cv_repeats = 2L, bandwidth_sample = 1000,
    backend = c("exact", "nystrom"), nystrom_rank = NULL, seed = 1) {
    backend <- match.arg(backend)
    if (length(cv_repeats) != 1L || !is.finite(cv_repeats) || cv_repeats < 1L || cv_repeats != as.integer(cv_repeats)) {
        stop("cv_repeats must be one positive integer.")
    }
    cv_repeats <- as.integer(cv_repeats)
    repeat_fold_id <- lapply(seq_len(cv_repeats), function(repeat_index) {
        make_z_stratified_folds(train_t, K = cv_folds, seed = seed + 100000L * (repeat_index - 1L))
    })
    splits <- vector("list", cv_repeats * cv_folds)
    split_index <- 0L
    for (repeat_index in seq_len(cv_repeats)) {
        for (fold in seq_len(cv_folds)) {
            split_index <- split_index + 1L
            validation_index <- which(repeat_fold_id[[repeat_index]] == fold)
            splits[[split_index]] <- list(repeat_id = repeat_index, fold = fold, training_index = which(repeat_fold_id[[repeat_index]] !=
                fold), validation_index = validation_index)
        }
    }
    coordinate_specs <- vector("list", length(splits))
    nystrom_landmark_specs <- vector("list", length(splits))
    for (split_index in seq_along(splits)) {
        split <- splits[[split_index]]
        repeat_offset <- 100000L * (split$repeat_id - 1L)
        training_index <- split$training_index
        coordinate_specs[[split_index]] <- prepare_rkhs_coordinates(train_t = train_t[training_index,
            , drop = FALSE], x_cols = x_cols, bandwidth_sample = bandwidth_sample, seed = seed + repeat_offset +
            1000L * split$fold)
        if (backend == "nystrom") {
            nystrom_landmark_specs[[split_index]] <- prepare_gamma_nystrom_landmarks(coordinate_spec = coordinate_specs[[split_index]],
                rank = nystrom_rank, seed = seed + repeat_offset + 2000L * split$fold)
        }
    }
    final_coordinate_spec <- prepare_rkhs_coordinates(train_t = train_t, x_cols = x_cols, bandwidth_sample = bandwidth_sample,
        seed = seed + 50000)
    final_nystrom_landmark_spec <- if (backend == "nystrom") {
        prepare_gamma_nystrom_landmarks(coordinate_spec = final_coordinate_spec, rank = nystrom_rank,
            seed = seed + 60000)
    }
    else {
        NULL
    }
    list(fold_id = repeat_fold_id[[1L]], repeat_fold_id = repeat_fold_id, splits = splits, cv_folds = as.integer(cv_folds),
        cv_repeats = cv_repeats, coordinate_specs = coordinate_specs, final_coordinate_spec = final_coordinate_spec,
        backend = backend, nystrom_rank = nystrom_rank, nystrom_landmark_specs = nystrom_landmark_specs,
        final_nystrom_landmark_spec = final_nystrom_landmark_spec)
}

prepare_gamma_nystrom_landmarks <- function(coordinate_spec, rank = NULL, seed = 1) {
    n_training <- nrow(coordinate_spec$transformed_train)
    resolved_rank <- resolve_gamma_nystrom_rank(rank, n_training)
    had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_random_seed) {
        saved_random_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
        if (had_random_seed) {
            assign(".Random.seed", saved_random_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)
    set.seed(seed)
    landmark_indices <- if (resolved_rank == n_training) {
        seq_len(n_training)
    }
    else {
        sort(sample(seq_len(n_training), resolved_rank, replace = FALSE))
    }
    list(rank_requested = resolved_rank, landmark_indices = landmark_indices, transformed_landmarks = coordinate_spec$transformed_train[landmark_indices,
        , drop = FALSE])
}

prepare_rkhs_coordinates <- function(train_t, x_cols = "X", bandwidth_sample = 1000, seed = 1) {
    non_numeric <- x_cols[!vapply(train_t[, x_cols, drop = FALSE], is.numeric, logical(1))]
    if (length(non_numeric) > 0) {
        stop("RKHS coordinates must be numeric: ", paste(non_numeric, collapse = ", "))
    }
    raw_coordinates <- cbind(Y = train_t$Y, as.matrix(train_t[, x_cols, drop = FALSE]))
    center <- colMeans(raw_coordinates)
    scale_value <- apply(raw_coordinates, 2, stats::sd)
    scale_value[!is.finite(scale_value) | scale_value < 1e-08] <- 1
    transform_coordinates <- function(newdata) {
        coordinates <- cbind(Y = newdata$Y, as.matrix(newdata[, x_cols, drop = FALSE]))
        sweep(sweep(coordinates, 2, center, "-"), 2, scale_value, "/")
    }
    transformed_train <- transform_coordinates(train_t)
    set.seed(seed)
    sample_index <- sample(seq_len(nrow(transformed_train)), min(bandwidth_sample, nrow(transformed_train)))
    distances <- as.numeric(stats::dist(transformed_train[sample_index, , drop = FALSE]))
    distances <- distances[is.finite(distances) & distances > 0]
    sigma <- if (length(distances) == 0) {
        1
    }
    else {
        stats::median(distances)
    }
    sigma <- max(sigma, 1e-04)
    list(transformed_train = transformed_train, transform = transform_coordinates, center = center, scale = scale_value,
        sigma = sigma)
}

project_gamma_to_compatible_interval <- function(gamma, nuisance, t_value, Lambda0, Lambda1, margin_fraction = 1e-04) {
    interval <- compute_compatible_gamma_interval(nuisance = nuisance, t_value = t_value, Lambda0 = Lambda0,
        Lambda1 = Lambda1)
    width <- interval$upper - interval$lower
    if (any(!is.finite(width)) || any(width <= 0)) {
        stop(sprintf("The compatible gamma interval is empty for treatment %d.", t_value))
    }
    interior_lower <- interval$lower + margin_fraction * width
    interior_upper <- interval$upper - margin_fraction * width
    projected <- pmin(interior_upper, pmax(interior_lower, gamma))
    list(gamma = projected, original = gamma, lower = interval$lower, upper = interval$upper, projected = abs(projected -
        gamma) > 10 * .Machine$double.eps)
}

resolve_gamma_lambda_grid <- function(lambda_grid, n_training, treatment) {
    if (length(n_training) != 1L || !is.finite(n_training) || n_training < 1L) {
        stop("n_training must be a positive integer.")
    }
    if (is.function(lambda_grid)) {
        lambda_grid <- lambda_grid(n_training = as.integer(n_training), treatment = treatment)
    }
    lambda_grid <- sort(unique(as.numeric(lambda_grid)))
    if (length(lambda_grid) == 0L || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
        stop(paste0("The resolved gamma lambda grid must contain positive finite ", "values."))
    }
    lambda_grid
}

resolve_gamma_nystrom_rank <- function(rank, n_training) {
    if (length(n_training) != 1L || !is.finite(n_training) || n_training < 1L) {
        stop("n_training must be a positive integer.")
    }
    n_training <- as.integer(n_training)
    if (is.function(rank)) {
        rank <- rank(n_training)
    }
    if (is.null(rank)) {
        rank <- ceiling(10 * sqrt(n_training))
    }
    if (length(rank) != 1L || !is.finite(rank) || rank < 1L || abs(rank - round(rank)) > sqrt(.Machine$double.eps)) {
        stop(paste0("gamma_nystrom_rank must be NULL, one positive integer, or a ", "function returning one positive integer from n_training."))
    }
    min(as.integer(round(rank)), n_training)
}

rkhs_kernel <- function(A, B, sigma, kernel_type = c("gaussian", "sobolev_matern32", "sobolev_matern52")) {
    kernel_type <- match.arg(kernel_type)
    if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
        stop("sigma must be one positive finite number.")
    }
    scaled_distance <- rkhs_pairwise_distance(A, B)/sigma
    if (kernel_type == "gaussian") {
        return(exp(-scaled_distance^2/2))
    }
    if (kernel_type == "sobolev_matern32") {
        argument <- sqrt(3) * scaled_distance
        return((1 + argument) * exp(-argument))
    }
    argument <- sqrt(5) * scaled_distance
    (1 + argument + argument^2/3) * exp(-argument)
}

rkhs_pairwise_distance <- function(A, B) {
    A <- as.matrix(A)
    B <- as.matrix(B)
    A2 <- rowSums(A^2)
    B2 <- rowSums(B^2)
    distance_squared <- outer(A2, B2, "+") - 2 * tcrossprod(A, B)
    distance_squared[distance_squared < 0] <- 0
    sqrt(distance_squared)
}

rkhs_validation_risk <- function(validation_t, nuisance_validation_t, gamma_raw, t_value) {
    w_tilde <- compute_riesz_weight(nuisance = nuisance_validation_t, t_value = t_value)
    mean(0.5 * as.numeric(validation_t$Z == 0) * gamma_raw^2 - as.numeric(validation_t$Z == 1) * w_tilde *
        gamma_raw)
}

select_rkhs_lambda <- function(train_t, nuisance_t, t_value, x_cols = "X", lambda_grid = 10^seq(-5, -1,
    length.out = 9), bandwidth_multipliers = 1/sqrt(2), cv_plan, selection_rule = c("minimum", "one_se"),
    bandwidth_sample = 1000, backend = c("exact", "nystrom"), nystrom_rank = NULL, nystrom_eigen_tolerance = 1e-10,
    nystrom_jitter = 1e-10, kernel_type = c("gaussian", "sobolev_matern32", "sobolev_matern52")) {
    backend <- match.arg(backend)
    kernel_type <- match.arg(kernel_type)
    selection_rule <- match.arg(selection_rule)
    lambda_grid <- sort(unique(as.numeric(lambda_grid)))
    if (length(lambda_grid) == 0 || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
        stop("lambda_grid must contain positive finite values.")
    }
    bandwidth_multipliers <- sort(unique(as.numeric(bandwidth_multipliers)))
    if (length(bandwidth_multipliers) == 0 || any(!is.finite(bandwidth_multipliers)) || any(bandwidth_multipliers <=
        0)) {
        stop("bandwidth_multipliers must contain positive finite values.")
    }
    candidate_grid <- expand.grid(lambda = lambda_grid, bandwidth_multiplier = bandwidth_multipliers,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    splits <- if (!is.null(cv_plan$splits)) {
        cv_plan$splits
    }
    else {
        lapply(seq_len(max(cv_plan$fold_id)), function(fold) {
            list(repeat_id = 1L, fold = fold, training_index = which(cv_plan$fold_id != fold), validation_index = which(cv_plan$fold_id ==
                fold))
        })
    }
    cv_resamples <- length(splits)
    fold_risk <- matrix(NA_real_, nrow = nrow(candidate_grid), ncol = cv_resamples)
    fold_effective_rank <- matrix(NA_integer_, nrow = length(bandwidth_multipliers), ncol = cv_resamples)
    for (resample in seq_along(splits)) {
        validation_index <- splits[[resample]]$validation_index
        training_index <- splits[[resample]]$training_index
        cv_train <- train_t[training_index, , drop = FALSE]
        cv_validation <- train_t[validation_index, , drop = FALSE]
        nuisance_cv_train <- subset_nuisance(nuisance_t, training_index)
        nuisance_cv_validation <- subset_nuisance(nuisance_t, validation_index)
        nystrom_feature_specs <- if (backend == "nystrom") {
            lapply(bandwidth_multipliers, function(multiplier) {
                make_gamma_nystrom_features(coordinate_spec = cv_plan$coordinate_specs[[resample]], landmark_spec = cv_plan$nystrom_landmark_specs[[resample]],
                  bandwidth_multiplier = multiplier, kernel_type = kernel_type, eigen_tolerance = nystrom_eigen_tolerance,
                  numerical_jitter = nystrom_jitter)
            })
        }
        else {
            NULL
        }
        if (backend == "nystrom") {
            fold_effective_rank[, resample] <- vapply(nystrom_feature_specs, function(specification) specification$effective_rank,
                integer(1))
        }
        for (candidate in seq_len(nrow(candidate_grid))) {
            bandwidth_index <- match(candidate_grid$bandwidth_multiplier[candidate], bandwidth_multipliers)
            fit <- fit_gamma_rkhs_reduced(train_t = cv_train, nuisance_t = nuisance_cv_train, t_value = t_value,
                x_cols = x_cols, lambda = candidate_grid$lambda[candidate], bandwidth_multiplier = candidate_grid$bandwidth_multiplier[candidate],
                coordinate_spec = cv_plan$coordinate_specs[[resample]], bandwidth_sample = bandwidth_sample,
                backend = backend, nystrom_rank = nystrom_rank, nystrom_landmark_spec = cv_plan$nystrom_landmark_specs[[resample]],
                nystrom_feature_spec = if (backend == "nystrom") {
                  nystrom_feature_specs[[bandwidth_index]]
                }
                else {
                  NULL
                }, nystrom_eigen_tolerance = nystrom_eigen_tolerance, nystrom_jitter = nystrom_jitter,
                kernel_type = kernel_type)
            raw_validation <- fit$predict_raw(cv_validation)
            fold_risk[candidate, resample] <- rkhs_validation_risk(validation_t = cv_validation, nuisance_validation_t = nuisance_cv_validation,
                gamma_raw = raw_validation, t_value = t_value)
        }
    }
    mean_risk <- rowMeans(fold_risk)
    standard_error <- apply(fold_risk, 1, stats::sd)/sqrt(cv_resamples)
    minimum_index <- which.min(mean_risk)
    selected_index <- if (selection_rule == "minimum") {
        minimum_index
    }
    else {
        threshold <- mean_risk[minimum_index] + standard_error[minimum_index]
        acceptable <- which(mean_risk <= threshold)
        ordering <- order(candidate_grid$lambda[acceptable], candidate_grid$bandwidth_multiplier[acceptable],
            decreasing = TRUE)
        acceptable[ordering[1]]
    }
    list(selected_lambda = candidate_grid$lambda[selected_index], selected_bandwidth_multiplier = candidate_grid$bandwidth_multiplier[selected_index],
        cv_table = data.frame(lambda = candidate_grid$lambda, bandwidth_multiplier = candidate_grid$bandwidth_multiplier,
            rkhs_backend = backend, rkhs_kernel = kernel_type, nystrom_effective_rank_min = if (backend ==
                "nystrom") {
                apply(fold_effective_rank, 1L, min)[match(candidate_grid$bandwidth_multiplier, bandwidth_multipliers)]
            } else {
                NA_integer_
            }, nystrom_effective_rank_max = if (backend == "nystrom") {
                apply(fold_effective_rank, 1L, max)[match(candidate_grid$bandwidth_multiplier, bandwidth_multipliers)]
            } else {
                NA_integer_
            }, cv_folds = if (is.null(cv_plan$cv_folds)) {
                max(cv_plan$fold_id)
            } else {
                cv_plan$cv_folds
            }, cv_repeats = if (is.null(cv_plan$cv_repeats)) {
                1L
            } else {
                cv_plan$cv_repeats
            }, cv_resamples = cv_resamples, mean_risk = mean_risk, standard_error = standard_error, selected = seq_len(nrow(candidate_grid)) ==
                selected_index), cv_resample_table = do.call(rbind, lapply(seq_len(nrow(candidate_grid)),
            function(candidate) {
                data.frame(lambda = candidate_grid$lambda[candidate], bandwidth_multiplier = candidate_grid$bandwidth_multiplier[candidate],
                  cv_repeat = vapply(splits, `[[`, integer(1L), "repeat_id"), fold = vapply(splits, `[[`,
                    integer(1L), "fold"), risk = as.numeric(fold_risk[candidate, ]), selected = candidate ==
                    selected_index, stringsAsFactors = FALSE)
            })), fold_risk = fold_risk)
}

solve_gamma_exact_reduced_system <- function(kernel, baseline_indicator, right_hand_side, ridge, eigen_tolerance = 1e-12) {
    kernel <- as.matrix(kernel)
    n_training <- nrow(kernel)
    if (n_training < 1L || ncol(kernel) != n_training || length(baseline_indicator) != n_training ||
        length(right_hand_side) != n_training) {
        stop("The exact gamma system has incompatible dimensions.")
    }
    if (length(ridge) != 1L || !is.finite(ridge) || ridge <= 0) {
        stop("ridge must be one positive finite number.")
    }
    if (length(eigen_tolerance) != 1L || !is.finite(eigen_tolerance) || eigen_tolerance <= 0) {
        stop("eigen_tolerance must be one positive finite number.")
    }
    baseline <- which(as.numeric(baseline_indicator) > 0.5)
    comparison <- setdiff(seq_len(n_training), baseline)
    alpha <- numeric(n_training)
    if (length(comparison) > 0L) {
        alpha[comparison] <- right_hand_side[comparison]/ridge
    }
    if (length(baseline) == 0L) {
        return(list(solution = alpha, method = "comparison_closed_form", minimum_eigenvalue = NA_real_,
            eigenvalue_cutoff = NA_real_))
    }
    baseline_kernel <- kernel[baseline, baseline, drop = FALSE]
    baseline_kernel <- (baseline_kernel + t(baseline_kernel))/2
    diag(baseline_kernel) <- diag(baseline_kernel) + ridge
    baseline_rhs <- right_hand_side[baseline]
    if (length(comparison) > 0L) {
        baseline_rhs <- baseline_rhs - as.vector(kernel[baseline, comparison, drop = FALSE] %*% alpha[comparison])
    }
    factor <- tryCatch(chol(baseline_kernel), error = function(error) NULL)
    if (!is.null(factor)) {
        alpha[baseline] <- backsolve(factor, forwardsolve(t(factor), baseline_rhs))
        return(list(solution = alpha, method = "block_cholesky", minimum_eigenvalue = NA_real_, eigenvalue_cutoff = NA_real_))
    }
    decomposition <- eigen(baseline_kernel, symmetric = TRUE)
    maximum_eigenvalue <- max(decomposition$values)
    cutoff <- eigen_tolerance * max(1, maximum_eigenvalue)
    stabilized_values <- pmax(decomposition$values, cutoff)
    alpha[baseline] <- as.vector(decomposition$vectors %*% (crossprod(decomposition$vectors, baseline_rhs)/stabilized_values))
    list(solution = alpha, method = "block_eigenvalue_floor", minimum_eigenvalue = min(decomposition$values),
        eigenvalue_cutoff = cutoff)
}

x5_adaptive_provider <- function(method) {
    stopifnot(method == "rkhs_tuned")
    cache <- new.env(parent = emptyenv())
    provider <- function(fold, treatment, train, nuisance_train, x_cols, seed, ...) {
        key <- paste(fold, treatment, sep = "_")
        signature <- list(data = train[, c("Y", x_cols, "T_obs", "Z")], nu = nuisance_train, seed = seed)
        if (exists(key, cache, inherits = FALSE)) {
            z <- get(key, cache)
            stopifnot(identical(z$signature, signature))
            return(z$fit)
        }
        fit <- fit_gamma_rkhs(train, nuisance_train, treatment, x_cols, lambda_grid = default_gamma_lambda_grid,
            bandwidth_multipliers = c(0.25, 0.5, 1, 2, 4), cv_folds = 3, cv_repeats = 1, selection_rule = "minimum",
            backend = "nystrom", nystrom_rank = NULL, seed = seed)
        assign(key, list(signature = signature, fit = fit, fold = fold, treatment = treatment), cache)
        fit
    }
    attr(provider, "gamma_cache") <- cache
    provider
}

