# Main-text reproduction dependencies; published formulas retained.
binary_available_library_candidates <- function(candidates) {
    candidates <- unique(as.character(candidates))
    allowed <- c("glm_spline", "gam", "glmnet", "ranger")
    if (length(candidates) < 1L || any(!candidates %in% allowed)) {
        stop("Binary nuisance-library candidates must be drawn from: ", paste(allowed, collapse = ", "),
            ".")
    }
    package_by_candidate <- c(gam = "mgcv", glmnet = "glmnet", ranger = "ranger")
    available <- vapply(candidates, function(candidate) {
        if (candidate == "glm_spline")
            return(TRUE)
        requireNamespace(package_by_candidate[[candidate]], quietly = TRUE)
    }, logical(1))
    candidates[available]
}

binary_check_indicator <- function(value, name) {
    observed <- sort(unique(stats::na.omit(value)))
    if (length(observed) == 0L || !all(observed %in% c(0, 1))) {
        stop(name, " must be coded 0/1.")
    }
}

binary_cross_entropy <- function(outcome, probability) {
    outcome <- as.numeric(outcome)
    probability <- pmin(1 - 1e-08, pmax(1e-08, as.numeric(probability)))
    if (length(outcome) != length(probability) || any(!is.finite(outcome)) || any(!is.finite(probability))) {
        stop("Invalid inputs to binary_cross_entropy().")
    }
    mean(-outcome * log(probability) - (1 - outcome) * log1p(-probability))
}

binary_cross_entropy_loss <- function(outcome, probability) {
    outcome <- as.numeric(outcome)
    probability <- pmin(1 - 1e-08, pmax(1e-08, as.numeric(probability)))
    if (length(outcome) != length(probability) || any(!is.finite(outcome)) || any(!is.finite(probability))) {
        stop("Invalid inputs to binary_cross_entropy_loss().")
    }
    mean(-outcome * log(probability) - (1 - outcome) * log1p(-probability))
}

binary_fit_regression_library_cv <- function(data, response, x_cols, spline_df, gam_basis_dimension = NULL,
    gam_min_basis = 5L, gam_max_basis = 20L, gam_growth_exponent = 1/3, seed = 1L, candidates = c("glm_spline",
        "gam", "glmnet", "ranger"), cv_folds = 3L, selection_rule = c("one_se", "minimum"), ranger_num_trees = 500L) {
    selection_rule <- match.arg(selection_rule)
    candidates_requested <- unique(as.character(candidates))
    candidates <- binary_available_library_candidates(candidates_requested)
    if (length(candidates) < 1L) {
        stop("No requested binary nuisance-library candidate is available.")
    }
    outcome <- as.numeric(data[[response]])
    binary_check_indicator(outcome, response)
    fit_candidate <- function(training, candidate, local_seed) {
        binary_fit_regression_with_learner(data = training, response = response, x_cols = x_cols, spline_df = spline_df,
            learner = candidate, gam_basis_dimension = gam_basis_dimension, gam_min_basis = gam_min_basis,
            gam_max_basis = gam_max_basis, gam_growth_exponent = gam_growth_exponent, ranger_num_trees = ranger_num_trees,
            seed = local_seed)
    }
    fold_id <- binary_make_response_folds(outcome, K = cv_folds, seed = seed)
    if (is.null(fold_id)) {
        selected <- if ("glm_spline" %in% candidates)
            "glm_spline"
        else candidates[1L]
        final_fit <- fit_candidate(data, selected, seed + 900000L)
        final_fit$learner <- "library_cv"
        final_fit$selected_learner <- selected
        final_fit$cv_mean_risk <- NA_real_
        final_fit$cv_standard_error <- NA_real_
        final_fit$cv_table <- data.frame(candidate = selected, mean_risk = NA_real_, standard_error = NA_real_,
            selected = TRUE, cv_available = FALSE, stringsAsFactors = FALSE)
        final_fit$library_configuration <- list(candidates_requested = candidates_requested, candidates_available = candidates,
            cv_folds = 0L, selection_rule = selection_rule, selected_learner = selected)
        return(final_fit)
    }
    folds <- sort(unique(fold_id))
    risk_rows <- list()
    counter <- 0L
    for (candidate_index in seq_along(candidates)) {
        candidate <- candidates[candidate_index]
        for (fold in folds) {
            validation <- fold_id == fold
            fitted <- tryCatch(fit_candidate(data[!validation, , drop = FALSE], candidate, seed + 10000L *
                candidate_index + fold), error = identity)
            counter <- counter + 1L
            if (inherits(fitted, "error")) {
                risk_rows[[counter]] <- data.frame(candidate = candidate, fold = fold, risk = Inf, error = conditionMessage(fitted),
                  stringsAsFactors = FALSE)
            }
            else {
                prediction <- fitted$predict(data[validation, , drop = FALSE])
                risk_rows[[counter]] <- data.frame(candidate = candidate, fold = fold, risk = binary_cross_entropy_loss(outcome[validation],
                  prediction), error = "", stringsAsFactors = FALSE)
            }
        }
    }
    fold_risk <- do.call(rbind, risk_rows)
    cv_table <- do.call(rbind, lapply(candidates, function(candidate) {
        local <- fold_risk[fold_risk$candidate == candidate, , drop = FALSE]
        finite <- is.finite(local$risk)
        data.frame(candidate = candidate, mean_risk = if (all(finite))
            mean(local$risk)
        else Inf, standard_error = if (all(finite) && nrow(local) > 1L) {
            stats::sd(local$risk)/sqrt(nrow(local))
        }
        else if (all(finite)) {
            0
        }
        else {
            NA_real_
        }, failed_folds = sum(!finite), error = paste(unique(local$error[nzchar(local$error)]), collapse = " | "),
            stringsAsFactors = FALSE)
    }))
    if (!any(is.finite(cv_table$mean_risk))) {
        stop("Every binary nuisance-library candidate failed cross-validation.")
    }
    minimum <- which.min(cv_table$mean_risk)
    choice <- minimum
    if (selection_rule == "one_se") {
        cutoff <- cv_table$mean_risk[minimum] + ifelse(is.finite(cv_table$standard_error[minimum]), cv_table$standard_error[minimum],
            0)
        eligible <- which(cv_table$mean_risk <= cutoff)
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
    final_fit$library_configuration <- list(candidates_requested = candidates_requested, candidates_available = candidates,
        cv_folds = length(folds), selection_rule = selection_rule, selected_learner = selected)
    final_fit
}

binary_fit_regression_with_learner <- function(data, response, x_cols, spline_df, learner = c("glm_spline",
    "gam", "ranger", "glmnet", "library_cv", "super_learner"), gam_basis_dimension = NULL, gam_min_basis = 5L,
    gam_max_basis = 20L, gam_growth_exponent = 1/3, library_candidates = c("glm_spline", "gam", "glmnet",
        "ranger"), library_cv_folds = 3L, library_selection_rule = c("one_se", "minimum"), ranger_num_trees = 500L,
    super_learner_candidates = binary_default_super_learner_candidates(), super_learner_cv_folds = 3L,
    super_learner_ranger_num_trees = 500L, super_learner_hal_max_degree = 2L, seed = 1L) {
    learner <- match.arg(learner)
    library_selection_rule <- match.arg(library_selection_rule)
    outcome <- as.numeric(data[[response]])
    binary_check_indicator(outcome, response)
    if (length(unique(outcome)) == 1L) {
        probability <- (sum(outcome) + 0.5)/(length(outcome) + 1)
        return(list(predict = function(newdata) rep(probability, nrow(newdata)), converged = TRUE, rank = 1L,
            constant = TRUE, learner = learner, selected_learner = if (learner == "library_cv") "constant" else learner,
            cv_mean_risk = NA_real_, cv_standard_error = NA_real_, cv_table = data.frame()))
    }
    if (learner == "library_cv") {
        return(binary_fit_regression_library_cv(data = data, response = response, x_cols = x_cols, spline_df = spline_df,
            gam_basis_dimension = gam_basis_dimension, gam_min_basis = gam_min_basis, gam_max_basis = gam_max_basis,
            gam_growth_exponent = gam_growth_exponent, seed = seed, candidates = library_candidates,
            cv_folds = library_cv_folds, selection_rule = library_selection_rule, ranger_num_trees = ranger_num_trees))
    }
    if (learner == "super_learner") {
        return(binary_fit_regression_super_learner(data = data, response = response, x_cols = x_cols,
            candidates = super_learner_candidates, cv_folds = super_learner_cv_folds, ranger_num_trees = super_learner_ranger_num_trees,
            hal_max_degree = super_learner_hal_max_degree, seed = seed))
    }
    if (learner == "gam") {
        if (!requireNamespace("mgcv", quietly = TRUE)) {
            stop("The binary GAM nuisance learner requires package 'mgcv'.")
        }
        selected_basis_dimension <- binary_resolve_gam_basis_dimension(data = data, x_cols = x_cols,
            basis_dimension = gam_basis_dimension, minimum = gam_min_basis, maximum = gam_max_basis,
            growth_exponent = gam_growth_exponent)
        fit <- suppressWarnings(mgcv::gam(binary_make_gam_formula(response, data, x_cols, selected_basis_dimension),
            data = data, family = stats::binomial(), method = "REML", select = TRUE))
        return(list(predict = function(newdata) {
            as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
        }, converged = isTRUE(fit$converged), rank = fit$rank, constant = FALSE, learner = learner, basis_dimension = selected_basis_dimension,
            effective_df = sum(fit$edf)))
    }
    if (learner == "ranger") {
        if (!requireNamespace("ranger", quietly = TRUE)) {
            stop("The binary ranger nuisance learner requires package 'ranger'.")
        }
        ranger_data <- data[, c(response, x_cols), drop = FALSE]
        ranger_data[[response]] <- factor(outcome, levels = c(0, 1))
        fit <- ranger::ranger(formula = stats::reformulate(x_cols, response = response), data = ranger_data,
            probability = TRUE, num.trees = as.integer(ranger_num_trees), mtry = max(1L, floor(sqrt(length(x_cols)))),
            min.node.size = 20L, seed = as.integer(seed), respect.unordered.factors = "order")
        return(list(predict = function(newdata) {
            prediction <- predict(fit, data = newdata)$predictions
            if (is.null(dim(prediction))) return(as.numeric(prediction))
            if ("1" %in% colnames(prediction)) {
                return(as.numeric(prediction[, "1"]))
            }
            as.numeric(prediction[, ncol(prediction)])
        }, converged = TRUE, rank = as.integer(length(x_cols)), constant = FALSE, learner = learner))
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
        single_predictor_continuous <- ncol(design) == 1L && length(unique(design[, 1L])) > 2L
        if (ncol(design) == 1L) {
            auxiliary <- if (single_predictor_continuous) {
                design[, 1L]^2
            }
            else {
                rep(0, nrow(design))
            }
            design <- cbind(design, .glmnet_auxiliary = auxiliary)
        }
        minimum_class_size <- min(tabulate(outcome + 1L, nbins = 2L))
        cv_folds <- min(5L, minimum_class_size)
        if (cv_folds < 2L) {
            stop("The glmnet nuisance fit requires both outcome classes.")
        }
        had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (had_seed)
            old_seed <- get(".Random.seed", envir = .GlobalEnv)
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
        fit <- suppressWarnings(glmnet::cv.glmnet(x = design, y = outcome, family = "binomial", alpha = 1,
            foldid = fold_id, type.measure = "deviance", standardize = TRUE, intercept = TRUE, parallel = FALSE))
        coefficients <- stats::coef(fit, s = "lambda.1se")
        active_count <- sum(as.numeric(coefficients) != 0) - 1L
        return(list(predict = function(newdata) {
            new_design <- as.matrix(newdata[, x_cols, drop = FALSE])
            storage.mode(new_design) <- "double"
            if (ncol(new_design) == 1L) {
                auxiliary <- if (single_predictor_continuous) {
                  new_design[, 1L]^2
                } else {
                  rep(0, nrow(new_design))
                }
                new_design <- cbind(new_design, .glmnet_auxiliary = auxiliary)
            }
            as.numeric(stats::predict(fit, newx = new_design, s = "lambda.1se", type = "response"))
        }, converged = TRUE, rank = as.integer(max(0L, active_count)), constant = FALSE, learner = learner,
            selected_lambda = unname(fit$lambda.1se), minimum_lambda = unname(fit$lambda.min)))
    }
    fit <- suppressWarnings(stats::glm(binary_make_spline_formula(response, data, x_cols, spline_df),
        data = data, family = stats::binomial(), control = stats::glm.control(maxit = 100L)))
    list(predict = function(newdata) {
        as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
    }, converged = isTRUE(fit$converged), rank = fit$rank, constant = FALSE, learner = learner)
}

binary_make_gam_formula <- function(response, data, x_cols, basis_dimension) {
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

binary_make_interaction_gam_formula <- function(response, data, x_cols, basis_dimension, interaction_basis_dimension = 3L) {
    unique_counts <- vapply(x_cols, function(variable) {
        length(unique(stats::na.omit(data[[variable]])))
    }, integer(1))
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
        }, character(1))
    }
    stats::as.formula(paste0("`", response, "` ~ ", paste(c(main_terms, interaction_terms), collapse = " + ")))
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
    if (K < 2L)
        return(NULL)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed)
        old_seed <- get(".Random.seed", envir = .GlobalEnv)
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

binary_make_spline_formula <- function(response, data, x_cols, spline_df) {
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

binary_pairwise_interaction_design <- function(data) {
    data <- as.data.frame(data)
    if (!ncol(data))
        stop("At least one predictor is required.")
    design <- as.matrix(data)
    storage.mode(design) <- "double"
    if (any(!is.finite(design))) {
        stop("The pairwise-interaction design contains non-finite values.")
    }
    if (ncol(design) == 1L) {
        squared <- matrix(design[, 1L]^2, ncol = 1L, dimnames = list(NULL, paste0(colnames(design), "_squared")))
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
    }, character(1))
    cbind(design, interactions)
}

binary_resolve_gam_basis_dimension <- function(data, x_cols, basis_dimension = NULL, minimum = 5L, maximum = 20L,
    growth_exponent = 1/3) {
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
        0 || growth_exponent > 1) {
        stop("Invalid adaptive GAM basis controls.")
    }
    observations_per_covariate <- nrow(data)/max(1L, length(x_cols))
    proposed <- ceiling(observations_per_covariate^growth_exponent)
    as.integer(min(maximum, max(minimum, proposed)))
}

calibrate_learned_witness_maximum <- function(signals, metadata, bootstrap_replications = 999L, confidence_level = 0.95,
    seed = 1L, cluster_id = NULL) {
    n <- nrow(signals)
    clustered <- !is.null(cluster_id)
    if (is.null(cluster_id))
        cluster_id <- seq_len(n)
    if (length(cluster_id) != n || anyNA(cluster_id)) {
        stop("cluster_id must be complete and have one value per signal row.")
    }
    cluster_id <- as.character(cluster_id)
    G <- length(unique(cluster_id))
    if (G < 2L)
        stop("At least two independent clusters are required.")
    estimate <- colMeans(signals)
    centered_individual <- sweep(signals, 2L, estimate, "-")
    cluster_signal <- rowsum(centered_individual, cluster_id, reorder = FALSE) * (G/n)
    standard_deviation <- apply(cluster_signal, 2L, stats::sd)
    keep <- is.finite(estimate) & is.finite(standard_deviation) & standard_deviation > 1e-08
    if (!any(keep))
        stop("No nondegenerate learned-witness moments remain.")
    signals <- signals[, keep, drop = FALSE]
    metadata <- metadata[keep, , drop = FALSE]
    estimate <- estimate[keep]
    standard_deviation <- standard_deviation[keep]
    cluster_signal <- cluster_signal[, keep, drop = FALSE]
    z_value <- sqrt(G) * estimate/standard_deviation
    lambda_key <- interaction(metadata$Lambda0, metadata$Lambda1, drop = TRUE)
    lambda_levels <- levels(lambda_key)
    observed <- vapply(lambda_levels, function(level) max(z_value[lambda_key == level]), numeric(1))
    observed_cvm <- vapply(lambda_levels, function(level) mean(pmax(z_value[lambda_key == level], 0)^2),
        numeric(1))
    centered <- sweep(cluster_signal, 2L, colMeans(cluster_signal), "-")
    centered <- sweep(centered, 2L, standard_deviation, "/")
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed)
        old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
        if (had_seed) {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)
    set.seed(seed)
    bootstrap <- matrix(NA_real_, bootstrap_replications, length(lambda_levels), dimnames = list(NULL,
        lambda_levels))
    bootstrap_cvm <- bootstrap
    for (replication in seq_len(bootstrap_replications)) {
        draw <- as.vector(crossprod(stats::rnorm(G), centered))/sqrt(G)
        bootstrap[replication, ] <- vapply(lambda_levels, function(level) max(draw[lambda_key == level]),
            numeric(1))
        bootstrap_cvm[replication, ] <- vapply(lambda_levels, function(level) mean(pmax(draw[lambda_key ==
            level], 0)^2), numeric(1))
    }
    alpha <- 1 - confidence_level
    rows <- lapply(seq_along(lambda_levels), function(index) {
        local_metadata <- metadata[lambda_key == lambda_levels[index], ]
        local_bootstrap <- bootstrap[, index]
        local_cvm_bootstrap <- bootstrap_cvm[, index]
        data.frame(Lambda0 = local_metadata$Lambda0[1L], Lambda1 = local_metadata$Lambda1[1L], statistic = observed[index],
            critical_value = as.numeric(stats::quantile(local_bootstrap, confidence_level, type = 8,
                names = FALSE)), p_value = (1 + sum(local_bootstrap >= observed[index]))/(bootstrap_replications +
                1), reject = FALSE, cvm_statistic = observed_cvm[index], cvm_critical_value = as.numeric(stats::quantile(local_cvm_bootstrap,
                confidence_level, type = 8, names = FALSE)), cvm_p_value = (1 + sum(local_cvm_bootstrap >=
                observed_cvm[index]))/(bootstrap_replications + 1), cvm_reject = FALSE, moments_tested = nrow(local_metadata),
            stringsAsFactors = FALSE)
    })
    summary <- do.call(rbind, rows)
    summary$reject <- summary$p_value <= alpha
    summary$cvm_reject <- summary$cvm_p_value <= alpha
    summary <- summary[order(summary$Lambda0, summary$Lambda1), ]
    summary$ordered_reject <- FALSE
    summary$cvm_ordered_reject <- FALSE
    continue <- TRUE
    for (index in seq_len(nrow(summary))) {
        if (continue && summary$reject[index]) {
            summary$ordered_reject[index] <- TRUE
        }
        else {
            continue <- FALSE
        }
    }
    continue <- TRUE
    for (index in seq_len(nrow(summary))) {
        if (continue && summary$cvm_reject[index]) {
            summary$cvm_ordered_reject[index] <- TRUE
        }
        else {
            continue <- FALSE
        }
    }
    simultaneous_bootstrap <- apply(bootstrap, 1L, max)
    simultaneous_cvm_bootstrap <- apply(bootstrap_cvm, 1L, max)
    summary$grid_simultaneous_p_value <- vapply(summary$statistic, function(value) (1 + sum(simultaneous_bootstrap >=
        value))/(bootstrap_replications + 1), numeric(1))
    summary$grid_simultaneous_reject <- summary$grid_simultaneous_p_value <= alpha
    summary$cvm_grid_simultaneous_p_value <- vapply(summary$cvm_statistic, function(value) {
        (1 + sum(simultaneous_cvm_bootstrap >= value))/(bootstrap_replications + 1)
    }, numeric(1))
    summary$cvm_grid_simultaneous_reject <- summary$cvm_grid_simultaneous_p_value <= alpha
    moment_detail <- cbind(metadata, estimate = estimate, standard_error = standard_deviation/sqrt(G),
        z_value = z_value)
    moment_detail <- moment_detail[order(moment_detail$Lambda0, -moment_detail$z_value), ]
    list(summary = summary, moment_detail = moment_detail, bootstrap = bootstrap, cluster_count = G,
        clustered_inference = clustered)
}

check_binary <- function(x, name) {
    values <- sort(unique(stats::na.omit(x)))
    if (!all(values %in% c(0, 1))) {
        stop(sprintf("%s must be coded 0/1.", name))
    }
}

check_learned_witness_dependencies <- function() {
    required <- c("fit_nuisance_models", "fit_gamma_rkhs_fixed", "make_tz_stratified_folds", "fit_series_conditional_allocations",
        "compute_dual_msm_moment_primitives", "make_arm_lambda_signal_blocks")
    missing <- required[!vapply(required, exists, logical(1), mode = "function")]
    if (length(missing) > 0L) {
        stop("Source R/sharp_bounds_plugin.R and ", "R/conservative_compatibility_test.R first. Missing: ",
            paste(missing, collapse = ", "), call. = FALSE)
    }
    invisible(TRUE)
}

clip_positive <- function(x, floor = 1e-06, cap = Inf) {
    x[!is.finite(x)] <- floor
    x <- pmax(floor, x)
    if (is.finite(cap)) {
        x <- pmin(cap, x)
    }
    x
}

clip_probability <- function(p, epsilon = 0.05) {
    pmin(1 - epsilon, pmax(epsilon, p))
}

compute_dual_msm_moment_primitives <- function(p0, p1, Lambda0, Lambda1) {
    n_lower <- p0 + (1 - p0)/Lambda0
    d_lower <- p1 + (1 - p1) * Lambda1
    gamma_lower <- n_lower/d_lower
    d_gamma_lower_p0 <- (1 - 1/Lambda0)/d_lower
    d_gamma_lower_p1 <- n_lower * (Lambda1 - 1)/d_lower^2
    n_upper <- p0 + Lambda0 * (1 - p0)
    d_upper <- p1 + (1 - p1)/Lambda1
    gamma_upper <- n_upper/d_upper
    d_gamma_upper_p0 <- (1 - Lambda0)/d_upper
    d_gamma_upper_p1 <- -n_upper * (1 - 1/Lambda1)/d_upper^2
    denominator <- 1 - p0
    c_lower <- d_upper/denominator
    c_upper <- d_lower/denominator
    b <- -p0/denominator
    list(gamma_lower = gamma_lower, gamma_upper = gamma_upper, d_gamma_lower_p0 = d_gamma_lower_p0, d_gamma_lower_p1 = d_gamma_lower_p1,
        d_gamma_upper_p0 = d_gamma_upper_p0, d_gamma_upper_p1 = d_gamma_upper_p1, c_lower = c_lower,
        c_upper = c_upper, b = b, d_c_lower_p0 = d_upper/denominator^2, d_c_lower_p1 = (1 - 1/Lambda1)/denominator,
        d_c_upper_p0 = d_lower/denominator^2, d_c_upper_p1 = (1 - Lambda1)/denominator, d_b_p0 = -1/denominator^2)
}

compute_learned_family_base <- function(family, primitive, m0, m1, Lambda) {
    switch(family, pointwise_lower = primitive$gamma_lower * m0 - m1, pointwise_upper = m1 - primitive$gamma_upper *
        m0, normalization_lower = primitive$c_lower * m1 + (primitive$b - 1/Lambda) * m0 - (1 - 1/Lambda),
        normalization_upper = -primitive$c_upper * m1 + (Lambda - primitive$b) * m0 - (Lambda - 1), stop("Unknown learned-witness family: ",
            family))
}

compute_learned_witness_scores <- function(gamma, pi0, pi1, treatment, Lambda) {
    p0 <- if (treatment == 1L)
        pi0
    else 1 - pi0
    p1 <- if (treatment == 1L)
        pi1
    else 1 - pi1
    primitive <- compute_dual_msm_moment_primitives(p0, p1, Lambda, Lambda)
    data.frame(pointwise_lower = primitive$gamma_lower - gamma, pointwise_upper = gamma - primitive$gamma_upper,
        normalization_lower = primitive$c_lower * gamma + primitive$b - 1/Lambda, normalization_upper = Lambda -
            primitive$c_upper * gamma - primitive$b, stringsAsFactors = FALSE)
}

d_u <- function(u) {
    dtnormtc(u, mu = 0, s = 1)
}

dtnormtc <- function(x, mu = 0, s = 1) {
    standardized <- (x - mu)/s
    normalizing_constant <- stats::pnorm(2) - stats::pnorm(-2)
    ifelse(standardized >= -2 & standardized <= 2, stats::dnorm(standardized)/(s * normalizing_constant),
        0)
}

fit_fold_learned_witness_signals <- function(training, validation, probability_training, probability_validation,
    gamma_training, gamma_validation, treatment, Lambda, x_cols, allocation_spline_df, witness_tail_fractions,
    localization_tail_fractions, witness_shapes) {
    score_training <- compute_learned_witness_scores(gamma_training, probability_training$pi0, probability_training$pi1,
        treatment, Lambda)
    score_validation <- compute_learned_witness_scores(gamma_validation, probability_validation$pi0,
        probability_validation$pi1, treatment, Lambda)
    families <- names(score_training)
    response_training <- response_validation <- list()
    response_family <- character()
    for (family in families) {
        family_response_count <- length(response_training)
        score_local <- score_training[[family]]
        active_local <- training$T_obs == treatment & is.finite(score_local)
        cutoffs <- make_witness_cutoffs(score_local, active = active_local, upper_tail_fractions = witness_tail_fractions)
        if ("indicator" %in% witness_shapes) {
            for (name in names(cutoffs)) {
                key <- paste(family, "indicator", name, sep = "__")
                response_training[[key]] <- as.numeric(score_training[[family]] > cutoffs[name])
                response_validation[[key]] <- as.numeric(score_validation[[family]] > cutoffs[name])
                response_family[key] <- family
            }
        }
        positive_training <- pmax(score_training[[family]], 0)
        positive_validation <- pmax(score_validation[[family]], 0)
        scale <- max(positive_training[active_local], na.rm = TRUE)
        if (is.finite(scale) && scale > 1e-10) {
            for (shape in intersect(witness_shapes, c("linear", "sqrt"))) {
                transform <- if (shape == "linear")
                  identity
                else sqrt
                key <- paste(family, shape, "positive", sep = "__")
                response_training[[key]] <- pmin(transform(positive_training/scale), 1)
                response_validation[[key]] <- pmin(transform(positive_validation/scale), 1)
                response_family[key] <- family
            }
        }
        if (length(response_training) == family_response_count) {
            key <- paste(family, "fallback", "zero", sep = "__")
            response_training[[key]] <- as.numeric(score_training[[family]] > 0)
            response_validation[[key]] <- as.numeric(score_validation[[family]] > 0)
            response_family[key] <- family
        }
    }
    response_training <- do.call(cbind, response_training)
    response_validation <- do.call(cbind, response_validation)
    storage.mode(response_training) <- "double"
    storage.mode(response_validation) <- "double"
    mean_fit <- fit_witness_allocation_means(training, validation, response_training, treatment, x_cols,
        allocation_spline_df)
    p0_training <- if (treatment == 1L) {
        probability_training$pi0
    }
    else 1 - probability_training$pi0
    p1_training <- if (treatment == 1L) {
        probability_training$pi1
    }
    else 1 - probability_training$pi1
    p0_validation <- if (treatment == 1L) {
        probability_validation$pi0
    }
    else 1 - probability_validation$pi0
    p1_validation <- if (treatment == 1L) {
        probability_validation$pi1
    }
    else 1 - probability_validation$pi1
    primitive_training <- compute_dual_msm_moment_primitives(p0_training, p1_training, Lambda, Lambda)
    primitive_validation <- compute_dual_msm_moment_primitives(p0_validation, p1_validation, Lambda,
        Lambda)
    output_signals <- output_metadata <- diagnostics <- list()
    for (family in families) {
        columns <- which(response_family == family)
        local_mean_fit <- lapply(mean_fit, function(value) {
            value[, columns, drop = FALSE]
        })
        selected <- select_learned_witness(family = family, response_training = response_training[, columns,
            drop = FALSE], response_validation = response_validation[, columns, drop = FALSE], mean_fit = local_mean_fit,
            primitive_training = primitive_training, primitive_validation = primitive_validation, Lambda = Lambda,
            localization_tail_fractions = localization_tail_fractions)
        block <- make_arm_lambda_signal_blocks(data = validation, H = selected$H_validation, R = selected$R_validation,
            m0 = selected$m0_validation, m1 = selected$m1_validation, e = probability_validation$e, pi0 = probability_validation$pi0,
            pi1 = probability_validation$pi1, t_value = treatment, Lambda0 = Lambda, Lambda1 = Lambda)
        keep <- block$metadata$family == family
        if (sum(keep) != 1L)
            stop("A learned witness did not map to one moment.")
        output_signals[[family]] <- block$signals[, keep, drop = FALSE]
        local_metadata <- block$metadata[keep, , drop = FALSE]
        local_metadata$r <- selected$r_name
        local_metadata$h <- paste0("learned_", selected$h_name)
        output_metadata[[family]] <- local_metadata
        diagnostics[[family]] <- data.frame(treatment = treatment, Lambda = Lambda, family = family,
            selected_r = selected$r_name, selected_h = selected$h_name, selection_criterion = selected$criterion,
            training_plugin_mean = selected$base_training_mean, training_plugin_standard_deviation = selected$base_training_standard_deviation,
            training_h_fraction = selected$h_training_fraction, validation_r_fraction = mean(selected$R_validation),
            validation_h_fraction = mean(selected$H_validation), stringsAsFactors = FALSE)
    }
    list(signals = do.call(cbind, output_signals), metadata = do.call(rbind, output_metadata), diagnostics = do.call(rbind,
        diagnostics))
}

fit_series_conditional_allocations <- function(training, validation, response, x_cols, spline_df = 5L,
    probability_epsilon = 1e-04) {
    if (nrow(training) != nrow(response)) {
        stop("The allocation response must align with training rows.")
    }
    output <- matrix(NA_real_, nrow(validation), ncol(response))
    colnames(output) <- colnames(response)
    constant <- apply(response, 2L, function(value) {
        max(value) - min(value) <= 10 * .Machine$double.eps
    })
    if (any(constant)) {
        output[, constant] <- matrix(rep(response[1L, constant], each = nrow(validation)), nrow = nrow(validation))
    }
    if (any(!constant)) {
        variable_columns <- x_cols[vapply(training[, x_cols, drop = FALSE], function(value) length(unique(value)) >=
            4L, logical(1))]
        binary_columns <- setdiff(x_cols, variable_columns)
        training_continuous <- if (length(variable_columns)) {
            make_spline_design(data = training, x_cols = variable_columns, spline_df = spline_df)
        }
        else NULL
        make_design <- function(newdata) {
            continuous <- if (is.null(training_continuous)) {
                matrix(1, nrow(newdata), 1L, dimnames = list(NULL, "intercept"))
            }
            else training_continuous$predict_matrix(newdata)
            if (!length(binary_columns))
                return(continuous)
            cbind(continuous, as.matrix(newdata[, binary_columns, drop = FALSE]))
        }
        design <- list(matrix = make_design(training), predict_matrix = make_design)
        fit <- stats::lm.fit(x = design$matrix, y = response[, !constant, drop = FALSE])
        coefficients <- fit$coefficients
        coefficients[!is.finite(coefficients)] <- 0
        prediction <- design$predict_matrix(validation) %*% coefficients
        output[, !constant] <- pmin(1 - probability_epsilon, pmax(probability_epsilon, prediction))
    }
    output
}

fit_witness_allocation_means <- function(training, validation, response_training, treatment, x_cols,
    spline_df) {
    output <- list()
    for (z_value in 0:1) {
        local <- training$T_obs == treatment & training$Z == z_value
        if (sum(local) <= spline_df * length(x_cols) + 3L) {
            stop("Too few observations for a learned-witness allocation regression.")
        }
        output[[paste0("train_z", z_value)]] <- fit_series_conditional_allocations(training = training[local,
            , drop = FALSE], validation = training, response = response_training[local, , drop = FALSE],
            x_cols = x_cols, spline_df = spline_df)
        output[[paste0("validation_z", z_value)]] <- fit_series_conditional_allocations(training = training[local,
            , drop = FALSE], validation = validation, response = response_training[local, , drop = FALSE],
            x_cols = x_cols, spline_df = spline_df)
    }
    output
}

gauss_legendre_rule <- function(n_nodes = 80, lower = -2, upper = 2) {
    if (n_nodes < 2) {
        stop("n_nodes must be at least 2.")
    }
    index <- seq_len(n_nodes - 1)
    beta <- index/sqrt(4 * index^2 - 1)
    jacobi <- matrix(0, nrow = n_nodes, ncol = n_nodes)
    jacobi[cbind(index, index + 1)] <- beta
    jacobi[cbind(index + 1, index)] <- beta
    eig <- eigen(jacobi, symmetric = TRUE)
    ordering <- order(eig$values)
    nodes_standard <- eig$values[ordering]
    weights_standard <- 2 * eig$vectors[1, ordering]^2
    list(nodes = (upper - lower) * nodes_standard/2 + (upper + lower)/2, weights = (upper - lower) *
        weights_standard/2)
}

learned_witness_probability_fit <- function(training, validation, x_cols, learner, positivity_epsilon,
    treatment_spline_df, gam_basis_dimension, integrated_confounder_strength, integrated_quadrature_nodes,
    super_learner_candidates = NULL, super_learner_cv_folds = 3L, super_learner_ranger_num_trees = 500L,
    super_learner_hal_max_degree = 2L, library_seed = 1L) {
    if (learner == "oracle") {
        required <- c("e_true", "pi0_true", "pi1_true")
        if (!all(required %in% names(training)) || !all(required %in% names(validation))) {
            stop("Oracle learned-witness testing requires true probability columns.")
        }
        return(list(training = list(e = training$e_true, pi0 = training$pi0_true, pi1 = training$pi1_true),
            validation = list(e = validation$e_true, pi0 = validation$pi0_true, pi1 = validation$pi1_true),
            fit = NULL))
    }
    if (learner == "glm_main") {
        formula_for <- function(response) stats::reformulate(x_cols, response)
        instrument_fit <- stats::glm(formula_for("Z"), data = training, family = stats::binomial())
        treatment_fit_z0 <- stats::glm(formula_for("T_obs"), data = training[training$Z == 0, , drop = FALSE],
            family = stats::binomial())
        treatment_fit_z1 <- stats::glm(formula_for("T_obs"), data = training[training$Z == 1, , drop = FALSE],
            family = stats::binomial())
        predict_probability <- function(newdata) list(e = clip_probability(1 - stats::predict(instrument_fit,
            newdata, type = "response"), positivity_epsilon), pi0 = clip_probability(stats::predict(treatment_fit_z0,
            newdata, type = "response"), positivity_epsilon), pi1 = clip_probability(stats::predict(treatment_fit_z1,
            newdata, type = "response"), positivity_epsilon))
        return(list(training = predict_probability(training), validation = predict_probability(validation),
            fit = list(instrument = instrument_fit, treatment_z0 = treatment_fit_z0, treatment_z1 = treatment_fit_z1)))
    }
    if (learner %in% c("binary_glm_spline", "binary_gam", "binary_glmnet", "binary_library_cv", "binary_super_learner")) {
        if (!exists("binary_fit_regression_with_learner", mode = "function")) {
            stop("Source R/binary_outcome_bounds.R for binary nuisance learners.")
        }
        binary_learner <- sub("^binary_", "", learner)
        fit_one <- function(local_data, response, seed_offset) {
            binary_fit_regression_with_learner(local_data, response, x_cols, if (binary_learner == "gam")
                gam_basis_dimension
            else treatment_spline_df, learner = binary_learner, super_learner_candidates = super_learner_candidates,
                super_learner_cv_folds = super_learner_cv_folds, super_learner_ranger_num_trees = super_learner_ranger_num_trees,
                super_learner_hal_max_degree = super_learner_hal_max_degree, seed = as.integer(library_seed) +
                  seed_offset)
        }
        instrument_fit <- fit_one(training, "Z", 1L)
        treatment_fit_z0 <- fit_one(training[training$Z == 0, , drop = FALSE], "T_obs", 2L)
        treatment_fit_z1 <- fit_one(training[training$Z == 1, , drop = FALSE], "T_obs", 3L)
        predict_probability <- function(newdata) list(e = clip_probability(1 - instrument_fit$predict(newdata),
            positivity_epsilon), pi0 = clip_probability(treatment_fit_z0$predict(newdata), positivity_epsilon),
            pi1 = clip_probability(treatment_fit_z1$predict(newdata), positivity_epsilon))
        return(list(training = predict_probability(training), validation = predict_probability(validation),
            fit = list(instrument = instrument_fit, treatment_z0 = treatment_fit_z0, treatment_z1 = treatment_fit_z1)))
    }
    fit <- fit_nuisance_models(train = training, x_cols = x_cols, positivity_epsilon = positivity_epsilon,
        learner = learner, treatment_spline_df = treatment_spline_df, gam_basis_dimension = gam_basis_dimension,
        integrated_confounder_strength = integrated_confounder_strength, integrated_quadrature_nodes = integrated_quadrature_nodes,
        super_learner_candidates = super_learner_candidates, super_learner_cv_folds = super_learner_cv_folds,
        super_learner_ranger_num_trees = super_learner_ranger_num_trees, super_learner_hal_max_degree = super_learner_hal_max_degree,
        library_seed = as.integer(library_seed))
    list(training = fit$predict(training), validation = fit$predict(validation), fit = fit)
}

make_additive_spline_formula <- function(response, predictors, spline_df = 5, data = NULL) {
    terms <- vapply(predictors, function(name) {
        use_linear <- !is.null(data) && length(unique(data[[name]][is.finite(data[[name]])])) <= spline_df +
            1L
        if (use_linear)
            name
        else sprintf("splines::ns(%s, df=%d)", name, spline_df)
    }, character(1))
    stats::as.formula(paste(response, "~", paste(terms, collapse = " + ")))
}

make_arm_lambda_signal_blocks <- function(data, H, R, m0, m1, e, pi0, pi1, t_value, Lambda0, Lambda1) {
    n <- nrow(data)
    p0 <- if (t_value == 1L)
        pi0
    else 1 - pi0
    p1 <- if (t_value == 1L)
        pi1
    else 1 - pi1
    primitive <- compute_dual_msm_moment_primitives(p0, p1, Lambda0, Lambda1)
    residual_p0 <- as.numeric(data$T_obs == t_value) - p0
    residual_p1 <- as.numeric(data$T_obs == t_value) - p1
    score_p0 <- as.numeric(data$Z == 0)/e * residual_p0
    score_p1 <- as.numeric(data$Z == 1)/(1 - e) * residual_p1
    score_m0 <- as.numeric(data$T_obs == t_value & data$Z == 0)/(e * p0) * (R - m0)
    score_m1 <- as.numeric(data$T_obs == t_value & data$Z == 1)/((1 - e) * p1) * (R - m1)
    a <- 1/Lambda0
    build_family <- function(base, dm0, dm1, dp0, dp1, family) {
        unweighted <- base + dm0 * score_m0 + dm1 * score_m1 + dp0 * score_p0 + dp1 * score_p1
        blocks <- lapply(seq_len(ncol(unweighted)), function(index) {
            H * unweighted[, index]
        })
        matrix <- do.call(cbind, blocks)
        metadata <- expand.grid(h = colnames(H), r = colnames(R), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        metadata <- metadata[order(match(metadata$r, colnames(R)), match(metadata$h, colnames(H))), ,
            drop = FALSE]
        metadata$family <- family
        list(signals = matrix, metadata = metadata)
    }
    families <- list(pointwise_lower = build_family(base = primitive$gamma_lower * m0 - m1, dm0 = primitive$gamma_lower,
        dm1 = -1, dp0 = primitive$d_gamma_lower_p0 * m0, dp1 = primitive$d_gamma_lower_p1 * m0, family = "pointwise_lower"),
        pointwise_upper = build_family(base = m1 - primitive$gamma_upper * m0, dm0 = -primitive$gamma_upper,
            dm1 = 1, dp0 = -primitive$d_gamma_upper_p0 * m0, dp1 = -primitive$d_gamma_upper_p1 * m0,
            family = "pointwise_upper"), normalization_lower = build_family(base = primitive$c_lower *
            m1 + (primitive$b - a) * m0 - (1 - a), dm0 = primitive$b - a, dm1 = primitive$c_lower, dp0 = primitive$d_c_lower_p0 *
            m1 + primitive$d_b_p0 * m0, dp1 = primitive$d_c_lower_p1 * m1, family = "normalization_lower"),
        normalization_upper = build_family(base = -primitive$c_upper * m1 + (Lambda0 - primitive$b) *
            m0 - (Lambda0 - 1), dm0 = Lambda0 - primitive$b, dm1 = -primitive$c_upper, dp0 = -primitive$d_c_upper_p0 *
            m1 - primitive$d_b_p0 * m0, dp1 = -primitive$d_c_upper_p1 * m1, family = "normalization_upper"))
    signals <- do.call(cbind, lapply(families, `[[`, "signals"))
    metadata <- do.call(rbind, lapply(families, `[[`, "metadata"))
    metadata$treatment <- t_value
    metadata$Lambda0 <- Lambda0
    metadata$Lambda1 <- Lambda1
    rownames(metadata) <- NULL
    list(signals = signals, metadata = metadata)
}

make_localization_cutoffs <- function(base, family, upper_tail_fractions = c(0.1, 0.2, 0.4)) {
    if (!grepl("^normalization_", family))
        return(c(all = -Inf))
    local <- base[is.finite(base)]
    cutoffs <- c(positive = 0, stats::setNames(as.numeric(stats::quantile(local, probs = 1 - upper_tail_fractions,
        type = 8, names = FALSE)), paste0("top", format(100 * upper_tail_fractions, trim = TRUE))), all = -Inf)
    cutoffs[!duplicated(round(cutoffs, 12))]
}

make_spline_design <- function(data, x_cols = "X", spline_df = 5) {
    non_numeric <- x_cols[!vapply(data[, x_cols, drop = FALSE], is.numeric, logical(1))]
    if (length(non_numeric) > 0) {
        stop("All x_cols must be numeric for the spline learners: ", paste(non_numeric, collapse = ", "))
    }
    basis_specs <- vector("list", length(x_cols))
    names(basis_specs) <- x_cols
    training_parts <- list(intercept = rep(1, nrow(data)))
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
        new_parts <- list(intercept = rep(1, nrow(newdata)))
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

make_tz_cluster_stratified_folds <- function(data, cluster_id, K = 5L, seed = 1L) {
    n <- nrow(data)
    K <- as.integer(K)
    if (length(cluster_id) != n || anyNA(cluster_id)) {
        stop("cluster_id must be complete and have one value per observation.")
    }
    if (!is.finite(K) || K < 2L)
        stop("K must be at least two.")
    cluster_id <- as.character(cluster_id)
    cluster_rows <- split(seq_len(n), cluster_id)
    if (length(cluster_rows) < K) {
        stop("There must be at least K independent clusters.")
    }
    cluster_summary <- lapply(cluster_rows, function(index) {
        local_z <- unique(data$Z[index])
        if (length(local_z) != 1L) {
            stop("Every cluster must have one common instrument value.")
        }
        c(Z = local_z, size = length(index), treated = sum(data$T_obs[index] == 1))
    })
    cluster_summary <- as.data.frame(do.call(rbind, cluster_summary))
    cluster_summary$cluster <- rownames(cluster_summary)
    stratum <- interaction(cluster_summary$Z, cluster_summary$size, cluster_summary$treated, drop = TRUE)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed)
        old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
        if (had_seed) {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)
    set.seed(seed)
    cluster_fold <- integer(nrow(cluster_summary))
    for (level in levels(stratum)) {
        index <- which(stratum == level)
        index <- sample(index, length(index), replace = FALSE)
        cluster_fold[index] <- rep(seq_len(K), length.out = length(index))
    }
    fold_id <- cluster_fold[match(cluster_id, cluster_summary$cluster)]
    if (anyNA(fold_id) || any(fold_id < 1L | fold_id > K)) {
        stop("Cluster-stratified fold construction failed.")
    }
    if (any(vapply(split(fold_id, cluster_id), function(value) length(unique(value)) != 1L, logical(1)))) {
        stop("At least one cluster was split across folds.")
    }
    cell <- table(fold_id, data$T_obs, data$Z)
    if (any(cell == 0L)) {
        stop("A cluster-stratified validation fold has an empty T-by-Z cell.")
    }
    as.integer(fold_id)
}

make_tz_stratified_folds <- function(data, K = 5, seed = 1) {
    check_binary(data$T_obs, "T_obs")
    check_binary(data$Z, "Z")
    if (K < 2) {
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

make_witness_cutoffs <- function(score, active, upper_tail_fractions = c(0.05, 0.1, 0.2, 0.4)) {
    local <- score[active & is.finite(score)]
    if (length(local) < 20L)
        stop("Too few observations to learn a witness.")
    cutoffs <- c(zero = 0, stats::setNames(as.numeric(stats::quantile(local, probs = 1 - upper_tail_fractions,
        type = 8, names = FALSE)), paste0("top", format(100 * upper_tail_fractions, trim = TRUE))))
    cutoffs[!duplicated(round(cutoffs, 12))]
}

make_z_stratified_folds <- function(data, K = 3, seed = 1) {
    if (K < 2) {
        stop("K must be at least 2.")
    }
    set.seed(seed)
    fold_id <- integer(nrow(data))
    for (z_value in c(0, 1)) {
        index <- which(data$Z == z_value)
        if (length(index) < K) {
            stop(sprintf("Only %d observations with Z=%d, fewer than K=%d.", length(index), z_value,
                K))
        }
        fold_id[index] <- sample(rep(seq_len(K), length.out = length(index)))
    }
    fold_id
}

select_learned_witness <- function(family, response_training, response_validation, mean_fit, primitive_training,
    primitive_validation, Lambda, localization_tail_fractions) {
    candidates <- list()
    counter <- 0L
    for (column in seq_len(ncol(response_training))) {
        base_training <- compute_learned_family_base(family, primitive_training, mean_fit$train_z0[,
            column], mean_fit$train_z1[, column], Lambda)
        base_validation <- compute_learned_family_base(family, primitive_validation, mean_fit$validation_z0[,
            column], mean_fit$validation_z1[, column], Lambda)
        h_cutoffs <- make_localization_cutoffs(base_training, family, upper_tail_fractions = localization_tail_fractions)
        for (h_name in names(h_cutoffs)) {
            h_training <- as.numeric(base_training > h_cutoffs[h_name])
            active_fraction <- mean(h_training)
            if (active_fraction < 0.025)
                next
            plug_in <- h_training * base_training
            standard_deviation <- stats::sd(plug_in)
            criterion <- if (is.finite(standard_deviation) && standard_deviation > 1e-10) {
                sqrt(length(plug_in)) * mean(plug_in)/standard_deviation
            }
            else {
                -Inf
            }
            counter <- counter + 1L
            candidates[[counter]] <- list(column = column, r_name = colnames(response_training)[column],
                h_name = h_name, h_cutoff = unname(h_cutoffs[h_name]), criterion = criterion, h_validation = as.numeric(base_validation >
                  h_cutoffs[h_name]), base_training_mean = mean(plug_in), base_training_standard_deviation = standard_deviation,
                h_training_fraction = active_fraction)
        }
    }
    if (length(candidates) == 0L)
        stop("No learned witness candidate remained.")
    criterion <- vapply(candidates, `[[`, numeric(1), "criterion")
    selected <- candidates[[which.max(criterion)]]
    selected$R_validation <- response_validation[, selected$column, drop = FALSE]
    selected$m0_validation <- mean_fit$validation_z0[, selected$column, drop = FALSE]
    selected$m1_validation <- mean_fit$validation_z1[, selected$column, drop = FALSE]
    selected$H_validation <- matrix(selected$h_validation, ncol = 1L, dimnames = list(NULL, paste0("learned_",
        family)))
    selected
}

test_learned_witness_dual_msm_compatibility <- function(data, lambda_grid, x_cols = "X", folds = 5L,
    primitive_learner = c("integrated_parametric", "glm_spline", "glm_main", "gam", "ranger", "library_cv",
        "super_learner", "binary_glm_spline", "binary_gam", "binary_glmnet", "binary_library_cv", "binary_super_learner",
        "oracle"), gamma_learner = c("rkhs", "rkhs_cv", "binary_glm_spline", "binary_gam", "binary_glmnet",
        "binary_library_cv", "binary_super_learner", "discrete_multinomial_spline", "oracle"), positivity_epsilon = 0.05,
    treatment_spline_df = 5L, gam_basis_dimension = 6L, allocation_spline_df = 5L, integrated_confounder_strength = 1,
    integrated_quadrature_nodes = 40L, super_learner_candidates = NULL, super_learner_cv_folds = 3L,
    super_learner_ranger_num_trees = 500L, super_learner_hal_max_degree = 2L, gamma_lambda_by_t = c(`0` = 0.03162278,
        `1` = 0.01), gamma_bandwidth_multiplier = 4, gamma_cv_lambda_grid = function(n_training, treatment) {
        10^seq(-1, 0.5, by = 0.5)/sqrt(n_training)
    }, gamma_cv_bandwidth_multipliers = c(1, 2, 4, 8), gamma_cv_folds = 3L, gamma_cv_repeats = 2L, gamma_cv_selection_rule = c("one_se",
        "minimum"), gamma_nystrom_rank = 80L, gamma_kernel = "gaussian", gamma_discrete_support = NULL,
    gamma_discrete_df_grid = 1:3, gamma_discrete_cv_folds = 3L, gamma_discrete_pmf_epsilon = 1e-05, gamma_discrete_decay = 1e-04,
    witness_tail_fractions = c(0.05, 0.1, 0.2, 0.4), localization_tail_fractions = c(0.1, 0.2, 0.4),
    witness_shapes = c("indicator", "linear", "sqrt"), bootstrap_replications = 999L, confidence_level = 0.95,
    seed = 1L, cluster_id = NULL, progress = FALSE) {
    check_learned_witness_dependencies()
    validate_conservative_test_data(data, x_cols)
    primitive_learner <- match.arg(primitive_learner)
    gamma_learner <- match.arg(gamma_learner)
    gamma_cv_selection_rule <- match.arg(gamma_cv_selection_rule)
    witness_shapes <- unique(as.character(witness_shapes))
    if (!length(witness_shapes) || !all(witness_shapes %in% c("indicator", "linear", "sqrt"))) {
        stop("witness_shapes must use indicator, linear, and/or sqrt.")
    }
    lambda_grid <- sort(unique(as.numeric(lambda_grid)))
    if (!length(lambda_grid) || any(!is.finite(lambda_grid)) || any(lambda_grid < 1)) {
        stop("lambda_grid must contain finite values of at least one.")
    }
    if (gamma_learner == "oracle") {
        required <- c("gamma0_oracle", "gamma1_oracle")
        if (!all(required %in% names(data))) {
            stop("Oracle gamma witnesses require gamma0_oracle and gamma1_oracle.")
        }
    }
    clustered_inference <- !is.null(cluster_id)
    if (is.null(cluster_id)) {
        cluster_id <- seq_len(nrow(data))
    }
    else if (length(cluster_id) != nrow(data) || anyNA(cluster_id)) {
        stop("cluster_id must be complete and have one value per data row.")
    }
    fold_id <- if (clustered_inference) {
        make_tz_cluster_stratified_folds(data, cluster_id = cluster_id, K = folds, seed = seed)
    }
    else {
        make_tz_stratified_folds(data, K = folds, seed = seed)
    }
    signal_rows <- vector("list", folds)
    selected_signal_rows <- vector("list", folds)
    aggregate_signal_rows <- vector("list", folds)
    component_selected_signal_rows <- vector("list", folds)
    component_aggregate_signal_rows <- vector("list", folds)
    metadata_rows <- vector("list", folds)
    diagnostic_rows <- vector("list", folds)
    gamma_diagnostics <- list()
    gamma_cv_diagnostics <- list()
    gamma_cv_resample_diagnostics <- list()
    probability_diagnostics <- list()
    gamma_counter <- 0L
    gamma_cv_counter <- 0L
    for (fold in seq_len(folds)) {
        if (isTRUE(progress)) {
            cat("Starting learned-witness outer fold ", fold, "/", folds, ".\n", sep = "")
        }
        training_index <- which(fold_id != fold)
        validation_index <- which(fold_id == fold)
        training <- data[training_index, , drop = FALSE]
        validation <- data[validation_index, , drop = FALSE]
        probability <- learned_witness_probability_fit(training, validation, x_cols, primitive_learner,
            positivity_epsilon, treatment_spline_df, gam_basis_dimension, integrated_confounder_strength,
            integrated_quadrature_nodes, super_learner_candidates = super_learner_candidates, super_learner_cv_folds = super_learner_cv_folds,
            super_learner_ranger_num_trees = super_learner_ranger_num_trees, super_learner_hal_max_degree = super_learner_hal_max_degree,
            library_seed = seed + 10000L * fold)
        if (!is.null(probability$fit) && is.data.frame(probability$fit$diagnostics)) {
            local_probability_diagnostics <- probability$fit$diagnostics
            local_probability_diagnostics$outer_fold <- fold
            probability_diagnostics[[length(probability_diagnostics) + 1L]] <- local_probability_diagnostics
        }
        fold_signals <- fold_metadata <- fold_diagnostics <- list()
        counter <- 0L
        for (treatment in 0:1) {
            n_treatment_training <- sum(training$T_obs == treatment)
            if (gamma_learner == "oracle") {
                gamma_training <- training[[paste0("gamma", treatment, "_oracle")]]
                gamma_validation <- validation[[paste0("gamma", treatment, "_oracle")]]
                gamma_fit <- NULL
            }
            else if (gamma_learner == "discrete_multinomial_spline") {
                gamma_fit <- fit_learned_witness_discrete_gamma(training = training, validation = validation,
                  treatment = treatment, x_cols = x_cols, support = if (is.null(gamma_discrete_support)) {
                    sort(unique(data$Y))
                  }
                  else {
                    gamma_discrete_support
                  }, df_grid = gamma_discrete_df_grid, cv_folds = gamma_discrete_cv_folds, pmf_epsilon = gamma_discrete_pmf_epsilon,
                  multinomial_decay = gamma_discrete_decay, seed = seed + 10000L * fold + 1000L * treatment)
                gamma_training <- gamma_fit$training
                gamma_validation <- gamma_fit$validation
            }
            else if (grepl("^binary_", gamma_learner)) {
                gamma_fit <- fit_learned_witness_binary_gamma(training, validation, treatment, x_cols,
                  learner = sub("^binary_", "", gamma_learner), spline_df = if (gamma_learner == "binary_gam") {
                    gam_basis_dimension
                  }
                  else treatment_spline_df, super_learner_candidates = super_learner_candidates, super_learner_cv_folds = super_learner_cv_folds,
                  super_learner_ranger_num_trees = super_learner_ranger_num_trees, super_learner_hal_max_degree = super_learner_hal_max_degree,
                  seed = seed + 10000L * fold + 1000L * treatment)
                gamma_training <- gamma_fit$training
                gamma_validation <- gamma_fit$validation
            }
            else {
                resolved_nystrom_rank <- if (is.function(gamma_nystrom_rank)) {
                  gamma_nystrom_rank(n_training = n_treatment_training, treatment = treatment)
                }
                else if (is.null(gamma_nystrom_rank)) {
                  NULL
                }
                else {
                  gamma_nystrom_rank
                }
                if (!is.null(resolved_nystrom_rank)) {
                  if (length(resolved_nystrom_rank) != 1L || !is.finite(resolved_nystrom_rank) || resolved_nystrom_rank <
                    1) {
                    stop("The resolved gamma Nystrom rank must be positive.")
                  }
                  resolved_nystrom_rank <- min(as.integer(ceiling(resolved_nystrom_rank)), n_treatment_training)
                }
                if (gamma_learner == "rkhs_cv") {
                  gamma_fit <- fit_gamma_rkhs(train = training, nuisance_train = probability$training,
                    t_value = treatment, x_cols = x_cols, lambda_grid = gamma_cv_lambda_grid, bandwidth_multipliers = gamma_cv_bandwidth_multipliers,
                    cv_folds = gamma_cv_folds, cv_repeats = gamma_cv_repeats, selection_rule = gamma_cv_selection_rule,
                    backend = "nystrom", nystrom_rank = resolved_nystrom_rank, kernel_type = gamma_kernel,
                    seed = seed + 10000L * fold + 1000L * treatment)
                  gamma_cv_counter <- gamma_cv_counter + 1L
                  gamma_cv_diagnostics[[gamma_cv_counter]] <- transform(gamma_fit$cv_table, outer_fold = fold,
                    treatment = treatment, selection_rule = gamma_cv_selection_rule)
                  gamma_cv_resample_diagnostics[[gamma_cv_counter]] <- transform(gamma_fit$cv_resample_table,
                    outer_fold = fold, treatment = treatment, selection_rule = gamma_cv_selection_rule)
                }
                else {
                  gamma_lambda <- if (is.function(gamma_lambda_by_t)) {
                    gamma_lambda_by_t(n_training = n_treatment_training, treatment = treatment)
                  }
                  else {
                    gamma_lambda_by_t[as.character(treatment)]
                  }
                  if (length(gamma_lambda) != 1L || !is.finite(gamma_lambda) || gamma_lambda <= 0) {
                    stop("The resolved fixed gamma RKHS penalty must be one positive number.")
                  }
                  gamma_fit <- fit_gamma_rkhs_fixed(train = training, nuisance_train = probability$training,
                    t_value = treatment, x_cols = x_cols, lambda = gamma_lambda, bandwidth_multiplier = gamma_bandwidth_multiplier,
                    backend = "nystrom", nystrom_rank = resolved_nystrom_rank, kernel_type = gamma_kernel,
                    seed = seed + 10000L * fold + 1000L * treatment)
                }
                gamma_training <- gamma_fit$predict(training)$normalized
                gamma_validation <- gamma_fit$predict(validation)$normalized
            }
            gamma_counter <- gamma_counter + 1L
            gamma_diagnostics[[gamma_counter]] <- data.frame(fold = fold, treatment = treatment, learner = gamma_learner,
                lambda = if (is.null(gamma_fit) || is.null(gamma_fit$selected_lambda)) {
                  NA_real_
                }
                else gamma_fit$selected_lambda, lambda_constant = if (is.null(gamma_fit) || is.null(gamma_fit$selected_lambda)) {
                  NA_real_
                }
                else gamma_fit$selected_lambda * sqrt(n_treatment_training), bandwidth_multiplier = if (is.null(gamma_fit) ||
                  is.null(gamma_fit$selected_bandwidth_multiplier)) {
                  NA_real_
                }
                else gamma_fit$selected_bandwidth_multiplier, selected_bandwidth = if (is.null(gamma_fit) ||
                  is.null(gamma_fit$selected_bandwidth)) {
                  NA_real_
                }
                else gamma_fit$selected_bandwidth, nystrom_rank_requested = if (is.null(gamma_fit) ||
                  is.null(gamma_fit$nystrom_rank_requested)) {
                  NA_integer_
                }
                else gamma_fit$nystrom_rank_requested, nystrom_effective_rank = if (is.null(gamma_fit) ||
                  is.null(gamma_fit$nystrom_effective_rank)) {
                  NA_integer_
                }
                else gamma_fit$nystrom_effective_rank, selected_df_z0 = if (is.null(gamma_fit) || is.null(gamma_fit$selected_df_z0)) {
                  NA_integer_
                }
                else gamma_fit$selected_df_z0, selected_df_z1 = if (is.null(gamma_fit) || is.null(gamma_fit$selected_df_z1)) {
                  NA_integer_
                }
                else gamma_fit$selected_df_z1, discrete_pmf_epsilon = if (is.null(gamma_fit) || is.null(gamma_fit$pmf_epsilon)) {
                  NA_real_
                }
                else gamma_fit$pmf_epsilon, training_mean = mean(gamma_training), validation_mean = mean(gamma_validation),
                training_oracle_bias = if (paste0("gamma", treatment, "_oracle") %in% names(training)) {
                  active <- training$T_obs == treatment
                  mean(gamma_training[active] - training[[paste0("gamma", treatment, "_oracle")]][active])
                }
                else NA_real_, training_oracle_rmse = if (paste0("gamma", treatment, "_oracle") %in%
                  names(training)) {
                  active <- training$T_obs == treatment
                  sqrt(mean((gamma_training[active] - training[[paste0("gamma", treatment, "_oracle")]][active])^2))
                }
                else NA_real_, validation_oracle_bias = if (paste0("gamma", treatment, "_oracle") %in%
                  names(validation)) {
                  active <- validation$T_obs == treatment
                  mean(gamma_validation[active] - validation[[paste0("gamma", treatment, "_oracle")]][active])
                }
                else NA_real_, validation_oracle_rmse = if (paste0("gamma", treatment, "_oracle") %in%
                  names(validation)) {
                  active <- validation$T_obs == treatment
                  sqrt(mean((gamma_validation[active] - validation[[paste0("gamma", treatment, "_oracle")]][active])^2))
                }
                else NA_real_, stringsAsFactors = FALSE)
            if (isTRUE(progress)) {
                cat("Completed gamma nuisance for fold ", fold, ", treatment ", treatment, ".\n", sep = "")
            }
            for (Lambda in lambda_grid) {
                counter <- counter + 1L
                local <- fit_fold_learned_witness_signals(training, validation, probability$training,
                  probability$validation, gamma_training, gamma_validation, treatment, Lambda, x_cols,
                  allocation_spline_df, witness_tail_fractions, localization_tail_fractions, witness_shapes)
                fold_signals[[counter]] <- local$signals
                fold_metadata[[counter]] <- local$metadata
                local$diagnostics$fold <- fold
                fold_diagnostics[[counter]] <- local$diagnostics
            }
        }
        local_signals <- do.call(cbind, fold_signals)
        local_metadata <- do.call(rbind, fold_metadata)
        moment_key <- paste(local_metadata$Lambda0, local_metadata$treatment, local_metadata$family,
            sep = "::")
        ordering <- order(moment_key)
        local_signals <- local_signals[, ordering, drop = FALSE]
        local_metadata <- local_metadata[ordering, , drop = FALSE]
        local_diagnostics <- do.call(rbind, fold_diagnostics)
        if (fold == 1L) {
            reference_key <- moment_key[ordering]
            metadata_reference <- local_metadata
            metadata_reference$r <- "cross_fitted_learned"
            metadata_reference$h <- "cross_fitted_learned"
        }
        else if (!identical(moment_key[ordering], reference_key)) {
            stop("Learned-witness moment columns differ across folds.")
        }
        diagnostic_key <- paste(local_diagnostics$Lambda, local_diagnostics$treatment, local_diagnostics$family,
            sep = "::")
        metadata_key <- paste(local_metadata$Lambda0, local_metadata$treatment, local_metadata$family,
            sep = "::")
        criterion <- local_diagnostics$selection_criterion[match(metadata_key, diagnostic_key)]
        training_standard_deviation <- local_diagnostics$training_plugin_standard_deviation[match(metadata_key,
            diagnostic_key)]
        if (anyNA(criterion))
            stop("Could not match learned aggregation weights.")
        if (anyNA(training_standard_deviation)) {
            stop("Could not match learned aggregation standard deviations.")
        }
        selected_columns <- vapply(lambda_grid, function(Lambda) {
            index <- which(local_metadata$Lambda0 == Lambda)
            index[which.max(criterion[index])]
        }, integer(1))
        selected_matrix <- local_signals[, selected_columns, drop = FALSE]
        colnames(selected_matrix) <- paste0("Lambda_", lambda_grid)
        selected_key <- metadata_key[selected_columns]
        local_diagnostics$selected_across_families <- diagnostic_key %in% selected_key
        aggregate_columns <- lapply(lambda_grid, function(Lambda) {
            index <- which(local_metadata$Lambda0 == Lambda)
            weight <- pmax(criterion[index], 0)/pmax(training_standard_deviation[index], 1e-10)
            if (!any(weight > 0))
                weight <- rep(1, length(index))
            weight <- weight/sum(weight)
            as.vector(local_signals[, index, drop = FALSE] %*% weight)
        })
        aggregate_matrix <- do.call(cbind, aggregate_columns)
        colnames(aggregate_matrix) <- paste0("Lambda_", lambda_grid)
        component_names <- c("pointwise", "normalization")
        component_selected_columns <- list()
        component_aggregate_columns <- list()
        component_counter <- 0L
        for (Lambda in lambda_grid) {
            for (component in component_names) {
                component_counter <- component_counter + 1L
                index <- which(local_metadata$Lambda0 == Lambda & startsWith(local_metadata$family, paste0(component,
                  "_")))
                if (!length(index)) {
                  stop("No learned-witness moments remained for component ", component)
                }
                selected_index <- index[which.max(criterion[index])]
                component_selected_columns[[component_counter]] <- local_signals[, selected_index]
                weight <- pmax(criterion[index], 0)/pmax(training_standard_deviation[index], 1e-10)
                if (!any(weight > 0)) {
                  weight <- numeric(length(index))
                  weight[which.max(criterion[index])] <- 1
                }
                else {
                  weight <- weight/sum(weight)
                }
                component_aggregate_columns[[component_counter]] <- as.vector(local_signals[, index,
                  drop = FALSE] %*% weight)
            }
        }
        component_column_names <- unlist(lapply(lambda_grid, function(Lambda) {
            paste0("Lambda_", Lambda, "__", component_names)
        }))
        component_selected_matrix <- do.call(cbind, component_selected_columns)
        component_aggregate_matrix <- do.call(cbind, component_aggregate_columns)
        colnames(component_selected_matrix) <- component_column_names
        colnames(component_aggregate_matrix) <- component_column_names
        signal_rows[[fold]] <- data.frame(row_id = validation_index, local_signals, check.names = FALSE)
        selected_signal_rows[[fold]] <- data.frame(row_id = validation_index, selected_matrix, check.names = FALSE)
        aggregate_signal_rows[[fold]] <- data.frame(row_id = validation_index, aggregate_matrix, check.names = FALSE)
        component_selected_signal_rows[[fold]] <- data.frame(row_id = validation_index, component_selected_matrix,
            check.names = FALSE)
        component_aggregate_signal_rows[[fold]] <- data.frame(row_id = validation_index, component_aggregate_matrix,
            check.names = FALSE)
        diagnostic_rows[[fold]] <- local_diagnostics
        if (isTRUE(progress)) {
            cat("Completed learned-witness outer fold ", fold, "/", folds, ".\n", sep = "")
        }
    }
    stacked <- do.call(rbind, signal_rows)
    stacked <- stacked[order(stacked$row_id), , drop = FALSE]
    signals <- as.matrix(stacked[, -1L, drop = FALSE])
    selected_stacked <- do.call(rbind, selected_signal_rows)
    selected_stacked <- selected_stacked[order(selected_stacked$row_id), , drop = FALSE]
    selected_signals <- as.matrix(selected_stacked[, -1L, drop = FALSE])
    aggregate_stacked <- do.call(rbind, aggregate_signal_rows)
    aggregate_stacked <- aggregate_stacked[order(aggregate_stacked$row_id), , drop = FALSE]
    aggregate_signals <- as.matrix(aggregate_stacked[, -1L, drop = FALSE])
    component_selected_stacked <- do.call(rbind, component_selected_signal_rows)
    component_selected_stacked <- component_selected_stacked[order(component_selected_stacked$row_id),
        , drop = FALSE]
    component_selected_signals <- as.matrix(component_selected_stacked[, -1L, drop = FALSE])
    component_aggregate_stacked <- do.call(rbind, component_aggregate_signal_rows)
    component_aggregate_stacked <- component_aggregate_stacked[order(component_aggregate_stacked$row_id),
        , drop = FALSE]
    component_aggregate_signals <- as.matrix(component_aggregate_stacked[, -1L, drop = FALSE])
    calibration <- calibrate_learned_witness_maximum(signals, metadata_reference, bootstrap_replications,
        confidence_level, seed + 900000L, cluster_id = cluster_id)
    selected_metadata <- data.frame(h = "cross_fitted_learned", r = "cross_fitted_learned", family = "training_selected_family",
        treatment = -1L, Lambda0 = lambda_grid, Lambda1 = lambda_grid, stringsAsFactors = FALSE)
    selected_calibration <- calibrate_learned_witness_maximum(selected_signals, selected_metadata, bootstrap_replications,
        confidence_level, seed + 905000L, cluster_id = cluster_id)
    selected_summary <- selected_calibration$summary
    names(selected_summary)[names(selected_summary) == "statistic"] <- "selected_statistic"
    names(selected_summary)[names(selected_summary) == "critical_value"] <- "selected_critical_value"
    names(selected_summary)[names(selected_summary) == "p_value"] <- "selected_p_value"
    names(selected_summary)[names(selected_summary) == "reject"] <- "selected_reject"
    names(selected_summary)[names(selected_summary) == "ordered_reject"] <- "selected_ordered_reject"
    names(selected_summary)[names(selected_summary) == "grid_simultaneous_p_value"] <- "selected_grid_simultaneous_p_value"
    names(selected_summary)[names(selected_summary) == "grid_simultaneous_reject"] <- "selected_grid_simultaneous_reject"
    selected_summary$moments_tested <- NULL
    selected_summary <- selected_summary[, !startsWith(names(selected_summary), "cvm_"), drop = FALSE]
    aggregate_metadata <- data.frame(h = "cross_fitted_learned", r = "cross_fitted_learned", family = "directed_aggregate",
        treatment = -1L, Lambda0 = lambda_grid, Lambda1 = lambda_grid, stringsAsFactors = FALSE)
    aggregate_calibration <- calibrate_learned_witness_maximum(aggregate_signals, aggregate_metadata,
        bootstrap_replications, confidence_level, seed + 910000L, cluster_id = cluster_id)
    aggregate_summary <- aggregate_calibration$summary
    names(aggregate_summary)[names(aggregate_summary) == "statistic"] <- "aggregate_statistic"
    names(aggregate_summary)[names(aggregate_summary) == "critical_value"] <- "aggregate_critical_value"
    names(aggregate_summary)[names(aggregate_summary) == "p_value"] <- "aggregate_p_value"
    names(aggregate_summary)[names(aggregate_summary) == "reject"] <- "aggregate_reject"
    names(aggregate_summary)[names(aggregate_summary) == "ordered_reject"] <- "aggregate_ordered_reject"
    names(aggregate_summary)[names(aggregate_summary) == "grid_simultaneous_p_value"] <- "aggregate_grid_simultaneous_p_value"
    names(aggregate_summary)[names(aggregate_summary) == "grid_simultaneous_reject"] <- "aggregate_grid_simultaneous_reject"
    aggregate_summary$moments_tested <- NULL
    aggregate_summary <- aggregate_summary[, !startsWith(names(aggregate_summary), "cvm_"), drop = FALSE]
    component_metadata <- do.call(rbind, lapply(lambda_grid, function(Lambda) {
        data.frame(h = "cross_fitted_learned", r = "cross_fitted_learned", family = c("pointwise_component",
            "normalization_component"), treatment = -1L, Lambda0 = Lambda, Lambda1 = Lambda, stringsAsFactors = FALSE)
    }))
    component_selected_calibration <- calibrate_learned_witness_maximum(component_selected_signals, component_metadata,
        bootstrap_replications, confidence_level, seed + 915000L, cluster_id = cluster_id)
    component_selected_summary <- component_selected_calibration$summary
    names(component_selected_summary)[names(component_selected_summary) == "statistic"] <- "component_selected_statistic"
    names(component_selected_summary)[names(component_selected_summary) == "critical_value"] <- "component_selected_critical_value"
    names(component_selected_summary)[names(component_selected_summary) == "p_value"] <- "component_selected_p_value"
    names(component_selected_summary)[names(component_selected_summary) == "reject"] <- "component_selected_reject"
    names(component_selected_summary)[names(component_selected_summary) == "ordered_reject"] <- "component_selected_ordered_reject"
    names(component_selected_summary)[names(component_selected_summary) == "grid_simultaneous_p_value"] <- "component_selected_grid_simultaneous_p_value"
    names(component_selected_summary)[names(component_selected_summary) == "grid_simultaneous_reject"] <- "component_selected_grid_simultaneous_reject"
    component_selected_summary$moments_tested <- NULL
    component_selected_summary <- component_selected_summary[, !startsWith(names(component_selected_summary),
        "cvm_"), drop = FALSE]
    component_aggregate_calibration <- calibrate_learned_witness_maximum(component_aggregate_signals,
        component_metadata, bootstrap_replications, confidence_level, seed + 920000L, cluster_id = cluster_id)
    component_aggregate_summary <- component_aggregate_calibration$summary
    names(component_aggregate_summary)[names(component_aggregate_summary) == "statistic"] <- "component_aggregate_statistic"
    names(component_aggregate_summary)[names(component_aggregate_summary) == "critical_value"] <- "component_aggregate_critical_value"
    names(component_aggregate_summary)[names(component_aggregate_summary) == "p_value"] <- "component_aggregate_p_value"
    names(component_aggregate_summary)[names(component_aggregate_summary) == "reject"] <- "component_aggregate_reject"
    names(component_aggregate_summary)[names(component_aggregate_summary) == "ordered_reject"] <- "component_aggregate_ordered_reject"
    names(component_aggregate_summary)[names(component_aggregate_summary) == "grid_simultaneous_p_value"] <- "component_aggregate_grid_simultaneous_p_value"
    names(component_aggregate_summary)[names(component_aggregate_summary) == "grid_simultaneous_reject"] <- "component_aggregate_grid_simultaneous_reject"
    component_aggregate_summary$moments_tested <- NULL
    component_aggregate_summary <- component_aggregate_summary[, !startsWith(names(component_aggregate_summary),
        "cvm_"), drop = FALSE]
    calibration$summary <- merge(calibration$summary, selected_summary, by = c("Lambda0", "Lambda1"),
        sort = TRUE)
    calibration$summary <- merge(calibration$summary, aggregate_summary, by = c("Lambda0", "Lambda1"),
        sort = TRUE)
    calibration$summary <- merge(calibration$summary, component_selected_summary, by = c("Lambda0", "Lambda1"),
        sort = TRUE)
    calibration$summary <- merge(calibration$summary, component_aggregate_summary, by = c("Lambda0",
        "Lambda1"), sort = TRUE)
    rejected <- calibration$summary$Lambda0[calibration$summary$ordered_reject]
    cvm_rejected <- calibration$summary$Lambda0[calibration$summary$cvm_ordered_reject]
    selected_rejected <- calibration$summary$Lambda0[calibration$summary$selected_ordered_reject]
    component_selected_rejected <- calibration$summary$Lambda0[calibration$summary$component_selected_ordered_reject]
    component_aggregate_rejected <- calibration$summary$Lambda0[calibration$summary$component_aggregate_ordered_reject]
    list(summary = calibration$summary, moment_detail = calibration$moment_detail, selected_moment_detail = selected_calibration$moment_detail,
        component_selected_moment_detail = component_selected_calibration$moment_detail, component_aggregate_moment_detail = component_aggregate_calibration$moment_detail,
        largest_ordered_refuted_grid_value = if (length(rejected)) {
            max(rejected)
        } else NA_real_, largest_cvm_ordered_refuted_grid_value = if (length(cvm_rejected)) {
            max(cvm_rejected)
        } else NA_real_, largest_selected_ordered_refuted_grid_value = if (length(selected_rejected)) max(selected_rejected) else NA_real_,
        largest_component_selected_ordered_refuted_grid_value = if (length(component_selected_rejected)) {
            max(component_selected_rejected)
        } else NA_real_, largest_component_aggregate_ordered_refuted_grid_value = if (length(component_aggregate_rejected)) {
            max(component_aggregate_rejected)
        } else NA_real_, witness_diagnostics = do.call(rbind, diagnostic_rows), aggregate_moment_detail = aggregate_calibration$moment_detail,
        gamma_diagnostics = do.call(rbind, gamma_diagnostics), gamma_cv = if (length(gamma_cv_diagnostics)) {
            do.call(rbind, gamma_cv_diagnostics)
        } else data.frame(), gamma_cv_resamples = if (length(gamma_cv_resample_diagnostics)) {
            do.call(rbind, gamma_cv_resample_diagnostics)
        } else data.frame(), probability_diagnostics = if (length(probability_diagnostics)) {
            do.call(rbind, probability_diagnostics)
        } else data.frame(), fold_id = fold_id, settings = list(lambda_grid = lambda_grid, folds = folds,
            primitive_learner = primitive_learner, super_learner_candidates = super_learner_candidates,
            super_learner_cv_folds = super_learner_cv_folds, super_learner_ranger_num_trees = super_learner_ranger_num_trees,
            gamma_learner = gamma_learner, gamma_cv_selection_rule = gamma_cv_selection_rule, gamma_cv_folds = gamma_cv_folds,
            gamma_cv_repeats = gamma_cv_repeats, gamma_discrete_support = gamma_discrete_support, gamma_discrete_df_grid = gamma_discrete_df_grid,
            gamma_discrete_cv_folds = gamma_discrete_cv_folds, gamma_discrete_pmf_epsilon = gamma_discrete_pmf_epsilon,
            gamma_discrete_decay = gamma_discrete_decay, bootstrap_replications = bootstrap_replications,
            confidence_level = confidence_level, clustered_inference = clustered_inference, cluster_count = length(unique(cluster_id)),
            progress = isTRUE(progress), inferential_scope = paste("Honest cross-fitted learned-witness test of necessary compatibility",
                "conditions; not a uniform adaptive test over all witnesses.")))
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
        stop("Y and the test covariates must be finite numeric variables. Invalid: ", paste(numeric_columns[invalid],
            collapse = ", "), call. = FALSE)
    }
    invisible(TRUE)
}

