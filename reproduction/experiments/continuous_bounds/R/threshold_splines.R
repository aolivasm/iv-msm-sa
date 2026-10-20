# Main-text reproduction dependencies; published formulas retained.
restricted_x5_contract <- function() list(spline_knots = c(-1.2, -0.4, 0.4, 1.2), boundary = c(-2, 2),
    q_penalty = "1/n_training; squared coefficients in full-rank training-whitened features", target = "unpenalized fixed-class pseudo-bound, NOT normalized-direct plug-in")

restricted_x5_features <- function(data, threshold_class, contract = restricted_x5_contract()) {
    z <- contract
    x_cols <- paste0("X", 1:5)
    stopifnot(all(x_cols %in% names(data)))
    if (threshold_class == "fixed_spline5") {
        return(do.call(cbind, lapply(x_cols, function(k) splines::ns(data[[k]], knots = z$spline_knots,
            Boundary.knots = z$boundary, intercept = FALSE))))
    }
    stop("Only the paper spline class is supported.")
}

restricted_x5_fit_q <- function(train, y, level, weight, threshold_class, n_training, contract = restricted_x5_contract()) {
    f <- restricted_x5_features(train, threshold_class, contract)
    center <- colMeans(f)
    centered <- sweep(f, 2, center, "-")
    s <- svd(centered, nu = 0)
    if (length(s$d) != ncol(f) || min(s$d) < max(s$d) * 1e-12)
        stop("Threshold design is not full rank; class will not be silently reduced.")
    rotation <- sweep(s$v, 2, s$d/sqrt(nrow(f)), "/")
    features <- centered %*% rotation
    penalty <- 1/as.integer(n_training)
    fit <- fit_weighted_quantile_ridge_admm(features, y, level, weight, penalty, rho = 0.01, maxit = 12000L,
        abs_tol = 1e-05, rel_tol = 1e-04)
    if (fit$convergence != 0L)
        fit <- fit_weighted_quantile_ridge_admm(features, y, level, weight, penalty, rho = 0.01, maxit = 60000L,
            abs_tol = 1e-05, rel_tol = 1e-04)
    raw <- as.vector(rotation %*% fit$coefficients[-1])
    coefficients <- c(fit$coefficients[1] - sum(center * raw), raw)
    list(coefficients = coefficients, threshold_class = threshold_class, threshold_contract = contract,
        convergence = fit$convergence, iterations = fit$iterations, penalty = penalty, primal_residual = fit$primal_residual,
        dual_residual = fit$dual_residual, objective = fit$objective, penalized_objective = fit$penalized_objective)
}

restricted_x5_predict_q <- function(model, newdata) {
    contract <- if (is.null(model$threshold_contract)) {
        restricted_x5_contract()
    }
    else model$threshold_contract
    as.vector(cbind(1, restricted_x5_features(newdata, model$threshold_class, contract)) %*% model$coefficients)
}

restricted_x5_provider <- function(threshold_class, store, contract = restricted_x5_contract()) {
    force(threshold_class)
    force(store)
    force(contract)
    function(mode, fold, treatment, ...) {
        force(mode)
        force(fold)
        force(treatment)
        list(thresholds = function(train, validation, nuisance_train, envelope_train, tau_train, ...) {
            use <- train$T_obs == treatment & train$Z == 0
            weights <- compute_target_odds(nuisance_train, treatment)[use] * envelope_train$D[use]/nuisance_train$e[use]
            lo <- restricted_x5_fit_q(train[use, , drop = FALSE], train$Y[use], 1 - tau_train[use], weights,
                threshold_class, nrow(train), contract)
            hi <- restricted_x5_fit_q(train[use, , drop = FALSE], train$Y[use], tau_train[use], weights,
                threshold_class, nrow(train), contract)
            key <- paste(mode, fold, treatment, sep = "_")
            assign(key, list(mode = mode, fold = fold, treatment = treatment, lower = lo, upper = hi),
                envir = store)
            lower_predict <- function(newdata) restricted_x5_predict_q(lo, newdata)
            upper_predict <- function(newdata) restricted_x5_predict_q(hi, newdata)
            list(q_lower_train = lower_predict(train), q_upper_train = upper_predict(train), q_lower_validation = lower_predict(validation),
                q_upper_validation = upper_predict(validation), lower_predict = lower_predict, upper_predict = upper_predict)
        })
    }
}

x5s_threshold_provider <- function(threshold_class, store) {
    e <- new.env(parent = globalenv())
    sys.source("R/threshold_splines.R", e)
    e$restricted_x5_fit_q <- function(train, y, level, weight, threshold_class, n_training, contract = e$restricted_x5_contract()) {
        f <- e$restricted_x5_features(train, threshold_class, contract)
        center <- colMeans(f)
        centered <- sweep(f, 2, center, "-")
        s <- svd(centered, nu = 0)
        if (length(s$d) != ncol(f) || min(s$d) < max(s$d) * 1e-12)
            stop("Threshold basis rank deficient.")
        rotation <- sweep(s$v, 2, s$d/sqrt(nrow(f)), "/")
        features <- centered %*% rotation
        penalty <- 1/as.integer(n_training)
        attempts <- list()
        for (rho in c(0.01, 0.1, 1, 0.001)) {
            fit <- fit_weighted_quantile_ridge_admm(features, y, level, weight, penalty, rho = rho, maxit = 60000L,
                abs_tol = 1e-05, rel_tol = 1e-04)
            attempts[[length(attempts) + 1]] <- list(rho = rho, convergence = fit$convergence, iterations = fit$iterations,
                primal_residual = fit$primal_residual, dual_residual = fit$dual_residual)
            if (fit$convergence == 0L)
                break
        }
        if (fit$convergence != 0L)
            stop("Threshold optimization failed certification after rho retries.")
        raw <- drop(rotation %*% fit$coefficients[-1])
        list(coefficients = c(fit$coefficients[1] - sum(center * raw), raw), threshold_class = threshold_class,
            threshold_contract = contract, convergence = fit$convergence, iterations = fit$iterations,
            penalty = penalty, primal_residual = fit$primal_residual, dual_residual = fit$dual_residual,
            objective = fit$objective, penalized_objective = fit$penalized_objective, attempts = attempts)
    }
    e$restricted_x5_provider(threshold_class, store)
}

