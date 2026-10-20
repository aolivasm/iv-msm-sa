# Main-text reproduction dependencies; published formulas retained.
fingerprint <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

generate <- function(t, n = config$n) {
    set.seed(t$seed)
    x <- matrix(qnorm(pnorm(-2) + runif(n * 5) * (pnorm(2) - pnorm(-2))), ncol = 5)
    colnames(x) <- paste0("X", 1:5)
    z <- rbinom(n, 1, 0.5)
    g <- config$settings[t$setting, ]
    h <- drop(x %*% c(0.14, 0.12, 0.1, 0.08, 0.06))
    tr <- rbinom(n, 1, plogis(h + g$a * (z - 0.5)))
    shape <- 1 + 0.6 * sin(x[, 1] * x[, 2]) + 0.4 * sin(x[, 3] * x[, 4])
    p <- plogis(-0.4 + 1.5 * tr + drop(x %*% c(0.2, -0.25, 0.3, 0.15, -0.2)) + (2 * tr - 1) * g$b * shape *
        z)
    data.frame(Y = rbinom(n, 1, p), T_obs = tr, Z = z, x)
}

new_session <- function() {
    engine <- new.env(parent = globalenv())
    sys.source("R/compatibility_test.R", engine)
    probability_original <- engine$learned_witness_probability_fit
    cache <- new.env(parent = emptyenv())
    engine$learned_witness_probability_fit <- function(...) {
        args <- list(...)
        key <- as.character(args$library_seed)
        hash <- fingerprint(args)
        if (exists(key, cache, inherits = FALSE)) {
            z <- get(key, cache)
            stopifnot(identical(z$hash, hash))
            return(z$value)
        }
        z <- do.call(probability_original, args)
        diagnostics <- data.frame(component = c("e", "pi0", "pi1"), training_clip_rate = vapply(c("e",
            "pi0", "pi1"), function(k) mean(z$training[[k]] <= config$positivity_epsilon | z$training[[k]] >=
            1 - config$positivity_epsilon), numeric(1)))
        z$fit <- list(diagnostics = diagnostics)
        assign(key, list(hash = hash, value = z), cache)
        z
    }
    forbidden <- function(...) stop("Forbidden: compatibility projection or q fitting in formal testing.")
    engine$project_gamma_to_compatible_interval <- forbidden
    engine$fit_learned_witness_binary_gamma <- function(...) stop("Classification gamma is forbidden.")
    engine$fit_learned_witness_discrete_gamma <- function(...) stop("PMF-ratio gamma is not this experiment.")
    previous <- NULL
    run <- function(data, method, seed, progress = TRUE) {
        stopifnot(method %in% config$methods)
        selection <- list()
        original_gamma <- fit_gamma_rkhs
        engine$fit_gamma_rkhs <- function(train, nuisance_train, t_value, x_cols, seed, ...) {
            fit <- {
                original_gamma(train = train, nuisance_train = nuisance_train, t_value = t_value, x_cols = x_cols,
                  seed = seed, lambda_grid = default_gamma_lambda_grid, bandwidth_multipliers = 1/sqrt(2),
                  cv_folds = 3L, cv_repeats = 2L, selection_rule = "minimum", backend = "nystrom", nystrom_rank = NULL)
            }
            record <- fit[intersect(names(fit), c("selected_lambda", "selected_level", "feature_count",
                "base_bandwidth", "selected_bandwidth", "selected_bandwidth_multiplier", "cv_table",
                "cv_resample_table", "rkhs_backend", "normalization_constant", "training_negative_rate",
                "clip_rate", "gradient"))]
            record$treatment <- t_value
            record$seed <- seed
            record$n_arm <- sum(train$T_obs == t_value)
            record$selected_C <- fit$selected_lambda * sqrt(record$n_arm)
            selection[[length(selection) + 1L]] <<- record
            fit
        }
        fit <- engine$test_learned_witness_dual_msm_compatibility(data = data, lambda_grid = config$Lambda,
            x_cols = paste0("X", 1:5), folds = config$folds, primitive_learner = "binary_super_learner",
            gamma_learner = "rkhs_cv", positivity_epsilon = config$positivity_epsilon, super_learner_candidates = config$primitives,
            super_learner_cv_folds = 3L, super_learner_ranger_num_trees = 500L, gamma_cv_lambda_grid = default_gamma_lambda_grid,
            gamma_cv_bandwidth_multipliers = 1/sqrt(2), gamma_cv_folds = 3L, gamma_cv_repeats = 2L, gamma_cv_selection_rule = "minimum",
            gamma_nystrom_rank = NULL, allocation_spline_df = config$allocation_df, bootstrap_replications = config$bootstrap,
            seed = seed, progress = progress)
        hashes <- vapply(sort(ls(cache)), function(key) {
            z <- get(key, cache)$value
            fingerprint(z[c("training", "validation")])
        }, character(1))
        paired <- list(data = fingerprint(data), folds = fingerprint(fit$fold_id), primitives = hashes)
        if (!is.null(previous))
            stopifnot(identical(previous, paired))
        previous <<- paired
        stopifnot(length(selection) == 2L * config$folds, length(hashes) == config$folds)
        fit$settings$gamma_learner <- method
        fit$gamma_diagnostics$learner <- method
        fit$settings$gamma_normalization <- "Positive prediction, training baseline-stratum mean normalization only; no compatibility projection or binary PMF reconstruction"
        list(summary = fit$summary, gamma_diagnostics = fit$gamma_diagnostics, gamma_selection = selection,
            probability_diagnostics = fit$probability_diagnostics, settings = fit$settings, paired = paired,
            fold_sizes = tabulate(fit$fold_id, config$folds), witness_candidate_counts = as.data.frame(table(fit$witness_diagnostics$family)),
            projection_applied = FALSE, q_fits = 0L)
    }
    list(run = run, engine = engine)
}

