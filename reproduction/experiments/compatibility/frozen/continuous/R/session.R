# Main-text reproduction dependencies; published formulas retained.
new_session <- function() {
    engine <- new.env(parent = globalenv())
    sys.source("R/compatibility_test.R", engine)
    original <- engine$learned_witness_probability_fit
    cache <- new.env(parent = emptyenv())
    engine$learned_witness_probability_fit <- function(...) {
        args <- list(...)
        key <- as.character(args$library_seed)
        if (exists(key, cache, inherits = FALSE)) {
            z <- get(key, cache)
            stopifnot(identical(args, z$args))
            return(z$fit)
        }
        f <- do.call(original, args)
        assign(key, list(args = args, fit = f), cache)
        f
    }
    engine$project_gamma_to_compatible_interval <- function(...) stop("Forbidden compatibility projection in testing.")
    last_folds <- NULL
    run <- function(d, method, seed, B = config$bootstrap) {
        provider <- x5_adaptive_provider(method)
        engine$fit_gamma_rkhs <- function(train, nuisance_train, t_value, x_cols, seed, ...) provider(fold = seed,
            treatment = t_value, train = train, nuisance_train = nuisance_train, x_cols = x_cols, seed = seed)
        f <- engine$test_learned_witness_dual_msm_compatibility(d, config$Lambda, x_cols = paste0("X",
            1:5), folds = config$folds, primitive_learner = "super_learner", gamma_learner = "rkhs_cv",
            super_learner_candidates = config$primitives, super_learner_cv_folds = 3, super_learner_ranger_num_trees = 500,
            gamma_cv_selection_rule = "minimum", gamma_cv_folds = 3, gamma_cv_repeats = 1, gamma_cv_bandwidth_multipliers = c(0.25,
                0.5, 1, 2, 4), gamma_cv_lambda_grid = default_gamma_lambda_grid, gamma_nystrom_rank = NULL,
            allocation_spline_df = config$allocation_df, bootstrap_replications = B, seed = seed)
        if (!is.null(last_folds))
            stopifnot(identical(last_folds, f$fold_id))
        last_folds <<- f$fold_id
        f$settings$gamma_learner <- method
        f$settings$gamma_provider <- "Frozen x5_adaptive_provider; actual CV selections retained separately"
        selections <- lapply(as.list(attr(provider, "gamma_cache")), function(z) {
            f <- z$fit
            f[intersect(names(f), c("selected_lambda", "selected_level", "selected_bandwidth_multiplier",
                "base_bandwidth", "selected_bandwidth", "feature_count", "rkhs_backend", "cv_table",
                "cv_resample_table", "training_negative_rate", "normalization_constant", "clip_rate",
                "gradient"))]
        })
        list(summary = f$summary, gamma_diagnostics = f$gamma_diagnostics, gamma_selection = selections,
            probability_diagnostics = f$probability_diagnostics, settings = f$settings, witness_diagnostics = f$witness_diagnostics,
            fold_sizes = tabulate(f$fold_id, config$folds))
    }
    list(run = run, engine = engine)
}

