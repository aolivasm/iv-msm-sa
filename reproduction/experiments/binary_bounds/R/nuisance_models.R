# Main-text reproduction dependencies; published formulas retained.
binary_arm_nuisance_from_predictions <- function(predictions, t_value) {
    if (t_value == 1) {
        list(r0 = predictions$pi0, r1 = predictions$pi1, a0 = predictions$a_t1_z0, a1 = predictions$a_t1_z1)
    }
    else {
        list(r0 = 1 - predictions$pi0, r1 = 1 - predictions$pi1, a0 = predictions$a_t0_z0, a1 = predictions$a_t0_z1)
    }
}

binary_default_super_learner_candidates <- function() {
    c("glm_main", "glm_interactions", "glmnet_interactions", "ranger", "xgboost")
}

binary_fit_regression_super_learner <- function(data, response, x_cols, candidates = binary_default_super_learner_candidates(),
    cv_folds = 3L, ranger_num_trees = 500L, hal_max_degree = 2L, seed = 1L) {
    if (!requireNamespace("SuperLearner", quietly = TRUE)) {
        stop("The binary primitive Super Learner requires 'SuperLearner'.")
    }
    candidate_map <- binary_super_learner_candidate_map()
    candidates <- unique(as.character(candidates))
    invalid <- setdiff(candidates, names(candidate_map))
    if (!length(candidates) || length(invalid)) {
        stop("Unknown binary Super Learner candidate(s): ", paste(invalid, collapse = ", "), ".")
    }
    package_by_candidate <- c(glmnet = "glmnet", glmnet_interactions = "glmnet", gam = "mgcv", gam_interactions = "mgcv",
        ranger = "ranger", xgboost = "xgboost", hal = "hal9001")
    missing_packages <- unique(vapply(candidates, function(candidate) {
        package <- unname(package_by_candidate[candidate])
        if (is.na(package) || requireNamespace(package, quietly = TRUE)) {
            ""
        }
        else {
            package
        }
    }, character(1)))
    missing_packages <- missing_packages[nzchar(missing_packages)]
    if (length(missing_packages)) {
        stop("The requested binary Super Learner library is incomplete; install: ", paste(missing_packages,
            collapse = ", "), ".")
    }
    outcome <- as.numeric(data[[response]])
    binary_check_indicator(outcome, response)
    if (length(unique(outcome)) == 1L) {
        probability <- (sum(outcome) + 0.5)/(length(outcome) + 1)
        return(list(predict = function(newdata) rep(probability, nrow(newdata)), converged = TRUE, rank = 1L,
            constant = TRUE, learner = "super_learner", selected_learner = "constant", cv_mean_risk = NA_real_,
            cv_standard_error = NA_real_, cv_table = data.frame(), ensemble_weights = setNames(1, "constant")))
    }
    cv_folds <- as.integer(cv_folds)
    minimum_class_size <- min(tabulate(outcome + 1L, nbins = 2L))
    cv_folds <- min(cv_folds, minimum_class_size)
    if (cv_folds < 2L)
        stop("Super Learner requires both outcome classes.")
    design <- binary_super_learner_numeric_design(data, x_cols)
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
    All <- SuperLearner::All
    SL.glm <- SuperLearner::SL.glm
    SL.glm.interaction <- SuperLearner::SL.glm.interaction
    SL.glmnet <- SuperLearner::SL.glmnet
    base_binary_ranger <- SL.binary_ranger
    base_binary_xgboost <- SL.binary_xgboost
    base_binary_hal <- SL.binary_hal
    SL.binary_ranger <- function(Y, X, newX, family, obsWeights, id, ...) {
        base_binary_ranger(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, id = id,
            num_trees = ranger_num_trees)
    }
    SL.binary_hal <- function(Y, X, newX, family, obsWeights, id, ...) {
        base_binary_hal(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, id = id,
            max_degree = hal_max_degree)
    }
    SL.binary_xgboost <- function(Y, X, newX, family, obsWeights, id, ...) {
        base_binary_xgboost(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, id = id)
    }
    fit_environment <- environment()
    fitted <- suppressWarnings(SuperLearner::SuperLearner(Y = outcome, X = as.data.frame(design$training),
        family = stats::binomial(), SL.library = unname(candidate_map[candidates]), method = binary_super_learner_method(),
        verbose = FALSE, control = list(saveFitLibrary = TRUE, trimLogit = 1e-05), cvControl = list(V = cv_folds,
            stratifyCV = TRUE, shuffle = TRUE), env = fit_environment))
    coefficients <- as.numeric(fitted$coef)
    names(coefficients) <- candidates
    optimizer_status <- fitted$method$optimizer$status
    converged <- is.null(optimizer_status) || optimizer_status %in% 1:4
    cv_table <- data.frame(candidate = candidates, wrapper = unname(candidate_map[candidates]), mean_risk = as.numeric(fitted$cvRisk),
        coefficient = coefficients, selected = coefficients > 1e-08, cv_available = TRUE, stringsAsFactors = FALSE)
    predict_probability <- function(newdata) {
        new_design <- as.data.frame(design$transform(newdata))
        prediction <- stats::predict(fitted, newdata = new_design, onlySL = TRUE)
        as.numeric(prediction$pred)
    }
    list(predict = predict_probability, converged = converged, rank = as.integer(sum(coefficients > 1e-08)),
        constant = FALSE, learner = "super_learner", selected_learner = "convex_stack", cv_mean_risk = binary_cross_entropy_loss(outcome,
            as.numeric(fitted$Z %*% coefficients)), cv_standard_error = NA_real_, cv_table = cv_table,
        fold_risk = data.frame(), ensemble_weights = coefficients, super_learner_fit = fitted, library_configuration = list(candidates_requested = candidates,
            candidates_available = candidates, cv_folds = cv_folds, meta_learner = "simplex_probability_mixture_bernoulli_log_loss",
            ranger_num_trees = as.integer(ranger_num_trees), xgboost_num_trees = 250L, xgboost_max_depth = 3L,
            xgboost_shrinkage = 0.05, xgboost_minimum_child_weight = 10L, xgboost_threads = 1L, hal_max_degree = as.integer(hal_max_degree),
            hal_num_knots_rule = "max(5, ceiling(n_training^(1/3)))"))
}

binary_recycle_arm_nuisance <- function(nuisance) {
    required <- c("r0", "r1", "a0", "a1")
    if (!all(required %in% names(nuisance))) {
        stop("Arm nuisance must contain r0, r1, a0, and a1.")
    }
    target_length <- max(vapply(nuisance[required], length, integer(1)))
    output <- lapply(nuisance[required], rep_len, length.out = target_length)
    if (any(!is.finite(unlist(output, use.names = FALSE)))) {
        stop("Arm nuisance contains non-finite values.")
    }
    output
}

binary_super_learner_candidate_map <- function() {
    c(glm_main = "SL.glm", glm_spline = "SL.binary_glm_spline", glm_interactions = "SL.glm.interaction",
        glmnet = "SL.glmnet", glmnet_interactions = "SL.binary_glmnet_interactions", gam = "SL.binary_mgcv_additive",
        gam_interactions = "SL.binary_mgcv_interactions", ranger = "SL.binary_ranger", xgboost = "SL.binary_xgboost",
        hal = "SL.binary_hal")
}

binary_super_learner_method <- function() {
    if (!requireNamespace("nloptr", quietly = TRUE)) {
        stop("The convex log-loss Super Learner requires package 'nloptr'.")
    }
    compute_coef <- function(Z, Y, libraryNames, verbose, obsWeights, control, errorsInLibrary = NULL,
        ...) {
        if (is.null(dim(Z)))
            Z <- matrix(Z, ncol = 1L)
        Z <- as.matrix(Z)
        probability <- Z
        probability[] <- pmin(1 - 1e-08, pmax(1e-08, as.numeric(Z)))
        weights <- as.numeric(obsWeights)
        weights <- weights/sum(weights)
        individual_risk <- vapply(seq_len(ncol(probability)), function(index) {
            binary_cross_entropy_loss(Y, probability[, index])
        }, numeric(1))
        names(individual_risk) <- libraryNames
        unavailable <- apply(!is.finite(probability), 2L, any)
        if (!is.null(errorsInLibrary))
            unavailable <- unavailable | errorsInLibrary
        available <- which(!unavailable)
        if (!length(available))
            stop("Every Super Learner candidate failed.")
        local_probability <- probability[, available, drop = FALSE]
        objective <- function(beta) {
            combined <- pmin(1 - 1e-08, pmax(1e-08, as.numeric(local_probability %*% beta)))
            loss <- -Y * log(combined) - (1 - Y) * log1p(-combined)
            gradient_factor <- weights * ((combined - Y)/(combined * (1 - combined)))
            list(objective = sum(weights * loss), gradient = as.numeric(crossprod(local_probability,
                gradient_factor)))
        }
        start <- rep(1/length(available), length(available))
        optimizer <- nloptr::nloptr(x0 = start, eval_f = objective, lb = rep(0, length(available)), ub = rep(1,
            length(available)), eval_g_eq = function(beta) sum(beta) - 1, eval_jac_g_eq = function(beta) rep(1,
            length(beta)), opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-09, maxeval = 1000L,
            print_level = 0L))
        coefficient <- rep(0, ncol(probability))
        coefficient[available] <- optimizer$solution
        coefficient[coefficient < 1e-08] <- 0
        coefficient <- coefficient/sum(coefficient)
        names(coefficient) <- libraryNames
        if (verbose) {
            message("Convex log-loss optimizer status: ", optimizer$status, " (", optimizer$message,
                ")")
        }
        list(cvRisk = individual_risk, coef = coefficient, optimizer = optimizer)
    }
    compute_prediction <- function(predY, coef, control, ...) {
        if (is.null(dim(predY)))
            predY <- matrix(predY, ncol = 1L)
        prediction <- as.numeric(as.matrix(predY) %*% as.numeric(coef))
        pmin(1 - 1e-08, pmax(1e-08, prediction))
    }
    list(require = "nloptr", computeCoef = compute_coef, computePred = compute_prediction)
}

binary_super_learner_numeric_design <- function(data, x_cols) {
    predictor_data <- data[, x_cols, drop = FALSE]
    formula <- stats::reformulate(x_cols, response = NULL, intercept = FALSE)
    terms <- stats::terms(formula, data = predictor_data)
    frame <- stats::model.frame(terms, data = predictor_data, na.action = stats::na.fail, drop.unused.levels = TRUE)
    xlevels <- stats::.getXlevels(terms, frame)
    contrasts <- lapply(frame[vapply(frame, is.factor, logical(1))], contrasts)
    design <- stats::model.matrix(terms, data = frame, contrasts.arg = contrasts)
    storage.mode(design) <- "double"
    if (!ncol(design) || any(!is.finite(design))) {
        stop("The Super Learner predictor design is empty or non-finite.")
    }
    training_columns <- colnames(design)
    transform <- function(newdata) {
        new_frame <- stats::model.frame(terms, data = newdata[, x_cols, drop = FALSE], na.action = stats::na.fail,
            xlev = xlevels)
        new_design <- stats::model.matrix(terms, data = new_frame, contrasts.arg = contrasts)
        missing <- setdiff(training_columns, colnames(new_design))
        if (length(missing)) {
            new_design <- cbind(new_design, matrix(0, nrow = nrow(new_design), ncol = length(missing),
                dimnames = list(NULL, missing)))
        }
        extra <- setdiff(colnames(new_design), training_columns)
        if (length(extra)) {
            stop("New Super Learner data generated unknown design columns: ", paste(extra, collapse = ", "),
                ".")
        }
        new_design <- new_design[, training_columns, drop = FALSE]
        storage.mode(new_design) <- "double"
        if (any(!is.finite(new_design))) {
            stop("New Super Learner predictor values are non-finite.")
        }
        new_design
    }
    list(training = design, transform = transform, columns = training_columns)
}

evaluate_binary_bounds_from_nuisances <- function(data, nuisance_fit, Lambda0 = 2, Lambda1 = Lambda0,
    epsilon = 0.01, confidence_level = 0.95, compatibility_action = c("project", "none"), compatibility_tolerance = 1e-12) {
    compatibility_action <- match.arg(compatibility_action)
    predictions <- nuisance_fit$predictions
    if (!is.finite(confidence_level) || confidence_level <= 0 || confidence_level >= 1) {
        stop("confidence_level must lie in (0,1).")
    }
    hard <- smooth <- signals <- unstabilized_hard <- unstabilized_smooth <- stabilization <- list()
    iv_hard <- iv_smooth <- iv_signals <- list()
    for (t_value in c(0, 1)) {
        key <- as.character(t_value)
        raw_arm <- binary_arm_nuisance_from_predictions(predictions, t_value)
        iv_hard[[key]] <- binary_iv_arm_endpoint(raw_arm, epsilon = 0, compute_gradient = FALSE)
        iv_smooth[[key]] <- binary_iv_arm_endpoint(raw_arm, epsilon = epsilon, compute_gradient = TRUE)
        iv_signals[[key]] <- list(lower = binary_arm_one_step_signal(data, predictions, raw_arm, iv_smooth[[key]],
            t_value, "lower"), upper = binary_arm_one_step_signal(data, predictions, raw_arm, iv_smooth[[key]],
            t_value, "upper"))
        unstabilized_hard[[key]] <- binary_arm_endpoint(raw_arm, Lambda0, Lambda1, epsilon = 0, compute_gradient = FALSE)
        unstabilized_smooth[[key]] <- binary_arm_endpoint(raw_arm, Lambda0, Lambda1, epsilon = epsilon,
            compute_gradient = FALSE)
        stabilization[[key]] <- binary_project_arm_outcome_probability(raw_arm, Lambda0 = Lambda0, Lambda1 = Lambda1,
            action = compatibility_action, tolerance = compatibility_tolerance)
        arm <- stabilization[[key]]$nuisance
        hard[[key]] <- binary_arm_endpoint(arm, Lambda0, Lambda1, epsilon = 0, compute_gradient = FALSE)
        smooth[[key]] <- binary_arm_endpoint(arm, Lambda0, Lambda1, epsilon = epsilon, compute_gradient = TRUE)
        signals[[key]] <- list(lower = binary_arm_one_step_signal(data, predictions, arm, smooth[[key]],
            t_value, "lower"), upper = binary_arm_one_step_signal(data, predictions, arm, smooth[[key]],
            t_value, "upper"))
    }
    endpoint_vector <- function(container, use_signal = FALSE) {
        get_value <- function(treatment, direction) {
            value <- container[[as.character(treatment)]][[direction]]
            if (use_signal)
                mean(value)
            else mean(value)
        }
        mu1_lower <- get_value(1, "lower")
        mu1_upper <- get_value(1, "upper")
        mu0_lower <- get_value(0, "lower")
        mu0_upper <- get_value(0, "upper")
        c(mu1_lower = mu1_lower, mu1_upper = mu1_upper, mu0_lower = mu0_lower, mu0_upper = mu0_upper,
            ate_lower = mu1_lower - mu0_upper, ate_upper = mu1_upper - mu0_lower)
    }
    hard_values <- endpoint_vector(hard)
    smooth_values <- endpoint_vector(smooth)
    one_step_values <- endpoint_vector(signals, use_signal = TRUE)
    iv_hard_values <- endpoint_vector(iv_hard)
    iv_smooth_values <- endpoint_vector(iv_smooth)
    iv_one_step_values <- endpoint_vector(iv_signals, use_signal = TRUE)
    ate_lower_signal <- signals[["1"]]$lower - signals[["0"]]$upper
    ate_upper_signal <- signals[["1"]]$upper - signals[["0"]]$lower
    standard_errors <- c(mu1_lower = stats::sd(signals[["1"]]$lower)/sqrt(nrow(data)), mu1_upper = stats::sd(signals[["1"]]$upper)/sqrt(nrow(data)),
        mu0_lower = stats::sd(signals[["0"]]$lower)/sqrt(nrow(data)), mu0_upper = stats::sd(signals[["0"]]$upper)/sqrt(nrow(data)),
        ate_lower = stats::sd(ate_lower_signal)/sqrt(nrow(data)), ate_upper = stats::sd(ate_upper_signal)/sqrt(nrow(data)))
    iv_ate_lower_signal <- iv_signals[["1"]]$lower - iv_signals[["0"]]$upper
    iv_ate_upper_signal <- iv_signals[["1"]]$upper - iv_signals[["0"]]$lower
    iv_standard_errors <- c(mu1_lower = stats::sd(iv_signals[["1"]]$lower)/sqrt(nrow(data)), mu1_upper = stats::sd(iv_signals[["1"]]$upper)/sqrt(nrow(data)),
        mu0_lower = stats::sd(iv_signals[["0"]]$lower)/sqrt(nrow(data)), mu0_upper = stats::sd(iv_signals[["0"]]$upper)/sqrt(nrow(data)),
        ate_lower = stats::sd(iv_ate_lower_signal)/sqrt(nrow(data)), ate_upper = stats::sd(iv_ate_upper_signal)/sqrt(nrow(data)))
    long_rows <- function(method, values, se = rep(NA_real_, length(values))) {
        data.frame(method = method, endpoint = names(values), estimate = unname(values), standard_error = unname(se[names(values)]),
            stringsAsFactors = FALSE)
    }
    estimates <- rbind(long_rows("dual_hard_plugin", hard_values), long_rows("dual_smooth_plugin", smooth_values),
        long_rows("dual_smooth_one_step", one_step_values, standard_errors), long_rows("iv_only_smooth_plugin",
            iv_smooth_values), long_rows("iv_only_smooth_one_step", iv_one_step_values, iv_standard_errors))
    comparator_wide <- binary_comparator_plugin_values(predictions, Lambda0)
    comparator_long <- do.call(rbind, lapply(seq_len(nrow(comparator_wide)), function(i) {
        values <- unlist(comparator_wide[i, setdiff(names(comparator_wide), "method")])
        names(values) <- setdiff(names(comparator_wide), "method")
        long_rows(comparator_wide$method[i], values)
    }))
    estimates <- rbind(estimates, comparator_long)
    critical_two_sided <- stats::qnorm(1 - (1 - confidence_level)/2)
    one_step_rows <- estimates$method == "dual_smooth_one_step"
    intervals <- estimates[one_step_rows, , drop = FALSE]
    intervals$confidence_level <- confidence_level
    intervals$ci_lower <- intervals$estimate - critical_two_sided * intervals$standard_error
    intervals$ci_upper <- intervals$estimate + critical_two_sided * intervals$standard_error
    ate_confidence_set <- data.frame(confidence_level = confidence_level, lower = one_step_values["ate_lower"] -
        critical_two_sided * standard_errors["ate_lower"], upper = one_step_values["ate_upper"] + critical_two_sided *
        standard_errors["ate_upper"], stringsAsFactors = FALSE)
    iv_one_step_rows <- estimates$method == "iv_only_smooth_one_step"
    iv_endpoint_intervals <- estimates[iv_one_step_rows, , drop = FALSE]
    iv_endpoint_intervals$confidence_level <- confidence_level
    iv_endpoint_intervals$ci_lower <- iv_endpoint_intervals$estimate - critical_two_sided * iv_endpoint_intervals$standard_error
    iv_endpoint_intervals$ci_upper <- iv_endpoint_intervals$estimate + critical_two_sided * iv_endpoint_intervals$standard_error
    iv_ate_confidence_set <- data.frame(confidence_level = confidence_level, lower = iv_one_step_values["ate_lower"] -
        critical_two_sided * iv_standard_errors["ate_lower"], upper = iv_one_step_values["ate_upper"] +
        critical_two_sided * iv_standard_errors["ate_upper"], stringsAsFactors = FALSE)
    diagnostics <- do.call(rbind, lapply(c(0, 1), function(t_value) {
        key <- as.character(t_value)
        raw_hard <- unstabilized_hard[[key]]
        raw_smooth <- unstabilized_smooth[[key]]
        post_hard <- hard[[key]]
        post_smooth <- smooth[[key]]
        hard_gamma_magnitude <- pmax(raw_hard$compatible_gamma_lower - raw_hard$gamma_one, raw_hard$gamma_one -
            raw_hard$compatible_gamma_upper, raw_hard$compatible_gamma_lower - raw_hard$gamma_zero, raw_hard$gamma_zero -
            raw_hard$compatible_gamma_upper, 0)
        hard_envelope_magnitude <- pmax(raw_hard$L_one - raw_hard$U_one, raw_hard$L_zero - raw_hard$U_zero,
            0)
        hard_theta_magnitude <- pmax(raw_hard$theta_lower - raw_hard$theta_upper, 0)
        smooth_theta_magnitude <- pmax(raw_smooth$theta_lower - raw_smooth$theta_upper, 0)
        diagnostic_tolerance <- 1e-08
        post_gamma_violation <- post_hard$gamma_violation
        post_envelope_crossing <- post_hard$envelope_crossing
        post_hard_theta_crossing <- post_hard$theta_crossing
        post_smooth_theta_crossing <- post_smooth$theta_crossing
        data.frame(treatment = t_value, compatibility_action = compatibility_action, compatibility_projection_required_rate = mean(stabilization[[key]]$required),
            compatibility_projection_rate = if (compatibility_action == "project") {
                mean(stabilization[[key]]$required)
            }
            else 0, mean_compatibility_projection_magnitude = mean(stabilization[[key]]$magnitude), maximum_compatibility_projection_magnitude = max(stabilization[[key]]$magnitude),
            hard_gamma_violation_rate = mean(raw_hard$gamma_violation), hard_envelope_crossing_rate = mean(raw_hard$envelope_crossing),
            hard_theta_crossing_rate = mean(raw_hard$theta_crossing), smooth_theta_crossing_rate = mean(raw_smooth$theta_crossing),
            post_hard_gamma_violation_rate = mean(post_gamma_violation), post_hard_envelope_crossing_rate = mean(post_envelope_crossing),
            post_hard_theta_crossing_rate = mean(post_hard_theta_crossing), post_smooth_theta_crossing_rate = mean(post_smooth_theta_crossing),
            diagnostic_tolerance = diagnostic_tolerance, hard_gamma_violation_rate_above_tolerance = mean(hard_gamma_magnitude >
                diagnostic_tolerance), hard_envelope_crossing_rate_above_tolerance = mean(hard_envelope_magnitude >
                diagnostic_tolerance), hard_theta_crossing_rate_above_tolerance = mean(hard_theta_magnitude >
                diagnostic_tolerance), smooth_theta_crossing_rate_above_tolerance = mean(smooth_theta_magnitude >
                diagnostic_tolerance), mean_hard_gamma_violation_magnitude = mean(hard_gamma_magnitude),
            maximum_hard_gamma_violation_magnitude = max(hard_gamma_magnitude), mean_hard_envelope_crossing_magnitude = mean(hard_envelope_magnitude),
            maximum_hard_envelope_crossing_magnitude = max(hard_envelope_magnitude), mean_hard_theta_crossing_magnitude = mean(hard_theta_magnitude),
            maximum_hard_theta_crossing_magnitude = max(hard_theta_magnitude), mean_smooth_theta_crossing_magnitude = mean(smooth_theta_magnitude),
            maximum_smooth_theta_crossing_magnitude = max(smooth_theta_magnitude), stringsAsFactors = FALSE)
    }))
    list(estimates = estimates, endpoint_intervals = intervals, ate_confidence_set = ate_confidence_set,
        iv_endpoint_intervals = iv_endpoint_intervals, iv_ate_confidence_set = iv_ate_confidence_set,
        diagnostics = diagnostics, hard = hard, smooth = smooth, unstabilized_hard = unstabilized_hard,
        unstabilized_smooth = unstabilized_smooth, stabilization = stabilization, signals = signals,
        iv_hard = iv_hard, iv_smooth = iv_smooth, iv_signals = iv_signals, Lambda0 = Lambda0, Lambda1 = Lambda1,
        epsilon = epsilon, compatibility_action = compatibility_action, compatibility_tolerance = compatibility_tolerance)
}

fit_binary_cross_fitted_direct_rkhs_nuisances <- function(data, x_cols, K = 5L, fold_id = NULL, seed = 1L,
    overlap_epsilon = 0.05, outcome_probability_epsilon = 1e-04, spline_df = 5L, gam_min_basis = 5L,
    gam_max_basis = 20L, gam_growth_exponent = 1/3, primitive_learner = c("super_learner", "library_cv",
        "gam", "ranger", "glm_spline", "glmnet"), library_candidates = c("glm_spline", "gam", "glmnet",
        "ranger"), library_cv_folds = 3L, library_selection_rule = c("one_se", "minimum"), library_ranger_num_trees = 500L,
    super_learner_candidates = binary_default_super_learner_candidates(), super_learner_cv_folds = 3L,
    super_learner_ranger_num_trees = 500L, super_learner_hal_max_degree = 2L, rkhs_lambda_by_t = NULL,
    rkhs_lambda_grid = default_gamma_lambda_grid, rkhs_bandwidth_multipliers = 1/sqrt(2), rkhs_cv_folds = 3L,
    rkhs_cv_repeats = 2L, rkhs_selection_rule = c("minimum", "one_se"), rkhs_bandwidth_multiplier = 1/sqrt(2),
    rkhs_gamma_floor = 1e-06, rkhs_gamma_cap = 100, rkhs_backend = c("nystrom", "exact"), rkhs_nystrom_rank = NULL,
    rkhs_kernel = c("gaussian", "sobolev_matern32", "sobolev_matern52"), rkhs_bandwidth_sample = 1000L) {
    binary_rkhs_gamma_requirements()
    primitive_learner <- match.arg(primitive_learner)
    library_selection_rule <- match.arg(library_selection_rule)
    rkhs_selection_rule <- match.arg(rkhs_selection_rule)
    rkhs_backend <- match.arg(rkhs_backend)
    rkhs_kernel <- match.arg(rkhs_kernel)
    gamma_tuning_mode <- if (is.null(rkhs_lambda_by_t))
        "joint_cv"
    else "fixed"
    if (gamma_tuning_mode == "joint_cv" && (length(rkhs_cv_folds) != 1L || !is.finite(rkhs_cv_folds) ||
        rkhs_cv_folds < 2L || rkhs_cv_folds != as.integer(rkhs_cv_folds))) {
        stop("rkhs_cv_folds must be one integer of at least two.")
    }
    if (gamma_tuning_mode == "joint_cv" && (length(rkhs_cv_repeats) != 1L || !is.finite(rkhs_cv_repeats) ||
        rkhs_cv_repeats < 1L || rkhs_cv_repeats != as.integer(rkhs_cv_repeats))) {
        stop("rkhs_cv_repeats must be one positive integer.")
    }
    required <- c("Y", "T_obs", "Z", x_cols)
    missing <- setdiff(required, names(data))
    if (length(missing))
        stop("Missing columns: ", paste(missing, collapse = ", "))
    if (is.null(fold_id)) {
        fold_id <- binary_make_tz_folds(data, K = K, seed = seed)
    }
    if (length(fold_id) != nrow(data))
        stop("fold_id has the wrong length.")
    folds <- sort(unique(fold_id))
    if (length(folds) < 2L)
        stop("At least two outer folds are required.")
    prediction_names <- c("e", "pi0", "pi1", "a_t0_z0", "a_t0_z1", "a_t1_z0", "a_t1_z1")
    predictions <- as.data.frame(setNames(replicate(length(prediction_names), rep(NA_real_, nrow(data)),
        simplify = FALSE), prediction_names))
    fit_diagnostics <- gamma_diagnostics <- reconstruction_diagnostics <- list()
    cv_tables <- gamma_cv_tables <- gamma_cv_resample_tables <- list()
    fit_counter <- gamma_counter <- reconstruction_counter <- 0L
    fit_component <- function(training, response, local_seed) {
        binary_rkhs_gamma_fit_component(training = training, response = response, x_cols = x_cols, seed = local_seed,
            spline_df = spline_df, gam_min_basis = gam_min_basis, gam_max_basis = gam_max_basis, gam_growth_exponent = gam_growth_exponent,
            library_candidates = library_candidates, library_cv_folds = library_cv_folds, library_selection_rule = library_selection_rule,
            library_ranger_num_trees = library_ranger_num_trees, learner = primitive_learner, super_learner_candidates = super_learner_candidates,
            super_learner_cv_folds = super_learner_cv_folds, super_learner_ranger_num_trees = super_learner_ranger_num_trees,
            super_learner_hal_max_degree = super_learner_hal_max_degree)
    }
    clip_design <- function(value) {
        binary_clip_probability(value, overlap_epsilon)
    }
    clip_outcome <- function(value) {
        binary_clip_probability(value, outcome_probability_epsilon)
    }
    for (fold in folds) {
        training <- data[fold_id != fold, , drop = FALSE]
        validation <- data[fold_id == fold, , drop = FALSE]
        validation_index <- which(fold_id == fold)
        fits <- list(e = fit_component(training, "Z", seed + 1000L * fold + 1L), pi0 = fit_component(training[training$Z ==
            0, , drop = FALSE], "T_obs", seed + 1000L * fold + 2L), pi1 = fit_component(training[training$Z ==
            1, , drop = FALSE], "T_obs", seed + 1000L * fold + 3L))
        for (treatment in 0:1) {
            baseline_name <- sprintf("a_t%d_z0", treatment)
            cell <- training$T_obs == treatment & training$Z == 0
            fits[[baseline_name]] <- fit_component(training[cell, , drop = FALSE], "Y", seed + 1000L *
                fold + 100L * treatment + 4L)
        }
        raw_training <- list(e = 1 - fits$e$predict(training), pi0 = fits$pi0$predict(training), pi1 = fits$pi1$predict(training))
        raw_validation <- list(e = 1 - fits$e$predict(validation), pi0 = fits$pi0$predict(validation),
            pi1 = fits$pi1$predict(validation))
        for (treatment in 0:1) {
            baseline_name <- sprintf("a_t%d_z0", treatment)
            raw_training[[baseline_name]] <- fits[[baseline_name]]$predict(training)
            raw_validation[[baseline_name]] <- fits[[baseline_name]]$predict(validation)
        }
        estimated_training <- data.frame(e = clip_design(raw_training$e), pi0 = clip_design(raw_training$pi0),
            pi1 = clip_design(raw_training$pi1), stringsAsFactors = FALSE)
        estimated_validation <- data.frame(e = clip_design(raw_validation$e), pi0 = clip_design(raw_validation$pi0),
            pi1 = clip_design(raw_validation$pi1), stringsAsFactors = FALSE)
        for (treatment in 0:1) {
            baseline_name <- sprintf("a_t%d_z0", treatment)
            estimated_training[[baseline_name]] <- clip_outcome(raw_training[[baseline_name]])
            estimated_validation[[baseline_name]] <- clip_outcome(raw_validation[[baseline_name]])
        }
        predictions[validation_index, c("e", "pi0", "pi1")] <- estimated_validation[c("e", "pi0", "pi1")]
        for (treatment in 0:1) {
            baseline_name <- sprintf("a_t%d_z0", treatment)
            comparison_name <- sprintf("a_t%d_z1", treatment)
            n_training_t <- sum(training$T_obs == treatment)
            gamma_seed <- seed + 100000L * fold + 1000L * treatment
            gamma_fit <- if (gamma_tuning_mode == "joint_cv") {
                fit_gamma_rkhs(train = training, nuisance_train = estimated_training, t_value = treatment,
                  x_cols = x_cols, lambda_grid = rkhs_lambda_grid, bandwidth_multipliers = rkhs_bandwidth_multipliers,
                  cv_folds = rkhs_cv_folds, cv_repeats = rkhs_cv_repeats, selection_rule = rkhs_selection_rule,
                  gamma_floor = rkhs_gamma_floor, gamma_cap = rkhs_gamma_cap, bandwidth_sample = rkhs_bandwidth_sample,
                  backend = rkhs_backend, nystrom_rank = rkhs_nystrom_rank, kernel_type = rkhs_kernel,
                  seed = gamma_seed)
            }
            else {
                selected_lambda <- binary_resolve_rkhs_lambda(rkhs_lambda_by_t = rkhs_lambda_by_t, treatment = treatment,
                  n_training = n_training_t)
                fit_gamma_rkhs_fixed(train = training, nuisance_train = estimated_training, t_value = treatment,
                  x_cols = x_cols, lambda = selected_lambda, bandwidth_multiplier = rkhs_bandwidth_multiplier,
                  gamma_floor = rkhs_gamma_floor, gamma_cap = rkhs_gamma_cap, bandwidth_sample = rkhs_bandwidth_sample,
                  backend = rkhs_backend, nystrom_rank = rkhs_nystrom_rank, kernel_type = rkhs_kernel,
                  seed = gamma_seed)
            }
            if (!is.null(gamma_fit$cv_table) && nrow(gamma_fit$cv_table)) {
                table <- gamma_fit$cv_table
                table$outer_fold <- fold
                table$treatment <- treatment
                table$n_training_t <- n_training_t
                table$lambda_sqrt_n_constant <- table$lambda * sqrt(n_training_t)
                table$selection_rule <- rkhs_selection_rule
                gamma_cv_tables[[length(gamma_cv_tables) + 1L]] <- table
            }
            if (!is.null(gamma_fit$cv_resample_table) && nrow(gamma_fit$cv_resample_table)) {
                resample_table <- gamma_fit$cv_resample_table
                resample_table$outer_fold <- fold
                resample_table$treatment <- treatment
                resample_table$n_training_t <- n_training_t
                resample_table$lambda_sqrt_n_constant <- resample_table$lambda * sqrt(n_training_t)
                resample_table$selection_rule <- rkhs_selection_rule
                gamma_cv_resample_tables[[length(gamma_cv_resample_tables) + 1L]] <- resample_table
            }
            new_zero <- new_one <- validation
            new_zero$Y <- 0
            new_one$Y <- 1
            prediction_zero <- gamma_fit$predict(new_zero)
            prediction_one <- gamma_fit$predict(new_one)
            reconstructed <- binary_rkhs_gamma_reconstruct(a0 = estimated_validation[[baseline_name]],
                gamma_zero = prediction_zero$normalized, gamma_one = prediction_one$normalized, probability_epsilon = outcome_probability_epsilon,
                gamma_floor = rkhs_gamma_floor)
            predictions[[baseline_name]][validation_index] <- estimated_validation[[baseline_name]]
            predictions[[comparison_name]][validation_index] <- reconstructed$a1
            gamma_counter <- gamma_counter + 1L
            gamma_diagnostics[[gamma_counter]] <- data.frame(fold = fold, treatment = treatment, n_training_t = n_training_t,
                tuning_mode = gamma_tuning_mode, cv_folds = if (gamma_tuning_mode == "joint_cv") {
                  as.integer(rkhs_cv_folds)
                }
                else NA_integer_, cv_repeats = if (gamma_tuning_mode == "joint_cv") {
                  as.integer(rkhs_cv_repeats)
                }
                else NA_integer_, selected_lambda = gamma_fit$selected_lambda, selected_lambda_sqrt_n_constant = gamma_fit$selected_lambda *
                  sqrt(n_training_t), selected_bandwidth_multiplier = gamma_fit$selected_bandwidth_multiplier,
                base_bandwidth = gamma_fit$base_bandwidth, selected_bandwidth = gamma_fit$selected_bandwidth,
                normalization_constant = gamma_fit$normalization_constant, training_negative_rate = gamma_fit$training_negative_rate,
                rkhs_backend = gamma_fit$rkhs_backend, rkhs_kernel = gamma_fit$rkhs_kernel, nystrom_rank_requested = gamma_fit$nystrom_rank_requested,
                nystrom_landmark_count = gamma_fit$nystrom_landmark_count, nystrom_effective_rank = gamma_fit$nystrom_effective_rank,
                stringsAsFactors = FALSE)
            reconstruction_counter <- reconstruction_counter + 1L
            reconstruction_diagnostics[[reconstruction_counter]] <- data.frame(fold = fold, treatment = treatment,
                mean_absolute_pre_normalization_error = mean(abs(reconstructed$normalizer - 1)), maximum_absolute_pre_normalization_error = max(abs(reconstructed$normalizer -
                  1)), lower_probability_clipping_rate = reconstructed$lower_probability_clipping_rate,
                upper_probability_clipping_rate = reconstructed$upper_probability_clipping_rate, stringsAsFactors = FALSE)
        }
        component_epsilon <- c(e = overlap_epsilon, pi0 = overlap_epsilon, pi1 = overlap_epsilon, a_t0_z0 = outcome_probability_epsilon,
            a_t1_z0 = outcome_probability_epsilon)
        for (component in names(fits)) {
            fit_counter <- fit_counter + 1L
            fit <- fits[[component]]
            fit_diagnostics[[fit_counter]] <- data.frame(fold = fold, component = component, learner = fit$learner,
                selected_learner = if (is.null(fit$selected_learner)) {
                  fit$learner
                }
                else {
                  fit$selected_learner
                }, converged = isTRUE(fit$converged), rank = if (is.null(fit$rank))
                  NA_integer_
                else as.integer(fit$rank), constant_fit = isTRUE(fit$constant), basis_dimension = if (is.null(fit$basis_dimension)) {
                  NA_integer_
                }
                else {
                  as.integer(fit$basis_dimension)
                }, effective_df = if (is.null(fit$effective_df)) {
                  NA_real_
                }
                else {
                  as.numeric(fit$effective_df)
                }, cv_mean_risk = if (is.null(fit$cv_mean_risk)) {
                  NA_real_
                }
                else {
                  as.numeric(fit$cv_mean_risk)
                }, truncation_epsilon = unname(component_epsilon[component]), lower_truncation_rate = mean(raw_validation[[component]] <
                  component_epsilon[[component]]), upper_truncation_rate = mean(raw_validation[[component]] >
                  1 - component_epsilon[[component]]), stringsAsFactors = FALSE)
            if (!is.null(fit$cv_table) && nrow(fit$cv_table)) {
                table <- fit$cv_table
                table$outer_fold <- fold
                table$component <- component
                cv_tables[[length(cv_tables) + 1L]] <- table
            }
        }
    }
    if (any(!is.finite(as.matrix(predictions)))) {
        stop("At least one direct-RKHS binary nuisance prediction is non-finite.")
    }
    list(predictions = predictions, fold_id = fold_id, learner = "direct_rkhs_gamma", diagnostics = do.call(rbind,
        fit_diagnostics), cv_table = if (length(cv_tables)) do.call(rbind, cv_tables) else data.frame(),
        gamma_cv_table = if (length(gamma_cv_tables)) {
            do.call(rbind, gamma_cv_tables)
        } else {
            data.frame()
        }, gamma_cv_resample_table = if (length(gamma_cv_resample_tables)) {
            do.call(rbind, gamma_cv_resample_tables)
        } else {
            data.frame()
        }, gamma_diagnostics = do.call(rbind, gamma_diagnostics), reconstruction_diagnostics = do.call(rbind,
            reconstruction_diagnostics), configuration = list(K = length(folds), primitive_learner = primitive_learner,
            overlap_epsilon = overlap_epsilon, outcome_probability_epsilon = outcome_probability_epsilon,
            library_candidates = library_candidates, library_cv_folds = as.integer(library_cv_folds),
            library_selection_rule = library_selection_rule, super_learner_candidates = super_learner_candidates,
            super_learner_cv_folds = as.integer(super_learner_cv_folds), super_learner_ranger_num_trees = as.integer(super_learner_ranger_num_trees),
            super_learner_hal_max_degree = as.integer(super_learner_hal_max_degree), rkhs_tuning_mode = gamma_tuning_mode,
            rkhs_lambda_specification = if (gamma_tuning_mode == "joint_cv") {
                if (identical(rkhs_lambda_grid, default_gamma_lambda_grid)) {
                  "C/sqrt(n_t)"
                } else if (is.function(rkhs_lambda_grid)) {
                  "user_function_of_treatment_arm_training_size"
                } else {
                  "user_candidate_grid"
                }
            } else if (is.function(rkhs_lambda_by_t)) {
                "fixed_function_of_treatment_arm_training_size"
            } else {
                rkhs_lambda_by_t
            }, rkhs_lambda_grid_at_n1 = if (gamma_tuning_mode == "joint_cv") {
                resolve_gamma_lambda_grid(rkhs_lambda_grid, n_training = 1L, treatment = 0L)
            } else {
                NULL
            }, rkhs_bandwidth_multipliers = if (gamma_tuning_mode == "joint_cv") {
                rkhs_bandwidth_multipliers
            } else {
                rkhs_bandwidth_multiplier
            }, rkhs_cv_folds = if (gamma_tuning_mode == "joint_cv") {
                as.integer(rkhs_cv_folds)
            } else {
                NA_integer_
            }, rkhs_cv_repeats = if (gamma_tuning_mode == "joint_cv") {
                as.integer(rkhs_cv_repeats)
            } else {
                NA_integer_
            }, rkhs_selection_rule = if (gamma_tuning_mode == "joint_cv") {
                rkhs_selection_rule
            } else {
                NA_character_
            }, rkhs_gamma_floor = rkhs_gamma_floor, rkhs_gamma_cap = rkhs_gamma_cap, rkhs_backend = rkhs_backend,
            rkhs_nystrom_rank = rkhs_nystrom_rank, rkhs_kernel = rkhs_kernel))
}

observed_dgp_nuisances <- function(x, specification) {
    features <- observed_dgp_features(x, specification)
    treatment_intercept <- if (is.null(specification$treatment_intercept))
        -0.15
    else specification$treatment_intercept
    probability_z1 <- stats::plogis(features$instrument)
    pi0 <- stats::plogis(treatment_intercept + features$treatment)
    pi1 <- stats::plogis(treatment_intercept + features$treatment + specification$instrument_log_odds_strength)
    list(e = 1 - probability_z1, probability_z1 = probability_z1, pi0 = pi0, pi1 = pi1, outcome_feature = features$outcome,
        tilt_shape = features$tilt_shape)
}

predict.SL.binary_glm_spline <- function(object, newdata, ...) {
    as.numeric(stats::predict(object$object, newdata = as.data.frame(newdata), type = "response"))
}

predict.SL.binary_glmnet_interactions <- function(object, newdata, ...) {
    design <- binary_pairwise_interaction_design(newdata)
    as.numeric(stats::predict(object$object, newx = design, type = "response", s = "lambda.1se"))
}

predict.SL.binary_hal <- function(object, newdata, ...) {
    as.numeric(stats::predict(object$object, new_data = as.matrix(newdata)))
}

predict.SL.binary_mgcv_additive <- function(object, newdata, ...) {
    as.numeric(stats::predict(object$object, newdata = as.data.frame(newdata), type = "response"))
}

predict.SL.binary_mgcv_interactions <- function(object, newdata, ...) {
    as.numeric(stats::predict(object$object, newdata = as.data.frame(newdata), type = "response"))
}

predict.SL.binary_xgboost <- function(object, newdata, ...) {
    design <- as.matrix(newdata)
    storage.mode(design) <- "double"
    as.numeric(stats::predict(object$object, newdata = design))
}

SL.binary_glm_spline <- function(Y, X, newX, family, obsWeights, spline_df = 3L, ...) {
    training <- as.data.frame(X)
    training$.outcome <- as.numeric(Y)
    formula <- binary_make_spline_formula(".outcome", training, setdiff(names(training), ".outcome"),
        spline_df)
    environment(formula) <- environment()
    fit <- suppressWarnings(stats::glm(formula, data = training, family = family, weights = obsWeights,
        control = stats::glm.control(maxit = 100L)))
    prediction <- as.numeric(stats::predict(fit, newdata = as.data.frame(newX), type = "response"))
    output <- list(pred = prediction, fit = list(object = fit))
    class(output$fit) <- "SL.binary_glm_spline"
    output
}

SL.binary_glmnet_interactions <- function(Y, X, newX, family, obsWeights, id, nfolds = 5L, nlambda = 75L,
    ...) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("SL.binary_glmnet_interactions requires package 'glmnet'.")
    }
    design <- binary_pairwise_interaction_design(X)
    new_design <- binary_pairwise_interaction_design(newX)
    minimum_class_size <- min(tabulate(as.integer(Y) + 1L, nbins = 2L))
    local_folds <- min(as.integer(nfolds), minimum_class_size)
    if (local_folds < 4L) {
        stop("At least four observations per outcome class are required by glmnet.")
    }
    fit <- suppressWarnings(glmnet::cv.glmnet(x = design, y = as.numeric(Y), weights = obsWeights, family = family$family,
        alpha = 1, nfolds = local_folds, nlambda = as.integer(nlambda), type.measure = "deviance", standardize = TRUE,
        intercept = TRUE, parallel = FALSE))
    prediction <- as.numeric(stats::predict(fit, newx = new_design, type = "response", s = "lambda.1se"))
    output <- list(pred = prediction, fit = list(object = fit, use_one_se = TRUE))
    class(output$fit) <- "SL.binary_glmnet_interactions"
    output
}

SL.binary_hal <- function(Y, X, newX, family, obsWeights, id, max_degree = 2L, num_knots = NULL, ...) {
    if (!requireNamespace("hal9001", quietly = TRUE)) {
        stop("SL.binary_hal requires package 'hal9001'.")
    }
    if (is.null(num_knots)) {
        num_knots <- max(5L, ceiling(nrow(X)^(1/3)))
    }
    hal9001::SL.hal9001(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, id = id,
        max_degree = as.integer(max_degree), smoothness_orders = 1, num_knots = as.integer(num_knots))
}

SL.binary_mgcv_additive <- function(Y, X, newX, family, obsWeights, basis_dimension = NULL, ...) {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("SL.binary_mgcv_additive requires package 'mgcv'.")
    }
    training <- as.data.frame(X)
    training$.outcome <- as.numeric(Y)
    predictors <- setdiff(names(training), ".outcome")
    local_basis <- binary_resolve_gam_basis_dimension(training, predictors, basis_dimension = basis_dimension,
        minimum = 5L, maximum = 20L, growth_exponent = 1/3)
    formula <- binary_make_gam_formula(".outcome", training, predictors, local_basis)
    environment(formula) <- environment()
    fit <- suppressWarnings(mgcv::gam(formula, data = training, family = family, weights = obsWeights,
        method = "REML", select = TRUE))
    prediction <- as.numeric(stats::predict(fit, newdata = as.data.frame(newX), type = "response"))
    output <- list(pred = prediction, fit = list(object = fit))
    class(output$fit) <- "SL.binary_mgcv_additive"
    output
}

SL.binary_mgcv_interactions <- function(Y, X, newX, family, obsWeights, basis_dimension = NULL, interaction_basis_dimension = 3L,
    ...) {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("SL.binary_mgcv_interactions requires package 'mgcv'.")
    }
    training <- as.data.frame(X)
    training$.outcome <- as.numeric(Y)
    predictors <- setdiff(names(training), ".outcome")
    local_basis <- binary_resolve_gam_basis_dimension(training, predictors, basis_dimension = basis_dimension,
        minimum = 5L, maximum = 20L, growth_exponent = 1/3)
    formula <- binary_make_interaction_gam_formula(".outcome", training, predictors, local_basis, interaction_basis_dimension = interaction_basis_dimension)
    environment(formula) <- environment()
    fit <- suppressWarnings(mgcv::gam(formula, data = training, family = family, weights = obsWeights,
        method = "REML", select = TRUE))
    prediction <- as.numeric(stats::predict(fit, newdata = as.data.frame(newX), type = "response"))
    output <- list(pred = prediction, fit = list(object = fit))
    class(output$fit) <- "SL.binary_mgcv_interactions"
    output
}

SL.binary_ranger <- function(Y, X, newX, family, obsWeights, id, num_trees = 500L, min_node_size = 20L,
    ...) {
    SuperLearner::SL.ranger(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, id = id,
        num.trees = as.integer(num_trees), min.node.size = as.integer(min_node_size), num.threads = 1L,
        verbose = FALSE)
}

SL.binary_xgboost <- function(Y, X, newX, family, obsWeights, id, ntrees = 250L, max_depth = 3L, shrinkage = 0.05,
    minobspernode = 10L, ...) {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop("SL.binary_xgboost requires package 'xgboost'.")
    }
    if (!identical(family$family, "binomial")) {
        stop("SL.binary_xgboost currently supports binomial outcomes only.")
    }
    training_design <- as.matrix(X)
    validation_design <- as.matrix(newX)
    storage.mode(training_design) <- "double"
    storage.mode(validation_design) <- "double"
    training_matrix <- xgboost::xgb.DMatrix(data = training_design, label = as.numeric(Y), weight = as.numeric(obsWeights))
    fit <- xgboost::xgb.train(params = list(objective = "binary:logistic", eval_metric = "logloss", max_depth = as.integer(max_depth),
        eta = as.numeric(shrinkage), min_child_weight = as.integer(minobspernode), tree_method = "hist",
        nthread = 1L), data = training_matrix, nrounds = as.integer(ntrees), verbose = 0)
    prediction <- as.numeric(stats::predict(fit, newdata = validation_design))
    output <- list(pred = prediction, fit = list(object = fit))
    class(output$fit) <- "SL.binary_xgboost"
    output
}

subset_nuisance <- function(nuisance, index) {
    list(e = nuisance$e[index], pi0 = nuisance$pi0[index], pi1 = nuisance$pi1[index])
}

