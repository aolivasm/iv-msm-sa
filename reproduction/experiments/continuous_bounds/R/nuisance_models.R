# Main-text reproduction dependencies; published formulas retained.
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

compute_nuisance_diagnostics_simple <- function(data) {
    accuracy <- function(estimated, truth) {
        correlation <- if (stats::sd(estimated) > 0 && stats::sd(truth) > 0) {
            stats::cor(estimated, truth)
        }
        else {
            NA_real_
        }
        data.frame(bias = mean(estimated - truth), rmse = sqrt(mean((estimated - truth)^2)), mae = mean(abs(estimated -
            truth)), correlation = correlation, mean_estimated = mean(estimated), mean_true = mean(truth))
    }
    estimated_nuisance <- list(e = data$e_hat, pi0 = data$pi0_hat, pi1 = data$pi1_hat)
    if (!all(c("e_true", "pi0_true", "pi1_true") %in% names(data))) {
        unavailable <- data.frame(bias = NA_real_, rmse = NA_real_, mae = NA_real_, correlation = NA_real_,
            mean_estimated = NA_real_, mean_true = NA_real_)
        return(list(e = unavailable, pi0 = unavailable, pi1 = unavailable, w1 = unavailable, w0 = unavailable))
    }
    oracle_nuisance <- list(e = data$e_true, pi0 = data$pi0_true, pi1 = data$pi1_true)
    w1_hat <- compute_riesz_weight(estimated_nuisance, t_value = 1)
    w1_true <- compute_riesz_weight(oracle_nuisance, t_value = 1)
    w0_hat <- compute_riesz_weight(estimated_nuisance, t_value = 0)
    w0_true <- compute_riesz_weight(oracle_nuisance, t_value = 0)
    list(e = accuracy(data$e_hat, data$e_true), pi0 = accuracy(data$pi0_hat, data$pi0_true), pi1 = accuracy(data$pi1_hat,
        data$pi1_true), w1 = accuracy(w1_hat, w1_true), w0 = accuracy(w0_hat, w0_true))
}

fit_nuisance_models <- function(train, x_cols = "X", positivity_epsilon = 0.05, learner = c("gam", "ranger",
    "glm_spline", "library_cv", "super_learner", "integrated_parametric"), treatment_spline_df = 5L,
    gam_basis_dimension = NULL, gam_minimum_basis_dimension = 5L, gam_maximum_basis_dimension = 20L,
    gam_basis_growth_exponent = 1/3, gam_method = "REML", integrated_confounder_strength = 1, integrated_quadrature_nodes = 40L,
    integrated_maxit = 500L, library_candidates = c("glm_spline", "gam", "ranger"), library_cv_folds = 3L,
    library_selection_rule = c("one_se", "minimum"), library_ranger_num_trees = 500L, library_seed = 1L,
    super_learner_candidates = NULL, super_learner_cv_folds = 3L, super_learner_ranger_num_trees = 500L,
    super_learner_hal_max_degree = 2L) {
    learner <- match.arg(learner)
    library_selection_rule <- match.arg(library_selection_rule)
    if (learner == "glm_spline") {
        return(fit_nuisance_models_glm(train = train, x_cols = x_cols, positivity_epsilon = positivity_epsilon,
            treatment_spline_df = treatment_spline_df))
    }
    if (learner == "gam") {
        return(fit_nuisance_models_gam(train = train, x_cols = x_cols, positivity_epsilon = positivity_epsilon,
            basis_dimension = gam_basis_dimension, minimum_basis_dimension = gam_minimum_basis_dimension,
            maximum_basis_dimension = gam_maximum_basis_dimension, basis_growth_exponent = gam_basis_growth_exponent,
            method = gam_method))
    }
    if (learner == "ranger") {
        return(fit_nuisance_models_ranger(train = train, x_cols = x_cols, positivity_epsilon = positivity_epsilon))
    }
    if (learner == "library_cv") {
        return(fit_nuisance_models_library_cv(train = train, x_cols = x_cols, positivity_epsilon = positivity_epsilon,
            candidates = library_candidates, cv_folds = library_cv_folds, selection_rule = library_selection_rule,
            treatment_spline_df = treatment_spline_df, gam_basis_dimension = gam_basis_dimension, gam_minimum_basis_dimension = gam_minimum_basis_dimension,
            gam_maximum_basis_dimension = gam_maximum_basis_dimension, gam_basis_growth_exponent = gam_basis_growth_exponent,
            gam_method = gam_method, ranger_num_trees = library_ranger_num_trees, seed = library_seed))
    }
    if (learner == "super_learner") {
        return(fit_nuisance_models_super_learner(train = train, x_cols = x_cols, positivity_epsilon = positivity_epsilon,
            candidates = super_learner_candidates, cv_folds = super_learner_cv_folds, ranger_num_trees = super_learner_ranger_num_trees,
            hal_max_degree = super_learner_hal_max_degree, seed = library_seed))
    }
    fit_nuisance_models_integrated_parametric(train = train, x_cols = x_cols, positivity_epsilon = positivity_epsilon,
        confounder_strength = integrated_confounder_strength, quadrature_nodes = integrated_quadrature_nodes,
        maxit = integrated_maxit)
}

fit_nuisance_models_gam <- function(train, x_cols = "X", positivity_epsilon = 0.05, basis_dimension = NULL,
    minimum_basis_dimension = 5L, maximum_basis_dimension = 20L, basis_growth_exponent = 1/3, method = "REML") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("The optional nuisance_learner='gam' requires the recommended mgcv package.")
    }
    non_numeric <- x_cols[!vapply(train[, x_cols, drop = FALSE], is.numeric, logical(1))]
    if (length(non_numeric) > 0L) {
        stop("GAM nuisance predictors must be numeric: ", paste(non_numeric, collapse = ", "))
    }
    if (length(method) != 1L || !nzchar(method)) {
        stop("nuisance_gam_method must be one nonempty string.")
    }
    basis_for <- function(data) {
        resolve_nuisance_gam_basis_dimension(data = data, x_cols = x_cols, basis_dimension = basis_dimension,
            minimum = minimum_basis_dimension, maximum = maximum_basis_dimension, growth_exponent = basis_growth_exponent)
    }
    make_formula <- function(response, data, selected_basis_dimension) {
        smooth_terms <- vapply(x_cols, function(variable) {
            unique_count <- length(unique(data[[variable]][is.finite(data[[variable]])]))
            if (unique_count <= 2L) {
                return(sprintf("`%s`", variable))
            }
            local_k <- min(selected_basis_dimension, unique_count - 1L)
            local_k <- max(3L, local_k)
            sprintf("s(`%s`, bs='cr', k=%d)", variable, local_k)
        }, character(1))
        formula <- stats::as.formula(paste0("`", response, "` ~ ", paste(smooth_terms, collapse = " + ")))
        environment(formula) <- asNamespace("mgcv")
        formula
    }
    e_basis_dimension <- basis_for(train)
    z0_training <- train[train$Z == 0, , drop = FALSE]
    z1_training <- train[train$Z == 1, , drop = FALSE]
    pi0_basis_dimension <- basis_for(z0_training)
    pi1_basis_dimension <- basis_for(z1_training)
    instrument_fit <- mgcv::gam(formula = make_formula("Z", train, e_basis_dimension), data = train,
        family = stats::binomial(), method = method, select = TRUE)
    treatment_fit_z0 <- mgcv::gam(formula = make_formula("T_obs", z0_training, pi0_basis_dimension),
        data = z0_training, family = stats::binomial(), method = method, select = TRUE)
    treatment_fit_z1 <- mgcv::gam(formula = make_formula("T_obs", z1_training, pi1_basis_dimension),
        data = z1_training, family = stats::binomial(), method = method, select = TRUE)
    predict_nuisance <- function(newdata) {
        p_z1 <- stats::predict(instrument_fit, newdata = newdata, type = "response")
        pi0 <- stats::predict(treatment_fit_z0, newdata = newdata, type = "response")
        pi1 <- stats::predict(treatment_fit_z1, newdata = newdata, type = "response")
        list(e = clip_probability(1 - p_z1, positivity_epsilon), pi0 = clip_probability(pi0, positivity_epsilon),
            pi1 = clip_probability(pi1, positivity_epsilon))
    }
    gam_diagnostics <- function(fit, component, selected_basis_dimension) {
        data.frame(component = component, learner = "gam", converged = isTRUE(fit$converged), model_rank = fit$rank,
            effective_df = sum(fit$edf), basis_dimension = selected_basis_dimension, stringsAsFactors = FALSE)
    }
    list(predict = predict_nuisance, instrument_fit = instrument_fit, treatment_fit_z0 = treatment_fit_z0,
        treatment_fit_z1 = treatment_fit_z1, learner = "gam", diagnostics = rbind(gam_diagnostics(instrument_fit,
            "e", e_basis_dimension), gam_diagnostics(treatment_fit_z0, "pi0", pi0_basis_dimension), gam_diagnostics(treatment_fit_z1,
            "pi1", pi1_basis_dimension)), configuration = list(basis_selection = if (is.null(basis_dimension)) {
            "adaptive_training_stratum_size"
        } else {
            "fixed"
        }, minimum_basis_dimension = as.integer(minimum_basis_dimension), maximum_basis_dimension = as.integer(maximum_basis_dimension),
            basis_growth_exponent = basis_growth_exponent, selected_basis_dimensions = c(e = e_basis_dimension,
                pi0 = pi0_basis_dimension, pi1 = pi1_basis_dimension)))
}

fit_nuisance_models_glm <- function(train, x_cols = "X", positivity_epsilon = 0.05, treatment_spline_df = 5) {
    instrument_formula <- stats::reformulate(x_cols, response = "Z")
    instrument_fit <- stats::glm(formula = instrument_formula, data = train, family = stats::binomial(),
        control = stats::glm.control(maxit = 100))
    treatment_formula <- make_additive_spline_formula(response = "T_obs", predictors = x_cols, spline_df = treatment_spline_df,
        data = train)
    treatment_fit_z0 <- stats::glm(formula = treatment_formula, data = train[train$Z == 0, , drop = FALSE],
        family = stats::binomial(), control = stats::glm.control(maxit = 100))
    treatment_fit_z1 <- stats::glm(formula = treatment_formula, data = train[train$Z == 1, , drop = FALSE],
        family = stats::binomial(), control = stats::glm.control(maxit = 100))
    predict_nuisance <- function(newdata) {
        p_z1 <- stats::predict(instrument_fit, newdata = newdata, type = "response")
        pi0 <- stats::predict(treatment_fit_z0, newdata = newdata, type = "response")
        pi1 <- stats::predict(treatment_fit_z1, newdata = newdata, type = "response")
        list(e = clip_probability(1 - p_z1, positivity_epsilon), pi0 = clip_probability(pi0, positivity_epsilon),
            pi1 = clip_probability(pi1, positivity_epsilon))
    }
    list(predict = predict_nuisance, instrument_fit = instrument_fit, treatment_fit_z0 = treatment_fit_z0,
        treatment_fit_z1 = treatment_fit_z1, learner = "glm_spline", diagnostics = data.frame(component = c("e",
            "pi0", "pi1"), learner = "glm_spline", converged = c(isTRUE(instrument_fit$converged), isTRUE(treatment_fit_z0$converged),
            isTRUE(treatment_fit_z1$converged)), model_rank = c(instrument_fit$rank, treatment_fit_z0$rank,
            treatment_fit_z1$rank), effective_df = c(instrument_fit$rank, treatment_fit_z0$rank, treatment_fit_z1$rank),
            stringsAsFactors = FALSE))
}

fit_nuisance_models_integrated_parametric <- function(train, x_cols = "X", positivity_epsilon = 0.05,
    confounder_strength = 1, quadrature_nodes = 40L, maxit = 500L) {
    non_numeric <- x_cols[!vapply(train[, x_cols, drop = FALSE], is.numeric, logical(1))]
    if (length(non_numeric) > 0L) {
        stop("Integrated-parametric nuisance predictors must be numeric: ", paste(non_numeric, collapse = ", "))
    }
    if (length(confounder_strength) != 1L || !is.finite(confounder_strength)) {
        stop("confounder_strength must be one finite number.")
    }
    if (length(quadrature_nodes) != 1L || !is.finite(quadrature_nodes) || quadrature_nodes < 2L || abs(quadrature_nodes -
        round(quadrature_nodes)) > sqrt(.Machine$double.eps)) {
        stop("nuisance_integrated_nodes must be one integer of at least 2.")
    }
    if (length(maxit) != 1L || !is.finite(maxit) || maxit < 1L || abs(maxit - round(maxit)) > sqrt(.Machine$double.eps)) {
        stop("nuisance_integrated_maxit must be one positive integer.")
    }
    instrument_formula <- stats::reformulate(x_cols, response = "Z")
    instrument_fit <- stats::glm(formula = instrument_formula, data = train, family = stats::binomial(),
        control = stats::glm.control(maxit = 100))
    rule <- gauss_legendre_rule(n_nodes = as.integer(quadrature_nodes), lower = -2, upper = 2)
    integration_weights <- rule$weights * d_u(rule$nodes)
    integration_weights <- integration_weights/sum(integration_weights)
    treatment <- as.numeric(train$T_obs)
    check_binary(treatment, "T_obs")
    treatment_design <- cbind(`(Intercept)` = 1, as.matrix(train[, x_cols, drop = FALSE]), Z = as.numeric(train$Z))
    storage.mode(treatment_design) <- "double"
    marginal_probability <- function(coefficients, design) {
        linear_predictor <- as.vector(design %*% coefficients)
        probability_matrix <- stats::plogis(outer(linear_predictor, confounder_strength * rule$nodes,
            "+"))
        as.vector(probability_matrix %*% integration_weights)
    }
    objective <- function(coefficients) {
        probability <- clip_probability(marginal_probability(coefficients, treatment_design), 1e-10)
        -sum(treatment * log(probability) + (1 - treatment) * log1p(-probability))
    }
    gradient <- function(coefficients) {
        linear_predictor <- as.vector(treatment_design %*% coefficients)
        conditional_probability <- stats::plogis(outer(linear_predictor, confounder_strength * rule$nodes,
            "+"))
        probability <- as.vector(conditional_probability %*% integration_weights)
        probability_derivative <- as.vector((conditional_probability * (1 - conditional_probability)) %*%
            integration_weights)
        probability <- clip_probability(probability, 1e-10)
        score_by_row <- (probability - treatment) * probability_derivative/(probability * (1 - probability))
        as.vector(crossprod(treatment_design, score_by_row))
    }
    initial_formula <- stats::reformulate(c(x_cols, "Z"), response = "T_obs")
    initial_fit <- suppressWarnings(stats::glm(formula = initial_formula, data = train, family = stats::binomial(),
        control = stats::glm.control(maxit = 100)))
    initial_coefficients <- stats::coef(initial_fit)
    initial_coefficients[!is.finite(initial_coefficients)] <- 0
    if (length(initial_coefficients) != ncol(treatment_design)) {
        initial_coefficients <- rep(0, ncol(treatment_design))
    }
    names(initial_coefficients) <- colnames(treatment_design)
    treatment_fit <- stats::optim(par = initial_coefficients, fn = objective, gr = gradient, method = "BFGS",
        control = list(maxit = as.integer(maxit), reltol = 1e-10))
    names(treatment_fit$par) <- colnames(treatment_design)
    make_treatment_design <- function(newdata, z_value) {
        design <- cbind(`(Intercept)` = 1, as.matrix(newdata[, x_cols, drop = FALSE]), Z = rep(z_value,
            nrow(newdata)))
        storage.mode(design) <- "double"
        design
    }
    predict_nuisance <- function(newdata) {
        p_z1 <- stats::predict(instrument_fit, newdata = newdata, type = "response")
        pi0 <- marginal_probability(treatment_fit$par, make_treatment_design(newdata, 0))
        pi1 <- marginal_probability(treatment_fit$par, make_treatment_design(newdata, 1))
        list(e = clip_probability(1 - p_z1, positivity_epsilon), pi0 = clip_probability(pi0, positivity_epsilon),
            pi1 = clip_probability(pi1, positivity_epsilon))
    }
    treatment_converged <- identical(treatment_fit$convergence, 0L)
    treatment_rank <- ncol(treatment_design)
    list(predict = predict_nuisance, instrument_fit = instrument_fit, treatment_fit = treatment_fit,
        treatment_fit_initial_glm = initial_fit, learner = "integrated_parametric", diagnostics = data.frame(component = c("e",
            "pi0", "pi1"), learner = "integrated_parametric", converged = c(isTRUE(instrument_fit$converged),
            treatment_converged, treatment_converged), model_rank = c(instrument_fit$rank, treatment_rank,
            treatment_rank), effective_df = c(instrument_fit$rank, treatment_rank, treatment_rank), n_training = rep(nrow(train),
            3L), stringsAsFactors = FALSE), configuration = list(confounder_strength = confounder_strength,
            quadrature_nodes = as.integer(quadrature_nodes), optimization_convergence = treatment_fit$convergence,
            optimization_message = treatment_fit$message, estimated_coefficients = treatment_fit$par))
}

fit_nuisance_models_library_cv <- function(train, x_cols = "X", positivity_epsilon = 0.05, candidates = c("glm_spline",
    "gam", "ranger"), cv_folds = 3L, selection_rule = c("one_se", "minimum"), treatment_spline_df = 5L,
    gam_basis_dimension = NULL, gam_minimum_basis_dimension = 5L, gam_maximum_basis_dimension = 20L,
    gam_basis_growth_exponent = 1/3, gam_method = "REML", ranger_num_trees = 500L, seed = 1L) {
    selection_rule <- match.arg(selection_rule)
    candidates <- unique(as.character(candidates))
    allowed <- c("glm_spline", "gam", "ranger")
    if (length(candidates) < 1L || any(!candidates %in% allowed)) {
        stop("Primitive nuisance library candidates must be drawn from: ", paste(allowed, collapse = ", "),
            ".")
    }
    cv_folds <- as.integer(cv_folds)
    if (length(cv_folds) != 1L || is.na(cv_folds) || cv_folds < 2L) {
        stop("Primitive nuisance library CV requires at least two folds.")
    }
    fit_candidate <- function(data, candidate, local_seed) {
        switch(candidate, glm_spline = fit_nuisance_models_glm(train = data, x_cols = x_cols, positivity_epsilon = positivity_epsilon,
            treatment_spline_df = treatment_spline_df), gam = fit_nuisance_models_gam(train = data, x_cols = x_cols,
            positivity_epsilon = positivity_epsilon, basis_dimension = gam_basis_dimension, minimum_basis_dimension = gam_minimum_basis_dimension,
            maximum_basis_dimension = gam_maximum_basis_dimension, basis_growth_exponent = gam_basis_growth_exponent,
            method = gam_method), ranger = fit_nuisance_models_ranger(train = data, x_cols = x_cols,
            positivity_epsilon = positivity_epsilon, num_trees = ranger_num_trees, seed = local_seed))
    }
    fold_id <- make_tz_stratified_folds(train, K = cv_folds, seed = as.integer(seed))
    risk_rows <- list()
    counter <- 1L
    for (candidate_index in seq_along(candidates)) {
        candidate <- candidates[candidate_index]
        for (fold in seq_len(cv_folds)) {
            validation <- fold_id == fold
            fit <- fit_candidate(train[!validation, , drop = FALSE], candidate, as.integer(seed + 10000L *
                candidate_index + fold))
            prediction <- fit$predict(train[validation, , drop = FALSE])
            local_data <- train[validation, , drop = FALSE]
            component_specifications <- list(e = list(index = rep(TRUE, nrow(local_data)), outcome = 1 -
                local_data$Z, prediction = prediction$e), pi0 = list(index = local_data$Z == 0, outcome = local_data$T_obs,
                prediction = prediction$pi0), pi1 = list(index = local_data$Z == 1, outcome = local_data$T_obs,
                prediction = prediction$pi1))
            for (component in names(component_specifications)) {
                specification <- component_specifications[[component]]
                index <- specification$index
                if (!any(index)) {
                  stop("A primitive nuisance CV fold contains an empty stratum.")
                }
                risk_rows[[counter]] <- data.frame(component = component, candidate = candidate, fold = fold,
                  risk = binary_cross_entropy(specification$outcome[index], specification$prediction[index]),
                  observations = sum(index), stringsAsFactors = FALSE)
                counter <- counter + 1L
            }
        }
    }
    fold_risk <- do.call(rbind, risk_rows)
    groups <- split(fold_risk, interaction(fold_risk$component, fold_risk$candidate, drop = TRUE, lex.order = TRUE))
    cv_table <- do.call(rbind, lapply(groups, function(group) {
        data.frame(component = group$component[1L], candidate = group$candidate[1L], mean_risk = mean(group$risk),
            standard_error = if (nrow(group) > 1L) {
                stats::sd(group$risk)/sqrt(nrow(group))
            }
            else {
                0
            }, stringsAsFactors = FALSE)
    }))
    rownames(cv_table) <- NULL
    complexity <- c(glm_spline = 1L, gam = 2L, ranger = 3L)
    selected <- stats::setNames(character(3L), c("e", "pi0", "pi1"))
    cv_table$selected <- FALSE
    for (component in names(selected)) {
        local <- which(cv_table$component == component)
        minimum <- local[which.min(cv_table$mean_risk[local])]
        choice <- minimum
        if (selection_rule == "one_se") {
            cutoff <- cv_table$mean_risk[minimum] + cv_table$standard_error[minimum]
            eligible <- local[cv_table$mean_risk[local] <= cutoff]
            choice <- eligible[which.min(complexity[cv_table$candidate[eligible]])]
        }
        cv_table$selected[choice] <- TRUE
        selected[component] <- cv_table$candidate[choice]
    }
    selected_candidates <- unique(unname(selected))
    final_fits <- stats::setNames(lapply(seq_along(selected_candidates), function(index) {
        fit_candidate(train, selected_candidates[index], as.integer(seed + 900000L + index))
    }), selected_candidates)
    predict_nuisance <- function(newdata) {
        candidate_predictions <- lapply(final_fits, function(fit) {
            fit$predict(newdata)
        })
        list(e = candidate_predictions[[selected[["e"]]]]$e, pi0 = candidate_predictions[[selected[["pi0"]]]]$pi0,
            pi1 = candidate_predictions[[selected[["pi1"]]]]$pi1)
    }
    selected_diagnostics <- do.call(rbind, lapply(names(selected), function(component) {
        candidate <- selected[[component]]
        diagnostic <- final_fits[[candidate]]$diagnostics
        diagnostic <- diagnostic[diagnostic$component == component, , drop = FALSE]
        local_cv <- cv_table[cv_table$component == component & cv_table$candidate == candidate, , drop = FALSE]
        value_or <- function(name, default) {
            if (name %in% names(diagnostic))
                diagnostic[[name]][1L]
            else default
        }
        data.frame(component = component, learner = "library_cv", selected_learner = candidate, converged = value_or("converged",
            TRUE), model_rank = value_or("model_rank", NA_real_), effective_df = value_or("effective_df",
            NA_real_), n_training = value_or("n_training", if (component == "e")
            nrow(train)
        else sum(train$Z == if (component == "pi0") 0 else 1)), cv_mean_risk = local_cv$mean_risk[1L],
            cv_standard_error = local_cv$standard_error[1L], stringsAsFactors = FALSE)
    }))
    rownames(selected_diagnostics) <- NULL
    list(predict = predict_nuisance, fits = final_fits, learner = "library_cv", selected_learners = selected,
        diagnostics = selected_diagnostics, cv_table = cv_table, configuration = list(candidates = candidates,
            cv_folds = cv_folds, selection_rule = selection_rule, selected_learners = selected))
}

fit_nuisance_models_ranger <- function(train, x_cols = "X", positivity_epsilon = 0.05, num_trees = 500L,
    min_node_size = NULL, seed = 1L) {
    if (!requireNamespace("ranger", quietly = TRUE)) {
        stop("The optional nuisance_learner='ranger' requires the ranger package.")
    }
    non_numeric <- x_cols[!vapply(train[, x_cols, drop = FALSE], is.numeric, logical(1))]
    if (length(non_numeric) > 0L) {
        stop("Ranger nuisance predictors must be numeric: ", paste(non_numeric, collapse = ", "))
    }
    num_trees <- as.integer(num_trees)
    if (length(num_trees) != 1L || is.na(num_trees) || num_trees < 1L) {
        stop("num_trees must be one positive integer.")
    }
    fit_probability_forest <- function(data, response, local_seed) {
        outcome <- as.numeric(data[[response]])
        check_binary(outcome, response)
        if (length(unique(outcome)) == 1L) {
            probability <- (sum(outcome) + 0.5)/(length(outcome) + 1)
            return(list(predict = function(newdata) rep(probability, nrow(newdata)), fit = NULL, min_node_size = NA_integer_))
        }
        local_min_node_size <- if (is.null(min_node_size)) {
            max(5L, as.integer(floor(nrow(data)^(1/4))))
        }
        else {
            as.integer(min_node_size)
        }
        if (length(local_min_node_size) != 1L || is.na(local_min_node_size) || local_min_node_size <
            1L) {
            stop("min_node_size must be NULL or one positive integer.")
        }
        ranger_data <- data[, c(response, x_cols), drop = FALSE]
        ranger_data[[response]] <- factor(outcome, levels = c(0, 1))
        fit <- ranger::ranger(formula = stats::reformulate(x_cols, response = response), data = ranger_data,
            probability = TRUE, num.trees = num_trees, mtry = max(1L, floor(sqrt(length(x_cols)))), min.node.size = local_min_node_size,
            seed = as.integer(local_seed), respect.unordered.factors = "order", num.threads = 1L)
        list(predict = function(newdata) {
            prediction <- predict(fit, data = newdata)$predictions
            if (is.null(dim(prediction))) return(as.numeric(prediction))
            if ("1" %in% colnames(prediction)) {
                return(as.numeric(prediction[, "1"]))
            }
            as.numeric(prediction[, ncol(prediction)])
        }, fit = fit, min_node_size = local_min_node_size)
    }
    z0_training <- train[train$Z == 0, , drop = FALSE]
    z1_training <- train[train$Z == 1, , drop = FALSE]
    instrument_fit <- fit_probability_forest(train, "Z", seed)
    treatment_fit_z0 <- fit_probability_forest(z0_training, "T_obs", seed + 1L)
    treatment_fit_z1 <- fit_probability_forest(z1_training, "T_obs", seed + 2L)
    predict_nuisance <- function(newdata) {
        list(e = clip_probability(1 - instrument_fit$predict(newdata), positivity_epsilon), pi0 = clip_probability(treatment_fit_z0$predict(newdata),
            positivity_epsilon), pi1 = clip_probability(treatment_fit_z1$predict(newdata), positivity_epsilon))
    }
    list(predict = predict_nuisance, instrument_fit = instrument_fit$fit, treatment_fit_z0 = treatment_fit_z0$fit,
        treatment_fit_z1 = treatment_fit_z1$fit, learner = "ranger", diagnostics = data.frame(component = c("e",
            "pi0", "pi1"), learner = "ranger", converged = TRUE, model_rank = rep(length(x_cols), 3L),
            effective_df = rep(NA_real_, 3L), n_training = c(nrow(train), nrow(z0_training), nrow(z1_training)),
            min_node_size = c(instrument_fit$min_node_size, treatment_fit_z0$min_node_size, treatment_fit_z1$min_node_size),
            num_trees = rep(num_trees, 3L), stringsAsFactors = FALSE))
}

fit_nuisance_models_super_learner <- function(train, x_cols = "X", positivity_epsilon = 0.05, candidates = NULL,
    cv_folds = 3L, ranger_num_trees = 500L, hal_max_degree = 2L, seed = 1L) {
    required <- c("binary_fit_regression_super_learner", "binary_default_super_learner_candidates")
    missing <- required[!vapply(required, exists, logical(1), mode = "function")]
    if (length(missing)) {
        stop("Source R/binary_outcome_bounds.R before requesting the primitive ", "Super Learner; missing: ",
            paste(missing, collapse = ", "), ".")
    }
    if (is.null(candidates)) {
        candidates <- binary_default_super_learner_candidates()
    }
    candidates <- unique(as.character(candidates))
    fit_component <- function(data, response, local_seed) {
        binary_fit_regression_super_learner(data = data, response = response, x_cols = x_cols, candidates = candidates,
            cv_folds = cv_folds, ranger_num_trees = ranger_num_trees, hal_max_degree = hal_max_degree,
            seed = local_seed)
    }
    z0_training <- train[train$Z == 0, , drop = FALSE]
    z1_training <- train[train$Z == 1, , drop = FALSE]
    fits <- list(e = fit_component(train, "Z", seed + 1L), pi0 = fit_component(z0_training, "T_obs",
        seed + 2L), pi1 = fit_component(z1_training, "T_obs", seed + 3L))
    training_sizes <- c(e = nrow(train), pi0 = nrow(z0_training), pi1 = nrow(z1_training))
    diagnostics <- do.call(rbind, lapply(names(fits), function(component) {
        fit <- fits[[component]]
        row <- data.frame(component = component, learner = "super_learner", selected_learner = fit$selected_learner,
            converged = isTRUE(fit$converged), model_rank = as.integer(fit$rank), effective_df = NA_real_,
            n_training = as.integer(training_sizes[[component]]), cv_mean_risk = as.numeric(fit$cv_mean_risk),
            stringsAsFactors = FALSE)
        weights <- fit$ensemble_weights
        risk_table <- fit$cv_table
        for (candidate in candidates) {
            row[[paste0("weight_", candidate)]] <- if (!is.null(weights) && candidate %in% names(weights))
                as.numeric(weights[[candidate]])
            else 0
            local_risk <- if (!is.null(risk_table) && nrow(risk_table) && candidate %in% risk_table$candidate) {
                risk_table$mean_risk[match(candidate, risk_table$candidate)]
            }
            else {
                NA_real_
            }
            row[[paste0("risk_", candidate)]] <- as.numeric(local_risk)
        }
        row
    }))
    rownames(diagnostics) <- NULL
    cv_tables <- do.call(rbind, lapply(names(fits), function(component) {
        table <- fits[[component]]$cv_table
        if (is.null(table) || !nrow(table))
            return(NULL)
        table$component <- component
        table
    }))
    predict_nuisance <- function(newdata) {
        list(e = clip_probability(1 - fits$e$predict(newdata), positivity_epsilon), pi0 = clip_probability(fits$pi0$predict(newdata),
            positivity_epsilon), pi1 = clip_probability(fits$pi1$predict(newdata), positivity_epsilon))
    }
    list(predict = predict_nuisance, fits = fits, learner = "super_learner", diagnostics = diagnostics,
        cv_table = cv_tables, configuration = list(candidates = candidates, cv_folds = as.integer(cv_folds),
            ranger_num_trees = as.integer(ranger_num_trees), hal_max_degree = as.integer(hal_max_degree),
            metalearner = "simplex_probability_mixture_bernoulli_log_loss"))
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

resolve_nuisance_gam_basis_dimension <- function(data, x_cols, basis_dimension = NULL, minimum = 5L,
    maximum = 20L, growth_exponent = 1/3) {
    if (!is.null(basis_dimension)) {
        basis_dimension <- as.integer(basis_dimension)
        if (length(basis_dimension) != 1L || is.na(basis_dimension) || basis_dimension < 3L) {
            stop("nuisance_gam_k must be NULL or one integer of at least 3.")
        }
        return(basis_dimension)
    }
    minimum <- as.integer(minimum)
    maximum <- as.integer(maximum)
    if (length(minimum) != 1L || is.na(minimum) || minimum < 3L || length(maximum) != 1L || is.na(maximum) ||
        maximum < minimum || length(growth_exponent) != 1L || !is.finite(growth_exponent) || growth_exponent <=
        0 || growth_exponent > 1) {
        stop("Invalid adaptive nuisance-GAM basis controls.")
    }
    observations_per_covariate <- nrow(data)/max(1L, length(x_cols))
    proposed <- ceiling(observations_per_covariate^growth_exponent)
    as.integer(min(maximum, max(minimum, proposed)))
}

sfp_contract <- function() list(primitive_library = c("glm_main", "glm_interactions", "glmnet_interactions",
    "ranger", "xgboost"))

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

