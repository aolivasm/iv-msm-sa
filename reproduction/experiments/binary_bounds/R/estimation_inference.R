# Main-text reproduction dependencies; published formulas retained.
binary_arm_endpoint <- function(nuisance, Lambda0, Lambda1, epsilon = 0, compute_gradient = FALSE) {
    if (!is.finite(Lambda0) || !is.finite(Lambda1) || Lambda0 < 1 || Lambda1 < 1) {
        stop("Lambda0 and Lambda1 must be finite and at least one.")
    }
    if (!is.finite(epsilon) || epsilon < 0 || (epsilon > 0 && epsilon >= 1/Lambda0)) {
        stop("epsilon must be nonnegative and smaller than 1/Lambda0.")
    }
    nuisance <- binary_recycle_arm_nuisance(nuisance)
    r0 <- nuisance$r0
    r1 <- nuisance$r1
    a0 <- nuisance$a0
    a1 <- nuisance$a1
    n <- length(r0)
    if (any(r0 <= 0 | r0 >= 1 | r1 <= 0 | r1 >= 1 | a0 <= 0 | a0 >= 1 | a1 <= 0 | a1 >= 1)) {
        stop("Arm propensities and outcome regressions must lie strictly in (0,1).")
    }
    gradient_names <- c("r0", "r1", "a0", "a1")
    zero_gradient <- matrix(0, nrow = n, ncol = length(gradient_names), dimnames = list(NULL, gradient_names))
    gamma_one <- a1/a0
    gamma_zero <- (1 - a1)/(1 - a0)
    gamma_one_gradient <- zero_gradient
    gamma_one_gradient[, "a0"] <- -gamma_one/a0
    gamma_one_gradient[, "a1"] <- 1/a0
    gamma_zero_gradient <- zero_gradient
    gamma_zero_gradient[, "a0"] <- gamma_zero/(1 - a0)
    gamma_zero_gradient[, "a1"] <- -1/(1 - a0)
    construct_envelope <- function(gamma, gamma_gradient) {
        candidate <- function(c_value) {
            h_value <- r1 + (1 - r1) * c_value
            value <- (gamma * h_value - r0)/(1 - r0)
            derivative <- gamma_gradient * (h_value/(1 - r0))
            derivative[, "r0"] <- derivative[, "r0"] + (gamma * h_value - 1)/(1 - r0)^2
            derivative[, "r1"] <- derivative[, "r1"] + gamma * (1 - c_value)/(1 - r0)
            list(value = value, derivative = derivative)
        }
        if (Lambda0 == 1 && Lambda1 == 1) {
            return(list(L = rep(1, n), U = rep(1, n), derivative_L = zero_gradient, derivative_U = zero_gradient,
                candidate_lower = rep(1, n), candidate_upper = rep(1, n)))
        }
        candidate_lower <- candidate(1/Lambda1)
        candidate_upper <- candidate(Lambda1)
        lower <- binary_smooth_max_below(rep(1/Lambda0, n), candidate_lower$value, epsilon)
        upper <- binary_smooth_min_above(rep(Lambda0, n), candidate_upper$value, epsilon)
        list(L = lower$value, U = upper$value, derivative_L = candidate_lower$derivative * lower$weight_second,
            derivative_U = candidate_upper$derivative * upper$weight_second, candidate_lower = candidate_lower$value,
            candidate_upper = candidate_upper$value)
    }
    envelope_one <- construct_envelope(gamma_one, gamma_one_gradient)
    envelope_zero <- construct_envelope(gamma_zero, gamma_zero_gradient)
    lower_first <- a0 * envelope_one$L
    lower_second <- 1 - (1 - a0) * envelope_zero$U
    upper_first <- a0 * envelope_one$U
    upper_second <- 1 - (1 - a0) * envelope_zero$L
    lower_pair <- binary_smooth_max_below(lower_first, lower_second, epsilon)
    upper_pair <- binary_smooth_min_above(upper_first, upper_second, epsilon)
    derivative_lower_first <- envelope_one$derivative_L * a0
    derivative_lower_first[, "a0"] <- derivative_lower_first[, "a0"] + envelope_one$L
    derivative_lower_second <- -envelope_zero$derivative_U * (1 - a0)
    derivative_lower_second[, "a0"] <- derivative_lower_second[, "a0"] + envelope_zero$U
    derivative_upper_first <- envelope_one$derivative_U * a0
    derivative_upper_first[, "a0"] <- derivative_upper_first[, "a0"] + envelope_one$U
    derivative_upper_second <- -envelope_zero$derivative_L * (1 - a0)
    derivative_upper_second[, "a0"] <- derivative_upper_second[, "a0"] + envelope_zero$L
    derivative_theta_lower <- derivative_lower_first * lower_pair$weight_first + derivative_lower_second *
        lower_pair$weight_second
    derivative_theta_upper <- derivative_upper_first * upper_pair$weight_first + derivative_upper_second *
        upper_pair$weight_second
    make_mean <- function(theta, derivative_theta) {
        value <- r0 * a0 + (1 - r0) * theta
        derivative <- derivative_theta * (1 - r0)
        derivative[, "r0"] <- derivative[, "r0"] + a0 - theta
        derivative[, "a0"] <- derivative[, "a0"] + r0
        list(value = value, derivative = derivative)
    }
    lower_mean <- make_mean(lower_pair$value, derivative_theta_lower)
    upper_mean <- make_mean(upper_pair$value, derivative_theta_upper)
    interval <- binary_gamma_compatible_interval(r0, r1, Lambda0, Lambda1)
    comparison_tolerance <- 100 * .Machine$double.eps
    lower_tolerance <- comparison_tolerance * pmax(1, abs(interval$lower))
    upper_tolerance <- comparison_tolerance * pmax(1, abs(interval$upper))
    gamma_violation_one <- gamma_one < interval$lower - lower_tolerance | gamma_one > interval$upper +
        upper_tolerance
    gamma_violation_zero <- gamma_zero < interval$lower - lower_tolerance | gamma_zero > interval$upper +
        upper_tolerance
    output <- list(lower = lower_mean$value, upper = upper_mean$value, theta_lower = lower_pair$value,
        theta_upper = upper_pair$value, gamma_one = gamma_one, gamma_zero = gamma_zero, compatible_gamma_lower = interval$lower,
        compatible_gamma_upper = interval$upper, L_one = envelope_one$L, U_one = envelope_one$U, L_zero = envelope_zero$L,
        U_zero = envelope_zero$U, candidate_lower_one = envelope_one$candidate_lower, candidate_upper_one = envelope_one$candidate_upper,
        candidate_lower_zero = envelope_zero$candidate_lower, candidate_upper_zero = envelope_zero$candidate_upper,
        gamma_violation = gamma_violation_one | gamma_violation_zero, envelope_crossing = envelope_one$L >
            envelope_one$U + comparison_tolerance * pmax(1, abs(envelope_one$U)) | envelope_zero$L >
            envelope_zero$U + comparison_tolerance * pmax(1, abs(envelope_zero$U)), theta_crossing = lower_pair$value >
            upper_pair$value + comparison_tolerance * pmax(1, abs(upper_pair$value)), epsilon = epsilon)
    if (compute_gradient) {
        output$gradient_lower <- lower_mean$derivative
        output$gradient_upper <- upper_mean$derivative
    }
    output
}

binary_arm_one_step_signal <- function(data, predictions, arm_nuisance, endpoint, t_value, direction) {
    direction <- match.arg(direction, c("lower", "upper"))
    gradient <- if (direction == "lower") {
        endpoint$gradient_lower
    }
    else {
        endpoint$gradient_upper
    }
    value <- if (direction == "lower")
        endpoint$lower
    else endpoint$upper
    observed_arm <- as.numeric(data$T_obs == t_value)
    signal <- value
    for (z_value in c(0, 1)) {
        z_probability <- if (z_value == 0)
            predictions$e
        else 1 - predictions$e
        z_indicator <- as.numeric(data$Z == z_value)
        r_name <- paste0("r", z_value)
        a_name <- paste0("a", z_value)
        r_value <- arm_nuisance[[r_name]]
        a_value <- arm_nuisance[[a_name]]
        signal <- signal + gradient[, r_name] * z_indicator/z_probability * (observed_arm - r_value)
        signal <- signal + gradient[, a_name] * z_indicator * observed_arm/(z_probability * r_value) *
            (data$Y - a_value)
    }
    signal
}

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

binary_clip_probability <- function(value, epsilon = 0.05) {
    pmin(1 - epsilon, pmax(epsilon, value))
}

binary_comparator_plugin_values <- function(predictions, Lambda) {
    e <- predictions$e
    output <- list()
    for (t_value in c(0, 1)) {
        arm <- binary_arm_nuisance_from_predictions(predictions, t_value)
        marginal_r <- e * arm$r0 + (1 - e) * arm$r1
        marginal_a <- (e * arm$r0 * arm$a0 + (1 - e) * arm$r1 * arm$a1)/marginal_r
        tan_x <- binary_standard_msm_arm(marginal_r, marginal_a, Lambda)
        tan_z0 <- binary_standard_msm_arm(arm$r0, arm$a0, Lambda)
        tan_z1 <- binary_standard_msm_arm(arm$r1, arm$a1, Lambda)
        tan_xz <- list(lower = e * tan_z0$lower + (1 - e) * tan_z1$lower, upper = e * tan_z0$upper +
            (1 - e) * tan_z1$upper)
        iv_only <- list(lower = pmax(arm$r0 * arm$a0, arm$r1 * arm$a1), upper = pmin(arm$r0 * arm$a0 +
            1 - arm$r0, arm$r1 * arm$a1 + 1 - arm$r1))
        output[[as.character(t_value)]] <- list(tan_x = tan_x, tan_xz = tan_xz, iv_only = iv_only)
    }
    make_row <- function(method, component) {
        mu1_lower <- mean(output[["1"]][[component]]$lower)
        mu1_upper <- mean(output[["1"]][[component]]$upper)
        mu0_lower <- mean(output[["0"]][[component]]$lower)
        mu0_upper <- mean(output[["0"]][[component]]$upper)
        data.frame(method = method, mu1_lower = mu1_lower, mu1_upper = mu1_upper, mu0_lower = mu0_lower,
            mu0_upper = mu0_upper, ate_lower = mu1_lower - mu0_upper, ate_upper = mu1_upper - mu0_lower,
            stringsAsFactors = FALSE)
    }
    rbind(make_row("tan_x_plugin", "tan_x"), make_row("tan_xz_plugin", "tan_xz"), make_row("iv_only_plugin",
        "iv_only"))
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

binary_iv_arm_endpoint <- function(nuisance, epsilon = 0, compute_gradient = FALSE) {
    if (!is.finite(epsilon) || epsilon < 0) {
        stop("epsilon must be one nonnegative finite number.")
    }
    nuisance <- binary_recycle_arm_nuisance(nuisance)
    r0 <- nuisance$r0
    r1 <- nuisance$r1
    a0 <- nuisance$a0
    a1 <- nuisance$a1
    n <- length(r0)
    if (any(r0 <= 0 | r0 >= 1 | r1 <= 0 | r1 >= 1 | a0 <= 0 | a0 >= 1 | a1 <= 0 | a1 >= 1)) {
        stop("Arm propensities and outcome regressions must lie strictly in (0,1).")
    }
    lower0 <- r0 * a0
    lower1 <- r1 * a1
    upper0 <- 1 - r0 * (1 - a0)
    upper1 <- 1 - r1 * (1 - a1)
    lower <- binary_smooth_max_below(lower0, lower1, epsilon)
    upper <- binary_smooth_min_above(upper0, upper1, epsilon)
    output <- list(lower = lower$value, upper = upper$value, hard_lower = pmax(lower0, lower1), hard_upper = pmin(upper0,
        upper1), epsilon = epsilon)
    if (compute_gradient) {
        gradient_names <- c("r0", "r1", "a0", "a1")
        gradient_lower <- gradient_upper <- matrix(0, nrow = n, ncol = length(gradient_names), dimnames = list(NULL,
            gradient_names))
        gradient_lower[, "r0"] <- lower$weight_first * a0
        gradient_lower[, "a0"] <- lower$weight_first * r0
        gradient_lower[, "r1"] <- lower$weight_second * a1
        gradient_lower[, "a1"] <- lower$weight_second * r1
        gradient_upper[, "r0"] <- upper$weight_first * (a0 - 1)
        gradient_upper[, "a0"] <- upper$weight_first * r0
        gradient_upper[, "r1"] <- upper$weight_second * (a1 - 1)
        gradient_upper[, "a1"] <- upper$weight_second * r1
        output$gradient_lower <- gradient_lower
        output$gradient_upper <- gradient_upper
    }
    output
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

binary_make_tz_folds <- function(data, K = 5L, seed = 1L) {
    binary_check_indicator(data$T_obs, "T_obs")
    binary_check_indicator(data$Z, "Z")
    K <- as.integer(K)
    if (K < 2L)
        stop("K must be at least two.")
    set.seed(seed)
    strata <- interaction(data$T_obs, data$Z, drop = TRUE)
    fold_id <- integer(nrow(data))
    for (stratum in levels(strata)) {
        index <- which(strata == stratum)
        if (length(index) < K) {
            stop("A (T,Z) stratum contains fewer observations than folds.")
        }
        fold_id[index] <- sample(rep(seq_len(K), length.out = length(index)))
    }
    fold_id
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

binary_project_arm_outcome_probability <- function(arm_nuisance, Lambda0, Lambda1 = Lambda0, action = c("project",
    "none"), tolerance = 1e-12) {
    action <- match.arg(action)
    arm_nuisance <- binary_recycle_arm_nuisance(arm_nuisance)
    interval <- binary_gamma_compatible_interval(arm_nuisance$r0, arm_nuisance$r1, Lambda0, Lambda1)
    lower <- pmax(interval$lower * arm_nuisance$a0, 1 - interval$upper * (1 - arm_nuisance$a0))
    upper <- pmin(interval$upper * arm_nuisance$a0, 1 - interval$lower * (1 - arm_nuisance$a0))
    if (any(lower > upper + tolerance)) {
        stop("The binary compatible interval for a1 is empty; no midpoint repair ", "was applied.")
    }
    upper <- pmax(upper, lower)
    projected_a1 <- pmin(upper, pmax(lower, arm_nuisance$a1))
    magnitude <- abs(projected_a1 - arm_nuisance$a1)
    required <- magnitude > tolerance
    stabilized <- arm_nuisance
    if (action == "project")
        stabilized$a1 <- projected_a1
    list(nuisance = stabilized, original_a1 = arm_nuisance$a1, projected_a1 = projected_a1, lower = lower,
        upper = upper, required = required, magnitude = magnitude, action = action, tolerance = tolerance)
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

binary_smooth_max_below <- function(first, second, epsilon) {
    if (epsilon == 0) {
        tie <- abs(first - second) <= 10 * .Machine$double.eps
        return(list(value = pmax(first, second), weight_first = ifelse(tie, 0.5, as.numeric(first > second)),
            weight_second = ifelse(tie, 0.5, as.numeric(second > first))))
    }
    bandwidth <- epsilon/log(2)
    maximum <- pmax(first, second)
    exponential_first <- exp((first - maximum)/bandwidth)
    exponential_second <- exp((second - maximum)/bandwidth)
    denominator <- exponential_first + exponential_second
    list(value = maximum + bandwidth * log(denominator/2), weight_first = exponential_first/denominator,
        weight_second = exponential_second/denominator)
}

binary_smooth_min_above <- function(first, second, epsilon) {
    if (epsilon == 0) {
        tie <- abs(first - second) <= 10 * .Machine$double.eps
        return(list(value = pmin(first, second), weight_first = ifelse(tie, 0.5, as.numeric(first < second)),
            weight_second = ifelse(tie, 0.5, as.numeric(second < first))))
    }
    bandwidth <- epsilon/log(2)
    minimum <- pmin(first, second)
    exponential_first <- exp(-(first - minimum)/bandwidth)
    exponential_second <- exp(-(second - minimum)/bandwidth)
    denominator <- exponential_first + exponential_second
    list(value = minimum - bandwidth * log(denominator/2), weight_first = exponential_first/denominator,
        weight_second = exponential_second/denominator)
}

binary_standard_msm_arm <- function(r, a, Lambda) {
    theta_lower <- pmax(a/Lambda, 1 - (1 - a) * Lambda)
    theta_upper <- pmin(a * Lambda, 1 - (1 - a)/Lambda)
    list(lower = r * a + (1 - r) * theta_lower, upper = r * a + (1 - r) * theta_upper)
}

clip_positive <- function(x, floor = 1e-06, cap = Inf) {
    x[!is.finite(x)] <- floor
    x <- pmax(floor, x)
    if (is.finite(cap)) {
        x <- pmin(cap, x)
    }
    x
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

observed_binary_probabilities <- function(x, treatment, specification) {
    treatment <- as.integer(treatment)
    if (!treatment %in% 0:1)
        stop("treatment must be zero or one.")
    features <- observed_dgp_features(x, specification)
    baseline <- stats::plogis(specification$binary_outcome_intercept + specification$binary_observed_treatment_log_odds_shift *
        treatment + features$outcome)
    comparison <- stats::plogis(stats::qlogis(baseline) + ifelse(treatment == 1L, 1, -1) * specification$binary_tilt_scale *
        features$tilt_shape)
    list(a0 = baseline, a1 = comparison)
}

observed_dgp_features <- function(x, specification) {
    x <- observed_dgp_validate_x(x, specification)
    if (isTRUE(specification$flat_compatibility_design)) {
        x1 <- x[, 1L]
        return(list(instrument = 0.2 * x1, treatment = rep(0, length(x1)), outcome = 0.65 * x1, tilt_shape = rep(1,
            length(x1))))
    }
    if (specification$dimension == 1L) {
        x1 <- x[, 1L]
        return(list(instrument = 0.35 * x1, treatment = 0.7 * x1, outcome = 0.65 * x1 + 0.15 * sin(pi *
            x1/2), tilt_shape = 0.75 + 0.25 * sin(pi * x1/2)^2))
    }
    x1 <- x[, 1L]
    x2 <- x[, 2L]
    x3 <- x[, 3L]
    x4 <- x[, 4L]
    x5 <- x[, 5L]
    x5_treatment_coefficient <- if (is.null(specification$x5_treatment_coefficient))
        0
    else specification$x5_treatment_coefficient
    treatment_feature <- 0.42 * x1 - 0.3 * x2 + 0.22 * x3 + 0.16 * x1 * x2 + 0.14 * sin(x4) + x5_treatment_coefficient *
        x5
    tilt_shape <- 0.75 + 0.25 * sin((x1 + x2 - x3)/sqrt(3))^2
    if (isTRUE(specification$flat_compatibility_treatment)) {
        treatment_feature <- rep(0, length(x1))
    }
    treatment_feature_scale <- if (is.null(specification$treatment_feature_scale))
        1
    else specification$treatment_feature_scale
    treatment_feature <- treatment_feature_scale * treatment_feature
    if (isTRUE(specification$constant_compatibility_tilt_shape)) {
        tilt_shape <- rep(1, length(x1))
    }
    list(instrument = 0.28 * x1 - 0.2 * x2 + 0.12 * x3, treatment = treatment_feature, outcome = 0.45 *
        x1 + 0.3 * x2 - 0.22 * x3 + 0.16 * x4 * x5 + 0.12 * (x1^2 - 1), tilt_shape = tilt_shape)
}

observed_dgp_validate_x <- function(x, specification) {
    columns <- observed_dgp_x_columns(specification$dimension)
    if (is.data.frame(x)) {
        missing <- setdiff(columns, names(x))
        if (length(missing)) {
            stop("Missing observed-DGP covariates: ", paste(missing, collapse = ", "))
        }
        matrix_value <- as.matrix(x[, columns, drop = FALSE])
    }
    else {
        matrix_value <- as.matrix(x)
        if (specification$dimension == 1L && ncol(matrix_value) != 1L) {
            matrix_value <- matrix(as.numeric(x), ncol = 1L)
        }
    }
    storage.mode(matrix_value) <- "double"
    if (ncol(matrix_value) != specification$dimension || any(!is.finite(matrix_value))) {
        stop("x has the wrong dimension or contains non-finite values.")
    }
    colnames(matrix_value) <- columns
    matrix_value
}

observed_dgp_x_columns <- function(dimension) {
    if (dimension == 1L)
        "X"
    else paste0("X", seq_len(dimension))
}

