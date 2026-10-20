estimate_unconfounded_aipw <- function(data, nuisance_fit, confidence_level) {
    predictions <- nuisance_fit$predictions
    signals <- list()
    for (treatment in c(0, 1)) {
        treatment_probability <- if (treatment == 1) {
            ifelse(data$Z == 0, predictions$pi0, predictions$pi1)
        }
        else {
            ifelse(data$Z == 0, 1 - predictions$pi0, 1 - predictions$pi1)
        }
        conditional_mean <- if (treatment == 1) {
            ifelse(data$Z == 0, predictions$a_t1_z0, predictions$a_t1_z1)
        }
        else {
            ifelse(data$Z == 0, predictions$a_t0_z0, predictions$a_t0_z1)
        }
        signals[[as.character(treatment)]] <- conditional_mean + as.numeric(data$T_obs == treatment)/treatment_probability *
            (data$Y - conditional_mean)
    }
    ate_signal <- signals[["1"]] - signals[["0"]]
    estimates <- c(mu1 = mean(signals[["1"]]), mu0 = mean(signals[["0"]]), ate = mean(ate_signal))
    standard_errors <- c(mu1 = stats::sd(signals[["1"]])/sqrt(nrow(data)), mu0 = stats::sd(signals[["0"]])/sqrt(nrow(data)),
        ate = stats::sd(ate_signal)/sqrt(nrow(data)))
    critical <- stats::qnorm(1 - (1 - confidence_level)/2)
    data.frame(estimand = names(estimates), estimate = unname(estimates), standard_error = unname(standard_errors),
        confidence_level = confidence_level, ci_lower = unname(estimates - critical * standard_errors),
        ci_upper = unname(estimates + critical * standard_errors), stringsAsFactors = FALSE)
}
