# Internal deterministic check of population identities, not a user estimator example.
# Population identities only: no fitting, standard errors, or tests.
if (!file.exists('R/binary_continuum_band.R'))
  stop('Run this script from the repository root.')
source('R/binary_continuum_band.R')
evaluate <- function(Lambda) {
  treated <- binary_flat_arm_bounds(.30, .70, .45, .60, Lambda)
  untreated <- binary_flat_arm_bounds(.70, .30, .34, .40, Lambda)
  compatible <- treated$lower <= treated$upper &&
    untreated$lower <= untreated$upper
  data.frame(Lambda = Lambda, compatible = compatible,
    mu1_lower = treated$lower, mu1_upper = treated$upper,
    mu0_lower = untreated$lower, mu0_upper = untreated$upper,
    ate_lower = if (compatible) treated$lower - untreated$upper else NA_real_,
    ate_upper = if (compatible) treated$upper - untreated$lower else NA_real_)
}
answer <- do.call(rbind, lapply(c(1, 1.375, 1.6, 2, 3), evaluate))
stopifnot(!answer$compatible[1])
at16 <- answer[answer$Lambda == 1.6, ]
stopifnot(at16$compatible, abs(at16$ate_lower - .1313) < 1e-10,
          abs(at16$ate_upper - .292625) < 1e-10)
print(answer, row.names = FALSE, digits = 6)
cat('\nCrossed mean intervals indicate population incompatibility; no repair is applied.\n')
cat('ATE bounds are omitted for incompatible specifications. These are not confidence intervals.\n')
