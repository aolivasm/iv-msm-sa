# Numerical pathwise-derivative helpers for the dual-MSM EIF.
#
# These functions are diagnostics. Source R/sharp_bounds_plugin.R before this
# file so the envelope and EIF primitive constructors are available.

normalize_audit_density <- function(density, weights) {
  mass <- sum(density * weights)
  if (!is.finite(mass) || mass <= 0) {
    stop("The audit density has invalid mass.")
  }
  density / mass
}

audit_density_tilt <- function(density, score, epsilon, weights) {
  normalize_audit_density(
    density * exp(epsilon * score),
    weights
  )
}

audit_nuisance_from_arm_probability <- function(
    t_value,
    r0,
    r1,
    n
) {
  if (!t_value %in% c(0, 1)) stop("t_value must be zero or one.")
  if (any(c(r0, r1) <= 0) || any(c(r0, r1) >= 1)) {
    stop("Arm probabilities must lie strictly between zero and one.")
  }
  list(
    e = rep(0.5, n),
    pi0 = rep(if (t_value == 1) r0 else 1 - r0, n),
    pi1 = rep(if (t_value == 1) r1 else 1 - r1, n)
  )
}

evaluate_eif_arm_pathwise_endpoint <- function(
    y,
    weights,
    f0,
    f1,
    r0,
    r1,
    t_value,
    Lambda0 = 5,
    Lambda1 = 5,
    direction = c("lower", "upper"),
    return_gradient = FALSE
) {
  direction <- match.arg(direction)
  nuisance <- audit_nuisance_from_arm_probability(
    t_value = t_value,
    r0 = r0,
    r1 = r1,
    n = length(y)
  )
  gamma <- f1 / f0
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = gamma,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    gamma_compatibility_action = "stop"
  )

  probability_mass0 <- f0 * weights
  gap_mass <- probability_mass0 * envelope$D
  target_gap_mass <- if (direction == "lower") {
    1 - sum(probability_mass0 * envelope$L)
  } else {
    sum(probability_mass0 * envelope$U) - 1
  }
  target_gap_mass <- min(sum(gap_mass), max(0, target_gap_mass))
  cumulative <- cumsum(gap_mass)
  boundary <- which(cumulative >= target_gap_mass)[1L]
  if (is.na(boundary)) boundary <- length(y)
  before <- if (boundary == 1L) 0 else cumulative[boundary - 1L]
  fraction <- if (gap_mass[boundary] > 0) {
    (target_gap_mass - before) / gap_mass[boundary]
  } else {
    0
  }
  fraction <- min(1, max(0, fraction))
  low_tail <- as.numeric(seq_along(y) < boundary)
  low_tail[boundary] <- fraction
  allocate_upper <- if (direction == "lower") low_tail else 1 - low_tail
  lambda <- allocate_upper * envelope$U +
    (1 - allocate_upper) * envelope$L
  q <- y[boundary]
  outcome_mean0 <- sum(probability_mass0 * y)
  theta <- sum(probability_mass0 * lambda * y)
  value <- r0 * outcome_mean0 + (1 - r0) * theta

  if (!return_gradient) return(value)

  primitive <- compute_eif_row_primitives(
    data = data.frame(Y = y),
    nuisance = nuisance,
    envelope = envelope,
    t_value = t_value,
    Lambda0 = Lambda0,
    Lambda1 = Lambda1,
    q_lower = rep(q, length(y)),
    q_upper = rep(q, length(y))
  )[[direction]]

  list(
    value = value,
    r0 = outcome_mean0 - theta +
      (1 - r0) * sum(probability_mass0 * primitive$D_r0),
    r1 = (1 - r0) * sum(probability_mass0 * primitive$D_r1),
    f0 = primitive$H0,
    f1 = primitive$H1,
    q = q,
    lambda_normalization_error =
      sum(probability_mass0 * lambda) - 1
  )
}

run_eif_treatment_arm_pathwise_audit <- function(
    y_nodes = 40001L,
    epsilon = 1e-5,
    Lambda0 = 5,
    Lambda1 = 5
) {
  required <- c(
    "compute_projected_envelopes",
    "compute_eif_row_primitives"
  )
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing) > 0L) {
    stop(
      "Source R/sharp_bounds_plugin.R before running the EIF audit. Missing: ",
      paste(missing, collapse = ", ")
    )
  }
  y_nodes <- as.integer(y_nodes)
  if (y_nodes < 2001L) stop("y_nodes must be at least 2001.")
  if (!is.finite(epsilon) || epsilon <= 0) {
    stop("epsilon must be positive.")
  }

  y <- seq(-2, 2, length.out = y_nodes)
  step <- y[2L] - y[1L]
  weights <- rep(step, y_nodes)
  weights[c(1L, y_nodes)] <- step / 2

  rows <- list()
  counter <- 1L
  for (t_value in c(0, 1)) {
    if (t_value == 0) {
      f0 <- stats::dnorm(y, mean = -0.20, sd = 0.90)
      f1 <- stats::dnorm(y, mean = 0.05, sd = 0.95)
      r0 <- 0.55
      r1 <- 0.35
    } else {
      f0 <- stats::dnorm(y, mean = -0.15, sd = 0.90)
      f1 <- stats::dnorm(y, mean = 0.20, sd = 0.95)
      r0 <- 0.45
      r1 <- 0.65
    }
    f0 <- normalize_audit_density(f0, weights)
    f1 <- normalize_audit_density(f1, weights)
    score0 <- y + 0.3 * y^2
    score0 <- score0 - sum(weights * f0 * score0)
    score1 <- sin(y) - 0.2 * y
    score1 <- score1 - sum(weights * f1 * score1)
    pi_sign <- if (t_value == 1) 1 else -1

    for (direction in c("lower", "upper")) {
      gradient <- evaluate_eif_arm_pathwise_endpoint(
        y = y,
        weights = weights,
        f0 = f0,
        f1 = f1,
        r0 = r0,
        r1 = r1,
        t_value = t_value,
        Lambda0 = Lambda0,
        Lambda1 = Lambda1,
        direction = direction,
        return_gradient = TRUE
      )

      perturbations <- list(
        arm_probability_z0 = list(
          plus = list(f0 = f0, f1 = f1, r0 = r0 + epsilon, r1 = r1),
          minus = list(f0 = f0, f1 = f1, r0 = r0 - epsilon, r1 = r1),
          analytic = gradient$r0
        ),
        arm_probability_z1 = list(
          plus = list(f0 = f0, f1 = f1, r0 = r0, r1 = r1 + epsilon),
          minus = list(f0 = f0, f1 = f1, r0 = r0, r1 = r1 - epsilon),
          analytic = gradient$r1
        ),
        treatment_probability_z0 = list(
          plus = list(
            f0 = f0, f1 = f1,
            r0 = r0 + pi_sign * epsilon, r1 = r1
          ),
          minus = list(
            f0 = f0, f1 = f1,
            r0 = r0 - pi_sign * epsilon, r1 = r1
          ),
          analytic = pi_sign * gradient$r0
        ),
        treatment_probability_z1 = list(
          plus = list(
            f0 = f0, f1 = f1,
            r0 = r0, r1 = r1 + pi_sign * epsilon
          ),
          minus = list(
            f0 = f0, f1 = f1,
            r0 = r0, r1 = r1 - pi_sign * epsilon
          ),
          analytic = pi_sign * gradient$r1
        ),
        density_z0 = list(
          plus = list(
            f0 = audit_density_tilt(f0, score0, epsilon, weights),
            f1 = f1, r0 = r0, r1 = r1
          ),
          minus = list(
            f0 = audit_density_tilt(f0, score0, -epsilon, weights),
            f1 = f1, r0 = r0, r1 = r1
          ),
          analytic = sum(weights * f0 * gradient$f0 * score0)
        ),
        density_z1 = list(
          plus = list(
            f0 = f0,
            f1 = audit_density_tilt(f1, score1, epsilon, weights),
            r0 = r0, r1 = r1
          ),
          minus = list(
            f0 = f0,
            f1 = audit_density_tilt(f1, score1, -epsilon, weights),
            r0 = r0, r1 = r1
          ),
          analytic = sum(weights * f1 * gradient$f1 * score1)
        )
      )

      for (submodel in names(perturbations)) {
        perturbation <- perturbations[[submodel]]
        common <- list(
          y = y,
          weights = weights,
          t_value = t_value,
          Lambda0 = Lambda0,
          Lambda1 = Lambda1,
          direction = direction
        )
        numerical <- (
          do.call(
            evaluate_eif_arm_pathwise_endpoint,
            c(common, perturbation$plus)
          ) -
            do.call(
              evaluate_eif_arm_pathwise_endpoint,
              c(common, perturbation$minus)
            )
        ) / (2 * epsilon)
        rows[[counter]] <- data.frame(
          treatment = t_value,
          direction = direction,
          submodel = submodel,
          analytic = perturbation$analytic,
          numerical = numerical,
          difference = perturbation$analytic - numerical,
          absolute_difference = abs(perturbation$analytic - numerical),
          lambda_normalization_error =
            gradient$lambda_normalization_error,
          threshold = gradient$q,
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }
  }

  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}
