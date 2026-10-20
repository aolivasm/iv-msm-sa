# Oracle observed-data experiment for pointwise envelope nonregularity.
#
# The law has six exact envelope-tie covariate strata and one regular
# first-stage stratum.  In every stratum,
#
#   P(T = 1 | X, Z = 0) = 0.10.
#
# In the tie strata P(T = 1 | X, Z = 1) is also 0.10.  In the regular
# stratum it is 0.90.  With 90% population mass on the tie strata and 10% on
# the regular stratum, the marginal first stage is 0.08 while every treatment
# probability and its complement is at least 0.10.
#
# Outcomes follow
#
#   Y = treatment_effect * T + E,   E ~ Uniform(0, 1),
#
# independently of Z conditional on (T, X), so gamma_t(Y, X) = 1.  The
# experiment holds the instrument propensity, covariate probabilities, and
# outcome law at their known population values and estimates only pi_0(X) and
# pi_1(X) in the exact-tie strata.  This is the analogue of an oracle
# exceptional-law experiment: it isolates the max/min envelope derivative
# instead of allowing unrelated nuisance noise to conceal it.
#
# The primary target is either endpoint of the mu1 sensitivity interval.  The
# hard and shifted-log-mean-exp envelopes agree exactly at a population tie.
# Consequently, smoothing changes neither endpoint in the 90% tie portion of
# the population.  The remaining target displacement is reported explicitly,
# and the smoothed confidence interval is translated to the sharp target
# without changing its length.  Coverage after this translation cannot be
# attributed to widening the confidence interval.

make_oracle_envelope_tie_spec <- function(
    sensitivity_lambda = 3,
    tie_strata = 6L,
    tie_mass = 0.9,
    probability_z1 = 0.5,
    treatment_probability_z0 = 0.1,
    regular_treatment_probability_z1 = 0.9,
    treatment_effect = 0.75
) {
  sensitivity_lambda <- as.numeric(sensitivity_lambda)
  tie_strata <- as.integer(tie_strata)
  numeric_values <- c(
    sensitivity_lambda, tie_mass, probability_z1,
    treatment_probability_z0,
    regular_treatment_probability_z1,
    treatment_effect
  )
  if (
    length(tie_strata) != 1L || is.na(tie_strata) || tie_strata < 1L ||
      any(!is.finite(numeric_values)) || sensitivity_lambda <= 1 ||
      tie_mass <= 0 || tie_mass >= 1 ||
      probability_z1 <= 0 || probability_z1 >= 1 ||
      treatment_probability_z0 <= 0 ||
      treatment_probability_z0 >= 0.5 ||
      regular_treatment_probability_z1 <= treatment_probability_z0 ||
      regular_treatment_probability_z1 >= 1
  ) {
    stop("Invalid oracle envelope-tie specification.")
  }

  x_probabilities <- c(
    rep(tie_mass / tie_strata, tie_strata),
    1 - tie_mass
  )
  pi0 <- rep(treatment_probability_z0, tie_strata + 1L)
  pi1 <- c(
    rep(treatment_probability_z0, tie_strata),
    regular_treatment_probability_z1
  )
  overlap_constant <- min(
    probability_z1, 1 - probability_z1,
    pi0, 1 - pi0, pi1, 1 - pi1
  )

  specification <- list(
    name = "oracle_observed_data_exact_envelope_ties",
    sensitivity_lambda = sensitivity_lambda,
    tie_strata = tie_strata,
    regular_stratum = tie_strata + 1L,
    x_values = seq_len(tie_strata + 1L),
    x_probabilities = x_probabilities,
    tie_mass = tie_mass,
    probability_z1 = probability_z1,
    pi0 = pi0,
    pi1 = pi1,
    treatment_effect = treatment_effect,
    overlap_constant = overlap_constant,
    marginal_treatment_probability_z0 = sum(x_probabilities * pi0),
    marginal_treatment_probability_z1 = sum(x_probabilities * pi1),
    marginal_first_stage = sum(x_probabilities * (pi1 - pi0)),
    gamma_t = 1
  )
  class(specification) <- c(
    "oracle_envelope_tie_specification", "list"
  )
  specification
}

oracle_envelope_tie_x_index <- function(x, specification) {
  x <- as.integer(x)
  index <- match(x, specification$x_values)
  if (anyNA(index)) {
    stop("X contains a value outside the oracle envelope-tie support.")
  }
  index
}

generate_oracle_envelope_tie_data <- function(
    n,
    specification = make_oracle_envelope_tie_spec(),
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  x <- sample(
    specification$x_values,
    size = n,
    replace = TRUE,
    prob = specification$x_probabilities
  )
  index <- oracle_envelope_tie_x_index(x, specification)
  z <- stats::rbinom(n, 1L, specification$probability_z1)
  pi0 <- specification$pi0[index]
  pi1 <- specification$pi1[index]
  probability_t <- ifelse(z == 1L, pi1, pi0)
  treatment <- stats::rbinom(n, 1L, probability_t)
  residual <- stats::runif(n)
  y <- specification$treatment_effect * treatment + residual

  output <- data.frame(
    Y = y,
    T_obs = treatment,
    Z = z,
    X = x,
    e_true = rep(1 - specification$probability_z1, n),
    pi0_true = pi0,
    pi1_true = pi1,
    gamma0_oracle = rep(1, n),
    gamma1_oracle = rep(1, n),
    stringsAsFactors = FALSE
  )
  for (level in seq_len(specification$tie_strata)) {
    output[[paste0("X_tie_", level)]] <- as.numeric(x == level)
  }
  output$X_regular <- as.numeric(x == specification$regular_stratum)
  output
}

make_oracle_envelope_tie_nuisance_provider <- function(
    specification = make_oracle_envelope_tie_spec()
) {
  cache <- new.env(parent = emptyenv())

  function(
      fold,
      train,
      x_cols,
      seed,
      positivity_epsilon = 0.05,
      ...
  ) {
    key <- paste0("fold_", as.integer(fold))
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    required <- c("X", "Z", "T_obs")
    if (!all(required %in% names(train))) {
      stop("The saturated nuisance provider requires X, Z, and T_obs.")
    }
    levels <- specification$x_values

    cell_probability <- function(response, subset = rep(TRUE, nrow(train))) {
      raw <- count <- numeric(length(levels))
      for (index in seq_along(levels)) {
        local <- subset & train$X == levels[index]
        count[index] <- sum(local)
        if (count[index] == 0L) {
          stop("An outer-fold saturated nuisance cell is empty.")
        }
        raw[index] <- mean(response[local])
      }
      stabilized <- pmin(
        pmax(raw, positivity_epsilon),
        1 - positivity_epsilon
      )
      list(
        raw = raw,
        stabilized = stabilized,
        count = count,
        stabilization_rate = mean(stabilized != raw)
      )
    }

    e_fit <- cell_probability(as.numeric(train$Z == 0))
    pi0_fit <- cell_probability(train$T_obs, train$Z == 0)
    pi1_fit <- cell_probability(train$T_obs, train$Z == 1)

    predict_nuisance <- function(newdata) {
      index <- match(as.integer(newdata$X), levels)
      if (anyNA(index)) {
        stop("Prediction data contain an unsupported X value.")
      }
      list(
        e = e_fit$stabilized[index],
        pi0 = pi0_fit$stabilized[index],
        pi1 = pi1_fit$stabilized[index]
      )
    }

    output <- list(
      predict = predict_nuisance,
      learner = "saturated_cell",
      diagnostics = data.frame(
        component = c("e", "pi0", "pi1"),
        learner = "saturated_cell",
        converged = TRUE,
        model_rank = rep(length(levels), 3L),
        effective_df = rep(length(levels), 3L),
        n_training = c(
          nrow(train), sum(train$Z == 0), sum(train$Z == 1)
        ),
        minimum_cell_count = c(
          min(e_fit$count), min(pi0_fit$count), min(pi1_fit$count)
        ),
        prediction_stabilization_rate = c(
          e_fit$stabilization_rate,
          pi0_fit$stabilization_rate,
          pi1_fit$stabilization_rate
        ),
        stringsAsFactors = FALSE
      ),
      fitted_probabilities = list(
        e = e_fit,
        pi0 = pi0_fit,
        pi1 = pi1_fit
      )
    )
    assign(key, output, envir = cache)
    output
  }
}

oracle_envelope_tie_evaluate_mu1 <- function(
    pi0,
    pi1,
    endpoint = c("lower", "upper"),
    sensitivity_lambda = 3,
    envelope_smoothing_epsilon = 0,
    treatment_effect = 0.75
) {
  endpoint <- match.arg(endpoint)
  pi0 <- as.numeric(pi0)
  pi1 <- as.numeric(pi1)
  pi1 <- rep_len(pi1, length(pi0))
  Lambda <- as.numeric(sensitivity_lambda)
  epsilon <- as.numeric(envelope_smoothing_epsilon)
  values <- c(pi0, pi1, Lambda, epsilon, treatment_effect)
  if (
    any(!is.finite(values)) || any(pi0 <= 0) || any(pi0 >= 1) ||
      any(pi1 <= 0) || any(pi1 >= 1) || Lambda <= 1 ||
      length(epsilon) != 1L || epsilon < 0 || epsilon >= 1 / Lambda
  ) {
    stop("Invalid inputs to oracle_envelope_tie_evaluate_mu1().")
  }

  h_lower <- pi1 + (1 - pi1) / Lambda
  h_upper <- pi1 + (1 - pi1) * Lambda
  candidate_lower <- (h_lower - pi0) / (1 - pi0)
  candidate_upper <- (h_upper - pi0) / (1 - pi0)

  derivative_candidate_lower_pi0 <-
    (h_lower - 1) / (1 - pi0)^2
  derivative_candidate_lower_pi1 <-
    (1 - 1 / Lambda) / (1 - pi0)
  derivative_candidate_upper_pi0 <-
    (h_upper - 1) / (1 - pi0)^2
  derivative_candidate_upper_pi1 <-
    (1 - Lambda) / (1 - pi0)

  if (epsilon == 0) {
    lower_weight <- as.numeric(candidate_lower > 1 / Lambda)
    upper_weight <- as.numeric(candidate_upper < Lambda)
    L <- pmax(1 / Lambda, candidate_lower)
    U <- pmin(Lambda, candidate_upper)
  } else {
    bandwidth <- epsilon / log(2)
    lower_weight <- stats::plogis(
      (candidate_lower - 1 / Lambda) / bandwidth
    )
    upper_weight <- stats::plogis(
      (Lambda - candidate_upper) / bandwidth
    )

    lower_maximum <- pmax(1 / Lambda, candidate_lower)
    L <- lower_maximum + bandwidth * log(
      (
        exp((1 / Lambda - lower_maximum) / bandwidth) +
          exp((candidate_lower - lower_maximum) / bandwidth)
      ) / 2
    )
    upper_minimum <- pmin(Lambda, candidate_upper)
    U <- upper_minimum - bandwidth * log(
      (
        exp(-(Lambda - upper_minimum) / bandwidth) +
          exp(-(candidate_upper - upper_minimum) / bandwidth)
      ) / 2
    )
  }

  derivative_L_pi0 <- lower_weight * derivative_candidate_lower_pi0
  derivative_L_pi1 <- lower_weight * derivative_candidate_lower_pi1
  derivative_U_pi0 <- upper_weight * derivative_candidate_upper_pi0
  derivative_U_pi1 <- upper_weight * derivative_candidate_upper_pi1

  if (any(!is.finite(L)) || any(!is.finite(U)) || any(L >= U)) {
    stop("The requested nuisance values produce invalid envelopes.")
  }

  outcome_mean <- treatment_effect + 0.5
  if (endpoint == "lower") {
    quantile_residual <- (1 - L) / (U - L)
    lower_partial_moment <-
      treatment_effect * quantile_residual + quantile_residual^2 / 2
    tilted_mean <-
      L * outcome_mean + (U - L) * lower_partial_moment

    derivative_tilted_mean <- function(derivative_L, derivative_U) {
      derivative_quantile <- (
        -derivative_L * (U - L) -
          (1 - L) * (derivative_U - derivative_L)
      ) / (U - L)^2
      derivative_L * outcome_mean +
        (derivative_U - derivative_L) * lower_partial_moment +
        (U - L) *
          (treatment_effect + quantile_residual) * derivative_quantile
    }
  } else {
    quantile_residual <- (U - 1) / (U - L)
    lower_partial_moment <-
      treatment_effect * quantile_residual + quantile_residual^2 / 2
    tilted_mean <-
      U * outcome_mean + (L - U) * lower_partial_moment

    derivative_tilted_mean <- function(derivative_L, derivative_U) {
      derivative_quantile <- (
        derivative_U * (U - L) -
          (U - 1) * (derivative_U - derivative_L)
      ) / (U - L)^2
      derivative_U * outcome_mean +
        (derivative_L - derivative_U) * lower_partial_moment +
        (L - U) *
          (treatment_effect + quantile_residual) * derivative_quantile
    }
  }

  derivative_tilted_pi0 <- derivative_tilted_mean(
    derivative_L_pi0, derivative_U_pi0
  )
  derivative_tilted_pi1 <- derivative_tilted_mean(
    derivative_L_pi1, derivative_U_pi1
  )

  estimate <- pi0 * outcome_mean + (1 - pi0) * tilted_mean
  derivative_pi0 <- outcome_mean - tilted_mean +
    (1 - pi0) * derivative_tilted_pi0
  derivative_pi1 <- (1 - pi0) * derivative_tilted_pi1

  list(
    estimate = estimate,
    derivative_pi0 = derivative_pi0,
    derivative_pi1 = derivative_pi1,
    L = L,
    U = U,
    candidate_lower = candidate_lower,
    candidate_upper = candidate_upper,
    lower_derivative_weight = lower_weight,
    upper_derivative_weight = upper_weight,
    quantile = treatment_effect + quantile_residual,
    tilted_mean = tilted_mean,
    endpoint = endpoint,
    envelope_smoothing_epsilon = epsilon
  )
}

oracle_envelope_tie_population <- function(
    specification = make_oracle_envelope_tie_spec(),
    envelope_smoothing_epsilon = 0
) {
  endpoints <- lapply(c("lower", "upper"), function(endpoint) {
    local <- oracle_envelope_tie_evaluate_mu1(
      pi0 = specification$pi0,
      pi1 = specification$pi1,
      endpoint = endpoint,
      sensitivity_lambda = specification$sensitivity_lambda,
      envelope_smoothing_epsilon = envelope_smoothing_epsilon,
      treatment_effect = specification$treatment_effect
    )
    sum(specification$x_probabilities * local$estimate)
  })
  names(endpoints) <- c("mu1_lower", "mu1_upper")
  data.frame(
    envelope_smoothing_epsilon = envelope_smoothing_epsilon,
    mu1_lower = endpoints$mu1_lower,
    mu1_upper = endpoints$mu1_upper,
    mu1_width = endpoints$mu1_upper - endpoints$mu1_lower,
    mu1_truth = specification$treatment_effect + 0.5,
    stringsAsFactors = FALSE
  )
}

oracle_envelope_tie_population_grid <- function(
    specification = make_oracle_envelope_tie_spec(),
    epsilon_grid = c(0, 0.005, 0.01, 0.05, 0.1)
) {
  epsilon_grid <- sort(unique(as.numeric(epsilon_grid)))
  if (
    length(epsilon_grid) < 2L || any(!is.finite(epsilon_grid)) ||
      any(epsilon_grid < 0) ||
      any(epsilon_grid >= 1 / specification$sensitivity_lambda) ||
      !any(epsilon_grid == 0)
  ) {
    stop("epsilon_grid must include zero and valid nonnegative values.")
  }
  output <- do.call(rbind, lapply(epsilon_grid, function(epsilon) {
    oracle_envelope_tie_population(
      specification,
      envelope_smoothing_epsilon = epsilon
    )
  }))
  hard <- output[output$envelope_smoothing_epsilon == 0, , drop = FALSE]
  output$lower_target_shift_from_hard <- output$mu1_lower - hard$mu1_lower
  output$upper_target_shift_from_hard <- output$mu1_upper - hard$mu1_upper
  output$width_increase_from_hard <- output$mu1_width - hard$mu1_width
  rownames(output) <- NULL
  output
}

oracle_envelope_tie_draw_sufficient_statistics <- function(
    n,
    replications,
    specification,
    seed
) {
  n <- as.integer(n)
  replications <- as.integer(replications)
  if (
    length(n) != 1L || is.na(n) || n < 200L ||
      length(replications) != 1L || is.na(replications) || replications < 1L
  ) {
    stop("Invalid sufficient-statistic simulation dimensions.")
  }
  set.seed(as.integer(seed))
  joint_probability <- as.vector(outer(
    specification$x_probabilities,
    c(1 - specification$probability_z1, specification$probability_z1)
  ))
  cell_count <- t(stats::rmultinom(
    replications,
    size = n,
    prob = joint_probability
  ))
  tie <- seq_len(specification$tie_strata)
  z1_offset <- length(specification$x_values)
  n_z0 <- cell_count[, tie, drop = FALSE]
  n_z1 <- cell_count[, z1_offset + tie, drop = FALSE]
  if (any(n_z0 == 0L) || any(n_z1 == 0L)) {
    stop(
      "At least one simulated tie-stratum instrument cell was empty; increase n."
    )
  }
  p_tie <- specification$pi0[1L]
  treated_z0 <- matrix(
    stats::rbinom(length(n_z0), size = c(n_z0), prob = p_tie),
    nrow = replications,
    ncol = specification$tie_strata
  )
  treated_z1 <- matrix(
    stats::rbinom(length(n_z1), size = c(n_z1), prob = p_tie),
    nrow = replications,
    ncol = specification$tie_strata
  )
  list(
    n_z0 = n_z0,
    n_z1 = n_z1,
    treated_z0 = treated_z0,
    treated_z1 = treated_z1
  )
}

run_oracle_envelope_tie_simulation <- function(
    sample_sizes = c(1500L, 3000L, 6000L),
    replications = 20000L,
    epsilon_grid = c(0, 0.005, 0.01, 0.05, 0.1),
    specification = make_oracle_envelope_tie_spec(),
    seed_start = 94026L,
    batch_size = 10000L,
    confidence_level = 0.95,
    numerical_probability_clip = 1e-8
) {
  sample_sizes <- as.integer(sample_sizes)
  replications <- as.integer(replications)
  batch_size <- as.integer(batch_size)
  epsilon_grid <- sort(unique(as.numeric(epsilon_grid)))
  if (
    length(sample_sizes) < 1L || anyNA(sample_sizes) ||
      any(sample_sizes < 200L) || replications < 1L || batch_size < 1L ||
      any(!is.finite(epsilon_grid)) || !any(epsilon_grid == 0) ||
      any(epsilon_grid < 0) ||
      any(epsilon_grid >= 1 / specification$sensitivity_lambda) ||
      !is.finite(confidence_level) || confidence_level <= 0 ||
      confidence_level >= 1 ||
      !is.finite(numerical_probability_clip) ||
      numerical_probability_clip <= 0 ||
      numerical_probability_clip >= specification$overlap_constant
  ) {
    stop("Invalid oracle envelope-tie simulation settings.")
  }

  population <- oracle_envelope_tie_population_grid(
    specification, epsilon_grid
  )
  hard_population <- population[
    population$envelope_smoothing_epsilon == 0, , drop = FALSE
  ]
  critical_value <- stats::qnorm(1 - (1 - confidence_level) / 2)
  endpoint_names <- c("mu1_lower", "mu1_upper")
  detail <- list()
  detail_counter <- 1L

  for (sample_index in seq_along(sample_sizes)) {
    n <- sample_sizes[sample_index]
    completed <- 0L
    batch_index <- 0L
    while (completed < replications) {
      local_replications <- min(batch_size, replications - completed)
      batch_index <- batch_index + 1L
      sufficient <- oracle_envelope_tie_draw_sufficient_statistics(
        n = n,
        replications = local_replications,
        specification = specification,
        seed = seed_start + 100000L * sample_index + batch_index - 1L
      )
      raw_pi0 <- sufficient$treated_z0 / sufficient$n_z0
      raw_pi1 <- sufficient$treated_z1 / sufficient$n_z1
      clipped_pi0 <- pmin(
        pmax(raw_pi0, numerical_probability_clip),
        1 - numerical_probability_clip
      )
      clipped_pi1 <- pmin(
        pmax(raw_pi1, numerical_probability_clip),
        1 - numerical_probability_clip
      )
      stabilization_rate <- rowMeans(
        clipped_pi0 != raw_pi0 | clipped_pi1 != raw_pi1
      )
      tie_branch_z1_rate <- rowMeans(clipped_pi1 > clipped_pi0)
      empirical_exact_tie_rate <- rowMeans(clipped_pi1 == clipped_pi0)
      variance_pi0 <-
        clipped_pi0 * (1 - clipped_pi0) / sufficient$n_z0
      variance_pi1 <-
        clipped_pi1 * (1 - clipped_pi1) / sufficient$n_z1
      tie_weight <- specification$tie_mass / specification$tie_strata
      regular_index <- specification$regular_stratum

      for (epsilon in epsilon_grid) {
        local_population <- population[
          population$envelope_smoothing_epsilon == epsilon,
          , drop = FALSE
        ]
        for (endpoint in c("lower", "upper")) {
          endpoint_name <- paste0("mu1_", endpoint)
          evaluated <- oracle_envelope_tie_evaluate_mu1(
            pi0 = c(clipped_pi0),
            pi1 = c(clipped_pi1),
            endpoint = endpoint,
            sensitivity_lambda = specification$sensitivity_lambda,
            envelope_smoothing_epsilon = epsilon,
            treatment_effect = specification$treatment_effect
          )
          tie_estimate <- matrix(
            evaluated$estimate,
            nrow = local_replications,
            ncol = specification$tie_strata
          )
          gradient_pi0 <- matrix(
            evaluated$derivative_pi0,
            nrow = local_replications,
            ncol = specification$tie_strata
          )
          gradient_pi1 <- matrix(
            evaluated$derivative_pi1,
            nrow = local_replications,
            ncol = specification$tie_strata
          )
          regular <- oracle_envelope_tie_evaluate_mu1(
            pi0 = specification$pi0[regular_index],
            pi1 = specification$pi1[regular_index],
            endpoint = endpoint,
            sensitivity_lambda = specification$sensitivity_lambda,
            envelope_smoothing_epsilon = epsilon,
            treatment_effect = specification$treatment_effect
          )
          estimate <- tie_weight * rowSums(tie_estimate) +
            (1 - specification$tie_mass) * regular$estimate
          standard_error <- tie_weight * sqrt(rowSums(
            gradient_pi0^2 * variance_pi0 +
              gradient_pi1^2 * variance_pi1
          ))
          confidence_lower <- estimate - critical_value * standard_error
          confidence_upper <- estimate + critical_value * standard_error
          target <- local_population[[endpoint_name]][1L]
          sharp_target <- hard_population[[endpoint_name]][1L]
          translation <- sharp_target - target
          translated_lower <- confidence_lower + translation
          translated_upper <- confidence_upper + translation

          detail[[detail_counter]] <- data.frame(
            replication = completed + seq_len(local_replications),
            sample_size = n,
            epsilon = epsilon,
            endpoint = endpoint_name,
            population_target = target,
            sharp_population_target = sharp_target,
            target_shift_from_sharp = target - sharp_target,
            estimate = estimate,
            standard_error = standard_error,
            confidence_lower = confidence_lower,
            confidence_upper = confidence_upper,
            centered_coverage = confidence_lower <= target &
              confidence_upper >= target,
            translated_confidence_lower = translated_lower,
            translated_confidence_upper = translated_upper,
            translated_sharp_coverage = translated_lower <= sharp_target &
              translated_upper >= sharp_target,
            confidence_length = confidence_upper - confidence_lower,
            tie_branch_z1_rate = tie_branch_z1_rate,
            empirical_exact_tie_rate = empirical_exact_tie_rate,
            probability_stabilization_rate = stabilization_rate,
            minimum_tie_cell_count = apply(
              pmin(sufficient$n_z0, sufficient$n_z1), 1L, min
            ),
            stringsAsFactors = FALSE
          )
          detail_counter <- detail_counter + 1L
        }
      }
      completed <- completed + local_replications
    }
  }

  detail <- do.call(rbind, detail)
  rownames(detail) <- NULL
  groups <- split(
    detail,
    interaction(
      detail$sample_size, detail$epsilon, detail$endpoint,
      drop = TRUE
    )
  )
  summary <- do.call(rbind, lapply(groups, function(group) {
    z_score <- (group$estimate - group$population_target) /
      group$standard_error
    coverage <- mean(group$centered_coverage)
    data.frame(
      sample_size = group$sample_size[1L],
      epsilon = group$epsilon[1L],
      endpoint = group$endpoint[1L],
      replications = nrow(group),
      population_target = group$population_target[1L],
      sharp_population_target = group$sharp_population_target[1L],
      target_shift_from_sharp = group$target_shift_from_sharp[1L],
      mean_estimate = mean(group$estimate),
      bias = mean(group$estimate - group$population_target),
      empirical_standard_deviation = stats::sd(group$estimate),
      mean_standard_error = mean(group$standard_error),
      z_mean = mean(z_score),
      z_standard_deviation = stats::sd(z_score),
      centered_coverage = coverage,
      centered_coverage_mcse = sqrt(coverage * (1 - coverage) / nrow(group)),
      translated_sharp_coverage = mean(group$translated_sharp_coverage),
      mean_confidence_length = mean(group$confidence_length),
      mean_tie_branch_z1_rate = mean(group$tie_branch_z1_rate),
      mean_empirical_exact_tie_rate = mean(group$empirical_exact_tie_rate),
      mean_probability_stabilization_rate = mean(
        group$probability_stabilization_rate
      ),
      any_probability_stabilization_rate = mean(
        group$probability_stabilization_rate > 0
      ),
      maximum_probability_stabilization_rate = max(
        group$probability_stabilization_rate
      ),
      minimum_tie_cell_count = min(group$minimum_tie_cell_count),
      stringsAsFactors = FALSE
    )
  }))
  summary <- summary[order(
    summary$sample_size, summary$endpoint, summary$epsilon
  ), , drop = FALSE]
  rownames(summary) <- NULL

  list(
    specification = specification,
    configuration = list(
      sample_sizes = sample_sizes,
      replications = replications,
      epsilon_grid = epsilon_grid,
      seed_start = seed_start,
      batch_size = batch_size,
      confidence_level = confidence_level,
      numerical_probability_clip = numerical_probability_clip,
      isolated_nuisances = c("pi0_tie_strata", "pi1_tie_strata"),
      oracle_nuisances = c(
        "P(X)", "P(Z|X)", "outcome_law", "gamma_t",
        "regular_stratum_treatment_probabilities"
      )
    ),
    population_targets = population,
    detail = detail,
    summary = summary
  )
}
