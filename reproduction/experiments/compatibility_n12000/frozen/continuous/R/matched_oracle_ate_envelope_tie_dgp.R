# Matched oracle-envelope experiment for the final ATE envelope-tie DGP.
#
# This module uses exactly the observed-data law defined by
# make_fully_estimated_ate_envelope_tie_spec(), but removes estimation noise
# from every component except the treatment probabilities in the positive-mass
# tie strata.  Those probabilities are estimated by their empirical X,Z cell
# proportions and then propagated through either the hard or smoothed
# envelopes.  The X distribution, instrument propensity, outcome law,
# density ratios, regular-stratum treatment probabilities, thresholds, and
# outcome integrations are held at their population values.

matched_oracle_ate_require_dependencies <- function() {
  required <- c(
    "compute_projected_envelopes",
    "compute_eif_envelope_derivatives",
    "fully_estimated_ate_envelope_tie_residual_mean",
    "fully_estimated_ate_envelope_tie_residual_quantile",
    "fully_estimated_ate_envelope_tie_residual_partial_mean",
    "fully_estimated_ate_envelope_tie_population_grid"
  )
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing)) {
    stop(
      "Source R/sharp_bounds_plugin.R and ",
      "R/fully_estimated_ate_envelope_tie_dgp.R first. Missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

matched_oracle_ate_recycle <- function(value, length_out, name) {
  value <- as.numeric(value)
  if (!length(value) || any(!is.finite(value)) ||
      !(length(value) %in% c(1L, length_out))) {
    stop(name, " must be finite and have length one or the common length.")
  }
  rep_len(value, length_out)
}

matched_oracle_ate_evaluate_arm_endpoint <- function(
    pi0,
    pi1,
    treatment,
    endpoint = c("lower", "upper"),
    location = 0,
    specification,
    envelope_smoothing_epsilon = 0
) {
  matched_oracle_ate_require_dependencies()
  endpoint <- match.arg(endpoint)
  treatment <- as.integer(treatment)
  if (length(treatment) != 1L || is.na(treatment) ||
      !treatment %in% 0:1) {
    stop("treatment must equal zero or one.")
  }
  common_length <- max(length(pi0), length(pi1), length(location))
  pi0 <- matched_oracle_ate_recycle(pi0, common_length, "pi0")
  pi1 <- matched_oracle_ate_recycle(pi1, common_length, "pi1")
  location <- matched_oracle_ate_recycle(
    location, common_length, "location"
  )
  if (any(pi0 <= 0) || any(pi0 >= 1) ||
      any(pi1 <= 0) || any(pi1 >= 1)) {
    stop("pi0 and pi1 must lie strictly between zero and one.")
  }
  epsilon <- as.numeric(envelope_smoothing_epsilon)
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon < 0 ||
      epsilon >= 1 / specification$sensitivity_lambda) {
    stop("Invalid envelope_smoothing_epsilon.")
  }

  nuisance <- list(
    e = rep(1 - specification$probability_z1, common_length),
    pi0 = pi0,
    pi1 = pi1
  )
  gamma <- rep(1, common_length)
  envelope <- compute_projected_envelopes(
    nuisance = nuisance,
    gamma = gamma,
    t_value = treatment,
    Lambda0 = specification$sensitivity_lambda,
    Lambda1 = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = epsilon,
    gamma_compatibility_action = "stop"
  )
  derivative <- compute_eif_envelope_derivatives(
    nuisance = nuisance,
    gamma = gamma,
    t_value = treatment,
    Lambda0 = specification$sensitivity_lambda,
    Lambda1 = specification$sensitivity_lambda,
    envelope_smoothing_epsilon = epsilon
  )
  L <- envelope$L
  U <- envelope$U
  if (any(!is.finite(L)) || any(!is.finite(U)) || any(L >= U)) {
    stop("Estimated tie probabilities produced invalid envelopes.")
  }

  # compute_eif_envelope_derivatives() differentiates with respect to the
  # treatment-arm probabilities r_z=P(T=t|X,Z=z).  Convert those derivatives
  # to derivatives with respect to pi_z=P(T=1|X,Z=z).
  pi_sign <- if (treatment == 1L) 1 else -1
  derivative_L_pi0 <- pi_sign * derivative$lower$r0
  derivative_L_pi1 <- pi_sign * derivative$lower$r1
  derivative_U_pi0 <- pi_sign * derivative$upper$r0
  derivative_U_pi1 <- pi_sign * derivative$upper$r1

  residual_mean <- fully_estimated_ate_envelope_tie_residual_mean(
    specification
  )
  outcome_mean <- location + specification$treatment_effect * treatment +
    residual_mean
  outcome_location <- location + specification$treatment_effect * treatment

  if (endpoint == "lower") {
    cutoff_probability <- (1 - L) / (U - L)
    lower_partial_moment <- outcome_location * cutoff_probability +
      fully_estimated_ate_envelope_tie_residual_partial_mean(
        cutoff_probability, specification
      )
    tilted_mean <- L * outcome_mean +
      (U - L) * lower_partial_moment

    derivative_tilted_mean <- function(derivative_L, derivative_U) {
      derivative_probability <- (
        -derivative_L * (U - L) -
          (1 - L) * (derivative_U - derivative_L)
      ) / (U - L)^2
      cutoff_residual <- fully_estimated_ate_envelope_tie_residual_quantile(
        cutoff_probability, specification
      )
      derivative_L * outcome_mean +
        (derivative_U - derivative_L) * lower_partial_moment +
        (U - L) * (outcome_location + cutoff_residual) *
          derivative_probability
    }
  } else {
    cutoff_probability <- (U - 1) / (U - L)
    lower_partial_moment <- outcome_location * cutoff_probability +
      fully_estimated_ate_envelope_tie_residual_partial_mean(
        cutoff_probability, specification
      )
    tilted_mean <- U * outcome_mean +
      (L - U) * lower_partial_moment

    derivative_tilted_mean <- function(derivative_L, derivative_U) {
      derivative_probability <- (
        derivative_U * (U - L) -
          (U - 1) * (derivative_U - derivative_L)
      ) / (U - L)^2
      cutoff_residual <- fully_estimated_ate_envelope_tie_residual_quantile(
        cutoff_probability, specification
      )
      derivative_U * outcome_mean +
        (derivative_L - derivative_U) * lower_partial_moment +
        (L - U) * (outcome_location + cutoff_residual) *
          derivative_probability
    }
  }

  derivative_tilted_pi0 <- derivative_tilted_mean(
    derivative_L_pi0, derivative_U_pi0
  )
  derivative_tilted_pi1 <- derivative_tilted_mean(
    derivative_L_pi1, derivative_U_pi1
  )
  r0 <- if (treatment == 1L) pi0 else 1 - pi0
  estimate <- r0 * outcome_mean + (1 - r0) * tilted_mean
  derivative_pi0 <- pi_sign * (outcome_mean - tilted_mean) +
    (1 - r0) * derivative_tilted_pi0
  derivative_pi1 <- (1 - r0) * derivative_tilted_pi1

  list(
    estimate = estimate,
    derivative_pi0 = derivative_pi0,
    derivative_pi1 = derivative_pi1,
    L = L,
    U = U,
    candidate_lower = envelope$candidate_lower,
    candidate_upper = envelope$candidate_upper,
    lower_derivative_weight = derivative$lower_derivative_weight,
    upper_derivative_weight = derivative$upper_derivative_weight,
    cutoff_probability = cutoff_probability,
    tilted_mean = tilted_mean,
    treatment = treatment,
    endpoint = endpoint,
    envelope_smoothing_epsilon = epsilon
  )
}

matched_oracle_ate_draw_sufficient_statistics <- function(
    n,
    replications,
    specification,
    seed
) {
  n <- as.integer(n)
  replications <- as.integer(replications)
  if (length(n) != 1L || is.na(n) || n < 200L ||
      length(replications) != 1L || is.na(replications) ||
      replications < 1L) {
    stop("Invalid sufficient-statistic simulation dimensions.")
  }
  set.seed(as.integer(seed))
  joint_probability <- as.vector(outer(
    specification$x_probabilities,
    c(1 - specification$probability_z1, specification$probability_z1)
  ))
  cell_count <- t(stats::rmultinom(
    replications, size = n, prob = joint_probability
  ))
  tie <- seq_len(specification$tie_strata)
  z1_offset <- length(specification$x_values)
  n_z0 <- cell_count[, tie, drop = FALSE]
  n_z1 <- cell_count[, z1_offset + tie, drop = FALSE]
  if (any(n_z0 == 0L) || any(n_z1 == 0L)) {
    stop("A simulated tie-stratum instrument cell was empty; increase n.")
  }
  probability_z0 <- rep(specification$pi0[tie], each = replications)
  probability_z1 <- rep(specification$pi1[tie], each = replications)
  treated_z0 <- matrix(
    stats::rbinom(length(n_z0), size = c(n_z0), prob = probability_z0),
    nrow = replications,
    ncol = specification$tie_strata
  )
  treated_z1 <- matrix(
    stats::rbinom(length(n_z1), size = c(n_z1), prob = probability_z1),
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

matched_oracle_ate_summarize_endpoints <- function(detail) {
  if (!nrow(detail)) return(data.frame())
  groups <- split(
    detail,
    interaction(
      detail$sample_size, detail$epsilon, detail$endpoint,
      drop = TRUE, lex.order = TRUE
    )
  )
  output <- do.call(rbind, lapply(groups, function(group) {
    z_score <- (group$estimate - group$population_target) /
      group$standard_error
    data.frame(
      sample_size = group$sample_size[1L],
      epsilon = group$epsilon[1L],
      endpoint = group$endpoint[1L],
      replications = nrow(group),
      population_target = group$population_target[1L],
      sharp_population_target = group$sharp_population_target[1L],
      target_shift_from_sharp = group$target_shift_from_sharp[1L],
      bias = mean(group$estimate - group$population_target),
      empirical_standard_deviation = stats::sd(group$estimate),
      mean_standard_error = mean(group$standard_error),
      z_mean = mean(z_score),
      z_standard_deviation = stats::sd(z_score),
      own_target_one_sided_97_5_coverage = mean(
        group$own_target_one_sided_covered
      ),
      sharp_target_one_sided_97_5_coverage = mean(
        group$sharp_target_one_sided_covered
      ),
      translated_sharp_one_sided_97_5_coverage = mean(
        group$translated_sharp_one_sided_covered
      ),
      own_target_two_sided_95_coverage = mean(
        group$own_target_two_sided_covered
      ),
      stringsAsFactors = FALSE
    )
  }))
  output[order(output$sample_size, output$epsilon, output$endpoint), ,
         drop = FALSE]
}

matched_oracle_ate_summarize_intervals <- function(detail) {
  if (!nrow(detail)) return(data.frame())
  groups <- split(
    detail,
    interaction(
      detail$sample_size, detail$epsilon,
      drop = TRUE, lex.order = TRUE
    )
  )
  output <- do.call(rbind, lapply(groups, function(group) {
    data.frame(
      sample_size = group$sample_size[1L],
      epsilon = group$epsilon[1L],
      replications = nrow(group),
      population_lower = group$population_lower[1L],
      population_upper = group$population_upper[1L],
      sharp_population_lower = group$sharp_population_lower[1L],
      sharp_population_upper = group$sharp_population_upper[1L],
      mean_point_width = mean(group$point_width),
      mean_outward_width = mean(group$outward_width),
      mean_translated_outward_width = mean(
        group$translated_outward_width
      ),
      own_target_outward_95_coverage = mean(
        group$own_target_outward_covered
      ),
      sharp_target_outward_95_coverage = mean(
        group$sharp_target_outward_covered
      ),
      translated_sharp_outward_95_coverage = mean(
        group$translated_sharp_outward_covered
      ),
      stringsAsFactors = FALSE
    )
  }))
  output[order(output$sample_size, output$epsilon), , drop = FALSE]
}

run_matched_oracle_ate_envelope_tie_simulation <- function(
    sample_sizes = c(1500L, 3000L, 6000L),
    replications = 20000L,
    epsilon_grid = c(0, 0.01),
    specification,
    seed_start = 264000L,
    batch_size = 5000L,
    confidence_level = 0.95,
    numerical_probability_clip = 1e-8
) {
  matched_oracle_ate_require_dependencies()
  sample_sizes <- as.integer(sample_sizes)
  replications <- as.integer(replications)
  batch_size <- as.integer(batch_size)
  epsilon_grid <- sort(unique(as.numeric(epsilon_grid)))
  if (!length(sample_sizes) || anyNA(sample_sizes) ||
      any(sample_sizes < 200L) || replications < 1L || batch_size < 1L ||
      any(!is.finite(epsilon_grid)) || !any(epsilon_grid == 0) ||
      any(epsilon_grid < 0) ||
      any(epsilon_grid >= 1 / specification$sensitivity_lambda) ||
      !is.finite(confidence_level) || confidence_level <= 0 ||
      confidence_level >= 1 || !is.finite(numerical_probability_clip) ||
      numerical_probability_clip <= 0 ||
      numerical_probability_clip >= specification$overlap_constant) {
    stop("Invalid matched oracle ATE simulation settings.")
  }

  population <- fully_estimated_ate_envelope_tie_population_grid(
    specification = specification,
    epsilon_grid = epsilon_grid
  )
  hard_population <- population[
    population$envelope_smoothing_epsilon == 0, , drop = FALSE
  ]
  critical_value <- stats::qnorm(1 - (1 - confidence_level) / 2)
  tie <- seq_len(specification$tie_strata)
  tie_weight <- specification$tie_mass / specification$tie_strata
  regular <- specification$regular_stratum
  endpoint_detail <- interval_detail <- diagnostics <- list()
  endpoint_counter <- interval_counter <- diagnostic_counter <- 1L

  for (sample_index in seq_along(sample_sizes)) {
    sample_size <- sample_sizes[sample_index]
    completed <- 0L
    batch_index <- 0L
    while (completed < replications) {
      local_replications <- min(batch_size, replications - completed)
      batch_index <- batch_index + 1L
      sufficient <- matched_oracle_ate_draw_sufficient_statistics(
        n = sample_size,
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
      variance_pi0 <- clipped_pi0 * (1 - clipped_pi0) / sufficient$n_z0
      variance_pi1 <- clipped_pi1 * (1 - clipped_pi1) / sufficient$n_z1
      stabilization_rate <- rowMeans(
        clipped_pi0 != raw_pi0 | clipped_pi1 != raw_pi1
      )
      branch_z1_rate <- rowMeans(clipped_pi1 > clipped_pi0)
      empirical_tie_rate <- rowMeans(clipped_pi1 == clipped_pi0)
      minimum_cell_count <- apply(
        pmin(sufficient$n_z0, sufficient$n_z1), 1L, min
      )
      replication_index <- completed + seq_len(local_replications)

      flattened_pi0 <- c(clipped_pi0)
      flattened_pi1 <- c(clipped_pi1)
      flattened_location <- rep(
        specification$outcome_x_shifts[tie],
        each = local_replications
      )

      aggregate_arm_endpoint <- function(treatment, endpoint, epsilon) {
        evaluated <- matched_oracle_ate_evaluate_arm_endpoint(
          pi0 = flattened_pi0,
          pi1 = flattened_pi1,
          treatment = treatment,
          endpoint = endpoint,
          location = flattened_location,
          specification = specification,
          envelope_smoothing_epsilon = epsilon
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
        regular_evaluated <- matched_oracle_ate_evaluate_arm_endpoint(
          pi0 = specification$pi0[regular],
          pi1 = specification$pi1[regular],
          treatment = treatment,
          endpoint = endpoint,
          location = specification$outcome_x_shifts[regular],
          specification = specification,
          envelope_smoothing_epsilon = epsilon
        )
        list(
          estimate = tie_weight * rowSums(tie_estimate) +
            (1 - specification$tie_mass) * regular_evaluated$estimate,
          gradient_pi0 = tie_weight * gradient_pi0,
          gradient_pi1 = tie_weight * gradient_pi1
        )
      }

      for (epsilon in epsilon_grid) {
        arm <- list(
          mu0_lower = aggregate_arm_endpoint(0L, "lower", epsilon),
          mu0_upper = aggregate_arm_endpoint(0L, "upper", epsilon),
          mu1_lower = aggregate_arm_endpoint(1L, "lower", epsilon),
          mu1_upper = aggregate_arm_endpoint(1L, "upper", epsilon)
        )
        ate <- list(
          lower = list(
            estimate = arm$mu1_lower$estimate - arm$mu0_upper$estimate,
            gradient_pi0 = arm$mu1_lower$gradient_pi0 -
              arm$mu0_upper$gradient_pi0,
            gradient_pi1 = arm$mu1_lower$gradient_pi1 -
              arm$mu0_upper$gradient_pi1
          ),
          upper = list(
            estimate = arm$mu1_upper$estimate - arm$mu0_lower$estimate,
            gradient_pi0 = arm$mu1_upper$gradient_pi0 -
              arm$mu0_lower$gradient_pi0,
            gradient_pi1 = arm$mu1_upper$gradient_pi1 -
              arm$mu0_lower$gradient_pi1
          )
        )
        target <- population[
          population$envelope_smoothing_epsilon == epsilon,
          , drop = FALSE
        ]
        endpoint_rows <- list()
        for (endpoint in c("lower", "upper")) {
          endpoint_name <- paste0("ate_", endpoint)
          local <- ate[[endpoint]]
          standard_error <- sqrt(rowSums(
            local$gradient_pi0^2 * variance_pi0 +
              local$gradient_pi1^2 * variance_pi1
          ))
          confidence_lower <- local$estimate -
            critical_value * standard_error
          confidence_upper <- local$estimate +
            critical_value * standard_error
          population_target <- target[[endpoint_name]][1L]
          sharp_target <- hard_population[[endpoint_name]][1L]
          translation <- sharp_target - population_target
          translated_lower <- confidence_lower + translation
          translated_upper <- confidence_upper + translation
          lower_endpoint <- identical(endpoint, "lower")
          endpoint_rows[[endpoint]] <- data.frame(
            replication = replication_index,
            sample_size = sample_size,
            epsilon = epsilon,
            endpoint = endpoint_name,
            population_target = population_target,
            sharp_population_target = sharp_target,
            target_shift_from_sharp = population_target - sharp_target,
            estimate = local$estimate,
            standard_error = standard_error,
            confidence_lower = confidence_lower,
            confidence_upper = confidence_upper,
            translated_confidence_lower = translated_lower,
            translated_confidence_upper = translated_upper,
            own_target_one_sided_covered = if (lower_endpoint) {
              confidence_lower <= population_target
            } else {
              confidence_upper >= population_target
            },
            sharp_target_one_sided_covered = if (lower_endpoint) {
              confidence_lower <= sharp_target
            } else {
              confidence_upper >= sharp_target
            },
            translated_sharp_one_sided_covered = if (lower_endpoint) {
              translated_lower <= sharp_target
            } else {
              translated_upper >= sharp_target
            },
            own_target_two_sided_covered =
              confidence_lower <= population_target &
              confidence_upper >= population_target,
            stringsAsFactors = FALSE
          )
          endpoint_detail[[endpoint_counter]] <- endpoint_rows[[endpoint]]
          endpoint_counter <- endpoint_counter + 1L
        }

        lower_row <- endpoint_rows$lower
        upper_row <- endpoint_rows$upper
        interval_detail[[interval_counter]] <- data.frame(
          replication = replication_index,
          sample_size = sample_size,
          epsilon = epsilon,
          population_lower = lower_row$population_target,
          population_upper = upper_row$population_target,
          sharp_population_lower = lower_row$sharp_population_target,
          sharp_population_upper = upper_row$sharp_population_target,
          estimate_lower = lower_row$estimate,
          estimate_upper = upper_row$estimate,
          point_width = upper_row$estimate - lower_row$estimate,
          confidence_lower = lower_row$confidence_lower,
          confidence_upper = upper_row$confidence_upper,
          outward_width = upper_row$confidence_upper -
            lower_row$confidence_lower,
          translated_confidence_lower =
            lower_row$translated_confidence_lower,
          translated_confidence_upper =
            upper_row$translated_confidence_upper,
          translated_outward_width =
            upper_row$translated_confidence_upper -
            lower_row$translated_confidence_lower,
          own_target_outward_covered =
            lower_row$own_target_one_sided_covered &
            upper_row$own_target_one_sided_covered,
          sharp_target_outward_covered =
            lower_row$sharp_target_one_sided_covered &
            upper_row$sharp_target_one_sided_covered,
          translated_sharp_outward_covered =
            lower_row$translated_sharp_one_sided_covered &
            upper_row$translated_sharp_one_sided_covered,
          stringsAsFactors = FALSE
        )
        interval_counter <- interval_counter + 1L

        diagnostics[[diagnostic_counter]] <- data.frame(
          replication = replication_index,
          sample_size = sample_size,
          epsilon = epsilon,
          tie_branch_z1_rate = branch_z1_rate,
          empirical_exact_tie_rate = empirical_tie_rate,
          probability_stabilization_rate = stabilization_rate,
          minimum_tie_cell_count = minimum_cell_count,
          stringsAsFactors = FALSE
        )
        diagnostic_counter <- diagnostic_counter + 1L
      }
      completed <- completed + local_replications
    }
  }

  endpoint_detail <- do.call(rbind, endpoint_detail)
  interval_detail <- do.call(rbind, interval_detail)
  diagnostics <- do.call(rbind, diagnostics)
  rownames(endpoint_detail) <- NULL
  rownames(interval_detail) <- NULL
  rownames(diagnostics) <- NULL
  list(
    specification = specification,
    population_targets = population,
    endpoint_detail = endpoint_detail,
    endpoint_summary = matched_oracle_ate_summarize_endpoints(
      endpoint_detail
    ),
    interval_detail = interval_detail,
    interval_summary = matched_oracle_ate_summarize_intervals(
      interval_detail
    ),
    diagnostics = diagnostics,
    configuration = list(
      sample_sizes = sample_sizes,
      replications = replications,
      epsilon_grid = epsilon_grid,
      seed_start = as.integer(seed_start),
      batch_size = batch_size,
      confidence_level = confidence_level,
      critical_value = critical_value,
      numerical_probability_clip = numerical_probability_clip
    )
  )
}
