# Five-dimensional extension of the bounded continuous-outcome DGP.
#
# Every measured covariate contributes to the instrument, treatment, and
# outcome through the standardized index
#
#   X_index = sum_j w_j X_j,  ||w||_2 = 1.
#
# The default equal weights preserve approximately the same variance of the
# measured-confounder signal as the original one-dimensional DGP.  Estimation
# nevertheless receives X1,...,Xd separately; X_index is retained only for
# oracle calculations and must not be included among the fitted predictors.

if (!exists("generate_iv_data", mode = "function")) {
  if (!file.exists("R/sharp_bounds_plugin.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/sharp_bounds_plugin.R")
}

multidimensional_x_columns <- function(dimension = 5L) {
  dimension <- as.integer(dimension)
  if (length(dimension) != 1L || is.na(dimension) || dimension < 1L) {
    stop("dimension must be one positive integer.")
  }
  paste0("X", seq_len(dimension))
}

resolve_multidimensional_x_weights <- function(
    dimension = 5L,
    x_weights = NULL
) {
  dimension <- as.integer(dimension)
  if (length(dimension) != 1L || is.na(dimension) || dimension < 1L) {
    stop("dimension must be one positive integer.")
  }
  if (is.null(x_weights)) {
    x_weights <- rep(1 / sqrt(dimension), dimension)
  }
  x_weights <- as.numeric(x_weights)
  if (length(x_weights) != dimension || any(!is.finite(x_weights))) {
    stop("x_weights must contain one finite value per covariate.")
  }
  norm <- sqrt(sum(x_weights^2))
  if (!is.finite(norm) || norm <= 0) {
    stop("x_weights must have positive Euclidean norm.")
  }
  x_weights / norm
}

multidimensional_x_index <- function(
    data,
    x_cols = multidimensional_x_columns(5L),
    x_weights = NULL
) {
  missing <- setdiff(x_cols, names(data))
  if (length(missing) > 0L) {
    stop("Missing covariates: ", paste(missing, collapse = ", "))
  }
  weights <- resolve_multidimensional_x_weights(
    length(x_cols), x_weights
  )
  design <- as.matrix(data[, x_cols, drop = FALSE])
  storage.mode(design) <- "double"
  if (any(!is.finite(design))) {
    stop("The multidimensional covariate matrix contains non-finite values.")
  }
  as.vector(design %*% weights)
}

get_multidimensional_outcome_support <- function(
    treatment_effect = 2,
    dimension = 5L,
    x_weights = NULL
) {
  if (length(treatment_effect) != 1L || !is.finite(treatment_effect)) {
    stop("treatment_effect must be one finite number.")
  }
  weights <- resolve_multidimensional_x_weights(dimension, x_weights)
  # Each X_j, U, and outcome error lies in [-2,2].
  x_index_bound <- 2 * sum(abs(weights))
  outcome_half_width <- 1.5 * x_index_bound + 1.5 * 2 + 2
  list(
    `0` = c(-outcome_half_width, outcome_half_width),
    `1` = treatment_effect + c(-outcome_half_width, outcome_half_width)
  )
}

generate_iv_data_multidimensional <- function(
    n = 5000,
    dimension = 5L,
    x_weights = NULL,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  x_cols <- multidimensional_x_columns(dimension)
  weights <- resolve_multidimensional_x_weights(dimension, x_weights)
  if (!is.null(seed)) set.seed(seed)

  x_matrix <- matrix(
    rnormtc(n * length(x_cols), mu = 0, s = 1),
    nrow = n,
    ncol = length(x_cols),
    dimnames = list(NULL, x_cols)
  )
  x_index <- as.vector(x_matrix %*% weights)
  U <- rnormtc(n, mu = 0, s = 1)
  prob_Z <- stats::plogis(0.5 * x_index)
  Z <- stats::rbinom(n, 1, prob_Z)
  prob_T <- stats::plogis(
    treatment_intercept + x_index + instrument_strength * Z +
      confounder_strength * U
  )
  T_obs <- stats::rbinom(n, 1, prob_T)

  mean_Y0 <- 1.5 * x_index + 1.5 * U
  mean_Y1 <- treatment_effect + 1.5 * x_index + 1.5 * U
  Y0 <- rnormtc(n, mu = mean_Y0, s = 1)
  Y1 <- rnormtc(n, mu = mean_Y1, s = 1)
  Y <- ifelse(T_obs == 1, Y1, Y0)

  cbind(
    data.frame(
      Y = Y,
      T_obs = T_obs,
      Z = Z,
      stringsAsFactors = FALSE
    ),
    as.data.frame(x_matrix),
    data.frame(
      X_index = x_index,
      U = U,
      Y0 = Y0,
      Y1 = Y1,
      prob_Z = prob_Z,
      prob_T = prob_T,
      stringsAsFactors = FALSE
    )
  )
}

compute_oracle_quantities_multidimensional <- function(
    data,
    x_cols = multidimensional_x_columns(5L),
    x_weights = NULL,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    n_nodes = 80L,
    batch_size = 1000L
) {
  x_index <- multidimensional_x_index(data, x_cols, x_weights)
  compute_oracle_quantities(
    data = data.frame(Y = data$Y, X = x_index),
    confounder_strength = confounder_strength,
    treatment_intercept = treatment_intercept,
    instrument_strength = instrument_strength,
    treatment_effect = treatment_effect,
    n_nodes = n_nodes,
    batch_size = batch_size
  )
}

add_oracle_quantities_multidimensional <- function(
    data,
    x_cols = multidimensional_x_columns(5L),
    x_weights = NULL,
    confounder_strength = 1,
    treatment_intercept = 0,
    instrument_strength = 1,
    treatment_effect = 2,
    n_nodes = 80L,
    batch_size = 1000L
) {
  oracle <- compute_oracle_quantities_multidimensional(
    data = data,
    x_cols = x_cols,
    x_weights = x_weights,
    confounder_strength = confounder_strength,
    treatment_intercept = treatment_intercept,
    instrument_strength = instrument_strength,
    treatment_effect = treatment_effect,
    n_nodes = n_nodes,
    batch_size = batch_size
  )
  cbind(data, oracle)
}

make_multidimensional_score_rule <- function(
    dimension = 5L,
    x_weights = NULL,
    n_nodes = 160L,
    reference_draws = 1000000L,
    seed = 47001L
) {
  dimension <- as.integer(dimension)
  n_nodes <- as.integer(n_nodes)
  reference_draws <- as.integer(reference_draws)
  if (n_nodes < 2L || reference_draws < 1000L) {
    stop("n_nodes must be at least 2 and reference_draws at least 1000.")
  }
  weights <- resolve_multidimensional_x_weights(dimension, x_weights)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  score <- numeric(reference_draws)
  for (index in seq_len(dimension)) {
    score <- score + weights[index] * rnormtc(
      reference_draws, mu = 0, s = 1
    )
  }
  probabilities <- (seq_len(n_nodes) - 0.5) / n_nodes
  nodes <- as.numeric(stats::quantile(
    score, probabilities, names = FALSE, type = 8
  ))
  list(
    nodes = nodes,
    weights = rep(1 / n_nodes, n_nodes),
    reference_draws = reference_draws,
    score_mean = mean(score),
    score_sd = stats::sd(score),
    score_minimum = min(score),
    score_maximum = max(score),
    dimension = dimension,
    x_weights = weights,
    seed = as.integer(seed)
  )
}
