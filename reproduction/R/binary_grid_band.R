# Simultaneous outward calibration from cross-fitted binary ATE signals.
# Rows must be independent analysis units, in identical order across columns.
# This is NOT a nuisance fitter, compatibility test, or clustered/weighted SE.
# Validity requires the paper's uniform expansion/stability conditions.
binary_grid_band <- function(grid, lower_signal, upper_signal,
                             level = .95, draws = 10000L, seed = 1L) {
  lo <- as.matrix(lower_signal); hi <- as.matrix(upper_signal)
  stopifnot(is.numeric(grid), length(grid) >= 2L, all(is.finite(grid)),
    min(grid) >= 1, all(diff(grid) > 0), identical(dim(lo), dim(hi)),
    nrow(lo) >= 3L, ncol(lo) == length(grid),
    all(is.finite(lo)), all(is.finite(hi)),
    length(level) == 1L, is.finite(level), level > 0, level < 1,
    length(draws) == 1L, is.finite(draws), draws >= 100,
    draws == as.integer(draws), length(seed) == 1L, is.finite(seed))
  # Interleave lower/upper signals to match the historical OHIE calibration.
  S <- matrix(NA_real_, nrow(lo), 2L * ncol(lo))
  S[, seq(1L, ncol(S), 2L)] <- lo
  S[, seq(2L, ncol(S), 2L)] <- hi
  standard_error <- apply(S, 2L, sd) / sqrt(nrow(S))
  if (any(standard_error <= 0)) stop('Zero-variance signals require separate treatment.')
  signs <- rep(c(1, -1), ncol(lo))
  covariance <- cor(S) * outer(signs, signs) * (nrow(S) - 1) / nrow(S)
  eig <- eigen(covariance, symmetric = TRUE)
  if (min(eig$values) < -1e-8) stop('Invalid multiplier covariance.')
  A <- sweep(eig$vectors, 2L, sqrt(pmax(eig$values, 0)), '*')
  had_seed <- exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get('.Random.seed', envir = .GlobalEnv)
  on.exit(if (had_seed) assign('.Random.seed', old_seed, envir = .GlobalEnv)
          else if (exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE))
            rm('.Random.seed', envir = .GlobalEnv))
  set.seed(seed)
  maxima <- numeric(draws)
  for (start in seq.int(1L, draws, by = 2000L)) {
    ix <- start:min(start + 1999L, draws)
    z <- A %*% matrix(rnorm(ncol(S) * length(ix)), nrow = ncol(S))
    maxima[ix] <- apply(z, 2L, max)
  }
  critical <- unname(quantile(maxima, level, type = 8))
  out <- data.frame(Lambda = grid, lower = colMeans(lo), upper = colMeans(hi),
    se_lower = standard_error[seq(1L, ncol(S), 2L)],
    se_upper = standard_error[seq(2L, ncol(S), 2L)])
  out$simultaneous_lower <- out$lower - critical * out$se_lower
  out$simultaneous_upper <- out$upper + critical * out$se_upper
  out$empty <- out$simultaneous_lower > out$simultaneous_upper
  list(table = out, critical = critical, covariance = covariance,
       level = level, draws = draws, seed = seed,
       scope = 'Grid-only outward band; iid units; no projection or interval repair.')
}
