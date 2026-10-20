#!/usr/bin/env Rscript
# Deterministic population MSM targets for the completed binary X5 stress
# comparators. No data simulation, nuisance fitting, or source mutation.
#
# Run from the project root:
# Rscript scripts/analysis_binary_x5_comparator_population_v1.R [project_root] [comparator_summary.rds]
# Existing output is never overwritten. Source this file to access the
# calculation functions without running integration or writing an artifact.

binary_x5_msm_conditional <- function(r, p, Lambda) {
  stopifnot(Lambda >= 1, all(r > 0 & r < 1), all(p > 0 & p < 1))
  list(
    lower = r * p + (1 - r) * pmax(p / Lambda, 1 - Lambda * (1 - p)),
    upper = r * p + (1 - r) * pmin(Lambda * p, 1 - (1 - p) / Lambda)
  )
}

binary_x5_msm_quadrature <- function(order) {
  stopifnot(length(order) == 1L, order >= 2L, order == as.integer(order))
  i <- seq_len(order - 1L)
  jacobi <- matrix(0, order, order)
  jacobi[cbind(i, i + 1L)] <- i / sqrt(4 * i^2 - 1)
  eig <- eigen(jacobi + t(jacobi), symmetric = TRUE)
  u <- (eig$values + 1) / 2
  weights1 <- eig$vectors[1L, ]^2
  index <- as.matrix(expand.grid(rep(list(seq_len(order)), 5L)))
  x <- matrix(stats::qnorm(stats::pnorm(-2) +
                            u[index] * (stats::pnorm(2) - stats::pnorm(-2))),
              ncol = 5L)
  weights <- Reduce(`*`, lapply(seq_len(5L), function(j) weights1[index[, j]]))
  stopifnot(all(weights > 0), abs(sum(weights) - 1) < 1e-12)
  list(
    h = drop(x %*% c(.14, .12, .10, .08, .06)),
    m = drop(x %*% c(.20, -.25, .30, .15, -.20)),
    shape = 1 + .6 * sin(x[, 1L] * x[, 2L]) + .4 * sin(x[, 3L] * x[, 4L]),
    weights = weights, order = order, points = nrow(x),
    weight_sum_error = abs(sum(weights) - 1)
  )
}

binary_x5_msm_integrate <- function(order, configuration) {
  rule <- binary_x5_msm_quadrature(order)
  integrate <- function(value) sum(rule$weights * value)
  Lambda_grid <- sort(unique(c(1, configuration$Lambda)))
  rows <- checks <- list()
  for (setting in seq_len(nrow(configuration$settings))) {
    design <- configuration$settings[setting, ]
    for (Lambda in Lambda_grid) {
      arms <- list()
      for (t in 0:1) {
        r0 <- stats::plogis(rule$h - design$a / 2)
        r1 <- stats::plogis(rule$h + design$a / 2)
        if (t == 0L) {
          r0 <- 1 - r0
          r1 <- 1 - r1
        }
        p0 <- stats::plogis(-.4 + 1.5 * t + rule$m)
        p1 <- stats::plogis(-.4 + 1.5 * t + rule$m +
                             (2 * t - 1) * design$b * rule$shape)
        # X-only observed outcome probabilities mix over Z CONDITIONAL on
        # treatment. Z is marginally equiprobable, but P(Z|T=t,X) is not .5.
        r <- (r0 + r1) / 2
        p <- (r0 * p0 + r1 * p1) / (r0 + r1)
        x_bounds <- binary_x5_msm_conditional(r, p, Lambda)
        z0_bounds <- binary_x5_msm_conditional(r0, p0, Lambda)
        z1_bounds <- binary_x5_msm_conditional(r1, p1, Lambda)
        xz_bounds <- list(lower = (z0_bounds$lower + z1_bounds$lower) / 2,
                          upper = (z0_bounds$upper + z1_bounds$upper) / 2)
        arms[[as.character(t)]] <- list(
          x = vapply(x_bounds, integrate, numeric(1)),
          xz = vapply(xz_bounds, integrate, numeric(1))
        )
        checks[[length(checks) + 1L]] <- data.frame(
          order = order, setting = setting, treatment = t, Lambda = Lambda,
          treatment_mixture_error = max(abs(r * p - (r0 * p0 + r1 * p1) / 2)),
          weight_normalization_error = rule$weight_sum_error,
          lambda_one_collapse_error = if (Lambda == 1) {
            max(abs(x_bounds$lower - p), abs(x_bounds$upper - p),
                abs(xz_bounds$lower - (p0 + p1) / 2),
                abs(xz_bounds$upper - (p0 + p1) / 2))
          } else 0,
          conditional_crossing = max(0, x_bounds$lower - x_bounds$upper,
                                    xz_bounds$lower - xz_bounds$upper)
        )
      }
      for (adjustment in c("x", "xz")) {
        a0 <- arms[["0"]][[adjustment]]
        a1 <- arms[["1"]][[adjustment]]
        rows[[length(rows) + 1L]] <- data.frame(
          setting = setting, a = design$a, b = design$b,
          Lambda = Lambda, adjustment = adjustment, order = order,
          mu0_lower = unname(a0["lower"]), mu0_upper = unname(a0["upper"]),
          mu1_lower = unname(a1["lower"]), mu1_upper = unname(a1["upper"]),
          lower = unname(a1["lower"] - a0["upper"]),
          upper = unname(a1["upper"] - a0["lower"])
        )
      }
    }
  }
  list(targets = do.call(rbind, rows), checks = do.call(rbind, checks),
       quadrature = list(order = order, points = rule$points,
                         rule = "Tensor Gauss-Legendre on five truncated-normal probability coordinates"))
}

binary_x5_msm_population <- function(project_root = getwd(), comparator_summary = NULL) {
  project_root <- normalizePath(project_root, mustWork = TRUE)
  source_paths <- c(
    iv_design = "experiments/binary_bounds/scripts/common.R",
    comparator_design = "experiments/binary_comparators/R/design.R",
    iv_population = "experiments/binary_bounds/targets/population.rds"
  )
  absolute_paths <- file.path(project_root, source_paths)
  if (!all(file.exists(absolute_paths)))
    stop("Generate binary population targets with paper.R experiment binary_bounds targets first.")
  observed_hash <- unname(tools::md5sum(absolute_paths))
  # Hashes record locally generated target provenance; serialized outputs are
  # intentionally not distributed. Generator parity is checked below.
  iv_population <- readRDS(absolute_paths[3L])
  configuration <- iv_population$config
  stopifnot(configuration$n == 6000L, configuration$reps == 500L,
            identical(configuration$Lambda, c(1.8, 2, 2.5, 3)),
            identical(configuration$settings$a, c(.55, 1.1, 1.65, 1.1)),
            identical(configuration$settings$b, c(.075, .075, .075, .15)))
  extract_generator <- function(path) {
    expressions <- as.list(parse(path))
    hits <- vapply(expressions, function(expr) {
      is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        identical(expr[[2L]], as.name("generate"))
    }, logical(1))
    stopifnot(sum(hits) == 1L)
    expressions[[which(hits)]][[3L]]
  }
  stopifnot(identical(extract_generator(absolute_paths[1L]),
                      extract_generator(absolute_paths[2L])))
  uploaded <- NULL
  if (!is.null(comparator_summary)) {
    comparator_summary <- normalizePath(comparator_summary, mustWork = TRUE)
    uploaded <- readRDS(comparator_summary)
    stopifnot(identical(uploaded$signature$config, configuration),
              identical(unname(uploaded$signature$manifest["R/design.R"]), observed_hash[2L]))
  }
  started <- proc.time()[["elapsed"]]
  coarse <- binary_x5_msm_integrate(11L, configuration)
  fine <- binary_x5_msm_integrate(15L, configuration)
  columns <- c("mu0_lower", "mu0_upper", "mu1_lower", "mu1_upper", "lower", "upper")
  keys <- c("setting", "Lambda", "adjustment")
  stopifnot(identical(coarse$targets[keys], fine$targets[keys]))
  refinement <- fine$targets[keys]
  errors <- abs(as.matrix(fine$targets[columns]) - as.matrix(coarse$targets[columns]))
  refinement$max_endpoint_error <- apply(errors, 1L, max)
  refinement$max_ate_error <- apply(errors[, c("lower", "upper"), drop = FALSE], 1L, max)
  targets <- subset(fine$targets, Lambda != 1)
  nesting <- merge(targets[c(keys, "lower", "upper")], iv_population$targets,
                   by = c("setting", "Lambda"), suffixes = c("_msm", "_iv"))
  nesting$lower_gap <- nesting$lower_iv - nesting$lower_msm
  nesting$upper_gap <- nesting$upper_msm - nesting$upper_iv
  checks <- rbind(coarse$checks, fine$checks)
  max_refinement <- max(refinement$max_endpoint_error)
  max_ate_refinement <- max(refinement$max_ate_error)
  tolerance <- 1e-5
  maximum_algebraic_error <- max(as.matrix(checks[c(
    "treatment_mixture_error", "weight_normalization_error",
    "lambda_one_collapse_error", "conditional_crossing")]))
  nesting_violation <- max(0, -nesting$lower_gap, -nesting$upper_gap)
  validated <- max_refinement < tolerance && maximum_algebraic_error < 1e-12 &&
    nesting_violation < tolerance && all(is.finite(as.matrix(targets[columns])))
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    version = "binary_x5_comparator_population_v1",
    validated = validated,
    validation_status = if (validated) "validated_numerically" else "not_validated",
    targets = targets,
    no_unmeasured_confounding_benchmarks = subset(fine$targets, Lambda == 1),
    coarse = coarse$targets,
    refinement = refinement,
    max_refinement = max_refinement,
    max_ate_refinement = max_ate_refinement,
    validation_tolerance = tolerance,
    checks = checks,
    maximum_algebraic_error = maximum_algebraic_error,
    iv_nesting = nesting,
    maximum_iv_nesting_violation = nesting_violation,
    iv_population_refinement = iv_population$max_refinement,
    configuration = configuration,
    formulas = list(
      x_treatment = "r=(r0+r1)/2",
      x_outcome = "p=(r0*p0+r1*p1)/(r0+r1)",
      xz_averaging = "Equal P(Z=z|X)=1/2 weights AFTER conditional MSM extrema",
      lower = "r*p+(1-r)*max(p/Lambda,1-Lambda*(1-p))",
      upper = "r*p+(1-r)*min(Lambda*p,1-(1-p)/Lambda)",
      ate = "lower=mu1_lower-mu0_upper; upper=mu1_upper-mu0_lower"
    ),
    provenance = list(
      sources = data.frame(role = names(source_paths), path = absolute_paths,
                            md5 = observed_hash, row.names = NULL),
      parsed_generators_identical = TRUE,
      comparator_summary = comparator_summary,
      comparator_summary_md5 = if (!is.null(comparator_summary))
        unname(tools::md5sum(comparator_summary)) else NULL,
      comparator_implementation = if (!is.null(uploaded)) uploaded$signature$implementation else NULL,
      analysis_script_md5 = unname(tools::md5sum(file.path(project_root,
        "scripts/analysis_binary_x5_comparator_population_v1.R"))),
      R = as.character(getRversion()),
      computed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    integration = list(coarse = coarse$quadrature, fine = fine$quadrature,
                       runtime_seconds = elapsed, processes = 1L),
    caveat = paste("Numerically validated population MSM bounds, separately for X and XZ.",
      "Quadrature refinement is not a certified integration error bound.",
      "Lambda=1 values are no-unmeasured-confounding benchmarks, not asserted causal means.",
      "No Monte Carlo generation, nuisance fitting, or classifier gamma estimation.")
  )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  root <- if (length(args)) args[1L] else getwd()
  comparator_summary <- if (length(args) > 1L) args[2L] else NULL
  output <- file.path(root, "results/binary_x5_comparator_population_v1/targets.rds")
  if (file.exists(output)) stop("Output already exists; refusing to overwrite: ", output)
  result <- binary_x5_msm_population(root, comparator_summary)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("unpublished_", tmpdir = dirname(output))
  saveRDS(result, temporary)
  stopifnot(identical(readRDS(temporary), result))
  if (!file.rename(temporary, output)) stop("Could not publish population target artifact.")
  cat("Status:", result$validation_status, "\nOutput:", output,
      "\nRuntime seconds:", result$integration$runtime_seconds,
      "\nMaximum endpoint refinement:", result$max_refinement,
      "\nMaximum ATE refinement:", result$max_ate_refinement,
      "\nMaximum algebraic error:", result$maximum_algebraic_error,
      "\nMaximum IV nesting violation:", result$maximum_iv_nesting_violation, "\n")
  print(result$targets[c("setting", "Lambda", "adjustment", "lower", "upper")], row.names = FALSE)
  if (!result$validated) quit(status = 1L)
}
