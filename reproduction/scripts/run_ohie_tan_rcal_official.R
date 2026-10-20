#!/usr/bin/env Rscript

# Tan's official regularized calibrated sensitivity estimator for the Oregon
# Health Insurance Experiment. This runner calls RCALsa::sa.ate.regu.cv
# directly, uses the main-effect/two-way-interaction library and the 0.8%
# sparsity filter from Tan's empirical application, and saves each lower or
# upper endpoint as soon as it completes. The reported point estimate is the
# package's augmented IPW estimator with relaxation (the `est` row).
# Official package archive:
# https://sites.sas.rutgers.edu/ztan/Publication/RCALsa_1.0.zip
#
# Examples:
#   RCALSA_LIBRARY=/path/to/rcalsa-library \
#     Rscript scripts/run_ohie_tan_rcal_official.R
#
#   IVMSM_OHIE_TAN_ADJUSTMENT=xz IVMSM_OHIE_TAN_LAMBDAS=1.6 \
#     RCALSA_LIBRARY=/path/to/rcalsa-library \
#     Rscript scripts/run_ohie_tan_rcal_official.R

covariate_specification <- tolower(trimws(Sys.getenv(
  "IVMSM_OHIE_COVARIATE_SPEC", unset = "expanded"
)))
if (!covariate_specification %in% c("expanded", "augmented")) {
  stop("IVMSM_OHIE_COVARIATE_SPEC must be expanded or augmented.")
}
input_file <- trimws(Sys.getenv(
  "IVMSM_OHIE_INPUT_FILE",
  unset = if (covariate_specification == "augmented") {
    "data/ohie_application_augmented.rds"
  } else {
    "data/ohie_application.rds"
  }
))
if (!file.exists(input_file)) {
  stop(
    "Prepare the requested OHIE application data and execute this script ",
    "from the repository root. Missing: ", input_file,
    call. = FALSE
  )
}

read_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(
    name, unset = if (default) "true" else "false"
  )))
  if (!value %in% c("true", "false")) {
    stop(name, " must be true or false.", call. = FALSE)
  }
  identical(value, "true")
}

lambda_tag <- function(lambda) {
  value <- format(lambda, scientific = FALSE, trim = TRUE, digits = 12L)
  value <- sub("[.]0+$", "", value)
  value <- sub("([.][0-9]*?)0+$", "\\1", value)
  gsub(".", "p", value, fixed = TRUE)
}

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "-"), tmpdir = dirname(path)
  )
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically save ", path, ".", call. = FALSE)
  }
  invisible(path)
}

atomic_write_csv <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "-"), tmpdir = dirname(path)
  )
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(object, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically save ", path, ".", call. = FALSE)
  }
  invisible(path)
}

load_rcalsa <- function() {
  extra_libraries <- unique(c(
    Sys.getenv("RCALSA_LIBRARY", unset = ""),
    "/private/tmp/rcalsa-official-lib",
    "/tmp/RCALsa-library"
  ))
  extra_libraries <- extra_libraries[
    nzchar(extra_libraries) & dir.exists(extra_libraries)
  ]
  if (length(extra_libraries)) {
    .libPaths(c(extra_libraries, .libPaths()))
  }
  if (!requireNamespace("RCALsa", quietly = TRUE)) {
    stop(
      "RCALsa is required. Set RCALSA_LIBRARY to the library containing ",
      "the authors' RCALsa package and its RCAL and lpSolveAPI dependencies.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

build_application_dictionary <- function(
    main,
    minimum_nonzero_fraction = 0.008
) {
  main <- as.matrix(main)
  storage.mode(main) <- "double"
  original_names <- colnames(main)
  minimum_nonzero <- ceiling(minimum_nonzero_fraction * nrow(main))
  pairs <- if (ncol(main) > 1L) {
    utils::combn(ncol(main), 2L)
  } else {
    matrix(integer(), 2L, 0L)
  }
  interactions <- if (ncol(pairs)) {
    value <- vapply(seq_len(ncol(pairs)), function(index) {
      main[, pairs[1L, index]] * main[, pairs[2L, index]]
    }, numeric(nrow(main)))
    if (is.null(dim(value))) value <- matrix(value, ncol = 1L)
    colnames(value) <- paste0(
      original_names[pairs[1L, ]], ":", original_names[pairs[2L, ]]
    )
    retain <- colSums(value != 0) >= minimum_nonzero
    value[, retain, drop = FALSE]
  } else {
    matrix(numeric(), nrow(main), 0L)
  }
  # Tan's application library uses all main effects and two-way interactions;
  # it does not add squared terms. Sparse interactions are removed using the
  # same 0.8% nonzero-frequency rule used in the paper's application.
  full <- cbind(main, interactions)
  centers <- colMeans(full)
  scales <- apply(full, 2L, stats::sd)
  valid <- is.finite(scales) & scales > 0
  full <- sweep(full[, valid, drop = FALSE], 2L, centers[valid], "-")
  full <- sweep(full, 2L, scales[valid], "/")
  if (any(!is.finite(full))) {
    stop("Tan's working-model dictionary contains nonfinite values.")
  }
  # The added baseline indicators create redundant interaction columns. The
  # RCALsa optimization requires a full-rank working design, so remove only
  # columns that add no linear span after an intercept and all main effects
  # have entered. With the original 11-column specification this removes
  # nothing; with the compact augmented specification it removes five
  # redundant interactions. This is numerical de-duplication, not variable
  # selection.
  design_rank <- qr(
    cbind(`(Intercept)` = 1, full), tol = 1e-10, LAPACK = FALSE
  )
  independent <- sort(design_rank$pivot[seq_len(design_rank$rank)])
  retained_columns <- independent[independent > 1L] - 1L
  removed_columns <- setdiff(seq_len(ncol(full)), retained_columns)
  if (any(seq_len(ncol(main)) %in% removed_columns)) {
    stop("Rank screening attempted to remove an augmented main effect.")
  }
  removed_names <- colnames(full)[removed_columns]
  full <- full[, retained_columns, drop = FALSE]
  output <- list(
    x = full,
    main_terms = ncol(main),
    candidate_interactions = ncol(pairs),
    retained_interactions = ncol(full) - ncol(main),
    retained_terms = ncol(full),
    minimum_interaction_nonzero = minimum_nonzero,
    minimum_interaction_nonzero_fraction = minimum_nonzero_fraction
  )
  if (length(removed_columns)) {
    output$frequency_retained_interactions <- ncol(interactions)
    output$removed_linear_dependencies <- length(removed_columns)
    output$removed_linear_dependency_names <- removed_names
  }
  output
}

extract_diagnostics <- function(fit, lambda, bound) {
  components <- list(
    ps_t0 = fit$ps[[1L]], ps_t1 = fit$ps[[2L]],
    or_t0 = fit$or[[1L]], or_t1 = fit$or[[2L]],
    qr_t0 = fit$qr[[1L]], qr_t1 = fit$qr[[2L]]
  )
  do.call(rbind, lapply(names(components), function(component) {
    value <- components[[component]]
    selected_penalty <- if (!is.null(value$sel.rho)) {
      value$sel.rho[1L]
    } else if (!is.null(value$pen)) {
      value$pen[1L]
    } else {
      NA_real_
    }
    selected_nonzero <- if (!is.null(value$sel.nz)) {
      value$sel.nz[1L]
    } else if (!is.null(value$coef)) {
      sum(abs(value$coef[-1L]) > 0)
    } else {
      NA_real_
    }
    data.frame(
      Lambda = lambda,
      bound = bound,
      component = component,
      selected_penalty = selected_penalty,
      selected_nonzero = selected_nonzero,
      stringsAsFactors = FALSE
    )
  }))
}

valid_endpoint <- function(record, lambda, bound, adjustment, contract) {
  is.list(record) && identical(record$schema_version, 1L) &&
    isTRUE(all.equal(record$Lambda, lambda, tolerance = 1e-12)) &&
    identical(record$bound, bound) &&
    identical(record$adjustment, adjustment) &&
    identical(record$contract, contract) &&
    is.finite(record$estimate) && is.finite(record$standard_error) &&
    record$standard_error >= 0
}

read_endpoint <- function(path, lambda, bound, adjustment, contract) {
  if (!file.exists(path)) return(NULL)
  record <- tryCatch(readRDS(path), error = function(error) NULL)
  if (!valid_endpoint(record, lambda, bound, adjustment, contract)) NULL else record
}

collect_results <- function(root, lambdas, adjustment, contract) {
  endpoint_rows <- list()
  diagnostic_rows <- list()
  bound_rows <- list()
  critical <- stats::qnorm(0.975)
  for (lambda in lambdas) {
    directory <- file.path(root, paste0("lambda_", lambda_tag(lambda)))
    lower <- read_endpoint(
      file.path(directory, "endpoint_lower.rds"),
      lambda, "lower", adjustment, contract
    )
    upper <- read_endpoint(
      file.path(directory, "endpoint_upper.rds"),
      lambda, "upper", adjustment, contract
    )
    for (record in list(lower, upper)) {
      if (is.null(record)) next
      endpoint_rows[[length(endpoint_rows) + 1L]] <- data.frame(
        Lambda = record$Lambda,
        bound = record$bound,
        estimate = record$estimate,
        standard_error = record$standard_error,
        elapsed_seconds = record$elapsed_seconds,
        stringsAsFactors = FALSE
      )
      diagnostic_rows[[length(diagnostic_rows) + 1L]] <- record$diagnostics
    }
    if (!is.null(lower) && !is.null(upper)) {
      bound_rows[[length(bound_rows) + 1L]] <- data.frame(
        method = "Tan official RCAL",
        adjustment = adjustment,
        Lambda = lambda,
        ate_lower = lower$estimate,
        ate_upper = upper$estimate,
        se_lower = lower$standard_error,
        se_upper = upper$standard_error,
        confidence_lower = lower$estimate - critical * lower$standard_error,
        confidence_upper = upper$estimate + critical * upper$standard_error,
        stringsAsFactors = FALSE
      )
    }
  }
  bind <- function(rows) {
    if (!length(rows)) return(data.frame())
    output <- do.call(rbind, rows)
    rownames(output) <- NULL
    output
  }
  bounds <- bind(bound_rows)
  if (nrow(bounds)) {
    bounds$confidence_width <- bounds$confidence_upper - bounds$confidence_lower
    bounds <- bounds[order(bounds$Lambda), , drop = FALSE]
  }
  list(
    bounds = bounds,
    endpoints = bind(endpoint_rows),
    diagnostics = bind(diagnostic_rows)
  )
}

load_rcalsa()
prepared <- readRDS(input_file)
if (covariate_specification == "augmented") {
  data <- prepared$ed_household_index_overlap_augmented
  x_cols <- prepared$ed_overlap_augmented_x_cols
} else {
  data <- prepared$ed_household_index_overlap
  x_cols <- prepared$ed_overlap_expanded_x_cols
}
stopifnot(
  is.data.frame(data), nrow(data) == 22170L,
  all(c("Y", "T_obs", "Z", x_cols) %in% names(data)),
  !anyNA(data[, c("Y", "T_obs", "Z", x_cols)]),
  all(data$Y %in% 0:1), all(data$T_obs %in% 0:1)
)

adjustment <- tolower(Sys.getenv("IVMSM_OHIE_TAN_ADJUSTMENT", unset = "x"))
if (!adjustment %in% c("x", "xz")) {
  stop("IVMSM_OHIE_TAN_ADJUSTMENT must be x or xz.", call. = FALSE)
}
lambda_text <- Sys.getenv(
  "IVMSM_OHIE_TAN_LAMBDAS",
  unset = "1,1.2,1.25,1.4,1.5,1.6,1.75,1.8,2,2.5,3,3.5,4,4.5,5"
)
lambdas <- suppressWarnings(as.numeric(strsplit(
  lambda_text, ",", fixed = TRUE
)[[1L]]))
if (!length(lambdas) || any(!is.finite(lambdas)) || any(lambdas < 1)) {
  stop("All sensitivity parameters must be finite and at least one.")
}
lambdas <- unique(lambdas)
folds <- suppressWarnings(as.integer(Sys.getenv(
  "IVMSM_OHIE_TAN_FOLDS", unset = "5"
)))
if (length(folds) != 1L || is.na(folds) || folds < 2L) {
  stop("IVMSM_OHIE_TAN_FOLDS must be an integer of at least two.")
}
force <- read_flag("IVMSM_OHIE_TAN_FORCE")
seed <- 20260814L

raw_main <- as.matrix(data[, x_cols, drop = FALSE])
if (adjustment == "xz") {
  raw_main <- cbind(raw_main, instrument_Z = data$Z)
}
dictionary <- build_application_dictionary(raw_main)
x <- dictionary$x
default_version <- if (covariate_specification == "expanded") {
  paste0("ohie_ed_tan_rcal_official_", adjustment, "_v2")
} else {
  paste0("ohie_ed_tan_rcal_official_augmented_", adjustment, "_v1")
}
version <- trimws(Sys.getenv(
  "IVMSM_OHIE_TAN_VERSION", unset = default_version
))
if (!nzchar(version)) stop("IVMSM_OHIE_TAN_VERSION cannot be empty.")
output_root <- file.path("results", "ohie_application", version)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

contract <- c(
  list(data_md5 = unname(tools::md5sum(input_file))),
  if (covariate_specification == "augmented") {
    list(covariate_specification = covariate_specification)
  },
  list(
    adjustment = adjustment,
    x_cols = colnames(raw_main),
    observations = nrow(data),
    dictionary = dictionary[setdiff(names(dictionary), "x")],
    folds = rep(folds, 3L),
    nlam = rep(25L, 3L),
    lam_fac = rep(2^(-1 / 4), 3L),
    propensity_loss = "cal",
    outcome_loss = "ml",
    quantile_option = "qu",
    seed = seed,
    RCALsa_version = as.character(utils::packageVersion("RCALsa"))
  )
)
atomic_save_rds(contract, file.path(output_root, "analysis_contract.rds"))
atomic_write_csv(data.frame(
  adjustment = adjustment,
  observations = nrow(data),
  main_terms = dictionary$main_terms,
  candidate_interactions = dictionary$candidate_interactions,
  retained_interactions = dictionary$retained_interactions,
  retained_terms = dictionary$retained_terms,
  removed_linear_dependencies = if (
    is.null(dictionary$removed_linear_dependencies)
  ) 0L else dictionary$removed_linear_dependencies,
  folds = folds,
  cv_penalties = 25L,
  penalty_ratio = 2^(-1 / 4),
  RCALsa_version = as.character(utils::packageVersion("RCALsa")),
  stringsAsFactors = FALSE
), file.path(output_root, "configuration.csv"))

any_failure <- FALSE
for (lambda in lambdas) {
  directory <- file.path(output_root, paste0("lambda_", lambda_tag(lambda)))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  for (bound in c("lower", "upper")) {
    endpoint_path <- file.path(directory, paste0("endpoint_", bound, ".rds"))
    failure_path <- file.path(directory, paste0("failure_", bound, ".rds"))
    existing <- if (force) NULL else read_endpoint(
      endpoint_path, lambda, bound, adjustment, contract
    )
    if (!is.null(existing)) {
      if (file.exists(failure_path)) unlink(failure_path)
      cat(
        "Reusing Tan RCAL adjustment=", adjustment,
        ", Lambda=", lambda, ", endpoint=", bound, ".\n", sep = ""
      )
      next
    }
    cat(
      "Starting Tan RCAL adjustment=", adjustment,
      ", Lambda=", lambda, ", endpoint=", bound, ".\n", sep = ""
    )
    started <- Sys.time()
    set.seed(seed)
    fit <- tryCatch(
      RCALsa::sa.ate.regu.cv(
        fold = rep(folds, 3L),
        y = data$Y,
        tr = data$T_obs,
        x = x,
        gam = lambda,
        upper = identical(bound, "upper"),
        nlam = rep(25L, 3L),
        lam.fac = rep(2^(-1 / 4), 3L),
        ploss = "cal",
        yloss = "ml",
        option = "qu"
      ),
      error = identity
    )
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (inherits(fit, "error")) {
      any_failure <- TRUE
      atomic_save_rds(list(
        Lambda = lambda,
        bound = bound,
        adjustment = adjustment,
        message = conditionMessage(fit),
        elapsed_seconds = elapsed,
        contract = contract,
        session_info = utils::sessionInfo()
      ), failure_path)
      cat("FAILED: ", conditionMessage(fit), "\n", sep = "")
      next
    }
    estimate <- unname(fit$est["est", "diff"])
    standard_error <- sqrt(unname(fit$est["var", "diff"]))
    if (!is.finite(estimate) || !is.finite(standard_error) ||
        standard_error < 0) {
      stop("RCALsa returned an invalid endpoint or standard error.")
    }
    record <- list(
      schema_version = 1L,
      Lambda = lambda,
      bound = bound,
      adjustment = adjustment,
      estimate = estimate,
      standard_error = standard_error,
      elapsed_seconds = elapsed,
      completed_at = as.character(Sys.time()),
      diagnostics = extract_diagnostics(fit, lambda, bound),
      contract = contract,
      fit = fit
    )
    atomic_save_rds(record, endpoint_path)
    if (file.exists(failure_path)) unlink(failure_path)
    cat(
      "Completed Tan RCAL adjustment=", adjustment,
      ", Lambda=", lambda, ", endpoint=", bound,
      " (", round(elapsed, 1), " seconds).\n", sep = ""
    )
  }
  collected <- collect_results(output_root, lambdas, adjustment, contract)
  if (nrow(collected$bounds)) {
    atomic_write_csv(
      collected$bounds, file.path(output_root, "tan_rcal_bounds.csv")
    )
  }
  if (nrow(collected$endpoints)) {
    atomic_write_csv(
      collected$endpoints, file.path(output_root, "endpoint_progress.csv")
    )
  }
  if (nrow(collected$diagnostics)) {
    atomic_write_csv(
      collected$diagnostics, file.path(output_root, "tuning_diagnostics.csv")
    )
  }
}

collected <- collect_results(output_root, lambdas, adjustment, contract)
if (nrow(collected$bounds)) print(collected$bounds, row.names = FALSE)
cat("Results saved in: ", output_root, "\n", sep = "")
if (any_failure) {
  stop("At least one official RCALsa endpoint failed; inspect and rerun.")
}
