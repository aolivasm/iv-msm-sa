# Adapter for the authors' released Dorn--Guo--Kallus DVDS implementation.
#
# This file intentionally does not reimplement the DVDS estimating equations.
# It loads simplified_dvds.R from the authors' replication repository and
# calls dvds.binary() or dvds.continuous() from that file.  The nuisance fits
# below reproduce the GRF fitting pattern in the authors' dvds_sim.R script:
# one probability forest for the treatment propensity and separate arm-wise
# probability/regression forests for the outcome means.  In-arm predictions
# are GRF out-of-bag predictions; predictions for the other arm are ordinary
# forest predictions.
#
# Authors' repository:
#   https://github.com/CausalML/DoublyValidDoublySharpSensitivity
# Audited commit:
#   26cbbf6e90dcd0648c837c0d58568375291dd522

dvds_official_audited_commit <-
  "26cbbf6e90dcd0648c837c0d58568375291dd522"

dvds_official_required_packages <- function() {
  c("dplyr", "progress", "ggplot2", "grf")
}

dvds_official_check_packages <- function() {
  missing <- dvds_official_required_packages()[!vapply(
    dvds_official_required_packages(), requireNamespace, logical(1),
    quietly = TRUE
  )]
  if (length(missing)) {
    stop(
      "The authors' DVDS implementation requires: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

dvds_official_source_candidates <- function() {
  unique(c(
    trimws(Sys.getenv("IVMSM_DVDS_OFFICIAL_SOURCE", unset = "")),
    file.path("external", "DoublyValidDoublySharpSensitivity"),
    file.path("results", "official_sources",
              "DoublyValidDoublySharpSensitivity"),
    "/tmp/iv-msm-dvds-official-audit"
  ))
}

dvds_official_locate_source <- function(source_directory = NULL) {
  candidates <- if (is.null(source_directory)) {
    dvds_official_source_candidates()
  } else {
    as.character(source_directory)
  }
  candidates <- candidates[nzchar(candidates)]
  valid <- candidates[file.exists(file.path(candidates, "simplified_dvds.R"))]
  if (!length(valid)) {
    stop(
      "Could not find the authors' simplified_dvds.R. Clone the official ",
      "repository and set IVMSM_DVDS_OFFICIAL_SOURCE to its directory.",
      call. = FALSE
    )
  }
  normalizePath(valid[1L], mustWork = TRUE)
}

dvds_official_git_commit <- function(source_directory) {
  value <- tryCatch(
    suppressWarnings(system2(
      "git",
      c("-C", shQuote(source_directory), "rev-parse", "HEAD"),
      stdout = TRUE, stderr = FALSE
    )),
    error = function(error) character()
  )
  if (length(value) == 1L && grepl("^[0-9a-f]{40}$", value)) value else NA_character_
}

dvds_official_load <- function(
    source_directory = NULL,
    required_commit = dvds_official_audited_commit,
    strict_commit = TRUE
) {
  dvds_official_check_packages()
  source_directory <- dvds_official_locate_source(source_directory)
  commit <- dvds_official_git_commit(source_directory)
  if (isTRUE(strict_commit) &&
      (is.na(commit) || !identical(commit, required_commit))) {
    stop(
      "The DVDS source is not at the audited commit ", required_commit,
      "; observed commit: ", ifelse(is.na(commit), "unavailable", commit),
      ". Set IVMSM_DVDS_OFFICIAL_STRICT_COMMIT=false only after auditing ",
      "the new source.",
      call. = FALSE
    )
  }
  authors <- new.env(parent = globalenv())
  sys.source(
    file.path(source_directory, "simplified_dvds.R"),
    envir = authors,
    chdir = TRUE
  )
  required_functions <- c("dvds.binary", "dvds.continuous")
  missing <- required_functions[!vapply(
    required_functions, exists, logical(1),
    envir = authors, mode = "function", inherits = FALSE
  )]
  if (length(missing)) {
    stop(
      "The authors' source did not define: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  list(
    environment = authors,
    source_directory = source_directory,
    source_file = file.path(source_directory, "simplified_dvds.R"),
    git_commit = commit,
    required_commit = required_commit,
    strict_commit = isTRUE(strict_commit)
  )
}

dvds_official_validate_binary <- function(value, name) {
  observed <- sort(unique(as.numeric(value)))
  if (anyNA(observed) || !length(observed) ||
      !all(observed %in% c(0, 1))) {
    stop(name, " must be coded 0/1.", call. = FALSE)
  }
  invisible(TRUE)
}

dvds_official_prediction_matrix <- function(prediction) {
  value <- prediction$predictions
  # The propensity line in the released dvds_sim.R uses $prediction while all
  # other lines use the current GRF field $predictions.  Accept both solely to
  # make the adapter work with the GRF version used by the authors.
  if (is.null(value)) value <- prediction$prediction
  if (is.null(value)) {
    stop("GRF did not return a prediction component.", call. = FALSE)
  }
  value
}

dvds_official_probability_one <- function(prediction) {
  value <- dvds_official_prediction_matrix(prediction)
  if (is.null(dim(value))) return(as.numeric(value))
  if ("1" %in% colnames(value)) return(as.numeric(value[, "1"]))
  if (ncol(value) == 2L) return(as.numeric(value[, 2L]))
  stop("GRF did not return the probability of class 1.", call. = FALSE)
}

dvds_official_numeric_prediction <- function(prediction) {
  value <- dvds_official_prediction_matrix(prediction)
  if (!is.null(dim(value))) value <- value[, 1L]
  as.numeric(value)
}

dvds_official_fit_nuisances <- function(x, treatment, outcome, outcome_type) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  treatment <- as.integer(treatment)
  outcome <- as.numeric(outcome)
  outcome_type <- match.arg(outcome_type, c("continuous", "binary"))
  dvds_official_validate_binary(treatment, "treatment")
  if (outcome_type == "binary") {
    dvds_official_validate_binary(outcome, "outcome")
  }
  if (nrow(x) != length(treatment) || length(outcome) != length(treatment) ||
      any(!is.finite(x)) || any(!is.finite(outcome))) {
    stop("The DVDS analysis variables have incompatible dimensions or values.")
  }
  if (any(tabulate(treatment + 1L, nbins = 2L) < 2L)) {
    stop("Each treatment arm needs at least two observations for GRF.")
  }

  propensity_fit <- grf::probability_forest(
    X = x, Y = factor(treatment, levels = c(0, 1))
  )
  propensity <- dvds_official_probability_one(predict(propensity_fit))

  mu <- matrix(NA_real_, nrow(x), 2L,
               dimnames = list(NULL, c("mu0", "mu1")))
  outcome_fits <- vector("list", 2L)
  for (arm in 0:1) {
    in_arm <- treatment == arm
    if (outcome_type == "binary") {
      if (length(unique(outcome[in_arm])) < 2L) {
        stop(
          "The authors' probability-forest nuisance fit requires both ",
          "outcome classes in treatment arm ", arm, "."
        )
      }
      fit <- grf::probability_forest(
        X = x[in_arm, , drop = FALSE],
        Y = factor(outcome[in_arm], levels = c(0, 1))
      )
      predict_value <- function(newdata = NULL) {
        prediction <- if (is.null(newdata)) {
          predict(fit)
        } else {
          predict(fit, newdata)
        }
        dvds_official_probability_one(prediction)
      }
    } else {
      fit <- grf::regression_forest(
        X = x[in_arm, , drop = FALSE], Y = outcome[in_arm]
      )
      predict_value <- function(newdata = NULL) {
        prediction <- if (is.null(newdata)) {
          predict(fit)
        } else {
          predict(fit, newdata)
        }
        dvds_official_numeric_prediction(prediction)
      }
    }
    mu[in_arm, arm + 1L] <- predict_value()
    mu[!in_arm, arm + 1L] <- predict_value(x[!in_arm, , drop = FALSE])
    outcome_fits[[arm + 1L]] <- fit
  }
  if (any(!is.finite(propensity)) || any(!is.finite(mu))) {
    stop("The official GRF nuisance fits produced non-finite predictions.")
  }
  list(
    propensity = propensity,
    mu0 = mu[, "mu0"],
    mu1 = mu[, "mu1"],
    propensity_fit = propensity_fit,
    outcome_fits = outcome_fits,
    diagnostics = data.frame(
      component = c("propensity", "mu0", "mu1"),
      minimum = c(min(propensity), min(mu[, "mu0"]), min(mu[, "mu1"])),
      maximum = c(max(propensity), max(mu[, "mu0"]), max(mu[, "mu1"])),
      mean = c(mean(propensity), mean(mu[, "mu0"]), mean(mu[, "mu1"])),
      sd = c(stats::sd(propensity), stats::sd(mu[, "mu0"]),
             stats::sd(mu[, "mu1"])),
      stringsAsFactors = FALSE
    )
  )
}

dvds_official_bounds_to_endpoints <- function(bounds, method, adjustment) {
  critical <- stats::qnorm(0.975)
  lower_se <- (bounds$lower - bounds$lower.CI) / critical
  upper_se <- (bounds$upper.CI - bounds$upper) / critical
  output <- rbind(
    data.frame(
      Lambda = bounds$Lambda, endpoint = "ate_lower",
      estimate = bounds$lower, standard_error = lower_se,
      outward_confidence_endpoint = bounds$lower.CI
    ),
    data.frame(
      Lambda = bounds$Lambda, endpoint = "ate_upper",
      estimate = bounds$upper, standard_error = upper_se,
      outward_confidence_endpoint = bounds$upper.CI
    )
  )
  output$method <- method
  output$adjustment <- adjustment
  output <- output[order(output$Lambda, output$endpoint), , drop = FALSE]
  rownames(output) <- NULL
  output
}

estimate_dvds_official_authors <- function(
    data,
    x_cols,
    outcome_type = c("continuous", "binary"),
    Lambdas = c(2, 2.5, 3, 4, 5),
    treatment_col = "T_obs",
    outcome_col = "Y",
    adjustment = c("x", "xz"),
    authors_source = NULL,
    source_directory = NULL,
    strict_commit = TRUE,
    keep_fits = FALSE
) {
  outcome_type <- match.arg(outcome_type)
  adjustment <- match.arg(adjustment)
  Lambdas <- sort(unique(as.numeric(Lambdas)))
  if (!length(Lambdas) || any(!is.finite(Lambdas)) || any(Lambdas < 1)) {
    stop("Lambdas must be finite and at least one.")
  }
  required <- unique(c(x_cols, treatment_col, outcome_col))
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Missing DVDS analysis columns: ", paste(missing, collapse = ", "))
  }
  authors_source <- if (is.null(authors_source)) {
    dvds_official_load(
      source_directory = source_directory,
      strict_commit = strict_commit
    )
  } else authors_source
  x <- as.matrix(data[, x_cols, drop = FALSE])
  treatment <- as.integer(data[[treatment_col]])
  outcome <- as.numeric(data[[outcome_col]])
  nuisance <- dvds_official_fit_nuisances(
    x, treatment, outcome, outcome_type
  )

  method <- paste0("dorn_guo_kallus_dvds_official_authors_", adjustment)
  authors <- authors_source$environment
  if (outcome_type == "continuous") {
    # This is one direct call to the authors' vector-Lambda function.  Its ATE
    # columns are correct for each Lambda; the released arm-summary columns
    # average across Lambda and are therefore deliberately not reported here.
    bounds <- get("dvds.continuous", envir = authors, inherits = FALSE)(
      X = x, Y = outcome, Z = treatment,
      mu0 = nuisance$mu0, mu1 = nuisance$mu1,
      e = nuisance$propensity, Lambdas = Lambdas, alpha = 0.05
    )
  } else {
    binary_function <- get("dvds.binary", envir = authors, inherits = FALSE)
    rows <- lapply(Lambdas, function(Lambda) {
      binary_function(
        Y = outcome, Z = treatment,
        mu0 = nuisance$mu0, mu1 = nuisance$mu1,
        e = nuisance$propensity, Lambda = Lambda, alpha = 0.05
      )
    })
    bounds <- do.call(rbind, rows)
    rownames(bounds) <- NULL
  }
  required_bounds <- c("Lambda", "lower", "upper", "lower.CI", "upper.CI")
  if (!all(required_bounds %in% names(bounds)) ||
      any(!is.finite(as.matrix(bounds[, required_bounds]))) ||
      any(bounds$lower > bounds$upper)) {
    stop("The authors' DVDS function returned invalid ATE bounds.")
  }
  output <- list(
    bounds = bounds,
    endpoint_detail = dvds_official_bounds_to_endpoints(
      bounds, method, adjustment
    ),
    diagnostics = nuisance$diagnostics,
    provenance = data.frame(
      implementation = "authors released simplified_dvds.R",
      repository = paste0(
        "https://github.com/CausalML/",
        "DoublyValidDoublySharpSensitivity"
      ),
      source_directory = authors_source$source_directory,
      source_file = authors_source$source_file,
      git_commit = authors_source$git_commit,
      audited_commit = authors_source$required_commit,
      strict_commit = authors_source$strict_commit,
      nuisance_recipe = "authors dvds_sim.R GRF recipe",
      propensity_truncation = "none",
      stringsAsFactors = FALSE
    )
  )
  if (isTRUE(keep_fits)) output$nuisance_fits <- nuisance
  output
}
