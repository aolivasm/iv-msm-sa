# Resource-only wrapper for the frozen official DVDS nuisance adapter.
# Source the official adapter first, then explicitly install the returned copy:
#   dvds_official_fit_nuisances <- single_thread_dvds_nuisances(
#     dvds_official_fit_nuisances
#   )
# No source file, package namespace, forest tuning, or RNG seed is modified.

single_thread_dvds_nuisances <- function(fit_function) {
  if (!is.function(fit_function) || is.primitive(fit_function)) {
    stop("fit_function must be the sourced DVDS nuisance function.")
  }
  original_body <- body(fit_function)
  original_environment <- environment(fit_function)
  forest_heads <- list(
    probability_forest = quote(grf::probability_forest),
    regression_forest = quote(grf::regression_forest)
  )
  counts <- c(probability_forest = 0L, regression_forest = 0L)

  rewrite <- function(node, remove = FALSE) {
    if (!is.call(node)) return(node)
    matched <- names(forest_heads)[vapply(
      forest_heads, function(head) identical(node[[1L]], head), logical(1)
    )]
    if (length(matched)) {
      if (isTRUE(remove)) {
        if (!identical(node[["num.threads"]], 1L)) {
          stop("Unexpected thread argument while verifying the wrapper.")
        }
        node[["num.threads"]] <- NULL
      } else {
        if ("num.threads" %in% names(node)) {
          stop("The input adapter already sets num.threads; inspect its source.")
        }
        counts[[matched]] <<- counts[[matched]] + 1L
        node[["num.threads"]] <- 1L
      }
    }
    # In-place language-object replacement preserves source attributes and
    # argument names. Missing arguments (e.g. x[, column]) must remain missing.
    for (index in seq_along(node)) {
      if (!identical(node[[index]], quote(expr = ))) {
        node[index] <- list(rewrite(node[[index]], remove = remove))
      }
    }
    node
  }

  threaded_body <- rewrite(original_body)
  expected <- c(probability_forest = 2L, regression_forest = 1L)
  if (!identical(counts, expected)) {
    stop("The frozen DVDS nuisance adapter must contain exactly three forest calls.")
  }
  if (!identical(rewrite(threaded_body, remove = TRUE), original_body)) {
    stop("AST verification failed: changes extend beyond num.threads arguments.")
  }

  local_environment <- new.env(parent = original_environment)
  local_environment$predict <- function(object, ...) {
    arguments <- list(...)
    arguments$num.threads <- 1L
    do.call(stats::predict, c(list(object = object), arguments))
  }
  result <- fit_function
  body(result) <- threaded_body
  environment(result) <- local_environment
  attr(result, "single_thread_dvds") <- list(
    num.threads = 1L,
    forest_calls = counts,
    forest_ast_changes_only = TRUE,
    prediction_thread_override = TRUE
  )
  stopifnot(
    identical(formals(result), formals(fit_function)),
    identical(body(fit_function), original_body),
    identical(environment(fit_function), original_environment)
  )
  result
}

# Alternate descriptive entry point; also returns a function, never assigns it.
configure_dvds_single_thread <- single_thread_dvds_nuisances
