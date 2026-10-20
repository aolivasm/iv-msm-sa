# GPFS account quota, not shared-filesystem free space. The parser is separate
# so small tests can exercise failure cases without querying the cluster.
n12_quota_bytes <- function(value) {
  if (length(value) != 1L || is.na(value) ||
      !grepl("^[0-9]+([.][0-9]+)?[KMGTPE]?$", value))
    stop("Unrecognized GPFS quota value: ", paste(value, collapse = ", "))
  unit <- sub("^[0-9.]+", "", value)
  number <- as.numeric(sub("[KMGTPE]$", "", value))
  result <- number * 1024^(match(unit, c("", "K", "M", "G", "T", "P", "E")) - 1L)
  if (!is.finite(result) || result < 0) stop("Invalid GPFS quota value: ", value)
  result
}

n12_parse_quota <- function(raw, filesystem = "mmfs1") {
  if (!is.character(raw) || !length(raw)) stop("Empty GPFS quota output.")
  status <- attr(raw, "status")
  if (!is.null(status) && !identical(as.integer(status), 0L))
    stop("GPFS quota query failed with status ", status, ".")
  rows <- grep("\\bUSR\\b", raw, value = TRUE)
  if (length(rows) != 1L) stop("Expected exactly one GPFS user quota row.")
  fields <- strsplit(trimws(rows), "[[:space:]]+")[[1L]]
  j <- match("USR", fields)
  # LPC includes a Fileset column (for example sas1-homeFS); older output
  # omits it. Match USR, accepting only these two documented row layouts.
  if (is.na(j) || !j %in% c(2L, 3L) || fields[1L] != filesystem ||
      length(fields) < j + 4L ||
      (j == 3L && !grepl("^[A-Za-z0-9_.-]+$", fields[2L])))
    stop("Unrecognized filesystem or incomplete GPFS account quota row.")
  values <- vapply(fields[j + 1:4], n12_quota_bytes, numeric(1))
  names(values) <- c("used_bytes", "soft_bytes", "hard_bytes", "in_doubt_bytes")
  limits <- values[c("soft_bytes", "hard_bytes")]
  limits <- limits[limits > 0]
  if (!length(limits)) stop("Cannot determine a finite GPFS account quota.")
  limit <- min(limits)
  available <- unname(limit - values["used_bytes"] - values["in_doubt_bytes"])
  if (!is.finite(available)) stop("Cannot determine GPFS account headroom.")
  c(as.list(values), list(limit_bytes = unname(limit), available_bytes = available))
}

n12_check_storage <- function(required_bytes = 8 * 1024^3) {
  if (length(required_bytes) != 1L || !is.numeric(required_bytes) ||
      !is.finite(required_bytes) || required_bytes < 0)
    stop("required_bytes must be one finite nonnegative number.")
  required_bytes <- max(required_bytes, 8 * 1024^3)
  tool <- "/usr/lpp/mmfs/bin/mmlsquota"
  if (!file.exists(tool) || file.access(tool, 1L) != 0L)
    stop("GPFS quota tool unavailable; run this check on LPC. No quota bypass is provided.")
  user <- unname(Sys.info()[["user"]])
  if (length(user) != 1L || is.na(user) || !grepl("^[A-Za-z0-9_.-]+$", user))
    stop("Cannot determine the current account for GPFS quota checking.")
  raw <- system2(tool, c("-u", user, "--block-size", "auto", "mmfs1"),
                 stdout = TRUE, stderr = TRUE)
  cat(paste(raw, collapse = "\n"), "\n")
  result <- n12_parse_quota(raw)
  result$required_bytes <- required_bytes
  cat(sprintf("GPFS account quota: %.3f GiB; used: %.3f GiB; in_doubt: %.3f GiB.\n",
              result$limit_bytes / 1024^3, result$used_bytes / 1024^3,
              result$in_doubt_bytes / 1024^3))
  cat(sprintf("Conservative headroom: %.3f GiB; required headroom: %.3f GiB.\n",
              result$available_bytes / 1024^3, required_bytes / 1024^3))
  if (result$available_bytes < required_bytes)
    stop("Insufficient GPFS account quota after in_doubt allocation; jobs are not released.")
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Usage: Rscript scripts/storage.R [required_bytes]")
  required <- if (length(args)) suppressWarnings(as.numeric(args[1L])) else 8 * 1024^3
  n12_check_storage(required)
}
