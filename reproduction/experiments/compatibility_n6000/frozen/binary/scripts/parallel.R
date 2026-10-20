#!/usr/bin/env Rscript
# Independent Rscript workers, never forked estimator state. Source-only common.R.
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1,
           VECLIB_MAXIMUM_THREADS = 1, RCPP_PARALLEL_NUM_THREADS = 1,
           BLIS_NUM_THREADS = 1, NUMEXPR_NUM_THREADS = 1)
args <- commandArgs(TRUE)
stage <- if (length(args)) args[[1L]] else "run"
if (!stage %in% c("run", "archive")) stop("Use run or archive.")
worker_text <- Sys.getenv("BINARY_COMPAT_WORKERS", "8")
if (!grepl("^[1-8]$", worker_text)) stop("BINARY_COMPAT_WORKERS must be 1 through 8.")
workers <- as.integer(worker_text)
launch_dir <- Sys.getenv("BINARY_COMPAT_LAUNCH_DIR")
lock_dir <- Sys.getenv("BINARY_COMPAT_LOCK_DIR")
token <- Sys.getenv("BINARY_COMPAT_LAUNCH_TOKEN")
if (!nzchar(launch_dir) || !dir.exists(launch_dir) || !nzchar(lock_dir) ||
    !nzchar(token) || !file.exists(file.path(lock_dir, "token")) ||
    !identical(readLines(file.path(lock_dir, "token"), warn = FALSE), token)) {
  stop("Use bash scripts/launch_mac.sh; an owned launcher lock is required.")
}
drained_file <- file.path(launch_dir, "queue_drained.flag")
if (stage == "run") {
  if (file.exists(file.path(launch_dir, "progress.tsv"))) {
    stop("This invocation already has runner records; use the launcher for a fresh locked resume.")
  }
  writeLines("No worker has been launched by this invocation.", drained_file)
}
source("scripts/common.R")
if (!is.data.frame(tasks) || !"id" %in% names(tasks) || !nrow(tasks) ||
    anyNA(tasks$id) || anyDuplicated(tasks$id) ||
    any(tasks$id < 1 | tasks$id != as.integer(tasks$id)) ||
    !length(gate_ids) || anyNA(gate_ids) || anyDuplicated(gate_ids) ||
    !all(gate_ids %in% tasks$id)) stop("Invalid task or gate inventory.")
if (!is.character(config$result_root) || length(config$result_root) != 1L ||
    !grepl("^results/[A-Za-z0-9_-]+$", config$result_root)) {
  stop("config$result_root must be a simple bundle-relative results directory.")
}
result_root <- config$result_root
rscript <- file.path(R.home("bin"), "Rscript")
say <- function(...) {cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n"); flush.console()}
last_worker_line <- function(path) {
  tryCatch({
    size <- file.info(path)$size
    if (!is.finite(size) || size < 1) return("(no task output yet)")
    con <- file(path, "rb"); on.exit(close(con))
    seek(con, where = max(0, size - 4096), origin = "start")
    bytes <- readBin(con, "raw", n = 4096L)
    lines <- trimws(strsplit(rawToChar(bytes), "[\r\n]+")[[1L]])
    lines <- lines[nzchar(lines)]
    if (length(lines)) tail(lines, 1L) else "(no task output yet)"
  }, error = function(e) "(task log temporarily unavailable)")
}
free_space_check <- function() {
  output <- system2("df", c("-Pk", shQuote(getwd())), stdout = TRUE, stderr = TRUE)
  value <- suppressWarnings(as.numeric(strsplit(trimws(tail(output, 1L)), "[[:space:]]+")[[1L]][4L]))
  if (!is.null(attr(output, "status")) || !is.finite(value) || value < 2097152) {
    stop("Cannot confirm 2 GiB of available space; checkpoints and logs retained.")
  }
}

if (stage == "archive") {
  required <- c(file.path(result_root, c("summary.rds", "signature.rds")),
                "scripts", "R", "provenance", "README.txt", "manifest.rds",
                "study_config.rds")
  if (!all(file.exists(required))) stop("Collector/source artifacts missing; no archive created.")
  free_space_check()
  archive_dir <- file.path("archives", paste0("binary_compat_", token))
  if (!dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Archive directory already exists; refusing to overwrite.")
  }
  archive_file <- file.path(archive_dir, paste0("binary_x5_compatibility_", token, ".tar.gz"))
  launch_relative <- substring(normalizePath(launch_dir), nchar(normalizePath(getwd())) + 2L)
  if (!startsWith(launch_relative, "logs/launch_")) stop("Launch logs are outside this bundle.")
  # Full append-only logs stay local. This archive includes progress and 4 KiB tails.
  selected_logs <- list.files(launch_relative, full.names = TRUE)
  selected_logs <- selected_logs[basename(selected_logs) != "workers"]
  members <- c(required, selected_logs)
  optional <- file.path(result_root, "gate_report.rds")
  if (file.exists(optional)) members <- c(members, optional)
  Sys.setenv(COPYFILE_DISABLE = "1")
  status <- system2("tar", c("-czf", shQuote(archive_file), shQuote(members)))
  if (!identical(as.integer(status), 0L)) stop("Archive creation failed; partial archive retained.")
  status <- system2("gzip", c("-t", shQuote(archive_file)))
  if (!identical(as.integer(status), 0L)) stop("Archive integrity check failed; archive retained for inspection.")
  size <- file.info(archive_file)$size / 1024^2
  say("COMPLETE:", normalizePath(archive_file), sprintf("(%.2f MiB)", size))
  if (size > 20) say("Archive exceeds the 20 MiB compactness target; inspect before transferring.")
  say("Replication checkpoints and full worker logs remain local; compact tails are included.")
  quit(status = 0L)
}

worker_dir <- file.path(launch_dir, "workers")
tail_dir <- file.path(launch_dir, "worker_log_tails")
dir.create(worker_dir, showWarnings = FALSE)
dir.create(tail_dir, showWarnings = FALSE)
progress_file <- file.path(launch_dir, "progress.tsv")
writeLines("phase\tid\texit_code\telapsed_seconds\tlog", progress_file)
writeLines(c("Full worker logs remain in this launch's workers directory and are not archived.",
             "worker_log_tails contains at most the last 4096 bytes per task, not complete logs.",
             "No observation-level data, model fits, or replication checkpoints enter the archive."),
           file.path(launch_dir, "ARCHIVE_LOG_NOTES.txt"))
active <- list()
uncertain <- FALSE
completed <- data.frame(phase = character(), id = integer(), code = integer(), seconds = numeric())
is_alive <- function(pid) {
  identical(as.integer(suppressWarnings(system2("/bin/kill", c("-0", pid),
                                                 stdout = FALSE, stderr = FALSE))), 0L)
}
start_one <- function(id, phase) {
  stem <- file.path(worker_dir, sprintf("task_%05d", id))
  paths <- paste0(stem, c(".log", ".pid", ".status"))
  if (any(file.exists(paths))) stop("Worker artifacts already exist; refusing to overwrite: ", stem)
  log <- paths[[1L]]; pid <- paths[[2L]]; status <- paths[[3L]]
  # The R child PID is recorded so an interrupted controller cannot mistake an
  # orphan worker for an idle experiment. Exit status publication is atomic.
  command <- paste0(
    shQuote(rscript), " scripts/run.R task ", as.integer(id),
    " >> ", shQuote(log), " 2>&1 &\nworker=$!\n",
    "printf '%s\\n' \"$worker\" > ", shQuote(paste0(pid, ".tmp")), "\n",
    "mv ", shQuote(paste0(pid, ".tmp")), " ", shQuote(pid), "\n",
    "wait \"$worker\"\ncode=$?\n",
    "printf '%s\\n' \"$code\" > ", shQuote(paste0(status, ".tmp")), "\n",
    "mv ", shQuote(paste0(status, ".tmp")), " ", shQuote(status), "\nexit 0\n")
  if (file.exists(drained_file) && !file.remove(drained_file)) {
    stop("Could not clear this invocation's drainage marker before dispatch.")
  }
  code <- tryCatch(system2("/bin/sh", c("-c", shQuote(command)), wait = FALSE,
                           stdout = FALSE, stderr = log), error = function(e) {
    uncertain <<- TRUE
    stop("Could not confirm launch of task ", id, ": ", conditionMessage(e))
  })
  if (!identical(as.integer(code), 0L)) {
    uncertain <<- TRUE
    stop("Could not confirm launch of task ", id)
  }
  list(id = id, phase = phase, log = log, pid = pid, status = status,
       started = Sys.time(), dead_since = NULL)
}
retain_tail <- function(log) {
  if (!file.exists(log)) return(invisible(NULL))
  con <- file(log, "rb"); on.exit(close(con))
  seek(con, where = max(0, file.info(log)$size - 4096), origin = "start")
  bytes <- readBin(con, "raw", n = 4096L)
  writeBin(bytes, file.path(tail_dir, basename(log)))
}
poll <- function() {
  if (!length(active)) return(invisible(NULL))
  for (key in names(active)) {
    item <- active[[key]]
    code <- NA_integer_
    if (file.exists(item$status)) {
      value <- readLines(item$status, warn = FALSE)
      if (length(value) != 1L || !grepl("^[0-9]+$", value)) {
        uncertain <<- TRUE; code <- 255L
      } else code <- as.integer(value)
    } else if (file.exists(item$pid)) {
      pid <- readLines(item$pid, warn = FALSE)
      if (length(pid) != 1L || !grepl("^[1-9][0-9]*$", pid)) {
        uncertain <<- TRUE; code <- 255L
      } else if (!is_alive(pid)) {
        if (is.null(item$dead_since)) item$dead_since <- Sys.time()
        # Allow the supervising shell time to publish its exit status.
        if (as.numeric(difftime(Sys.time(), item$dead_since, units = "secs")) > 5) code <- 255L
      }
    } else if (as.numeric(difftime(Sys.time(), item$started, units = "secs")) > 30) {
      uncertain <<- TRUE; code <- 255L
    }
    active[[key]] <<- item
    if (!is.na(code)) {
      seconds <- as.numeric(difftime(Sys.time(), item$started, units = "secs"))
      retain_tail(item$log)
      cat(item$phase, item$id, code, round(seconds, 3), item$log,
          sep = "\t", file = progress_file, append = TRUE); cat("\n", file = progress_file, append = TRUE)
      completed <<- rbind(completed, data.frame(phase = item$phase, id = item$id, code = code, seconds = seconds))
      active[[key]] <<- NULL
      say(item$phase, "task", item$id, "exit", code, sprintf("after %.1f s;", seconds),
          nrow(completed), "/", nrow(tasks), "datasets finished in this invocation.")
    }
  }
  invisible(NULL)
}
drain <- function() {
  if (length(active)) say("No further dispatch; waiting for", length(active), "active worker(s).")
  suspendInterrupts({
    while (length(active)) {poll(); if (length(active)) Sys.sleep(1)}
    if (!uncertain) writeLines("All launched R workers have exited.", drained_file)
  })
}
main <- function() {
  on.exit(drain(), add = TRUE)
  for (phase in c("gate", "full")) {
    free_space_check()
    ids <- if (phase == "gate") gate_ids else setdiff(tasks$id, gate_ids)
    limit <- if (phase == "gate") min(2L, workers) else workers
    next_index <- 1L
    last_report <- Sys.time()
    say("Stage", phase, ":", length(ids), "datasets;", limit, "independent workers.")
    while (next_index <= length(ids) || length(active)) {
      poll()
      if (any(completed$code != 0L)) {
        drain()
        stop("Task failure: pending tasks were not launched; failed checkpoints and logs retained. Inspect before retrying.")
      }
      while (next_index <= length(ids) && length(active) < limit) {
        id <- ids[[next_index]]
        # Do not allow a keyboard interrupt between spawn and registration.
        suspendInterrupts(active[[as.character(id)]] <<- start_one(id, phase))
        next_index <- next_index + 1L
      }
      if (as.numeric(difftime(Sys.time(), last_report, units = "secs")) >= 60) {
        timings <- completed$seconds[completed$code == 0L & completed$seconds >= 5]
        remaining <- length(ids) - next_index + 1L + length(active)
        eta <- if (length(timings)) sprintf("rough remaining %.1f hours", remaining * median(timings) / limit / 3600) else "ETA unavailable until actual tasks finish"
        say(phase, ":", length(active), "active;", remaining, "remaining;", eta,
            "(resume skips and setting differences affect timing).")
        if (phase == "gate") {
          for (item in active) say("gate task", item$id, "latest:", last_worker_line(item$log))
        }
        last_report <- Sys.time()
      }
      if (length(active)) Sys.sleep(1)
    }
    if (phase == "gate") {
      code <- system2(rscript, c("scripts/run.R", "check_gate"))
      if (!identical(as.integer(code), 0L)) stop("Gate validation failed; full phase blocked.")
      say("GATE PASSED. Full dispatch uses the same sealed n6000 design and", workers, "workers.")
    }
  }
}
tryCatch(main(), interrupt = function(e) {
  say("Interrupted; active workers were drained, checkpoints retained.")
  quit(status = 130L)
})
