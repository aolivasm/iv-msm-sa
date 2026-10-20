# Create local configuration/integrity metadata, never simulation results.
prepare_experiments <- function(root = getwd()) {
  root <- normalizePath(root, mustWork = TRUE)
  configs <- dget(file.path(root, 'scripts/experiment_settings.R'))
  for (name in names(configs)) {
    path <- file.path(root, name)
    if (file.exists(path)) stopifnot(identical(readRDS(path), configs[[name]]))
    else saveRDS(configs[[name]], path)
  }
  manifests <- list.files(file.path(root, 'experiments'), pattern = '^run[.]R$',
                          recursive = TRUE, full.names = TRUE)
  bundles <- unique(dirname(dirname(manifests)))
  bundles <- bundles[order(nchar(bundles), decreasing = TRUE)]
  for (bundle in bundles) {
    files <- sort(list.files(bundle, pattern = '[.]R$|^study_config[.]rds$',
                            recursive = TRUE, full.names = TRUE))
    hashes <- tools::md5sum(files)
    names(hashes) <- substring(files, nchar(bundle) + 2L)
    path <- file.path(bundle, 'manifest.rds')
    if (file.exists(path)) {
      if (!identical(readRDS(path), hashes))
        stop('Source changed since preparation: ', bundle,
             '. Use a fresh checkout for a new run; do not mix checkpoints.')
    } else saveRDS(hashes, path)
  }
  invisible(TRUE)
}
if (sys.nframe() == 0L) {
  prepare_experiments()
  cat('Local configuration and source manifests prepared; no fits run.\n')
}
