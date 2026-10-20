# Read-only dependency check. Never installs packages or accesses study data.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop('Usage: check_environment.R [minimal|fitting]')
profile <- if (length(args)) args[[1L]] else 'minimal'
if (!profile %in% c('minimal', 'fitting')) stop('Unknown profile: ', profile)
if (getRversion() < '4.4.0') stop('Use R 4.4 or later for these checks.')
cat(R.version.string, '\nProfile:', profile, '\n')
packages <- if (profile == 'minimal') character() else
  c('mgcv', 'glmnet', 'ranger', 'xgboost', 'SuperLearner', 'nnls',
    'nloptr', 'hal9001', 'haven', 'digest')
available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
if (length(packages)) {
  versions <- vapply(packages, function(p) {
    if (available[[p]]) as.character(packageVersion(p)) else 'MISSING'
  }, character(1))
  print(data.frame(package = packages, version = versions), row.names = FALSE)
}
if (any(!available)) stop('Missing packages: ', paste(packages[!available], collapse = ', '))
cat('Environment check passed. This does not verify frozen-run version parity.\n')
