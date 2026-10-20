# Explicit allowlist: never discover/run every simulation or test automatically.
if (length(commandArgs(trailingOnly = TRUE))) stop('This script takes no arguments.')
files <- c('tests/check_binary_population_identities.R', 'tests/test_binary_continuum_band.R',
           'tests/test_binary_grid_band.R')
if (!all(file.exists(files))) stop('Run from the reproduction directory; smoke-test files are missing.')
rscript <- file.path(R.home('bin'), 'Rscript')
for (file in files) {
  cat('\nRunning:', file, '\n')
  status <- system2(rscript, c('--vanilla', shQuote(file)))
  if (status != 0L) stop('Smoke test failed: ', file, ' (exit ', status, ')')
}
cat('\nSmoke tests passed. Full fitted simulations and application inference were not run.\n')
