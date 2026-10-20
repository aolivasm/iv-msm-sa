#!/usr/bin/env Rscript
# RKHS-only bands. Run from the release root after fitting OHIE nuisances.
# Saves aggregate estimates/diagnostics only, never person-level signals.
source('R/binary_outcome_bounds.R')
source('R/binary_grid_band.R')
source('R/binary_continuum_band.R')

ohie_rkhs_bands <- function(data, checkpoint, grid=(120:200)/100,
                            draws=200000L, seed=20260920L) {
  stopifnot(identical(checkpoint$contract$person_id,data$person_id),
            !anyDuplicated(data$person_id), !anyDuplicated(data$household_id),
            length(data$person_id)==nrow(data), length(data$household_id)==nrow(data),
            length(grid)>=2L, all(diff(grid)>0), all(grid>=1))
  lower <- upper <- matrix(NA_real_,nrow(data),length(grid))
  diagnostics <- vector('list',length(grid))
  for (j in seq_along(grid)) {
    f <- evaluate_binary_bounds_from_nuisances(data,checkpoint$nuisance_fit,
      Lambda0=grid[j],Lambda1=grid[j],epsilon=0,compatibility_action='project')
    lower[,j] <- f$signals[['1']]$lower-f$signals[['0']]$upper
    upper[,j] <- f$signals[['1']]$upper-f$signals[['0']]$lower
    diagnostics[[j]] <- transform(f$diagnostics,Lambda=grid[j])
  }
  band <- binary_grid_band(grid,lower,upper,draws=draws,seed=seed)
  band$table$pointwise_lower <- band$table$lower-qnorm(.975)*band$table$se_lower
  band$table$pointwise_upper <- band$table$upper+qnorm(.975)*band$table$se_upper
  band$diagnostics <- do.call(rbind,diagnostics)
  band$interpretation <- paste('Separate 95% outward band for the binary ATE.',
    'Independent household-index adults; requires model and uniform regularity conditions.',
    'Projection diagnostics are retained and are not compatibility tests.',
    'Use binary_continuum_band for conservative between-grid extension.')
  band
}

if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if (length(args)!=3L) stop(paste('Usage: Rscript scripts/run_ohie_rkhs_bands.R',
    'PREPARED_DATA.rds NUISANCE_CHECKPOINT.rds NEW_OUTPUT_DIRECTORY'))
  if (dir.exists(args[3]) || file.exists(args[3])) stop('Output already exists; refusing overwrite.')
  d <- readRDS(args[1])$ed_household_index_overlap
  checkpoint <- readRDS(args[2])
  ans <- ohie_rkhs_bands(d,checkpoint)
  dir.create(args[3],recursive=TRUE)
  saveRDS(ans,file.path(args[3],'rkhs_grid_band.rds'))
  write.csv(ans$table,file.path(args[3],'rkhs_grid_band.csv'),row.names=FALSE)
  write.csv(ans$diagnostics,file.path(args[3],'diagnostics.csv'),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(args[3],'sessionInfo.txt'))
  cat('Saved aggregate RKHS band and diagnostics. No person-level output saved.\n')
}
