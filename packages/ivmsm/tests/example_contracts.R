library(ivmsm)
# Example data-flow checks only: no nuisance fitting or multiplier simulation.
files <- list.files(system.file('examples',package='ivmsm'),pattern='[.]R$',full.names=TRUE)
invisible(lapply(files,parse))
e <- new.env(parent=globalenv()); fits <- evaluations <- calibrations <- 0L
e$estimate_binary_sensitivity_band <- estimate_binary_sensitivity_band
environment(e$estimate_binary_sensitivity_band) <- e
e$print <- function(x,...)invisible(x)
e$plot <- function(x,...)invisible(x)
e$estimate_binary_dual_bounds <- function(data,...) {
  fits <<- fits+1L
  stopifnot(all(c('Y','T_obs','Z','X')%in%names(data)))
  list(nuisance_fit=list(marker='same training fits'),compatibility_action='project',compatibility_tolerance=1e-12)
}
e$evaluate_binary_bounds_from_nuisances <- function(data,nuisance_fit,Lambda0,Lambda1,epsilon,...) {
  evaluations <<- evaluations+1L
  stopifnot(identical(nuisance_fit$marker,'same training fits'),Lambda0==Lambda1,epsilon==0)
  v <- seq_len(nrow(data))/nrow(data)
  list(signals=list('0'=list(lower=v,upper=v+.2),
                   '1'=list(lower=v+Lambda0,upper=v+Lambda0+.4)),
       diagnostics=data.frame(checked=TRUE))
}
e$binary_grid_band <- function(grid,lower_signal,upper_signal,...) {
  calibrations <<- calibrations+1L
  stopifnot(nrow(lower_signal)==2000L,ncol(lower_signal)==21L,
    max(abs(lower_signal-matrix(rep(grid-.2,each=2000),2000)))<1e-12,
    max(abs(upper_signal-matrix(rep(grid+.4,each=2000),2000)))<1e-12)
  list(table=data.frame(Lambda=grid,lower=grid-.2,upper=grid+.4,
    simultaneous_lower=grid-.3,simultaneous_upper=grid+.5))
}
e$binary_continuum_band <- function(...)data.frame(Lambda=1.2,lower=-.1,upper=.5,empty=FALSE)
invisible(capture.output(sys.source(system.file('examples/binary_bands.R',package='ivmsm'),e)))
stopifnot(fits==1L,evaluations==21L,calibrations==1L,length(e$fit$diagnostics)==21L)
for(name in c('finite_iv.R','ordered_outcome.R','binary_bands.R')) {
  text <- readLines(system.file('examples',name,package='ivmsm'))
  stopifnot(!any(grepl('finite_iv_bounds\\(|discrete_arm_endpoint\\(|noise <-',text)))
}
cat('Installed examples parsed; fitted-band data flow passed with mocked learners and calibration.\n')
