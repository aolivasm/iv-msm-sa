library(ivmsm)
source(system.file('examples/continuous_fit.R',package='ivmsm'))
stopifnot(inherits(f,'ivmsm_fit'),!is.null(f$endpoints))
if (requireNamespace('nnet',quietly=TRUE)) {
  source(system.file('examples/ordered_outcome.R',package='ivmsm'))
  stopifnot(is.list(fit),is.data.frame(fit$diagnostics))
}
