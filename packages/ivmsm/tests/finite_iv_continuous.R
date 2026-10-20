library(ivmsm)
set.seed(51)
d <- expand.grid(received=0:1, lottery=c('low','middle','high'), repeat_id=1:60)
d$age <- runif(nrow(d),-1,1)
d$response <- .5*d$received+.2*d$age+rnorm(nrow(d),sd=.4)
fast <- list(gamma_fixed_lambda=.1,gamma_rank=8,spline_df=2,
             threshold_maxit=10000L)
f <- suppressWarnings(estimate_finite_iv_continuous(d,'response','received','lottery',
   'age',Lambda=2,reference='low',folds=2,control=fast,seed=4))
stopifnot(nrow(f$endpoints)==6L,nrow(f$outward_sets)==3L,
  all(is.finite(f$signals)),length(f$fold_id)==nrow(d),
  f$instrument_levels[1]=='low',length(f$warnings)>=0,
  max(abs(f$signals[,'ate_lower']-(f$signals[,'mu1_lower']-f$signals[,'mu0_upper'])))<1e-12,
  max(abs(f$direct_plugin-f$pseudo_plugin-colMeans(f$direct_signals-f$pseudo_signals)))<1e-12)
for(i in seq_len(nrow(f$diagnostics))) {
  row <- f$diagnostics[i,]; k <- row$fold; key <- row$endpoint
  stopifnot(abs(row$direct_minus_pseudo-
    mean(f$direct_signals[f$fold_id==k,key]-f$pseudo_signals[f$fold_id==k,key]))<1e-10)
}
# Out-of-fold outcome perturbation cannot change that fold's fitted thresholds.
d2 <- d; d2$response[f$fold_id==1] <- d2$response[f$fold_id==1]+.1
g <- suppressWarnings(estimate_finite_iv_continuous(d2,'response','received','lottery',
  'age',Lambda=2,reference='low',folds=2,fold_id=f$fold_id,control=fast,seed=4))
for(key in paste0('threshold_1_',c('mu0_lower','mu0_upper','mu1_lower','mu1_upper')))
  stopifnot(identical(f$tuning[[key]],g$tuning[[key]]))
# No covariates and user-supplied fixed spline basis both work.
base <- suppressWarnings(estimate_finite_iv_continuous(d,'response','received','lottery',
  Lambda=2,folds=2,control=fast,seed=4,keep_signals=FALSE))
stopifnot(is.null(base$signals),base$threshold_dimension==1L)
fixed_spline <- function(x) cbind(1,splines::ns(x$age,knots=0,Boundary.knots=c(-1,1)))
sp <- suppressWarnings(estimate_finite_iv_continuous(d,'response','received','lottery',
  'age',Lambda=c(low=2,middle=2.2,high=2.4),reference='middle',folds=2,
  threshold_basis=fixed_spline,control=fast,seed=4))
stopifnot(sp$instrument_levels[1]=='middle',sp$Lambda[1]==2.2,sp$threshold_dimension==3L)
bad <- try(estimate_finite_iv_continuous(d,'response','received','lottery','age',
  folds=2,threshold_basis=function(x)cbind(1,x$age-mean(x$age)),control=fast),silent=TRUE)
stopifnot(inherits(bad,'try-error'))
