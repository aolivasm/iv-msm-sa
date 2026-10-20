library(ivmsm)
set.seed(36)
d <- expand.grid(t=0:1,z=c('a','b','c'),rep=1:40)
d$x <- runif(nrow(d),-1,1); d$y <- .3*d$t+.2*d$x+rnorm(nrow(d),sd=.3)
f <- suppressWarnings(estimate_finite_iv_continuous(d,'y','t','z','x',folds=2,
  control=list(gamma_C=c(.1,1),gamma_rank=6,spline_df=2),seed=18))
g <- f$tuning[grepl('^gamma_',names(f$tuning))]
stopifnot(length(g)==8L)
for(v in g) {
  stopifnot(v$selected_bandwidth_multiplier==1/sqrt(2),
            abs(v$selected_bandwidth-v$base_bandwidth/sqrt(2))<1e-12,
            is.finite(v$selected_C),nrow(v$cv_resample_table)>0)
}
