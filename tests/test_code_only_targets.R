# Small deterministic integration fixture. No estimators or random datasets.
old<-getwd()
setwd('experiments/binary_bounds');source('scripts/common.R');setwd(old)
b<-list(config=config)
expressions<-as.list(parse('scripts/compute_paper_targets.R'))
for(expr in expressions) {
  if(is.call(expr) && identical(expr[[1]],as.name('local')))eval(expr)
  if(is.call(expr) && identical(expr[[1]],as.name('<-')) &&
     identical(expr[[2]],as.name('continuous_msm')))eval(expr)
}
d<-continuous_msm(64L,0L)
stopifnot(nrow(d)==16L,all(is.finite(d$lower)),all(d$lower<d$upper))
for(j in 1:4){z<-d[d$setting==j,];stopifnot(all(diff(z$lower)<0),all(diff(z$upper)>0))}
cat('PASS: continuous population integration fixture and sensitivity nesting.\n')
