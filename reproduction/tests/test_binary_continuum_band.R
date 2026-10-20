source('R/binary_continuum_band.R')
source('R/binary_outcome_bounds.R')
set.seed(20260920)
nu<-list(r0=runif(300,.1,.9),r1=runif(300,.1,.9),a0=runif(300,.1,.9),a1=runif(300,.1,.9))
# Algebra holds even if the resulting interval is empty; do not project nuisances.
for(L in c(1.01,1.2,1.6,2,3)) {
  f<-do.call(binary_flat_arm_bounds,c(nu,list(Lambda=L)))
  g<-binary_arm_endpoint(nu,Lambda0=L,Lambda1=L,epsilon=0)
  stopifnot(max(abs(f$lower-g$lower))<1e-12,max(abs(f$upper-g$upper))<1e-12)
}
grid<-seq(1.2,2,length.out=9);at<-seq(1.2,2,length.out=1001)
# A compatible, common-stratum observed law, heterogeneous across rows.
nu$r1<-nu$r0;nu$a1<-nu$a0
curve<-function(L) {
  f<-do.call(binary_flat_arm_bounds,c(nu,list(Lambda=L)))
  z<-nu;z$r0<-1-z$r0;z$r1<-1-z$r1;z$a0<-z$a1<-.25
  g<-do.call(binary_flat_arm_bounds,c(z,list(Lambda=L)))
  c(lower=mean(f$lower)-mean(g$upper),upper=mean(f$upper)-mean(g$lower))
}
truth<-t(vapply(at,curve,numeric(2)));gtruth<-t(vapply(grid,curve,numeric(2)))
lower<-gtruth[,1]-runif(length(grid),0,.05);upper<-gtruth[,2]+runif(length(grid),0,.05)
band<-binary_continuum_band(grid,lower,upper,at)
stopifnot(all(band$lower<=truth[,1]+1e-12),all(band$upper>=truth[,2]-1e-12),
  all(diff(band$lower)<=1e-12),all(diff(band$upper)>=-1e-12),
  all(abs(diff(band$lower))<=2*diff(at)+1e-12),all(abs(diff(band$upper))<=2*diff(at)+1e-12))
for(i in seq_along(at)) {
  right<-which(grid>=at[i])[1]
  stopifnot(band$lower[i]>=lower[right]-1e-12,band$upper[i]<=upper[right]+1e-12)
}
bad<-binary_continuum_band(c(1.2,2),c(.5,.5),c(-.5,-.5))
stopifnot(all(bad$empty)) # No silent midpoint repair of contradictory limits.
cat('Flat-envelope identity, continuum containment, nesting, Lipschitz and empty-set checks passed.\n')
