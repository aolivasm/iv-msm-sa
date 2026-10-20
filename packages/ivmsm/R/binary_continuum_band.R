# Deterministic extension of simultaneous grid confidence limits.
# It does not create/repair the grid's statistical coverage guarantee.
binary_continuum_band <- function(grid, lower, upper, at=grid, lipschitz=2) {
  stopifnot(length(grid)>=2L,length(lower)==length(grid),length(upper)==length(grid),
    all(is.finite(c(grid,lower,upper,at,lipschitz))),all(diff(grid)>0),
    min(grid)>=1,all(at>=min(grid)&at<=max(grid)),length(lipschitz)==1L,lipschitz>=2)
  # 2 is a proved universal bound for binary IV-MSM ATE endpoints along Lambda0=Lambda1.
  # Do not substitute a smaller estimated constant without an appropriate guarantee.
  penalty <- lipschitz*pmax(outer(at,grid,'-'),0)
  lo <- apply(sweep(-penalty,2,lower,'+'),1,max)
  hi <- apply(sweep(penalty,2,upper,'+'),1,min)
  data.frame(Lambda=at,lower=lo,upper=hi,empty=lo>hi)
}

binary_flat_arm_bounds <- function(r0,r1,a0,a1,Lambda) {
  stopifnot(length(Lambda)==1,Lambda>=1,is.finite(Lambda))
  lower<-cbind(r0*a0+(1-r0)*a0/Lambda,
    r0*a0+(1-r0)*(1-(1-a0)*Lambda),
    r1*a1+(1-r1)*a1/Lambda,
    r1*a1+(1-r1)*(1-(1-a1)*Lambda))
  upper<-cbind(r0*a0+(1-r0)*a0*Lambda,
    r0*a0+(1-r0)*(1-(1-a0)/Lambda),
    r1*a1+(1-r1)*a1*Lambda,
    r1*a1+(1-r1)*(1-(1-a1)/Lambda))
  list(lower=apply(lower,1,max),upper=apply(upper,1,min))
}
