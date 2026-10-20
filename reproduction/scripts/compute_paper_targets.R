#!/usr/bin/env Rscript
# Deterministic numerical integration, not simulated datasets or fitted results.
# First generate the IV-MSM targets using paper.R experiment ... targets.
out <- 'results/population_targets.rds'
if (file.exists(out)) stop('Refusing to overwrite targets: ', out)
rows <- list()
add <- function(d, outcome, display) {
  rows[[length(rows)+1L]] <<- data.frame(outcome=outcome, setting=d$setting,
    Lambda=d$Lambda, display=display, target_lower=d$lower, target_upper=d$upper)
}
for (name in c('continuous_L2_L3','continuous_L18_L25')) {
  d <- readRDS(file.path('results/work',name,'results/population_targets.rds'))$targets
  add(subset(d,class=='fixed_spline5'),'Continuous','IV-MSM')
}
b <- readRDS('results/work/binary_bounds/targets/population.rds')
add(b$targets,'Binary','IV-MSM')
source('scripts/analysis_binary_x5_comparator_population_v1.R')
bc <- binary_x5_msm_integrate(11L,b$config)$targets
bf <- binary_x5_msm_integrate(15L,b$config)$targets
stopifnot(max(abs(bc$lower-bf$lower),abs(bc$upper-bf$upper)) < 1e-5)
for (method in c('Tan','Dorn')) add(subset(bf,adjustment=='x' & Lambda!=1),'Binary',method)

# Continuous observed outcome is theta'X + .75*T + a truncated-normal residual.
# In the X-adjusted MSM, mix residual distributions using P(Z | T=t,X).
local({
  old<-getwd();on.exit(setwd(old));setwd('experiments/continuous_bounds')
  source('R/x5_simple_population_targets.R',local=globalenv())
})
continuous_msm <- function(n, shift) {
  x <- x5p_rule(n,shift);xm<-as.matrix(x)
  h<-drop(xm%*%c(.14,.12,.10,.08,.06))
  location<-drop(xm%*%c(.20,-.25,.30,.15,-.20))
  shape<-1+.6*sin(x$X1*x$X2)+.4*sin(x$X3*x$X4)
  settings<-b$config$settings;answer<-list()
  for(j in 1:4)for(L in c(1.8,2,2.5,3)) {
    arms<-lapply(0:1,function(t) {
      r0<-plogis(h-settings$a[j]/2);r1<-plogis(h+settings$a[j]/2)
      if(t==0){r0<-1-r0;r1<-1-r1}
      pt<-(r0+r1)/2;w0<-r0/(r0+r1);eta<-(2*t-1)*settings$b[j]*shape
      moments<-function(q,eta) {
        den<-pnorm(2-eta)-pnorm(-2-eta)
        mass<-(pnorm(q-eta)-pnorm(-2-eta))/den
        first<-eta*mass+(dnorm(-2-eta)-dnorm(q-eta))/den
        list(mass=mass,first=first)
      }
      mix<-function(q){a<-moments(q,0);z<-moments(q,eta)
        list(mass=w0*a$mass+(1-w0)*z$mass,first=w0*a$first+(1-w0)*z$first)}
      mean_r<-mix(2)$first
      quant<-function(prob){lo<-rep(-2,n);hi<-rep(2,n)
        for(k in 1:45){mid<-(lo+hi)/2;below<-mix(mid)$mass<prob
          lo[below]<-mid[below];hi[!below]<-mid[!below]};(lo+hi)/2}
      tau<-1/(L+1)
      low<-mean_r/L+(L-1/L)*mix(quant(tau))$first
      high<-mean_r/L+(L-1/L)*(mean_r-mix(quant(1-tau))$first)
      c(lower=mean(location+.75*t+pt*mean_r+(1-pt)*low),
        upper=mean(location+.75*t+pt*mean_r+(1-pt)*high))
    })
    answer[[length(answer)+1L]]<-data.frame(setting=j,Lambda=L,
      lower=unname(arms[[2]]['lower']-arms[[1]]['upper']),
      upper=unname(arms[[2]]['upper']-arms[[1]]['lower']))
  }
  do.call(rbind,answer)
}
coarse<-continuous_msm(8192L,1729L);fine<-continuous_msm(16384L,9821L)
error<-max(abs(coarse$lower-fine$lower),abs(coarse$upper-fine$upper))
if(!is.finite(error)||error>=5e-4)stop('Continuous MSM target refinement failed: ',error)
for(method in c('Tan','Dorn'))add(fine,'Continuous',method)
targets<-do.call(rbind,rows)
stopifnot(!anyDuplicated(targets[c('outcome','setting','Lambda','display')]))
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
saveRDS(targets,out)
cat('Population targets computed locally; continuous MSM refinement:',error,'\n')
