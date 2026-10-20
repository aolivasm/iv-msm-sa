library(ivmsm)
ns <- asNamespace('ivmsm')
bounds <- get('.biv_bounds',ns); score <- get('.biv_signal',ns)
project <- get('.biv_project',ns)
# Deterministic algebra checks, including category-specific sensitivity values.
r <- matrix(c(.2,.4,.6,.7,.3,.5),2,3)
p <- matrix(c(.3,.35,.4,.45,.5,.55),2,3); L <- c(1.5,2,3)
b <- bounds(r,p,L)
for(side in c('lower','upper')) {
 v <- b[[side]]
 for(i in 1:2) {
  j <- v$j[i]; h <- 1e-7; a <- r; a[i,j] <- a[i,j]+h
  stopifnot(abs((bounds(a,p,L)[[side]]$value[i]-v$value[i])/h-v$dr[i])<1e-6)
  a <- p; a[i,j] <- a[i,j]+h
  stopifnot(abs((bounds(r,a,L)[[side]]$value[i]-v$value[i])/h-v$dp[i])<1e-6)
 }
}
pop <- finite_iv_bounds(r,array(c(1-p,p),c(2,3,2)),
  array(c(1-p,p),c(2,3,2)),support=0:1,Lambda=L)
stopifnot(max(abs(pop$arms$mu1$endpoints-cbind(b$lower$value,b$upper$value)))<1e-12)
rr <- matrix(.5,2,3); pp <- matrix(rep(c(.1,.9,.5),each=2),2,3)
pr <- project(rr,pp,rep(2,3),'project')
stopifnot(pr$rate==1,pr$projection_rate==1,
 all(bounds(rr,pr$p,rep(2,3))$lower$value<=bounds(rr,pr$p,rep(2,3))$upper$value+1e-9),
 identical(project(rr,pp,rep(2,3),'none')$p,pp),
 inherits(try(project(rr,pp,rep(2,3),'stop'),silent=TRUE),'try-error'))

# Exercise the public orchestration with deterministic mocked nuisance fits.
# No fitted simulation or learner invocation takes place.
env <- new.env(parent=ns); seen <- list()
env$.fiv_primitives <- function(tr,J,x,c,seed) {
 seen[[length(seen)+1L]] <<- tr[[x]]
 predict <- function(d)list(e=matrix(1/J,nrow(d),J),pi=matrix(.5,nrow(d),J),
  instrument_clip_rate=0,treatment_clip_rate=0)
 list(predict=predict,tuning=list(mock=TRUE))
}
env$.biv_base <- function(tr,va,t,x,ctl,seed) {
 stopifnot(!any(tr[[x]]%in%va[[x]]))
 list(p=rep(.4+.1*t,nrow(va)),tuning=list(mock=TRUE))
}
env$.biv_gamma <- function(pd,pn,t,x,ctl,seed)
 list(predict=function(d)list(normalized=rep(1,nrow(d))),selected_lambda=.1)
f <- estimate_finite_iv_binary; environment(f) <- env
d <- expand.grid(person=1:12,T=0:1,Z=c('low','middle','high'))
d$X <- seq_len(nrow(d)); d$Y <- d$person%%2; folds <- (d$person%%2)+1
set.seed(17); before <- .Random.seed
fit <- f(d,'Y','T','Z','X',reference='low',folds=2,fold_id=folds,Lambda=2,
 control=list(gamma_fixed_lambda=.1))
stopifnot(identical(before,.Random.seed),nrow(fit$endpoints)==6,nrow(fit$outward_sets)==3,
 nrow(fit$diagnostics)==4,all(fit$diagnostics$projection_rate==0),length(seen)==2,
 identical(fit$instrument_levels,c('low','high','middle')),
 max(abs(fit$endpoints$estimate-colMeans(fit$signals)))<1e-12,
 max(abs(fit$endpoints$standard_error-apply(fit$signals,2,sd)/sqrt(nrow(d))))<1e-12)
expected <- matrix(NA_real_,nrow(d),4)
dd <- data.frame(Y=d$Y,T_obs=d$T,Z=match(d$Z,fit$instrument_levels))
for(t in 0:1) expected[,2*t+1:2] <- score(dd,t,matrix(1/3,nrow(d),3),
 matrix(.5,nrow(d),3),matrix(.4+.1*t,nrow(d),3),rep(2,3))
stopifnot(max(abs(expected-fit$signals[,1:4]))<1e-12)
named <- f(d,'Y','T','Z','X',reference='low',folds=2,fold_id=folds,
 Lambda=c(middle=3,low=2,high=4),control=list(gamma_fixed_lambda=.1),keep_signals=FALSE)
stopifnot(identical(named$Lambda,c(low=2,high=4,middle=3)),is.null(named$signals))
bad <- function(...)inherits(try(f(...),silent=TRUE),'try-error')
stopifnot(bad(d,'Y','T','Z','X',control=list(typo=1)),
 bad(d,'Y','T','Z','X',fold_id=rep(1,nrow(d))),
 bad(transform(d,Y=2),'Y','T','Z','X'),
 bad(d,'Y','T','Z','X',Lambda=c(foo=2,bar=2,baz=2)))
cat('Binary finite-IV formulas, derivatives, projection and mocked cross-fitting passed; no learners run.\n')

# Verify RKHS tuning dispatch without fitting a kernel model.
ge <- new.env(parent=ns); calls <- list()
capture <- function(...) {calls[[length(calls)+1L]] <<- list(...);list(selected_lambda=.1)}
ge$fit_gamma_rkhs <- capture; ge$fit_gamma_rkhs_fixed <- capture
g <- get('.biv_gamma',ns); environment(g) <- ge
ctl <- get('.fiv_control',ns)(list(),3)
pd <- data.frame(Y=c(0,1,0,1),T_obs=c(0,0,1,1),Z=c(0,1,0,1))
invisible(g(pd,list(),0,character(),ctl,1))
stopifnot(calls[[1]]$bandwidth_multipliers==1/sqrt(2),
 calls[[1]]$cv_repeats==2,calls[[1]]$cv_folds==3,
 calls[[1]]$selection_rule=='minimum',
 max(abs(calls[[1]]$lambda_grid(100)-ctl$gamma_C/10))<1e-12)
ctl$gamma_fixed_lambda <- .2
invisible(g(pd,list(),0,character(),ctl,1))
stopifnot(calls[[2]]$lambda==.2,calls[[2]]$bandwidth_multiplier==1/sqrt(2))
cat('Fixed M/2 bandwidth and repeated Riesz-CV dispatch verified with mocks.\n')
