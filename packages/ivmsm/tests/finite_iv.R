library(ivmsm)
pi <- matrix(c(.3,.7),1,2)
p0 <- array(c(.66,.60,.34,.40),c(1,2,2))
p1 <- array(c(.55,.40,.45,.60),c(1,2,2))
v <- finite_iv_bounds(pi,p0,p1,0:1,1.6)
a <- binary_flat_arm_bounds(.3,.7,.45,.60,1.6)
b <- binary_flat_arm_bounds(.7,.3,.34,.40,1.6)
stopifnot(v$compatible,abs(v$estimates['ate_lower']-(a$lower-b$upper))<1e-12,
          abs(v$estimates['ate_upper']-(a$upper-b$lower))<1e-12)
bad <- finite_iv_bounds(pi,p0,p1,0:1,1)
stopifnot(!bad$compatible,all(is.na(bad$estimates)))
# Duplicating a category does not change the intersection.
idx <- c(1,2,2)
w <- finite_iv_bounds(pi[,idx,drop=FALSE],p0[,idx,,drop=FALSE],p1[,idx,,drop=FALSE],0:1,1.6)
stopifnot(max(abs(v$estimates-w$estimates))<1e-12)
# Multivalued support agrees with existing binary-IV water-filling implementation.
p <- array(rep(c(.2,.3,.5),each=2),c(1,2,3))
f <- finite_iv_bounds(pi,p,p,0:2,2)
old <- discrete_arm_endpoint(list(r0=.3,r1=.7,p0=matrix(c(.2,.3,.5),1),
  p1=matrix(c(.2,.3,.5),1)),support=0:2,Lambda0=2)
stopifnot(abs(f$estimates['mu1_lower']-old$exact_lower)<1e-10,
          abs(f$estimates['mu1_upper']-old$exact_upper)<1e-10)
d <- expand.grid(Y=0:2,T=0:1,Z=letters[1:4],rep=1:3)
stopifnot(fit_finite_iv_bounds(d,'Y','T','Z',0:2)$compatible)
stopifnot(inherits(try(finite_iv_bounds(pi,p0,p1,1:0),silent=TRUE),'try-error'))
