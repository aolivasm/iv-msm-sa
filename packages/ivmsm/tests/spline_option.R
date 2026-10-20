library(ivmsm)
make <- get('.fiv_spline_basis',asNamespace('ivmsm'))
x <- data.frame(X=seq(-.9,.9,length.out=20),W=seq(1,3,length.out=20),B=rep(0:1,10))
spec <- list(X=list(knots=0,boundary=c(-1,1)),W=list(knots=2,boundary=c(0,4)),B=list(linear=TRUE))
b <- make(names(x),spec)
expected <- cbind(1,splines::ns(x$X,knots=0,Boundary.knots=c(-1,1)),
 splines::ns(x$W,knots=2,Boundary.knots=c(0,4)),x$B)
stopifnot(isTRUE(all.equal(unname(b(x)),unname(expected))),
 max(abs(b(x)[1:8,]-b(x[1:8,])))<1e-12)
r <- make(names(x),reference=x)
stopifnot(max(abs(r(x)[1:8,]-r(x[1:8,])))<1e-12,
 inherits(try(make(names(x)),silent=TRUE),'try-error'),
 inherits(try(make('X',list(X=list(knots=2,boundary=c(-1,1)))),silent=TRUE),'try-error'))
cat('Spline option: additive basis, fixed rowwise map, reference freezing and validation passed; no fits.\n')
