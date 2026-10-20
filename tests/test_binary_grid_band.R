source('R/binary_grid_band.R')
source('R/binary_continuum_band.R')
set.seed(12)
lo <- matrix(rnorm(300), 100, 3)
hi <- lo + matrix(rnorm(300), 100, 3) + 1
saved <- .Random.seed
fit <- binary_grid_band(c(1.2, 1.6, 2), lo, hi, draws = 500, seed = 17)
stopifnot(identical(saved, .Random.seed),
          identical(fit, binary_grid_band(c(1.2, 1.6, 2), lo, hi, draws = 500, seed = 17)))
# Independent direct multiplier construction verifies signs and covariance.
S <- cbind(lo[,1], hi[,1], lo[,2], hi[,2], lo[,3], hi[,3])
V <- cor(S) * outer(rep(c(1,-1),3), rep(c(1,-1),3)) * 99/100
e <- eigen(V, symmetric = TRUE)
A <- e$vectors %*% diag(sqrt(pmax(e$values, 0)))
set.seed(17)
critical <- unname(quantile(apply(A %*% matrix(rnorm(6*500),6),2,max),.95,type=8))
stopifnot(abs(fit$critical-critical) < 1e-12,
          max(abs(fit$table$se_lower-apply(lo,2,sd)/10)) < 1e-12,
          all(fit$table$simultaneous_lower < colMeans(lo)),
          all(fit$table$simultaneous_upper > colMeans(hi)))
bad <- try(binary_grid_band(c(1.2,1.6,2),lo*0,hi,draws=500), silent=TRUE)
stopifnot(inherits(bad,'try-error'))
cat('Binary grid calibration, signs, standard errors, reproducibility, and zero-variance guard passed.\n')
