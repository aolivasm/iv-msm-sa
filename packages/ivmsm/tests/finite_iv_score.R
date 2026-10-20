library(ivmsm)
component <- getFromNamespace('.fiv_components','ivmsm')
score <- getFromNamespace('.fiv_score','ivmsm')
# Enumerated observed law: verify analytic fixed-threshold score against tilts.
J <- 3L; support <- c(-1,-.2,.6,1.5); K <- length(support)
cells <- expand.grid(y=seq_len(K),t=0:1,z=seq_len(J))
mass <- exp(.13*cells$z+.2*cells$t-.3*support[cells$y]^2+
            .09*cells$z*support[cells$y]+.1*cells$t*support[cells$y])
mass <- mass/sum(mass); L <- c(2,2.3,2.6); q <- .17
evaluate <- function(m,side) {
  e <- r <- numeric(J); p <- matrix(0,J,K)
  for(j in seq_len(J)) {
    zj <- cells$z==j; tj <- zj & cells$t==1
    e[j] <- sum(m[zj]); r[j] <- sum(m[tj])/e[j]
    p[j,] <- m[tj]/sum(m[tj])
  }
  rows <- nrow(cells)
  rr <- matrix(rep(r,each=rows),rows,J); ee <- matrix(rep(e,each=rows),rows,J)
  gam <- t(p[,cells$y,drop=FALSE])/p[1,cells$y]
  cc <- component(support[cells$y],rep(q,rows),rr,gam,L,side)
  Hbar <- Cbar <- numeric(J)
  for(j in seq_len(J)) {
    ix <- cells$z==j & cells$t==1
    Hbar[j] <- sum(p[j,]*cc$H[ix,j]); Cbar[j] <- sum(p[j,]*cc$C[ix,j])
  }
  psi <- q+sum(Hbar)
  signal <- score(rep(q,rows),matrix(rep(Hbar,each=rows),rows),
    matrix(rep(Cbar,each=rows),rows),rr,ee,cells$z,cells$t,1,cc$H)
  list(psi=psi,signal=signal)
}
for(side in c('lower','upper')) {
  a <- evaluate(mass,side)
  stopifnot(abs(sum(mass*a$signal)-a$psi)<1e-12)
  for(b in 1:5) {
    h <- sin(seq_along(mass)*b); eps <- 1e-6
    plus <- mass*exp(eps*h); plus <- plus/sum(plus)
    minus <- mass*exp(-eps*h); minus <- minus/sum(minus)
    numerical <- (evaluate(plus,side)$psi-evaluate(minus,side)$psi)/(2*eps)
    analytic <- sum(mass*(a$signal-a$psi)*h)
    stopifnot(abs(numerical-analytic)<1e-7)
  }
}
