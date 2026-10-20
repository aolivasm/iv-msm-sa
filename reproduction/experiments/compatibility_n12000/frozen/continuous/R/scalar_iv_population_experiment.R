# Standalone population calculations; base R only. No changes to production DGPs.
# X ~ N(0,1) truncated to [-2,2], Z independent Bernoulli(1/2).
# logit P(T=1|X,Z) = .5 X + a (Z-.5).
# Y|T=t,X=x,Z=z = .65 x + .75 t + R,
# R ~ N(z*(2t-1)*b,1) truncated to [-2,2].
# This specifies P, not a unique full-data Q or a true population ATE.

siv_normalizer <- function(eta) pnorm(2-eta)-pnorm(-2-eta)
siv_moments <- function(lo=-2, hi=2, eta=0) {
  if (hi <= lo) return(c(mass=0, moment=0))
  z <- siv_normalizer(eta)
  mass <- (pnorm(hi-eta)-pnorm(lo-eta))/z
  c(mass=mass, moment=eta*mass+(dnorm(lo-eta)-dnorm(hi-eta))/z)
}
siv_gamma <- function(r, eta) {
  exp(eta*r-eta^2/2)*siv_normalizer(0)/siv_normalizer(eta)
}
siv_gamma_inverse <- function(value, eta) {
  (log(value)+eta^2/2-log(siv_normalizer(0)/siv_normalizer(eta)))/eta
}
siv_probabilities <- function(x, a, t=1) {
  pi <- plogis(.5*x+c(-a/2,a/2))
  if(t==1) pi else 1-pi
}

# Integrate max/min(k,c*gamma+d) against f0 analytically, including the
# residual first moment. Split at the exact branch crossing, if present.
siv_envelope_integral <- function(c, d, k, eta, maximum, lo=-2, hi=2) {
  if(hi <= lo) return(c(mass=0,moment=0))
  cuts <- c(lo,hi)
  gstar <- (k-d)/c
  if(eta != 0 && gstar > 0) {
    root <- siv_gamma_inverse(gstar,eta)
    if(root > lo && root < hi) cuts <- c(lo,root,hi)
  }
  ans <- c(mass=0,moment=0)
  for(j in seq_len(length(cuts)-1L)) {
    left <- cuts[j]; right <- cuts[j+1L]
    affine <- c*siv_gamma((left+right)/2,eta)+d
    use_affine <- if(maximum) affine > k else affine < k
    m0 <- siv_moments(left,right,0)
    ans <- ans + if(use_affine) c*siv_moments(left,right,eta)+d*m0 else k*m0
  }
  ans
}
siv_envelopes <- function(x,a,b,t,Lambda) {
  p <- siv_probabilities(x,a,t); eta <- (2*t-1)*b
  d <- -p[1]/(1-p[1])
  cl <- ((1-p[2])/Lambda+p[2])/(1-p[1])
  cu <- ((1-p[2])*Lambda+p[2])/(1-p[1])
  lower <- function(lo=-2,hi=2) siv_envelope_integral(cl,d,1/Lambda,eta,TRUE,lo,hi)
  upper <- function(lo=-2,hi=2) siv_envelope_integral(cu,d,Lambda,eta,FALSE,lo,hi)
  g <- siv_gamma(c(-2,2),eta)
  point_gap <- max(pmax(1/Lambda,cl*g+d)-pmin(Lambda,cu*g+d))
  list(p=p,eta=eta,cl=cl,cu=cu,d=d,lower=lower,upper=upper,
       slack=c(pointwise=point_gap,lower_mean=lower()[1]-1,upper_mean=1-upper()[1]))
}

# Global searches are numerical, not formal interval-arithmetic certificates:
# dense support grid INCLUDING endpoints, followed by local maximization.
siv_support_max <- function(fun,n_grid=81L) {
  x <- seq(-2,2,length.out=n_grid)
  y <- vapply(x,fun,numeric(1))
  peaks <- which(y[2:(n_grid-1)] >= y[1:(n_grid-2)] &
                 y[2:(n_grid-1)] >= y[3:n_grid])+1L
  if(diff(range(y)) < 1e-13) peaks <- integer()
  best <- max(y)
  for(i in peaks) best <- max(best,optimize(fun,c(x[i-1],x[i+1]),
                                           maximum=TRUE,tol=1e-9)$objective)
  best
}
siv_compatibility <- function(a,b,Lambda,n_grid=81L) {
  unlist(lapply(0:1,function(t) sapply(seq_len(3),function(j)
    siv_support_max(function(x)
      unname(siv_envelopes(x,a,b,t,Lambda)$slack[j]),n_grid))))
}
siv_iv_only_violation <- function(a,b,n_grid=81L) {
  # Marginal IV independence is feasible iff integral max(p0*f0,p1*f1)<=1
  # in both arms. This does not impose joint potential-outcome independence.
  max(vapply(0:1,function(t) siv_support_max(function(x) {
    p <- siv_probabilities(x,a,t)
    siv_envelope_integral(p[2],0,p[1],(2*t-1)*b,TRUE)[1]-1
  },n_grid),numeric(1)))
}
siv_frontier <- function(a,b,n_grid=81L) {
  if(siv_iv_only_violation(a,b,n_grid) > 1e-9) return(Inf)
  if(b==0) return(1)
  f <- function(L) max(siv_compatibility(a,b,L,n_grid))
  upper <- 2
  while(f(upper)>0 && upper<1024) upper <- 2*upper
  if(f(upper)>0) return(Inf)
  uniroot(f,c(1,upper),tol=1e-8)$root
}
siv_x_rule <- function(n=80L) {
  j <- seq_len(n-1L); v <- j/sqrt(4*j*j-1)
  mat <- matrix(0,n,n); mat[cbind(j,j+1)] <- v; mat[cbind(j+1,j)] <- v
  e <- eigen(mat,symmetric=TRUE); idx <- order(e$values)
  x <- 2*e$values[idx]
  w <- 4*e$vectors[1,idx]^2*dnorm(x)/siv_normalizer(0)
  list(x=x,w=w/sum(w))
}

siv_extreme <- function(env) {
  L <- env$lower(); U <- env$upper(); remaining <- 1-L[1]
  if(max(env$slack)>2e-8) return(c(lower=NA,upper=NA,norm_error=NA))
  if(U[1]-L[1] < 1e-12) return(c(lower=L[2],upper=L[2],norm_error=abs(L[1]-1)))
  mass_to <- function(q) env$upper(hi=q)[1]-env$lower(hi=q)[1]
  quant <- function(target) {
    if(target<=1e-12) return(-2)
    if(target>=U[1]-L[1]-1e-12) return(2)
    uniroot(function(q) mass_to(q)-target,c(-2,2),tol=1e-11)$root
  }
  qlo <- quant(remaining); qhi <- quant(U[1]-L[1]-remaining)
  # Evaluate U below qlo / L above qlo; L below qhi / U above qhi.
  lo <- env$upper(hi=qlo)+env$lower(lo=qlo)
  hi <- env$lower(hi=qhi)+env$upper(lo=qhi)
  c(lower=unname(lo[2]),upper=unname(hi[2]),norm_error=max(abs(c(lo[1],hi[1])-1)))
}
siv_msm_residual <- function(eta,Lambda) {
  m <- siv_moments(eta=eta)[2]
  if(Lambda==1) return(c(lower=m,upper=m))
  tau <- 1/(Lambda+1)
  quant <- function(p) eta+qnorm(pnorm(-2-eta)+p*siv_normalizer(eta))
  c(lower=m/Lambda+(Lambda-1/Lambda)*siv_moments(hi=quant(tau),eta=eta)[2],
    upper=m/Lambda+(Lambda-1/Lambda)*siv_moments(lo=quant(1-tau),eta=eta)[2])
}

# Conventional MSM adjusting only for X. The observed arm distribution is
# a mixture over Z with weights P(Z=z|T=t,X), not P(Z=z|X).
siv_msm_x_population <- function(a,b,Lambda,n_x=80L) {
  rule <- siv_x_rule(n_x)
  mu <- matrix(0,2,2,dimnames=list(c('mu0','mu1'),c('lower','upper')))
  for(i in seq_along(rule$x)) for(t in 0:1) {
    x <- rule$x[i]; p <- siv_probabilities(x,a,t)
    pt <- mean(p); w0 <- p[1]/sum(p); eta <- (2*t-1)*b
    mix <- function(lo=-2,hi=2) w0*siv_moments(lo,hi,0)+(1-w0)*siv_moments(lo,hi,eta)
    m <- mix()[2]
    if(Lambda==1) ex <- rep(m,2) else {
      tau <- 1/(Lambda+1)
      quant <- function(prob) uniroot(function(q) mix(hi=q)[1]-prob,
                                     c(-2,2),tol=1e-11)$root
      ex <- m/Lambda+(Lambda-1/Lambda)*c(mix(hi=quant(tau))[2],mix(lo=quant(1-tau))[2])
    }
    mu[t+1,] <- mu[t+1,]+rule$w[i]*(.65*x+.75*t+pt*m+(1-pt)*ex)
  }
  list(mu=mu,ate=c(lower=unname(mu[2,1]-mu[1,2]),upper=unname(mu[2,2]-mu[1,1])))
}
siv_margin_mass <- function(c,d,k,eta,s) {
  if(eta==0) return(as.numeric(abs(c+d-k)<=s))
  low <- (k-d-s)/c; high <- (k-d+s)/c
  if(high<=0) return(0)
  rr <- c(if(low<=0) if(eta>0) -Inf else Inf else siv_gamma_inverse(low,eta),
          siv_gamma_inverse(high,eta))
  siv_moments(max(-2,min(rr)),min(2,max(rr)))[1]
}
siv_population <- function(a,b,Lambda,compatible,n_x=80L) {
  rule <- siv_x_rule(n_x); thresholds <- c(1e-4,1e-3,.01,.05,.1)
  mu <- matrix(0,2,4,dimnames=list(c('mu0','mu1'),c('iv_lower','iv_upper','msm_lower','msm_upper')))
  margin <- array(0,c(2,2,length(thresholds)))
  norm_error <- 0; t0_mass <- numeric(2); observed_sum <- numeric(2)
  for(i in seq_along(rule$x)) for(t in 0:1) {
    x <- rule$x[i]; w <- rule$w[i]; location <- .65*x+.75*t
    env <- siv_envelopes(x,a,b,t,Lambda); p <- env$p; eta <- env$eta
    m0 <- siv_moments()[2]; m1 <- siv_moments(eta=eta)[2]
    observed_sum[t+1] <- observed_sum[t+1]+w*.5*(p[1]*(location+m0)+p[2]*(location+m1))
    msm <- .5*(p[1]*m0+(1-p[1])*siv_msm_residual(0,Lambda)+
               p[2]*m1+(1-p[2])*siv_msm_residual(eta,Lambda))+location
    mu[t+1,3:4] <- mu[t+1,3:4]+w*msm
    if(compatible) {
      ex <- siv_extreme(env)
      mu[t+1,1:2] <- mu[t+1,1:2]+w*(location+p[1]*m0+(1-p[1])*ex[1:2])
      norm_error <- max(norm_error,ex[3])
    }
    t0_mass[t+1] <- t0_mass[t+1]+w*p[1]
    for(j in seq_along(thresholds)) {
      margin[t+1,1,j] <- margin[t+1,1,j]+w*p[1]*siv_margin_mass(env$cl,env$d,1/Lambda,eta,thresholds[j])
      margin[t+1,2,j] <- margin[t+1,2,j]+w*p[1]*siv_margin_mass(env$cu,env$d,Lambda,eta,thresholds[j])
    }
  }
  if(!compatible) mu[,1:2] <- NA_real_
  ate <- c(iv_lower=mu[2,1]-mu[1,2],iv_upper=mu[2,2]-mu[1,1],
           msm_lower=mu[2,3]-mu[1,4],msm_upper=mu[2,4]-mu[1,3])
  names(ate) <- c('iv_lower','iv_upper','msm_lower','msm_upper')
  margins <- expand.grid(treatment=0:1,branch=c('lower','upper'),s=thresholds)
  margins$probability <- as.vector(sweep(margin,1,t0_mass,'/'))
  list(ate=ate,mu=mu,margin=margins,norm_error=if(compatible) norm_error else NA_real_,
       atu=2*mu[2,1:2]-2*sum(observed_sum))
}
siv_design_metrics <- function(a,b,n_x=160L) {
  r <- siv_x_rule(n_x); p0 <- plogis(.5*r$x-a/2);p1 <- plogis(.5*r$x+a/2)
  pt <- (p0+p1)/2
  cross <- if(b==0) 0 else siv_gamma_inverse(1,b)
  tv <- if(b==0) 0 else siv_moments(lo=cross,eta=b)[1]-siv_moments(lo=cross,eta=0)[1]
  chi2 <- exp(b*b)*siv_normalizer(2*b)*siv_normalizer(0)/siv_normalizer(b)^2-1
  er <- siv_moments(eta=b)[2]
  data.frame(a=a,b=b,first_stage=sum(r$w*(p1-p0)),
    partial_R2=sum(r$w*.25*(p1-p0)^2)/sum(r$w*pt*(1-pt)),
    treatment_prevalence=sum(r$w*pt),min_conditional_TZ_cell=.5*plogis(-1-a/2),
    min_marginal_TZ_cell=.5*sum(r$w*p0),TV=unname(tv),chi_squared=chi2,
    max_abs_log_gamma=max(abs(log(siv_gamma(c(-2,2),b)))),
    reduced_form=sum(r$w*(.75*(p1-p0)+(2*p1-1)*er)))
}
