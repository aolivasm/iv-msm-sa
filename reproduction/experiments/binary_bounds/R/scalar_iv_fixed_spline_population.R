# Population fixed-class pseudo-bounds for the new scalar-X factorial DGP.
# Run from the repository root. No fitted nuisance functions or regularization.
source('R/scalar_iv_population_experiment.R')

siv_spline_basis <- function(x) {
  cbind(1, splines::ns(x, knots=c(-1.2,-.4,.4,1.2),
                     Boundary.knots=c(-2,2), intercept=FALSE))
}

siv_spline_arm <- function(a,b,Lambda,t,n_x=80L) {
  rule <- siv_x_rule(n_x); x <- rule$x
  raw <- siv_spline_basis(x)
  # Invertible whitening changes coordinates, not the fixed function class.
  transform <- solve(chol(crossprod(raw,rule$w*raw)))
  design <- raw %*% transform
  env <- lapply(x,function(xx) siv_envelopes(xx,a,b,t,Lambda))
  stopifnot(max(vapply(env,function(e) max(e$slack),numeric(1)))<1e-8)
  loc <- .65*x+.75*t
  p <- vapply(env,function(e) e$p[1],numeric(1))
  lmom <- t(vapply(env,function(e) e$lower(),numeric(2)))
  umom <- t(vapply(env,function(e) e$upper(),numeric(2)))
  total <- umom[,1]-lmom[,1]
  total_y <- umom[,2]-lmom[,2]+loc*total
  LY <- lmom[,2]+loc*lmom[,1]
  ow <- rule$w*(1-p)
  values <- function(q) {
    cut <- pmax(-2,pmin(2,q-loc))
    v <- t(vapply(seq_along(x),function(i)
      env[[i]]$upper(hi=cut[i])-env[[i]]$lower(hi=cut[i]),numeric(2)))
    list(mass=v[,1],moment=v[,2]+loc*v[,1])
  }
  output <- list()
  for(bound in c('lower','upper')) {
    target <- if(bound=='lower') 1-lmom[,1] else umom[,1]-1
    point_q <- loc+vapply(seq_along(x),function(i) {
      uniroot(function(r) env[[i]]$upper(hi=r)[1]-
                env[[i]]$lower(hi=r)[1]-target[i],c(-2,2),tol=1e-11)$root
    },numeric(1))
    objective <- function(beta) {
      q <- as.vector(design%*%beta); v <- values(q)
      sum(ow*(q*(v$mass-target)-v$moment))
    }
    gradient <- function(beta) {
      q <- as.vector(design%*%beta)
      as.vector(crossprod(design,ow*(values(q)$mass-target)))
    }
    initial <- lm.wfit(design,point_q,ow*total)$coefficients
    opt <- optim(initial,objective,gradient,method='BFGS',
                 control=list(maxit=2000,reltol=1e-14))
    # Newton polishing with backtracking: analytic curvature of the integrated
    # check loss. Outside outcome support its curvature is zero.
    beta <- opt$par
    for(iter in seq_len(30)) {
      g <- gradient(beta)
      if(max(abs(g))<1e-11) break
      q <- as.vector(design%*%beta); r <- q-loc
      curvature <- vapply(seq_along(x),function(i) {
        e <- env[[i]]; gamma <- siv_gamma(r[i],e$eta)
        gap <- min(Lambda,e$cu*gamma+e$d)-max(1/Lambda,e$cl*gamma+e$d)
        if(abs(r[i])>=2) 0 else gap*dnorm(r[i])/siv_normalizer(0)
      },numeric(1))
      h <- crossprod(design,(ow*curvature)*design)
      step <- solve(h,g); scale <- 1; old <- objective(beta)
      while(scale>1e-8 && objective(beta-scale*step)>old+1e-14) scale <- scale/2
      beta <- beta-scale*step
    }
    endpoint <- function(q) {
      v <- values(q)
      missing <- if(bound=='lower') LY+v$moment+q*(1-lmom[,1]-v$mass) else
        LY+total_y-v$moment+q*(1-umom[,1]+v$mass)
      sum(rule$w*(p*loc+(1-p)*missing))
    }
    q <- as.vector(design%*%beta); v <- values(q)
    normalization <- if(bound=='lower') lmom[,1]+v$mass-1 else umom[,1]-v$mass-1
    output[[bound]] <- list(endpoint=paste0('mu',t,'_',bound),
      target=endpoint(q),unrestricted=endpoint(point_q),
      maximum_gradient=max(abs(gradient(beta))),
      max_conditional_normalization_error=max(abs(normalization)),
      mean_abs_normalization_error=sum(rule$w*abs(normalization)),
      direct_minus_pseudo=sum(ow*q*normalization),
      coefficients=as.vector(transform%*%beta),optimizer_code=opt$convergence)
  }
  output
}

siv_spline_population <- function(a,b,Lambda,n_x=80L) {
  arms <- c(siv_spline_arm(a,b,Lambda,0,n_x),siv_spline_arm(a,b,Lambda,1,n_x))
  detail <- do.call(rbind,lapply(arms,function(v)
    as.data.frame(v[setdiff(names(v),'coefficients')],stringsAsFactors=FALSE)))
  get <- function(t,b,col='target') detail[detail$endpoint==paste0('mu',t,'_',b),col]
  ate <- c(lower=get(1,'lower')-get(0,'upper'),upper=get(1,'upper')-get(0,'lower'))
  list(detail=detail,ate=ate,models=arms)
}
