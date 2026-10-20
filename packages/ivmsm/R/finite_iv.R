# General finite-IV interfaces; no experiment-specific DGP or file paths.
.iv_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)))
    stop(name, ' must be a finite numeric matrix.')
  if (nrow(x) < 1L || ncol(x) < 2L) stop(name, ' must have at least two IV columns.')
  x
}
.iv_lambda <- function(Lambda, J) {
  if (!is.numeric(Lambda) || !length(Lambda) %in% c(1L, J) ||
      any(!is.finite(Lambda)) || any(Lambda < 1))
    stop('Lambda must be finite, at least one, and scalar or length J.')
  rep(Lambda, length.out = J)
}

finite_iv_bounds <- function(pi, p0, p1, support, Lambda = 2, weights = NULL) {
  pi <- .iv_matrix(pi, 'pi'); n <- nrow(pi); J <- ncol(pi)
  if (any(pi <= 0 | pi >= 1)) stop('Treatment probabilities must be strictly interior.')
  if (!is.numeric(support) || length(support) < 2L ||
      any(!is.finite(support)) || any(diff(support) <= 0))
    stop('support must be a strictly increasing finite numeric vector.')
  K <- length(support); L <- .iv_lambda(Lambda, J)
  for (p in list(p0, p1)) {
    if (!is.numeric(p) || !identical(dim(p), c(n, J, K)) ||
        any(!is.finite(p)) || any(p < 0 | p > 1) ||
        max(abs(apply(p, c(1, 2), sum) - 1)) > 1e-8)
      stop('p0 and p1 must be n by J by length(support) probability arrays summing to one.')
  }
  if (is.null(weights)) weights <- rep(1/n, n)
  if (!is.numeric(weights) || length(weights) != n || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) stop('Invalid averaging weights.')
  weights <- weights / sum(weights)
  arm <- function(r, p) {
    A <- r + sweep(1-r, 2, L, '/')
    B <- r + sweep(1-r, 2, L, '*')
    lower <- upper <- matrix(NA_real_, n, K)
    for (k in seq_len(K)) {
      pk <- matrix(p[,,k], n, J)
      lower[,k] <- apply(A * pk, 1, max)
      upper[,k] <- apply(B * pk, 1, min)
    }
    endpoints <- matrix(NA_real_, n, 2L, dimnames = list(NULL, c('lower','upper')))
    feasible <- logical(n)
    for (i in seq_len(n)) {
      a <- discrete_waterfill(lower[i,], upper[i,], support, 'lower')
      b <- discrete_waterfill(lower[i,], upper[i,], support, 'upper')
      feasible[i] <- a$feasible && b$feasible
      if (feasible[i]) endpoints[i,] <- c(a$value, b$value)
    }
    list(endpoints = endpoints, feasible = feasible, lower_mass = lower,
         upper_mass = upper, diagnostics = data.frame(
           pointwise_crossing = apply(pmax(lower-upper, 0), 1, max),
           lower_mass_excess = pmax(rowSums(lower)-1, 0),
           upper_mass_deficit = pmax(1-rowSums(upper), 0)))
  }
  a0 <- arm(1-pi, p0); a1 <- arm(pi, p1)
  conditional <- cbind(mu0_lower=a0$endpoints[,1], mu0_upper=a0$endpoints[,2],
    mu1_lower=a1$endpoints[,1], mu1_upper=a1$endpoints[,2],
    ate_lower=a1$endpoints[,1]-a0$endpoints[,2],
    ate_upper=a1$endpoints[,2]-a0$endpoints[,1])
  compatible <- all(a0$feasible & a1$feasible)
  # Never average only the feasible subset or silently change the population.
  estimates <- if (compatible) colSums(conditional * weights) else
    setNames(rep(NA_real_, ncol(conditional)), colnames(conditional))
  list(estimates=estimates, conditional=conditional, compatible=compatible,
       arms=list(mu0=a0, mu1=a1), Lambda=L, support=support, weights=weights,
       inference='Population/plugin bounds only; no standard errors or formal test.',
       projection_rate=0)
}

fit_finite_iv_bounds <- function(data, outcome, treatment, instrument,
                                 support, Lambda = 2, strata = character()) {
  cols <- c(outcome, treatment, instrument, strata)
  if (!is.data.frame(data) || !all(cols %in% names(data)) || !nrow(data) ||
      anyNA(data[,cols,drop=FALSE])) stop('Supply complete data and valid column names.')
  y <- data[[outcome]]; t <- data[[treatment]]
  if (!is.numeric(y) || !all(y %in% support) || !all(t %in% c(0,1)))
    stop('Outcomes must lie in support and treatment must be coded 0/1.')
  z <- factor(data[[instrument]]); J <- nlevels(z)
  if (J < 2L) stop('At least two observed instrument categories are required.')
  g <- if (!length(strata)) factor(rep('all', nrow(data))) else
    do.call(interaction, c(lapply(data[strata], factor), list(drop=TRUE, lex.order=TRUE)))
  n <- nlevels(g); K <- length(support)
  pi <- matrix(NA_real_, n, J); p0 <- p1 <- array(0, c(n,J,K))
  counts <- array(0L, c(n,J,2L))
  for (i in seq_len(n)) for (j in seq_len(J)) {
    cell <- as.integer(g)==i & as.integer(z)==j
    counts[i,j,] <- c(sum(cell & t==0), sum(cell & t==1))
    if (any(counts[i,j,] == 0)) stop('Every stratum/IV cell needs both treatments; no smoothing is applied.')
    pi[i,j] <- mean(t[cell])
    p0[i,j,] <- tabulate(match(y[cell & t==0],support), nbins=K)/counts[i,j,1]
    p1[i,j,] <- tabulate(match(y[cell & t==1],support), nbins=K)/counts[i,j,2]
  }
  ans <- finite_iv_bounds(pi,p0,p1,support,Lambda,as.numeric(table(g)))
  ans$instrument_levels <- levels(z); ans$stratum_levels <- levels(g)
  ans$cell_counts <- counts
  ans$nuisances <- list(pi=pi,p0=p0,p1=p1)
  ans$inference <- 'Empirical cell plugin only; no cross-fitting, regularization, or confidence intervals.'
  ans
}

finite_iv_envelopes <- function(r, gamma, Lambda = 2, reference = 1L) {
  r <- .iv_matrix(r,'r'); gamma <- .iv_matrix(gamma,'gamma')
  if (!identical(dim(r),dim(gamma)) || any(r<=0 | r>=1) || any(gamma<=0) ||
      length(reference)!=1L || !reference %in% seq_len(ncol(r)) ||
      max(abs(gamma[,reference]-1))>1e-8) stop('Invalid arm propensities, density ratios, or reference.')
  L <- .iv_lambda(Lambda,ncol(r))
  lower <- (r+sweep(1-r,2,L,'/'))*gamma
  upper <- (r+sweep(1-r,2,L,'*'))*gamma
  a <- apply(lower,1,max); b <- apply(upper,1,min)
  list(lower=a,upper=b,crossed=a>b,lower_candidates=lower,upper_candidates=upper,
       reference=reference,projection_rate=0,
       scope='Common potential-outcome density relative to reference observed arm density; not latent-selection envelopes. Normalization feasibility must also be checked.')
}
