# Covariate-assisted exact Balke--Pearl bounds, Levis et al. (2025),
# equations (1)--(3), with cross-fitted selected-branch one-step signals.
# Input predictions must be out-of-fold; this file does not fit nuisances.
# Cell order: 00.0, 01.0, 10.0, 11.0, 00.1, 01.1, 10.1, 11.1 (Y,T.Z).
bp_coefficients <- function() {
  lower <- rbind(c(1,0,0,0,0,0,0,1), c(0,0,0,1,1,0,0,0),
    c(0,0,0,0,0,-1,-1,0), c(0,-1,-1,0,0,0,0,0),
    c(0,-1,-1,1,0,0,-1,-1), c(0,0,-1,-1,0,-1,-1,1),
    c(-1,-1,0,0,1,-1,-1,0), c(1,-1,-1,0,-1,-1,0,0))
  upper <- rbind(c(0,0,-1,0,0,-1,0,0), c(0,-1,0,0,0,0,-1,0),
    c(0,0,0,0,1,0,0,1), c(1,0,0,1,0,0,0,0),
    c(1,-1,0,1,1,1,0,0), c(1,1,0,0,1,-1,0,1),
    c(0,0,1,1,1,0,-1,1), c(1,0,-1,1,0,0,1,1))
  list(lower=lower, upper=upper, lower_offset=c(-1,-1,rep(0,6)),
       upper_offset=c(1,1,rep(0,6)))
}

bp_cell_predictions <- function(p) {
  out <- do.call(cbind, lapply(0:1, function(z) {
    r <- p[[paste0('pi',z)]]
    a0 <- p[[paste0('a_t0_z',z)]]; a1 <- p[[paste0('a_t1_z',z)]]
    cbind((1-r)*(1-a0), r*(1-a1), (1-r)*a0, r*a1)
  }))
  colnames(out) <- c('00.0','01.0','10.0','11.0','00.1','01.1','10.1','11.1')
  stopifnot(all(is.finite(out)), all(out >= 0), all(out <= 1),
    max(abs(rowSums(out[,1:4])-1)) < 1e-12,
    max(abs(rowSums(out[,5:8])-1)) < 1e-12)
  out
}

bp_evaluate <- function(data, predictions) {
  n <- nrow(data); p <- bp_cell_predictions(predictions)
  stopifnot(nrow(p)==n, all(data$Y %in% 0:1), all(data$T_obs %in% 0:1),
    all(data$Z %in% 0:1), all(predictions$e>0 & predictions$e<1))
  residual <- p * 0
  for (z in 0:1) for (y in 0:1) for (t in 0:1) {
    j <- 4*z + 2*y + t + 1
    ez <- if(z==0) predictions$e else 1-predictions$e
    residual[,j] <- (data$Z==z)/ez * ((data$Y==y & data$T_obs==t)-p[,j])
  }
  cf <- bp_coefficients(); values <- signals <- winners <- gaps <- list()
  for (side in c('lower','upper')) {
    branches <- sweep(p %*% t(cf[[side]]),2,cf[[paste0(side,'_offset')]],'+')
    # Deterministic first branch at ties; ties are recorded, not smoothed away.
    w <- max.col(if(side=='lower') branches else -branches, ties.method='first')
    values[[side]] <- branches[cbind(seq_len(n),w)]
    signals[[side]] <- values[[side]] + rowSums(residual * cf[[side]][w,,drop=FALSE])
    winners[[side]] <- w
    gaps[[side]] <- apply(branches,1,function(v) {
      s <- sort(v, decreasing=side=='lower'); abs(s[1]-s[2])
    })
  }
  signal <- do.call(cbind,signals)
  estimate <- colMeans(signal); se <- apply(signal,2,sd)/sqrt(n)
  # Bonferroni outward interval, matching existing IV-only reporting.
  ci <- estimate + c(-1,1)*qnorm(.975)*se
  list(estimate=estimate, se=se, outward=ci, plugin=sapply(values,mean),
    signals=signal, covariance=cov(signal)/n, winners=winners, gaps=gaps,
    cells=p, fitted_crossing_rate=mean(values$lower>values$upper),
    diagnostics=do.call(rbind,lapply(names(gaps),function(s) data.frame(
      endpoint=s, tolerance=c(1e-12,.001,.005,.01),
      fraction=sapply(c(1e-12,.001,.005,.01),function(a) mean(gaps[[s]]<=a))))))
}
