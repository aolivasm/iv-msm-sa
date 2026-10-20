# Display methods preserve all original numerical results and diagnostics.
.ivmsm_result <- function(x,kind) {
  x$analysis_type <- kind
  for(nm in c('ate_confidence_set','iv_ate_confidence_set'))
    if(!is.null(x[[nm]])) rownames(x[[nm]]) <- 'ATE'
  class(x) <- c('ivmsm_fit','list'); x
}
summary.ivmsm_fit <- function(object,...) {
  ep <- if(!is.null(object$endpoints)) object$endpoints else
    if(!is.null(object$endpoint_summary)) object$endpoint_summary else object$endpoint_intervals
  ep <- ep[,intersect(c('endpoint','estimate','standard_error'),names(ep)),drop=FALSE]
  sets <- object$outward_sets
  if(is.null(sets)) sets <- cbind(parameter='ATE',object$ate_confidence_set)
  sets <- sets[,intersect(c('parameter','confidence_level','lower','upper'),names(sets)),drop=FALSE]
  d <- object$diagnostics; rows <- list()
  labels <- c(compatibility_projection_rate='Adjusted fitted values (%)',
    post_smooth_theta_crossing_rate='Remaining bound inconsistencies (%)',
    projection_rate='Adjusted fitted values (%)',
    incompatibility_rate='Initial bound inconsistencies (%)',
    heldout_gamma_projection='Adjusted fitted ratios (%)',
    conditional_mean_projection='Adjusted conditional means (%)',
    smooth_any_failure_rate='Bound inconsistencies (%)')
  for(nm in intersect(names(labels),names(d))) {
    for(i in seq_len(nrow(d))) {
      group <- paste(c(if('fold'%in%names(d))paste('Fold',d$fold[i]),
        if('treatment'%in%names(d))paste('treatment',d$treatment[i]),
        if('endpoint'%in%names(d))d$endpoint[i]),collapse=', ')
      rows[[length(rows)+1L]] <- data.frame(group=group,check=unname(labels[nm]),percent=100*d[[nm]][i])
    }
  }
  diag <- if(length(rows)) do.call(rbind,rows) else data.frame()
  # Keep display short without averaging unequal-sized validation folds.
  maximum <- nrow(diag)>8L
  if(maximum) diag <- aggregate(percent~check,data=diag,FUN=max)
  structure(list(type=object$analysis_type,Lambda=object$Lambda,
    Lambda0=object$Lambda0,Lambda1=object$Lambda1,
    epsilon=object$epsilon,envelope_epsilon=object$envelope_epsilon,
    threshold_epsilon=object$threshold_epsilon,endpoints=ep,outward_sets=sets,
    diagnostics=diag,diagnostic_maximum=maximum,
    warnings=object$warnings),class='summary.ivmsm_fit')
}
print.summary.ivmsm_fit <- function(x,digits=4,diagnostics=FALSE,...) {
  cat('IV-MSM: ',x$type,'\n',sep='')
  if(!is.null(x$Lambda0)) cat('Sensitivity by IV stratum: Lambda_0 =',x$Lambda0,
    '; Lambda_1 =',x$Lambda1,'\n') else if(!is.null(x$Lambda)) {
      cat('Sensitivity by IV category:\n'); print(x$Lambda,digits=digits)
    }
  if(!is.null(x$epsilon)) cat('Envelope smoothing epsilon =',x$epsilon,'(0 = unsmoothed)\n')
  if(!is.null(x$envelope_epsilon)) cat('Envelope smoothing =',x$envelope_epsilon,
    '; threshold smoothing =',x$threshold_epsilon,'\n')
  ep <- x$endpoints
  pretty <- c(mu0_lower='Mean under no treatment: lower',mu0_upper='Mean under no treatment: upper',
    mu1_lower='Mean under treatment: lower',mu1_upper='Mean under treatment: upper',
    ate_lower='ATE: lower',ate_upper='ATE: upper')
  if(nrow(ep)) ep$endpoint <- ifelse(ep$endpoint%in%names(pretty),pretty[ep$endpoint],ep$endpoint)
  names(ep) <- c('Endpoint','Estimate','SE')[seq_len(ncol(ep))]
  cat('\nOne-step endpoint estimates\n');print(ep,row.names=FALSE,digits=digits)
  cs <- x$outward_sets
  names(cs) <- c(parameter='Parameter',confidence_level='Level',lower='Lower',upper='Upper')[names(cs)]
  cat('\nOutward confidence sets (pointwise)\n');print(cs,row.names=FALSE,digits=digits)
  if(any(cs$Lower>cs$Upper,na.rm=TRUE))
    cat('Warning: some lower confidence limits exceed upper limits; inspect the fitted results.\n')
  if(isTRUE(diagnostics) && nrow(x$diagnostics)) {
    cat('\nNumerical checks',if(x$diagnostic_maximum)' (maximum across displayed groups)' else '', '\n',sep='')
    print(x$diagnostics,row.names=FALSE,digits=digits)
  }
  if(isTRUE(diagnostics)) {
    cat('\nAdjustments are numerical safeguards, not compatibility tests or a guarantee of coverage.\n')
    cat('Full details: fit$diagnostics. Inference requires the stated regularity conditions.\n')
  }
  invisible(x)
}
print.ivmsm_fit <- function(x,...) {print(summary(x),...);invisible(x)}

mock_iv_data <- function(outcome=c('binary','ordered','continuous'),iv_categories=2L) {
  outcome <- match.arg(outcome)
  if(length(iv_categories)!=1L||!is.numeric(iv_categories)||!is.finite(iv_categories)||
     iv_categories<2||iv_categories!=floor(iv_categories)) stop('iv_categories must be an integer >= 2.')
  d <- expand.grid(person=1:60,T_obs=0:1,Z=seq_len(iv_categories)-1L)
  d$X <- (d$person-30.5)/30
  d$Y <- switch(outcome,binary=as.integer((d$person+d$T_obs)%%4<2),
    ordered=as.integer((d$person+d$T_obs)%%3),
    continuous=.5*d$T_obs+.2*d$X+.4*sin(d$person*1.7))
  d[,c('Y','T_obs','Z','X')]
}

estimate_binary_sensitivity_band <- function(data,x_cols='X',Lambda_range=c(1.2,2),
    grid_points=21L,confidence_level=.95,draws=10000L,seed=1L,...) {
  if(length(Lambda_range)!=2L||any(!is.finite(Lambda_range))||
     Lambda_range[1]<1||Lambda_range[2]<=Lambda_range[1]) stop('Supply an increasing sensitivity interval starting at >=1.')
  if(length(grid_points)!=1L||!is.finite(grid_points)||grid_points<2||grid_points!=floor(grid_points)) stop('Invalid grid_points.')
  dots <- list(...)
  if(any(c('Lambda0','Lambda1','epsilon','data','x_cols','confidence_level','seed')%in%names(dots)))
    stop('This band follows Lambda_0 = Lambda_1 with unsmoothed envelopes; do not override these arguments.')
  fit <- do.call(estimate_binary_dual_bounds,c(list(data=data,x_cols=x_cols,
    Lambda0=Lambda_range[1],Lambda1=Lambda_range[1],epsilon=0,
    confidence_level=confidence_level,seed=seed),dots))
  grid <- seq(Lambda_range[1],Lambda_range[2],length.out=grid_points)
  lo <- hi <- matrix(NA_real_,nrow(data),grid_points); diagnostic <- vector('list',grid_points)
  for(j in seq_along(grid)) {
    f <- evaluate_binary_bounds_from_nuisances(data,fit$nuisance_fit,
      Lambda0=grid[j],Lambda1=grid[j],epsilon=0,confidence_level=confidence_level,
      compatibility_action=fit$compatibility_action,compatibility_tolerance=fit$compatibility_tolerance)
    lo[,j] <- f$signals[['1']]$lower-f$signals[['0']]$upper
    hi[,j] <- f$signals[['1']]$upper-f$signals[['0']]$lower
    diagnostic[[j]] <- f$diagnostics
  }
  band <- binary_grid_band(grid,lo,hi,level=confidence_level,draws=draws,seed=seed)
  band$continuum <- binary_continuum_band(grid,band$table$simultaneous_lower,
    band$table$simultaneous_upper,at=seq(grid[1],tail(grid,1),length.out=201))
  band$diagnostics <- diagnostic;band$Lambda_range <- Lambda_range
  band$confidence_level <- confidence_level
  class(band) <- c('ivmsm_band','list'); band
}
summary.ivmsm_band <- function(object,...) object
print.ivmsm_band <- function(x,digits=4,...) {
  cat('IV-MSM: binary outcome, binary IV\nCommon sensitivity path: Lambda_0 = Lambda_1 in [',
    x$Lambda_range[1],', ',x$Lambda_range[2],']\n',sep='')
  cat(100*x$confidence_level,'% simultaneous outward ATE band; unsmoothed envelopes\n',sep='')
  cat('Calibrated jointly at',nrow(x$table),'sensitivity values; extended over the interval.\n')
  cat('Band limits evaluated at',nrow(x$continuum),'locations: fit$continuum. Use plot(fit) to view the band.\n')
  if(any(x$table$empty)||any(x$continuum$empty))
    cat('Warning: some lower confidence limits exceed upper limits; inspect the fitted results.\n')
  invisible(x)
}

plot.ivmsm_band <- function(x,xlab=expression(Lambda[0]==Lambda[1]),
    ylab='ATE',main='Simultaneous outward confidence band',...) {
  d <- x$continuum
  graphics::plot(d$Lambda,d$lower,type='n',ylim=range(c(d$lower,d$upper)),
    xlab=xlab,ylab=ylab,main=main,...)
  # Shade only adjacent valid intervals; never repair or fill crossed limits.
  for(i in seq_len(nrow(d)-1L)) if(!any(d$empty[c(i,i+1L)])) {
    j <- c(i,i+1L)
    graphics::polygon(c(d$Lambda[j],rev(d$Lambda[j])),
      c(d$lower[j],rev(d$upper[j])),col='#DCE9F4',border=NA)
  }
  graphics::lines(d$Lambda,d$lower,col='#204A78',lwd=2)
  graphics::lines(d$Lambda,d$upper,col='#204A78',lwd=2)
  graphics::abline(h=0,lty=3,col='gray50')
  if(any(d$empty)) warning('Crossed confidence limits are shown without repair or shading.',call.=FALSE)
  invisible(x)
}
