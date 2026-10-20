# General data interface for the binary nested-IV selected-branch estimator.
# No experiment generator, oracle, target file, or threshold optimization.
.biv_bounds <- function(r,p,L) {
  stopifnot(identical(dim(r),dim(p)),all(is.finite(c(r,p))),
    all(r>0 & r<1),all(p>=0 & p<=1))
  J <- ncol(r); n <- nrow(r); L <- .iv_lambda(L,J)
  A <- r+sweep(1-r,2,L,'/'); B <- r+sweep(1-r,2,L,'*')
  pick <- function(v,lower) {
    k <- max.col(if(lower) v else -v,ties.method='first')
    j <- (k-1L)%%J+1L; ix <- cbind(seq_len(n),j); second <- k>J
    if(lower) {
      dr <- ifelse(second,-(1-L[j])*(1-p[ix]),(1-1/L[j])*p[ix])
      dp <- ifelse(second,B[ix],A[ix])
    } else {
      dr <- ifelse(second,-(1-1/L[j])*(1-p[ix]),(1-L[j])*p[ix])
      dp <- ifelse(second,A[ix],B[ix])
    }
    value <- v[cbind(seq_len(n),k)]
    list(value=value,j=j,dr=dr,dp=dp,tie=rowSums(abs(v-value)<1e-10)>1L)
  }
  list(lower=pick(cbind(A*p,1-B*(1-p)),TRUE),
       upper=pick(cbind(B*p,1-A*(1-p)),FALSE))
}
.biv_project <- function(r,p,L,action) {
  v <- .biv_bounds(r,p,L); gap <- pmax(v$lower$value-v$upper$value,0)
  bad <- gap>1e-10; alpha <- rep(1,nrow(p))
  if(any(bad) && action=='stop') stop('Estimated binary bounds cross; compatibility_action=stop.')
  if(any(bad) && action=='project') {
    rr <- r[bad,,drop=FALSE]; pp <- p[bad,,drop=FALSE]; base <- pp[,1]
    lo <- rep(0,nrow(pp)); hi <- rep(1,nrow(pp))
    for(it in 1:45) {
      m <- (lo+hi)/2; b <- .biv_bounds(rr,base+(pp-base)*m,L)
      ok <- b$lower$value<=b$upper$value; lo[ok] <- m[ok]; hi[!ok] <- m[!ok]
    }
    alpha[bad] <- lo*(1-1e-8)
  }
  list(p=p[,1]+(p-p[,1])*alpha,rate=mean(bad),max_gap=max(gap),
       projection_rate=mean(alpha<1),max_shrinkage=max(1-alpha))
}
.biv_signal <- function(d,t,e,r,p,L) {
  b <- .biv_bounds(r,p,L)
  one <- function(v) {
    ix <- cbind(seq_len(nrow(d)),v$j); I <- as.numeric(d$Z==v$j)
    D <- as.numeric(d$T_obs==t)
    v$value+I/e[ix]*(v$dr*(D-r[ix])+v$dp*D/r[ix]*(d$Y-p[ix]))
  }
  cbind(lower=one(b$lower),upper=one(b$upper))
}
.biv_base <- function(tr,va,t,x,ctl,seed) {
  use <- tr$Z==1L & tr$T_obs==t; d <- tr[use,,drop=FALSE]
  if(!length(x)) return(list(p=rep(mean(d$Y),nrow(va)),tuning=list(intercept=TRUE)))
  f <- binary_fit_regression_with_learner(d,'Y',x,spline_df=ctl$spline_df,
    learner=ctl$primitive_learner,library_cv_folds=ctl$primitive_cv_folds,
    super_learner_candidates=ctl$super_learner_candidates,
    super_learner_cv_folds=ctl$primitive_cv_folds,seed=seed)
  list(p=f$predict(va),tuning=.fiv_compact(f))
}
.biv_gamma <- function(pd,pn,t,x,ctl,seed) {
  args <- list(train=pd,nuisance_train=pn,t_value=t,x_cols=x,
    backend=ctl$gamma_backend,nystrom_rank=ctl$gamma_rank,
    bandwidth_sample=ctl$gamma_bandwidth_sample,seed=seed)
  f <- if(is.null(ctl$gamma_fixed_lambda)) do.call(fit_gamma_rkhs,c(args,list(
    lambda_grid=function(n_training,treatment=NULL) ctl$gamma_C/sqrt(n_training),
    bandwidth_multipliers=1/sqrt(2),cv_folds=ctl$gamma_cv_folds,
    cv_repeats=ctl$gamma_cv_repeats,selection_rule='minimum'))) else
    do.call(fit_gamma_rkhs_fixed,c(args,list(lambda=ctl$gamma_fixed_lambda,
      bandwidth_multiplier=1/sqrt(2))))
  f$selected_C <- f$selected_lambda*sqrt(sum(pd$T_obs==t)); f
}

estimate_finite_iv_binary <- function(data,outcome,treatment,instrument,
    covariates=character(),Lambda=2,reference=NULL,folds=5L,fold_id=NULL,
    control=list(),confidence_level=.95,seed=1L,keep_signals=TRUE) {
  if(!all(vapply(list(outcome,treatment,instrument),function(v)
      is.character(v)&&length(v)==1L&&!is.na(v)&&nzchar(v),logical(1))) ||
     !is.character(covariates)) stop('Supply column names.')
  cols <- c(outcome,treatment,instrument,covariates)
  if(!is.data.frame(data)||!nrow(data)||anyDuplicated(cols)||
     !all(cols%in%names(data))||anyNA(data[,cols,drop=FALSE]))
    stop('Supply complete data and distinct column names.')
  if(!is.numeric(data[[outcome]])||!all(data[[outcome]]%in%0:1)||
     !is.numeric(data[[treatment]])||!all(data[[treatment]]%in%0:1))
    stop('Outcome and treatment must be numeric 0/1.')
  if(!all(vapply(data[covariates],function(v)is.numeric(v)&&all(is.finite(v)),logical(1))))
    stop('Covariates must be finite numeric columns; encode categorical covariates explicitly.')
  labels <- sort(unique(as.character(data[[instrument]]))); J <- length(labels); n <- nrow(data)
  if(J<2L) stop('At least two observed instrument categories required.')
  if(!is.null(names(Lambda))) {
    if(anyDuplicated(names(Lambda))||length(Lambda)!=J||!setequal(names(Lambda),labels))
      stop('Named Lambda must match the observed instrument labels.')
    Lambda <- Lambda[labels]
  }
  L <- .iv_lambda(Lambda,J)
  if(is.null(reference)) reference <- labels[1L]
  if(length(reference)!=1L||is.na(reference)||!as.character(reference)%in%labels)
    stop('reference must be an observed instrument label.')
  order <- c(match(as.character(reference),labels),which(labels!=as.character(reference)))
  labels <- labels[order]; L <- L[order]
  for(v in list(folds,seed)) if(!is.numeric(v)||length(v)!=1L||!is.finite(v)||v!=floor(v)) stop('Invalid folds or seed.')
  if(folds<2||folds>n||seed<0||seed>.Machine$integer.max-10000*folds-1000*J) stop('Invalid folds or seed.')
  if(length(confidence_level)!=1L||!is.finite(confidence_level)||confidence_level<=0||confidence_level>=1)
    stop('Invalid confidence level.')
  if(!is.logical(keep_signals)||length(keep_signals)!=1L||is.na(keep_signals)) stop('Invalid keep_signals.')
  allowed <- c('primitive_learner','spline_df','super_learner_candidates','primitive_cv_folds',
    'probability_floor','instrument_floor','gamma_C','gamma_cv_folds','gamma_cv_repeats',
    'gamma_backend','gamma_rank','gamma_fixed_lambda','gamma_bandwidth_sample',
    'compatibility_action','projection_warning')
  if(!is.list(control)||any(!names(control)%in%allowed)) stop('Unknown binary estimator control.')
  action <- if(is.null(control$compatibility_action)) 'project' else
    match.arg(control$compatibility_action,c('project','none','stop'))
  control$compatibility_action <- if(action=='none') 'project' else action
  ctl <- .fiv_control(control,J); ctl$compatibility_action <- action
  had_seed <- exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) old_seed <- get('.Random.seed',envir=.GlobalEnv)
  on.exit(if(had_seed) assign('.Random.seed',old_seed,envir=.GlobalEnv) else
    if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv))
  set.seed(seed)
  d <- data.frame(Y=data[[outcome]],T_obs=data[[treatment]],Z=match(as.character(data[[instrument]]),labels))
  x <- if(length(covariates)) paste0('.x',seq_along(covariates)) else character()
  d[x] <- data[covariates]
  if(is.null(fold_id)) {
    fold_id <- integer(n)
    for(j in seq_len(J)) for(t in 0:1) {
      ii <- which(d$Z==j & d$T_obs==t)
      if(length(ii)<folds) stop('Each instrument/treatment cell needs at least folds observations.')
      fold_id[ii] <- sample(rep(seq_len(folds),length.out=length(ii)))
    }
  }
  if(!is.numeric(fold_id)||length(fold_id)!=n||anyNA(fold_id)||
     !setequal(unique(fold_id),seq_len(folds))) stop('Invalid fold_id.')
  S <- plugin <- matrix(NA_real_,n,4,dimnames=list(NULL,c('mu0_lower','mu0_upper','mu1_lower','mu1_upper')))
  diagnostics <- tuning <- list(); warnings <- character()
  result <- withCallingHandlers({
    for(k in seq_len(folds)) {
      tr <- d[fold_id!=k,,drop=FALSE]; va <- d[fold_id==k,,drop=FALSE]
      counts <- table(factor(tr$Z,levels=seq_len(J)),factor(tr$T_obs,levels=0:1))
      minimum <- if(is.null(ctl$gamma_fixed_lambda)) ctl$gamma_cv_folds else 2L
      if(any(counts<minimum)) stop('Insufficient training observations in an instrument/treatment cell.')
      prim <- .fiv_primitives(tr,J,x,ctl,seed+10000*k)
      nt <- prim$predict(tr); nv <- prim$predict(va)
      tuning[[paste0('primitive_',k)]] <- prim$tuning
      for(t in 0:1) {
        base <- .biv_base(tr,va,t,x,ctl,seed+10000*k+100*t)
        if(length(base$p)!=nrow(va)||any(!is.finite(base$p))) stop('Invalid baseline outcome predictions.')
        b <- pmin(1-1e-4,pmax(1e-4,base$p)); p <- matrix(b,nrow(va),J)
        tuning[[paste('baseline',k,t,sep='_')]] <- base$tuning
        norm <- clipping <- numeric()
        for(j in 2:J) {
          pair <- tr$Z%in%c(1L,j); pd <- tr[pair,,drop=FALSE]; pd$Z <- as.integer(pd$Z==j)
          pn <- list(e=nt$e[pair,1]/(nt$e[pair,1]+nt$e[pair,j]),pi0=nt$pi[pair,1],pi1=nt$pi[pair,j])
          gf <- .biv_gamma(pd,pn,t,x,ctl,seed+10000*k+1000*t+j)
          zero <- one <- va; zero$Y <- 0; one$Y <- 1
          rc <- binary_rkhs_gamma_reconstruct(b,gf$predict(zero)$normalized,gf$predict(one)$normalized)
          p[,j] <- rc$a1; norm <- c(norm,abs(rc$normalizer-1))
          clipping <- c(clipping,rc$lower_probability_clipping_rate+rc$upper_probability_clipping_rate)
          tuning[[paste('gamma',k,t,j,sep='_')]] <- .fiv_compact(gf)
        }
        r <- if(t==1) nv$pi else 1-nv$pi
        proj <- .biv_project(r,p,L,action); v <- .biv_bounds(r,proj$p,L)
        keys <- paste0('mu',t,c('_lower','_upper'))
        S[fold_id==k,keys] <- .biv_signal(va,t,nv$e,r,proj$p,L)
        plugin[fold_id==k,keys] <- cbind(v$lower$value,v$upper$value)
        diagnostics[[length(diagnostics)+1L]] <- data.frame(fold=k,treatment=t,
          incompatibility_rate=proj$rate,maximum_crossing=proj$max_gap,
          projection_rate=proj$projection_rate,maximum_shrinkage=proj$max_shrinkage,
          baseline_probability_clip=mean(b!=base$p),reconstructed_probability_clip=mean(clipping),
          conditional_normalization_mae=mean(norm),lower_tie=mean(v$lower$tie),upper_tie=mean(v$upper$tie),
          instrument_clip=nv$instrument_clip_rate,treatment_clip=nv$treatment_clip_rate)
      }
    }
    add <- function(s)cbind(s,ate_lower=s[,'mu1_lower']-s[,'mu0_upper'],ate_upper=s[,'mu1_upper']-s[,'mu0_lower'])
    S <- add(S); plugin <- add(plugin)
    if(any(!is.finite(c(S,plugin)))) stop('Nonfinite fitted signals.')
    est <- colMeans(S); se <- apply(S,2,sd)/sqrt(n); z <- qnorm(1-(1-confidence_level)/2)
    endpoints <- data.frame(endpoint=names(est),estimate=unname(est),standard_error=unname(se),
      interval_lower=unname(est-z*se),interval_upper=unname(est+z*se))
    sets <- do.call(rbind,lapply(c('mu0','mu1','ate'),function(a) {
      lo <- paste0(a,'_lower'); hi <- paste0(a,'_upper')
      data.frame(parameter=a,lower=unname(est[lo]-z*se[lo]),upper=unname(est[hi]+z*se[hi]),
        point_lower=unname(est[lo]),point_upper=unname(est[hi]),confidence_level=confidence_level)
    }))
    sets$empty <- sets$lower>sets$upper; diag <- do.call(rbind,diagnostics)
    if(any(diag$incompatibility_rate>0)) warning('Estimated incompatibility recorded; inspect diagnostics. This is not a formal test.')
    if(any(diag$projection_rate>ctl$projection_warning)) warning('Substantial compatibility stabilization; inference is not certified by projection.')
    list(endpoints=endpoints,outward_sets=sets,plugin=colMeans(plugin),diagnostics=diag,
      tuning=tuning,fold_id=fold_id,instrument_levels=labels,Lambda=setNames(L,labels),
      control=ctl[intersect(names(ctl),allowed)],n=n,seed=seed,
      signals=if(keep_signals) S else NULL,
      inference=list(target='Binary finite-IV sharp endpoints; selected-branch one-step estimator.',
        conditions='Requires nuisance convergence rates, branch-margin/no-tie and influence-function stabilization conditions; projections do not ensure validity.',
        coverage='Pointwise in Lambda and separately for each parameter, not simultaneous bands.',
        units='Independent observations; no clustered or sampling-weighted inference.',
        diagnostics='Clipping, normalization and projection are numerical diagnostics, not formal compatibility tests.'))
  },warning=function(w) warnings <<- c(warnings,conditionMessage(w)))
  result$warnings <- unique(warnings)
  .ivmsm_result(result,'binary outcome, finite-category IV')
}
