# Cross-fitted general-J extension of ref_fiv_run. No DGP, target, or file paths.
.fiv_control <- function(control, J) {
  defaults <- list(primitive_learner='glm_spline', spline_df=3L,
    mean_learner='spline_fixed', mean_control=NULL,
    super_learner_candidates=binary_default_super_learner_candidates(),
    primitive_cv_folds=3L, probability_floor=.05,
    instrument_floor=min(.005,.5/J), mean_margin=1e-4,
    gamma_C=10^seq(-1.5,1,.5), gamma_cv_folds=3L, gamma_cv_repeats=2L,
    gamma_backend='nystrom', gamma_rank=NULL, gamma_fixed_lambda=NULL,
    gamma_bandwidth_sample=1000L, threshold_maxit=60000L,
    compatibility_action='project', projection_warning=.05)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!nzchar(names(control))) || anyDuplicated(names(control)))))
    stop('control must be a uniquely named list.')
  unknown <- setdiff(names(control),names(defaults))
  if (length(unknown)) stop('Unknown control: ',paste(unknown,collapse=', '))
  defaults[names(control)] <- control; c <- defaults
  c$primitive_learner <- match.arg(c$primitive_learner,c('glm_spline','super_learner','library_cv','gam','ranger','glmnet'))
  c$mean_learner <- match.arg(c$mean_learner,c('spline_fixed','spline_cv','gam','ranger','library_cv'))
  c$gamma_backend <- match.arg(c$gamma_backend,c('nystrom','exact'))
  c$compatibility_action <- match.arg(c$compatibility_action,c('project','stop'))
  for (key in c('spline_df','primitive_cv_folds','gamma_cv_folds','gamma_cv_repeats','gamma_bandwidth_sample','threshold_maxit')) {
    x <- c[[key]]
    if (length(x)!=1L || !is.finite(x) || x!=as.integer(x) || x<1) stop('Invalid control: ',key)
  }
  if(c$gamma_cv_folds<2 || c$primitive_cv_folds<2) stop('CV requires at least two folds.')
  for(key in c('probability_floor','instrument_floor','mean_margin','projection_warning'))
    if(length(c[[key]])!=1L || !is.finite(c[[key]]) || c[[key]]<=0 || c[[key]]>=1)
      stop('Invalid control: ',key)
  if(c$probability_floor>=.5 || c$instrument_floor>=1/J) stop('Invalid probability floors.')
  if(!is.numeric(c$gamma_C) || !length(c$gamma_C) || any(!is.finite(c$gamma_C)) || any(c$gamma_C<=0)) stop('Invalid gamma_C.')
  if(!is.null(c$gamma_fixed_lambda) && (length(c$gamma_fixed_lambda)!=1L ||
      !is.finite(c$gamma_fixed_lambda) || c$gamma_fixed_lambda<=0)) stop('Invalid fixed gamma penalty.')
  if(!is.null(c$gamma_rank) && (length(c$gamma_rank)!=1L || !is.finite(c$gamma_rank) ||
      c$gamma_rank<1 || c$gamma_rank!=as.integer(c$gamma_rank))) stop('Invalid gamma rank.')
  c
}
.fiv_compact <- function(fit) {
  fields <- c('selected_lambda','selected_C','base_bandwidth','selected_bandwidth',
    'selected_bandwidth_multiplier','normalization_constant','training_negative_rate',
    'rkhs_backend','nystrom_rank_requested','nystrom_landmark_count','nystrom_effective_rank',
    'cv_table','cv_resample_table','ensemble_weights','converged','rank','configuration',
    'selected_learners','library_configuration')
  fit[intersect(fields,names(fit))]
}
.fiv_primitives <- function(tr,J,x,c,seed) {
  binary <- function(d,y,s) {
    if (!length(x)) {
      p <- mean(y)
      return(list(predict=function(nd) rep(p,nrow(nd))))
    }
    d$.response <- y
    binary_fit_regression_with_learner(d,'.response',x,spline_df=c$spline_df,
      learner=c$primitive_learner,library_cv_folds=c$primitive_cv_folds,
      super_learner_candidates=c$super_learner_candidates,
      super_learner_cv_folds=c$primitive_cv_folds,seed=s)
  }
  hazards <- lapply(seq_len(J-1L),function(j) {
    use <- tr$Z>=j; binary(tr[use,,drop=FALSE],as.integer(tr$Z[use]==j),seed+j)
  })
  propensity <- lapply(seq_len(J),function(j) {
    use <- tr$Z==j; binary(tr[use,,drop=FALSE],tr$T_obs[use],seed+100L+j)
  })
  predict <- function(d) {
    e <- matrix(0,nrow(d),J); remaining <- rep(1,nrow(d))
    for(j in seq_len(J-1L)) {
      h <- pmin(1,pmax(0,hazards[[j]]$predict(d)))
      e[,j] <- remaining*h; remaining <- remaining*(1-h)
    }
    e[,J] <- remaining; mix <- rep(0,nrow(d))
    for(j in seq_len(J)) {
      low <- e[,j]<c$instrument_floor
      mix[low] <- pmax(mix[low],(c$instrument_floor-e[low,j])/(1/J-e[low,j]))
    }
    e <- e*(1-mix)+mix/J
    rawpi <- do.call(cbind,lapply(propensity,function(f) f$predict(d)))
    pi <- rawpi
    pi[] <- pmin(1-c$probability_floor,pmax(c$probability_floor,rawpi))
    if(any(!is.finite(c(e,pi)))) stop('Nonfinite primitive predictions.')
    list(e=e,pi=pi,instrument_clip_rate=mean(mix>0),
         treatment_clip_rate=mean(abs(pi-rawpi)>1e-12))
  }
  list(predict=predict,tuning=list(instrument=lapply(hazards,.fiv_compact),
                                  treatment=lapply(propensity,.fiv_compact)))
}
.fiv_means <- function(d,values,x,c,seed) {
  values <- as.matrix(values)
  if(!length(x)) {
    m <- colMeans(values)
    return(list(predict=function(nd) matrix(rep(m,each=nrow(nd)),nrow(nd)),tuning=list(intercept=TRUE)))
  }
  names <- paste0('.response',seq_len(ncol(values)))
  d[names] <- as.data.frame(values)
  fit <- fit_eif_conditional_means(d,names,x,learner=c$mean_learner,
    spline_df=c$spline_df,control=c$mean_control,seed=seed)
  list(predict=fit$predict,tuning=.fiv_compact(fit))
}
.fiv_project <- function(g,r,L,action) {
  A <- r+sweep(1-r,2,L,'/'); B <- r+sweep(1-r,2,L,'*')
  alpha <- rep(1,nrow(r))
  for(j in seq_len(ncol(r))) for(k in seq_len(ncol(r))) {
    slope <- A[,j]*(g[,j]-1)-B[,k]*(g[,k]-1); use <- slope>0
    alpha[use] <- pmin(alpha[use],(B[use,k]-A[use,j])/slope[use])
  }
  alpha <- pmax(0,pmin(1,alpha)); bad <- alpha<1-1e-10
  if(action=='stop' && any(bad)) stop('Estimated density envelopes cross; compatibility_action=stop.')
  alpha[bad] <- alpha[bad]*(1-1e-8)
  list(gamma=1+(g-1)*alpha,alpha=alpha)
}
.fiv_components <- function(y,q,r,g,L,direction) {
  A <- r+sweep(1-r,2,L,'/'); B <- r+sweep(1-r,2,L,'*')
  lo <- A*g; up <- B*g
  il <- max.col(lo,ties.method='first'); iu <- max.col(-up,ties.method='first')
  idx <- seq_along(y); lower <- lo[cbind(idx,il)]; upper <- up[cbind(idx,iu)]
  high <- if(direction=='lower') y<=q else y>q
  j <- ifelse(high,iu,il); ix <- cbind(idx,j)
  coef <- ifelse(high,B[ix],A[ix]); deriv <- ifelse(high,1-L[j],1-1/L[j])
  H <- C <- matrix(0,length(y),ncol(r))
  H[ix] <- (y-q)*coef; C[ix] <- (y-q)*deriv
  list(lower=lower,upper=upper,rho=ifelse(high,upper,lower),H=H,C=C,
    lower_tie=rowSums(abs(lo-lower)<1e-10)>1L,
    upper_tie=rowSums(abs(up-upper)<1e-10)>1L)
}
.fiv_score <- function(q,Hbar,Cbar,r,e,z,tobs,t,H) {
  ix <- cbind(seq_along(q),z); D <- as.numeric(tobs==t)
  q+rowSums(Hbar)+(Cbar[ix]*(D-r[ix])+D/r[ix]*(H[ix]-Hbar[ix]))/e[ix]
}
.fiv_threshold <- function(F,y,level,weight,n_training,maxit) {
  if(any(!is.finite(c(F,y,level,weight))) || any(weight<0) ||
     any(level<=0 | level>=1)) stop('Invalid threshold regression inputs; inspect envelope diagnostics.')
  if(sum(weight)==0) return(list(coefficient=c(median(y),rep(0,ncol(F)-1L)),
    iterations=0L,convergence=0L,primal=0,dual=0,flat_objective=TRUE))
  if(qr(F)$rank<ncol(F)) stop('Threshold basis rank deficient in a training reference arm; class not changed.')
  # The fixed basis includes an intercept. Whiten only its nonconstant columns.
  features <- F[,-1L,drop=FALSE]
  if(ncol(features)) {
    center <- colMeans(features); centered <- sweep(features,2,center)
    s <- svd(centered,nu=0)
    if(min(s$d)<max(s$d)*1e-12) stop('Ill-conditioned threshold basis.')
    rotation <- sweep(s$v,2,s$d/sqrt(nrow(F)),'/')
    features <- centered %*% rotation
  }
  for(rho in c(.01,.1,1,.001)) {
    f <- fit_weighted_quantile_ridge_admm(features,y,level,weight,1/n_training,
                                         rho=rho,maxit=maxit,abs_tol=1e-5,rel_tol=1e-4)
    if(f$convergence==0L) break
  }
  if(f$convergence!=0L) stop('Threshold optimizer did not converge.')
  coef <- f$coefficients
  if(ncol(features)) {
    raw <- drop(rotation %*% coef[-1L]); coef <- c(coef[1L]-sum(center*raw),raw)
  }
  list(coefficient=coef,iterations=f$iterations,convergence=f$convergence,
       primal=f$primal_residual,dual=f$dual_residual)
}

.fiv_spline_basis <- function(covariates,spec=NULL,reference=NULL) {
  if(!length(covariates)||anyDuplicated(covariates)) stop('Supply distinct spline covariate names.')
  if(!is.null(spec)&&!is.null(reference)) stop('Supply spline_spec or spline_reference, not both.')
  if(!is.null(reference)) {
    if(!is.data.frame(reference)||!all(covariates%in%names(reference)))
      stop('spline_reference must contain the covariates; use independent reference data.')
    spec <- setNames(lapply(covariates,function(nm) {
      v <- reference[[nm]]
      if(!is.numeric(v)||!length(v)||any(!is.finite(v))||length(unique(v))<2)
        stop('Reference covariates must be finite numeric and nonconstant.')
      if(all(v%in%c(0,1))) return(list(linear=TRUE))
      boundary <- range(v)
      knots <- unique(as.numeric(stats::quantile(v,c(.25,.5,.75),names=FALSE)))
      list(knots=knots[knots>boundary[1]&knots<boundary[2]],boundary=boundary)
    }),covariates)
  }
  if(!is.list(spec)||is.null(names(spec))||anyDuplicated(names(spec))||!setequal(names(spec),covariates))
    stop('For splines, supply a named spline_spec or independent spline_reference data.')
  spec <- spec[covariates]
  for(s in spec) {
    if(!is.list(s)) stop('Each spline specification must be a list.')
    if(isTRUE(s$linear)) next
    b <- s$boundary; k <- s$knots
    if(!is.numeric(b)||length(b)!=2L||any(!is.finite(b))||b[1]>=b[2]||
       !is.numeric(k)||any(!is.finite(k))||anyDuplicated(k)||is.unsorted(k)||any(k<=b[1]|k>=b[2]))
      stop('Supply increasing boundary values and distinct interior knots (numeric(0) is allowed).')
  }
  force(spec);force(covariates)
  function(x) {
    if(!all(covariates%in%names(x))) stop('Missing spline covariates.')
    pieces <- lapply(covariates,function(nm) {
      v <- x[[nm]];s <- spec[[nm]]
      if(!is.numeric(v)||any(!is.finite(v))) stop('Spline covariates must be finite numeric.')
      if(isTRUE(s$linear)) return(matrix(v,ncol=1))
      splines::ns(v,knots=s$knots,Boundary.knots=s$boundary)
    })
    do.call(cbind,c(list(rep(1,nrow(x))),pieces))
  }
}

estimate_finite_iv_continuous <- function(data,outcome,treatment,instrument,
    covariates=character(),Lambda=2,reference=NULL,folds=5L,fold_id=NULL,
    threshold_basis=NULL,control=list(),confidence_level=.95,seed=1L,
    keep_signals=TRUE,threshold_class=c('linear','spline'),
    spline_spec=NULL,spline_reference=NULL) {
  threshold_class <- match.arg(threshold_class)
  if(threshold_class=='spline') {
    if(!is.null(threshold_basis)) stop('Use threshold_class or threshold_basis, not both.')
    threshold_basis <- .fiv_spline_basis(covariates,spline_spec,spline_reference)
  } else if(!is.null(spline_spec)||!is.null(spline_reference))
    stop('Spline controls require threshold_class="spline".')
  if(!all(vapply(list(outcome,treatment,instrument),function(v)
       is.character(v)&&length(v)==1L&&!is.na(v)&&nzchar(v),logical(1))) ||
     !is.character(covariates)) stop('Supply outcome, treatment, instrument, and covariates as column names.')
  cols <- c(outcome,treatment,instrument,covariates)
  if(!is.data.frame(data) || !nrow(data) || !all(cols %in% names(data)) ||
     anyDuplicated(cols) || anyNA(data[,cols,drop=FALSE])) stop('Supply complete data and distinct column names.')
  if(!is.numeric(data[[outcome]]) || any(!is.finite(data[[outcome]])) ||
     !all(data[[treatment]] %in% c(0,1))) stop('Outcome must be finite numeric; treatment must be 0/1.')
  if(length(covariates) && !all(vapply(data[covariates],function(v) is.numeric(v)&&all(is.finite(v)),logical(1))))
    stop('Covariates must be finite numeric columns; encode categorical covariates explicitly.')
  z <- factor(data[[instrument]]); labels <- levels(z); J <- length(labels); n <- nrow(data)
  if(J<2L) stop('At least two IV categories required.')
  if(!is.null(names(Lambda))) {
    if(length(Lambda)!=J || anyDuplicated(names(Lambda)) || !setequal(names(Lambda),labels))
      stop('Named Lambda must have exactly one value for each instrument label.')
    Lambda <- Lambda[labels]
  }
  L <- .iv_lambda(Lambda,J)
  if(is.null(reference)) reference <- labels[1L]
  if(length(reference)!=1L || !as.character(reference) %in% labels) stop('reference must be an observed IV label.')
  order <- c(match(as.character(reference),labels),which(labels!=as.character(reference)))
  labels <- labels[order]; L <- L[order]
  if(length(folds)!=1 || !is.finite(folds) || folds<2 || folds!=as.integer(folds)) stop('Invalid folds.')
  if(length(confidence_level)!=1 || !is.finite(confidence_level) || confidence_level<=0 || confidence_level>=1) stop('Invalid confidence level.')
  if(length(seed)!=1 || !is.finite(seed) || seed<0 || seed!=as.integer(seed) || seed>.Machine$integer.max/100) stop('Invalid seed.')
  if(!is.logical(keep_signals) || length(keep_signals)!=1 || is.na(keep_signals)) stop('Invalid keep_signals.')
  ctl <- .fiv_control(control,J)
  had_seed <- exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) old_seed <- get('.Random.seed',envir=.GlobalEnv)
  on.exit(if(had_seed) assign('.Random.seed',old_seed,envir=.GlobalEnv) else
    if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv))
  set.seed(seed)
  d <- data.frame(Y=data[[outcome]],T_obs=as.integer(data[[treatment]]),Z=match(as.character(z),labels))
  x <- paste0('.x',seq_along(covariates)); if(!length(covariates)) x <- character()
  d[x] <- data[covariates]
  if(is.null(fold_id)) {
    fold_id <- integer(n)
    for(j in seq_len(J)) for(t in 0:1) {
      ids <- which(d$Z==j & d$T_obs==t)
      if(length(ids)<folds) stop('Each IV/treatment cell needs at least folds observations.')
      fold_id[ids] <- sample(rep(seq_len(folds),length.out=length(ids)))
    }
  }
  if(!is.numeric(fold_id) || length(fold_id)!=n || anyNA(fold_id) || !setequal(unique(fold_id),seq_len(folds))) stop('Invalid fold_id.')
  # Caller basis is evaluated on covariates only; it must be prespecified and rowwise.
  X <- data[covariates]
  F <- if(is.null(threshold_basis)) cbind(1,as.matrix(X)) else {
    if(!is.function(threshold_basis)) stop('threshold_basis must be a function.')
    as.matrix(threshold_basis(X))
  }
  if(!is.numeric(F) || nrow(F)!=n || ncol(F)<1 || any(!is.finite(F)) || any(abs(F[,1]-1)>1e-12))
    stop('Threshold basis must return a finite n-row numeric matrix with first column one.')
  keys <- c('mu0_lower','mu0_upper','mu1_lower','mu1_upper')
  signals <- direct <- pseudo <- matrix(NA_real_,n,4,dimnames=list(NULL,keys))
  diagnostics <- tuning <- list(); warnings <- character()
  result <- withCallingHandlers({
    for(k in seq_len(folds)) {
      use <- fold_id!=k; tr <- d[use,,drop=FALSE]; va <- d[!use,,drop=FALSE]
      if(!is.null(threshold_basis)) for(part in list(use,!use)) {
        local_basis <- as.matrix(threshold_basis(X[part,,drop=FALSE]))
        if(!identical(dim(local_basis),dim(F[part,,drop=FALSE])) ||
           any(!is.finite(local_basis)) || max(abs(local_basis-F[part,,drop=FALSE]))>1e-10)
          stop('threshold_basis must be a fixed rowwise map, not a data-adaptive basis.')
      }
      counts <- table(factor(tr$Z,levels=seq_len(J)),factor(tr$T_obs,levels=0:1))
      minimum <- if(is.null(ctl$gamma_fixed_lambda)) ctl$gamma_cv_folds else 2L
      if(any(counts<minimum)) stop('Insufficient training cell observations for gamma fitting in fold ',k)
      prim <- .fiv_primitives(tr,J,x,ctl,seed+10000*k)
      nt <- prim$predict(tr); nv <- prim$predict(va)
      tuning[[paste0('primitive_',k)]] <- prim$tuning
      for(t in 0:1) {
        rt <- if(t==1) nt$pi else 1-nt$pi; rv <- if(t==1) nv$pi else 1-nv$pi
        gt <- matrix(1,nrow(tr),J); gv <- matrix(1,nrow(va),J)
        for(j in 2:J) {
          pair <- tr$Z %in% c(1L,j); pd <- tr[pair,,drop=FALSE]; pd$Z <- as.integer(pd$Z==j)
          pn <- list(e=nt$e[pair,1]/(nt$e[pair,1]+nt$e[pair,j]),pi0=nt$pi[pair,1],pi1=nt$pi[pair,j])
          args <- list(train=pd,nuisance_train=pn,t_value=t,x_cols=x,
            backend=ctl$gamma_backend,nystrom_rank=ctl$gamma_rank,
            bandwidth_sample=ctl$gamma_bandwidth_sample,seed=seed+10000*k+1000*t+j)
          gf <- if(is.null(ctl$gamma_fixed_lambda)) do.call(fit_gamma_rkhs,c(args,list(
            lambda_grid=function(n_training,treatment=NULL) ctl$gamma_C/sqrt(n_training),
            bandwidth_multipliers=1/sqrt(2),cv_folds=ctl$gamma_cv_folds,
            cv_repeats=ctl$gamma_cv_repeats,selection_rule='minimum'))) else
            do.call(fit_gamma_rkhs_fixed,c(args,list(lambda=ctl$gamma_fixed_lambda,bandwidth_multiplier=1/sqrt(2))))
          gf$selected_C <- gf$selected_lambda*sqrt(sum(pd$T_obs==t))
          gt[tr$T_obs==t,j] <- gf$predict(tr[tr$T_obs==t,,drop=FALSE])$normalized
          if(any(va$T_obs==t)) gv[va$T_obs==t,j] <- gf$predict(va[va$T_obs==t,,drop=FALSE])$normalized
          tuning[[paste('gamma',k,t,j,sep='_')]] <- .fiv_compact(gf)
        }
        pt <- .fiv_project(gt,rt,L,ctl$compatibility_action)
        pv <- .fiv_project(gv,rv,L,ctl$compatibility_action)
        env <- .fiv_components(tr$Y,rep(0,nrow(tr)),rt,pt$gamma,L,'lower')
        base <- tr$T_obs==t & tr$Z==1L
        mf <- .fiv_means(tr[base,,drop=FALSE],cbind(env$lower[base],env$upper[base]),x,ctl,seed+7000*k+t)
        pred <- as.matrix(mf$predict(tr))
        lb <- pmin(1-ctl$mean_margin,pmax(1e-6,pred[,1])); ub <- pmax(1+ctl$mean_margin,pred[,2])
        tuning[[paste('envelope',k,t,sep='_')]] <- mf$tuning
        for(side in c('lower','upper')) {
          key <- paste0('mu',t,'_',side)
          level <- if(side=='lower') (1-lb)/(ub-lb) else (ub-1)/(ub-lb)
          q <- .fiv_threshold(F[use,,drop=FALSE][base,,drop=FALSE],tr$Y[base],level[base],
            (env$upper-env$lower)[base]/(nt$e[base,1]*rt[base,1]),nrow(tr),ctl$threshold_maxit)
          qt <- drop(F[use,,drop=FALSE]%*%q$coefficient); qv <- drop(F[!use,,drop=FALSE]%*%q$coefficient)
          ct <- .fiv_components(tr$Y,qt,rt,pt$gamma,L,side)
          cv <- .fiv_components(va$Y,qv,rv,pv$gamma,L,side)
          Hbar <- Cbar <- matrix(0,nrow(va),J); mt <- list()
          for(j in seq_len(J)) {
            arm <- tr$T_obs==t & tr$Z==j
            f <- .fiv_means(tr[arm,,drop=FALSE],cbind(ct$H[arm,j],ct$C[arm,j]),x,ctl,seed+30000*k+100*t+j)
            p <- as.matrix(f$predict(va)); Hbar[,j] <- p[,1]; Cbar[,j] <- p[,2]; mt[[j]] <- f$tuning
          }
          signals[!use,key] <- .fiv_score(qv,Hbar,Cbar,rv,nv$e,va$Z,va$T_obs,t,cv$H)
          w <- as.numeric(va$Z==1L & va$T_obs==t)/(nv$e[,1]*rv[,1])
          direct[!use,key] <- w*cv$rho*va$Y
          pseudo[!use,key] <- qv+w*cv$rho*(va$Y-qv)
          tuning[[paste('threshold',k,key,sep='_')]] <- list(threshold=q,conditional_means=mt)
          a <- tr$T_obs==t; b <- va$T_obs==t
          diagnostics[[length(diagnostics)+1L]] <- data.frame(fold=k,endpoint=key,
            training_gamma_projection=mean(pt$alpha[a]<1-1e-10),
            heldout_gamma_projection=if(any(b)) mean(pv$alpha[b]<1-1e-10) else NA_real_,
            maximum_gamma_shrinkage=max(1-pt$alpha[a],1-pv$alpha[b]),
            conditional_mean_projection=mean(abs(lb-pred[,1])>1e-12 | abs(ub-pred[,2])>1e-12),
            normalization_residual=mean(w*cv$rho)-1,
            direct_minus_pseudo=mean(w*cv$rho*qv-qv),
            lower_branch_tie=if(any(b)) mean(cv$lower_tie[b]) else NA_real_,
            upper_branch_tie=if(any(b)) mean(cv$upper_tie[b]) else NA_real_,
            instrument_clip=nv$instrument_clip_rate,treatment_clip=nv$treatment_clip_rate)
        }
      }
    }
    if(any(!is.finite(c(signals,direct,pseudo)))) stop('Nonfinite fitted signals.')
    add_ate <- function(s) cbind(s,ate_lower=s[,'mu1_lower']-s[,'mu0_upper'],ate_upper=s[,'mu1_upper']-s[,'mu0_lower'])
    signals <- add_ate(signals); direct <- add_ate(direct); pseudo <- add_ate(pseudo)
    est <- colMeans(signals); se <- apply(signals,2,sd)/sqrt(n); zcrit <- qnorm(1-(1-confidence_level)/2)
    endpoints <- data.frame(endpoint=names(est),estimate=unname(est),standard_error=unname(se),
      interval_lower=unname(est-zcrit*se),interval_upper=unname(est+zcrit*se))
    sets <- do.call(rbind,lapply(c('mu0','mu1','ate'),function(p) {
      a <- paste0(p,'_lower'); b <- paste0(p,'_upper')
      data.frame(parameter=p,lower=unname(est[a]-zcrit*se[a]),upper=unname(est[b]+zcrit*se[b]),
        point_lower=unname(est[a]),point_upper=unname(est[b]),confidence_level=confidence_level)
    }))
    sets$empty <- sets$lower>sets$upper
    diag <- do.call(rbind,diagnostics)
    if(any(as.matrix(diag[c('training_gamma_projection','heldout_gamma_projection','conditional_mean_projection')])>ctl$projection_warning,na.rm=TRUE))
      warning('Substantial stabilization recorded; inspect diagnostics. This is not a formal compatibility rejection.')
    list(endpoints=endpoints,outward_sets=sets,direct_plugin=colMeans(direct),
      pseudo_plugin=colMeans(pseudo),diagnostics=diag,tuning=tuning,fold_id=fold_id,
      instrument_levels=labels,Lambda=setNames(L,labels),control=ctl,
      threshold_class=if(threshold_class=='spline') 'fixed additive natural cubic splines including intercept' else if(is.null(threshold_basis)) 'fixed linear including intercept' else 'user-supplied fixed basis including intercept',
      threshold_dimension=ncol(F),n=n,seed=seed,
      signals=if(keep_signals) signals else NULL,
      direct_signals=if(keep_signals) direct else NULL,pseudo_signals=if(keep_signals) pseudo else NULL,
      inference=list(target='restricted-class pseudo-bound endpoints',
        outward='Pointwise in Lambda; requires nuisance rates, margin and influence-function stabilization conditions. Not simultaneous across parameters.',
        endpoint='Two-sided intervals additionally require negligible optimization error and appropriate regularity at the optimizer.',
        units='Independent observations only; no clustered or sampling-weighted inference.',
        diagnostics='Projection and clipping are numerical diagnostics, not formal tests or certification of validity.'))
  },warning=function(w) {warnings <<- c(warnings,conditionMessage(w))})
  result$warnings <- unique(warnings)
  .ivmsm_result(result,'continuous outcome, finite-category IV; fixed threshold class')
}
