# Isolated scalar-X factorial pilot. Source from repository/package root.
source('R/sharp_bounds_plugin.R')
source('R/binary_outcome_bounds.R') # Shared primitive Super Learner implementation.
source('R/restricted_scalar_study.R')
source('R/scalar_iv_fixed_spline_population.R')

sfp_contract <- function() list(version='scalar_factorial_spline_lpc_v1',
  a=c(0,.55,1.1,1.65),b=c(0,.025,.075,.15),n=3000L,repetitions=50L,
  Lambda=c(2,3),folds=5L,threshold_class='fixed_spline5',
  planned_second_class='fixed_rbf8',gamma_backend='exact',
  gamma_C=10^seq(-1.5,1,by=.5),gamma_bandwidth_multiplier=1/sqrt(2),
  gamma_cv_folds=3L,gamma_cv_repeats=2L,
  modes=c('oracle_all','oracle_gamma','rkhs'),
  primitive_library=c('glm_main','glm_interactions','glmnet_interactions','ranger','xgboost'),
  target='unpenalized fixed-spline pseudo-bound',inference=FALSE)
sfp_tasks <- function() {
  z <- sfp_contract()
  g <- expand.grid(a=z$a,b=z$b,KEEP.OUT.ATTRS=FALSE)
  g$design <- seq_len(nrow(g))
  d <- g[rep(seq_len(nrow(g)),each=z$repetitions),]
  d$replication <- rep(seq_len(z$repetitions),nrow(g))
  d$seed <- 270000L+1000L*d$design+d$replication
  d$task <- seq_len(nrow(d));rownames(d) <- NULL;d
}
sfp_generate <- function(n,a,b,seed) {
  set.seed(seed)
  X <- qnorm(pnorm(-2)+runif(n)*siv_normalizer(0))
  Z <- rbinom(n,1,.5)
  pi0 <- plogis(.5*X-a/2);pi1 <- plogis(.5*X+a/2)
  T_obs <- rbinom(n,1,ifelse(Z==1,pi1,pi0))
  eta <- Z*(2*T_obs-1)*b
  r <- eta+qnorm(pnorm(-2-eta)+runif(n)*siv_normalizer(eta))
  Y <- .65*X+.75*T_obs+r
  data.frame(Y=Y,T_obs=T_obs,Z=Z,X=X,e_true=.5,pi0_true=pi0,pi1_true=pi1,
    gamma0_oracle=siv_gamma(Y-.65*X,-b),
    gamma1_oracle=siv_gamma(Y-.65*X-.75,b))
}
sfp_signature <- function() {
  f <- c('R/sharp_bounds_plugin.R','R/binary_outcome_bounds.R','R/restricted_scalar_study.R',
    'R/scalar_iv_population_experiment.R','R/scalar_iv_fixed_spline_population.R',
    'R/scalar_factorial_pilot.R','scripts/run_scalar_factorial_pilot.R',
    'results/scalar_iv_population/fixed_spline_pilot_v1/arm_diagnostics.csv')
  stopifnot(all(file.exists(f)))
  list(contract=sfp_contract(),source_md5=tools::md5sum(f),R=as.character(getRversion()),
    packages=sapply(c('SuperLearner','glmnet','ranger','xgboost','nnls','nloptr'),function(p)
      if(requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else NA_character_))
}
sfp_atomic <- function(object,path) {
  tmp <- tempfile('checkpoint_',dirname(path));saveRDS(object,tmp)
  if(!file.rename(tmp,path)) stop('Could not save checkpoint: ',path)
}
sfp_preflight <- function() {
  packages <- c('SuperLearner','glmnet','ranger','xgboost','nnls','nloptr')
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop('Missing R packages: ',paste(missing,collapse=', '))
  stopifnot(nrow(sfp_tasks())==800)
  d <- sfp_generate(200,1.1,.15,2)
  stopifnot(all(abs(d$Y-.65*d$X-.75*d$T_obs)<=2),all(abs(d$X)<=2))
  target <- read.csv('results/scalar_iv_population/fixed_spline_pilot_v1/arm_diagnostics.csv')
  stopifnot(nrow(target)==128,all(is.finite(target$target)))
  cat('Preflight passed: 800 tasks; exact RKHS; repeated 2 x 3-fold gamma CV.\n')
}
sfp_run <- function(index,smoke=FALSE) {
  sfp_preflight();z <- sfp_contract();tasks <- sfp_tasks()
  stopifnot(length(index)==1,index %in% tasks$task)
  task <- tasks[index,];n <- if(smoke) 600L else z$n
  root <- file.path('results',z$version,if(smoke) 'smoke' else 'pilot')
  path <- file.path(root,sprintf('task_%04d',index))
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  # Exclusive task lock; never delete another process's lock automatically.
  lock <- file.path(path,'running.lock')
  if(!dir.create(lock,showWarnings=FALSE)) stop('Task already locked: ',lock)
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  signature <- sfp_signature()
  meta <- file.path(path,'signature.rds')
  if(file.exists(meta)) stopifnot(identical(readRDS(meta),signature)) else sfp_atomic(signature,meta)
  data <- sfp_generate(n,task$a,task$b,task$seed)
  folds <- make_tz_stratified_folds(data,K=z$folds,seed=task$seed)
  # Full current tuning rule is retained in smoke; only n changes.
  gamma <- make_cached_gamma_rkhs_provider(lambda_grid=default_gamma_lambda_grid,
    bandwidth_multipliers=z$gamma_bandwidth_multiplier,cv_folds=z$gamma_cv_folds,
    cv_repeats=z$gamma_cv_repeats,selection_rule='minimum',backend='exact')
  primitive <- new.env(parent=emptyenv())
  provider <- function(fold,train,x_cols,seed,...) {
    key <- paste0('fold',fold)
    if(exists(key,primitive,inherits=FALSE)) return(get(key,primitive))
    value <- fit_nuisance_models(train,x_cols,positivity_epsilon=.05,
      learner='super_learner',treatment_spline_df=5,
      super_learner_candidates=z$primitive_library,super_learner_cv_folds=3,
      super_learner_ranger_num_trees=500,library_seed=seed)
    assign(key,value,primitive);value
  }
  targets <- read.csv('results/scalar_iv_population/fixed_spline_pilot_v1/arm_diagnostics.csv')
  statuses <- list()
  for(L in z$Lambda) {
    file <- file.path(path,paste0('Lambda',L,'.rds'))
    if(file.exists(file)) {
      old <- readRDS(file)
      stopifnot(identical(old$signature,signature))
      if(isTRUE(old$success)) {statuses[[as.character(L)]] <- old$status;next}
    }
    started <- proc.time()[3]
    result <- tryCatch({
      qstore <- new.env(parent=emptyenv())
      fit <- estimate_sharp_bounds_plugin(data,
        Lambda0_T1=L,Lambda1_T1=L,Lambda0_T0=L,Lambda1_T0=L,
        x_cols='X',K=z$folds,fold_id=folds,nuisance_learner='super_learner',
        nuisance_fit_provider_factory=provider,gamma_fit_provider_factory=gamma,
        gamma_rkhs_backend='exact',rkhs_cv_repeats=2,
        threshold_learner='external',
        bound_component_provider_factory=restricted_scalar_provider('fixed_spline5',qstore),
        conditional_mean_spline_df=5,estimated_compatibility_action='project',
        max_gamma_projection_rate=.10,estimation_modes=z$modes,
        compute_oracle_eif=FALSE,compute_estimated_eif=FALSE,seed=task$seed)
      # The shared estimator's latent-DGP truth defaults DO NOT apply here.
      fit$truth <- NULL
      pop <- targets[targets$a==task$a & targets$b==task$b & targets$Lambda==L,]
      target <- setNames(pop$target,pop$endpoint)
      target <- c(target,ate_lower=target[['mu1_lower']]-target[['mu0_upper']],
        ate_upper=target[['mu1_upper']]-target[['mu0_lower']])
      endpoints <- do.call(rbind,lapply(seq_len(nrow(fit$pseudo_summary)),function(i) {
        s <- fit$pseudo_summary[i,];direct <- fit$summary[fit$summary$mode==s$mode,]
        data.frame(mode=s$mode,endpoint=names(target),target=unname(target),
          estimate=as.numeric(s[1,names(target)]),
          direct_estimate=as.numeric(direct[1,names(target)]))
      }))
      endpoints$error <- endpoints$estimate-endpoints$target
      endpoints$direct_minus_pseudo <- endpoints$direct_estimate-endpoints$estimate
      stopifnot(nrow(endpoints)==18,all(is.finite(endpoints$estimate)),
        setequal(unique(endpoints$mode),z$modes))
      qm <- as.list(qstore)
      qdiag <- do.call(rbind,lapply(qm,function(m) do.call(rbind,lapply(c('lower','upper'),function(bd) {
        q <- m[[bd]]
        data.frame(mode=m$mode,fold=m$fold,treatment=m$treatment,bound=bd,
          convergence=q$convergence,duality_gap=q$duality_gap,
          tolerance=q$certificate_tolerance,penalty=q$penalty)
      }))))
      stopifnot(nrow(qdiag)==60)
      gamma_cv <- do.call(rbind,lapply(as.list(attr(gamma,'gamma_cache')),function(entry) {
        cv <- entry$fit$cv_resample_table
        cv$outer_fold <- entry$fold;cv$treatment <- entry$treatment
        cv$n_training_arm <- sum(entry$training_signature$T_obs==entry$treatment)
        cv$C <- cv$lambda*sqrt(cv$n_training_arm)
        cv
      }))
      stopifnot(nrow(gamma_cv)==360,setequal(unique(gamma_cv$cv_repeat),1:2))
      list(success=TRUE,error='',endpoints=endpoints,fit=fit,qdiagnostics=qdiag,
        gamma_cv_resamples=gamma_cv,
        qmodels=qm,threshold_converged=all(qdiag$convergence==0))
    },error=function(e) list(success=FALSE,error=conditionMessage(e),threshold_converged=FALSE))
    result$signature <- signature
    result$status <- cbind(task,n=n,Lambda=L,smoke=smoke,success=result$success,
      threshold_converged=result$threshold_converged,
      elapsed_seconds=unname(proc.time()[3]-started),error=result$error)
    sfp_atomic(result,file)
    statuses[[as.character(L)]] <- result$status
    print(result$status);flush.console()
  }
  status <- do.call(rbind,statuses)
  sfp_atomic(status,file.path(path,'status.rds'))
  capture.output(sessionInfo(),file=file.path(path,'session_info.txt'))
  if(!all(status$success)) stop('Task has failed configurations; checkpoints preserved.')
  invisible(status)
}
sfp_collect <- function(smoke=FALSE) {
  root <- file.path('results',sfp_contract()$version,if(smoke) 'smoke' else 'pilot')
  files <- list.files(root,pattern='^Lambda[23][.]rds$',recursive=TRUE,full.names=TRUE)
  if(!length(files)) stop('No checkpoints found.')
  statuses <- endpoints <- qdiag <- folds <- gamma <- cv <- list();signature <- NULL
  tag_rows <- function(tag,d) {
    if(is.null(d) || !nrow(d)) return(data.frame())
    tags <- tag[rep(1,nrow(d)),,drop=FALSE];rownames(tags) <- NULL
    rownames(d) <- NULL;cbind(tags,d)
  }
  for(file in files) {
    r <- readRDS(file)
    if(is.null(signature)) signature <- r$signature else stopifnot(identical(signature,r$signature))
    statuses[[file]] <- r$status
    if(!r$success) next
    tag <- r$status[c('task','a','b','n','seed','replication','Lambda')]
    endpoints[[file]] <- tag_rows(tag,r$endpoints)
    qdiag[[file]] <- tag_rows(tag,r$qdiagnostics)
    folds[[file]] <- tag_rows(tag,r$fit$fold_diagnostics)
    gamma[[file]] <- tag_rows(tag,r$fit$gamma_accuracy_by_fold)
    cv[[file]] <- tag_rows(tag,r$gamma_cv_resamples)
  }
  status <- restricted_scalar_bind(statuses);e <- restricted_scalar_bind(endpoints)
  write.csv(status,file.path(root,'status.csv'),row.names=FALSE)
  write.csv(e,file.path(root,'endpoints.csv'),row.names=FALSE)
  write.csv(restricted_scalar_bind(qdiag),file.path(root,'threshold_diagnostics.csv'),row.names=FALSE)
  write.csv(restricted_scalar_bind(folds),file.path(root,'fold_diagnostics.csv'),row.names=FALSE)
  write.csv(restricted_scalar_bind(gamma),file.path(root,'gamma_accuracy.csv'),row.names=FALSE)
  write.csv(restricted_scalar_bind(cv),file.path(root,'gamma_cv_resamples.csv'),row.names=FALSE)
  if(nrow(e)) {
    keys <- c('a','b','n','Lambda','mode','endpoint')
    summary <- do.call(rbind,lapply(split(e,interaction(e[keys],drop=TRUE)),function(g)
      cbind(g[1,keys],replications=nrow(g),target=g$target[1],bias=mean(g$error),
        rmse=sqrt(mean(g$error^2)),empirical_sd=sd(g$estimate),
        mean_direct_minus_pseudo=mean(g$direct_minus_pseudo))))
    write.csv(summary,file.path(root,'endpoint_summary.csv'),row.names=FALSE)
    wide <- reshape(e[c('task','a','b','n','Lambda','mode','endpoint','estimate','target')],
      idvar=c('task','a','b','n','Lambda','mode'),timevar='endpoint',direction='wide')
    wide$ate_width <- wide$estimate.ate_upper-wide$estimate.ate_lower
    wide$target_width <- wide$target.ate_upper-wide$target.ate_lower
    wide$ate_endpoints_cross <- wide$ate_width<0
    write.csv(wide,file.path(root,'widths.csv'),row.names=FALSE)
  }
  cat(nrow(status),'configurations;',sum(status$success),'successful;',
    sum(!status$success),'failed;',sum(!status$threshold_converged),'uncertified threshold fits.\n')
  cat('Expected configurations:',if(smoke) 6 else 1600,
      '; each dataset has two Lambda configurations.\n')
}
