# Isolated simple/stress X5 pilot. Never changes scalar or production defaults.
source('R/x5_simple_gamma.R')
source('R/scalar_log_spline_compatibility.R')

x5p_shape <- function(d,design) {
  stopifnot(design %in% c('simple','stress'))
  if(design=='simple')return(rep(1,nrow(d)))
  1+.6*sin(d$X1*d$X2)+.4*sin(d$X3*d$X4)
}
x5p_generate <- function(n,a,b,seed,design) {
  # Same random draws, X/Z/T and marginal outcome location in both designs.
  z <- x5s_generate(n,a,b,seed)
  # Recover the exact uniform residual draws without consuming additional RNG.
  loc <- drop(as.matrix(z[,paste0('X',1:5)])%*%x5s_contract()$theta)
  eta0 <- z$Z*(2*z$T_obs-1)*b
  r0 <- z$Y-loc-.75*z$T_obs
  u <- (pnorm(r0-eta0)-pnorm(-2-eta0))/siv_normalizer(eta0)
  shape <- x5p_shape(z,design)
  eta <- eta0*shape
  r <- eta+qnorm(pnorm(-2-eta)+u*siv_normalizer(eta))
  z$Y <- loc+.75*z$T_obs+r
  z$gamma0_oracle <- siv_gamma(z$Y-loc,-b*shape)
  z$gamma1_oracle <- siv_gamma(z$Y-loc-.75,b*shape)
  z
}

# Main natural splines and Y-by-X spline interactions, not a full tensor.
# Joint minimum CV over C and interaction penalty rho, including additive limit.
x5p_log_engine <- function(rho=1) {
  e <- x5s_gamma_engine('log_interaction');features <- e$sg_features
  e$sg_features <- function(d,type,spec=NULL) {
    f<-features(d,type,spec)
    if(is.infinite(rho))f$F<-f$F[,1:19,drop=FALSE]
    f
  }
  pen <- if(is.infinite(rho))c(0,rep(1,18)) else c(0,rep(1,18),rep(rho,45))
  e$sg_objective <- function(beta,F,w0,w1,lambda,type)
    mean(w0*exp(drop(F%*%beta))-w1*drop(F%*%beta))+lambda*sum(pen*beta^2)/2
  e$sg_gradient <- function(beta,F,w0,w1,lambda,type)
    drop(crossprod(F,w0*exp(drop(F%*%beta))-w1))/nrow(F)+lambda*pen*beta
  e
}
x5p_log_fit <- function(train,nu,t,seed,smoke=FALSE) {
  i<-which(train$T_obs==t);d<-train[i,];nu<-lapply(nu,function(v)v[i])
  constants<-if(smoke)c(.1,1) else 10^seq(-1.5,1,.5)
  multipliers<-if(smoke)c(1,Inf) else c(1,10,100,Inf)
  rows<-list();recoveries<-list()
  for(rep in 1:2) {
    folds<-make_z_stratified_folds(d,K=3,seed=seed+100000L*(rep-1))
    for(k in 1:3)for(rho in multipliers) {
      e<-x5p_log_engine(rho);tr<-which(folds!=k);va<-which(folds==k)
      for(C in constants) {
        f<-e$sg_fit_fixed(d[tr,],lapply(nu,function(v)v[tr]),t,'loglinear',C/sqrt(nrow(d)))
        g<-f$predict(d[va,])$normalized
        w0<-as.numeric(d$Z[va]==0)
        w1<-as.numeric(d$Z[va]==1)*compute_riesz_weight(lapply(nu,function(v)v[va]),t)
        rows[[length(rows)+1]]<-data.frame(repeat_id=rep,fold=k,C=C,multiplier=rho,
          lambda=C/sqrt(nrow(d)),risk=mean(w0*g-w1*log(g)))
      }
      recoveries<-c(recoveries,e$recovery)
    }
  }
  detail<-do.call(rbind,rows);cv<-aggregate(risk~C+multiplier+lambda,detail,mean)
  stopifnot(all(is.finite(cv$risk)))
  cv$selected<-seq_len(nrow(cv))==which.min(cv$risk);best<-cv[cv$selected,]
  e<-x5p_log_engine(best$multiplier)
  fit<-e$sg_fit_fixed(d,nu,t,'loglinear',best$lambda)
  fit$selected_lambda<-best$lambda;fit$selected_multiplier<-best$multiplier
  fit$cv_table<-cv;fit$cv_resample_table<-detail
  fit$numerical_recovery<-c(recoveries,e$recovery)
  fit
}

x5p_session <- function(smoke=FALSE) {
  e<-new.env(parent=globalenv());sys.source('R/learned_witness_compatibility_test.R',e)
  original<-e$learned_witness_probability_fit
  cache<-new.env(parent=emptyenv())
  e$learned_witness_probability_fit<-function(...) {
    args<-list(...);key<-as.character(args$library_seed)
    if(exists(key,cache,inherits=FALSE)) {
      old<-get(key,cache);stopifnot(identical(args,old$args));return(old$fit)
    }
    f<-do.call(original,args);assign(key,list(args=args,fit=f),cache);f
  }
  run<-function(d,method,L,seed,a,b,design) {
    stopifnot(method %in% c('rkhs','log_interaction','oracle_gamma','oracle_all'))
    records<-list()
    e$fit_gamma_rkhs<-function(train,nuisance_train,t_value,x_cols,seed,...) {
      args<-list(...)
      f<-if(method=='rkhs')fit_gamma_rkhs(train,nuisance_train,t_value,x_cols,
        lambda_grid=if(smoke)function(n_training,treatment)c(.1,1)/sqrt(n_training) else default_gamma_lambda_grid,
        bandwidth_multipliers=1/sqrt(2),cv_folds=3,cv_repeats=2,selection_rule='minimum',
        backend='nystrom',nystrom_rank=if(smoke)40L else ceiling(10*sqrt(sum(train$T_obs==t_value))),seed=seed)
      else x5p_log_fit(train,nuisance_train,t_value,seed,smoke)
      # Independent validation draws, not training or CV observations.
      val<-x5p_generate(if(smoke)1000L else 4000L,a,b,seed+700000L,design)
      vi<-val$T_obs==t_value & val$Z==0
      pred<-f$predict(val[vi,]);truth<-val[[paste0('gamma',t_value,'_oracle')]][vi]
      ri<-train$T_obs==t_value & train$Z==0
      records[[length(records)+1]]<<-list(treatment=t_value,seed=seed,
        rmse=sqrt(mean((pred$normalized-truth)^2)),mae=mean(abs(pred$normalized-truth)),
        training_normalization=mean(f$predict(train[ri,])$normalized),
        validation_clip_rate=mean(pred$raw<1e-6|pred$raw>100),
        validation_normalization=mean(pred$normalized),validation_n=sum(vi),
        selected_lambda=f$selected_lambda,selected_multiplier=f$selected_multiplier,
        bandwidth_multiplier=f$selected_bandwidth_multiplier,base_bandwidth=f$base_bandwidth,
        selected_bandwidth=f$selected_bandwidth,rank=f$nystrom_effective_rank,
        cv=f$cv_table,cv_resamples=f$cv_resample_table,recovery=f$numerical_recovery)
      f
    }
    oracle<-method %in% c('oracle_gamma','oracle_all')
    f<-e$test_learned_witness_dual_msm_compatibility(d,L,x_cols=paste0('X',1:5),
      folds=if(smoke)3L else 5L,primitive_learner=if(method=='oracle_all')'oracle' else 'super_learner',
      gamma_learner=if(oracle)'oracle' else 'rkhs_cv',
      super_learner_candidates=if(smoke)c('glm_main','glm_interactions') else sfp_contract()$primitive_library,
      super_learner_cv_folds=3,super_learner_ranger_num_trees=if(smoke)50L else 500L,
      gamma_cv_folds=3,gamma_cv_repeats=2,gamma_cv_selection_rule='minimum',
      gamma_cv_bandwidth_multipliers=1/sqrt(2),allocation_spline_df=5,
      bootstrap_replications=if(smoke)99L else 399L,seed=seed)
    stopifnot(identical(f$fold_id,make_tz_stratified_folds(d,K=if(smoke)3L else 5L,seed=seed)))
    list(summary=f$summary,gamma_diagnostics=f$gamma_diagnostics,
      gamma_audit=records,moment_detail=f$moment_detail,settings=f$settings,
      witness_diagnostics=f$witness_diagnostics,probability_diagnostics=f$probability_diagnostics,
      primitive_cache_entries=length(ls(cache)),fold_id=f$fold_id)
  }
  list(run=run,cache=cache)
}
