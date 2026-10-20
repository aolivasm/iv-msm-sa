# Isolated adapter: no modification of the production compatibility test.
source('R/scalar_factorial_pilot.R')
source('R/scalar_gamma_sieve_audit.R')
source('R/conservative_compatibility_test.R')
source('R/scalar_interaction_onestep.R')

sclc_session <- function() {
  engine <- new.env(parent=globalenv())
  sys.source('R/learned_witness_compatibility_test.R',envir=engine)
  probability_original <- engine$learned_witness_probability_fit
  cache <- new.env(parent=emptyenv())
  engine$learned_witness_probability_fit <- function(...) {
    args <- list(...);key<-as.character(args$library_seed)
    if(exists(key,cache,inherits=FALSE)) {
      old<-get(key,cache);stopifnot(identical(args,old$args));return(old$fit)
    }
    fit<-do.call(probability_original,args)
    assign(key,list(args=args,fit=fit),cache);fit
  }
  run <- function(data,method,lambda_grid,seed,bootstrap=399L,folds=5L) {
    stopifnot(method %in% c('spline3','log_additive','log_interaction'))
    provider<-if(method=='log_interaction') sio_provider() else sls_provider(method)
    engine$fit_gamma_rkhs <- function(train,nuisance_train,t_value,x_cols,seed,...) {
      provider(fold=seed,treatment=t_value,train=train,nuisance_train=nuisance_train,
        x_cols=x_cols,seed=seed)
    }
    fit<-engine$test_learned_witness_dual_msm_compatibility(data,lambda_grid,
      x_cols='X',folds=folds,primitive_learner='super_learner',
      gamma_learner='rkhs_cv',
      super_learner_candidates=sfp_contract()$primitive_library,
      super_learner_cv_folds=3,super_learner_ranger_num_trees=500,
      gamma_cv_selection_rule='minimum',gamma_cv_folds=3,gamma_cv_repeats=2,
      gamma_cv_bandwidth_multipliers=1/sqrt(2),gamma_cv_lambda_grid=default_gamma_lambda_grid,
      allocation_spline_df=5,bootstrap_replications=bootstrap,seed=seed)
    fit$settings$gamma_learner<-method
    fit$settings$gamma_backend<-'experimental_direct_sieve'
    fit$gamma_fits<-as.list(attr(provider,'gamma_cache'))
    fit$gamma_diagnostics$learner<-method
    fit$primitive_cache_entries<-length(ls(cache))
    fit
  }
  list(run=run,engine=engine,cache=cache)
}

