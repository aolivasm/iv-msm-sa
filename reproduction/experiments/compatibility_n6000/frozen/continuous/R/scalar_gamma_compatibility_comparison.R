# Isolated adapter: no modification of the production compatibility test.
source('R/scalar_factorial_pilot.R')
source('R/scalar_gamma_sieve_audit.R')
source('R/conservative_compatibility_test.R')

sgc_session <- function() {
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
    stopifnot(method %in% c('rkhs','loglinear','riesz_spline','oracle_gamma'))
    engine$fit_gamma_rkhs <- function(train,nuisance_train,t_value,x_cols,seed,...) {
      stopifnot(identical(x_cols,'X'))
      if(method=='rkhs') return(fit_gamma_rkhs(train,nuisance_train,t_value,
        x_cols=x_cols,backend='exact',bandwidth_multipliers=1/sqrt(2),
        lambda_grid=default_gamma_lambda_grid,cv_folds=3,cv_repeats=2,
        selection_rule='minimum',seed=seed))
      sg_fit(train,nuisance_train,t_value,method,seed)
    }
    fit<-engine$test_learned_witness_dual_msm_compatibility(data,lambda_grid,
      x_cols='X',folds=folds,primitive_learner='super_learner',
      gamma_learner=if(method=='oracle_gamma') 'oracle' else 'rkhs_cv',
      super_learner_candidates=sfp_contract()$primitive_library,
      super_learner_cv_folds=3,super_learner_ranger_num_trees=500,
      gamma_cv_selection_rule='minimum',gamma_cv_folds=3,gamma_cv_repeats=2,
      gamma_cv_bandwidth_multipliers=1/sqrt(2),gamma_cv_lambda_grid=default_gamma_lambda_grid,
      allocation_spline_df=5,bootstrap_replications=bootstrap,seed=seed)
    fit$settings$gamma_learner<-method
    fit$settings$gamma_backend<-if(method=='rkhs') 'exact' else 'direct_or_oracle'
    fit$settings$oracle_scope<-'oracle_gamma retains estimated primitives and allocation regressions'
    fit$gamma_diagnostics$learner<-method
    fit$primitive_cache_entries<-length(ls(cache))
    fit
  }
  list(run=run,engine=engine,cache=cache)
}
