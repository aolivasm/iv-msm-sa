# Private threshold adapter: preserve class, penalty and tolerances; retry rho.
x5s_threshold_provider<-function(threshold_class,store){
 e<-new.env(parent=globalenv());sys.source('R/restricted_x5_study.R',e)
 e$restricted_x5_fit_q<-function(train,y,level,weight,threshold_class,n_training,
                                contract=e$restricted_x5_contract()){
  f<-e$restricted_x5_features(train,threshold_class,contract)
  center<-colMeans(f);centered<-sweep(f,2,center,'-');s<-svd(centered,nu=0)
  if(length(s$d)!=ncol(f)||min(s$d)<max(s$d)*1e-12)stop('Threshold basis rank deficient.')
  rotation<-sweep(s$v,2,s$d/sqrt(nrow(f)),'/');features<-centered%*%rotation
  penalty<-1/as.integer(n_training);attempts<-list()
  for(rho in c(.01,.1,1,.001)){
   fit<-fit_weighted_quantile_ridge_admm(features,y,level,weight,penalty,
     rho=rho,maxit=60000L,abs_tol=1e-5,rel_tol=1e-4)
   attempts[[length(attempts)+1]]<-list(rho=rho,convergence=fit$convergence,
     iterations=fit$iterations,primal_residual=fit$primal_residual,dual_residual=fit$dual_residual)
   if(fit$convergence==0L)break
  }
  if(fit$convergence!=0L)stop('Threshold optimization failed certification after rho retries.')
  raw<-drop(rotation%*%fit$coefficients[-1])
  list(coefficients=c(fit$coefficients[1]-sum(center*raw),raw),
    threshold_class=threshold_class,threshold_contract=contract,
    convergence=fit$convergence,iterations=fit$iterations,penalty=penalty,
    primal_residual=fit$primal_residual,dual_residual=fit$dual_residual,
    objective=fit$objective,penalized_objective=fit$penalized_objective,attempts=attempts)
 }
 e$restricted_x5_provider(threshold_class,store)
}
