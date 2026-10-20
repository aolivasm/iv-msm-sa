# Frozen threshold classes and matching population objectives for the X5 pilot.
# This file adds no overrides to the public estimator.
restricted_x5_contract <- function() {
  list(version = "x5_fixed_classes_n1500_3000_6000_r100_fullgrid_v1",
       design = "continuous_x5_nonlinear", sample_sizes = c(1500L,3000L,6000L),
       seeds = 170001:170100, Lambda = c(1.5,2,2.5,3,4,5),
       epsilon = c(0,0.01), classes = c("fixed_spline5","fixed_rbf60"),
       spline_knots = c(-1.2,-0.4,0.4,1.2), boundary = c(-2,2),
       rbf_rank = 60L, rbf_bandwidth = 2, rbf_skip = 20L,
       q_penalty = "1/n_training; squared coefficients in full-rank training-whitened features",
       target = "unpenalized fixed-class pseudo-bound, NOT normalized-direct plug-in",
       population_tolerance = 5e-4, population_gradient_tolerance = 1e-4,
       diagnostic_replications = 10L)
}
restricted_x5_cells <- function() {
  z <- restricted_x5_contract()
  expand.grid(threshold_class=z$classes,Lambda=z$Lambda,epsilon=z$epsilon,
              stringsAsFactors=FALSE)
}
restricted_x5_features <- function(
    data, threshold_class, contract = restricted_x5_contract()
) {
  z <- contract; x_cols <- paste0("X",1:5)
  stopifnot(all(x_cols %in% names(data)))
  if (threshold_class == "fixed_spline5") {
    return(do.call(cbind,lapply(x_cols,function(k)
      splines::ns(data[[k]],knots=z$spline_knots,Boundary.knots=z$boundary,
                  intercept=FALSE))))
  }
  stopifnot(threshold_class == "fixed_rbf60")
  centers <- make_fixed_rbf_centers(z$rbf_rank,5,lower=-2,upper=2,
                                    skip=z$rbf_skip,marginal="truncated_normal")
  make_threshold_fixed_rbf_features(data,x_cols,centers,z$rbf_bandwidth)$features
}
restricted_x5_bind <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z)>0,x)
  if (!length(x)) return(data.frame())
  cols <- unique(unlist(lapply(x,names)))
  x <- lapply(x,function(z) {
    for (k in setdiff(cols,names(z))) z[[k]] <- NA
    z[cols]
  })
  rownames_out <- do.call(rbind,x); rownames(rownames_out)<-NULL; rownames_out
}
restricted_x5_save <- function(object,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temporary <- tempfile(pattern=".checkpoint-",tmpdir=dirname(path))
  on.exit(unlink(temporary),add=TRUE)
  saveRDS(object,temporary,compress=TRUE)
  if (!file.rename(temporary,path)) stop("Atomic checkpoint rename failed: ",path)
  invisible(path)
}
# A storage schema, not a change to the statistical contract. Retain sufficient
# summaries for the paper and q coefficients for population-regret evaluation.
# Seeds and frozen source reproduce data/folds. Observation-level EIF vectors
# and full nuisance learners are not necessary for these reporting tasks.
restricted_x5_compact_cell <- function(object) {
  required <- c("success","contract","status","error")
  if(isTRUE(object$success)) required <- c(required,"endpoints","qmodels",
    "qdiagnostics","diagnostics","eif_diagnostics","nuisance_diagnostics",
    "gamma_diagnostics","threshold_converged","plugin_summary")
  if(!is.list(object) || !all(required %in% names(object)))
    stop("Incomplete checkpoint; refusing to discard any of its contents.")
  if(!is.null(object$fold_id)) {
    fold_sizes <- tabulate(object$fold_id,nbins=5L)
    if(!is.null(object$fold_sizes) && !identical(object$fold_sizes,fold_sizes))
      stop("Inconsistent fold sizes; refusing compaction.")
    object$fold_sizes <- fold_sizes
  }
  # Keep unknown fields rather than silently discarding unrecognized results.
  object[c("signals","fold_id")] <- NULL
  object$storage_schema <- "paper_minimal_v1"
  object
}
restricted_x5_fit_q <- function(
    train, y, level, weight, threshold_class, n_training,
    contract = restricted_x5_contract()
) {
  f <- restricted_x5_features(train,threshold_class,contract)
  center <- colMeans(f); centered <- sweep(f,2,center,"-")
  s <- svd(centered,nu=0)
  if (length(s$d)!=ncol(f) || min(s$d)<max(s$d)*1e-12)
    stop("Threshold design is not full rank; class will not be silently reduced.")
  rotation <- sweep(s$v,2,s$d/sqrt(nrow(f)),"/")
  features <- centered %*% rotation
  penalty <- 1/as.integer(n_training)
  fit <- fit_weighted_quantile_ridge_admm(features,y,level,weight,penalty,
    rho=.01,maxit=12000L,abs_tol=1e-5,rel_tol=1e-4)
  if(fit$convergence!=0L)
    fit <- fit_weighted_quantile_ridge_admm(features,y,level,weight,penalty,
      rho=.01,maxit=60000L,abs_tol=1e-5,rel_tol=1e-4)
  raw <- as.vector(rotation %*% fit$coefficients[-1])
  coefficients <- c(fit$coefficients[1]-sum(center*raw),raw)
  list(coefficients=coefficients,threshold_class=threshold_class,
       threshold_contract=contract,
       convergence=fit$convergence,iterations=fit$iterations,
       penalty=penalty,primal_residual=fit$primal_residual,
       dual_residual=fit$dual_residual,objective=fit$objective,
       penalized_objective=fit$penalized_objective)
}
restricted_x5_predict_q <- function(model,newdata) {
  contract <- if (is.null(model$threshold_contract)) {
    restricted_x5_contract()
  } else model$threshold_contract
  as.vector(cbind(1,restricted_x5_features(
    newdata,model$threshold_class,contract
  )) %*%
              model$coefficients)
}
restricted_x5_provider <- function(
    threshold_class, store, contract = restricted_x5_contract()
) {
  force(threshold_class); force(store); force(contract)
  function(mode,fold,treatment,...) {
    force(mode); force(fold); force(treatment)
    list(thresholds=function(train,validation,nuisance_train,envelope_train,
                             tau_train,...) {
      use <- train$T_obs==treatment & train$Z==0
      weights <- compute_target_odds(nuisance_train,treatment)[use] *
        envelope_train$D[use]/nuisance_train$e[use]
      lo <- restricted_x5_fit_q(train[use,,drop=FALSE],train$Y[use],
        1-tau_train[use],weights,threshold_class,nrow(train),contract)
      hi <- restricted_x5_fit_q(train[use,,drop=FALSE],train$Y[use],
        tau_train[use],weights,threshold_class,nrow(train),contract)
      key <- paste(mode,fold,treatment,sep="_")
      assign(key,list(mode=mode,fold=fold,treatment=treatment,
                      lower=lo,upper=hi),envir=store)
      lower_predict <- function(newdata) restricted_x5_predict_q(lo,newdata)
      upper_predict <- function(newdata) restricted_x5_predict_q(hi,newdata)
      list(q_lower_train=lower_predict(train),q_upper_train=upper_predict(train),
           q_lower_validation=lower_predict(validation),
           q_upper_validation=upper_predict(validation),
           lower_predict=lower_predict,upper_predict=upper_predict)
    })
  }
}
# Build cumulative integrals without materializing repeated (Y,X) data.
restricted_x5_population_arm <- function(
    cell, treatment, x_nodes, y_nodes,
    specification = make_observed_data_dgp_spec("continuous_x5_nonlinear"),
    contract = restricted_x5_contract()
) {
  rule <- observed_dgp_covariate_rule(specification,x_nodes_per_dimension=x_nodes)
  x <- as.data.frame(rule$x); nx <- nrow(x)
  design <- cbind(1,restricted_x5_features(
    x,cell$threshold_class,contract
  ))
  nuisance <- observed_dgp_nuisances(x,specification)
  location <- observed_dgp_outcome_location(x,treatment,specification)
  tilt <- observed_dgp_tilt(x,treatment,specification)
  normalizer <- (pnorm(2)-pnorm(-2))/(pnorm(2-tilt)-pnorm(-2-tilt))
  residual <- seq(-2,2,length.out=y_nodes)
  mass <- dnorm(residual); mass[c(1,y_nodes)] <- mass[c(1,y_nodes)]/2
  mass <- mass/sum(mass)
  cumulative <- cumulative_y <- matrix(0,nx,y_nodes+1L)
  Lbar <- Ubar <- LY <- UY <- EY <- numeric(nx)
  min_gap <- Inf
  for (j in seq_len(y_nodes)) {
    gamma <- exp(tilt*residual[j]-tilt^2/2)*normalizer
    envelope <- observed_dgp_envelopes(gamma,nuisance$pi0,nuisance$pi1,
                          treatment,cell$Lambda,cell$Lambda,cell$epsilon)
    if (any(gamma<envelope$gamma_lower-1e-7 |
            gamma>envelope$gamma_upper+1e-7))
      stop("Population compatibility failed; no projection is permitted.")
    gap <- envelope$U-envelope$L
    min_gap <- min(min_gap,gap)
    yy <- location+residual[j]
    cumulative[,j+1L] <- cumulative[,j]+mass[j]*gap
    cumulative_y[,j+1L] <- cumulative_y[,j]+mass[j]*gap*yy
    Lbar <- Lbar+mass[j]*envelope$L; Ubar<-Ubar+mass[j]*envelope$U
    LY<-LY+mass[j]*envelope$L*yy; UY<-UY+mass[j]*envelope$U*yy
    EY<-EY+mass[j]*yy
  }
  if (min_gap<0 || any(Lbar>1+1e-7 | Ubar<1-1e-7))
    stop("Population normalization feasibility failed.")
  total <- cumulative[,y_nodes+1L]; total_y <- cumulative_y[,y_nodes+1L]
  factual <- if(treatment==1L) nuisance$pi0 else 1-nuisance$pi0
  objective_weight <- rule$weights*(1-factual)
  index <- seq_len(nx)
  output <- list()
  for (bound in c("lower","upper")) {
    level <- if(bound=="lower") (1-Lbar)/total else (Ubar-1)/total
    if(any(!is.finite(level) | level<=0 | level>=1))
      stop("Population quantile levels outside (0,1).")
    values <- function(q) {
      j <- findInterval(q-location,residual)+1L
      below<-cumulative[cbind(index,j)]; y_below<-cumulative_y[cbind(index,j)]
      loss<-(1-level)*(q*below-y_below)+
        level*(total_y-y_below-q*(total-below))
      list(loss=loss,score=below-level*total,below=below,y_below=y_below)
    }
    pointwise_index <- vapply(index,function(i)
      which(cumulative[i,-1]>=level[i]*total[i])[1L],integer(1))
    pointwise_q <- location+residual[pointwise_index]
    initial <- lm.wfit(design,pointwise_q,pmax(objective_weight*total,1e-12))$coefficients
    if(any(!is.finite(initial))) stop("Population class basis is rank deficient.")
    objective <- function(beta) sum(objective_weight*values(as.vector(design%*%beta))$loss)
    gradient <- function(beta) as.vector(crossprod(design,
                       objective_weight*values(as.vector(design%*%beta))$score))
    numerical_ridge <- 1e-10
    opt <- optim(initial,function(b) objective(b)+numerical_ridge*sum(b[-1]^2)/2,
      function(b) gradient(b)+c(0,numerical_ridge*b[-1]),method="BFGS",
      control=list(maxit=4000L,reltol=1e-12))
    endpoint <- function(q) {
      v<-values(q)
      missing <- if(bound=="lower")
        q*(1-Lbar-v$below)+LY+v$y_below else
        q*(1-Ubar+v$below)+UY-v$y_below
      sum(rule$weights*(factual*EY+(1-factual)*missing))
    }
    q <- as.vector(design%*%opt$par)
    target<-endpoint(q); sharp<-endpoint(pointwise_q)
    direct_residual <- sum(objective_weight*q*
      if(bound=="lower") (Lbar+values(q)$below-1) else (Ubar-values(q)$below-1))
    output[[bound]]<-list(endpoint=paste0("mu",treatment,"_",bound),
      target=target,unrestricted=sharp,coefficients=opt$par,
      convergence=opt$convergence,maximum_gradient=max(abs(gradient(opt$par))),
      numerical_ridge=numerical_ridge,direct_minus_pseudo=direct_residual,
      x_nodes=x_nodes,y_nodes=y_nodes)
  }
  output
}
restricted_x5_population <- function(
    cell, x_nodes, y_nodes,
    specification = make_observed_data_dgp_spec("continuous_x5_nonlinear"),
    contract = restricted_x5_contract()
) {
  arms<-unlist(lapply(0:1,function(t)
    restricted_x5_population_arm(
      cell,t,x_nodes,y_nodes,specification,contract
    )),recursive=FALSE)
  detail<-do.call(rbind,lapply(arms,function(z)
    data.frame(endpoint=z$endpoint,target=z$target,unrestricted=z$unrestricted,
      convergence=z$convergence,maximum_gradient=z$maximum_gradient,
      direct_minus_pseudo=z$direct_minus_pseudo)))
  endpoint_value<-function(name,column) detail[detail$endpoint==name,column]
  for (b in c("lower","upper")) {
    opposite<-if(b=="lower") "upper" else "lower"
    detail<-rbind(detail,data.frame(endpoint=paste0("ate_",b),
      target=endpoint_value(paste0("mu1_",b),"target")-
        endpoint_value(paste0("mu0_",opposite),"target"),
      unrestricted=endpoint_value(paste0("mu1_",b),"unrestricted")-
        endpoint_value(paste0("mu0_",opposite),"unrestricted"),
      convergence=NA_integer_,maximum_gradient=NA_real_,direct_minus_pseudo=NA_real_))
  }
  list(detail=detail,models=arms,x_nodes=x_nodes,y_nodes=y_nodes)
}
