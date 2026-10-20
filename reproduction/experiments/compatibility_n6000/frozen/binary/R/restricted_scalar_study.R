# Scalar-only fixed-class study. Does not modify or source the running X5 driver.
restricted_scalar_contract <- function() {
  list(version = "scalar_fixed_classes_n750_1500_3000_6000_fullgrid_v1",
       design = "continuous_scalar_main", sample_sizes = c(750L,1500L,3000L,6000L),
       seeds = 170001:170500, Lambda = c(1.5,2,2.5,3,4,5),
       epsilon = c(0,0.01), classes = c("fixed_spline5","fixed_rbf8"),
       spline_knots = c(-1.2,-0.4,0.4,1.2), boundary = c(-2,2),
       rbf_rank = 8L, rbf_bandwidth = 1, rbf_skip = 20L,
       q_penalty = "1/n_training; squared coefficients in full-rank training-whitened features",
       q_solver = "smooth continuation with exact check-loss primal-dual gap <= 1e-7*(1+abs(primal))",
       target = "unpenalized fixed-class pseudo-bound, NOT normalized-direct plug-in",
       population_tolerance = 5e-4, population_gradient_tolerance = 1e-4,
       diagnostic_replications = 10L, folds = 5L,
       primitive_library = c("glm_main","glm_interactions","glmnet_interactions",
                             "ranger","xgboost"),
       gamma_penalty_constants = 10^seq(-1.5,1,by=.5),
       gamma_bandwidth_multipliers = c(.25,.5,1,2,4),
       gamma_cv_folds = 3L, gamma_selection_rule = "minimum",
       population_resolutions = list(c(80L,801L),c(160L,1601L),c(320L,3201L)))
}
restricted_scalar_tasks <- function(replications = 100L) {
  z <- restricted_scalar_contract()
  if (length(replications)!=1L || !is.finite(replications) ||
      replications!=as.integer(replications) || replications<1L ||
      replications>length(z$seeds)) stop("Replications must be an integer from 1 to 500.")
  expand.grid(sample_size=z$sample_sizes, seed=head(z$seeds,replications),
              KEEP.OUT.ATTRS=FALSE)
}
restricted_scalar_cells <- function() {
  z <- restricted_scalar_contract()
  expand.grid(threshold_class=z$classes,Lambda=z$Lambda,epsilon=z$epsilon,
              stringsAsFactors=FALSE)
}
restricted_scalar_features <- function(data, threshold_class) {
  z <- restricted_scalar_contract(); x_cols <- "X"
  stopifnot(all(x_cols %in% names(data)))
  if (threshold_class == "fixed_spline5") {
    return(do.call(cbind,lapply(x_cols,function(k)
      splines::ns(data[[k]],knots=z$spline_knots,Boundary.knots=z$boundary,
                  intercept=FALSE))))
  }
  stopifnot(threshold_class == "fixed_rbf8")
  centers <- make_fixed_rbf_centers(z$rbf_rank,1,lower=-2,upper=2,
                                    skip=z$rbf_skip,marginal="truncated_normal")
  make_threshold_fixed_rbf_features(data,x_cols,centers,z$rbf_bandwidth)$features
}
restricted_scalar_bind <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z)>0,x)
  if (!length(x)) return(data.frame())
  cols <- unique(unlist(lapply(x,names)))
  x <- lapply(x,function(z) {
    for (k in setdiff(cols,names(z))) z[[k]] <- NA
    z[cols]
  })
  rownames_out <- do.call(rbind,x); rownames(rownames_out)<-NULL; rownames_out
}
restricted_scalar_save <- function(object,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temporary <- tempfile(pattern=".checkpoint-",tmpdir=dirname(path))
  saveRDS(object,temporary,compress=TRUE)
  if (!file.rename(temporary,path)) stop("Atomic checkpoint rename failed: ",path)
  invisible(path)
}
restricted_scalar_fit_q <- function(train,y,level,weight,threshold_class,n_training) {
  f <- restricted_scalar_features(train,threshold_class)
  center <- colMeans(f); centered <- sweep(f,2,center,"-")
  s <- svd(centered,nu=0)
  if (length(s$d)!=ncol(f) || min(s$d)<max(s$d)*1e-12)
    stop("Threshold design is not full rank; class will not be silently reduced.")
  rotation <- sweep(s$v,2,s$d/sqrt(nrow(f)),"/")
  features <- centered %*% rotation
  penalty <- 1/as.integer(n_training)
  fit <- restricted_scalar_solve_quantile(features,y,level,weight,penalty)
  raw <- as.vector(rotation %*% fit$coefficients[-1])
  coefficients <- c(fit$coefficients[1]-sum(center*raw),raw)
  list(coefficients=coefficients,threshold_class=threshold_class,
       convergence=fit$convergence,iterations=fit$iterations,
       penalty=penalty,primal_residual=fit$primal_residual,
       dual_residual=fit$dual_residual,objective=fit$objective,
       penalized_objective=fit$penalized_objective,solver=fit$solver,
       duality_gap=fit$duality_gap,certificate_tolerance=fit$certificate_tolerance,
       stationarity_residual=fit$stationarity_residual,solver_exit_code=fit$solver_exit_code)
}

# Solve the SAME nonsmooth penalized check-loss problem as ADMM. Smooth losses
# are used only to generate candidate coefficients. Success is determined by a
# primal-dual certificate for the ORIGINAL, unsmoothed empirical objective.
# For D=(1,F), dual feasibility is sum(u)=0 and w*(tau-1)<=u<=w*tau;
# the dual objective is y'u-||F'u||^2/(2*lambda). Weak duality therefore bounds
# the remaining numerical optimization error, independently of optim's code.
restricted_scalar_solve_quantile <- function(features,y,level,weight,lambda) {
  features <- as.matrix(features); design <- cbind(1,features)
  n <- length(y); level <- rep_len(level,n)
  stopifnot(nrow(design)==n,length(weight)==n,all(is.finite(design)),
    all(is.finite(y)),all(is.finite(weight)),all(weight>=0),sum(weight)>0,
    all(level>0 & level<1),length(lambda)==1L,is.finite(lambda),lambda>0)
  w <- weight/sum(weight); lower <- w*(level-1); upper <- w*level
  beta <- lm.wfit(design,y,pmax(w,1e-12))$coefficients
  if(any(!is.finite(beta))) stop("Rank-deficient quantile initialization.")
  scale <- max(1,sd(y)); evaluations <- 0L; convergence <- 1L
  for(delta in scale*c(1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9)) {
    objective <- function(b) {
      r <- y-as.vector(design%*%b)
      sum(w*((level-.5)*r+.5*sqrt(r^2+delta^2)))+lambda*sum(b[-1L]^2)/2
    }
    gradient <- function(b) {
      r <- y-as.vector(design%*%b)
      u <- w*(level-.5+.5*r/sqrt(r^2+delta^2))
      -as.vector(crossprod(design,u))+c(0,lambda*b[-1L])
    }
    opt <- optim(beta,objective,gradient,method="BFGS",
                 control=list(maxit=2000L,reltol=1e-14))
    beta <- opt$par; evaluations <- evaluations+unname(opt$counts["gradient"])
    fitted <- as.vector(design%*%beta); r <- y-fitted
    u <- w*(level-.5+.5*r/sqrt(r^2+delta^2))
    # Euclidean projection onto the box plus intercept equality.
    left <- min(u-upper); right <- max(u-lower)
    for(j in 1:80) {
      shift <- (left+right)/2
      candidate <- pmax(lower,pmin(upper,u-shift))
      if(sum(candidate)>0) left <- shift else right <- shift
    }
    dual <- pmax(lower,pmin(upper,u-(left+right)/2))
    unpenalized <- sum(w*check_loss(y,fitted,level))
    primal <- unpenalized+lambda*sum(beta[-1L]^2)/2
    dual_objective <- sum(y*dual)-sum(crossprod(features,dual)^2)/(2*lambda)
    # Near a kink, tiny residual rounding can spoil the smooth derivative even
    # when the primal value is accurate. Recover the allowed subgradients at
    # near-zero residuals by bounded least squares, then enforce the intercept
    # equality exactly. Every candidate remains feasible for the original dual.
    for(active_tolerance in scale*c(1e-7,1e-6,1e-5,1e-4)) {
      active <- which(abs(r)<=active_tolerance & w>0)
      if(!length(active)) next
      candidate <- ifelse(r>=0,upper,lower)
      target_sum <- -sum(candidate[-active])
      if(target_sum<sum(lower[active])-1e-12 || target_sum>sum(upper[active])+1e-12) next
      a <- t(design[active,,drop=FALSE])
      rhs <- c(0,lambda*beta[-1L])-as.vector(crossprod(design[-active,,drop=FALSE],candidate[-active]))
      decomp <- svd(a); keep <- decomp$d>max(decomp$d)*1e-12
      start <- as.vector(decomp$v[,keep,drop=FALSE]%*%
        (as.vector(crossprod(decomp$u[,keep,drop=FALSE],rhs))/decomp$d[keep]))
      start <- pmax(lower[active],pmin(upper[active],start))
      fit_dual <- optim(start,function(u) sum((a%*%u-rhs)^2),
        function(u) as.vector(2*crossprod(a,a%*%u-rhs)),method="L-BFGS-B",
        lower=lower[active],upper=upper[active],
        control=list(maxit=1000L,factr=1,pgtol=1e-12))
      ua <- fit_dual$par
      left <- min(ua-upper[active]); right <- max(ua-lower[active])
      for(j in 1:80) {
        shift <- (left+right)/2
        if(sum(pmax(lower[active],pmin(upper[active],ua-shift)))>target_sum)
          left <- shift else right <- shift
      }
      candidate[active] <- pmax(lower[active],pmin(upper[active],ua-(left+right)/2))
      candidate_objective <- sum(y*candidate)-sum(crossprod(features,candidate)^2)/(2*lambda)
      if(abs(sum(candidate))<1e-10 && candidate_objective>dual_objective) {
        dual <- candidate; dual_objective <- candidate_objective
      }
    }
    gap <- primal-dual_objective
    feasibility <- max(abs(sum(dual)),max(lower-dual),max(dual-upper))
    stationarity <- sqrt(sum((-as.vector(crossprod(design,dual))+
                               c(0,lambda*beta[-1L]))^2))
    tolerance <- 1e-7*(1+abs(primal))
    if(is.finite(gap) && gap>=-1e-10 && gap<=tolerance && feasibility<1e-10) {
      convergence <- 0L; break
    }
  }
  list(coefficients=beta,fitted=fitted,convergence=convergence,iterations=evaluations,
    primal_residual=NA_real_,dual_residual=NA_real_,objective=unpenalized,
    penalized_objective=primal,duality_gap=gap,certificate_tolerance=tolerance,
    stationarity_residual=stationarity,solver_exit_code=opt$convergence,
    solver="certified_smooth_continuation",dual_variable=dual,
    dual_objective=dual_objective)
}
restricted_scalar_predict_q <- function(model,newdata) {
  as.vector(cbind(1,restricted_scalar_features(newdata,model$threshold_class)) %*%
              model$coefficients)
}
restricted_scalar_provider <- function(threshold_class,store) {
  force(threshold_class); force(store)
  function(mode,fold,treatment,...) {
    force(mode); force(fold); force(treatment)
    list(thresholds=function(train,validation,nuisance_train,envelope_train,
                             tau_train,...) {
      use <- train$T_obs==treatment & train$Z==0
      weights <- compute_target_odds(nuisance_train,treatment)[use] *
        envelope_train$D[use]/nuisance_train$e[use]
      lo <- restricted_scalar_fit_q(train[use,,drop=FALSE],train$Y[use],
        1-tau_train[use],weights,threshold_class,nrow(train))
      hi <- restricted_scalar_fit_q(train[use,,drop=FALSE],train$Y[use],
        tau_train[use],weights,threshold_class,nrow(train))
      key <- paste(mode,fold,treatment,sep="_")
      assign(key,list(mode=mode,fold=fold,treatment=treatment,
                      lower=lo,upper=hi),envir=store)
      lower_predict <- function(newdata) restricted_scalar_predict_q(lo,newdata)
      upper_predict <- function(newdata) restricted_scalar_predict_q(hi,newdata)
      list(q_lower_train=lower_predict(train),q_upper_train=upper_predict(train),
           q_lower_validation=lower_predict(validation),
           q_upper_validation=upper_predict(validation),
           lower_predict=lower_predict,upper_predict=upper_predict)
    })
  }
}

# Piecewise-linear interpolation of the envelope-gap density. Its cumulative
# mass and first moment are integrated analytically within each residual cell,
# so the objective and its displayed gradient are derivatives of the SAME
# numerical population functional. No fitted nuisances or empirical data enter.
restricted_scalar_population_arm <- function(cell,treatment,x_nodes,y_nodes) {
  specification <- make_observed_data_dgp_spec("continuous_scalar_main")
  rule <- observed_dgp_covariate_rule(specification,x_nodes=x_nodes)
  x <- as.data.frame(rule$x); nx <- nrow(x)
  raw_features <- restricted_scalar_features(x,cell$threshold_class)
  center <- colMeans(raw_features)
  centered <- sweep(raw_features,2,center,"-")
  s <- svd(centered,nu=0)
  if(min(s$d)<max(s$d)*1e-12) stop("Population class basis is rank deficient.")
  rotation <- sweep(s$v,2,s$d/sqrt(nx),"/")
  design <- cbind(1,centered %*% rotation)
  nuisance <- observed_dgp_nuisances(x,specification)
  location <- observed_dgp_outcome_location(x,treatment,specification)
  tilt <- observed_dgp_tilt(x,treatment,specification)
  normalizer <- (pnorm(2)-pnorm(-2))/(pnorm(2-tilt)-pnorm(-2-tilt))
  residual <- seq(-2,2,length.out=y_nodes); h <- residual[2]-residual[1]
  trap <- rep(h,y_nodes); trap[c(1,y_nodes)] <- h/2
  baseline_density <- dnorm(residual)/sum(trap*dnorm(residual))
  density <- cumulative <- cumulative_y <- matrix(0,nx,y_nodes)
  Lbar <- Ubar <- LY <- UY <- EY <- numeric(nx)
  for(j in seq_len(y_nodes)) {
    gamma <- exp(tilt*residual[j]-tilt^2/2)*normalizer
    envelope <- observed_dgp_envelopes(gamma,nuisance$pi0,nuisance$pi1,
      treatment,cell$Lambda,cell$Lambda,cell$epsilon)
    if(any(gamma<envelope$gamma_lower-1e-7 | gamma>envelope$gamma_upper+1e-7) ||
       any(envelope$U<envelope$L))
      stop("Population compatibility failed; no projection is permitted.")
    density[,j] <- (envelope$U-envelope$L)*baseline_density[j]
    yy <- location+residual[j]; mass <- trap[j]*baseline_density[j]
    Lbar <- Lbar+mass*envelope$L; Ubar <- Ubar+mass*envelope$U
    LY <- LY+mass*envelope$L*yy; UY <- UY+mass*envelope$U*yy
    EY <- EY+mass*yy
    if(j>1L) {
      d0 <- density[,j-1L]; slope <- (density[,j]-d0)/h
      dm <- h*(d0+density[,j])/2
      cumulative[,j] <- cumulative[,j-1L]+dm
      cumulative_y[,j] <- cumulative_y[,j-1L]+
        (location+residual[j-1L])*dm+d0*h^2/2+slope*h^3/3
    }
  }
  if(any(Lbar>1+1e-7 | Ubar<1-1e-7)) stop("Population normalization infeasible.")
  total <- cumulative[,y_nodes]; total_y <- cumulative_y[,y_nodes]
  # Match the first moment used by the continuous interpolated gap exactly.
  UY <- LY+total_y
  factual <- if(treatment==1L) nuisance$pi0 else 1-nuisance$pi0
  objective_weight <- rule$weights*(1-factual)
  index <- seq_len(nx)
  values <- function(q) {
    qr <- pmin(2,pmax(-2,q-location))
    j <- pmin(y_nodes-1L,pmax(1L,findInterval(qr,residual)))
    u <- qr-residual[j]; d0 <- density[cbind(index,j)]
    slope <- (density[cbind(index,j+1L)]-d0)/h
    dm <- d0*u+slope*u^2/2
    below <- cumulative[cbind(index,j)]+dm
    y_below <- cumulative_y[cbind(index,j)]+(location+residual[j])*dm+
      d0*u^2/2+slope*u^3/3
    list(below=below,y_below=y_below)
  }
  output <- list()
  for(bound in c("lower","upper")) {
    level <- if(bound=="lower") (1-Lbar)/total else (Ubar-1)/total
    if(any(!is.finite(level) | level<=0 | level>=1)) stop("Invalid population levels.")
    pointwise_q <- location+vapply(index,function(i) {
      target <- level[i]*total[i]
      j <- min(y_nodes-1L,max(1L,findInterval(target,cumulative[i,])))
      deficit <- target-cumulative[i,j]; d0 <- density[i,j]
      slope <- (density[i,j+1L]-d0)/h
      denominator <- d0+sqrt(max(0,d0^2+2*slope*deficit))
      u <- if(deficit<=0) 0 else if(denominator>0) 2*deficit/denominator else h
      residual[j]+min(h,max(0,u))
    },numeric(1))
    objective <- function(beta) {
      q <- as.vector(design%*%beta); v <- values(q)
      loss <- (1-level)*(q*v$below-v$y_below)+
        level*(total_y-v$y_below-q*(total-v$below))
      sum(objective_weight*loss)
    }
    gradient <- function(beta) as.vector(crossprod(design,
      objective_weight*(values(as.vector(design%*%beta))$below-level*total)))
    initial <- lm.wfit(design,pointwise_q,pmax(objective_weight*total,1e-12))$coefficients
    opt <- optim(initial,objective,gradient,method="BFGS",
                 control=list(maxit=4000L,reltol=1e-13))
    q <- as.vector(design%*%opt$par)
    endpoint <- function(q) {
      v <- values(q)
      missing <- if(bound=="lower") q*(1-Lbar-v$below)+LY+v$y_below else
        q*(1-Ubar+v$below)+UY-v$y_below
      sum(rule$weights*(factual*EY+(1-factual)*missing))
    }
    raw <- as.vector(rotation%*%opt$par[-1L])
    coefficients <- c(opt$par[1L]-sum(center*raw),raw)
    residual_value <- sum(objective_weight*q*if(bound=="lower")
      (Lbar+values(q)$below-1) else (Ubar-values(q)$below-1))
    output[[bound]] <- list(endpoint=paste0("mu",treatment,"_",bound),
      target=endpoint(q),unrestricted=endpoint(pointwise_q),coefficients=coefficients,
      convergence=opt$convergence,maximum_gradient=max(abs(gradient(opt$par))),
      numerical_ridge=0,direct_minus_pseudo=residual_value,
      x_nodes=x_nodes,y_nodes=y_nodes)
  }
  output
}
restricted_scalar_population <- function(cell,x_nodes,y_nodes) {
  arms <- unlist(lapply(0:1,function(t)
    restricted_scalar_population_arm(cell,t,x_nodes,y_nodes)),recursive=FALSE)
  detail <- do.call(rbind,lapply(arms,function(z)
    data.frame(endpoint=z$endpoint,target=z$target,unrestricted=z$unrestricted,
      convergence=z$convergence,maximum_gradient=z$maximum_gradient,
      direct_minus_pseudo=z$direct_minus_pseudo)))
  get_value <- function(name,column) detail[detail$endpoint==name,column]
  for(b in c("lower","upper")) {
    opposite <- if(b=="lower") "upper" else "lower"
    detail <- rbind(detail,data.frame(endpoint=paste0("ate_",b),
      target=get_value(paste0("mu1_",b),"target")-get_value(paste0("mu0_",opposite),"target"),
      unrestricted=get_value(paste0("mu1_",b),"unrestricted")-
        get_value(paste0("mu0_",opposite),"unrestricted"),
      convergence=NA_integer_,maximum_gradient=NA_real_,
      direct_minus_pseudo=get_value(paste0("mu1_",b),"direct_minus_pseudo")-
        get_value(paste0("mu0_",opposite),"direct_minus_pseudo")))
  }
  list(detail=detail,models=arms,x_nodes=x_nodes,y_nodes=y_nodes)
}
