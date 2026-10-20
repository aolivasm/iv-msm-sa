# Experimental diagnostic adapter; production testing code is unchanged.
source('R/scalar_gamma_compatibility_comparison.R')
source('R/learned_witness_compatibility_test.R')

# Reconstruct exactly the training-frozen candidate r functions of the engine.
scd_spec <- function(training, gamma, probability, treatment, Lambda, tails) {
  scores <- compute_learned_witness_scores(gamma, probability$pi0,
    probability$pi1, treatment, Lambda)
  active <- training$T_obs == treatment
  lapply(scores, function(s) list(cutoffs=make_witness_cutoffs(s, active, tails),
    scale=max(pmax(s[active],0))))
}

scd_response <- function(scores, keys, spec) {
  out <- lapply(keys,function(key) {
    bits <- strsplit(key,'__',fixed=TRUE)[[1]]
    s <- scores[[bits[1]]]; z <- spec[[bits[1]]]
    if(bits[2]=='indicator') return(as.numeric(s > z$cutoffs[bits[3]]))
    if(bits[2]=='fallback') return(as.numeric(s > 0))
    value <- pmax(s,0)/z$scale
    if(bits[2]=='sqrt') value <- sqrt(value)
    pmin(value,1)
  })
  out <- do.call(cbind,out); colnames(out)<-keys; out
}

# Midpoint integration over the truncated-normal CDF. These are numerical
# oracle conditional means, not fitted regressions. Resolution is recorded.
scd_means <- function(x,treatment,Lambda,b,keys,spec,nodes=256L,a=1.1) {
  ans <- lapply(0:1,function(z) {
    eta <- z*(2*treatment-1)*b
    residual <- eta+qnorm(pnorm(-2-eta)+((seq_len(nodes)-.5)/nodes)*siv_normalizer(eta))
    xx <- rep(x,each=nodes); rr <- rep(residual,length(x))
    scores <- compute_learned_witness_scores(siv_gamma(rr,(2*treatment-1)*b),
      plogis(.5*xx-a/2),plogis(.5*xx+a/2),treatment,Lambda)
    v <- scd_response(scores,keys,spec)
    out <- vapply(seq_along(keys),function(j)
      colMeans(matrix(v[,j],nrow=nodes)),numeric(length(x)))
    out<-matrix(out,nrow=length(x),dimnames=list(NULL,keys))
    # Integrate indicator witnesses analytically: every score is affine in
    # gamma, and gamma is monotone in the residual for nonzero b.
    g0<-siv_gamma(-2,(2*treatment-1)*b)
    g1<-siv_gamma(2,(2*treatment-1)*b)
    sl<-compute_learned_witness_scores(rep(g0,length(x)),plogis(.5*x-a/2),
      plogis(.5*x+a/2),treatment,Lambda)
    su<-compute_learned_witness_scores(rep(g1,length(x)),plogis(.5*x-a/2),
      plogis(.5*x+a/2),treatment,Lambda)
    for(j in seq_along(keys)) {
      bits<-strsplit(keys[j],'__',fixed=TRUE)[[1]]
      if(!bits[2] %in% c('indicator','fallback')) next
      cut<-if(bits[2]=='fallback') 0 else spec[[bits[1]]]$cutoffs[bits[3]]
      lo<-sl[[bits[1]]]; hi<-su[[bits[1]]]
      mass<-as.numeric(lo>cut & hi>cut)
      crossing<-which((lo>cut)!=(hi>cut))
      if(length(crossing)) {
        gstar<-g0+(cut-lo[crossing])/(hi[crossing]-lo[crossing])*(g1-g0)
        root<-pmax(-2,pmin(2,siv_gamma_inverse(gstar,(2*treatment-1)*b)))
        cdf<-(pnorm(root-eta)-pnorm(-2-eta))/siv_normalizer(eta)
        mass[crossing]<-ifelse(hi[crossing]>lo[crossing],1-cdf,cdf)
      }
      out[,j]<-mass
    }
    out
  })
  list(m0=ans[[1]],m1=ans[[2]])
}

scd_session <- function(b,nodes=256L,population_nodes=401L) {
  session <- sgc_session(); engine <- session$engine
  # Give primitive learner its own cache namespace.
  probability <- engine$learned_witness_probability_fit
  oracle_probability <- function(args) list(
    training=list(e=args$training$e_true,pi0=args$training$pi0_true,pi1=args$training$pi1_true),
    validation=list(e=args$validation$e_true,pi0=args$validation$pi0_true,pi1=args$validation$pi1_true),fit=NULL)
  method_now <- ''; context <- NULL; audit <- list(); call_index <- 0L
  engine$learned_witness_probability_fit <- function(...) {
    args<-list(...)
    if(method_now=='oracle_all') {
      args$training<-args[[1]];args$validation<-args[[2]]
      return(oracle_probability(args))
    }
    probability(...)
  }
  base_fold <- engine$fit_fold_learned_witness_signals
  base_allocation <- engine$fit_witness_allocation_means
  base_select <- engine$select_learned_witness
  engine$fit_fold_learned_witness_signals <- function(...) {
    args <- list(...);call_index<<-call_index+1L
    # Engine passes these arguments positionally.
    names(args)<-names(formals(base_fold))[seq_along(args)]
    context<<-args
    if(method_now=='oracle_all') context$spec<<-scd_spec(args$training,args$gamma_training,
      args$probability_training,args$treatment,args$Lambda,args$witness_tail_fractions)
    do.call(base_fold,args)
  }
  engine$fit_witness_allocation_means <- function(training,validation,response_training,
      treatment,x_cols,spline_df) {
    if(method_now!='oracle_all') return(base_allocation(training,validation,response_training,
      treatment,x_cols,spline_df))
    keys<-colnames(response_training)
    # Check reconstruction against the engine's actual training candidates.
    s<-compute_learned_witness_scores(context$gamma_training,
      context$probability_training$pi0,context$probability_training$pi1,treatment,context$Lambda)
    stopifnot(max(abs(scd_response(s,keys,context$spec)-response_training))<1e-12)
    tr<-scd_means(training$X,treatment,context$Lambda,b,keys,context$spec,nodes)
    va<-scd_means(validation$X,treatment,context$Lambda,b,keys,context$spec,nodes)
    list(train_z0=tr$m0,train_z1=tr$m1,validation_z0=va$m0,validation_z1=va$m1)
  }
  engine$select_learned_witness <- function(...) {
    args<-list(...); selected<-do.call(base_select,args)
    if(method_now=='oracle_all') {
      x<-qnorm(pnorm(-2)+((seq_len(population_nodes)-.5)/population_nodes)*siv_normalizer(0))
      t<-context$treatment; L<-context$Lambda
      p0<-plogis(.5*x-.55);p1<-plogis(.5*x+.55)
      if(t==0) {p0<-1-p0;p1<-1-p1}
      primitive<-compute_dual_msm_moment_primitives(p0,p1,L,L)
      h<-NULL
      rows<-lapply(c(nodes,2L*nodes),function(k) {
        m<-scd_means(x,t,L,b,selected$r_name,context$spec,k)
        base<-compute_learned_family_base(args$family,primitive,m$m0[,1],m$m1[,1],L)
        # Freeze h using the primary numerical oracle resolution.
        if(k==nodes) h<<-as.numeric(base>selected$h_cutoff)
        data.frame(call=call_index,treatment=t,Lambda=L,family=args$family,
          selected_r=selected$r_name,selected_h=selected$h_name,nodes=k,
          population_moment=mean(h*base),h_fraction=mean(h),
          mean_r_z0=mean(m$m0[,1]),mean_r_z1=mean(m$m1[,1]))
      })
      audit[[length(audit)+1L]]<<-do.call(rbind,rows)
    }
    selected
  }
  run<-function(data,method,grid,seed,bootstrap=399L) {
    stopifnot(method %in% c('loglinear','riesz_spline','oracle_gamma','oracle_all'))
    method_now<<-method; audit<<-list();call_index<<-0L
    f<-session$run(data,if(method=='oracle_all') 'oracle_gamma' else method,grid,seed,bootstrap)
    f$settings$diagnostic_method<-method
    f$settings$oracle_scope<-if(method=='oracle_all')
      'True gamma, e and pi; numerical oracle allocation means; training-selected witnesses.' else
      'Estimated allocation regressions and primitives; gamma as named.'
    f$population_witness_audit<-if(length(audit)) do.call(rbind,audit) else data.frame()
    f$settings$allocation_quadrature_nodes<-if(method=='oracle_all') nodes else NA_integer_
    f
  }
  list(run=run)
}
