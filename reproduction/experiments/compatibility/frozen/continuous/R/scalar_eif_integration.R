# Simulation-only conditional EIF integration for the ORIGINAL scalar DGP.
# Fitted gamma and thresholds are retained. No production estimator changes.
sei_rule <- function(n) {
 j<-seq_len(n-1);v<-j/sqrt(4*j*j-1);A<-matrix(0,n,n)
 A[cbind(j,j+1)]<-v;A[cbind(j+1,j)]<-v
 e<-eigen(A,symmetric=TRUE);o<-order(e$values)
 list(x=e$values[o],w=2*e$vectors[1,o]^2)
}
sei_factory <- function(nodes=80L) {
 rule<-sei_rule(nodes)
 function(mode,fold,treatment,gamma_fit,x_cols) {
  stopifnot(mode=='rkhs',identical(x_cols,'X'))
  function(validation,nuisance_validation,envelope_validation,t_value,Lambda0,Lambda1,
   q_lower_validation,q_upper_validation,baseline_outcomes,comparison_outcomes) {
   stopifnot(t_value==treatment)
   n<-nrow(validation)
   out0<-matrix(NA_real_,n,length(baseline_outcomes),dimnames=list(NULL,baseline_outcomes))
   out1<-matrix(NA_real_,n,length(comparison_outcomes),dimnames=list(NULL,comparison_outcomes))
   for(start in seq(1,n,by=20)) {
    ii<-start:min(n,start+19);N<-length(ii);cutlist<-vector('list',N)
    center<-.65*validation$X[ii]+.75*t_value
    for(j in seq_along(ii)) {
     cuts<-sort(c(-2,2,pmax(-2,pmin(2,c(q_lower_validation[ii[j]],q_upper_validation[ii[j]])-center[j]))))
     p0<-nuisance_validation$pi0[ii[j]];p1<-nuisance_validation$pi1[ii[j]]
     if(t_value==0){p0<-1-p0;p1<-1-p1}
     a0<--p0/(1-p0);cl<-((1-p1)/Lambda1+p1)/(1-p0);cu<-((1-p1)*Lambda1+p1)/(1-p0)
     levels<-c((1/Lambda0-a0)/cl,(Lambda0-a0)/cu)
     grid<-seq(-2,2,length.out=65)
     gf<-function(rr)gamma_fit$predict(data.frame(Y=rr+center[j],X=rep(validation$X[ii[j]],length(rr))))$normalized
     gv<-gf(grid)
     for(level in levels){v<-gv-level;idx<-which(v[-length(v)]*v[-1]<0)
      for(k in idx)cuts<-c(cuts,uniroot(function(rr)gf(rr)-level,c(grid[k],grid[k+1]),tol=1e-10)$root)}
     cutlist[[j]]<-sort(unique(cuts))
    }
    K<-max(lengths(cutlist)-1)*nodes;r<-w<-matrix(0,N,K)
    for(j in seq_along(ii)){
     cuts<-cutlist[[j]]
     for(k in seq_len(length(cuts)-1)){cols<-((k-1)*nodes+1):(k*nodes);half<-(cuts[k+1]-cuts[k])/2
      r[j,cols]<-(cuts[k+1]+cuts[k])/2+half*rule$x;w[j,cols]<-half*rule$w}
    }
    eta<-(2*t_value-1)*.15
    p0<-w*dnorm(r)/siv_normalizer(0);p1<-w*dnorm(r-eta)/siv_normalizer(eta)
    stopifnot(max(abs(rowSums(p0)-1))<1e-9,max(abs(rowSums(p1)-1))<1e-9)
    yy<-as.vector(t(r+center));xx<-rep(validation$X[ii],each=K)
    nu<-lapply(nuisance_validation,function(v)rep(v[ii],each=K))
    gam<-gamma_fit$predict(data.frame(Y=yy,X=xx))$normalized
    env<-compute_projected_envelopes(nu,gam,t_value,Lambda0,Lambda1,
     gamma_compatibility_action='project',gamma_projection_margin=1e-4)
    pr<-compute_eif_row_primitives(data.frame(Y=yy),nu,env,t_value,Lambda0,Lambda1,
     rep(q_lower_validation[ii],each=K),rep(q_upper_validation[ii],each=K))
    avg<-function(v,p)rowSums(p*matrix(v,nrow=N,byrow=TRUE))
    out0[ii,'.m0']<-avg(yy,p0)
    for(ep in c('lower','upper')){
     v<-pr[[ep]]
     out0[ii,paste0('.theta_',ep)]<-avg(v$lambda*yy,p0)
     out0[ii,paste0('.D_r0_',ep)]<-avg(v$D_r0,p0)
     out0[ii,paste0('.D_r1_',ep)]<-avg(v$D_r1,p0)
     out0[ii,paste0('.H0_',ep)]<-avg(v$H0,p0)
     out1[ii,paste0('.H1_',ep)]<-avg(v$H1,p1)
    }
   }
   list(baseline_prediction=out0,comparison_prediction=out1)
  }
 }
}
