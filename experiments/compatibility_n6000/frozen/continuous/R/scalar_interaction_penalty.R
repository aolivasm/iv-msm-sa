# Isolated direct density-ratio experiment; production defaults unchanged.
source('R/scalar_log_spline.R')
sip_engine <- function(multiplier=1) {
 stopifnot(length(multiplier)==1,multiplier>=1)
 e<-sls_engine(if(is.infinite(multiplier))'log_additive' else 'log_tensor')
 if(is.infinite(multiplier))return(e)
 penalty<-c(0,rep(1,6),rep(multiplier,9))
 e$sg_objective<-function(beta,F,w0,w1,lambda,type)
  mean(w0*exp(drop(F%*%beta))-w1*drop(F%*%beta))+lambda*sum(penalty*beta^2)/2
 e$sg_gradient<-function(beta,F,w0,w1,lambda,type)
  drop(crossprod(F,w0*exp(drop(F%*%beta))-w1))/nrow(F)+lambda*penalty*beta
 e
}
sip_fit <- function(train,nuisance_train,t,seed,
 constants=10^seq(-1.5,1,.5),multipliers=c(1,10,100,Inf)) {
 i<-which(train$T_obs==t);d<-train[i,];nu<-lapply(nuisance_train,function(v)v[i])
 rows<-list()
 for(rep in 1:2){
  folds<-make_z_stratified_folds(d,K=3,seed=seed+100000L*(rep-1))
  for(k in 1:3)for(mult in multipliers){
   e<-sip_engine(mult);tr<-which(folds!=k);va<-which(folds==k)
   for(C in constants){
    lam<-C/sqrt(nrow(d))
    f<-e$sg_fit_fixed(d[tr,],lapply(nu,function(v)v[tr]),t,'loglinear',lam)
    g<-f$predict(d[va,])$normalized
    w0<-as.numeric(d$Z[va]==0)
    w1<-as.numeric(d$Z[va]==1)*compute_riesz_weight(lapply(nu,function(v)v[va]),t)
    rows[[length(rows)+1]]<-data.frame(repeat_id=rep,fold=k,C=C,multiplier=mult,
     lambda=lam,risk=mean(w0*g-w1*log(g)))
   }
  }
 }
 detail<-do.call(rbind,rows)
 cv<-aggregate(risk~C+multiplier+lambda,detail,mean)
 cv$selected<-seq_len(nrow(cv))==which.min(cv$risk);best<-cv[cv$selected,]
 e<-sip_engine(best$multiplier)
 fit<-e$sg_fit_fixed(d,nu,t,'loglinear',best$lambda)
 fit$selected_lambda<-best$lambda;fit$selected_multiplier<-best$multiplier
 fit$cv_table<-cv;fit$cv_resample_table<-detail
 fit
}
