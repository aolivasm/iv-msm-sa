# Experimental growing sparse tensor spline Riesz learner.
# No RKHS/classification fit. No claim of verified one-step nuisance rates.
source('R/x5_growing_gamma.R')
x5ss_budget <- function(n_arm)max(2L,as.integer(floor(log2(n_arm)/2)-1L))
x5ss_features <- function(d,level,spec=NULL){
 v<-as.matrix(d[,c('Y',paste0('X',1:5))])
 if(is.null(spec))spec<-list(center=colMeans(v),scale=apply(v,2,sd))
 stopifnot(all(is.finite(spec$scale)),all(spec$scale>0),level>=1)
 # Bounded invertible coordinate map, learned using training data only.
 u<-.5+atan(sweep(sweep(v,2,spec$center),2,spec$scale,'/'))/pi
 F<-matrix(1,nrow(d),1);cost<-0L
 for(j in seq_len(ncol(v))){
  old<-F;oldcost<-cost;parts<-list(old);costs<-list(oldcost)
  for(k in seq_len(level)){
   B<-if(k==1)matrix(u[,j],ncol=1) else {
    m<-2^(k-2);centers<-(seq_len(m)-.5)/m
    vapply(centers,function(c)pmax(0,1-2*m*abs(u[,j]-c)),numeric(nrow(d)))
   }
   eligible<-which(oldcost+k<=level)
   if(length(eligible))for(h in seq_len(ncol(B))){
    parts[[length(parts)+1]]<-old[,eligible,drop=FALSE]*B[,h]
    costs[[length(costs)+1]]<-oldcost[eligible]+k
   }
  }
  F<-do.call(cbind,parts);cost<-unlist(costs,use.names=FALSE)
 }
 list(F=F,spec=spec,cost=cost)
}
x5ss_engine <- function(){
 e<-new.env(parent=environment(x5g_fit))
 e$x5g_features<-x5ss_features
 f<-x5g_fit;b<-as.list(body(f))
 hit<-which(vapply(b,function(z)is.call(z)&&identical(z[[1]],as.name('<-'))&&
   identical(z[[2]],as.name('maxlevel')),logical(1)))
 stopifnot(length(hit)==1)
 b[[hit]]<-quote(maxlevel<-x5ss_budget(nt));body(f)<-as.call(b);environment(f)<-e
 f
}
x5ss_fit <- function(train,nu,t,seed){
 fit<-x5ss_engine()(train,nu,t,seed)
 fit$normalization_constant<-fit$normalizer
 fit$training_negative_rate<-mean(fit$predict(train[train$T_obs==t,])$raw<=0)
 fit$rkhs_backend<-'growing_sparse_tensor_spline_riesz'
 for(k in c('selected_bandwidth_multiplier','base_bandwidth','selected_bandwidth',
 'nystrom_rank_requested','nystrom_landmark_count','nystrom_effective_rank',
 'nystrom_min_retained_eigenvalue','nystrom_max_eigenvalue'))fit[[k]]<-NA_real_
 fit
}
