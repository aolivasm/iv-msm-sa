# Nested direct-Riesz sieve prototype. Experimental, not production defaults.
source('R/x5_simple_gamma.R')
x5g_features<-function(d,level,spec=NULL){
 v<-as.matrix(d[,c('Y',paste0('X',1:5))]);stopifnot(level %in% 1:3)
 if(is.null(spec)){
  spec<-list(center=colMeans(v),scale=apply(v,2,sd))
  stopifnot(all(spec$scale>0));z<-sweep(sweep(v,2,spec$center,'-'),2,spec$scale,'/')
  spec$coarse<-spec$fine<-vector('list',6)
  for(j in 1:6){
   knots<-quantile(z[,j],c(.25,.5,.75),names=FALSE)
   spec$coarse[[j]]<-splines::ns(z[,j],knots=knots[2],Boundary.knots=range(z[,j]))
   spec$fine[[j]]<-splines::ns(z[,j],knots=knots,Boundary.knots=range(z[,j]))
  }
 }
 z<-sweep(sweep(v,2,spec$center,'-'),2,spec$scale,'/')
 B<-lapply(1:6,function(j)predict(if(level==1)spec$coarse[[j]] else spec$fine[[j]],z[,j]))
 tensor<-function(bs){ans<-bs[[1]];for(next_b in bs[-1])ans<-do.call(cbind,lapply(seq_len(ncol(ans)),function(k)ans[,k]*next_b));ans}
 parts<-c(list(rep(1,nrow(d))),B,lapply(combn(1:6,2,simplify=FALSE),function(js)tensor(B[js])))
 if(level==3){C<-lapply(1:6,function(j)predict(spec$coarse[[j]],z[,j]))
  parts<-c(parts,lapply(combn(1:6,3,simplify=FALSE),function(js)tensor(C[js])))}
 list(F=do.call(cbind,parts),spec=spec)
}
x5g_solve<-function(F,w0,w1,lambda){
 penalty<-c(0,rep(lambda,ncol(F)-1))
 H<-crossprod(F,F*w0)/nrow(F)+diag(penalty)
 rhs<-drop(crossprod(F,w1-w0))/nrow(F)
 beta<-drop(solve(H,rhs));gradient<-max(abs(drop(H%*%beta)-rhs))
 stopifnot(gradient<1e-8,all(is.finite(beta)))
 raw<-1+drop(F%*%beta);positive<-pmin(100,pmax(1e-6,raw));normalizer<-mean(positive[w0==1])
 list(beta=beta,normalizer=normalizer,gradient=gradient,clip_rate=mean(raw<1e-6|raw>100))
}
x5g_predict<-function(F,fit){raw<-1+drop(F%*%fit$beta);positive<-pmin(100,pmax(1e-6,raw))
 list(raw=raw,positive=positive,normalized=pmin(100,pmax(1e-6,positive/fit$normalizer)))}
x5g_fit<-function(train,nu,t,seed){
 use<-train$T_obs==t;d<-train[use,];nu<-lapply(nu,function(v)v[use]);nt<-nrow(d)
 # Predeclared sample-size expansion; 73, 265, 425 coordinates.
 maxlevel<-if(nt>=2500)3L else if(nt>=1000)2L else 1L
 levels<-seq_len(maxlevel);constants<-10^seq(-1.5,1,.5);rows<-list()
 for(rep in 1:2){
  folds<-make_z_stratified_folds(d,3,seed+100000L*(rep-1))
  for(k in 1:3){tr<-folds!=k;va<-!tr;w0<-as.numeric(d$Z[tr]==0)
   w1<-as.numeric(d$Z[tr]==1)*compute_riesz_weight(lapply(nu,function(v)v[tr]),t)
   v0<-as.numeric(d$Z[va]==0);v1<-as.numeric(d$Z[va]==1)*compute_riesz_weight(lapply(nu,function(v)v[va]),t)
   spec<-NULL
   for(level in levels){
    f<-x5g_features(d[tr,],level,spec);spec<-f$spec;fv<-x5g_features(d[va,],level,spec)$F
    for(C in constants){
     fit<-x5g_solve(f$F,w0,w1,C/sqrt(nt));g<-x5g_predict(fv,fit)$normalized
     rows[[length(rows)+1]]<-data.frame(repeat_id=rep,fold=k,level=level,dimension=ncol(f$F),C=C,lambda=C/sqrt(nt),risk=mean(v0*g^2/2-v1*g))
    }
   }
  }
 }
 detail<-do.call(rbind,rows);cv<-aggregate(risk~level+dimension+C+lambda,detail,mean)
 # Deterministic tie ordering, minimum mean validation risk.
 cv<-cv[order(cv$level,cv$C),];best<-which.min(cv$risk);cv$selected<-seq_len(nrow(cv))==best
 selected<-cv[best,];f<-x5g_features(d,selected$level)
 fit<-x5g_solve(f$F,as.numeric(d$Z==0),as.numeric(d$Z==1)*compute_riesz_weight(nu,t),selected$lambda)
 fit$predict<-local({ff<-fit;spec<-f$spec;level<-selected$level;function(newdata)x5g_predict(x5g_features(newdata,level,spec)$F,ff)})
 fit$selected_lambda<-selected$lambda;fit$selected_level<-selected$level;fit$feature_count<-selected$dimension
 fit$cv_table<-cv;fit$cv_resample_table<-detail;fit
}

x5g_generate<-function(n,a,b,seed,stress=FALSE){
 if(!stress)return(x5s_generate(n,a,b,seed))
 # Same random stream as x5s_generate: only the conditional outcome tilt changes.
 set.seed(seed);x<-matrix(qnorm(pnorm(-2)+runif(n*5)*siv_normalizer(0)),ncol=5)
 colnames(x)<-paste0('X',1:5);z<-rbinom(n,1,.5);c<-x5s_contract()
 h<-drop(x%*%c$beta);t<-rbinom(n,1,plogis(h+a*(z-.5)))
 shape<-1+.6*sin(x[,1]*x[,2])+.4*sin(x[,3]*x[,4])
 eta<-z*(2*t-1)*b*shape
 r<-eta+qnorm(pnorm(-2-eta)+runif(n)*siv_normalizer(eta));loc<-drop(x%*%c$theta)
 data.frame(Y=loc+.75*t+r,T_obs=t,Z=z,x,e_true=.5,pi0_true=plogis(h-a/2),pi1_true=plogis(h+a/2),
  gamma0_oracle=siv_gamma(loc+.75*t+r-loc,-b*shape),
  gamma1_oracle=siv_gamma(loc+.75*t+r-loc-.75,b*shape))
}
