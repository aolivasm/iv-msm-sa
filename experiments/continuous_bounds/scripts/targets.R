source('scripts/common.R')
source('R/x5_simple_population_targets.R')
c<-stress_config();out<-'results/population_targets.rds'
if(file.exists(out))stop('Targets already exist; refusing overwrite.')
shape<-function(x)1+.6*sin(x$X1*x$X2)+.4*sin(x$X3*x$X4)
arm<-x5p_arm;bdy<-as.list(body(arm))
hit<-which(vapply(bdy,function(z)identical(z,quote(eta <- (2*t-1)*b)),logical(1)))
stopifnot(length(hit)==1L);bdy[[hit]]<-quote(eta <- (2*t-1)*b*shape(x));body(arm)<-as.call(bdy)
slack<-function(x,a,b,L) {
 h<-drop(as.matrix(x)%*%x5s_contract()$beta);ans<-numeric()
 for(t in 0:1) {
  p0<-plogis(h-a/2);p1<-plogis(h+a/2);if(t==0){p0<-1-p0;p1<-1-p1}
  eta<-(2*t-1)*b*shape(x);d<--p0/(1-p0)
  cl<-((1-p1)/L+p1)/(1-p0);cu<-((1-p1)*L+p1)/(1-p0)
  lm<-x5p_integral(cl,d,1/L,eta,TRUE)[,1];um<-x5p_integral(cu,d,L,eta,FALSE)[,1]
  ans<-c(ans,max(lm-1),max(1-um),vapply(c(-2,2),function(r){g<-siv_gamma(r,eta)
    max(pmax(1/L,cl*g+d)-pmin(L,cu*g+d))},numeric(1)))
 }
 max(ans)
}
corner<-as.data.frame(expand.grid(rep(list(c(-2,2)),5)));names(corner)<-paste0('X',1:5)
interior<-data.frame(X1=1.8,X2=pi/3.6,X3=1.8,X4=pi/3.6,X5=1.8)
support<-rbind(x5p_rule(16384L,9821),corner,interior)
rows<-list();details<-list()
for(j in 1:4)for(L in c$Lambdas)for(cl in c$classes) {
 a<-c$settings$a[j];b<-c$settings$b[j];compatible<-stress_compatible(j,L)
 sampled_slack<-slack(support,a,b,L)
 if(!compatible) {
  witness_slack<-slack(interior,a,b,L);stopifnot(witness_slack>1e-4)
  rows[[length(rows)+1]]<-data.frame(setting=j,a=a,b=b,Lambda=L,class=cl,compatible=FALSE,
    lower=NA_real_,upper=NA_real_,refinement=NA_real_,mean_refinement=NA_real_,sampled_slack=sampled_slack)
  next
 }
 stopifnot(sampled_slack<0)
 values<-lapply(seq_along(c(4096L,8192L,16384L)),function(k) {
  n<-c(4096L,8192L,16384L)[k];x<-x5p_rule(n,c(0L,1729L,9821L)[k])
  v<-c(arm(x,a,b,L,0,cl),arm(x,a,b,L,1,cl))
  get<-function(ep)v[[which(vapply(v,function(z)z$endpoint==ep,logical(1)))]]$target
  stopifnot(max(vapply(v,`[[`,numeric(1),'gradient'))<1e-8,
    max(abs(vapply(v,`[[`,numeric(1),'direct_minus_pseudo')))<1e-6)
  list(n=n,models=v,ate=c(lower=get('mu1_lower')-get('mu0_upper'),upper=get('mu1_upper')-get('mu0_lower')))
 })
 ar<-do.call(rbind,lapply(values,`[[`,'ate'))
 means<-sapply(values,function(v)vapply(v$models,`[[`,numeric(1),'target'))
 gap<-max(apply(ar,2,function(x)diff(range(x))))
 mgap<-max(apply(means,1,function(x)diff(range(x))))
 stopifnot(gap<5e-4,mgap<5e-4,ar[3,1]<ar[3,2])
 rows[[length(rows)+1]]<-data.frame(setting=j,a=a,b=b,Lambda=L,class=cl,compatible=TRUE,
   lower=ar[3,1],upper=ar[3,2],refinement=gap,mean_refinement=mgap,sampled_slack=sampled_slack)
 details[[paste(j,L,cl)]]<-values
 cat('Validated population target:',j,L,cl,'range',gap,'\n');flush.console()
}
stress_atomic(list(targets=do.call(rbind,rows),details=details,config=c,
 witness=interior,witness_slack=slack(interior,1.1,.15,2),
 note='Compatible cells passed numerical support and integration checks, not a formal global certificate. No targets for incompatible cell.'),out,FALSE)
