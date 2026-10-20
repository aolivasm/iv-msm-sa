source('scripts/common.R')
integrate_order<-function(k){
 i<-1:(k-1);J<-matrix(0,k,k);J[cbind(i,i+1)]<-i/sqrt(4*i*i-1);J<-J+t(J)
 e<-eigen(J,symmetric=TRUE);u<-(e$values+1)/2;w<-e$vectors[1,]^2
 idx<-as.matrix(expand.grid(rep(list(seq_len(k)),5)))
 x<-matrix(qnorm(pnorm(-2)+u[idx]*(pnorm(2)-pnorm(-2))),ncol=5)
 wt<-apply(matrix(w[idx],ncol=5),1,prod);stopifnot(abs(sum(wt)-1)<1e-12)
 h<-drop(x%*%c(.14,.12,.10,.08,.06));m<-drop(x%*%c(.20,-.25,.30,.15,-.20))
 shape<-1+.6*sin(x[,1]*x[,2])+.4*sin(x[,3]*x[,4]);rows<-list()
 for(j in 1:4)for(L in config$Lambda){
  g<-config$settings[j,];arms<-lapply(0:1,function(t){
   r0<-plogis(h-g$a/2);r1<-plogis(h+g$a/2);if(t==0){r0<-1-r0;r1<-1-r1}
   p0<-plogis(-.4+1.5*t+m);p1<-plogis(-.4+1.5*t+m+(2*t-1)*g$b*shape)
   lo<-function(p,r)r*p+(1-r)*pmax(p/L,1-L*(1-p))
   hi<-function(p,r)r*p+(1-r)*pmin(L*p,1-(1-p)/L)
   lower<-pmax(lo(p0,r0),lo(p1,r1));upper<-pmin(hi(p0,r0),hi(p1,r1))
   b<-binary_arm_endpoint(list(r0=r0,r1=r1,a0=p0,a1=p1),Lambda0=L,Lambda1=L,epsilon=0,compute_gradient=FALSE)
   stopifnot(all(lower<=upper),!any(b$gamma_violation|b$envelope_crossing|b$theta_crossing),
    max(abs(lower-b$lower),abs(upper-b$upper))<1e-9)
   c(lower=sum(wt*lower),upper=sum(wt*upper))
  })
  rows[[length(rows)+1L]]<-data.frame(setting=j,Lambda=L,lower=arms[[2]]['lower']-arms[[1]]['upper'],upper=arms[[2]]['upper']-arms[[1]]['lower'])
 }
 do.call(rbind,rows)
}
coarse<-integrate_order(11);fine<-integrate_order(15)
err<-max(abs(as.matrix(coarse[,3:4])-as.matrix(fine[,3:4])))
cat('Maximum population refinement difference:',err,'\n')
stopifnot(err<1e-4)
atomic(list(targets=fine,coarse=coarse,max_refinement=err,orders=c(11,15),config=config,
 caveat='Tensor Gaussian quadrature refinement; frontier support maximization is numerical, not a certified global proof.'),'targets/population.rds')
