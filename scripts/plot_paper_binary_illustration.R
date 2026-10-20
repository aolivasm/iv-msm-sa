# Introductory population example; no estimated confidence limits or data fitting.
pi_z <- c(.30,.70); a_z <- c(.45,.60); b_z <- c(.34,.40)
ell <- function(a,L) pmax(a/L,1-L*(1-a))
upp <- function(a,L) pmin(L*a,1-(1-a)/L)
means <- function(L) list(
  lo1=pi_z*a_z+(1-pi_z)*ell(a_z,L),
  hi1=pi_z*a_z+(1-pi_z)*upp(a_z,L),
  lo0=(1-pi_z)*b_z+pi_z*ell(b_z,L),
  hi0=(1-pi_z)*b_z+pi_z*upp(b_z,L))
gap <- function(L0,L1) {
  m<-means(c(L0,L1)); min(min(m$hi1)-max(m$lo1),min(m$hi0)-max(m$lo0))
}
curves <- function(L) {
  m<-means(L); p<-mean(pi_z)
  a<-mean(pi_z*a_z)/p; b<-mean((1-pi_z)*b_z)/(1-p)
  c(iv_lower=max(m$lo1)-min(m$hi0),iv_upper=min(m$hi1)-max(m$lo0),
    msm_z_lower=mean(m$lo1)-mean(m$hi0),msm_z_upper=mean(m$hi1)-mean(m$lo0),
    msm_lower=p*a+(1-p)*ell(a,L)-(1-p)*b-p*upp(b,L),
    msm_upper=p*a+(1-p)*upp(a,L)-(1-p)*b-p*ell(b,L))
}
frontier <- uniroot(function(L)gap(L,L),c(1,2),tol=1e-12)$root
tipping <- uniroot(function(L)curves(L)['iv_lower'],c(frontier,3.5),tol=1e-12)$root
iv_only <- c(max(pi_z*a_z)-min((1-pi_z)*b_z+pi_z),
             min(pi_z*a_z+1-pi_z)-max((1-pi_z)*b_z))
stopifnot(abs(frontier-1.375)<1e-9,abs(tipping-189/59)<1e-9,
          max(abs(iv_only-c(-.118,.482)))<1e-12,
          max(abs(curves(1.6)[1:2]-c(.1313,.292625)))<1e-10)
xstart <- uniroot(function(x)gap(x,3.5),c(1,3.5),tol=1e-12)$root
xx <- seq(xstart,3.5,length.out=1801)
boundary <- vapply(xx,function(x) {
  if(gap(x,1)>=0)1 else if(abs(gap(x,3.5))<1e-10)3.5 else
    uniroot(function(y)gap(x,y),c(1,3.5),tol=1e-12)$root
},numeric(1))
L <- sort(unique(c(seq(1,3.5,length.out=2001),frontier,1.6,tipping)))
values <- t(vapply(L,curves,numeric(6)))
values[L<frontier,1:2] <- NA_real_
blue<-'#0072B2';orange<-'#D55E00';purple<-'#CC79A7'
draw <- function() {
  par(mfrow=c(1,2),mar=c(4.3,4.5,1.7,.7),mgp=c(2.7,.8,0),
      las=1,bty='o',cex=.95)
  plot(NA,xlim=c(1,3.5),ylim=c(1,3.5),xaxs='i',yaxs='i',
       xlab=expression(Lambda[0]),ylab=expression(Lambda[1]),axes=FALSE)
  rect(1,1,3.5,3.5,col='#FAE9ED',border=NA)
  polygon(c(xx,3.5,xstart),c(boundary,3.5,3.5),col='#DDF0FA',border=NA)
  lines(xx,boundary,col='gray40',lwd=1.4)
  axis(1,seq(1,3.5,.5));axis(2,seq(1,3.5,.5));box()
  abline(0,1,lty=3,col='gray30')
  points(frontier,frontier,pch=21,bg='white',cex=.85)
  text(frontier+.06,frontier+.04,'1.375',adj=0,cex=.8)
  points(1.6,1.6,pch=19,col=blue,cex=.9)
  text(1.66,1.60,'1.6',adj=0,col=blue,cex=.85)
  text(2.7,3.15,'compatible',col=blue,cex=.9)
  text(1.32,1.13,'incompatible',col='#9B3448',cex=.8)
  text(2.5,2.36,expression(Lambda[0]==Lambda[1]),srt=38,cex=.85)
  mtext('A',side=3,line=.45,adj=0,font=2)
  plot(NA,xlim=c(1,3.5),ylim=c(-.31,.54),xaxs='i',yaxs='i',
       xlab=expression(Lambda),ylab='ATE bound',axes=FALSE)
  rect(1,iv_only[1],3.5,iv_only[2],col='#F3F3F3',border=NA)
  rect(1,-.31,frontier,.54,col='#FAE9ED',border=NA)
  abline(h=iv_only,col='gray40',lty=4,lwd=1.2)
  abline(h=0,col='gray30',lty=3)
  abline(v=frontier,col='#9B3448',lty=3)
  for(k in 5:6)lines(L,values[,k],col=orange,lty=2,lwd=1.6)
  for(k in 3:4)lines(L,values[,k],col=purple,lty=5,lwd=1.6)
  for(k in 1:2)lines(L,values[,k],col=blue,lwd=2)
  axis(1,seq(1,3.5,.5));axis(2,seq(-.3,.5,.1),labels=sprintf('%.1f',seq(-.3,.5,.1)));box()
  v<-curves(1.6)
  points(rep(1.6,6),v,col=c(blue,blue,purple,purple,orange,orange),
         pch=c(19,19,2,2,1,1),cex=.75)
  text(1.64,.07,expression(Lambda==1.6),adj=0,cex=.85)
  segments(1.6,.085,1.6,.115,lwd=.7)
  points(tipping,0,pch=21,bg='white',col=blue,lwd=1.5,cex=.95)
  arrows(2.95,-.050,tipping-.015,-.007,length=.065,angle=22,col=blue,lwd=1)
  text(2.98,-.067,bquote('IV-MSM tipping point '*Lambda == .(round(tipping,3))),
       adj=1,col=blue,cex=.85)
  text(1.20,-.235,'incompatible',col='#9B3448',cex=.8)
  text(frontier+.04,-.283,expression(Lambda[min]==1.375),adj=0,col='#9B3448',cex=.85)
  legend(2.32,.27,c('IV-MSM','MSM, ignore Z','MSM, condition on Z','IV only'),
         col=c(blue,orange,purple,'gray40'),lty=c(1,2,5,4),lwd=c(2,1.6,1.6,1.2),
         bty='n',cex=.8)
  mtext('B',side=3,line=.45,adj=0,font=2)
}
path <- 'output/pdf/paper_figures_20261009/main_binary_illustration.pdf'
dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
pdf(path,width=11,height=4.9,family='Times',useDingbats=FALSE);draw();dev.off()
cat('Compatibility frontier:',frontier,'\nPopulation tipping point:',tipping,'\n')
