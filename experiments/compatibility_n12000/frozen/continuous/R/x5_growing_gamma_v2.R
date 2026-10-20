# New eligibility schedule only; v1 fits and files remain unchanged.
source('R/x5_growing_gamma.R')
x5g2_maxlevel<-function(n_arm){
 if(n_arm>=1000L)3L else if(n_arm>=700L)2L else 1L
}
# Clone the audited fitter, replacing exactly its schedule assignment.
x5g2_fit<-x5g_fit
b<-as.list(body(x5g2_fit))
hit<-which(vapply(b,function(z)is.call(z)&&identical(z[[1]],as.name('<-'))&&identical(z[[2]],as.name('maxlevel')),logical(1)))
stopifnot(length(hit)==1L)
b[[hit]]<-quote(maxlevel<-x5g2_maxlevel(nt))
body(x5g2_fit)<-as.call(b)
rm(b,hit)
