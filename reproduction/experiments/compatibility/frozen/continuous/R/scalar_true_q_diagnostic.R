# Simulation-only true-P fixed-class optimizer; not an unrestricted quantile.
stq_provider <- function(models,store) {
 force(models);force(store)
 function(mode,fold,treatment,...) {
  key<-paste(mode,fold,treatment,sep='_')
  lo<-models[[paste0('mu',treatment,'_lower')]]
  hi<-models[[paste0('mu',treatment,'_upper')]]
  list(thresholds=function(train,validation,...) {
   lower_predict<-function(newdata)drop(siv_spline_basis(newdata$X)%*%lo$coefficients)
   upper_predict<-function(newdata)drop(siv_spline_basis(newdata$X)%*%hi$coefficients)
   lo$convergence<-0L;hi$convergence<-0L
   assign(key,list(mode=mode,fold=fold,treatment=treatment,lower=lo,upper=hi),store)
   list(q_lower_train=lower_predict(train),q_upper_train=upper_predict(train),
    q_lower_validation=lower_predict(validation),q_upper_validation=upper_predict(validation),
    lower_predict=lower_predict,upper_predict=upper_predict)
  })
 }
}
