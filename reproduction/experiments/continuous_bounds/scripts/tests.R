source('scripts/common.R')
tmp<-tempfile('stress-tests-');dir.create(tmp)
c<-stress_config();tasks<-stress_tasks();s<-stress_spec();signature<-list(test=TRUE)
stopifnot(nrow(tasks)==2000,nrow(s)==8,all(tasks$n==6000),all(table(tasks$setting)==500),
 sum(!outer(1:4,c$Lambdas,Vectorize(stress_compatible)))==1)
t<-tasks[1501,];path<-file.path(tmp,'test.rds')
z<-stress_load(path,t,signature);stopifnot(length(z$records)==0)
make_record<-function(row,success=TRUE) {
 e<-data.frame(endpoint=c('ate_lower','ate_upper'),one_step_estimate=c(.2,.8),
  standard_error=c(.1,.1),confidence_lower=c(.004,.604),confidence_upper=c(.396,.996))
 compatible<-stress_compatible(t$setting,row$Lambda)
 v<-c(list(success=success,attempted=TRUE,compatible=compatible,task=t,signature=signature,
  method=row$method,class=row$class,Lambda=row$Lambda),
  if(success)stress_measure(e,compatible,data.frame(lower=.2,upper=.8)) else list(outward=NA))
 stress_check_record(v,t,row,signature);v
}
# Empty, partial, failed terminal, fully completed, and immutable resume.
stress_atomic(z,path);stopifnot(identical(stress_load(path,t,signature),z))
z$records[s$key[1]]<-list(make_record(s[1,,drop=FALSE]));stress_atomic(z,path)
stopifnot(length(stress_load(path,t,signature)$records)==1)
for(j in 2:nrow(s))z$records[s$key[j]]<-list(make_record(s[j,,drop=FALSE],success=j!=2))
stress_atomic(z,path);before<-tools::md5sum(path)
stopifnot(identical(stress_load(path,t,signature),z),identical(tools::md5sum(path),before))
# Simulated failed save must leave the previous checkpoint untouched.
badwriter<-stress_atomic;e<-new.env(parent=environment(stress_atomic))
e$saveRDS<-function(...)warning('simulated quota failure');environment(badwriter)<-e
stopifnot(inherits(try(badwriter(list(bad=TRUE),path),silent=TRUE),'try-error'),
 identical(tools::md5sum(path),before))
# Unreadable checkpoints and mismatched identities are rejected, never overwritten.
corrupt<-file.path(tmp,'empty.rds');file.create(corrupt)
stopifnot(inherits(try(stress_load(corrupt,t,signature),silent=TRUE),'try-error'))
bad<-z;bad$task$id<-999L;stress_atomic(bad,file.path(tmp,'wrong.rds'))
stopifnot(inherits(try(stress_load(file.path(tmp,'wrong.rds'),t,signature),silent=TRUE),'try-error'))
source('R/x5_growing_gamma.R')
for(j in 1:4) {
 a<-c$settings$a[j];b<-c$settings$b[j]
 base<-x5g_generate(200,a,b,820001,FALSE);stress<-x5g_generate(200,a,b,820001,TRUE)
 stopifnot(identical(base[c('T_obs','Z',paste0('X',1:5))],stress[c('T_obs','Z',paste0('X',1:5))]),
  !identical(base$Y,stress$Y))
}
cat('PASS: 2000-task mapping; 8 fits/task; incompatible coverage NA; empty/partial/failed/complete checkpoints; immutable resume; failed-write preservation; corrupt checkpoint rejection; paired DGP inputs.\n')
