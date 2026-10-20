# Compact study bookkeeping. Statistical routines are frozen from the main study.
source('paper_config.R')
stress_config <- function() paper_profile$config
stress_tasks <- function() {
 c<-stress_config();r<-c$repetitions;s<-nrow(c$settings)
 data.frame(setting=rep(seq_len(s),each=r),n=c$n,
 replication=rep(seq_len(r),s),id=seq_len(r*s),seed=c$seed_base+rep(seq_len(r),s))
}
stress_spec <- function() {
 c<-stress_config();s<-expand.grid(method=c$methods,Lambda=c$Lambdas,class=c$classes,stringsAsFactors=FALSE)
 s$key<-with(s,paste0(method,'_',class,'_L',Lambda));s
}
stress_compatible <- function(setting,Lambda) !(setting==4L && Lambda%in%paper_profile$incompatible_Lambda)
stress_atomic <- function(object,path,overwrite=TRUE) {
 if(!overwrite && file.exists(path))stop('Refusing overwrite: ',path)
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 tmp<-tempfile('unpublished_',dirname(path));on.exit(if(file.exists(tmp))unlink(tmp))
 withCallingHandlers(saveRDS(object,tmp,compress='gzip'),warning=function(w)stop(conditionMessage(w)))
 stopifnot(file.info(tmp)$size>0,identical(readRDS(tmp),object))
 if(!file.rename(tmp,path))stop('Cannot publish checkpoint: ',path)
 stopifnot(identical(readRDS(path),object));invisible(path)
}
stress_check_record <- function(v,t,s,signature) {
 stopifnot(is.list(v),identical(v$task,t),identical(v$signature,signature),
   identical(v$method,s$method),identical(v$class,s$class),identical(v$Lambda,s$Lambda),
   is.logical(v$success),length(v$success)==1L,!is.na(v$success),isTRUE(v$attempted),
   identical(v$compatible,stress_compatible(t$setting,s$Lambda)))
 if(!v$compatible)stopifnot(is.na(v$outward),
   if(isTRUE(v$success)) all(is.na(v$ate$target))&&all(is.na(v$ate$centered_coverage)) else TRUE)
 invisible(TRUE)
}
stress_load <- function(path,t,signature) {
 if(!file.exists(path))return(list(task=t,signature=signature,records=list(),gamma_selection=list()))
 z<-readRDS(path)
 stopifnot(identical(z$task,t),identical(z$signature,signature),is.list(z$records))
 if(length(z$records)) {
   stopifnot(!is.null(names(z$records)),!anyDuplicated(names(z$records)))
   s<-stress_spec();stopifnot(all(names(z$records)%in%s$key))
   for(k in names(z$records))stress_check_record(z$records[[k]],t,s[s$key==k,,drop=FALSE],signature)
 }
 z
}
stress_measure <- function(e,compatible,target=NULL) {
 ate<-subset(e,endpoint%in%c('ate_lower','ate_upper'))
 stopifnot(nrow(ate)==2,setequal(ate$endpoint,c('ate_lower','ate_upper')),
   all(is.finite(ate$one_step_estimate)),all(is.finite(ate$standard_error)),all(ate$standard_error>0))
 ate$target<-NA_real_
 if(compatible) {
   stopifnot(nrow(target)==1,all(is.finite(c(target$lower,target$upper))))
   ate$target<-ifelse(ate$endpoint=='ate_lower',target$lower,target$upper)
 }
 ate$error<-ate$one_step_estimate-ate$target
 ate$centered_coverage<-with(ate,confidence_lower<=target&confidence_upper>=target)
 lo<-subset(ate,endpoint=='ate_lower');hi<-subset(ate,endpoint=='ate_upper')
 list(ate=ate,outward=if(compatible)lo$confidence_lower<=lo$target&&hi$confidence_upper>=hi$target else NA,
   point_width=hi$one_step_estimate-lo$one_step_estimate,
   confidence_width=hi$confidence_upper-lo$confidence_lower,
   crossed_point=hi$one_step_estimate<lo$one_step_estimate,
   crossed_confidence=hi$confidence_upper<lo$confidence_lower)
}
