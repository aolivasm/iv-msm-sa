# Standalone paired stress comparator bookkeeping; no baseline writes.
source('R/x5_simple_study.R')
# Evaluate ONLY the exact generator from the frozen LPC source, not its learners.
ast<-parse('R/x5_growing_gamma_frozen.R')
definition<-Filter(function(v)is.call(v)&&identical(v[[1]],as.name('<-'))&&
 identical(v[[2]],as.name('x5g_generate')),as.list(ast))
stopifnot(length(definition)==1L);eval(definition[[1]])
config<-list(n=6000L,repetitions=500L,settings=data.frame(a=c(.55,1.1,1.65,1.1),b=c(.075,.075,.075,.15)),
 Lambdas=c(2,3),methods=c('tan','dorn'),adjustments=c('x','xz'),seed_base=820000L,stress=TRUE)
tasks<-data.frame(setting=rep(1:4,each=500),n=6000L,replication=rep(1:500,4),id=1:2000,seed=820000L+rep(1:500,4))
spec<-expand.grid(Lambda=config$Lambdas,method=config$methods,adjustment=config$adjustments,stringsAsFactors=FALSE)
spec$key<-with(spec,paste(method,adjustment,Lambda,sep='_'))
gate_ids<-c(1L,501L,1001L,1501L)
dictionary<-function(d,cols){
 x<-as.matrix(d[cols]);xx<-as.matrix(d[paste0('X',1:5)])^2
 colnames(xx)<-paste0('X',1:5,'_sq')
 pairs<-combn(ncol(x),2);inter<-sapply(seq_len(ncol(pairs)),function(i)x[,pairs[1,i]]*x[,pairs[2,i]])
 f<-scale(cbind(x,xx,inter));stopifnot(all(is.finite(f)),ncol(f)==if('Z'%in%cols)26L else 20L);f
}
atomic<-function(object,path){
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 tmp<-tempfile('unpublished_',dirname(path));on.exit(if(file.exists(tmp))unlink(tmp))
 withCallingHandlers(saveRDS(object,tmp,compress='gzip'),warning=function(w)stop(conditionMessage(w)))
 stopifnot(file.info(tmp)$size>0,identical(readRDS(tmp),object))
 if(!file.rename(tmp,path))stop('Cannot publish checkpoint: ',path)
 stopifnot(identical(readRDS(path),object));invisible(path)
}
check_record<-function(z,t,s,signature){
 stopifnot(identical(z$task,t),identical(z$signature,signature),
  identical(z$method,s$method),identical(z$adjustment,s$adjustment),identical(z$Lambda,s$Lambda),
  identical(z$iv_msm_compatible,!(t$setting==4L&&s$Lambda==2)))
 e<-z$endpoints
 stopifnot(nrow(e)==2,setequal(e$endpoint,c('ate_lower','ate_upper')),
  is.logical(e$success),!anyNA(e$success),identical(z$success,all(e$success)),
  all(is.finite(e$estimate[e$success])),all(is.finite(e$standard_error[e$success])),
  all(e$standard_error[e$success]>=0),all(nzchar(e$error[!e$success])))
 if(any(e$success)){
  q<-e[e$success,]
  stopifnot(max(abs(q$confidence_lower-(q$estimate-qnorm(.975)*q$standard_error)))<1e-10,
   max(abs(q$confidence_upper-(q$estimate+qnorm(.975)*q$standard_error)))<1e-10)
 }
 invisible(TRUE)
}
load_task<-function(path,t,signature){
 if(!file.exists(path))return(list(task=t,signature=signature,records=list()))
 z<-readRDS(path)
 stopifnot(identical(z$task,t),identical(z$signature,signature),is.list(z$records))
 if(length(z$records)){
  stopifnot(!is.null(names(z$records)),!anyDuplicated(names(z$records)),all(names(z$records)%in%spec$key))
  for(k in names(z$records))check_record(z$records[[k]],t,spec[spec$key==k,,drop=FALSE],signature)
 }
 z
}
run_task<-function(id,root,signature,fit_fun){
 stopifnot(id%in%tasks$id);t<-tasks[id,,drop=FALSE]
 path<-file.path(root,'runs',sprintf('task_%04d.rds',id))
 lock<-file.path(root,'runs',sprintf('task_%04d.lock',id))
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 if(!dir.create(lock,showWarnings=FALSE))stop('Task lock exists; inspect process state: ',lock)
 on.exit(unlink(lock,recursive=TRUE))
 z<-load_task(path,t,signature)
 if(length(z$records)==nrow(spec)){cat('ALREADY COMPLETE',id,'; unchanged\n');return(invisible(z))}
 st<-config$settings[t$setting,];d<-x5g_generate(t$n,st$a,st$b,t$seed,TRUE)
 for(j in seq_len(nrow(spec))){
  s<-spec[j,,drop=FALSE];if(s$key%in%names(z$records))next
  result<-fit_fun(d,t,s,signature)
  check_record(result,t,s,signature);z$records[[s$key]]<-result;atomic(z,path)
  cat(format(Sys.time()),'task',id,'setting',t$setting,'rep',t$replication,
   s$method,s$adjustment,'Lambda',s$Lambda,'success',result$success,'seconds',result$seconds,'\n');flush.console()
 }
 cat('TASK COMPLETE',id,'; verified',length(z$records),'interval fits\n');invisible(z)
}
