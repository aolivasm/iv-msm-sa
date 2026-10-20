Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,RCPP_PARALLEL_NUM_THREADS=1)
source('R/design.R')
source('R/dvds_official_authors.R')
# Explicitly constrain the adapter's namespace-qualified forest calls too.
thread_calls<-0L
single_thread<-function(x){
 if(!is.call(x))return(x)
 for(i in seq_along(x)[-1])if(is.call(x[[i]]))x[[i]]<-single_thread(x[[i]])
 if(any(vapply(list(quote(grf::probability_forest),quote(grf::regression_forest),quote(predict)),function(f)identical(x[[1]],f),logical(1)))){
  x$num.threads<-1L;thread_calls<<-thread_calls+1L
 }
 x
}
body(dvds_official_fit_nuisances)<-single_thread(body(dvds_official_fit_nuisances))
stopifnot(thread_calls==8L)
stopifnot(requireNamespace('RCALsa',quietly=TRUE),requireNamespace('grf',quietly=TRUE))
authors<-dvds_official_load()
stopifnot(identical(authors$git_commit,dvds_official_audited_commit))
for(nm in c('quantile_forest','regression_forest','probability_forest'))
 authors$environment[[nm]]<-local({f<-getExportedValue('grf',nm);function(X,...)f(as.matrix(X),...,num.threads=1)})
authors$environment$get_forest_weights<-function(forest,newdata=NULL,...)grf::get_forest_weights(forest,if(is.null(newdata))NULL else as.matrix(newdata),...,num.threads=1)
authors$environment$predict<-function(object,newdata=NULL,...)stats::predict(object,if(is.null(newdata))NULL else as.matrix(newdata),...,num.threads=1)
args<-commandArgs(TRUE);stage<-if(length(args))args[1] else 'preflight'
stopifnot(stage%in%c('preflight','task','gate_check','collect','smoke'))
m<-readRDS('manifest.rds');stopifnot(identical(unname(m),unname(tools::md5sum(names(m)))))
sig<-list(config=config,manifest=m,implementation='official binary Tan and binary DVDS',R=as.character(getRversion()),
 packages=sapply(c('RCALsa','grf'),function(p)as.character(packageVersion(p))),authors=authors$git_commit)
root<-'results/binary_x5_comparators_v1';meta<-file.path(root,'signature.rds')
if(file.exists(meta))stopifnot(identical(readRDS(meta),sig)) else {stopifnot(stage=='preflight');atomic(sig,meta)}
spec<-expand.grid(method=c('tan','dorn'),adjustment=config$adjustments,Lambda=config$Lambda,stringsAsFactors=FALSE)
spec$key<-paste(spec$method,spec$adjustment,spec$Lambda,sep='_')
dictionary<-function(d,cols){x<-as.matrix(d[cols]);sq<-as.matrix(d[paste0('X',1:5)])^2;p<-combn(ncol(x),2)
 f<-scale(cbind(x,sq,sapply(seq_len(ncol(p)),function(i)x[,p[1,i]]*x[,p[2,i]])));stopifnot(all(is.finite(f)));f}
check<-function(z,t){stopifnot(identical(z$signature,sig),identical(z$task,t),is.list(z$records),!anyDuplicated(names(z$records)),all(names(z$records)%in%spec$key))
 for(k in names(z$records)){e<-z$records[[k]]$endpoints;stopifnot(nrow(e)==2,identical(e$endpoint,c('ate_lower','ate_upper')),all(!is.na(e$success)),all(is.finite(e$estimate[e$success])),all(is.finite(e$standard_error[e$success])))}
}
if(stage=='preflight'){
 stopifnot(nrow(tasks)==2000,!anyDuplicated(tasks$seed),all(table(tasks$setting)==500))
 cat('PREFLIGHT PASSED: binary n6000; 2000 datasets, 16000 interval fits; 8 workers after gate.\n');quit(status=0)
}
if(stage%in%c('task','smoke')){
 id<-as.integer(args[2]);stopifnot(id%in%tasks$id);t<-tasks[id,,drop=FALSE]
 run<-function(){
  dir<-file.path(root,if(stage=='smoke')'smoke' else 'runs');dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  p<-file.path(dir,sprintf('task_%04d.rds',id));lock<-paste0(p,'.lock')
  if(!dir.create(lock,showWarnings=FALSE))stop('Task lock exists: ',lock)
  on.exit(unlink(lock,recursive=TRUE))
  z<-if(file.exists(p))readRDS(p) else list(signature=sig,task=t,records=list());check(z,t)
  if(length(z$records)==nrow(spec)){cat('ALREADY COMPLETE',id,'\n');return()}
  d<-generate(t,if(stage=='smoke')600L else 6000L)
  failed<-function(endpoint,msg)data.frame(endpoint=endpoint,success=FALSE,error=msg,estimate=NA_real_,standard_error=NA_real_)
  for(adj in config$adjustments){
   cols<-if(adj=='x')paste0('X',1:5) else c(paste0('X',1:5),'Z');design<-dictionary(d,cols);dcache<-NULL
   for(method in c('tan','dorn'))for(L in config$Lambda){
    key<-paste(method,adj,L,sep='_');if(key%in%names(z$records))next
    start<-proc.time()[3];w<-character()
    e<-withCallingHandlers({
     if(method=='tan')do.call(rbind,lapply(c('lower','upper'),function(bound){
      set.seed(t$seed+1000L+100L*match(bound,c('lower','upper'))+round(10*L))
      tryCatch({f<-RCALsa::saBin.ate.regu.cv(fold=rep(5L,2),y=d$Y,tr=d$T_obs,x=design,gam=L,upper=bound=='upper',nlam=rep(25L,2),lam.fac=rep(2^(-1/4),2))
       est<-unname(f$est['est','diff']);se<-sqrt(unname(f$est['var','diff']));stopifnot(is.finite(est),is.finite(se),se>=0)
       data.frame(endpoint=paste0('ate_',bound),success=TRUE,error='',estimate=est,standard_error=se)
      },error=function(e)failed(paste0('ate_',bound),conditionMessage(e)))
     })) else {
      if(is.null(dcache)){
       set.seed(t$seed+100000L)
       dcache<-tryCatch(estimate_dvds_official_authors(d,cols,'binary',config$Lambda,adjustment=adj,authors_source=authors),error=function(e)e)
      }
      if(inherits(dcache,'error'))failed(c('ate_lower','ate_upper'),conditionMessage(dcache)) else {
       v<-dcache$endpoint_detail;v<-v[v$Lambda==L,];v<-v[match(c('ate_lower','ate_upper'),v$endpoint),]
       stopifnot(nrow(v)==2,all(is.finite(v$estimate)),all(is.finite(v$standard_error)))
       data.frame(endpoint=v$endpoint,success=TRUE,error='',estimate=v$estimate,standard_error=v$standard_error)
      }
     }
    },warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')})
    e$confidence_lower<-e$estimate-qnorm(.975)*e$standard_error;e$confidence_upper<-e$estimate+qnorm(.975)*e$standard_error
    z$records[[key]]<-list(method=method,adjustment=adj,Lambda=L,endpoints=e,warnings=unique(w),seconds=unname(proc.time()[3]-start))
    check(z,t);atomic(z,p);cat('task',id,method,adj,'Lambda',L,'success',all(e$success),'\n');flush.console()
   }
  }
 };run();quit(status=0)
}
ids<-if(stage=='gate_check')c(1,501,1001,1501) else tasks$id
rows<-list()
for(id in ids){z<-readRDS(file.path(root,'runs',sprintf('task_%04d.rds',id)));check(z,tasks[id,,drop=FALSE]);stopifnot(setequal(names(z$records),spec$key))
 for(v in z$records)rows[[length(rows)+1L]]<-cbind(z$task,n=6000,method=v$method,adjustment=v$adjustment,Lambda=v$Lambda,v$endpoints)}
d<-do.call(rbind,rows)
if(stage=='gate_check'){stopifnot(all(d$success));cat('GATE PASSED: four n6000 datasets, 64 intervals.\n');quit(status=0)}
stopifnot(nrow(d)==nrow(tasks)*nrow(spec)*2L,!anyDuplicated(d[c('id','method','adjustment','Lambda','endpoint')]))
atomic(list(signature=sig,endpoints=d,failures=d[!d$success,]),file.path(root,'summary.rds'))
cat(sum(d$success),'of',nrow(d),'endpoints successful; failures retained. No coverage until MSM targets validated.\n')
