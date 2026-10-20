Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,
 VECLIB_MAXIMUM_THREADS=1,RCPP_PARALLEL_NUM_THREADS=1)
source('scripts/common.R')
manifest<-readRDS('manifest.rds')
stopifnot(identical(unname(manifest),unname(tools::md5sum(names(manifest)))))
source('R/dvds_official_authors.R');stopifnot(requireNamespace('RCALsa',quietly=TRUE))
if(paper_profile$single_thread_nuisance){
 source('scripts/grf_threads.R')
 dvds_official_fit_nuisances<-single_thread_dvds_nuisances(dvds_official_fit_nuisances)
}
authors<-dvds_official_load()
for(nm in c('quantile_forest','regression_forest','probability_forest'))
 authors$environment[[nm]]<-local({f<-getExportedValue('grf',nm);function(X,...)f(as.matrix(X),...,num.threads=1)})
authors$environment$get_forest_weights<-function(forest,newdata=NULL,...)
 grf::get_forest_weights(forest,if(is.null(newdata))NULL else as.matrix(newdata),...,num.threads=1)
authors$environment$predict<-function(object,newdata=NULL,...)
 stats::predict(object,if(is.null(newdata))NULL else as.matrix(newdata),...,num.threads=1)
versions<-sapply(c('RCALsa','grf'),function(p)as.character(packageVersion(p)))
baseline<-NULL # No private historical checkpoint required.
signature<-list(config=config,manifest=manifest,R=as.character(getRversion()),packages=versions,
 dvds_commit=authors$git_commit,baseline_signature=baseline,
 dictionary='Tan: standardized main effects, X squares, all pairwise interactions; no Z square')
args<-commandArgs(TRUE);stage<-if(length(args))args[1] else 'preflight'
stopifnot(stage%in%c('preflight','task','gate_check','collect'))
root<-paper_profile$result_root;dir.create(root,recursive=TRUE,showWarnings=FALSE)
meta<-file.path(root,'signature.rds')
if(file.exists(meta))stopifnot(identical(readRDS(meta),signature)) else {
 stopifnot(stage=='preflight');atomic(signature,meta)
}
fit_one<-function(d,t,s,signature){
 method<-s$method;adj<-s$adjustment;L<-s$Lambda
 cols<-if(adj=='x')paste0('X',1:5) else c(paste0('X',1:5),'Z')
 start<-proc.time()[3];warnings<-character();detail<-list()
 if(method=='tan'){
  design<-dictionary(d,cols)
  for(bound in c('lower','upper')){
   set.seed(t$seed+1000L+100L*match(bound,c('lower','upper'))+round(10*L))
   detail[[bound]]<-tryCatch(withCallingHandlers({
    f<-RCALsa::sa.ate.regu.cv(fold=rep(5L,3),y=d$Y,tr=d$T_obs,x=design,gam=L,
     upper=bound=='upper',nlam=rep(25L,3),lam.fac=rep(2^(-1/4),3),ploss='cal',yloss='gaus',option='qu')
    est<-unname(f$est['est','diff']);se<-sqrt(unname(f$est['var','diff']))
    stopifnot(is.finite(est),is.finite(se),se>=0)
    data.frame(endpoint=paste0('ate_',bound),success=TRUE,error='',estimate=est,standard_error=se)
   },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
    error=function(e)data.frame(endpoint=paste0('ate_',bound),success=FALSE,error=conditionMessage(e),estimate=NA_real_,standard_error=NA_real_))
  }
 } else {
  set.seed(t$seed+100000L)
  v<-tryCatch(withCallingHandlers({
   f<-estimate_dvds_official_authors(d,cols,'continuous',L,adjustment=adj,authors_source=authors)
   e<-f$endpoint_detail;stopifnot(nrow(e)==2,all(is.finite(e$estimate)),all(is.finite(e$standard_error)),all(e$standard_error>=0))
   data.frame(endpoint=e$endpoint,success=TRUE,error='',estimate=e$estimate,standard_error=e$standard_error)
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
   error=function(e)data.frame(endpoint=c('ate_lower','ate_upper'),success=FALSE,error=conditionMessage(e),estimate=NA_real_,standard_error=NA_real_))
  detail<-list(v)
 }
 e<-do.call(rbind,detail);e$confidence_lower<-e$estimate-qnorm(.975)*e$standard_error
 e$confidence_upper<-e$estimate+qnorm(.975)*e$standard_error
 st<-config$settings[t$setting,]
 list(signature=signature,task=t,a=st$a,b=st$b,method=method,adjustment=adj,Lambda=L,
  endpoints=e,success=all(e$success),warnings=unique(warnings),seconds=unname(proc.time()[3]-start),
  iv_msm_compatible=!(t$setting==4L&&L%in%paper_profile$incompatible_Lambda))
}
if(stage=='preflight'){
 for(id in gate_ids){t<-tasks[id,];st<-config$settings[t$setting,];d<-x5g_generate(t$n,st$a,st$b,t$seed,TRUE)
  stopifnot(nrow(d)==config$n,all(is.finite(dictionary(d,paste0('X',1:5)))),all(is.finite(dictionary(d,c(paste0('X',1:5),'Z')))))}
 cat('PREFLIGHT PASSED:',nrow(tasks),'datasets; n =',config$n,'; Lambda =',
     paste(config$Lambdas,collapse=','),';',nrow(tasks)*nrow(spec),'interval fits.\n')
 cat('Exact frozen LPC stress generator; baseline Tan/Dorn tuning and software; compact checkpoints.\n');quit(status=0)
}
if(stage=='task'){run_task(as.integer(args[2]),root,signature,fit_one);quit(status=0)}
ids<-if(stage=='gate_check')gate_ids else tasks$id
rows<-list();diagnostics<-list();k<-0L
for(id in ids){
 t<-tasks[id,,drop=FALSE];p<-file.path(root,'runs',sprintf('task_%04d.rds',id))
 z<-load_task(p,t,signature);stopifnot(length(z$records)==nrow(spec))
 for(key in spec$key){v<-z$records[[key]];k<-k+1L
  if(stage=='gate_check'&&!v$success)stop('Gate numerical failure: task ',id,' ',key,'; retained for inspection.')
  rows[[k]]<-cbind(t,a=v$a,b=v$b,method=v$method,adjustment=v$adjustment,Lambda=v$Lambda,
    iv_msm_compatible=v$iv_msm_compatible,v$endpoints)
  diagnostics[[k]]<-list(id=id,key=key,seconds=v$seconds,warnings=v$warnings)
 }
}
ep<-do.call(rbind,rows);rownames(ep)<-NULL
if(stage=='gate_check'){cat('GATE PASSED: four datasets; all 32 interval fits / 64 endpoints successful.\n');quit(status=0)}
stopifnot(nrow(ep)==32000,!anyDuplicated(ep[c('id','method','adjustment','Lambda','endpoint')]))
atomic(list(signature=signature,endpoints=ep,failures=ep[!ep$success,],diagnostics=diagnostics),file.path(root,'summary.rds'))
cat(sum(ep$success),'of',nrow(ep),'endpoints successful; failures retained.\n')
cat('Stress-specific conventional MSM targets must be validated before coverage; never reuse baseline or IV-MSM targets.\n')
