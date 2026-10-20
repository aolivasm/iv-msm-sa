# Storage and task contract for the focused two-outcome n=12000 study.
# This file contains no statistical estimator or DGP implementation.
n12_config <- function() list(
  version = 'x5_compatibility_n12000_v1', n = 12000L, reps = 500L,
  outcomes = c('continuous','binary'), methods = c('rkhs','growing_direct'),
  Lambda = c(1,1.1), folds = 5L, bootstrap = 399L, alpha = .05,
  settings = data.frame(setting = c(1:4,0L),
    a = c(.55,1.1,1.65,1.1,1.1), b = c(.075,.075,.075,.15,0)),
  primary = 'pointwise learned-witness maximum; reject',
  continuous_tuning = 'Frozen bandwidth multipliers .25,.5,1,2,4; one 3-fold inner CV partition',
  binary_tuning = 'Frozen kernel variance M/2; penalty-only two independent 3-fold inner CV partitions',
  growing_tuning = 'Frozen sample-size-dependent sparse direct Riesz sieve and selection rule',
  result_root = 'results/x5_compatibility_n12000_v1',
  max_checkpoint_bytes = 8*1024^2,
  note = 'No compatibility projection, classification/PMF-ratio gamma, threshold q, bounds, or comparators. Secondary grid-simultaneous outputs refer only to the two-value grid.')

n12_tasks <- function(config = n12_config()) {
  z <- do.call(rbind,lapply(config$outcomes,function(o) {
    x <- config$settings[rep(seq_len(nrow(config$settings)),each=config$reps),]
    x$outcome <- o; x$replication <- rep(seq_len(config$reps),nrow(config$settings))
    x$type <- ifelse(x$setting==0L,'null','stress')
    source_setting <- ifelse(x$setting==0L,5L,x$setting)
    x$seed <- if(o=='continuous') 820000L+x$replication else
      261000L+source_setting+100L*(x$replication-1L)
    x$parent_task_id <- if(o=='continuous')
      ifelse(x$setting==0L,4000L,2000L+(x$setting-1L)*500L)+x$replication else
      (source_setting-1L)*500L+x$replication
    x
  }))
  rownames(z) <- NULL; z$id <- seq_len(nrow(z)); z$n <- config$n
  z[c('id','outcome','setting','type','a','b','replication','seed','n','parent_task_id')]
}
n12_fingerprint <- function(x) digest::digest(x,algo='sha256',serialize=TRUE)
n12_atomic <- function(object,path,overwrite=TRUE,
                       writer=function(object,path)saveRDS(object,path,compress='gzip')) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  if(file.exists(path) && !overwrite) stop('Refusing to overwrite: ',path)
  tmp <- tempfile('unpublished_',tmpdir=dirname(path))
  on.exit(if(file.exists(tmp))unlink(tmp))
  withCallingHandlers(writer(object,tmp),warning=function(w)stop(conditionMessage(w)))
  if(!file.exists(tmp) || is.na(file.info(tmp)$size) || file.info(tmp)$size<=0)
    stop('Empty unpublished checkpoint: ',tmp)
  stopifnot(identical(readRDS(tmp),object))
  if(!overwrite && file.exists(path)) stop('Destination appeared: ',path)
  if(!file.rename(tmp,path)) stop('Atomic publication failed: ',path)
  stopifnot(identical(readRDS(path),object))
  invisible(path)
}
n12_check_record <- function(v,task,method,signature) {
  cfg <- n12_config()
  stopifnot(is.list(v),!anyDuplicated(names(v)),
    all(c('task','method','signature_hash','n_used','success','error','seconds','warnings')%in%names(v)),
    length(v$n_used)==1L,is.finite(v$n_used),length(v$seconds)==1L,
    is.finite(v$seconds),v$seconds>=0,is.character(v$warnings),!anyNA(v$warnings),
    identical(v$task,task),identical(v$method,method),
    identical(v$signature_hash,n12_fingerprint(signature)),v$n_used==task$n,
    is.logical(v$success),length(v$success)==1L,!is.na(v$success),
    length(v$error)==1L,is.character(v$error),!is.na(v$error))
  if(!v$success) {stopifnot(nzchar(v$error));return(invisible(TRUE))}
  s <- v$summary
  stopifnot(all(c('summary','gamma_selection','paired','fold_sizes','projection_applied','q_fits')%in%names(v)),
    is.data.frame(s),!anyDuplicated(names(s)),
    all(c('Lambda0','Lambda1','statistic','critical_value','p_value','reject','moments_tested')%in%names(s)),
    nrow(s)==length(cfg$Lambda),
    identical(as.numeric(s$Lambda0),cfg$Lambda),identical(as.numeric(s$Lambda1),cfg$Lambda),
    all(is.finite(s$statistic)),all(is.finite(s$critical_value)),
    all(is.finite(s$p_value)),all(s$p_value>=0 & s$p_value<=1),
    is.logical(s$reject),!anyNA(s$reject),all(s$moments_tested==8L),
    isFALSE(v$projection_applied),identical(as.integer(v$q_fits),0L),
    length(v$gamma_selection)==2L*cfg$folds,
    length(v$fold_sizes)==cfg$folds,all(v$fold_sizes>0),sum(v$fold_sizes)==task$n)
  valid_hash <- function(x)is.character(x)&&length(x)==1L&&!is.na(x)&&grepl('^[0-9a-f]{64}$',x)
  p <- v$paired
  stopifnot(is.list(p),valid_hash(p$data),valid_hash(p$folds),
    is.character(p$primitives),length(p$primitives)==cfg$folds,
    !is.null(names(p$primitives)),!anyDuplicated(names(p$primitives)),
    setequal(names(p$primitives),as.character(task$seed+10000L*seq_len(cfg$folds))),
    all(vapply(p$primitives,valid_hash,logical(1))))
  invisible(TRUE)
}
n12_load_task <- function(path,task,signature) {
  if(!file.exists(path))return(list(version=n12_config()$version,task=task,
    signature_hash=n12_fingerprint(signature),n_used=task$n,records=list()))
  z <- readRDS(path)
  stopifnot(is.list(z),!anyDuplicated(names(z)),
    all(c('version','task','signature_hash','n_used','records')%in%names(z)),
    length(z$n_used)==1L,is.finite(z$n_used),
    identical(z$version,n12_config()$version),identical(z$task,task),
    identical(z$signature_hash,n12_fingerprint(signature)),z$n_used==task$n,
    is.list(z$records),!anyDuplicated(names(z$records)),
    !length(z$records)||(!is.null(names(z$records))&&all(names(z$records)%in%n12_config()$methods)))
  for(m in names(z$records))n12_check_record(z$records[[m]],task,m,signature)
  good <- Filter(function(v)isTRUE(v$success),z$records)
  if(length(good)==2L)stopifnot(identical(good[[1L]]$paired,good[[2L]]$paired))
  z
}
n12_run_task <- function(id,signature,bundle_root=getwd(),
                         result_root='results/x5_compatibility_n12000_v1',
                         adapter_factory=n12_adapter_new,smoke=FALSE) {
  cfg <- n12_config(); tasks <- n12_tasks(cfg)
  stopifnot(length(id)==1L,!is.na(id),id%in%tasks$id)
  task <- tasks[id,,drop=FALSE]
  if(smoke)task$n <- 600L
  folder <- if(smoke)'smoke' else 'runs'
  path <- file.path(result_root,folder,sprintf('task_%04d.rds',id))
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  lock <- paste0(path,'.lock')
  if(!dir.create(lock,showWarnings=FALSE))stop('Task lock exists; inspect processes before any recovery: ',lock)
  on.exit(unlink(lock,recursive=TRUE))
  z <- n12_load_task(path,task,signature)
  if(setequal(names(z$records),cfg$methods)) {
    cat('ALREADY ATTEMPTED',id,'; checkpoint unchanged; successful',
        sum(vapply(z$records,function(v)isTRUE(v$success),logical(1))),'of 2\n')
    return(invisible(z))
  }
  prep <- tryCatch({
    a <- adapter_factory(bundle_root=bundle_root,outcome=task$outcome,Lambda=cfg$Lambda,n=task$n)
    list(adapter=a,data=a$generate(task,task$n))
  },error=function(e)e)
  for(method in cfg$methods) {
    if(method%in%names(z$records))next
    warnings <- character(); start <- proc.time()[3L]
    cat(format(Sys.time()),'START',id,task$outcome,task$type,'setting',task$setting,
        'rep',task$replication,'n',task$n,method,'Lambda 1,1.1\n');flush.console()
    fit <- tryCatch(withCallingHandlers({
      if(inherits(prep,'error'))stop(conditionMessage(prep))
      prep$adapter$run(prep$data,method,task$seed)
    },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
      error=function(e)e)
    v <- list(task=task,method=method,signature_hash=n12_fingerprint(signature),
      n_used=task$n,success=!inherits(fit,'error'),
      error=if(inherits(fit,'error'))conditionMessage(fit) else '',
      seconds=unname(proc.time()[3L]-start),warnings=unique(warnings))
    if(v$success) {
      # Whitelist small summaries, not fitted models or observation-level arrays.
      keep <- c('summary','gamma_diagnostics','gamma_selection','probability_diagnostics',
        'settings','paired','fold_sizes','witness_candidate_counts','projection_applied','q_fits')
      v <- c(v,fit[intersect(keep,names(fit))])
      v <- tryCatch({
        n12_check_record(v,task,method,signature)
        for(prior in z$records)if(isTRUE(prior$success))stopifnot(identical(prior$paired,v$paired))
        v
      },error=function(e){v$success<-FALSE;v$error<-paste('Output validation failed:',conditionMessage(e));v})
    }
    n12_check_record(v,task,method,signature)
    z$records[[method]] <- v
    n12_atomic(z,path,writer=function(object,path) {
      saveRDS(object,path,compress='gzip')
      if(file.info(path)$size>cfg$max_checkpoint_bytes)
        stop('Checkpoint exceeds compact-storage limit; original checkpoint preserved.')
    })
    cat(format(Sys.time()),'SAVED',id,task$outcome,method,'success',v$success,
        'seconds',v$seconds,v$error,'\n');flush.console()
  }
  invisible(z)
}

n12_manifest <- function(bundle_root=getwd()) {
  expected <- readRDS(file.path(bundle_root,'manifest.rds'))
  stopifnot(is.character(expected),length(expected)>0,!anyDuplicated(names(expected)),
    all(nzchar(names(expected))),!any(grepl('(^/|(^|/)[.][.](/|$))',names(expected))))
  observed <- tools::md5sum(file.path(bundle_root,names(expected)))
  bad <- names(expected)[is.na(observed)|unname(expected)!=unname(observed)]
  if(length(bad))stop('Sealed source mismatch: ',paste(bad,collapse=', '))
  expected
}
n12_signature <- function(bundle_root=getwd()) {
  packages <- c('SuperLearner','glmnet','ranger','xgboost','nloptr','nnls','digest')
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing))stop('Required R packages missing: ',paste(missing,collapse=', '))
  cfg <- n12_config()
  stopifnot(identical(readRDS(file.path(bundle_root,'study_config.rds')),cfg))
  list(config=cfg,manifest=n12_manifest(bundle_root),
    runtime=list(R=as.character(getRversion()),packages=setNames(vapply(packages,
      function(p)as.character(packageVersion(p)),character(1)),packages)))
}
