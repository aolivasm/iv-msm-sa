source('scripts/common.R')
stopifnot(nrow(tasks) == nrow(config$settings)*500L,
  !anyDuplicated(tasks$seed), all(table(tasks$setting)==500L),
  identical(config$primitives, c('glm_main','glm_interactions','glmnet_interactions','ranger','xgboost')),
  config$folds == 5L, config$bootstrap == 399L,
  identical(config$settings$a[1:4], c(.55,1.1,1.65,1.1)),
  identical(config$settings$b[1:4], c(.075,.075,.075,.15)))
# The first four designs reproduce the exact completed endpoint-study data.
old_env <- new.env(parent=globalenv())
old_env$config <- config; old_env$config$settings <- config$settings[1:4, , drop=FALSE]
old_generate <- eval(extract_generator('provenance/binary_stress_common.R'), old_env)
for(id in gate_ids[1:4]) {
  task <- tasks[id, , drop=FALSE]
  a <- generate(task,600L); b <- old_generate(task,600L)
  stopifnot(identical(a,b), all(a$Y %in% 0:1), all(a$T_obs %in% 0:1),
    all(a$Z %in% 0:1), all(abs(as.matrix(a[paste0('X',1:5)])) <= 2))
  fold <- make_tz_stratified_folds(a,5,task$seed)
  stopifnot(length(fold)==600L,all(table(fold,a$T_obs,a$Z)>0))
}
stopifnot(isTRUE(all.equal(default_gamma_lambda_grid(1000L)*sqrt(1000L),10^seq(-1.5,1,.5))))
# Oracle binary interval-intersection checks, not a fitted gamma classifier.
d <- generate(tasks[1,,drop=FALSE],600L)
truth <- binary_truth(d,1L)
for(t in 0:1) {
  q <- truth$outcome[[t+1L]]
  for(L in config$Lambda) {
    r0 <- if(t==1) truth$pi0 else 1-truth$pi0
    r1 <- if(t==1) truth$pi1 else 1-truth$pi1
    bounds <- function(p,r) cbind(lower=r*p+(1-r)*pmax(p/L,1-L*(1-p)),
                                  upper=r*p+(1-r)*pmin(L*p,1-(1-p)/L))
    e0<-bounds(q$p0,r0); e1<-bounds(q$p1,r1)
    stopifnot(all(is.finite(e0)),all(is.finite(e1)),
      all(e0[,1]<=e0[,2]+1e-12),all(e1[,1]<=e1[,2]+1e-12))
    if(L==1) stopifnot(max(abs(e0[,1]-q$p0),abs(e0[,2]-q$p0))<1e-12)
  }
}
# Exhaust the eight atoms of a binary (Z,T,Y) law, with exact population
# weights. These are oracle identities, not fitted PMF-ratio gamma learners,
# estimator runs, Monte Carlo draws, or a finite-sample validity claim.
oracle_engine <- new_session()$engine
check_binary_oracle_fixture <- function(q0, q1, Lambda, compatible,
                                        violation = NULL, r0 = .3, r1 = .6) {
  tolerance <- 1e-12
  f0 <- c(1-q0,q0); f1 <- c(1-q1,q1)
  gamma <- f1/f0
  interval <- function(q,r) c(lower=r*q+(1-r)*max(q/Lambda,1-Lambda*(1-q)),
                              upper=r*q+(1-r)*min(Lambda*q,1-(1-q)/Lambda))
  i0 <- interval(q0,r0); i1 <- interval(q1,r1)
  intersection <- c(lower=max(i0[1],i1[1]),upper=min(i0[2],i1[2]))
  stopifnot(identical(unname(intersection[1] <= intersection[2]+tolerance),compatible),
    abs(sum(f0*gamma)-1)<tolerance)
  for(t in 0:1) {
    # Exercise both the treatment-one and treatment-zero indexing branches.
    pi0 <- if(t==1L) r0 else 1-r0
    pi1 <- if(t==1L) r1 else 1-r1
    nuisance <- list(pi0=rep(pi0,2L),pi1=rep(pi1,2L))
    envelope <- compute_projected_envelopes(nuisance,gamma,t,Lambda,Lambda,
      gamma_compatibility_action='stop')
    scores <- oracle_engine$compute_learned_witness_scores(gamma,
      nuisance$pi0,nuisance$pi1,t,Lambda)
    stopifnot(!any(envelope$gamma_projected),identical(envelope$gamma_used,gamma),
      max(abs(scores$pointwise_lower-(envelope$gamma_compatible_lower-gamma)))<tolerance,
      max(abs(scores$pointwise_upper-(gamma-envelope$gamma_compatible_upper)))<tolerance)
    if(Lambda>1) stopifnot(
      max(abs(scores$normalization_lower-(envelope$candidate_lower-1/Lambda)))<tolerance,
      max(abs(scores$normalization_upper-(Lambda-envelope$candidate_upper)))<tolerance)
    expected <- c(pointwise_lower=sum(f0*pmax(scores$pointwise_lower,0)),
      pointwise_upper=sum(f0*pmax(scores$pointwise_upper,0)),
      normalization_lower=sum(f0*envelope$L)-1,
      normalization_upper=1-sum(f0*envelope$U))
    # Positive-part witnesses attain the exact normalization violations and
    # detect any binary-support pointwise crossing without projection.
    stopifnot(identical(all(expected<=tolerance),compatible))
    if(!is.null(violation)) {
      if(violation=='pointwise') stopifnot(max(expected[1:2])>tolerance,
                                          any(envelope$L>envelope$U+tolerance)) else {
        stopifnot(max(expected[1:2])<=tolerance,expected[[violation]]>tolerance,
          all(envelope$L<=envelope$U+tolerance))
      }
    }
    atoms <- expand.grid(Y=0:1,T_obs=0:1,Z=0:1)
    n <- nrow(atoms); e <- .4
    pi <- ifelse(atoms$Z==0,pi0,pi1)
    q <- ifelse(atoms$Z==0,q0,q1)
    weights <- ifelse(atoms$Z==0,e,1-e)*
      ifelse(atoms$T_obs==1,pi,1-pi)*ifelse(atoms$Y==1,q,1-q)
    stopifnot(abs(sum(weights)-1)<tolerance,all(weights>0))
    primitive <- compute_dual_msm_moment_primitives(r0,r1,Lambda,Lambda)
    for(family in names(scores)) {
      witness <- as.numeric(scores[[family]]>0)
      m0 <- sum(f0*witness); m1 <- sum(f1*witness)
      base <- oracle_engine$compute_learned_family_base(family,primitive,m0,m1,Lambda)
      block <- make_arm_lambda_signal_blocks(atoms,
        H=matrix(1,n,1,dimnames=list(NULL,'one')),
        R=matrix(witness[atoms$Y+1L],n,1,dimnames=list(NULL,'positive')),
        m0=matrix(m0,n,1),m1=matrix(m1,n,1),e=rep(e,n),
        pi0=rep(pi0,n),pi1=rep(pi1,n),t_value=t,Lambda0=Lambda,Lambda1=Lambda)
      keep <- block$metadata$family==family
      stopifnot(sum(keep)==1L,abs(base-expected[[family]])<tolerance,
        abs(sum(weights*block$signals[,keep])-expected[[family]])<tolerance)
    }
  }
  invisible(TRUE)
}
check_binary_oracle_fixture(.25,.25,1,TRUE)
check_binary_oracle_fixture(.25,.30,2,TRUE)
check_binary_oracle_fixture(.20,.80,1.1,FALSE,'pointwise')
check_binary_oracle_fixture(.32,.68,2,FALSE,'normalization_upper')
check_binary_oracle_fixture(.36,.63,2,FALSE,'normalization_lower',r0=.6,r1=.6)
labels<-population_labels()
stopifnot(all(grepl('incompatible',labels$label[labels$setting<=4 & labels$Lambda0%in%c(1,1.1)])))
if(nrow(config$settings)==5L) {
  stopifnot(config$settings$b[5]==0,config$settings$a[5]==1.1)
  null<-binary_truth(generate(tasks[2001,,drop=FALSE],600),5)
  stopifnot(all(vapply(null$outcome,function(z)identical(z$p0,z$p1),logical(1))))
}
# Isolated mocked persistence tests: no fitted estimator runs here.
saved_config <- config
testdir <- tempfile('binary_compat_unit_'); dir.create(testdir)
config$result_root <- file.path(testdir,'results')
sig <- list(test=TRUE,config=config)
calls <- 0L
fake <- function() list(run=function(data,method,seed) {
  calls <<- calls+1L
  list(summary=data.frame(Lambda0=config$Lambda,Lambda1=config$Lambda,
    p_value=rep(.5,5),reject=rep(FALSE,5)),
    gamma_selection=rep(list(list(test=TRUE)),10),
    paired=list(data=fingerprint(data),
      folds=fingerprint(rep(seq_len(5),length.out=nrow(data))),
      primitives=setNames(rep(fingerprint('test primitives'),5),
        as.character(seed+10000L*seq_len(5)))),
    fold_sizes=rep(nrow(data)/5,5),projection_applied=FALSE,q_fits=0L)
})
quiet_run <- function(...) {
  invisible(capture.output(value <- run_task(...)))
  value
}
first <- quiet_run(1L,sig,smoke=TRUE,session_factory=fake)
p<-file.path(config$result_root,'smoke','task_0001.rds')
before<-tools::md5sum(p)
again <- quiet_run(1L,sig,smoke=TRUE,session_factory=function()stop('Must not refit'))
stopifnot(calls==2L,identical(first,again),identical(before,tools::md5sum(p)))
lock<-paste0(p,'.lock');dir.create(lock)
locked<-try(run_task(1L,sig,smoke=TRUE,session_factory=fake),silent=TRUE)
stopifnot(inherits(locked,'try-error'));unlink(lock,recursive=TRUE)
failed_factory<-function()list(run=function(...)stop('intentional retained failure'))
failed<-quiet_run(2L,sig,smoke=TRUE,session_factory=failed_factory)
failed_again<-quiet_run(2L,sig,smoke=TRUE,session_factory=function()stop('Must not retry failure'))
stopifnot(identical(failed,failed_again),all(!vapply(failed$records,`[[`,logical(1),'success')))
malformed_factory<-function()list(run=function(data,method,seed) {
  fit<-fake()$run(data,method,seed)
  fit$summary$p_value[1]<-NA_real_
  fit
})
malformed<-quiet_run(3L,sig,smoke=TRUE,session_factory=malformed_factory)
malformed_path<-file.path(config$result_root,'smoke','task_0003.rds')
malformed_before<-tools::md5sum(malformed_path)
malformed_again<-quiet_run(3L,sig,smoke=TRUE,
  session_factory=function()stop('Must not retry malformed output'))
stopifnot(all(!vapply(malformed$records,`[[`,logical(1),'success')),
  all(vapply(malformed$records,function(v)startsWith(v$error,'Output validation failed:'),logical(1))),
  identical(readRDS(malformed_path),malformed),identical(malformed,malformed_again),
  identical(malformed_before,tools::md5sum(malformed_path)))
invalid<-first;invalid$records$rkhs$summary$Lambda0[1]<-1.01
stopifnot(inherits(try(check_record(invalid$records$rkhs,invalid$task,'rkhs',sig,600L),silent=TRUE),'try-error'))
config <- saved_config
# Only this freshly created temporary unit-test directory is removed.
unlink(testdir,recursive=TRUE)
cat('TESTS PASSED: frozen generator, grids, binary oracle intersections/envelopes/moments, no-overwrite resume, estimator/validation failure preservation and locks.\n')
