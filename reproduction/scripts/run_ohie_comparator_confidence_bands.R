#!/usr/bin/env Rscript
# Official comparator fits, preserved historical outputs, and joint grid calibration.
# Run from repository root. Checkpointed; safe to resume after interruption.
Sys.setenv(OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1,
           VECLIB_MAXIMUM_THREADS=1, MKL_NUM_THREADS=1)
.libPaths(c('results/official_r_library', .libPaths()))
source('R/dvds_official_authors.R')
out <- 'results/ohie_simultaneous_grid/official_comparators_dense_v1'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
grid <- (120:200)/100
prepared <- readRDS('data/ohie_application.rds')
d <- prepared$ed_household_index_overlap
xcols <- prepared$ed_overlap_expanded_x_cols
stopifnot(nrow(d)==22170L, !anyDuplicated(d$household_id))
atomic <- function(x,path) {
  tmp <- tempfile(tmpdir=dirname(path)); saveRDS(x,tmp)
  stopifnot(file.rename(tmp,path))
}
calibrate <- function(S,se) {
  n <- nrow(S); p <- ncol(S)
  # Use each implementation's own SE convention, including RCALsa's divisor n.
  C <- sweep(S,2,colMeans(S))
  V <- crossprod(C)/n^2/outer(se,se)
  sg <- rep(c(1,-1),p/2); V <- V*outer(sg,sg)
  e <- eigen(V,symmetric=TRUE); stopifnot(min(e$values)>-1e-8)
  A <- e$vectors %*% diag(sqrt(pmax(e$values,0)),nrow=p)
  B <- 200000L; set.seed(20260920L); mx <- numeric(B)
  for(start in seq(1L,B,by=2000L)) {
    ix <- start:min(start+1999L,B)
    mx[ix] <- apply(A %*% matrix(rnorm(p*length(ix)),p),2,max)
  }
  list(critical=unname(quantile(mx,.95,type=8)),draws=B,seed=20260920L,
       half_criticals=c(quantile(mx[1:(B/2)],.95,type=8),
                       quantile(mx[(B/2+1):B],.95,type=8)))
}
finish <- function(tab,S,name,provenance) {
  stopifnot(ncol(S)==2*nrow(tab),all(is.finite(S)))
  se <- as.vector(t(as.matrix(tab[,c('se_lower','se_upper')])))
  cal <- calibrate(S,se)
  tab$pointwise_lower <- tab$lower-qnorm(.975)*tab$se_lower
  tab$pointwise_upper <- tab$upper+qnorm(.975)*tab$se_upper
  tab$simultaneous_lower <- tab$lower-cal$critical*tab$se_lower
  tab$simultaneous_upper <- tab$upper+cal$critical*tab$se_upper
  tab$simultaneous_width <- tab$simultaneous_upper-tab$simultaneous_lower
  atomic(list(table=tab,signals=S,calibration=cal,person_id=d$person_id,
    provenance=provenance,scope=paste('Separate asymptotic 95% outward band per method;',
    'requires joint score approximation and nuisance regularity.',
    'Between grid points use right-endpoint limits by population nesting.')),
    file.path(out,paste0(name,'.rds')))
  message(name,' complete; critical=',round(cal$critical,5))
}

# Change only the returned object of a local copy: all official calculations remain intact.
authors <- dvds_official_load(strict_commit=TRUE)
official <- get('dvds.binary',authors$environment)
scorefun <- official; bd <- body(scorefun)
stopifnot(identical(bd[[length(bd)]][[1]],as.name('return')))
ret <- bd[[length(bd)]][[2]]
bd[[length(bd)]] <- substitute(list(bounds=B,signals=cbind(ATEmins,ATEplus)),list(B=ret))
body(scorefun) <- bd
# Limit forest fitting/prediction to one thread without changing the recipe.
singlethread <- function(expr) {
  if(!is.call(expr)) return(expr)
  expr <- as.call(lapply(as.list(expr),singlethread))
  nm <- paste(deparse(expr[[1]]),collapse='')
  if(nm %in% c('grf::probability_forest','predict')) expr$num.threads <- 1L
  expr
}
fitnuis <- dvds_official_fit_nuisances
body(fitnuis) <- singlethread(body(fitnuis))
old_dvds <- paste0('results/ohie_application/',
  'ohie_ed_dvds_official_authors_x_xz_lambda12_20_commit_26cbbf6e_v1')
for(adj in 'x') {
  name <- paste0('Dorn_',adj)
  if(file.exists(file.path(out,paste0(name,'.rds')))) next
  path <- file.path(out,paste0(name,'_nuisances.rds'))
  if(file.exists(path)) {
    cache <- readRDS(path); stopifnot(identical(cache$person_id,d$person_id))
    nu <- cache$nuisance
  } else {
    message('Regenerating official DVDS nuisances: ',adj)
    set.seed(20260814L+100000L*match(adj,c('x','xz')))
    cols <- if(adj=='x')xcols else c(xcols,'Z')
    nu <- fitnuis(as.matrix(d[,cols]),d$T_obs,d$Y,'binary')
    nu <- nu[c('propensity','mu0','mu1','diagnostics')]
    atomic(list(nuisance=nu,person_id=d$person_id,source=authors$git_commit),path)
  }
  S <- matrix(NA_real_,nrow(d),2*length(grid)); rows <- list()
  for(j in seq_along(grid)) {
    args <- list(Y=d$Y,Z=d$T_obs,mu0=nu$mu0,mu1=nu$mu1,e=nu$propensity,Lambda=grid[j])
    f <- do.call(scorefun,args); ref <- do.call(official,args)
    stopifnot(identical(f$bounds,ref))
    est <- colMeans(f$signals); se <- apply(f$signals,2,sd)/sqrt(nrow(d))
    stopifnot(max(abs(est-c(ref$lower,ref$upper)))<1e-12,
      abs(est[1]-qnorm(.975)*se[1]-ref$lower.CI)<1e-12,
      abs(est[2]+qnorm(.975)*se[2]-ref$upper.CI)<1e-12)
    S[,c(2*j-1,2*j)] <- f$signals
    rows[[j]] <- data.frame(Lambda=grid[j],lower=est[1],upper=est[2],se_lower=se[1],se_upper=se[2])
  }
  tab <- do.call(rbind,rows)
  ref <- old$fit$bounds; ix <- match(ref$Lambda,tab$Lambda)
  delta <- max(abs(as.matrix(tab[ix,c('lower','upper')])-as.matrix(ref[,c('lower','upper')])))
  message('DVDS ',adj,' historical endpoint max difference: ',format(delta,digits=10))
  # Do not silently combine a different regeneration with the historical analysis.
  stopifnot(delta<1e-8)
  finish(tab,S,name,list(official_commit=authors$git_commit,historical_max_difference=delta,
    source_md5=tools::md5sum(authors$source_file),recipe='Official GRF recipe, one thread'))
  rm(S,nu); gc()
}

# Reuse existing official Tan fits in a NEW directory, then run the unchanged runner.
for(adj in 'x') {
  name <- paste0('Tan_',adj)
  if(file.exists(file.path(out,paste0(name,'.rds')))) next
  old <- paste0('results/ohie_application/ohie_ed_tan_rcal_official_',adj,'_v2')
  version <- paste0('ohie_ed_tan_rcal_official_',adj,'_dense_bands_v1')
  dest <- file.path('results/ohie_application',version)
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  tag <- function(L)gsub('.', 'p', format(L,trim=TRUE,scientific=FALSE),fixed=TRUE)
  for(L in grid) for(side in c('lower','upper')) {
    rel <- file.path(paste0('lambda_',tag(L)),paste0('endpoint_',side,'.rds'))
    target <- file.path(dest,rel)
    if(file.exists(file.path(old,rel)) && !file.exists(target)) {
      dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
      stopifnot(file.copy(file.path(old,rel),target))
    }
  }
  Sys.setenv(RCALSA_LIBRARY=Sys.getenv('RCALSA_LIBRARY'),
    IVMSM_OHIE_COVARIATE_SPEC='expanded',IVMSM_OHIE_TAN_ADJUSTMENT=adj,
    IVMSM_OHIE_TAN_LAMBDAS=paste(grid,collapse=','),IVMSM_OHIE_TAN_VERSION=version,
    IVMSM_OHIE_TAN_FORCE='false')
  status <- system2(file.path(R.home('bin'),'Rscript'),
    'scripts/run_ohie_tan_rcal_official.R')
  if(status!=0)stop('Official Tan runner failed; completed endpoints retained.')
  S <- matrix(NA_real_,nrow(d),2*length(grid)); rows <- list()
  for(j in seq_along(grid)) {
    est <- se <- numeric(2)
    for(k in 1:2) {
      side <- c('lower','upper')[k]
      f <- readRDS(file.path(dest,paste0('lambda_',tag(grid[j])),paste0('endpoint_',side,'.rds')))
      s <- f$fit$influ[,2]-f$fit$influ[,1]
      # RCALsa's relaxed estimate differs from mean(influ); use its reported estimate.
      # Centering the official influence contributions preserves its variance exactly.
      v <- mean((s-f$fit$est['est0','diff'])^2)/nrow(d)
      stopifnot(abs(sqrt(v)-f$standard_error)<1e-12)
      S[,2*j-2+k] <- s; est[k] <- f$estimate; se[k] <- f$standard_error
    }
    rows[[j]] <- data.frame(Lambda=grid[j],lower=est[1],upper=est[2],se_lower=se[1],se_upper=se[2])
  }
  finish(do.call(rbind,rows),S,name,list(source=dest,package='RCALsa 1.0',
    estimator='sa.ate.regu.cv, est row (relaxed), official influ and var'))
  rm(S);gc()
}
message('Main-text X-adjusted comparator bands complete.')
