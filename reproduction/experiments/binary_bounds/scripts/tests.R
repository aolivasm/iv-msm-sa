# Run from a prepared binary workspace. No datasets, targets, or fits are generated.
source('scripts/common.R')
stopifnot(nrow(tasks)==2000L, !anyDuplicated(tasks$seed),
          all(table(tasks$setting)==500L),
          identical(config$methods,'rkhs'),
          identical(config$Lambda,c(1.8,2,2.5,3)))

# Resume must return before generating data or fitting nuisance functions.
generate <- function(...) stop('TEST GUARD: data generation is forbidden')
fit_binary_cross_fitted_direct_rkhs_nuisances <-
  function(...) stop('TEST GUARD: model fitting is forbidden')

local({
  root <- tempfile('binary-checkpoint-test-')
  dir.create(root)
  on.exit(unlink(root,recursive=TRUE))
  t <- tasks[1,,drop=FALSE]
  sig <- list(test=TRUE)
  p <- file.path(root,'task_0001.rds')
  expect_error <- function(expr) {
    stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),'error'))
  }
  record <- function(value) {
    list(task=t,signature=sig,
         records=setNames(rep(list(value),length(config$methods)),config$methods))
  }
  resume_unchanged <- function(z) {
    atomic(z,p)
    hash <- tools::md5sum(p)
    returned <- run_task(1,root,sig)
    stopifnot(identical(returned,z),identical(hash,tools::md5sum(p)),
              !dir.exists(paste0(p,'.lock')))
  }

  # Failures are retained for inspection, not silently refitted.
  failed <- record(list(success=FALSE,error='deliberate fixture failure'))
  resume_unchanged(failed)
  stopifnot(!readRDS(p)$records$rkhs$success)

  # Synthetic endpoint records exercise checkpoint validation, not estimation.
  successful <- record(list(success=TRUE,error='',fits=lapply(config$Lambda,function(L)
    list(Lambda0=L,endpoint_intervals=data.frame(
      endpoint=c('ate_lower','ate_upper'),estimate=c(-.1,.2),
      standard_error=c(.01,.02))))))
  resume_unchanged(successful)
  expect_error(check(successful,t,list(wrong=TRUE)))
  expect_error(check(successful,tasks[2,,drop=FALSE],sig))
  invalid <- successful
  names(invalid$records) <- 'unknown_learner'
  expect_error(check(invalid,t,sig))
  invalid <- successful
  invalid$records <- c(invalid$records,invalid$records)
  expect_error(check(invalid,t,sig))
  invalid <- successful
  invalid$records$rkhs$fits[[1]]$endpoint_intervals$standard_error[1] <- -1
  expect_error(check(invalid,t,sig))
  invalid <- successful
  invalid$records$rkhs$fits[[1]]$Lambda0 <- 99
  expect_error(check(invalid,t,sig))
})
cat('PASS: binary task mapping, single-learner fixtures, atomic resume, retained failures, invalid-checkpoint rejection; no simulations or fits.\n')
