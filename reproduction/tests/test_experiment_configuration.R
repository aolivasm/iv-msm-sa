# Configuration and dispatch checks only: no estimation or population integration.
plan <- system2(file.path(R.home('bin'),'Rscript'),
                c('scripts/paper.R','plan'),stdout=TRUE)
stopifnot(!any(grepl('continuous_L|comparators_L|compatibility_[0-9]',plan)),
          sum(grepl('^ *compatibility ',plan))==1L,
          sum(grepl('^ *continuous_bounds ',plan))==1L)
source('scripts/prepare_experiments.R')
local({
  original <- getwd()
  tmp <- tempfile('paper-configuration-')
  dir.create(tmp)
  on.exit({setwd(original);unlink(tmp,recursive=TRUE)})
  for(d in c('scripts','experiments','configurations'))
    stopifnot(file.copy(d,tmp,recursive=TRUE))
  prepare_experiments(tmp)
  prepare_experiments(tmp) # Safe repeat without changing checkpoints.
  settings <- paper_settings(tmp)
  stopifnot(length(settings$profiles)==9L,
            length(unique(vapply(settings$profiles,`[[`,character(1),'engine')))==5L)
  for(key in names(settings$profiles)) {
    p <- settings$profiles[[key]]
    setwd(paper_workspace(key,tmp))
    e <- new.env()
    source(if(p$engine=='binary_comparators')'R/design.R' else 'scripts/common.R',local=e)
    if(p$engine=='compatibility') {
      tasks <- e$n12_tasks()
      stopifnot(all(tasks$n==p$config$n),
                setequal(tasks$outcome,p$config$outcomes),
                setequal(tasks$setting,p$config$settings$setting))
    } else tasks <- if(p$engine=='continuous_bounds')e$stress_tasks() else e$tasks
    stopifnot(nrow(tasks)==p$tasks, !anyDuplicated(tasks$id),
              identical(as.integer(tasks$id),seq_len(p$tasks)))
    if(p$engine=='binary_bounds')
      stopifnot(system2(file.path(R.home('bin'),'Rscript'),
                        c('--vanilla','scripts/tests.R'))==0L)
    if(grepl('^continuous_|^comparators_',key)) {
      stopifnot(all(tasks$seed==820000L+tasks$replication),
                all(tasks$n==6000L),all(table(tasks$setting)==500L))
    }
  }
  setwd(tmp)
  cat('# changed configuration\n',file='configurations/main_experiments.R',append=TRUE)
  error <- tryCatch({prepare_experiments(tmp);NULL},error=identity)
  stopifnot(inherits(error,'error'),
            grepl('changed source/configuration',conditionMessage(error)))
})
cat('PASS: nine profiles, five engines, task contracts, idempotent preparation, changed-source rejection.\n')
