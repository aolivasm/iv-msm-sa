#!/usr/bin/env Rscript
# Default is a plan; no automatic study or scheduler launch.
source('scripts/prepare_experiments.R')
args<-commandArgs(TRUE);settings<-paper_settings();profiles<-settings$profiles
analyses<-unique(vapply(profiles,`[[`,character(1),'analysis'))
if(!length(args)||args[1]=='plan') {
  print(do.call(rbind,lapply(analyses,function(analysis){
    selected<-profiles[vapply(profiles,function(p)p$analysis==analysis,logical(1))]
    data.frame(analysis=analysis,
      task_IDs=paste0('1–',max(vapply(selected,`[[`,integer(1),'tasks'))))
  })),row.names=FALSE)
  cat('\nEach analysis includes the settings reported in the paper.\nRun from reproduction/:\n',
      'Rscript scripts/paper.R experiment ANALYSIS preflight\n',
      'Rscript scripts/paper.R experiment ANALYSIS task ID\n',
      'Rscript scripts/paper.R experiment ANALYSIS collect\n',sep='')
  quit(status=0)
}
stopifnot(length(args)>=3,args[1]=='experiment',
          args[3]%in%c('targets','preflight','task','collect'))
keys<-names(profiles)[vapply(profiles,function(p)p$analysis==args[2],logical(1))]
if(!length(keys))stop('Unknown analysis. Run paper.R plan.')
if(args[3]=='task') {
  if(length(args)!=4L||!grepl('^[1-9][0-9]*$',args[4]))stop('Supply one positive task ID.')
  keys<-keys[as.integer(args[4])<=vapply(profiles[keys],`[[`,integer(1),'tasks')]
  if(!length(keys))stop('Task ID exceeds the selected configuration(s).')
} else if(length(args)!=3L)stop('Unexpected arguments.')
root<-normalizePath('.');prepare_experiments(root,keys)
for(key in keys) {
  setwd(paper_workspace(key,root))
  if(args[3]=='targets') {
    script<-if(profiles[[key]]$engine=='binary_bounds')'scripts/population.R' else 'scripts/targets.R'
    if(!file.exists(script))stop('No separate population targets for ',key)
    dir.create('results',showWarnings=FALSE);dir.create('targets',showWarnings=FALSE)
    status<-system2(file.path(R.home('bin'),'Rscript'),script)
  } else status<-system2(file.path(R.home('bin'),'Rscript'),c('scripts/run.R',shQuote(args[-c(1,2)])))
  if(status!=0L)quit(status=status)
}
