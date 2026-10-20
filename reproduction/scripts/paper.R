#!/usr/bin/env Rscript
# Default is a plan; no automatic study or scheduler launch.
source('scripts/prepare_experiments.R')
args<-commandArgs(TRUE);settings<-paper_settings();profiles<-settings$profiles
if(!length(args)||args[1]=='plan') {
  print(do.call(rbind,lapply(names(profiles),function(key){
    p<-profiles[[key]];c<-p$config
    data.frame(configuration=key,analysis=p$analysis,n=c$n,
      Lambda=paste(if(is.null(c$Lambda))c$Lambdas else c$Lambda,collapse=','),
      tasks=p$tasks)
  })),row.names=FALSE)
  cat('\nSettings: configurations/main_experiments.R. Run from reproduction/:\n',
      'Rscript scripts/paper.R experiment ANALYSIS_OR_CONFIGURATION preflight\n',
      'Rscript scripts/paper.R experiment ANALYSIS_OR_CONFIGURATION task ID\n',
      'Rscript scripts/paper.R experiment ANALYSIS_OR_CONFIGURATION collect\n',sep='')
  quit(status=0)
}
stopifnot(length(args)>=3,args[1]=='experiment',
          args[3]%in%c('targets','preflight','task','collect'))
keys<-if(args[2]%in%names(profiles))args[2] else
  names(profiles)[vapply(profiles,function(p)p$analysis==args[2],logical(1))]
if(!length(keys))stop('Unknown analysis/configuration. Run paper.R plan.')
if(args[3]=='task') {
  if(length(args)!=4L||!grepl('^[1-9][0-9]*$',args[4]))stop('Supply one positive task ID.')
  if(any(as.integer(args[4])>vapply(profiles[keys],`[[`,integer(1),'tasks')))
    stop('Task ID exceeds a selected configuration. Select a specific configuration.')
} else if(length(args)!=3L)stop('Unexpected arguments.')
root<-normalizePath('.');prepare_experiments(root,keys)
for(key in keys) {
  setwd(paper_workspace(key,root));cat('Configuration:',key,'\n')
  if(args[3]=='targets') {
    script<-if(profiles[[key]]$engine=='binary_bounds')'scripts/population.R' else 'scripts/targets.R'
    if(!file.exists(script))stop('No separate population targets for ',key)
    dir.create('results',showWarnings=FALSE);dir.create('targets',showWarnings=FALSE)
    status<-system2(file.path(R.home('bin'),'Rscript'),script)
  } else status<-system2(file.path(R.home('bin'),'Rscript'),c('scripts/run.R',shQuote(args[-c(1,2)])))
  if(status!=0L)quit(status=status)
}
