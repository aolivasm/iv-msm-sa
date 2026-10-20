#!/usr/bin/env Rscript
# Explicit dispatcher. Default is a plan; no automatic study or scheduler launch.
args<-commandArgs(TRUE)
analyses<-c('continuous_bounds','continuous_bounds_additional_values','binary_bounds',
 'continuous_comparators','continuous_comparators_additional_values','binary_comparators',
 'compatibility_n6000','compatibility_n12000','binary_compatibility_n24000')
if(!length(args)||args[1]=='plan') {
 print(data.frame(analysis=analyses,tasks=c(rep(2000,6),5000,5000,500)),row.names=FALSE)
 cat('\nRun from the repository root:\n',
 'Rscript scripts/paper.R experiment ANALYSIS preflight\n',
 'Rscript scripts/paper.R experiment ANALYSIS task ID\n',
 'Rscript scripts/paper.R experiment ANALYSIS collect\n',
 '\nFigures from existing aggregate results (no fitting):\n',
 'Rscript scripts/render_paper_results.R\n',sep='')
 quit(status=0)
}
if(args[1]!='experiment'||length(args)<3||!args[2]%in%analyses||
 !args[3]%in%c('targets','preflight','task','collect'))stop('Use: paper.R plan, or experiment ANALYSIS targets|preflight|task ID|collect')
if(args[3]=='task') {
 if(length(args)!=4L||!grepl('^[1-9][0-9]*$',args[4]))stop('Supply one positive task ID.')
} else if(length(args)!=3L)stop('Unexpected arguments.')
source('scripts/prepare_experiments.R')
prepare_experiments()
setwd(file.path('experiments',args[2]))
if(args[3]=='targets') {
 script <- if(args[2]=='binary_bounds') 'scripts/population.R' else 'scripts/targets.R'
 if(!file.exists(script))stop('This analysis does not require separate population targets.')
 dir.create('results',showWarnings=FALSE);dir.create('targets',showWarnings=FALSE)
 status<-system2(file.path(R.home('bin'),'Rscript'),script)
 quit(status=status)
}
status<-system2(file.path(R.home('bin'),'Rscript'),c('scripts/run.R',shQuote(args[-c(1,2)])))
quit(status=status)
