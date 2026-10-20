#!/usr/bin/env Rscript
source('R/dvds_official_authors.R')
pk<-c('RCALsa','RCAL','lpSolveAPI','grf','dplyr','progress','ggplot2')
ok<-vapply(pk,requireNamespace,logical(1),quietly=TRUE)
if(!all(ok))stop('Missing comparator dependencies: ',paste(pk[!ok],collapse=', '),
 '. See docs/comparators.md.')
stopifnot(as.character(packageVersion('RCALsa'))=='1.0')
a<-dvds_official_load(strict_commit=TRUE)
print(data.frame(package=pk,version=vapply(pk,function(p)as.character(packageVersion(p)),character(1))),row.names=FALSE)
cat('Official DVDS commit verified: ',a$git_commit,'\nNo fits run.\n')
