#!/usr/bin/env Rscript
# Convert completed private application outputs into aggregate figure inputs.
args<-commandArgs(TRUE);out<-if(length(args))args[1] else 'results/paper_reproduction'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,'ohie_bands.rds')))stop('OHIE aggregate already exists; choose a new directory.')
source('R/binary_continuum_band.R')
paths<-c('IV-MSM'='results/ohie_rkhs_bands/rkhs_grid_band.rds',
 Tan='results/ohie_simultaneous_grid/official_comparators_dense_v1/Tan_x.rds',
 'Dorn-Guo-Kallus'='results/ohie_simultaneous_grid/official_comparators_dense_v1/Dorn_x.rds')
stopifnot(all(file.exists(paths)))
bands<-lapply(paths,function(p){
 z<-readRDS(p);tab<-z$table
 stopifnot(nrow(tab)==81L,max(abs(tab$Lambda-(120:200)/100))<1e-12,
  all(is.finite(tab$simultaneous_lower)),all(is.finite(tab$simultaneous_upper)))
 positive<-tab$simultaneous_lower>0
 tip<-if(any(positive))max(tab$Lambda[positive]+tab$simultaneous_lower[positive]/2) else 1.2
 list(grid=tab,critical=if(!is.null(z$critical))z$critical else z$calibration$critical,
  tipping=tip,at_1p6=binary_continuum_band(tab$Lambda,tab$simultaneous_lower,tab$simultaneous_upper,1.6))
})
iv<-readRDS('results/ohie_balke_pearl/exact_saved_nuisances_v1/results.rds')$summary[1,,drop=FALSE]
path<-paste0('results/ohie_learned_witness_compatibility/',
 'ed_household_index_overlap_expanded_gam_hal_direct_rkhs_half_median_minimum_r2_b999_lambda1_20_v1/compatibility_test.csv')
t<-read.csv(path);t<-subset(t,Lambda0%in%c(1,1.1,1.2))
stopifnot(nrow(t)==3L)
aipw<-read.csv('results/ohie_aipw_benchmark/summary.csv')
stopifnot(nrow(aipw)>0L)
saveRDS(list(bands=bands,iv_only=iv),file.path(out,'ohie_bands.rds'))
write.csv(aipw,file.path(out,'ohie_aipw.csv'),row.names=FALSE)
write.csv(t[,c('Lambda0','Lambda1','statistic','p_value','reject')],file.path(out,'ohie_compatibility.csv'),row.names=FALSE)
cat('Saved aggregate OHIE limits, benchmark, and tests; no participant data exported.\n')
