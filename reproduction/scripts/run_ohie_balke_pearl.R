#!/usr/bin/env Rscript
# Standalone runner from the project root. No nuisance refits or old-result writes.
source('R/balke_pearl_bounds.R')
source('R/binary_outcome_bounds.R')
out <- 'results/ohie_balke_pearl/exact_saved_nuisances_v1'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
src <- 'results/ohie_application/ohie_ed_expanded_gam_hal_gamma_fixed_median_minimum_rank10_r2_lambda12_20_v1/nuisance_checkpoint.rds'
cache <- readRDS(src)
dat <- readRDS('data/ohie_application.rds')$ed_household_index_overlap
stopifnot(identical(dat$person_id,cache$contract$person_id),
  identical(cache$nuisance_fit$fold_id,cache$contract$fold_id))
bp <- bp_evaluate(dat,cache$nuisance_fit$predictions)
rm <- evaluate_binary_bounds_from_nuisances(dat,cache$nuisance_fit,
  Lambda0=1.6,Lambda1=1.6,epsilon=0)
rm_signals <- cbind(lower=rm$iv_signals[['1']]$lower-rm$iv_signals[['0']]$upper,
  upper=rm$iv_signals[['1']]$upper-rm$iv_signals[['0']]$lower)
signal_difference <- max(abs(bp$signals-rm_signals))
cat('Maximum BP minus natural row-signal difference:',signal_difference,'\n')
tab <- subset(rm$estimates,method=='iv_only_smooth_one_step' &
  endpoint %in% c('ate_lower','ate_upper'))
tab <- tab[match(c('ate_lower','ate_upper'),tab$endpoint),]
print(tab)
summary <- data.frame(method=c('Robins-Manski','Balke-Pearl'),
  lower=c(tab$estimate[1],bp$estimate[1]),upper=c(tab$estimate[2],bp$estimate[2]),
  se_lower=c(tab$standard_error[1],bp$se[1]),
  se_upper=c(tab$standard_error[2],bp$se[2]),row.names=NULL)
summary$outward_lower <- summary$lower-qnorm(.975)*summary$se_lower
summary$outward_upper <- summary$upper+qnorm(.975)*summary$se_upper
summary$width <- summary$upper-summary$lower
saveRDS(list(summary=summary,balke_pearl=bp,rm_signal_difference=signal_difference,
  source=src,source_md5=tools::md5sum(src),
  contract=cache$contract,person_id=dat$person_id,
  note='Levis exact selected-branch one-step formula with existing nuisance learners, not a rerun of their simulation fitting pipeline. No smoothing or compatibility projection.'),
  file.path(out,'results.rds'))
write.table(summary,file.path(out,'summary.txt'),row.names=FALSE,quote=FALSE)
print(summary,digits=9); print(bp$diagnostics); print(bp$fitted_crossing_rate)
print(lapply(bp$winners,table))
