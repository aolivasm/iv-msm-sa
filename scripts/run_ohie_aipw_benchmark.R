#!/usr/bin/env Rscript
# The reported conventional AIPW benchmark; not an IV-identified causal effect.
# Reuses the original probability-nuisance recipe, not the later gamma fit.
source('R/binary_outcome_bounds.R')
source('R/ohie_aipw_benchmark.R')
prepared<-readRDS('data/ohie_application.rds')
d<-prepared$ed_household_index_overlap
stopifnot(nrow(d)==22170L,!anyDuplicated(d$household_id))
out<-'results/ohie_aipw_benchmark';dir.create(out,recursive=TRUE,showWarnings=FALSE)
path<-file.path(out,'nuisances.rds')
contract<-list(person_id=d$person_id,x_cols=prepared$ed_overlap_expanded_x_cols,seed=20260814L)
if(file.exists(path)) {
 z<-readRDS(path);stopifnot(identical(z$contract,contract));nf<-z$fit
} else {
 nf<-fit_binary_cross_fitted_nuisances(d,x_cols=contract$x_cols,
  fold_id=binary_make_tz_folds(d,K=5L,seed=contract$seed),spline_df=3L,
  learner='super_learner',overlap_epsilon=.05,outcome_probability_epsilon=1e-4,
  library_candidates=c('glm_spline','gam','glmnet','ranger'),library_cv_folds=3L,
  library_selection_rule='one_se',library_ranger_num_trees=500L,
  super_learner_candidates=c('glm_main','glm_interactions','glmnet_interactions','ranger','xgboost'),
  super_learner_cv_folds=3L,super_learner_ranger_num_trees=500L,seed=contract$seed)
 tmp<-tempfile(tmpdir=out);saveRDS(list(contract=contract,fit=nf),tmp)
 stopifnot(file.rename(tmp,path))
}
tab<-estimate_unconfounded_aipw(d,nf,.95)
write.csv(tab,file.path(out,'summary.csv'),row.names=FALSE)
print(tab,row.names=FALSE)
