# Main-text compatibility session; no historical simulation driver.
source('R/design.R')
source('R/nuisance_models.R');source('R/rkhs_gamma.R');source('R/compatibility_test.R')
config <- list(folds=5L,bootstrap=399L,allocation_df=5L,primitives=sfp_contract()$primitive_library,methods='rkhs_tuned')
source('R/session.R')
