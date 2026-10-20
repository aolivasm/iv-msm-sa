# Experimental comparison. Production test and previous bundles remain unchanged.
source('R/scalar_gamma_compatibility_comparison.R')
source('R/learned_witness_compatibility_test.R')

sw_probability <- function(train,validation,seed) learned_witness_probability_fit(
  training=train,validation=validation,x_cols='X',learner='super_learner',
  positivity_epsilon=.05,treatment_spline_df=5,gam_basis_dimension=6,
  integrated_confounder_strength=0,integrated_quadrature_nodes=21,
  super_learner_candidates=sfp_contract()$primitive_library,
  super_learner_cv_folds=3,super_learner_ranger_num_trees=500,library_seed=seed)

# All candidate functions/cutoffs and nuisance fits use development data.
# Only the internal selection rows are used to rank their one-step signals.
# Final-test rows are never inspected by the criterion.
sw_selector <- function(selection_data,selection_probability,selection_rows) {
  force(selection_data);force(selection_probability);force(selection_rows)
  function(family,response_training,response_validation,mean_fit,
           primitive_training,primitive_validation,Lambda,localization_tail_fractions) {
    candidates<-list()
    for(column in seq_len(ncol(response_training))) {
      bt<-compute_learned_family_base(family,primitive_training,
        mean_fit$train_z0[,column],mean_fit$train_z1[,column],Lambda)
      bv<-compute_learned_family_base(family,primitive_validation,
        mean_fit$validation_z0[,column],mean_fit$validation_z1[,column],Lambda)
      cuts<-make_localization_cutoffs(bt,family,localization_tail_fractions)
      for(h_name in names(cuts)) {
        ht<-as.numeric(bt>cuts[h_name]);hv<-as.numeric(bv>cuts[h_name])
        if(mean(ht)<.025) next
        block<-make_arm_lambda_signal_blocks(selection_data,
          H=matrix(hv[selection_rows],ncol=1,dimnames=list(NULL,'h')),
          R=response_validation[selection_rows,column,drop=FALSE],
          m0=mean_fit$validation_z0[selection_rows,column,drop=FALSE],
          m1=mean_fit$validation_z1[selection_rows,column,drop=FALSE],
          e=selection_probability$e,pi0=selection_probability$pi0,
          pi1=selection_probability$pi1,t_value=attr(selection_data,'treatment'),
          Lambda0=Lambda,Lambda1=Lambda)
        signal<-block$signals[,block$metadata$family==family]
        se<-sd(signal)/sqrt(length(signal))
        criterion<-if(is.finite(se)&&se>1e-10) mean(signal)/se else -Inf
        candidates[[length(candidates)+1L]]<-list(column=column,
          r_name=colnames(response_training)[column],h_name=h_name,
          h_cutoff=unname(cuts[h_name]),criterion=criterion,h_validation=hv,
          base_training_mean=mean(ht*bt),base_training_standard_deviation=sd(ht*bt),
          h_training_fraction=mean(ht),selection_mean=mean(signal),selection_se=se)
      }
    }
    stopifnot(length(candidates)>0)
    selected<-candidates[[which.max(vapply(candidates,`[[`,numeric(1),'criterion'))]]
    selected$R_validation<-response_validation[,selected$column,drop=FALSE]
    selected$m0_validation<-mean_fit$validation_z0[,selected$column,drop=FALSE]
    selected$m1_validation<-mean_fit$validation_z1[,selected$column,drop=FALSE]
    selected$H_validation<-matrix(selected$h_validation,ncol=1,
      dimnames=list(NULL,paste0('learned_',family)))
    selected
  }
}

sw_run <- function(data,variant,grid,seed,bootstrap=399L) {
  stopifnot(variant %in% c('current','revised_crossfit','revised_honest'))
  if(variant=='current') {
    f<-sgc_session()$run(data,'riesz_spline',grid,seed,bootstrap)
    f$settings$variant<-variant;return(f)
  }
  engine<-new.env(parent=globalenv())
  sys.source('R/learned_witness_compatibility_test.R',engine)
  honest<-variant=='revised_honest'
  outer<-make_tz_stratified_folds(data,K=if(honest) 2 else 5,seed=seed)
  signals<-metadata<-diagnostics<-gamma_diagnostics<-splits<-list()
  for(fold in if(honest) 1L else 1:5) {
    ti<-which(outer!=fold);vi<-which(outer==fold)
    inner<-make_tz_stratified_folds(data[ti,],K=2,seed=seed+200000L+fold)
    di<-ti[inner==1];si<-ti[inner==2]
    stopifnot(!length(intersect(di,si)),!length(intersect(c(di,si),vi)))
    dev<-data[di,];sel<-data[si,];test<-data[vi,]
    evaluation<-rbind(sel,test);sr<-seq_len(nrow(sel));vr<-nrow(sel)+seq_len(nrow(test))
    prob<-sw_probability(dev,evaluation,seed+10000L*fold)
    selection_prob<-lapply(prob$validation,function(v)v[sr])
    local_signals<-local_metadata<-list()
    for(t in 0:1) {
      g<-sg_fit(dev,prob$training,t,'riesz_spline',seed+10000L*fold+1000L*t)
      gt<-g$predict(dev)$normalized;gv<-g$predict(evaluation)$normalized
      selection_data<-sel;attr(selection_data,'treatment')<-t
      engine$select_learned_witness<-sw_selector(selection_data,selection_prob,sr)
      gamma_diagnostics[[length(gamma_diagnostics)+1]]<-list(fold=fold,treatment=t,
        selected_lambda=g$selected_lambda,cv=g$cv_table,cv_resamples=g$cv_resample_table,
        training_clip_rate=g$training_clip_rate,training_negative_rate=g$training_negative_rate)
      for(L in grid) {
        local<-engine$fit_fold_learned_witness_signals(dev,evaluation,
          prob$training,prob$validation,gt,gv,t,L,'X',5,c(.05,.1,.2,.4),
          c(.1,.2,.4),c('indicator','linear','sqrt'))
        j<-length(local_signals)+1L
        local_signals[[j]]<-local$signals[vr,,drop=FALSE]
        local_metadata[[j]]<-local$metadata
        # Engine fraction diagnostics refer to selection+test rows; label explicitly.
        diagnostics[[length(diagnostics)+1L]]<-transform(local$diagnostics,fold=fold,
          fraction_scope='selection_plus_test; not used for selection criterion')
      }
    }
    m<-do.call(rbind,local_metadata);s<-do.call(cbind,local_signals)
    key<-paste(m$Lambda0,m$treatment,m$family,sep='::');ord<-order(key)
    if(fold==1L) reference<-key[ord] else stopifnot(identical(reference,key[ord]))
    signals[[length(signals)+1L]]<-data.frame(row_id=vi,s[,ord,drop=FALSE])
    metadata[[length(metadata)+1L]]<-m[ord,]
    splits[[length(splits)+1L]]<-list(fold=fold,development=di,selection=si,test=vi)
  }
  stacked<-do.call(rbind,signals);stacked<-stacked[order(stacked$row_id),]
  stopifnot(!anyDuplicated(stacked$row_id))
  m<-metadata[[1]];m$r<-'training_selected';m$h<-'training_selected'
  cal<-calibrate_learned_witness_maximum(as.matrix(stacked[,-1,drop=FALSE]),m,
    bootstrap_replications=bootstrap,confidence_level=.95,seed=seed+900000L)
  list(summary=cal$summary,moment_detail=cal$moment_detail,
    witness_diagnostics=do.call(rbind,diagnostics),gamma_diagnostics=gamma_diagnostics,
    signals=as.matrix(stacked[,-1,drop=FALSE]),test_row_ids=stacked$row_id,splits=splits,
    settings=list(variant=variant,gamma='riesz_spline',n_total=nrow(data),
      n_test=nrow(stacked),bootstrap=bootstrap,inner_split='50/50 development/selection',
      allocation_spline_df=5,selection='internal held-out one-step mean / SE'))
}
