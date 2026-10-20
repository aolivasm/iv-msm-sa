# Artificial deterministic schema fixtures, not paper results or new simulations.
release<-normalizePath('.')
fixture<-tempfile('paper-collector-test-');dir.create(fixture)
dir.create(file.path(fixture,'results'),recursive=TRUE)
ep<-expand.grid(outcome=c('Continuous','Binary'),setting=1:4,
 Lambda=c(1.8,2,2.5,3),display=c('IV-MSM','Tan','Dorn'),replication=1:500,
 stringsAsFactors=FALSE)
ep<-subset(ep,!(outcome=='Continuous' & setting==4 & display=='IV-MSM' & Lambda%in%c(1.8,2.5)))
ep$n<-6000L;ep$lower<--.2;ep$upper<-.5
ep$target_lower<--.1;ep$target_upper<-.4
bad<-with(ep,outcome=='Continuous' & setting==4 & display=='IV-MSM' & Lambda==2)
ep$target_lower[bad]<-ep$target_upper[bad]<-NA_real_
saveRDS(unique(ep[c('outcome','setting','Lambda','display','target_lower','target_upper')]),
 file.path(fixture,'results/population_targets.rds'))
save_fixture<-function(name,path,obj){
 f<-file.path(fixture,'experiments',name,path);dir.create(dirname(f),recursive=TRUE,showWarnings=FALSE)
 saveRDS(obj,f)
}
for(name in c('continuous_bounds','continuous_bounds_additional_values')) {
 L<-if(name=='continuous_bounds')c(2,3)else c(1.8,2.5)
 d<-subset(ep,outcome=='Continuous'&display=='IV-MSM'&Lambda%in%L)
 d$method<-'rkhs_tuned';d$class<-'fixed_spline5';d$success<-TRUE
 d$compatible<-is.finite(d$target_lower);d$outward<-ifelse(d$compatible,TRUE,NA)
 d$confidence_width<-d$upper-d$lower+.2
 d$target_lower<-d$target_upper<-NULL
 result<-if(name=='continuous_bounds')'x5_stress_v1' else 'x5_stress_extra_lambda_v1'
 save_fixture(name,paste0('results/',result,'/summary.rds'),list(detail=d))
}
d<-subset(ep,outcome=='Binary'&display=='IV-MSM')
d$method<-'rkhs';d$estimate_lower<-d$lower;d$estimate_upper<-d$upper
d$outward<-TRUE;d$confidence_width<-d$upper-d$lower+.2;d$target_lower<-d$target_upper<-NULL
save_fixture('binary_bounds','results/binary_x5_stress_v1/summary.rds',list(status=data.frame(success=TRUE),detail=d))
for(name in c('continuous_comparators','continuous_comparators_additional_values','binary_comparators')) {
 outcome<-if(name=='binary_comparators')'Binary'else 'Continuous'
 L<-switch(name,continuous_comparators=c(2,3),continuous_comparators_additional_values=c(1.8,2.5),c(1.8,2,2.5,3))
 d<-ep[as.character(ep$outcome)==outcome&ep$display!='IV-MSM'&ep$Lambda%in%L,]
 d$method<-ifelse(d$display=='Tan','tan','dorn');d$adjustment<-'x';d$success<-TRUE
 lo<-transform(d,endpoint='ate_lower',estimate=lower,confidence_lower=lower-.1,confidence_upper=lower+.1)
 hi<-transform(d,endpoint='ate_upper',estimate=upper,confidence_lower=upper-.1,confidence_upper=upper+.1)
 result<-switch(name,continuous_comparators='x5_stress_comparators_v1',continuous_comparators_additional_values='x5_stress_extra_comparators_v1',binary_comparators='binary_x5_comparators_v1')
 save_fixture(name,paste0('results/',result,'/summary.rds'),list(endpoints=rbind(lo,hi)))
}
k<-expand.grid(n=c(6000,12000),outcome=c('continuous','binary'),setting=0:4,Lambda0=c(1,1.1))
k$method<-'rkhs';k$successful<-500L
for(n in c(6000,12000))save_fixture(paste0('compatibility_n',n),paste0('results/x5_compatibility_n',n,'_v1/summary.rds'),list(status=data.frame(success=TRUE),summary=k[k$n==n,]))
save_fixture('binary_compatibility_n24000','results/binary_s4_compatibility_n24000_v1/summary.rds',list(summary=data.frame(successful=500L)))
run<-function(destination){
 old<-getwd();on.exit(setwd(old));setwd(fixture)
 system2(file.path(R.home('bin'),'Rscript'),c(shQuote(file.path(release,'scripts/collect_paper_results.R')),destination),stdout=TRUE,stderr=TRUE)
}
ok<-run('complete');if(!is.null(attr(ok,'status')))stop(paste(ok,collapse='\n'))
z<-readRDS(file.path(fixture,'complete/bounds.rds'));stopifnot(nrow(z$main)==93L,nrow(z$endpoint_estimates)==47000L)
f<-file.path(fixture,'experiments/binary_bounds/results/binary_x5_stress_v1/summary.rds')
z<-readRDS(f);z$status$success<-FALSE;saveRDS(z,f)
bad<-suppressWarnings(run('must_fail'));stopifnot(!is.null(attr(bad,'status')),!dir.exists(file.path(fixture,'must_fail')))
cat('PASS: main collection schema and retained-failure rejection; no estimators called.\n')
unlink(fixture,recursive=TRUE)
