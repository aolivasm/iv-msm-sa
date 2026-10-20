#!/usr/bin/env Rscript
# Assemble completed main experiments. Never fits, retries, or omits failures.
args<-commandArgs(TRUE)
out<-if(length(args))args[1] else 'results/paper_reproduction'
if(dir.exists(out))stop('Choose a new output directory; refusing overwrite: ',out)
targets<-readRDS('results/population_targets.rds')
palette<-c('IV-MSM'='#0072B2',Tan='#D55E00',Dorn='#009E73')
read_experiment<-function(name,path) {
 f<-file.path('experiments',name,path)
 if(!file.exists(f))stop('Run and collect this analysis first: ',f)
 readRDS(f)
}
rows<-list()
add<-function(d,outcome,display) {
 d$outcome<-outcome;d$display<-display
 d<-merge(d,targets,
          by=c('outcome','setting','Lambda','display'),all.x=TRUE,sort=FALSE)
 stopifnot(all(is.finite(d$lower)),all(is.finite(d$upper)),all(d$n==6000))
 rows[[length(rows)+1L]]<<-d[,c('outcome','setting','Lambda','display','replication',
 'n','lower','upper','target_lower','target_upper','outward','confidence_width')]
}
for(name in c('continuous_bounds','continuous_bounds_additional_values')) {
 result<-if(name=='continuous_bounds')'x5_stress_v1' else 'x5_stress_extra_lambda_v1'
 d<-read_experiment(name,paste0('results/',result,'/summary.rds'))$detail
 d<-subset(d,method=='rkhs_tuned'&class=='fixed_spline5')
 stopifnot(all(d$success))
 # Only Lambda=2 is shown diagnostically in incompatible setting 4.
 d<-subset(d,compatible | (setting==4 & Lambda==2))
 add(d,'Continuous','IV-MSM')
}
z<-read_experiment('binary_bounds','results/binary_x5_stress_v1/summary.rds')
stopifnot(all(z$status$success));d<-subset(z$detail,method=='rkhs')
d$lower<-d$estimate_lower;d$upper<-d$estimate_upper
add(d,'Binary','IV-MSM')
for(name in c('continuous_comparators','continuous_comparators_additional_values','binary_comparators')) {
 result<-switch(name,continuous_comparators='x5_stress_comparators_v1',
 continuous_comparators_additional_values='x5_stress_extra_comparators_v1',
 binary_comparators='binary_x5_comparators_v1')
 z<-read_experiment(name,paste0('results/',result,'/summary.rds'))
 d<-subset(z$endpoints,adjustment=='x');stopifnot(all(d$success))
 keys<-c('setting','replication','method','Lambda','n')
 lo<-subset(d,endpoint=='ate_lower');hi<-subset(d,endpoint=='ate_upper')
 d<-merge(lo[,c(keys,'estimate','confidence_lower')],hi[,c(keys,'estimate','confidence_upper')],by=keys)
 names(d)[names(d)=='estimate.x']<-'lower';names(d)[names(d)=='estimate.y']<-'upper'
 d$outcome<-if(name=='binary_comparators')'Binary' else 'Continuous'
 d$display<-ifelse(d$method=='tan','Tan','Dorn')
 d<-merge(d,targets,by=c('outcome','setting','Lambda','display'),all.x=TRUE)
 stopifnot(all(is.finite(d$target_lower)),all(is.finite(d$target_upper)))
 d$outward<-d$confidence_lower<=d$target_lower & d$confidence_upper>=d$target_upper
 d$confidence_width<-d$confidence_upper-d$confidence_lower
 rows[[length(rows)+1L]]<-d[,names(rows[[1]])]
}
ep<-do.call(rbind,rows)
stopifnot(nrow(ep)==47000L,!anyDuplicated(ep[c('outcome','setting','Lambda','display','replication')]))
groups<-split(ep,interaction(ep[c('outcome','setting','Lambda','display')],drop=TRUE))
main<-do.call(rbind,lapply(groups,function(d) {
 stopifnot(nrow(d)==500L,setequal(d$replication,1:500))
 if(all(is.na(d$target_lower)))return(NULL)
 stopifnot(!anyNA(d$outward),all(is.finite(d$confidence_width)))
 ci<-binom.test(sum(d$outward),500)$conf.int
 data.frame(d[1,c('outcome','setting','Lambda','display')],R=500L,
 outward=mean(d$outward),confidence_width=mean(d$confidence_width),mc_lower=ci[1],mc_upper=ci[2])
}))
stopifnot(nrow(main)==93L)
compat<-do.call(rbind,lapply(c(6000,12000),function(n) {
 z<-read_experiment(paste0('compatibility_n',n),paste0('results/x5_compatibility_n',n,'_v1/summary.rds'))
 stopifnot(all(z$status$success));d<-subset(z$summary,method=='rkhs')
 d$outcome<-ifelse(d$outcome=='continuous','Continuous','Binary');d
}))
stopifnot(nrow(compat)==40L,all(compat$successful==500L))
n24<-read_experiment('binary_compatibility_n24000','results/binary_s4_compatibility_n24000_v1/summary.rds')$summary
stopifnot(n24$successful==500L)
dir.create(out,recursive=TRUE)
saveRDS(list(main=main,endpoint_estimates=ep,palette=palette),file.path(out,'bounds.rds'))
saveRDS(compat,file.path(out,'compatibility.rds'))
write.csv(n24,file.path(out,'binary_n24000.csv'),row.names=FALSE)
cat('Complete simulation results assembled without exclusions: ',out,'\n')
