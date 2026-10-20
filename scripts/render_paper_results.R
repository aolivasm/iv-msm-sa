#!/usr/bin/env Rscript
# Render aggregate results only. This script never calls an estimator.
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
args<-commandArgs(TRUE);input<-if(length(args))args[1] else 'results/paper_reproduction'
out<-if(length(args)>1)args[2] else 'output/paper'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
z<-readRDS(file.path(input,'bounds.rds'));main<-z$main;ep<-z$endpoint_estimates
stopifnot(nrow(main)==93L,all(main$R==500),nrow(ep)==47000L)
old<-list(application_palette=z$palette)
source('scripts/paper_plot_specification.R')
legend_theme<-theme(legend.position='bottom',legend.text=element_text(size=14),
 legend.key.width=grid::unit(1,'cm'))
for(nm in c('fig2','pA','pB')) {
 p<-get(nm)
 for(s in p$scales$scales)if(any(s$aesthetics%in%c('colour','fill','shape','linetype')))
  s$labels<-function(x)ifelse(x=='Dorn','Dorn-Guo-Kallus',x)
 assign(nm,p)
}
boxes<-fig2+labs(caption=NULL,y='ATE bound estimate')+legend_theme
coverage<-((pA+labs(title=NULL,tag='A'))/(pB+labs(title=NULL,tag='B')))+plot_layout(guides='collect')
coverage<-coverage&legend_theme
ggsave(file.path(out,'figure3_endpoint_estimates.pdf'),boxes,width=11.8,height=7.2)
ggsave(file.path(out,'figure4_coverage_width.pdf'),coverage,width=11.8,height=9.4)
d<-readRDS(file.path(input,'compatibility.rds'))
stopifnot(nrow(d)==40L,all(d$successful==500),all(d$rejected/d$successful==d$rejection))
d$outcome<-factor(d$outcome,c('Continuous','Binary'))
d$sample<-factor(d$n,c(6000,12000),c('6,000','12,000'))
d$design<-factor(d$setting,c(1:4,0),c('a = 0.55, b = 0.075','a = 1.10, b = 0.075',
 'a = 1.65, b = 0.075','a = 1.10, b = 0.15','a = 1.10, b = 0\n(control)'))
d$lambda<-factor(d$Lambda0)
p<-ggplot(d,aes(sample,rejection,group=lambda))+
 geom_hline(yintercept=c(.05,1),linetype='dotted',colour='grey40')+
 geom_line(aes(linetype=lambda),colour='#0072B2')+
 geom_errorbar(aes(ymin=mc_low,ymax=mc_high),width=.09,colour='#0072B2')+
 geom_point(aes(shape=lambda),size=2.7,colour='#0072B2')+facet_grid(outcome~design)+
 scale_linetype_manual(values=c('solid','dashed'),name=expression(Lambda))+
 scale_shape_manual(values=c(16,17),name=expression(Lambda))+
 scale_y_continuous(limits=c(0,1.025),breaks=c(0,.05,.25,.5,.75,1))+
 labs(x='Sample size',y='Rejection probability')+theme_bw(base_size=12)+
 theme(legend.position='bottom',panel.grid.minor=element_blank())
ggsave(file.path(out,'figure5_compatibility.pdf'),p,width=12,height=6.5)
write.csv(main,file.path(out,'numerical_results.csv'),row.names=FALSE)
write.csv(d,file.path(out,'compatibility_results.csv'),row.names=FALSE)

# The application input contains aggregate limits, never participant records.
f<-file.path(input,'ohie_bands.rds')
if(file.exists(f)) {
 source('R/binary_continuum_band.R');z<-readRDS(f)
 mesh<-seq(1.2,2,length.out=8001)
 bands<-lapply(z$bands,function(b)binary_continuum_band(b$grid$Lambda,
  b$grid$simultaneous_lower,b$grid$simultaneous_upper,mesh))
 stopifnot(length(bands)==3L,all(vapply(bands,function(b)!any(b$empty),logical(1))))
 cols<-c('#0072B2','#D55E00','#009E73');iv<-z$iv_only
 pdf(file.path(out,'figure6_ohie_bands.pdf'),width=10,height=4.8,useDingbats=FALSE)
 par(mfrow=c(1,2),mar=c(4.2,4.4,2,.8),oma=c(2,0,0,0),las=1,bty='l')
 plot(NA,xlim=c(1.2,2),ylim=c(-.28,.56),xlab=expression(Lambda),ylab='ATE (probability difference)',main='A. Simultaneous outward bands')
 abline(h=0,col='grey30');abline(h=c(iv$outward_lower,iv$outward_upper),col='#8C510A',lty=4)
 for(i in 1:3){lines(mesh,bands[[i]]$lower,col=cols[i],lty=i,lwd=1.8);lines(mesh,bands[[i]]$upper,col=cols[i],lty=i,lwd=1.8)}
 plot(NA,xlim=c(1.2,2),ylim=c(.08,.8),xlab=expression(Lambda),ylab='Width (probability difference)',main='B. Band widths')
 abline(h=iv$outward_upper-iv$outward_lower,col='#8C510A',lty=4)
 for(i in 1:3)lines(mesh,bands[[i]]$upper-bands[[i]]$lower,col=cols[i],lty=i,lwd=1.8)
 par(fig=c(0,1,0,1),mar=rep(0,4),oma=rep(0,4),new=TRUE);plot.new()
 legend('bottom',c(names(z$bands),'IV-only'),col=c(cols,'#8C510A'),lty=1:4,lwd=2,bty='n',horiz=TRUE)
 dev.off()
 report<-do.call(rbind,lapply(names(z$bands),function(name){b<-z$bands[[name]]
  data.frame(method=name,tipping=b$tipping,lower_1p6=b$at_1p6$lower,upper_1p6=b$at_1p6$upper)}))
 write.csv(report,file.path(out,'ohie_results.csv'),row.names=FALSE)
} else message('OHIE results not present; simulation figures only. Run collect_ohie_results.R to add them.')
for(name in c('binary_n24000.csv','ohie_compatibility.csv','ohie_aipw.csv')) {
 f<-file.path(input,name)
 if(file.exists(f))write.csv(read.csv(f),file.path(out,name),row.names=FALSE)
}
cat('Figures rendered from aggregate results; no simulations or fits run.\n')
