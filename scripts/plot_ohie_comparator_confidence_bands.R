#!/usr/bin/env Rscript
# Reads only completed, separately calibrated bands. Never interpolates missing fits.
root <- 'results/ohie_simultaneous_grid/official_comparators_dense_v1'
ivroot <- 'results/ohie_simultaneous_grid/exact_dense_lambda120_200_v1'
figdir <- 'output/figures'
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
stepcoords <- function(x,y) list(x=rep(x,each=2)[-length(rep(x,each=2))],y=rep(y,each=2)[-1])
stepdraw <- function(x,y,...) {v<-stepcoords(x,y);lines(v$x,v$y,...)}
colours <- c(RKHS='#2166AC',LogSpline='#7570B3',Tan='#D95F02',Dorn='#1B9E77')
labels <- c(RKHS='IV-MSM: RKHS',LogSpline='IV-MSM: interaction log-spline',Tan='Tan',Dorn='Dorn–Guo–Kallus')
readpanel <- function(adj) {
  paths <- c(RKHS=file.path(ivroot,'RKHS.rds'),LogSpline=file.path(ivroot,'LogSpline.rds'),
    Tan=file.path(root,paste0('Tan_',adj,'.rds')),Dorn=file.path(root,paste0('Dorn_',adj,'.rds')))
  lapply(paths[file.exists(paths)],readRDS)
}
panels <- lapply(c('x','xz'),readpanel)
complete <- all(vapply(panels,length,integer(1))==4)
suffix <- if(complete)'complete' else 'partial'
plotfile <- file.path(figdir,paste0('ohie_comparator_dense_bands_',suffix,'.png'))
png(plotfile,width=2400,height=1250,res=220)
par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),oma=c(4,0,0,0),mgp=c(2.7,.8,0),las=1,bty='l')
for(i in 1:2) {
  panel<-panels[[i]]
  plot(c(1.2,2),c(-.15,.45),type='n',xlab=expression(Lambda),ylab='ATE (probability difference)',
    main=c('Conventional MSM adjusting for X','Conventional MSM adjusting for (X, Z)')[i],xaxt='n')
  axis(1,at=seq(1.2,2,.2));abline(h=0,lwd=1.5,col='gray25');abline(v=1.6,lty=3,col='gray70')
  for(nm in names(panel)) {
    t<-panel[[nm]]$table
    for(side in c('lower','upper'))stepdraw(t$Lambda,t[[paste0('simultaneous_',side)]],
      col=colours[nm],lwd=2,lty=if(nm=='LogSpline')2 else 1)
  }
  legend('topleft',labels[names(panel)],col=colours[names(panel)],
    lty=ifelse(names(panel)=='LogSpline',2,1),lwd=2,bty='n',cex=.82)
}
mtext('95% simultaneous outward bands over 81 sensitivity values; each method calibrated separately.',outer=TRUE,side=1,line=1,cex=.85)
mtext(if(complete)'Steps use the right grid endpoint. IV-MSM curves are the same in both panels.' else
  'PARTIAL RESULTS: Tan dense-grid fits are still running; no Tan curves are shown.',outer=TRUE,side=1,line=2.5,cex=.85)
dev.off()
png(file.path(figdir,paste0('ohie_comparator_dense_widths_',suffix,'.png')),width=2400,height=1100,res=220)
par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),oma=c(2,0,0,0),mgp=c(2.7,.8,0),las=1,bty='l')
for(i in 1:2) {
  panel<-panels[[i]]
  plot(c(1.2,2),c(.07,.60),type='n',xlab=expression(Lambda),ylab='Simultaneous band width',
    main=c('Conventional MSM adjusting for X','Conventional MSM adjusting for (X, Z)')[i],xaxt='n')
  axis(1,at=seq(1.2,2,.2))
  for(nm in names(panel)) {
    t<-panel[[nm]]$table;stepdraw(t$Lambda,t$simultaneous_upper-t$simultaneous_lower,
      col=colours[nm],lwd=2,lty=if(nm=='LogSpline')2 else 1)
  }
  legend('topleft',labels[names(panel)],col=colours[names(panel)],lty=ifelse(names(panel)=='LogSpline',2,1),lwd=2,bty='n',cex=.82)
}
mtext(if(complete)'Separate 95% simultaneous coverage for each method.' else 'Partial comparison: Tan pending.',outer=TRUE,side=1,cex=.85)
dev.off()
summaries<-list()
for(i in 1:2)for(nm in names(panels[[i]])) {
  a<-panels[[i]][[nm]];t<-a$table[abs(a$table$Lambda-1.6)<1e-9,]
  summaries[[length(summaries)+1]]<-data.frame(method=nm,adjustment=c('x','xz')[i],
    critical=a$calibration$critical,Lambda=1.6,lower=t$simultaneous_lower,upper=t$simultaneous_upper,
    width=t$simultaneous_upper-t$simultaneous_lower)
}
capture.output(print(do.call(rbind,summaries),row.names=FALSE,digits=8),
  file=file.path(root,paste0('summary_',suffix,'.txt')))
message(plotfile)
