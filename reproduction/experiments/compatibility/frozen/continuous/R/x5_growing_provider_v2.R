source('R/x5_growing_provider.R')
source('R/x5_growing_gamma_v2.R')
x5g2_provider<-x5g_provider
body(x5g2_provider)<-do.call(substitute,list(body(x5g_provider),list(x5g_fit=as.name('x5g2_fit'))))
