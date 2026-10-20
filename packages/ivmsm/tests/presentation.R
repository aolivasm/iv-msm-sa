library(ivmsm)
wrap <- get('.ivmsm_result',asNamespace('ivmsm'))
fit <- wrap(list(Lambda0=2,Lambda1=3,epsilon=0,
 endpoint_intervals=data.frame(endpoint=c('ate_lower','ate_upper'),estimate=c(-.2,.4),standard_error=c(.1,.12)),
 ate_confidence_set=data.frame(confidence_level=.95,lower=-.4,upper=.64,row.names='ate_lower'),
 diagnostics=data.frame(treatment=0:1,compatibility_projection_rate=c(.08,0),post_smooth_theta_crossing_rate=c(.004,0))),
 'binary outcome, binary IV')
text <- capture.output(print(fit)); s <- summary(fit)
details <- capture.output(print(fit,diagnostics=TRUE))
summary_text <- capture.output(print(s))
stopifnot(rownames(fit$ate_confidence_set)=='ATE',inherits(s,'summary.ivmsm_fit'),
 any(grepl('Lambda_0 = 2',text)),any(grepl('Lambda_1 = 3',text)),
 !any(grepl('Numerical checks|Adjusted|safeguards',text)),
 !any(grepl('Numerical checks|Adjusted|safeguards',summary_text)),
 any(grepl('Adjusted fitted values',details)),!any(grepl('projection',text)),
 identical(s$endpoints,fit$endpoint_intervals),nrow(s$diagnostics)==4L)
for(type in c('binary','ordered','continuous'))for(J in c(2,3)) {
 d <- mock_iv_data(type,J)
 stopifnot(nrow(d)==120*J,identical(names(d),c('Y','T_obs','Z','X')),
  length(unique(d$Z))==J,!anyNA(d))
}
cat('Display, sensitivity labels, ATE row name and mock-data checks passed; no fits.\n')
band <- structure(list(Lambda_range=c(1.2,2),confidence_level=.95,
 table=data.frame(Lambda=seq(1.2,2,length.out=21),empty=FALSE),
 continuum=data.frame(Lambda=seq(1.2,2,length.out=201),lower=-.2,upper=.4,empty=FALSE)),
 class=c('ivmsm_band','list'))
btext <- capture.output(print(band))
stopifnot(any(grepl('21 sensitivity values',btext)),any(grepl('201 locations',btext)),
 !any(grepl('empty|Warning',btext)))
tmp <- tempfile(fileext='.pdf')
grDevices::pdf(tmp);plot(band);grDevices::dev.off()
stopifnot(file.info(tmp)$size>0);unlink(tmp)
band$continuum$empty[1] <- TRUE
stopifnot(any(grepl('Warning',capture.output(print(band)))))
cat('Band summary, plotting and crossing notification passed on deterministic fixtures.\n')
