source('scripts/common.R')
source('scripts/adapters.R')
cfg<-n12_config();tasks<-n12_tasks()
stopifnot(nrow(tasks)==500L,identical(tasks$id,1:500),all(tasks$n==24000L),
 all(tasks$outcome=='binary'),all(tasks$setting==4L),all(tasks$a==1.1),all(tasks$b==.15),
 all(tasks$seed==261004L+100L*(tasks$replication-1L)),identical(cfg$Lambda,1),identical(cfg$methods,'rkhs'))
n12_manifest()
contract<-n12_adapter_contract('binary')
stopifnot(contract$rkhs$cv_repeats==2L,contract$rkhs$cv_folds==3L,
 identical(contract$rkhs$bandwidth_multipliers,1/sqrt(2)),contract$bootstrap==399L)
# Compare generator output independently with the unchanged frozen parent.
a<-n12_adapter_new(getwd(),'binary',1,24000L)
p<-n12_adapter_load_parent(normalizePath('frozen/binary'))
x<-a$generate(tasks[1,,drop=FALSE],120L)
y<-n12_adapter_in_directory(normalizePath('frozen/binary'),function()p$generate(tasks[1,,drop=FALSE],120L))
stopifnot(identical(x,y))
# Mock fitted output only: exercise single-Lambda validation, resume and failures.
factory<-function(...) list(generate=function(task,n)task,run=function(data,method,seed) {
 stopifnot(method=='rkhs')
 list(summary=data.frame(Lambda0=1,Lambda1=1,statistic=1,critical_value=2,p_value=.5,reject=FALSE,moments_tested=8L),
 gamma_selection=rep(list(list()),10),paired=list(data=n12_fingerprint(data),
 folds=n12_fingerprint(1:5),primitives=setNames(rep(n12_fingerprint('p'),5),as.character(seed+10000L*(1:5)))),
 fold_sizes=rep(4800L,5),projection_applied=FALSE,q_fits=0L)
})
tmp<-tempfile('n24-test-');dir.create(tmp);sig<-list(test=TRUE)
z<-n12_run_task(1L,sig,result_root=tmp,adapter_factory=factory)
f<-file.path(tmp,'runs/task_0001.rds');hash<-tools::md5sum(f)
stopifnot(z$records$rkhs$success)
n12_run_task(1L,sig,result_root=tmp,adapter_factory=function(...)stop('must not refit'))
stopifnot(identical(hash,tools::md5sum(f)))
expect_error<-function(expr)stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),'error'))
expect_error(n12_load_task(f,tasks[2,,drop=FALSE],sig))
expect_error(n12_atomic(list(bad=TRUE),f,writer=function(object,path)stop('write failure')))
stopifnot(identical(hash,tools::md5sum(f)))
bad<-n12_run_task(2L,sig,result_root=tmp,adapter_factory=function(...)stop('mock failure'))
stopifnot(!bad$records$rkhs$success)
n12_run_task(2L,sig,result_root=tmp,adapter_factory=function(...)stop('must not retry'))
cat('TESTS PASSED: 500-task scope, frozen binary generator, tuning contract, single-Lambda records, unchanged resume, foreign-record rejection, failed-write protection, retained failures.\n')
