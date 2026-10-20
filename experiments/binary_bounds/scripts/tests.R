source('scripts/common.R')
stopifnot(nrow(tasks)==2000,!anyDuplicated(tasks$seed),all(table(tasks$setting)==500),
 identical(config$Lambda,c(1.8,2,2.5,3)),identical(generate(tasks[1,],50),generate(tasks[1,],50)))
pop<-readRDS('targets/population.rds');stopifnot(pop$max_refinement<1e-4,nrow(pop$targets)==16,
 all(pop$targets$lower<=pop$targets$upper),identical(config,pop$config))
front<-readRDS('targets/frontier_audit.rds')$summary
stopifnot(all(front$frontier<min(config$Lambda)),identical(front$a,config$settings$a),identical(front$b,config$settings$b))
root<-tempfile('binary_checkpoint_unit_');dir.create(root)
t<-tasks[1,,drop=FALSE];sig<-list(test=TRUE)
z<-list(task=t,signature=sig,records=setNames(rep(list(list(success=FALSE,error='test failure')),2),config$methods))
p<-file.path(root,'task_0001.rds');atomic(z,p);h<-tools::md5sum(p)
run_task(1,root,sig,n=600L);stopifnot(identical(h,tools::md5sum(p)))
stopifnot(inherits(try(check(z,t,list(wrong=TRUE)),silent=TRUE),'try-error'))
unlink(root,recursive=TRUE)
cat('UNIT TESTS PASSED: task mapping, unique seeds, targets, numerical compatibility audit, atomic checkpoint/resume.\n')
