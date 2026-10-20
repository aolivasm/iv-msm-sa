# Shared engines, isolated generated workspaces. Preparation never fits models.
paper_settings <- function(root=getwd())
  dget(file.path(root,'configurations/main_experiments.R'))
paper_workspace <- function(key,root=getwd())
  file.path(root,'results/work',key)
prepare_experiments <- function(root=getwd(),keys=NULL) {
  root<-normalizePath(root,mustWork=TRUE);settings<-paper_settings(root)
  if(is.null(keys))keys<-names(settings$profiles)
  stopifnot(all(keys%in%names(settings$profiles)))
  for(key in keys) {
    profile<-settings$profiles[[key]]
    src<-file.path(root,'experiments',profile$engine);dst<-paper_workspace(key,root)
    files<-sort(list.files(src,recursive=TRUE,full.names=TRUE));files<-files[!dir.exists(files)]
    fingerprints<-tools::md5sum(c(files,file.path(root,'configurations/main_experiments.R')))
    names(fingerprints)<-c(substring(files,nchar(src)+2L),'configuration')
    seal<-file.path(dst,'source_signature.rds')
    if(dir.exists(dst)) {
      if(!file.exists(seal)||!identical(readRDS(seal),fingerprints))
        stop('Incomplete workspace or changed source/configuration: ',dst,
             '. Preserve it and use a fresh checkout; do not mix checkpoints.')
      paths<-file.path(dst,names(fingerprints)[-length(fingerprints)])
      stopifnot(identical(unname(tools::md5sum(paths)),unname(fingerprints[-length(fingerprints)])))
      e<-new.env(parent=baseenv())
      sys.source(file.path(dst,'paper_config.R'),envir=e)
      stopifnot(identical(e$paper_profile,profile))
      if(profile$engine=='compatibility') {
        stopifnot(identical(readRDS(file.path(dst,'study_config.rds')),profile$config),
          identical(readRDS(file.path(dst,'frozen/binary/study_config.rds')),settings$frozen$binary))
      }
      next
    }
    dir.create(dst,recursive=TRUE)
    for(i in seq_along(files)) {
      target<-file.path(dst,substring(files[i],nchar(src)+2L))
      dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
      stopifnot(file.copy(files[i],target))
    }
    writeLines(c('paper_profile <-',capture.output(dput(profile))),file.path(dst,'paper_config.R'))
    if(profile$engine=='compatibility') {
      saveRDS(profile$config,file.path(dst,'study_config.rds'))
      saveRDS(settings$frozen$binary,file.path(dst,'frozen/binary/study_config.rds'))
    }
    bundles<-unique(dirname(dirname(list.files(dst,pattern='^run[.]R$',recursive=TRUE,full.names=TRUE))))
    for(bundle in bundles[order(nchar(bundles),decreasing=TRUE)]) {
      paths<-sort(list.files(bundle,pattern='[.]R$|^study_config[.]rds$',recursive=TRUE,full.names=TRUE))
      hashes<-tools::md5sum(paths);names(hashes)<-substring(paths,nchar(bundle)+2L)
      saveRDS(hashes,file.path(bundle,'manifest.rds'))
    }
    saveRDS(fingerprints,seal)
  }
  invisible(TRUE)
}
if(sys.nframe()==0L) {
  prepare_experiments()
  cat('Configured local workspaces prepared; no fits run.\n')
}
