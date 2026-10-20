# GPFS quota, not filesystem free capacity. Fail closed on unknown output.
cmd<-'/usr/lpp/mmfs/bin/mmlsquota'
if(!file.exists(cmd))stop('GPFS quota tool unavailable; verify account quota before submission.')
raw<-system2(cmd,c('-u',Sys.getenv('USER'),'--block-size','auto','mmfs1'),stdout=TRUE,stderr=TRUE)
cat(paste(raw,collapse='\n'),'\n')
if(!is.null(attr(raw,'status')))stop('Quota check failed.')
lines<-grep('\\bUSR\\b',raw,value=TRUE)
if(length(lines)!=1)stop('Expected exactly one user quota row; inspect manually.')
fields<-strsplit(trimws(lines),'[[:space:]]+')[[1]];j<-match('USR',fields)
bytes<-function(x) {
 if(!grepl('^[0-9.]+[KMGTPE]?$' ,x))stop('Unrecognized quota value: ',x)
 unit<-sub('^[0-9.]+','',x);n<-as.numeric(sub('[KMGTPE]$','',x))
 n*1024^(match(unit,c('','K','M','G','T','P','E'))-1)
}
used<-bytes(fields[j+1]);soft<-bytes(fields[j+2]);hard<-bytes(fields[j+3]);doubt<-bytes(fields[j+4])
limits<-c(soft,hard);limits<-limits[limits>0]
if(!length(limits))stop('Cannot determine finite account quota.')
headroom<-(min(limits)-used-doubt)/1024^3
cat('Conservative headroom after in-doubt allocation:',headroom,'GiB\n')
if(headroom<8)stop('Require at least 8 GiB headroom before this compact study.')
