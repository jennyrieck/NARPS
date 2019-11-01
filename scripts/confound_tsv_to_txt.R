narps.dir<-'./fmriprep/'

all.subjs<-grep('sub*',list.dirs(narps.dir,full.names=F,recursive=F),value = T)
runs<-c('run-01','run-02','run-03','run-04')

#confounds.to.keep<-c('X', 'Y', 'Z','RotX', 'RotY', 'RotZ','aCompCor00','aCompCor01','aCompCor02', 'aCompCor03','aCompCor04','aCompCor05', 'FramewiseDisplacement')
mpe<-c('X', 'Y', 'Z','RotX', 'RotY', 'RotZ')
acompcor<-c('aCompCor00','aCompCor01','aCompCor02', 'aCompCor03','aCompCor04','aCompCor05')

for(s in all.subjs){
  for(r in runs){
    
    these.confounds<-read.table(paste0(narps.dir,s,'/func/',
                                       s, '_task-MGT_', r, '_bold_confounds.tsv'),header=T,sep = '\t',na.strings = 'n/a')
    mpe.out<-these.confounds[,mpe]
    mpe.out[is.na(mpe.out)]<-0
    write.table(mpe.out,file = paste0(narps.dir,s,'/func/',s, '_', r, '_motion_parameters.txt'),sep='\t',row.names = F, col.names = F)
    
    acompcor.out<-these.confounds[,acompcor]
    acompcor.out[is.na(acompcor.out)]<-0
    write.table(acompcor.out,file = paste0(narps.dir,s,'/func/',s, '_', r, '_aCompCor.txt'),sep='\t',row.names = F, col.names = F)
    
    all.cov.out<-cbind(mpe.out,acompcor.out)
    write.table(acompcor.out,file = paste0(narps.dir,s,'/func/',s, '_', r, '_all_regressors.txt'),sep='\t',row.names = F, col.names = F)
    
  }
}