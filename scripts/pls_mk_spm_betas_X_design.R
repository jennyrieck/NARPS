## various changes made which need to be undone...

rm(list=ls())
library(oro.nifti)
library(TExPosition)
library(neuroim)
library(ours)

use.parametric <- T

pls.dir <- '../../pls/'

mask <- loadVolume('./masks/HarvardOxford-cort_subcort_cereb_binary_mask-thr0_2mm.nii')
mask.ind <- which(mask@.Data==1)

all.subs <- read.delim('../../event_tsvs_Nov16/participants.tsv', sep='\t',stringsAsFactors = F)

  first.level.dir<-'parametric_new_preproc'
  condition.labels<-c('gain_increase', 'loss_increase')
  n.cond<-length(condition.labels)

nii.dir<-paste0(pls.dir,'spm_con_first_level/',first.level.dir)

spm.betas<-matrix(0,dim(all.subs)[1]*n.cond,length(mask.ind))
aggregate.design <- matrix(NA,dim(all.subs)[1]*n.cond, 3)
colnames(aggregate.design) <- c("SUBJECT","GROUP","CONDITION")

for(s in 1:dim(all.subs)[1]){
  for(c in 1:n.cond){
    
    this.con <- loadVector(paste0(nii.dir,'/',all.subs$group[s],'_',all.subs$participant_id[s],'_con_',condition.labels[c],'.nii'),mask=mask)
    this.row <- (s-1)*n.cond + c
    spm.betas[this.row,] <- this.con@data
    
    aggregate.design[this.row,"SUBJECT"] <- all.subs$participant_id[s]
    aggregate.design[this.row,"GROUP"] <- all.subs$group[s]
    aggregate.design[this.row,"CONDITION"] <- condition.labels[c]
    
  }
  print(s)
}

rownames(aggregate.design) <- aggregate.design[,"SUBJECT"]

save(aggregate.design,file='../../rdata/aggregate.design_2mm.rda')
save(spm.betas,file='../../rdata/spm.betas_2mm.rda')