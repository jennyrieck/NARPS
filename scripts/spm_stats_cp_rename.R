### copy and rename parametric first level con files for easy use

stats.dir<-'../../spm/first_level/'
out.dir<-'../../pls/spm_con_first_level/parametric_new_preproc/'
all.subs<-read.delim('../../event_tsvs_Nov16/participants.tsv', sep='\t',stringsAsFactors = F)

for(s in 1:dim(all.subs)[1]){
  file.copy(paste0(stats.dir,all.subs$participant_id[s],'/parametric/con_0001.nii'),
            paste0(out.dir,all.subs$group[s],'_',all.subs$participant_id[s],'_con_gain_increase.nii'))
  file.copy(paste0(stats.dir,all.subs$participant_id[s],'/parametric/con_0003.nii'),
            paste0(out.dir,all.subs$group[s],'_',all.subs$participant_id[s],'_con_loss_increase.nii')) 
}

