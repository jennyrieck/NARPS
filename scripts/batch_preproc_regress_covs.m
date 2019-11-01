function batch_preproc_regress_covs(basedir,sub)

nRun = 1;
numIgnore = 0;
rmMean = 0;
detrendNum = 3;
wmmask = [];
vesmask = [];
csfmask = [];

for r = 1:4
    fn = sprintf('%s/fmriprep/%s/func/%s_task-MGT_run-0%d_bold_space-MNI152NLin2009cAsym_preproc.nii', basedir, sub, sub, r);
    fn_out = sprintf('%s/fmriprep/%s/func/%s_task-MGT_run-0%d_bold_space-MNI152NLin2009cAsym_preproc_regress_cov.nii', basedir, sub, sub, r);
    mpe = sprintf('%s/fmriprep/%s/func/%s_run-0%d_all_regressors.txt', basedir, sub, sub, r);
 
	if exist(fn, 'file')==2
		disp(fn)
	end
	if exist(mpe, 'file')==2   
		disp(mpe)
	end

	NoiseRegressingOut(fn,fn_out,nRun,numIgnore,rmMean,detrendNum,wmmask,csfmask,vesmask,mpe,0);
    
end

%%NoiseRegressingOut(Filename,Fileout,NRun,NumScantoIgnore,RemoveMean,TrendOrder,WMMaskName,CSFMaskName,VesselMaskName,MotionFile,Normalize_To_Std)
