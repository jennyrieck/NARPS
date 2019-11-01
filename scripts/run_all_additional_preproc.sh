basedir=$1
sub=$2

matlab -nodesktop -nosplash -r "batch_preproc_regress_covs(${basedir}, ${sub});"

for r in 1 2 3 4
do
	fslmaths ${basedir}/fmriprep/${sub}/func/${sub}_task-MGT_run-0${r}_bold_space-MNI152NLin2009cAsym_preproc_regress_cov.nii -kernel gauss 1.69851380042 -fmean ${basedir}/fmriprep/${sub}/func/${sub}_task-MGT_run-0${r}_bold_space-MNI152NLin2009cAsym_preproc_regress_cov_smooth4mm.nii
done


