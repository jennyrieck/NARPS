function run_spm_first_level(data_dir,sub)

%%% Run spm_fist_level_parametric.m
first_level_dir = 'first_level';

n_trs = 453;
n_runs=4;

preproc_prefix = '';
preproc_suffix = 'bold_space-MNI152NLin2009cAsym_preproc_regress_cov_smooth4mm';
mask = '../masks/HarvardOxford-cort_subcort_cereb_binary_mask-thr0_2mm.nii';
regressors=0;

spm_first_level_parametric(data_dir, data_dir, first_level_dir,sub,n_trs,n_runs,preproc_prefix, preproc_suffix,mask, regressors)
