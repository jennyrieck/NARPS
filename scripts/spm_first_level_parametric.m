function spm_first_level_parametric(narps_dir, write_dir, first_level_dir,sub,n_trs,n_runs,preproc_prefix, preproc_suffix,mask, regressors)
% created Jenny Rieck 7 Feb 2019 for use with #NARPs

% data_dir          :   base NARPs directory where 'fmriprep',
%                           'event_tsvs_Nov16', and 'spm' directories live
% write_dir         :   base directory to write out first level stats
%                           if unsure, just set write_dir = narps_dir
% first_level_dir   :   name for first_level_directory
% n_trs             :   number of TRs (should be consistent a across runs)
%                           for NARPS n_trs = 453
% n_runs            :   number of runs (should be consistent across people)
%                           for NARPs n_runs = 4
% sub               :   subject id (one person at a time)
% preproc_prefix    :   prefix (if any) for fmriprep preprocessed data
% preproc_suffix    :   suffix (if any) for fmriprep preprocessed data,
%                           do not include file extension (.nii, .nii.gz)
% mask              :   file name for first level mask to be applied
% regressors        :   1/0: include nuisance covariates (mpe, acompcor) in first level model?
%                           See coufound_tsv_to_txt.R to generate


if ~isdir(fullfile(write_dir, 'spm',first_level_dir))
    mkdir(fullfile(write_dir, 'spm',first_level_dir) )
end

if ~isdir(fullfile(write_dir, 'spm', first_level_dir, sub))
    mkdir(fullfile(write_dir, 'spm', first_level_dir, sub))
    mkdir(fullfile(write_dir, 'spm', first_level_dir, sub, 'parametric'));
end


matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(write_dir, 'spm',first_level_dir,sub, 'parametric')};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

for r = 1:n_runs
    soa = tdfread(fullfile(narps_dir,'event_tsvs_Nov16',sprintf('%s_task-MGT_run-0%d_events.tsv', sub, r)));
    for tr = 1:n_trs
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans{tr,1} = fullfile(narps_dir,'fmriprep',sub,'func',sprintf('%s%s_task-MGT_run-0%d_%s.nii,%d',preproc_prefix,sub,r,preproc_suffix,tr));
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.name = 'trial';
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.onset =  soa.onset;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.duration = soa.duration;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(1).name = 'gain';
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(1).param = soa.gain;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(1).poly = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(2).name = 'loss';
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(2).param = soa.loss;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(2).poly = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.orth = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
    if(regressors==1)
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {fullfile(narps_dir,'fmriprep',sub, 'func',sprintf('%s_run-0%d_spm_regressors.txt', sub, r))};
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
end

matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
matlabbatch{1}.spm.stats.fmri_spec.mask = {mask};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'gain increase';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0 1 0];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'gain decrease';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 -1 0];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'repl';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'loss increase';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'repl';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'loss decrease';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 0 -1];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'repl';
matlabbatch{3}.spm.stats.con.delete = 1;

save(fullfile(write_dir,'spm',first_level_dir,sprintf('%s_parametric.mat', sub)), 'matlabbatch');
spm_jobman('run',matlabbatch)
clear matlabbatch
