home_dir = '/project/3018051.01/ruggero/derivatives/first_level/task-rsa';
addpath('/project/3018051.01/ruggero/spm12');
subjects = {'sub-102', 'sub-201', 'sub-202'};
sessions = {'ses-mri01', 'ses-mri02'};  

for i = 1:length(subjects)
    subj = subjects{i};
    for j = 1:length(sessions)
        sess = sessions{j};
    
        clear matlabbatch

        stats_dir = sprintf('/project/3018051.01/ruggero/derivatives/first_level/task-rsa/%s/%s_SPM/SPM_stats/%s', subj, subj, sess);
        allCond_dir = sprintf('/project/3018051.01/ruggero/derivatives/first_level/task-rsa/%s/%s_SPM/%s_%s_all-conditions.mat', subj, subj, subj, sess);
        regres_dir = sprintf('/project/3018051.01/ruggero/derivatives/first_level/task-rsa/%s/%s_%s_task-rsa_run-concat_nuisance_regressors.txt', subj, subj, sess);


        matlabbatch{1}.spm.stats.fmri_spec.dir = {stats_dir};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.5;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

        % Construct the full path to the NIfTI folder
        scan_folder = fullfile(home_dir, ...
            subj, ...
            [subj '_SPM'], ...
            [subj '_input_spm'], ...
            sess);

        % Check and read NIfTI files
        scan_files = dir(fullfile(scan_folder, '*.nii'));

        % Build full paths + ',1'
        scan_paths = fullfile({scan_files.folder}, {scan_files.name});
        scan_paths = strcat(scan_paths, ',1');

        % ===== STEP 1: Model Specification =====
        fprintf('Running model specification for %s, %s\n', subj, sess);

        matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = scan_paths';
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {allCond_dir};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {regres_dir};
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
        % Run model specification and wait for it to complete
        spm_jobman('run', matlabbatch);

        % Check if SPM.mat was created
        if ~exist(fullfile(stats_dir, 'SPM.mat'), 'file')
            warning('SPM.mat not created for %s, %s. Skipping estimation and contrasts.', subj, sess);
            continue;
        end

        % ===== STEP 2: Model Estimation =====
        clear matlabbatch
        fprintf('Running model estimation for %s, %s\n', subj, sess);

        matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(stats_dir, 'SPM.mat')};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

        spm_jobman('run', matlabbatch);
        
        % ===== STEP 3: Contrast Manager =====
        clear matlabbatch

        matlabbatch{1}.spm.stats.con.spmmat = {fullfile(stats_dir, 'SPM.mat')};
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'All conditions';
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = [1 0 0 0
                                                                0 1 0 0
                                                                0 0 1 0
                                                                0 0 0 1];
        matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'All graphemes';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Colored graphemes';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 1 0 0];
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Trained graphemes';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 1 0];
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'oddball';
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 1];
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.delete = 0;
        
        fprintf('Running contrasts for %s, %s\n', subj, sess);
        
        spm_jobman('run', matlabbatch);

    end
end