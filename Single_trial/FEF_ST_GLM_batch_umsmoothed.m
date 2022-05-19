function FEF_ST_GLM_batch_umsmoothed(path, ID)
    %% Parameters
    nRun = 4;
    fwhm = 8;

    %% Loop over subJects
    for xSub = 1 : length(ID)
        subDir = fullfile(path, ID{xSub});
        % Images
        for xRun = 1 : nRun
            runDir  = fullfile(subDir, sprintf('Run%d', xRun));
            EPI_data_bsm{xRun, 1} = cellstr(spm_select('FPList', runDir, '^w2w.*.\.img'));
            Headmotion{xRun, 1}   = cellstr(spm_select('FPList', runDir, '^Headmotion.mat$'));               
        end
        
        % Onsets
        load(fullfile(subDir, 'Single_trial', 'TrialType.mat'));
        d_temp = dir(fullfile(subDir, 'Single_trial'));
        nTrial = length(Trial_Num_valid);
        for i = 1 : length(d_temp)
            if strcmp(d_temp(i).name, 'TrialType.mat')
                x = i+1;
            end
        end
        for xT = 1 : nTrial
            TrialDir{xT, 1} = fullfile(subDir, 'Single_trial', d_temp(x).name); 
            x = x+1;
            for xRun = 1 : nRun
                Onsets{xT, 1}{xRun,1} = cellstr(fullfile(TrialDir{xT, 1}, 'Onsets', ...
                            sprintf('Run%dOnsets.mat',xRun)));
            end
        end
        assignin('base', 'Headmotion', Headmotion)
        %% GLM loop over trials
       gcp;
       parfor xT = 1 : nTrial
            SPMDir = fullfile(TrialDir{xT, 1}, 'GLM');
            if ~exist(SPMDir, 'dir'); mkdir(SPMDir); end
            if ~exist(fullfile(SPMDir, 'con_0001.nii'), 'file')
                CV = Contrast_vector{xT,1};            
                FEF_ST_GLM_1st(SPMDir, nRun, EPI_data_bsm, Onsets{xT, 1}, Headmotion, CV);
            end
        end
   end
end

function FEF_ST_GLM_1st(SPMDir, nRun, smoothed_Data, Onsets, Headmotion, Contrast_vector)
    if exist(fullfile(SPMDir, 'SPM.mat'), 'file'); delete(fullfile(SPMDir, 'SPM.mat')); end
    %% Initialise SPM
    clear matlabbatch
    %% Specify 1st-level
    matlabbatch{1}.spm.stats.fmri_spec.dir = {SPMDir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    for xRun = 1 : nRun
        matlabbatch{1}.spm.stats.fmri_spec.sess(xRun).scans = smoothed_Data{xRun, 1};
        matlabbatch{1}.spm.stats.fmri_spec.sess(xRun).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(xRun).multi = Onsets{xRun, 1};
        matlabbatch{1}.spm.stats.fmri_spec.sess(xRun).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(xRun).multi_reg = Headmotion{xRun, 1};
        matlabbatch{1}.spm.stats.fmri_spec.sess(xRun).hpf = 256;
    end
    
    %% Model estimation
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    %% Contrast
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Current trial';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = Contrast_vector;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 1;
    
    %% Delete files
    matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Model estimation: Beta Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','beta'));
    matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('Model estimation: Analysis Mask', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','mask'));
    matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('Model estimation: ResMS Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','resms'));
    matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_move.files(4) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_move.files(5) = cfg_dep('Contrast Manager: All Stats Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spm'));
    matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
    %% Run batch
%     save('matlabbatch.mat', 'matlabbatch');
    spm('defaults', 'FMRI');
    spm_jobman('initcfg')
    spm_jobman('run', matlabbatch);
end