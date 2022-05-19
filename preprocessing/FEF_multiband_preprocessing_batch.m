function FEF_multiband_preprocessing_batch(path, ID, startSub, endSub)
    %% Parameters
    nRun = 4;
    nMotion = 6;
    fwhm = 8;
    
    %% Fieldmap Parameters
    FiledMap_TE = [4.92, 7.38];    % in ms
    echo_spacing = 0.65;       % in ms
    n_slices = 72;
    multiband_accel_factor = 6;
%     total_readout_time = (n_slices/multiband_accel_factor - 1) * echo_spacing;
    total_readout_time = (n_slices - 1) * echo_spacing;
       
    for xSub = startSub : endSub
        subDir = fullfile(path, ID{xSub},'scans');
       %% Load data    
       for xRun = 1 : nRun
            runDir = fullfile(subDir, sprintf('Run%d', xRun));
            runDir_SBRef = fullfile(subDir, sprintf('Run%d_SBRef', xRun));
            EPI_data{xSub,1}{xRun, 1} = cellstr(spm_select('FPList', runDir, '^*.\.img$'));
            SBRef_data{xSub,1}{xRun, 1} = cellstr(spm_select('FPList', runDir_SBRef, '^*.\.img$'));
       end
       T1Dir = fullfile(subDir, 'T1');
       t = cellstr(spm_select('FPList', T1Dir, '^*.\.img$'));
       T1_data{xSub,1} = t(1);     clear t;
       FiledMap_m_Dir = fullfile(subDir, 'FieldMap_magnitude'); 
       t = cellstr(spm_select('FPList', FiledMap_m_Dir, '^*.\.img$'));
       FM_m_data{xSub,1} = t(1);   clear t;
       FiledMap_pd_Dir = fullfile(subDir, 'FieldMap_phase_diff'); 
       t = cellstr(spm_select('FPList', FiledMap_pd_Dir, '^*.\.img$'));
       FM_pd_data{xSub,1} = t(1);   clear t;                      
    end

    gcp;
    parfor xSub = startSub : endSub
        subDir = fullfile(path, ID{xSub},'scans');
        my_multiband_preprocessing_main(subDir, nRun, fwhm, EPI_data{xSub,1}, SBRef_data{xSub,1}, ...
             FM_m_data{xSub,1}, FM_pd_data{xSub,1}, T1_data{xSub,1}, FiledMap_TE, total_readout_time);
    end
%     
    % Plot head motion
    for xSub = startSub : endSub
        subDir = fullfile(path, ID{xSub},'scans');
        plot_head_motion(nRun, subDir, nMotion)
    end
end

function my_multiband_preprocessing_main(subDir, nRun, fwhm, EPI_data, SBRef_data,  ...
         FM_m_data, FM_pd_data, T1_data, FiledMap_TE, total_readout_time)
        %% Preprocessing
%          cd(fullfile(subDir));
         % Realign and fieldmap correction
%          my_mutilband_preprocessing_1(nRun, EPI_data, SBRef_data,  ...
%              FM_m_data, FM_pd_data, FiledMap_TE, total_readout_time);
         my_fieldmap_correction(nRun, EPI_data, SBRef_data, ...
                 FM_m_data, FM_pd_data, FiledMap_TE, total_readout_time); 
             
        % Load filedmap corrected EPI and SBRef data
         for xRun = 1 : nRun
             runDir = fullfile(subDir, sprintf('Run%d', xRun));
             runDir_SBRef = fullfile(subDir, sprintf('Run%d_SBRef', xRun));
             EPI_data_f{xRun, 1} = cellstr(spm_select('FPList', runDir, '^f.*.\.img$'));
%              SBRef_data_f{xRun, 1} = cellstr(spm_select('FPList', runDir_SBRef, '^f.*.\.img$'));
             SBRef_data_f{xRun, 1} = cellstr(spm_select('FPList', runDir_SBRef, '^uu.*.\.img$'));
         end
         %  registe SBRef to T1, normalize coregisted EPI to MNI
          my_mutilband_preprocessing_2(nRun, EPI_data_f, SBRef_data_f,T1_data);
          
         % Normalize to both T1 template and T1 template within brainstem mask
         for xRun = 1 : nRun
             runDir = fullfile(subDir, sprintf('Run%d', xRun));
             runDir_SBRef = fullfile(subDir, sprintf('Run%d_SBRef', xRun));
             EPI_data_W{xRun, 1} = cellstr(spm_select('FPList', runDir, '^wfu.*.\.img$'));
             SBRef_data_W{xRun, 1} = cellstr(spm_select('FPList', runDir_SBRef, '^wuu.*.\.img$'));
         end
           T1Dir = fullfile(subDir, 'T1');
           t = cellstr(spm_select('FPList', T1Dir, '^w.*.\.img$'));
           T1_data_W = t(1);     clear t;
         normalize_add_brainstem_3(nRun, EPI_data_W, SBRef_data_W,T1_data_W);
         
         % Smoothing
         for xRun = 1 : nRun
             runDir = fullfile(subDir, sprintf('Run%d', xRun));
             runDir_SBRef = fullfile(subDir, sprintf('Run%d_SBRef', xRun));
             EPI_data_w{xRun, 1} = cellstr(spm_select('FPList', runDir, '^w.*.\.img$'));
             SBRef_data_w{xRun, 1} = cellstr(spm_select('FPList', runDir_SBRef, '^w.*.\.img$'));
         end
         my_smoothing(fwhm, nRun, EPI_data_w, SBRef_data_w);
end

function plot_head_motion(nRun, subDir, nMotion)
    R_all = [];
    for xRun = 1 : nRun
        RunDir = fullfile(subDir, sprintf('Run%d_SBRef', xRun));
        motion_file = spm_select('FPList', RunDir, '^*.\.txt$');
        t = textread(motion_file, '%f');
        nVol = length(t)/nMotion;
        for xVol = 1 :(nVol - 1)
            for i = 1 : nMotion
                R(xVol, i) = t(xVol * nMotion + i);
            end
        end
        R_all = [R_all; R];
        save(fullfile(RunDir, 'Headmotion.mat'), 'R');
        clear t R
    end
    %% Plot
    close all; figure; hold on
    subplot(2,1,1);
    plot(R_all(:,1),'b-'); hold on;
    plot(R_all(:,2),'g-'); hold on;
    plot(R_all(:,3),'r-'); hold on;
    xlabel('image'); ylabel('mm'); 
    legend('x translation','y translation', 'z translation',...
           'Orientation','horizontal','Location','SouthOutside');
    title('translation'); box('off');
    
    subplot(2,1,2);
    plot(R_all(:,4),'b-'); hold on;
    plot(R_all(:,5),'g-'); hold on;
    plot(R_all(:,6),'r-'); hold on;
    xlabel('image'); ylabel('degrees'); 
    legend('pitch','roll', 'yaw',...
           'Orientation','horizontal','Location','SouthOutside');
    title('rotation'); box('off');
    
    print('-dpdf',fullfile(subDir, 'Headmotion.pdf'));
end

function my_mutilband_preprocessing_1(nRun, EPI_data, SBRef_data, ...
             FM_m_data, FM_pd_data, FiledMap_TE, total_readout_time)
    clear matlabbatch
    
    %% Realign
    n = 1;
    for xRun = 1 : nRun
        matlabbatch{n}.spm.spatial.realignunwarp.data(1).scans = [SBRef_data{xRun,1}; EPI_data{xRun,1}];
        matlabbatch{n}.spm.spatial.realignunwarp.data(1).pmscan = '';
        matlabbatch{n}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
        matlabbatch{n}.spm.spatial.realignunwarp.eoptions.sep = 4;
        matlabbatch{n}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
        matlabbatch{n}.spm.spatial.realignunwarp.eoptions.rtm = 0;
        matlabbatch{n}.spm.spatial.realignunwarp.eoptions.einterp = 4;
        matlabbatch{n}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
        matlabbatch{n}.spm.spatial.realignunwarp.eoptions.weight = '';
        matlabbatch{n}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
        matlabbatch{n}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{n}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
        matlabbatch{n}.spm.spatial.realignunwarp.uweoptions.jm = 0;
        matlabbatch{n}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
        matlabbatch{n}.spm.spatial.realignunwarp.uweoptions.sot = [];
        matlabbatch{n}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
        matlabbatch{n}.spm.spatial.realignunwarp.uweoptions.rem = 1;
        matlabbatch{n}.spm.spatial.realignunwarp.uweoptions.noi = 5;
        matlabbatch{n}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
        matlabbatch{n}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
        matlabbatch{n}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
        matlabbatch{n}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
        matlabbatch{n}.spm.spatial.realignunwarp.uwroptions.mask = 1;
        matlabbatch{n}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
        n = n + 1;
    end
    
    %% FiledMap correction
    % Calculate VDM
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = FM_pd_data;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = FM_m_data;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = FiledMap_TE;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = total_readout_time;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 1;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = ...
        {'/Volumes/Data/matlab_tools/spm12/toolbox/FieldMap/T1.nii'};
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
    for xRun = 1 : nRun
        matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.session(xRun).epi = SBRef_data{xRun,1};
    end
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
    n = n + 1;
    
    % Apply VDM
    for xRun = 1 : nRun    
        matlabbatch{n}.spm.tools.fieldmap.applyvdm.data(xRun).scans(1) = ...
            cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', ...
            substruct('.','val', '{}',{xRun}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','sess', '()',{1}, '.','uwrfiles'));
        matlabbatch{n}.spm.tools.fieldmap.applyvdm.data(xRun).vdmfile(1) = ...
            cfg_dep(sprintf('Calculate VDM: Voxel displacement map (Subj 1, Session %d)', xRun), ...
            substruct('.','val', '{}',{nRun + 1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('()',{1}, '.','vdmfile', '{}',{xRun}));
    end
    matlabbatch{n}.spm.tools.fieldmap.applyvdm.roptions.pedir = 2;
    matlabbatch{n}.spm.tools.fieldmap.applyvdm.roptions.which = [2 1];
    matlabbatch{n}.spm.tools.fieldmap.applyvdm.roptions.rinterp = 4;
    matlabbatch{n}.spm.tools.fieldmap.applyvdm.roptions.wrap = [0 0 0];
    matlabbatch{n}.spm.tools.fieldmap.applyvdm.roptions.mask = 1;
    matlabbatch{n}.spm.tools.fieldmap.applyvdm.roptions.prefix = 'f';
    
    %% Run matlab batch
    spm('defaults', 'FMRI');
    spm_jobman('initcfg')
    spm_jobman('run', matlabbatch);
%  	save('matlabbatch.mat', 'matlabbatch');
end

function my_fieldmap_correction(nRun, EPI_data, SBRef_data, ...
             FM_m_data, FM_pd_data, FiledMap_TE, total_readout_time)
    n = 1;     
    %% FiledMap correction
    % Calculate VDM
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = FM_pd_data;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = FM_m_data;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = FiledMap_TE;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = total_readout_time;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 1;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = ...
        {'/Volumes/Data/matlab_tools/spm12/toolbox/FieldMap/T1.nii'};
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
    for xRun = 1 : nRun
        matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.session(xRun).epi = SBRef_data{xRun,1};
    end
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
    matlabbatch{n}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
    n = n + 1;
    
    % Apply VDM
    for xRun = 1 : nRun    
        matlabbatch{n}.spm.tools.fieldmap.applyvdm.data(xRun).scans = EPI_data{xRun,1};
        matlabbatch{n}.spm.tools.fieldmap.applyvdm.data(xRun).vdmfile(1) = ...
            cfg_dep(sprintf('Calculate VDM: Voxel displacement map (Subj 1, Session %d)', xRun), ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('()',{1}, '.','vdmfile', '{}',{xRun}));
    end
    matlabbatch{n}.spm.tools.fieldmap.applyvdm.roptions.pedir = 2;
    matlabbatch{n}.spm.tools.fieldmap.applyvdm.roptions.which = [2 1];
    matlabbatch{n}.spm.tools.fieldmap.applyvdm.roptions.rinterp = 4;
    matlabbatch{n}.spm.tools.fieldmap.applyvdm.roptions.wrap = [0 0 0];
    matlabbatch{n}.spm.tools.fieldmap.applyvdm.roptions.mask = 1;
    matlabbatch{n}.spm.tools.fieldmap.applyvdm.roptions.prefix = 'f';
    
    %% Run matlab batch
    spm('defaults', 'FMRI');
    spm_jobman('initcfg')
    spm_jobman('run', matlabbatch);
%  	save('matlabbatch.mat', 'matlabbatch');
end

function my_mutilband_preprocessing_2(nRun, EPI_data, SBRef_data, T1_data)
    clear matlabbatch
    n = 1;
    for xRun = 1 : nRun
        %% Coregistration
        matlabbatch{n}.spm.spatial.coreg.estimate.ref = SBRef_data{xRun, 1};
        matlabbatch{n}.spm.spatial.coreg.estimate.source = T1_data;
        matlabbatch{n}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{n}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{n}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{n}.spm.spatial.coreg.estimate.eoptions.tol = ...
            [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{n}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        n = n + 1;
        
        %% Segmentation
        matlabbatch{n}.spm.spatial.preproc.channel.vols(1) = ...
            cfg_dep('Coregister: Estimate: Coregistered Images', ...
            substruct('.','val', '{}',{(xRun - 1) * 3 + 1}, '.','val', '{}',{1},...
            '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
        matlabbatch{n}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{n}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{n}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{n}.spm.spatial.preproc.tissue(1).tpm = {'/Volumes/Data/matlab_tools/spm12/tpm/TPM.nii,1'};
        matlabbatch{n}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{n}.spm.spatial.preproc.tissue(1).native = [0 0];
        matlabbatch{n}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{n}.spm.spatial.preproc.tissue(2).tpm = {'/Volumes/Data/matlab_tools/spm12/tpm/TPM.nii,2'};
        matlabbatch{n}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{n}.spm.spatial.preproc.tissue(2).native = [0 0];
        matlabbatch{n}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{n}.spm.spatial.preproc.tissue(3).tpm = {'/Volumes/Data/matlab_tools/spm12/tpm/TPM.nii,3'};
        matlabbatch{n}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{n}.spm.spatial.preproc.tissue(3).native = [0 0];
        matlabbatch{n}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{n}.spm.spatial.preproc.tissue(4).tpm = {'/Volumes/Data/matlab_tools/spm12/tpm/TPM.nii,4'};
        matlabbatch{n}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{n}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{n}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{n}.spm.spatial.preproc.tissue(5).tpm = {'/Volumes/Data/matlab_tools/spm12/tpm/TPM.nii,5'};
        matlabbatch{n}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{n}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{n}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{n}.spm.spatial.preproc.tissue(6).tpm = {'/Volumes/Data/matlab_tools/spm12/tpm/TPM.nii,6'};
        matlabbatch{n}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{n}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{n}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{n}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{n}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{n}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{n}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{n}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{n}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{n}.spm.spatial.preproc.warp.write = [0 1];
        n = n + 1;
    
        %% Normalizeation: write
        matlabbatch{n}.spm.spatial.normalise.write.subj(1).def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        matlabbatch{n}.spm.spatial.normalise.write.subj(1).resample(1) = ...
            cfg_dep('Coregister: Estimate: Coregistered Images', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','cfiles'));
        
        matlabbatch{n}.spm.spatial.normalise.write.subj(2).def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{(xRun - 1) * 3 + 2}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        matlabbatch{n}.spm.spatial.normalise.write.subj(2).resample = SBRef_data{xRun, 1};
        
        matlabbatch{n}.spm.spatial.normalise.write.subj(3).def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{(xRun - 1) * 3 + 2}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        matlabbatch{n}.spm.spatial.normalise.write.subj(3).resample = EPI_data{xRun, 1};
        matlabbatch{n}.spm.spatial.normalise.write.woptions.bb = [-78, -112, -70; ...  
                                                                   78,   76,  85];
        matlabbatch{n}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{n}.spm.spatial.normalise.write.woptions.interp = 4;
        n = n + 1;
    end 
   
    %% Run matlab batch
    spm('defaults', 'FMRI');
    spm_jobman('initcfg')
    spm_jobman('run', matlabbatch);
%       save('matlabbatch.mat', 'matlabbatch');
end

function normalize_add_brainstem_3(nRun, EPI_data, SBRef_data, T1_data)
        clear matlabbatch
        for xRun = 1 : nRun
            % Normalize estimate
            matlabbatch{1}.spm.spatial.coreg.estimate.ref = SBRef_data{xRun, 1};
            matlabbatch{1}.spm.spatial.coreg.estimate.source = T1_data;
            matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            % Normalize write
            matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(1).source(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
            matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(1).wtsrc = '';
            matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(1).resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
            matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(2).source(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
            matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(2).wtsrc = '';
            matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(2).resample = SBRef_data{xRun, 1};
            matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(3).source(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
            matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(3).wtsrc = '';
            matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(3).resample = EPI_data{xRun, 1};
            matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.template = {
                                                                           '/Volumes/Data/matlab_tools/spm12/toolbox/OldNorm/T1.nii,1'
                                                                           '/Volumes/Data/matlab_tools/MNI152/T1_brainstem_masked.nii,1'
                                                                           };
            matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
            matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
            matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
            matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
            matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
            matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
            matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
            matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
            matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.bb = [-78 -112 -70
                                                                     78 76 85];
            matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.vox = [2 2 2];
            matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
            matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w2';
            
            % Run matlab batch
            spm('defaults', 'FMRI');
            spm_jobman('initcfg')
            spm_jobman('run', matlabbatch);
        end
end

function my_smoothing(fwhm, nRun, EPI_data, SBRef_data)
        clear matlabbatch
        spm('defaults', 'FMRI');
        spm_jobman('initcfg')
        for xRun = 1 : nRun
            matlabbatch{1}.spm.spatial.smooth.data = EPI_data{xRun, 1};
            matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = sprintf('s%d', fwhm);
            matlabbatch{2}.spm.spatial.smooth.data = SBRef_data{xRun, 1};
            matlabbatch{2}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
            matlabbatch{2}.spm.spatial.smooth.dtype = 0;
            matlabbatch{2}.spm.spatial.smooth.im = 0;
            matlabbatch{2}.spm.spatial.smooth.prefix = sprintf('s%d', fwhm);
            spm_jobman('run', matlabbatch);
        end
end