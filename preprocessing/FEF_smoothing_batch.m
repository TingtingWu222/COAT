function FEF_smoothing_batch(path, ID)
    %% Parameters
    nRun = 4;
    fwhm = 8;

    %% Initialise SPM

    %% Loop over subjects
    gcp;
    parfor xSub =  1 :nSub
        fprintf('Subject%s\n', ID{xSub});
        subDir = fullfile(path, ID{xSub});
        FEF_smoothing(subDir, fwhm, nRun);
    end
end

function FEF_smoothing(subDir, fwhm, nRun)
        clear matlabbatch
        spm('defaults', 'FMRI');
        spm_jobman('initcfg')
        for xRun = 1 : nRun
            RunDir = fullfile(subDir, sprintf('RUN%d',xRun));        
            matlabbatch{1}.spm.spatial.smooth.data = cellstr(spm_select('FPList', RunDir, '^wu.*.\.img'));
            matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = sprintf('s%d', fwhm);
            spm_jobman('run', matlabbatch);
        end
end