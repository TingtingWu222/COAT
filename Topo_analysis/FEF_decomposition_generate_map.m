function FEF_decomposition_generate_map(path, path_output, ID, xROI)
    %% Parameters
    ROI_name = {'FEF', 'IPS'};
    n_perm = 1000;

    % load ROI mask
    mask_name = fullfile(path, sprintf('%s_mask_cropped.nii', ROI_name{xROI}));
    [mask_box, xyz] = spm_read_vols(spm_vol(mask_name));
    mask_box(mask_box > 0) = 1;

    %% Loop over subjects
    for xSub = 1 : length(ID)
        fprintf('%s\n', ID{xSub});
        SubDir = fullfile(path_output, ID{xSub});
        if ~exist(fullfile(SubDir, ROI_name{xROI}),'dir')
             mkdir(fullfile(SubDir, ROI_name{xROI}, 'cc'));
        end
        % Load ROI image
        filename = fullfile(path, sprintf('%s_%s_Weight_cropped_original.nii', ...
            ROI_name{xROI}, ID{xSub}));
        v = spm_vol(filename);
        V = spm_read_vols(v);

        % Seperate positive and negative voxels
        above_zero = V; above_zero = above_zero * (-1); 
        above_zero(V < 0) = 10; above_zero(mask_box == 0) = 10;
        below_zero = V; below_zero(V > 0) = 10;  below_zero(mask_box == 0) = 10;

        % save origianl data
        filename_p = fullfile(SubDir,ROI_name{xROI},   sprintf('%s_%s_Weight_cropped_positive_r_0.nii', ...
                     ROI_name{xROI}, ID{xSub}));
        vNew = v; vNew.fname = filename_p; spm_write_vol(vNew, above_zero); clear vNew;    
        filename_n = fullfile(SubDir, ROI_name{xROI},  sprintf('%s_%s_Weight_cropped_negative_r_0.nii', ...
                     ROI_name{xROI}, ID{xSub}));
        vNew = v; vNew.fname = filename_n; spm_write_vol(vNew, below_zero); clear vNew;

        % randomize voxels
        for i = 1 : n_perm
            I = find(mask_box == 1);
            V_r = V;  V_r(I) = V_r(I(randperm(length(I))));
            % Seperate positive and negative voxels
            above_zero_r = V_r; above_zero_r = above_zero_r * (-1); 
            above_zero_r(V_r < 0) = 10; above_zero_r(mask_box == 0) = 10;
            below_zero_r = V_r; below_zero_r(V_r > 0) = 10;  below_zero_r(mask_box == 0) = 10;

            % save origianl data
            SubDir = fullfile(path_output, ID{xSub});
            if ~exist(SubDir, 'dir'); mkdir(SubDir); end
            filename_p = fullfile(SubDir, ROI_name{xROI}, sprintf('%s_%s_Weight_cropped_positive_r_%d.nii', ...
                         ROI_name{xROI}, ID{xSub},i));
            vNew = v;    vNew.fname = filename_p;
            spm_write_vol(vNew, above_zero_r); clear vNew;

            filename_n = fullfile(SubDir, ROI_name{xROI},  sprintf('%s_%s_Weight_cropped_negative_r_%d.nii', ...
                         ROI_name{xROI}, ID{xSub}, i));
            vNew = v;    vNew.fname = filename_n;
            spm_write_vol(vNew, above_zero_r); clear vNew;
        end

    end
end