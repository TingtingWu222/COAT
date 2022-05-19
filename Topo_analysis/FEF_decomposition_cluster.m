function FEF_decomposition_cluster(path, path_output, ID)

    %% Parameters
    ROI_name = {'FEF', 'IPS'};
    n_perm = 1000;

    xROI = 1;
    % load ROI mask
    mask_name = fullfile(path, sprintf('%s_mask_cropped.nii', ROI_name{xROI}));
    [mask_box, xyz] = spm_read_vols(spm_vol(mask_name));
    mask_box(mask_box > 0) = 1;

    %% Loop over subjects
    for xSub = 1 : length(ID)
        fprintf('%s\n', ID{xSub});
        SubDir = fullfile(path_output, ID{xSub}, ROI_name{xROI}, 'cc');
        % Load ROI image
        filename = fullfile(path, sprintf('%s_%s_Weight_cropped_original.nii', ...
            ROI_name{xROI}, ID{xSub}));
        v = spm_vol(filename);
        img = spm_read_vols(v);

        %% Loop over connected components
        for n = 0 : 10%n_perm
            % Connected components: positive
            filename = fullfile(SubDir, sprintf('%s_%s_positive_cc_r_%d.mat', ...
                ROI_name{xROI}, ID{xSub}, n));
            load(filename); clear filename;
            cc = double(cc);
            idx_background = cc(1);
            cc(mask_box == 0) = 0; cc(img <= 0) = 0;
            t = cc; t(t ~= idx_background) = 0;
            cc_t = bwlabeln(t,26);
            cc(cc == idx_background) = 0;
            a = unique(cc);
            for i = 2 : length(a)
                cc(cc == a(i)) = i - 1;
            end
            cc_t(cc_t~=0) = cc_t(cc_t~=0) + length(a) - 1; 
            cc_t(mask_box == 0) = 0;
            cc_p = cc + cc_t; 
            clear a t cc cc_t idx_background
            b = unique(cc_p);
            for i = 2 : length(b)
                k_positive(i-1,1) = sum(cc_p(:) == b(i));
            end
            n_cluster_positive(xSub,n+1) = length(b) - 1;
            cluster_size_positive(xSub, n+1) = mean(k_positive);
            clear k_positive cc cc_t t a b;

            % Connected components: negative
            filename = fullfile(SubDir, sprintf('%s_%s_negative_cc_r_%d.mat', ...
                ROI_name{xROI}, ID{xSub}, n));
            load(filename); clear filename;
            cc = double(cc);
            idx_background = cc(1);
            cc(mask_box == 0) = 0; cc(img >= 0) = 0;
            t = cc; t(t ~= idx_background) = 0;
            cc_t = bwlabeln(t,26);
            cc(cc == idx_background) = 0;
            a = unique(cc);
            for i = 2 : length(a)
                cc(cc == a(i)) = i - 1;
            end
            cc_t(cc_t~=0) = cc_t(cc_t~=0) + length(a) - 1; 
            cc_t(mask_box == 0) = 0;
            cc_n = cc + cc_t; 
            clear a t cc cc_t idx_background

            b = unique(cc_n);
            for i = 2 : length(b)
                k_negative(i-1,1) = sum(cc_n(:) == b(i));
            end
            n_cluster_negative(xSub,n+1) = length(b) - 1;
            cluster_size_negative(xSub, n+1) = mean(k_negative);
            clear k_negative cc cc_t t a b;
        end    
    end

    %% Group level
    n_cluster(:,1) = n_cluster_positive(:,1); 
    n_cluster(:,3) = n_cluster_negative(:,1); 
    n_cluster(:,2) = mean(n_cluster_positive(:,2:end),2);
    n_cluster(:,4) = mean(n_cluster_negative(:,2:end),2);

    cluster_size(:,1) = cluster_size_positive(:,1); 
    cluster_size(:,3) = cluster_size_negative(:,1);
    cluster_size(:,2) = mean(cluster_size_positive(:,2:end),2);
    cluster_size(:,4) = mean(cluster_size_negative(:,2:end),2);

    vol = cluster_size * 8;
    isotopic_size = nthroot(vol,3);
    save(fullfile(path_output, 'Topo_results.mat'), 'n_cluster', 'cluster_size', 'vol', 'isotopic_size');
end