function FEF_ST_cross_classification(path, ID, xROI, cross)
    % This script is for the cross-cueing classification
    % Inputs:
    % path: directtory of the ROI single-trial activation .mat files
    % ID: list of subject's ID
    % xROI: 1: FEF, 2: IPS, 3: FEF + FPN
    % cross: 1: Exogenous to Endogenous, 2: Endogenous to Exogenous
    %% Parameters
    mode = 'unsmoothed_EPI';
    ROI_list = {'FEF','IPS','FPN'};
    classname = {'covert','overt'}; 

    %% Loop over subjects
    for cross = 1 : 2
         %% output file name
         if cross == 1
            filename = fullfile(path, ['SVM_cross_Ex2En_Results_', ROI_list{xROI},'.mat']);
         elseif cross == 2
            filename = fullfile(path, ['SVM_cross_En2Ex_Results_', ROI_list{xROI},'.mat']);
         end
        for xSub = 1 : length(ID)
            fprintf('Sub%s\t',ID{xSub});
             % labels & exclusion flags
             load(fullfile(path, ID{xSub},[ROI_list{xROI},'_dat'], 'TrialType.mat'));
             SubDir = fullfile(path, ID{xSub},[ROI_list{xROI},'_dat'], mode);


             %% load data   
             if xROI == 1
               Trial_type = Trial_type_valid(Trial_FEF_outlier(:,1)+Trial_GM_outlier == 0);
               load(fullfile(SubDir, 'FEF_vector.mat'));
               t = (Trial_GM_outlier(Trial_FEF_outlier(:,1)==0, 1) == 0);
               Feature = zscore(FEF_v(t,:));    clear FEF_v   
            elseif xROI == 2
               Trial_type = Trial_type_valid(Trial_IPS_outlier(:,1)+Trial_GM_outlier == 0);
               load(fullfile(SubDir, 'IPS_vector.mat'));
               t = (Trial_GM_outlier(Trial_IPS_outlier(:,1)==0, 1)==0);
               Feature = zscore(IPS_v(t,:));    clear IPS_v   
            elseif xROI == 3
               Trial_type = Trial_type_valid(Trial_FPN_outlier(:,1)+Trial_GM_outlier == 0);
               load(fullfile(SubDir, 'FPN_vector.mat'));
               t = (Trial_GM_outlier(Trial_FPN_outlier(:,1)==0, 1)==0);
               Feature = zscore(FPN_v(t,:));    clear FPN_v   
             end
             y = cell(length(Trial_type),1);

             if cross == 1
                I_train = find(Trial_type == 1 | Trial_type == 3);
                I_test = find(Trial_type == 2 | Trial_type == 4);
                y(Trial_type <= 2) = {'covert'};  y(Trial_type >= 3) = {'overt'};
             elseif cross == 2
                I_train = find(Trial_type == 2 | Trial_type == 4);
                I_test = find(Trial_type == 1 | Trial_type == 3);
                y(Trial_type <= 2) = {'covert'};  y(Trial_type >= 3) = {'overt'};
             elseif cross == 3
                I_train = find(Trial_type == 1 | Trial_type == 2);
                I_test = find(Trial_type == 3 | Trial_type == 4);
                y(Trial_type == 1 | Trial_type == 3) = {'ex'};  y(Trial_type == 2 | Trial_type == 4) = {'en'};
             end

             %% Classification
             SVMModel = fitcsvm(Feature(I_train,:),y(I_train),...
               'KernelFunction','linear','Standardize',false,'Solver', 'L1QP', ...
               'ClassNames',classname);
             label = predict(SVMModel,Feature(I_test,:));
             ACC(xSub,cross) = sum(strcmp(y(I_test),label))/length(I_test);
             Precision_class1(xSub,cross) = sum(strcmp(y(I_test),classname{1}).* ...
                                      strcmp(label,classname{1}))/sum(strcmp(label,classname{1}));
             Recall_class1(xSub,cross) = sum(strcmp(y(I_test),classname{1}).* ...
                                 strcmp(label,classname{1}))/sum(strcmp(y(I_test),classname{1}));
             F1_score_class1(xSub,cross) = 2*(Recall_class1(xSub,cross) * Precision_class1(xSub,cross)) / ...
                                     (Recall_class1(xSub,cross) + Precision_class1(xSub,cross));

             Precision_class2(xSub,cross) = sum(strcmp(y(I_test),classname{2}).* ...
                                    strcmp(label,classname{2}))/sum(strcmp(label,classname{2}));
             Recall_class2(xSub,cross) = sum(strcmp(y(I_test),classname{2}).* ...
                                 strcmp(label,classname{2}))/sum(strcmp(y(I_test),classname{2}));
             F1_score_class2(xSub,cross) = 2*(Recall_class2(xSub,cross) * Precision_class2(xSub,cross)) / ...
                                    (Recall_class2(xSub,cross) + Precision_class2(xSub,cross));
             fprintf('ACC=%f\n',ACC(xSub,cross));
             clear SVMModel I_train I_test Trial_type
        end
    end
    %% Save results
    save(filename, 'ACC', 'Precision_class1','Recall_class1','F1_score_class1',...
         'Precision_class2','Recall_class2','F1_score_class2');
end
