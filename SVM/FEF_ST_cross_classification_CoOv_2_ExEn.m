function FEF_ST_cross_classification_CoOv_2_ExEn(path, ID, xROI)
% This script is for the cross "Covert vs. Overt" to "Enxogenous to Endogenous" classification
% Inputs:
% path: directtory of the ROI single-trial activation .mat files
% ID: list of subject's ID
% xROI: 1: FEF, 2: IPS, 3: FEF + FPN
% xM: 1: unsmoothed data,2: smoothed data

    %% Parameters
    mode = 'unsmoothed_EPI';
    classname = {'covert','overt'}; 
    ROI_list = {'FEF','IPS','FPN'};
    filename = fullfile(path, ['SVM_cross_CoOv_2_ExEn_Results_', ROI_list{xROI},'.mat']);
    n = 1000; % number of permutations
    
    %% Loop over subjects
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
         y = cell(length(Trial_type),1); y2 = cell(length(Trial_type),1);
         y(Trial_type <= 2) = {'covert'};  y(Trial_type >= 3) = {'overt'};
         y2(Trial_type == 1 | Trial_type == 3) = {'covert'};  % Ex
         y2(Trial_type == 2 | Trial_type == 4) = {'overt'};   % En

         %% Classification
         for i = 1 : n
             t = randperm(length(y));
             I_train = t(1:ceil(length(y)/2));
             I_test = t((ceil(length(y)/2)+1):end);
             SVMModel = fitcsvm(Feature(I_train,:), y(I_train,:), 'KernelFunction','linear','Standardize',false,'Solver', 'L1QP', 'ClassNames',classname);
             label = predict(SVMModel,Feature(I_test,:));
             ACC(xSub,i) = sum(strcmp(y2(I_test),label))/length(I_test);
             Precision_class1(xSub,i) = sum(strcmp(y2(I_test),classname{1}).* ...
                                      strcmp(label,classname{1}))/sum(strcmp(label,classname{1}));
             Recall_class1(xSub,i) = sum(strcmp(y2(I_test),classname{1}).* ...
                                 strcmp(label,classname{1}))/sum(strcmp(y(I_test),classname{1}));
             F1_score_class1(xSub,i) = 2*(Recall_class1(xSub,i) * Precision_class1(xSub,i)) / ...
                                     (Recall_class1(xSub,i) + Precision_class1(xSub,i));

             Precision_class2(xSub,i) = sum(strcmp(y2(I_test),classname{2}).* ...
                                    strcmp(label,classname{2}))/sum(strcmp(label,classname{2}));
             Recall_class2(xSub,i) = sum(strcmp(y2(I_test),classname{2}).* ...
                                 strcmp(label,classname{2}))/sum(strcmp(y(I_test),classname{2}));
             F1_score_class2(xSub,i) = 2*(Recall_class2(xSub,i) * Precision_class2(xSub,i)) / ...
                                    (Recall_class2(xSub,i) + Precision_class2(xSub,i));
          clear SVMModel  Trial_type I_train I_test t
         end
          fprintf('ACC=%f\n',mean(ACC(xSub,:)));
    end

    % %% Save results
    save(filename, 'ACC', 'Precision_class1','Recall_class1','F1_score_class1',...
         'Precision_class2','Recall_class2','F1_score_class2');
end


