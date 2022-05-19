function FEF_ST_classification(path, ID, xROI, xM)
% This script is for the overt vs. covert SVM classification
% Inputs:
% path: directtory of the ROI single-trial activation .mat files
% ID: list of subject's ID
% xROI: 1: FEF, 2: IPS, 3: FEF + FPN
% xM: 1: unsmoothed data,2: smoothed data
    clc; 

    %% Parameters
    mode = {'unsmoothed_EPI','smoothed_EPI'};
    ROI_list = {'FEF','IPS','FPN'};

    k = 10; %k-fold
    n_perm = 1000; % number of permutations

    %% Loop over subjects
    for xSub = 1 : length(ID)
        fprintf('Sub%s\t',ID{xSub});
        SubDir = fullfile(path, ID{xSub}, [ROI_list{xROI},'_dat'], mode{xM});

        %% load data
        fprintf('Load data\t');
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
         Trial_type = Trial_type_valid(Trial_FEF_outlier(:,xM) == 0 & Trial_GM_outlier == 0);
         y = cell(length(Trial_type),1); 
         y(Trial_type <= 2) = {'covert'};  y(Trial_type >= 3) = {'overt'};

         %% Classification
         SVM_Results = my_SVM_FEF(Feature, Target_coov, k, n_perm, 0);
         fprintf('Accuracy = %f\t',mean(SVM_Results.ACC(:)));

         %% Chance level baseline
         fprintf('Estimating chance level ...\n');
         SVM_Results_baseline = my_SVM_FEF(Feature, Target_coov, k, n_perm, 1);

         %% Save results
         filename = fullfile(SubDir, 'SVM_Results_voxelwise.mat');
         save(filename, 'SVM_Results','SVM_Results_baseline');
    end
end

function SVM_Results = my_SVM_FEF(Feature, Target_coov, k, nP, shuffle_target)
        classname = {'covert','overt'}; 
        y = Target_coov;
        
         
         %% Loop over permutations
         SVM_Results.weight = [];
         for n = 1 : nP
             if shuffle_target == 1
                 t = randperm(length(Target_coov));
                 y = Target_coov(t);                     
             end
             %% k-fold
             c = cvpartition(y,'KFold',k);
             for i = 1 : k                    
                 % Training SVM
                 SVMModel = fitcsvm(Feature(c1.test(i)==0,:),y_train(c.test(i)==0),...
                   'KernelFunction','linear','Standardize',false,'Solver', 'L1QP', ...
                   'ClassNames',classname);

                  % Training accuracy
                  label = predict(SVMModel,Feature(c.test(i)==0,:));
                  SVM_Results.ACC_train(n,i) = sum(strcmp(y(c.test(i)==0),label))/sum(c.test(i)==0);
                  clear label

                  % Testing performance
                  label = predict(SVMModel,Feature(c.test(i)==1,:));
                  SVM_Results.ACC(n,i) = sum(strcmp(y(c.test(i)==1),label))/sum(c.test(i)==1);

                  SVM_Results.Precision_class1(n,i) = sum(strcmp(y(c.test(i)==1),classname{1}).* ...
                                          strcmp(label,classname{1}))/sum(strcmp(label,classname{1}));
                  SVM_Results.Recall_class1(n,i) = sum(strcmp(y(c.test(i)==1),classname{1}).* ...
                                     strcmp(label,classname{1}))/sum(strcmp(y(c.test(i)==1),classname{1}));
                  SVM_Results.F1_score_class1(n,i) = 2*(Recall_class1(1,i) * Precision_class1(1,i)) / ...
                                         (Recall_class1(1,i) + Precision_class1(1,i));

                  SVM_Results.Precision_class2(n,i) = sum(strcmp(y(c.test(i)==1),classname{2}).* ...
                                        strcmp(label,classname{2}))/sum(strcmp(label,classname{2}));
                  SVM_Results.Recall_class2(n,i) = sum(strcmp(y(c.test(i)==1),classname{2}).* ...
                                     strcmp(label,classname{2}))/sum(strcmp(y(c.test(i)==1),classname{2}));
                  SVM_Results.F1_score_class2(n,i) = 2*(Recall_class2(1,i) * Precision_class2(1,i)) / ...
                                        (Recall_class2(1,i) + Precision_class2(1,i));
                  
                  w(i,:) = SVMModel.Beta'; 
                  clear SVMModel
             end
             SVM_Results.weight = [SVM_Results.weight; w];
         end
end