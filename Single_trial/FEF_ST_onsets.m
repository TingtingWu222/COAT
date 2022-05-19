function FEF_ST_onsets(path, ID)

    %% Loop over subject
    for xSub = 1 : nSub
        fprintf('Subject%s\n', ID{xSub});
        subDir = fullfile(path,ID{xSub});
        BehvDir = fullfile(subDir,'Behv');
        EyeTrackerDir = fullfile(subDir, 'EyeTracker');
        cd(subDir);
        filename = spm_select('List', BehvDir, '^*.\.xlsx');
        filename_eye = sprintf('Sub%s_EyeTracker_results.mat', ID{xSub});
        FEF_extract_ST_onsets(subDir, BehvDir, filename, EyeTrackerDir, filename_eye);
    end
end

function FEF_extract_ST_onsets(subDir, BehvDir, filename, EyeTrackerDir, filename_eye)
    %% load data
    [~,~,rawData] = xlsread(fullfile(BehvDir,filename));
    conditionNames = {'CoEn','CoEx','OvEn','OvEx'};
    shiftVal = 2;
    load(fullfile(EyeTrackerDir, filename_eye));
%     flag = 0;
    %% Orgnize the excel file
    Label = rawData(shiftVal,:);
    Block = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'Block')));
    Fix30sOnsets = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'Fixation15secondsStart.OnsetTime')));
    Trial_type = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'type')));
    CueOnsetTime = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'SlideCue.OnsetTime')));
    ACC = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'SlideTarget.ACC')));
    RT = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'SlideTarget.RT')));
    xy = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'xy1') : (findLabel(Label,'xy1') + 7)));
    for xT = 1 : length(xy)
        Loc(xT, 1) = find(xy(xT, :) == 1);
    end
    CueOnset = (CueOnsetTime - Fix30sOnsets) / 1000;
    nRun = length(unique(Block));

    %% Extract onsets for each trial
    block_idx = unique(Fix30sOnsets);
    t = 1;
    for xTrial = 1 : length(CueOnset)        
        temp = find(Trial_type == Trial_type(xTrial) & Block == Block(xTrial) & ...
                    find(CueOnset) ~= xTrial & ACC == 1 & eye_movement_outlier == 0, 1);
        if ACC(xTrial) == 1 && eye_movement_outlier(xTrial) == 0 && ~isempty(temp)
           if xTrial < 10
              TrialDir = sprintf('%s/Single_trial/Trial_00%d/Onsets',subDir,xTrial);
           elseif xTrial < 100
              TrialDir = sprintf('%s/Single_trial/Trial_0%d/Onsets',subDir,xTrial);
           else
              TrialDir = sprintf('%s/Single_trial/Trial_%d/Onsets',subDir,xTrial);
           end
           if ~exist(TrialDir,'dir'); mkdir(TrialDir); end
           Trial_type_valid(t,1) = Trial_type(xTrial);
           Trial_Loc_valid(t,1) = Loc(xTrial);
           Trial_Num_valid(t,1)  = xTrial;
           Contrast_vector{t, 1} = [];
            %% Make regressors
            for xRun = 1 : nRun
                xReg = 1;

                % Regressor for current trial
                temp1 = CueOnset(Fix30sOnsets == block_idx(xRun) & find(CueOnset) == xTrial);
                if ~isempty(temp1)
                    nRun_CurrentTrial = xRun;
                    onsets{xReg} = temp1;
                    names{xReg} = sprintf('Run%dCond%dCurrent_Trial',xRun,xReg);
                    durations{xReg } = 0;                    
                    xReg = xReg + 1;
                end

                % Regressor for trials with correct responses
                for xType = 1 : length(conditionNames)
                    temp = CueOnset(Trial_type == xType & Block == xRun & find(CueOnset) ~= xTrial & ...
                                            ACC == 1 & eye_movement_outlier == 0);
                    if ~isempty(temp)
                        onsets{xReg} = temp;
                        names{xReg}  = sprintf('Run%dCond%d%sCorr', xRun, xReg, conditionNames{xType} );
                        durations{xReg} = 0;
                        xReg = xReg + 1;
                    end
                end

                % Regressor for trials with incorrect responses
                for xType = 1 : length(conditionNames)
                    temp = CueOnset(Trial_type == xType & Block == xRun & ACC == 0);
                    if ~isempty(temp)
                        onsets{xReg} = temp;
                        names{xReg}  = sprintf('Run%dCond%d%sError', xRun, xReg, conditionNames{xType} );
                        durations{xReg} = 0;
                        xReg = xReg + 1;
                    end
                    clear temp
                end

                 % Regressor for trials with eye movement artifact
                for xType = 1 : length(conditionNames)
                    temp = CueOnset(Trial_type == xType & Block == xRun & ACC == 1 & eye_movement_outlier == 1);
                    if ~isempty(temp)
                        onsets{xReg} = temp;
                        names{xReg}  = sprintf('Run%dCond%d%sOutlier', xRun, xReg, conditionNames{xType} );
                        durations{xReg} = 0;
                        xReg = xReg + 1;
                    end
                    clear temp
                end

                % Save onset vectors
                Ofilename = sprintf('%s/Run%dOnsets.mat',TrialDir,xRun);
                save(Ofilename,'onsets','names','durations');
                
                weight = zeros(1, xReg - 1);
                if ~isempty(temp1)
                    weight(1) = 1;
                end
                Contrast_vector{t, 1} = [Contrast_vector{t, 1}, weight, zeros(1,6)];

                clear onsets names durations xReg RawData Label
            end
            t = t + 1;
        end
    end
    assignin('base', 'Contrast_vector',Contrast_vector);
    filename_nTrial = sprintf('%s/Single_trial/TrialType.mat',subDir);
    save(filename_nTrial,'Trial_type_valid', 'Trial_Loc_valid', 'Trial_Num_valid', 'Contrast_vector');
    clear Trial_type_valid Trial_Loc _valid Trial_Num_valid t Contrast_vector  
end

function [ flag ] = findLabel( Label,target)
%
% Function to return the column number in the label list which matchs the target
%
% Inputs:
% - Label: a cells structure contains all labels, exctracted from the
%          'shiftVal' row of the Excel file
% - target: a string, the name of the target label
%
% Output:
% - flag: the column number of the target label in original Excel file
%
    celllength = cellfun('length',Label); %Number of strings in each cell
    
    %find cells which contain target string
    targetcell = strfind(Label,target);  
    
    %transform targetcell from cell to matrix, empty cells are replaced
    %with 0
    temp = cellfun('isempty',targetcell);%find all empty cells
    for i = 1:length(temp)
        if temp(i) == 1
            targetcell{i} = 0;
        end
    end
    targetmat=cell2mat(targetcell);
    
    flag = 0;
    for i = 1:length(Label)
        if targetmat(i) == 1 && celllength(i) == length(target)%5
            flag = i;
        end
    end

end
