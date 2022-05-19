%% clear
clear all; clc; close all;

%% Directories
path = fullfile('/Volumes/Data/DataMelissa/COAT_RAW/FEF/Analyses_new');
ID = '21010';

%% Parameters
nRun = 3;
conditionNames = {'CoEn','CoEx','OvEn','OvEx'};
shiftVal = 2;
nSD = 3;
w = 64;

%% Plot parameters
res = [1600, 1200];
circle_loc_x = [0,  13, 18,  13,  0, -13, -18, -13]/100 * res(1);
circle_loc_y = [33, 23,  0, -23, -33, -23,  0,  23]/100 * res(2);
r = res(2) * 0.13/2;       
covert_ROI_r = res(2)/2/4;
overt_ROI_r = res(2)/2;
ang=0:0.2:2*pi; 
marker_color = {'b', 'b', 'r', 'r'};
line_type = {'bo-', 'bo--', 'ro-', 'ro--'};

%% Inputs
subDir = fullfile(path,ID);
BehvDir = fullfile(subDir,'Behv');
EyeTrackerDir = fullfile(subDir, 'EyeTracker');
%         cd(subDir);
filename = spm_select('List', BehvDir, '^*.\.xlsx');
filename_eye = sprintf('%s_GA.mat', ID);

%% read data
[~,~,rawData] = xlsread(fullfile(BehvDir,filename));
load(fullfile(EyeTrackerDir, filename_eye));
eye_movement_outlier = zeros(length(corr_par), 1);

%% Orgnize the excel file
Label = rawData(shiftVal,:);
Block = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'Block')));
Trial_type = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'type')));
ACC = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'SlideTarget.ACC')));
TrialCodeCue = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'TrialCodeCue')));
Block = Block(129:end);   ACC = ACC(129:end);
Trial_type = Trial_type(129:end);
TrialCodeCue = TrialCodeCue(129:end);

Cue_loc = mod(TrialCodeCue, 10);        
loc_n = unique(Cue_loc);
corr_par_x_raw = corr_par(:,3);
corr_par_y_raw = corr_par(:,4);

clear rawData

%% Calibration within each run, using moving window
for xRun = 1 : nRun
    corr_par_run = corr_par(Block == xRun, :);
    l = length(corr_par_run);
    for xTrial = 1 : l
        if xTrial > w/2 && xTrial < l - w/2
           tx = corr_par_run((xTrial - w/2) : (xTrial + w/2), 3); tx = tx(~isnan(tx));
           ty = corr_par_run((xTrial - w/2) : (xTrial + w/2), 4); ty = ty(~isnan(ty));               
        elseif xTrial <= w/2
           tx = corr_par_run(1 : (w/2 + xTrial -1), 3); tx = tx(~isnan(tx));
           ty = corr_par_run(1 : (w/2 + xTrial -1), 4); ty = ty(~isnan(ty));
        else
           tx = corr_par_run((xTrial - w/2) : l, 3); tx = tx(~isnan(tx));
           ty = corr_par_run((xTrial - w/2) : l, 4); ty = ty(~isnan(ty));               
        end            
           corr_par_x(xTrial, xRun) = corr_par_run(xTrial, 3) - mean(tx); clear tx
           corr_par_y(xTrial, xRun) = corr_par_run(xTrial, 4) - mean(ty); clear ty
    end
    clear corr_par_run
end
corr_par_x = reshape(corr_par_x, length(corr_par), 1);
corr_par_y = reshape(corr_par_y, length(corr_par), 1);

%% Remove trials with eye location outside 3SD of the mean of overt condition and recalibarate the data
overt_x_raw_mean = mean(corr_par_x);%(Trial_type == 3 | Trial_type == 4));
overt_x_raw_sd = std(corr_par_x);%(Trial_type == 3 | Trial_type == 4));
overt_y_raw_mean = mean(corr_par_x);%(Trial_type == 3 | Trial_type == 4));
overt_y_raw_sd = std(corr_par_y);%(Trial_type == 3 | Trial_type == 4));
for xTrial = 1 : length(corr_par_x)
    if corr_par_x(xTrial) > overt_x_raw_mean + nSD * overt_x_raw_sd || ...
       corr_par_x(xTrial) < overt_x_raw_mean - nSD * overt_x_raw_sd || ...
       corr_par_y(xTrial) > overt_y_raw_mean + nSD * overt_y_raw_sd || ...
       corr_par_y(xTrial) < overt_y_raw_mean - nSD * overt_y_raw_sd
       eye_movement_outlier(xTrial, 1) = 1;
    end                                
end
clear overt_x_raw overt_y_raw

corr_par_x_c = corr_par_x;  corr_par_x_c(eye_movement_outlier == 1) = nan;
corr_par_y_c = corr_par_y;  corr_par_y_c(eye_movement_outlier == 1) = nan;
for xRun = 1 : nRun
    corr_par_x_run = corr_par_x(Block == xRun, 1);
    corr_par_y_run = corr_par_y(Block == xRun, 1);
    corr_par_x_c_run = corr_par_x_c(Block == xRun, 1);
    corr_par_y_c_run = corr_par_y_c(Block == xRun, 1);
    l = length(corr_par_x_run);
    for xTrial = 1 : l
        if xTrial > w/2 && xTrial < l - w/2
           tx = corr_par_x_c_run((xTrial - w/2) : (xTrial + w/2)); tx = tx(~isnan(tx));
           ty = corr_par_y_c_run((xTrial - w/2) : (xTrial + w/2)); ty = ty(~isnan(ty));               
        elseif xTrial <= w/2
           tx = corr_par_x_c_run(1 : (w/2 + xTrial -1)); tx = tx(~isnan(tx));
           ty = corr_par_y_c_run(1 : (w/2 + xTrial -1)); ty = ty(~isnan(ty)); 
        else
           tx = corr_par_x_c_run((xTrial - w/2) : l); tx = tx(~isnan(tx));
           ty = corr_par_y_c_run((xTrial - w/2) : l); ty = ty(~isnan(ty));           
        end
       corr_par_x2(xTrial, xRun) = corr_par_x_run(xTrial) - mean(tx); clear tx
       corr_par_y2(xTrial, xRun) = corr_par_y_run(xTrial) - mean(ty); clear ty
    end
    clear corr_par_x_run corr_par_y_run corr_par_x_c_run corr_par_y_c_run
end
corr_par_x = reshape(corr_par_x2, length(corr_par), 1);
corr_par_y = reshape(corr_par_y2, length(corr_par), 1);

%% Remove covert trials with eye location outside the  ROI and
%  overt tirals with eye location within the covert ROI
for xTrial = 1 : length(corr_par_x)
    if isnan(corr_par_x(xTrial)) || isnan(corr_par_y(xTrial))
        eye_movement_outlier(xTrial, 1) = 1;
    end
    if Trial_type(xTrial) == 1 || Trial_type(xTrial) == 2
        if sqrt(corr_par_x(xTrial)^2 + corr_par_y(xTrial)^2) > covert_ROI_r
           eye_movement_outlier(xTrial, 1) = 1;
        end
    else
        if sqrt(corr_par_x(xTrial)^2 + corr_par_y(xTrial)^2) < covert_ROI_r || ...
           sqrt(corr_par_x(xTrial)^2 + corr_par_y(xTrial)^2) > overt_ROI_r 
           eye_movement_outlier(xTrial, 1) = 1;
        end
    end
end

%% Remove trials with eye movement excess 3SD of mean in each condition + location
for xType = 1 : length(conditionNames)
    for xLoc = 1 : length(loc_n)
        eye_loc_x{xType, xLoc} = corr_par_x(Trial_type == xType  & Cue_loc == xLoc & ACC == 1 & eye_movement_outlier == 0);
        eye_loc_y{xType, xLoc} = corr_par_y(Trial_type == xType  & Cue_loc == xLoc & ACC == 1 & eye_movement_outlier == 0);

        %% Remove trials with eye movement excess 3SD of mean                
        eye_loc_x_t_mean(xType, xLoc) = mean(eye_loc_x{xType, xLoc});
        eye_loc_x_t_sd(xType, xLoc)   = std(eye_loc_x{xType, xLoc});                
        eye_loc_y_t_mean(xType, xLoc) = mean(eye_loc_y{xType, xLoc});
        eye_loc_y_t_sd(xType, xLoc)   = std(eye_loc_y{xType, xLoc});

        eye_loc_x_upper(xType, xLoc) = eye_loc_x_t_mean(xType, xLoc) + nSD * eye_loc_x_t_sd(xType, xLoc);
        eye_loc_x_lower(xType, xLoc) = eye_loc_x_t_mean(xType, xLoc) - nSD * eye_loc_x_t_sd(xType, xLoc);
        eye_loc_y_upper(xType, xLoc) = eye_loc_y_t_mean(xType, xLoc) + nSD * eye_loc_y_t_sd(xType, xLoc);
        eye_loc_y_lower(xType, xLoc) = eye_loc_y_t_mean(xType, xLoc) - nSD * eye_loc_y_t_sd(xType, xLoc);
    end
end

for xTrial = 1 : length(corr_par_x)
    if corr_par_x(xTrial) > eye_loc_x_upper(Trial_type(xTrial), Cue_loc(xTrial)) || ...
       corr_par_x(xTrial) < eye_loc_x_lower(Trial_type(xTrial), Cue_loc(xTrial)) || ...
       corr_par_y(xTrial) > eye_loc_y_upper(Trial_type(xTrial), Cue_loc(xTrial)) || ...
       corr_par_y(xTrial) < eye_loc_y_lower(Trial_type(xTrial), Cue_loc(xTrial))
       eye_movement_outlier(xTrial, 1) = 1;
    end                                
end
clear eye_loc_x_upper eye_loc_x_lower eye_loc_y_upper eye_loc_y_lower eye_loc_x_t* eye_loc_y_t*


for xType = 1 : length(conditionNames)
    for xLoc = 1 : length(loc_n)
        eye_loc_x_clean{xType, xLoc} = corr_par_x(Trial_type == xType  & Cue_loc == xLoc ...
                                & ACC == 1 & eye_movement_outlier == 0);
        eye_loc_y_clean{xType, xLoc} = corr_par_y(Trial_type == xType  & Cue_loc == xLoc ...
                                & ACC == 1 & eye_movement_outlier == 0);
        eye_loc_x_mean(xType, xLoc) = mean(eye_loc_x_clean{xType, xLoc});
        eye_loc_y_mean(xType, xLoc) = mean(eye_loc_y_clean{xType, xLoc});
    end
end

%% Save
save(fullfile(EyeTrackerDir, sprintf('Sub%s_EyeTracker_results.mat', ID)),...
    'eye_loc_x_clean',  'eye_loc_x_mean', 'eye_loc_y_clean',  'eye_loc_y_mean',...
    'eye_movement_outlier');

 %% plot scatter of raw data
 close all; 
figure; hold on
for i = 1 : length(corr_par_x);
    if Trial_type(i) == 1 || Trial_type(i) == 2
        scatter(corr_par_x_raw(i), corr_par_y_raw(i), 'o', 'MarkerEdgeColor',[1 1 1],...
            'MarkerFaceColor',[0, 0, i/length(corr_par) * 0.8]); hold on;
    else
        scatter(corr_par_x_raw(i), corr_par_y_raw(i), 'o', 'MarkerEdgeColor',[1 1 1],...
            'MarkerFaceColor',[i/length(corr_par) * 0.8, 0, 0]); hold on;
    end
end
for xLoc = 1 : length(loc_n)
    xp=r*cos(ang);  yp=r*sin(ang);
    plot(circle_loc_x(xLoc)+xp, circle_loc_y(xLoc)+yp, 'k-', 'LineWidth',3); hold on;
end
xpc = covert_ROI_r *cos(ang); ypc = covert_ROI_r *sin (ang);
plot(xpc, ypc, 'k--', 'LineWidth', 1.5); hold on;
%         xpo = overt_ROI_r *cos(ang); ypo = overt_ROI_r *sin (ang);
%         plot(xpo, ypo, 'k--', 'LineWidth', 1.5); hold on;

xlim([min(-1 * max(abs(corr_par_x_raw)), -1 * res(1)/2), max(max(abs(corr_par_x_raw)), res(1)/2)]);
ylim([min(-1 * max(abs(corr_par_y_raw)), -1 * res(2)/2), max(max(abs(corr_par_y_raw)), res(2)/2)]);
drawaxis(gca, 'x', 0, 'movelabel', 1); drawaxis(gca, 'y', 0, 'movelabel', 1);
title(sprintf('Sub%s',ID)); box('off'); 
print('-dpdf', fullfile(EyeTrackerDir, sprintf('Sub%s_EyeTracker_scatter', ID)));


%% Plot scatter of clean data
figure; hold on
for xType = 1 : length(conditionNames)
    for xLoc = 1 : length(loc_n)
        x = eye_loc_x_clean{xType, xLoc};
        y = eye_loc_y_clean{xType, xLoc};
        scatter(x, y, 'o', 'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',marker_color{xType}); hold on;
    end
end
for xLoc = 1 : length(loc_n)
    xp=r*cos(ang);  yp=r*sin(ang);
    plot(circle_loc_x(xLoc)+xp, circle_loc_y(xLoc)+yp, 'k-', 'LineWidth',3); hold on;
end
xpc = covert_ROI_r *cos(ang); ypc = covert_ROI_r *sin (ang);
plot(xpc, ypc, 'k--', 'LineWidth', 1.5); hold on;
%         xpo = overt_ROI_r *cos(ang); ypo = overt_ROI_r *sin (ang);
%         plot(xpo, ypo, 'k--', 'LineWidth', 1.5); hold on;
t = (eye_movement_outlier == 1) + (ACC == 0);
text(600,500, sprintf('%d Trial removed', sum(t > 0)));

xlim([-1 * res(1)/2, res(1)/2]); ylim([-1 * res(2)/2, res(2)/2]);
drawaxis(gca, 'x', 0, 'movelabel', 1); drawaxis(gca, 'y', 0, 'movelabel', 1);
title(sprintf('Sub%s',ID)); box('off'); 
print('-dpdf', fullfile(EyeTrackerDir, sprintf('Sub%s_EyeTracker_Clean_scatter', ID)));

% Plot averaged eye location
figure; hold on
for i = 1 : length(conditionNames)
    x = eye_loc_x_mean(i,:); y = eye_loc_y_mean(i,:);
    x = [x, eye_loc_x_mean(i,1)]; y = [y, eye_loc_y_mean(i,1)];
    plot(x,y, line_type{i}, 'MarkerFaceColor', marker_color{i},'LineWidth',1.5); hold on
    clear x y
end

for xLoc = 1 : length(loc_n)
    xp=r*cos(ang);  yp=r*sin(ang);
    plot(circle_loc_x(xLoc)+xp, circle_loc_y(xLoc)+yp, 'k-', 'LineWidth',3); hold on;
end
xpc = covert_ROI_r *cos(ang); ypc = covert_ROI_r *sin (ang);
plot(xpc, ypc, 'k--', 'LineWidth', 1.5); hold on;
%         xpo = overt_ROI_r *cos(ang); ypo = overt_ROI_r *sin (ang);
%         plot(xpo, ypo, 'k--', 'LineWidth', 1.5); hold on;

xlim([-1 * res(1)/2, res(1)/2]); ylim([-1 * res(2)/2, res(2)/2]);
drawaxis(gca, 'x', 0, 'movelabel', 1); drawaxis(gca, 'y', 0, 'movelabel', 1);
legend(conditionNames); title(sprintf('Sub%s',ID));
box('off'); 
print('-dpdf', fullfile(EyeTrackerDir, sprintf('Sub%s_EyeTracker_results', ID)));
clear Trial_type ACC TrialCodeCue Cue_loc corr_par* eye_loc*
