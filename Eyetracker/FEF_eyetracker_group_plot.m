function FEF_eyetracker_group_plot(path, ID)
    %% Parameters
    nRun = 4;
    conditionNames = {'CoEn','CoEx','OvEn','OvEx'};
    nCond = length(conditionNames);
    nLoc = 8; 
    %% Plot parameters
    % res = [1600, 1200];
    res = [1280, 1024];
    % circle_loc_x = [0,  13, 18,  13,  0, -13, -18, -13]/100 * res(1);
    % circle_loc_y = [33, 23,  0, -23, -33, -23,  0,  23]/100 * res(2);
    circle_loc_x = [0,  13, 18,  13,  0, -13, -18, -13]/100 * res(1)* 1600/1280;
    circle_loc_y = [33, 23,  0, -23, -33, -23,  0,  23]/100 * res(2) ;%* 1200/1024;
    r = res(2) * 0.13/2;       
    covert_ROI_r = res(2)/2/4;
    overt_ROI_r = res(2)/2;
    marker_color = {'b', 'b', 'r', 'r'};
    line_type = {'b-', 'b--', 'r-', 'r--'};

    %% Loop over subject
    for xSub = 1 : nSub
             if xSub ~= 9 && xSub ~= 10 && xSub ~= 11 %&& xSub ~= 6 && xSub ~= 8  &&xSub ~= 2  
             subDir = fullfile(path,ID{xSub});
             BehvDir = fullfile(subDir,'Behv');
             EyeTrackerDir = fullfile(subDir, 'EyeTracker');
             load(fullfile(EyeTrackerDir, sprintf('Sub%s_EyeTracker_results.mat', ID{xSub})));
             for i = 1 : nCond
                 for j = 1 : nLoc
                     eye_loc_group{i, j}(xSub, 1) = eye_loc_x_mean(i, j);
                     eye_loc_group{i, j}(xSub, 2) = eye_loc_y_mean(i, j);
                 end
             end
             clear eye_loc_x_* eye_loc_y_*
         else
             for i = 1 : nCond
                 for j = 1 : nLoc
                     eye_loc_group{i, j}(xSub, :) = nan;
                 end          
             end
             end
    end

    % Mean x, y of each condition
    for i = 1 : nCond
        for j = 1 : nLoc
            eye_loc_group_valid{i, j} = eye_loc_group{i, j}(~isnan(eye_loc_group{i, j}(:,1)), :);
            eye_loc_x_group_mean(i,j) = mean(eye_loc_group_valid{i, j}(:, 1));
            eye_loc_y_group_mean(i,j) = mean(eye_loc_group_valid{i, j}(:, 2));
        end
    end

    for xLoc = 1 : nLoc
        overt_loc{1, xLoc}  = (eye_loc_group_valid{3, xLoc} + eye_loc_group_valid{4, xLoc})/2;
        overt_loc_m(1, xLoc) = mean(overt_loc{1, xLoc}(:, 1));
        overt_loc_m(2, xLoc) = mean(overt_loc{1, xLoc}(:, 2));
        overt_loc_sd(1, xLoc) = std(overt_loc{1, xLoc}(:, 1));
        overt_loc_sd(2, xLoc) = std(overt_loc{1, xLoc}(:, 2));
        covert_loc{1, xLoc}  = (eye_loc_group_valid{1, xLoc} + eye_loc_group_valid{2, xLoc})/2;
        covert_loc_m(1, xLoc) = mean(covert_loc{1, xLoc}(:, 1));
        covert_loc_m(2, xLoc) = mean(covert_loc{1, xLoc}(:, 2));
        covert_loc_sd(1, xLoc) = std(covert_loc{1, xLoc}(:, 1));
        covert_loc_sd(2, xLoc) = std(covert_loc{1, xLoc}(:, 2));
    end

    %% Plot averaged eye location
    interv = 0.01;
    figure; hold on

    ang=0:interv:(2*pi+interv);
    for xLoc = 1 : nLoc
        xp=r*cos(ang);  yp=r*sin(ang);
        plot(circle_loc_x(xLoc)+xp, circle_loc_y(xLoc)+yp, 'k-', 'LineWidth',3); hold on; 
        t = [overt_loc_m(1, xLoc) - overt_loc_sd(1, xLoc), overt_loc_m(2, xLoc) - overt_loc_sd(2, xLoc), ...
             overt_loc_sd(1, xLoc)*2,  overt_loc_sd(2, xLoc)*2];
        rectangle('Position',t, 'EdgeColor', 'none', 'FaceColor', [1 .5 .5]);
        t = [covert_loc_m(1, xLoc) - covert_loc_sd(1, xLoc), covert_loc_m(2, xLoc) - covert_loc_sd(2, xLoc), ...
            covert_loc_sd(1, xLoc)*2,  covert_loc_sd(2, xLoc)*2];
        rectangle('Position',t, 'EdgeColor', 'none', 'FaceColor', [.5 .5 1]);
    end

    for i = 3 : length(conditionNames)
         x = eye_loc_x_group_mean(i,:); y = eye_loc_y_group_mean(i,:);
        x = [x, eye_loc_x_group_mean(i,1)]; y = [y, eye_loc_y_group_mean(i,1)];
        plot(x,y, line_type{i},'LineWidth',1.5); hold on
        clear x y
    end
    xlim([-1 * res(1)/2, res(1)/2]); ylim([-1 * res(2)/2, res(2)/2]);
    drawaxis(gca, 'x', 0, 'movelabel', 1); drawaxis(gca, 'y', 0, 'movelabel', 1);
    legend(conditionNames); title(sprintf('Sub%s',ID{xSub}));
    box('off'); 
    print('-dpdf', fullfile(EyeTrackerDir, sprintf('Sub%s_EyeTracker_results', ID{xSub})));
end
