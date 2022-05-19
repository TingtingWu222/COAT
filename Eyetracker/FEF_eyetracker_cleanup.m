function FEF_eyetracker_cleanup(path, ID)
    %% Parameters
    nRun = 4;
    conditionNames = {'CoEn','CoEx','OvEn','OvEx'};
    shiftVal = 2;
    nSD = 2;

    %% Loop over subject
    for xSub = 1 : nSub
        fprintf('Subject%s\n', ID{xSub});
        subDir = fullfile(path,ID{xSub});
        BehvDir = fullfile(subDir,'Behv');
        EyeTrackerDir = fullfile(subDir, 'EyeTracker');
    %         cd(subDir);
        filename = spm_select('List', BehvDir, '^*.\.xlsx');
        filename_eye = sprintf('%s_GA.mat', ID{xSub});

        %% read .xlsx file
        [~,~,rawData] = xlsread(fullfile(BehvDir,filename));

        %% Orgnize the excel file
        Label = rawData(shiftVal,:);
        Trial_type = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'type')));
        ACC = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'SlideTarget.ACC')));
        TrialCodeCue = cell2mat(rawData(shiftVal+1:end,findLabel(Label,'TrialCodeCue')));
        Cue_loc = mod(TrialCodeCue, 10);        
        loc_n = unique(Cue_loc);

        %% Load eye tracker data
        load(fullfile(EyeTrackerDir, filename_eye));

        %% Calibration using moving window
        w = 64;
        for xTrial = 1 : length(corr_par)
            if xTrial > w/2 && xTrial < length(corr_par) - w/2
               t = corr_par((xTrial - w/2) : (xTrial + w/2), 3); t = t(~isnan(t));
               corr_par_x(xTrial,1) = corr_par(xTrial, 3) - mean(t); clear t
               t = corr_par((xTrial - w/2) : (xTrial + w/2), 4); t = t(~isnan(t));               
               corr_par_y(xTrial,1) = corr_par(xTrial, 4) - mean(t); clear t
            elseif xTrial <= w/2
               t = corr_par(1 : w, 3); t = t(~isnan(t));
               corr_par_x(xTrial,1) = corr_par(xTrial, 3) - mean(t); clear t
               t = corr_par(1 : w, 4); t = t(~isnan(t));
               corr_par_y(xTrial,1) = corr_par(xTrial, 4) - mean(t); clear t
            else
               t = corr_par((length(corr_par)- w) : length(corr_par), 3); t = t(~isnan(t));
               corr_par_x(xTrial,1) = corr_par(xTrial, 3) - mean(t); clear t
               t = corr_par((length(corr_par)- w) : length(corr_par), 4); t = t(~isnan(t));               
               corr_par_y(xTrial,1) = corr_par(xTrial, 4) - mean(t); clear t               
            end            
        end
        %% Eye movement corr
        for xType = 1 : length(conditionNames)
            for xLoc = 1 : length(loc_n)
    %                 eye_loc_x{xType, xLoc} = corr_par((Trial_type == xType  & Cue_loc == xLoc & ACC == 1 ...
    %                                          & ~isnan(corr_par(:,3)) & ~isnan(corr_par(:,4))), 3);
    %                 eye_loc_y{xType, xLoc} = corr_par((Trial_type == xType  & Cue_loc == xLoc & ACC == 1 ...
    %                                          & ~isnan(corr_par(:,3)) & ~isnan(corr_par(:,4))), 4);
                eye_loc_x{xType, xLoc} = corr_par_x(Trial_type == xType  & Cue_loc == xLoc & ACC == 1 ...
                                         & ~isnan(corr_par(:,3)) & ~isnan(corr_par(:,4)));
                eye_loc_y{xType, xLoc} = corr_par_y(Trial_type == xType  & Cue_loc == xLoc & ACC == 1 ...
                                         & ~isnan(corr_par(:,3)) & ~isnan(corr_par(:,4)));

                %% Remove trials with eye movement excess 3SD of mean                
                eye_loc_x_raw_mean(xType, xLoc) = mean(eye_loc_x{xType, xLoc});
                eye_loc_x_raw_sd(xType, xLoc)   = std(eye_loc_x{xType, xLoc});                
                eye_loc_y_raw_mean(xType, xLoc) = mean(eye_loc_y{xType, xLoc});
                eye_loc_y_raw_sd(xType, xLoc)   = std(eye_loc_y{xType, xLoc});

                eye_loc_x_upper(xType, xLoc) = eye_loc_x_raw_mean(xType, xLoc) + nSD * eye_loc_x_raw_sd(xType, xLoc);
                eye_loc_x_lower(xType, xLoc) = eye_loc_x_raw_mean(xType, xLoc) - nSD * eye_loc_x_raw_sd(xType, xLoc);
                eye_loc_y_upper(xType, xLoc) = eye_loc_y_raw_mean(xType, xLoc) + nSD * eye_loc_y_raw_sd(xType, xLoc);
                eye_loc_y_lower(xType, xLoc) = eye_loc_y_raw_mean(xType, xLoc) - nSD * eye_loc_y_raw_sd(xType, xLoc);

                x_t = eye_loc_x{xType, xLoc};
                y_t = eye_loc_y{xType, xLoc};
                v = (x_t >= eye_loc_x_lower(xType, xLoc) & x_t <= eye_loc_x_upper(xType, xLoc) & ... 
                     y_t >= eye_loc_y_lower(xType, xLoc) & y_t <= eye_loc_y_upper(xType, xLoc));
                eye_loc_x_clean{xType, xLoc} =  eye_loc_x{xType, xLoc}(v);
                eye_loc_y_clean{xType, xLoc} =  eye_loc_y{xType, xLoc}(v);
                clear x_t y_t

                eye_loc_x_clean_mean(xType, xLoc) = mean(eye_loc_x_clean{xType, xLoc});
                eye_loc_y_clean_mean(xType, xLoc) = mean(eye_loc_y_clean{xType, xLoc});

            end

            %% Move the original opint to the centroid of eye location
            eye_loc_x_cent = mean(mean(eye_loc_x_clean_mean));
            eye_loc_y_cent = mean(mean(eye_loc_y_clean_mean));
            eye_loc_x_mean = eye_loc_x_clean_mean - eye_loc_x_cent;
            eye_loc_y_mean = eye_loc_y_clean_mean - eye_loc_y_cent;

            save(fullfile(EyeTrackerDir, sprintf('Sub%s_EyeTracker_results.mat', ID{xSub})),...
                'eye_loc_*', 'eye_loc_y', 'eye_loc_x_mean', 'eye_loc_y_mean');

            %% plot scatter of raw data
            close all; 
            figure; hold on
            for i = 1 : length(corr_par_x);
                if Trial_type(i) == 1 || Trial_type(i) == 2
                    scatter(corr_par_x(i), corr_par_y(i), 'o', 'MarkerEdgeColor',[1 1 1],...
                        'MarkerFaceColor',[0, 0, i/length(corr_par) * 0.8]); hold on;
                else
                    scatter(corr_par_x(i), corr_par_y(i), 'o', 'MarkerEdgeColor',[1 1 1],...
                        'MarkerFaceColor',[i/length(corr_par) * 0.8, 0, 0]); hold on;
                end
            end
            circle_loc_x = [0,  13, 18,  13,  0, -13, -18, -13]/100 * 1600;
            circle_loc_y = [33, 23,  0, -23, -33, -23,  0,  23]/100 * 1200;
            r = 1200 * 0.13/2;        
            for xLoc = 1 : length(loc_n)
                ang=0:0.01:2*pi; 
                xp=r*cos(ang);
                yp=r*sin(ang);
                plot(circle_loc_x(xLoc)+xp, circle_loc_y(xLoc)+yp, 'k-', 'LineWidth',3); hold on;
            end
            xlim([min(-1 * max(abs(corr_par_x)), -800), max(max(abs(corr_par_x)),800)]);
            ylim([min(-1 * max(abs(corr_par_y)), -600), max(max(abs(corr_par_y)),600)]);
            drawaxis(gca, 'x', 0, 'movelabel', 1); drawaxis(gca, 'y', 0, 'movelabel', 1);
            title(sprintf('Sub%s',ID{xSub})); box('off'); 
            print('-dpdf', fullfile(EyeTrackerDir, sprintf('Sub%s_EyeTracker_scatter', ID{xSub})));


            %% Plot scatter of clean data
            figure; hold on
            marker_color = {'b', 'b', 'r', 'r'};
            for xType = 1 : length(conditionNames)
                for xLoc = 1 : length(loc_n)
                    x = eye_loc_x_clean{xType, xLoc} - eye_loc_x_cent;
                    y = eye_loc_y_clean{xType, xLoc} - eye_loc_x_cent;
                    scatter(x, y, 'o', 'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',marker_color{xType}); hold on;
                end
            end
            for xLoc = 1 : length(loc_n)
                ang=0:0.01:2*pi; 
                xp=r*cos(ang);
                yp=r*sin(ang);
                plot(circle_loc_x(xLoc)+xp, circle_loc_y(xLoc)+yp, 'k-', 'LineWidth',3); hold on;
            end
            xlim([min(-1 * max(abs(corr_par_x)), -800), max(max(abs(corr_par_x)),800)]);
            ylim([min(-1 * max(abs(corr_par_y)), -600), max(max(abs(corr_par_y)),600)]);
            drawaxis(gca, 'x', 0, 'movelabel', 1); drawaxis(gca, 'y', 0, 'movelabel', 1);
            title(sprintf('Sub%s',ID{xSub})); box('off'); 
            print('-dpdf', fullfile(EyeTrackerDir, sprintf('Sub%s_EyeTracker_Clean_scatter', ID{xSub})));

            %% Plot averaged eye location
            figure; hold on
            line_type = {'bo-', 'bo--', 'ro-', 'ro--'};
            for i = 1 : length(conditionNames)
                x = eye_loc_x_mean(i,:); y = eye_loc_y_mean(i,:);
                x = [x, eye_loc_x_mean(i,1)]; y = [y, eye_loc_y_mean(i,1)];
                plot(x,y, line_type{i}, 'MarkerFaceColor', marker_color{i},'LineWidth',1.5); hold on
                clear x y
            end

            for xLoc = 1 : length(loc_n)
                ang=0:0.01:2*pi; 
                xp=r*cos(ang);
                yp=r*sin(ang);
                plot(circle_loc_x(xLoc)+xp, circle_loc_y(xLoc)+yp, 'k-', 'LineWidth',3); hold on;
            end

            xlim([-800, 800]); ylim([-600, 600]);
            drawaxis(gca, 'x', 0, 'movelabel', 1); drawaxis(gca, 'y', 0, 'movelabel', 1);
            legend(conditionNames); title(sprintf('Sub%s',ID{xSub}));
            box('off'); 
            print('-dpdf', fullfile(EyeTrackerDir, sprintf('Sub%s_EyeTracker_results', ID{xSub})));
            clear Trial_type ACC TrialCodeCue Cue_loc coor_par* eye_loc*

         end
end

