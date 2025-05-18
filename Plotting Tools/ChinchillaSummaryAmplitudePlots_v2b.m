close all; clear; clc
%Load tables from experiments A, B, and C (Celia knows what these mean)
files = dir([cd filesep '*' filesep '*VOGResults.mat']);

tab = [];
for i = 1:length(files)
    load([files(i).folder filesep files(i).name],'all_results');
    for j = 1:height(all_results)
        if all_results{j,'CurrentAmp'} > 9
            all_results{j,'CurrentAmp'} = 5*round(all_results{j,'CurrentAmp'}/5);
        end
    end
    tab = [tab; all_results];
end

%Sort table
tab = sortrows(sortrows(sortrows(tab,'Frequency','ascend'),'CurrentAmp','ascend'),'Date','ascend');
tab_LA = tab(2:2:end,{'File','Date','Frequency','CurrentAmp','Cycles','StimAxis','MaxAxis_L','MaxAxis_R','MaxEyeMovement_L','MaxEyeMovement_R','MaxEyeMovement_L_sd','MaxEyeMovement_R_sd','MaxVel_LL','MaxVel_LL_sd','MaxVel_RL','MaxVel_RL_sd','MaxVel_LZ','MaxVel_LZ_sd','MaxVel_RZ','MaxVel_RZ_sd','MaxVel_LR','MaxVel_LR_sd','MaxVel_RR','MaxVel_RR_sd','MaxVel','MaxVel_sd','Align_L','Align_L_sd','Align_R','Align_R_sd','Align','Align_sd'});
tab_RP = tab(1:2:end,{'File','Date','Frequency','CurrentAmp','Cycles','StimAxis','MaxAxis_L','MaxAxis_R','MaxEyeMovement_L','MaxEyeMovement_R','MaxEyeMovement_L_sd','MaxEyeMovement_R_sd','MaxVel_LL','MaxVel_LL_sd','MaxVel_RL','MaxVel_RL_sd','MaxVel_LZ','MaxVel_LZ_sd','MaxVel_RZ','MaxVel_RZ_sd','MaxVel_LR','MaxVel_LR_sd','MaxVel_RR','MaxVel_RR_sd','MaxVel','MaxVel_sd','Align_L','Align_L_sd','Align_R','Align_R_sd','Align','Align_sd'});
curr_amps = unique(tab_LA.CurrentAmp(~isnan(tab_LA.CurrentAmp)));
%% Calculate misalignment 1D
plotfromday0 = 1;
results = calc_misalignment_chinch(tab_LA,plotfromday0);
misalign_LA_L = results.L_misalignment_table;
misalign_LA_R = results.R_misalignment_table;

freqs=1;

figure('Position', [100, 100, 500, 500], 'Color', 'white');
hold on
box on

x = unique(misalign_LA_R.CurrentAmp);
% Left eye
y = [results.MeanAngularDeviation_R(:,1)';
    results.MeanAngularDeviation_R(:,2)';
    results.MeanAngularDeviation_R(:,3)']';
y_sd = [results.StdAngularDeviation_R(:,1)';
    results.StdAngularDeviation_R(:,2)';
    results.StdAngularDeviation_R(:,3)']';



y_sd_high = y + y_sd; 
y_sd_low = y- y_sd;


b = bar(x,y);
b(1).FaceColor = '#b3cde3';
b(2).FaceColor = '#8c96c6';
b(3).FaceColor = '#8856a7';

hold on 
% Calculate the number of groups and number of bars in each group 
[ngroups,nbars] = size(y);
% Get the x coordinate of the bars 
xeb = nan(nbars, ngroups); 
for i = 1:nbars 
    xeb(i,:) = b(i).XEndPoints; 
end
% Plot the errorbars 
errorbar(xeb',y,y_sd,'k','linestyle','none'); 
hold off 

legend({'Day 0','Day 14','Day 19'},'Location','northwest', 'FontSize', 10)
title(['Left Anterior Stimulation at ',num2str(freqs(1)),'Hz'], 'FontSize', 14, 'FontWeight', 'bold')
xlim([0 max(curr_amps)+5])
%     ylim([0 max(tab_LA.MaxVel)+100])
xticks(curr_amps)
grid on
% sgtitle('Maximum velocity')
xlabel('Stimulation Amplitude (\muA)', 'FontSize', 12, 'FontWeight', 'bold')
if plotfromday0
    ylabel('Misalignment relative to Mean(Day 0) (deg)', 'FontSize', 12, 'FontWeight', 'bold')
else
    ylabel('Misalignment relative to LA (deg)', 'FontSize', 12, 'FontWeight', 'bold')
end

%% Calculate misalignment sphere plot
figure('Position', [100, 100, 500, 500], 'Color', 'white');
hold on
box on
% Draw a unit sphere
[X, Y, Z] = sphere(20);
s = surf(X, Y, Z, 'FaceAlpha', 1, 'EdgeAlpha', 0.15, 'FaceColor', [1 1 1]);
axis equal;
% grid on;
view(117.8,0);
% Add axis lines
line_length=5;
line([-line_length line_length], [0 0], [0 0], 'Color', 'g', 'LineStyle', '-', 'LineWidth', 3);
line([0 0], [-line_length line_length], [0 0], 'Color', 'b', 'LineStyle', '-', 'LineWidth', 3);
line([0 0], [0 0], [-line_length line_length], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 3);
for i = 1:height(tab_LA)
    if tab_LA.Date(i)==datetime(2025,02,06)
        col = '#b3cde3';
    elseif tab_LA.Date(i)==datetime(2025,02,21)
        col = '#8c96c6';
    elseif tab_LA.Date(i)==datetime(2025,02,26)
        col = '#8856a7';
    end
    mean_axis = mean(tab_LA.MaxAxis_R{i});
    mean_axis_norm = 1.05*mean_axis/norm(mean_axis);
    plot3(mean_axis_norm(1), mean_axis_norm(2), mean_axis_norm(3), 'o', ...
        'MarkerSize',7, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none');
end



%% Maximum L eye movement barplot
figure('Position', [100, 100, 500, 500], 'Color', 'white');
hold on
box on
freqs = 1;
curr_amps = unique(tab_LA.CurrentAmp(~isnan(tab_LA.CurrentAmp)));
exps = unique(tab_LA.Date);

x = unique(tab_LA.CurrentAmp);
% Left eye
y = [table2array(tab_LA(1:7,'MaxEyeMovement_R'))';
    table2array(tab_LA(8:14,'MaxEyeMovement_R'))';
    table2array(tab_LA(15:21,'MaxEyeMovement_R'))';]';
y_sd = [table2array(tab_LA(1:7,'MaxEyeMovement_R_sd'))';
    table2array(tab_LA(8:14,'MaxEyeMovement_R_sd'))';
    table2array(tab_LA(15:21,'MaxEyeMovement_R_sd'))']';
% Right eye
% y = [table2array(tab_LA(1:9,'MaxEyeMovement_R'))';
%     table2array(tab_LA([10:15 17:19],'MaxEyeMovement_R'))';
%     table2array(tab_LA([20:25 27:29],'MaxEyeMovement_R'))';]';
% y_sd = [table2array(tab_LA(1:9,'MaxEyeMovement_R_sd'))';
%     table2array(tab_LA([10:15 17:19],'MaxEyeMovement_R_sd'))';
%     table2array(tab_LA([20:25 27:29],'MaxEyeMovement_R_sd'))']';
% x(7) =[];


y_sd_high = y + y_sd; 
y_sd_low = y- y_sd;


b = bar(x,y);
b(1).FaceColor = '#b3cde3';
b(2).FaceColor = '#8c96c6';
b(3).FaceColor = '#8856a7';

hold on 
% Calculate the number of groups and number of bars in each group 
[ngroups,nbars] = size(y);
% Get the x coordinate of the bars 
xeb = nan(nbars, ngroups); 
for i = 1:nbars 
    xeb(i,:) = b(i).XEndPoints; 
end
% Plot the errorbars 
errorbar(xeb',y,y_sd,'k','linestyle','none'); 
hold off 

legend({'Day 0','Day 14','Day 19'},'Location','northwest', 'FontSize', 10)
title(['Left Anterior Stimulation at ',num2str(freqs(1)),'Hz'], 'FontSize', 14, 'FontWeight', 'bold')
xlim([0 max(curr_amps)+5])
%     ylim([0 max(tab_LA.MaxVel)+100])
xticks(curr_amps)
grid on
% sgtitle('Maximum velocity')
xlabel('Stimulation Amplitude (\muA)', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Maximum Eye Speed (dps)', 'FontSize', 12, 'FontWeight', 'bold')