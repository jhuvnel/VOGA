close all; clear; clc
warning off
%Load tables from experiments A, B, and C (Celia knows what these mean)
files = dir([cd filesep '*' filesep '*VOGResults.mat']);
files(contains(extractfield(files,'name'),'._')) = []; %remove un-analyzeable files


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
tab_LA = tab(2:2:end,{'File','Date','Frequency','CurrentAmp','Cycles','StimAxis','MaxAxis_R','MaxAxis_L','MaxEyeMovement_L','MaxEyeMovement_R','MaxEyeMovement_L_sd','MaxEyeMovement_R_sd','MaxVel_LL','MaxVel_LL_sd','MaxVel_RL','MaxVel_RL_sd','MaxVel_LZ','MaxVel_LZ_sd','MaxVel_RZ','MaxVel_RZ_sd','MaxVel_LR','MaxVel_LR_sd','MaxVel_RR','MaxVel_RR_sd','MaxVel','MaxVel_sd','Align_L','Align_L_sd','Align_R','Align_R_sd','Align','Align_sd'});
tab_RP = tab(1:2:end,{'File','Date','Frequency','CurrentAmp','Cycles','StimAxis','MaxAxis_R','MaxAxis_L','MaxEyeMovement_L','MaxEyeMovement_R','MaxEyeMovement_L_sd','MaxEyeMovement_R_sd','MaxVel_LL','MaxVel_LL_sd','MaxVel_RL','MaxVel_RL_sd','MaxVel_LZ','MaxVel_LZ_sd','MaxVel_RZ','MaxVel_RZ_sd','MaxVel_LR','MaxVel_LR_sd','MaxVel_RR','MaxVel_RR_sd','MaxVel','MaxVel_sd','Align_L','Align_L_sd','Align_R','Align_R_sd','Align','Align_sd'});

curr_amps = unique(tab_LA.CurrentAmp(~isnan(tab_LA.CurrentAmp)));
dates = unique(tab_LA.Date);


plotAmpFreqSweep(tab_LA,'R')
%%
plotFreqSweep(tab_LA,'R',40)

%%
plotAmpSweep(tab_LA,'R',10)

%%
plotMisalignmentSphere(tab_LA,'R',[],[])
%%

function [] = plotAmpFreqSweep(data,eye)
curr_amps = unique(data.CurrentAmp(~isnan(data.CurrentAmp)));
dates = unique(data.Date);
days = [1 15 20];
% colormap_jet = jet(length(curr_amps));
% colormap_jet = [
%     0.267, 0.004, 0.329;  % Dark Purple-Blue
%     0.282, 0.239, 0.545;  % Purple-Blue
%     0.253, 0.388, 0.620;  % Blue
%     0.208, 0.518, 0.584;  % Blue-Teal
%     0.161, 0.620, 0.522;  % Teal-Green
%     0.369, 0.678, 0.318;  % Green
%     0.639, 0.694, 0.200;  % Yellow-Green
%     0.867, 0.627, 0.188;  % Orange
%     0.706, 0.016, 0.149   % Dark Red
% ];
if dates(1) == datetime('06-Feb-2025') %Ch317
    colormap_jet = [
        0.0157 0.0157 0.3059;  % Dark Blue (jet start)
        %0.000, 0.255, 0.765;  % Blue
        0.000, 0.510, 0.980;  % Light Blue
        %0.000, 0.765, 0.765;  % Cyan (jet cyan)
        0.26, 1, 0.26;  % Green (jet green)
        0.765, 0.980, 0.000;  % Yellow-Green
        0.980, 0.765, 0.000;  % Yellow (jet yellow)
        0.980, 0.380, 0.000;  % Orange
        0.515, 0.000, 0.000   % Dark Red (jet end)
        ];
elseif dates(1) == datetime('22-Nov-2024') %Ch313
    data(data.CurrentAmp==15,:)=[];
    curr_amps = unique(data.CurrentAmp(~isnan(data.CurrentAmp)));
    colormap_jet = [
        0.0157 0.0157 0.3059;  % Dark Blue (jet start)
        0.000, 0.255, 0.765;  % Blue
        0.000, 0.510, 0.980;  % Light Blue
        0.000, 0.765, 0.765;  % Cyan (jet cyan)
        0.26, 1, 0.26;  % Green (jet green)
        0.765, 0.980, 0.000;  % Yellow-Green
        0.980, 0.765, 0.000;  % Yellow (jet yellow)
        0.980, 0.380, 0.000;  % Orange
        0.515, 0.000, 0.000   % Dark Red (jet end)
        ];
else
    colormap_jet = jet(length(curr_amps));
end

if strcmp(eye,'R')
    maxspeed_mean = 'MaxEyeMovement_R';
    maxspeed_std = 'MaxEyeMovement_R_sd';
else
end
% Create figure
figure('Units','inches','Position', [5 5 3.6 2.5]);
t=tiledlayout(1, length(dates), "TileSpacing","tight","Padding","tight");
% First pass to determine the maximum y-range needed
maxY = -Inf;
minY = Inf;

for d = 1:length(dates)
    currentDate = dates(d);
    dateData = data(data.Date == currentDate, :);

    if ~isempty(dateData)
        % Include standard deviation in y-range calculation
        maxY = max(maxY, max(table2array(dateData(:,maxspeed_mean)) + table2array(dateData(:,maxspeed_std))));
        minY = min(minY, min(table2array(dateData(:,maxspeed_mean)) - table2array(dateData(:,maxspeed_std))));
    end
end

% Define y-axis limits with a small buffer
yBuffer = 0.1 * (maxY - minY);
yLimits = [minY - yBuffer, maxY + yBuffer];

% Create the profile views
for d = 1:length(dates)
    nexttile
    currentDate = dates(d);
    dateData = data(data.Date == currentDate, :);

    % Create subplot for each date

    hold on;

    % Get unique amplitudes for this date
    date_amplitudes = unique(dateData.CurrentAmp);

    % Create a legend entries array
    legend_entries = cell(length(date_amplitudes), 1);

    % Plot a line for each amplitude
    for a = 1:length(date_amplitudes)
        amp = date_amplitudes(a);
        amp_data = dateData(dateData.CurrentAmp == amp, :);

        % Find the color index in the overall amplitude list
        color_idx = find(curr_amps == amp);

        % Sort by frequency
        [freq_sorted, sort_idx] = sort(amp_data.Frequency);
        eye_movement_sorted = amp_data(sort_idx,maxspeed_mean);
        sd_sorted = amp_data(sort_idx,maxspeed_std);
        eye_movement_sorted = table2array(eye_movement_sorted);
        sd_sorted = table2array(sd_sorted);

        % Plot the line with markers
        p(a) = plot(freq_sorted, eye_movement_sorted, '-o', 'Color', colormap_jet(color_idx,:), 'LineWidth', 1.8, 'MarkerSize',5);

        % Add error bars
        for i = 1:length(freq_sorted)
            x = freq_sorted(i);
            y = eye_movement_sorted(i);
            sd = sd_sorted(i);

            % Draw error bars (vertical lines)
            line([x x], [y-sd y+sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-', 'LineWidth', 1);

            % Draw caps on the error bars
            line([x-0.05 x+0.05], [y-sd y-sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-', 'LineWidth', 1);
            line([x-0.05 x+0.05], [y+sd y+sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-', 'LineWidth', 1);
        end

        % Prepare legend entry
        legend_entries{a} = [num2str(amp) ' uA'];
    end

    % Set consistent y-axis limits
    
    if dates(1) == datetime('06-Feb-2025') %Ch317
        ylim([-10 700])
    elseif dates(1) == datetime('22-Nov-2024') %Ch313
        ylim([-10 100])
    else
        ylim(yLimits);
    end
    xticks(freq_sorted);xlim([0 11])

    % Add labels and title
    if d==1
        ylabel('Peak eye speed (\circ/s)','FontSize',10)
        legend(p,legend_entries, 'Location', 'northwest','FontSize',8);
    elseif d==2
        yticklabels('')
        xlabel('Frequency (Hz)','FontSize',10)
    else
        yticklabels('')
    end
    grid on;

    hold off;

    ax = gca;
    ax.XAxis.FontSize = 8; 
    ax.YAxis.FontSize = 8; 
    
    xticklabels({'.5','1','2','3','8','10'})
    title(['Day ',num2str(days(d))],'FontSize',10)
    fontname(gcf,'Times')
    hold off;
end

end


function [] = plotFreqSweep(data,eye,currentamp)
col{1} = '#b3cde3';
col{2} = '#8c96c6';
col{3} = '#8856a7';
xoffset = 0.25*[-1 0 1];
data = data(data.CurrentAmp==currentamp,:);
curr_amps = unique(data.CurrentAmp(~isnan(data.CurrentAmp)));
dates = unique(data.Date);
if strcmp(eye,'R')
    maxspeed_mean = 'MaxEyeMovement_R';
    maxspeed_std = 'MaxEyeMovement_R_sd';
else
end
% Create figure
figure('Units','inches','Position', [5 5 1.8 2]);
% First pass to determine the maximum y-range needed
maxY = -Inf;
minY = Inf;

for d = 1:length(dates)
    currentDate = dates(d);
    dateData = data(data.Date == currentDate, :);

    if ~isempty(dateData)
        % Include standard deviation in y-range calculation
        maxY = max(maxY, max(table2array(dateData(:,maxspeed_mean)) + table2array(dateData(:,maxspeed_std))));
        minY = min(minY, min(table2array(dateData(:,maxspeed_mean)) - table2array(dateData(:,maxspeed_std))));
    end
end

% Define y-axis limits with a small buffer
yBuffer = 0.1 * (maxY - minY);
yLimits = [minY - yBuffer, maxY + yBuffer];

% Create the profile views
for d = 1:length(dates)
    currentDate = dates(d);
    dateData = data(data.Date == currentDate, :);

    % Create subplot for each date

    hold on;

    % Get unique amplitudes for this date
    date_amplitudes = unique(dateData.CurrentAmp);

    % Plot a line for each amplitude
    for a = 1:length(date_amplitudes)
        amp = date_amplitudes(a);
        amp_data = dateData(dateData.CurrentAmp == amp, :);

        % Find the color index in the overall amplitude list
        color_idx = find(curr_amps == amp);

        % Sort by frequency
        [freq_sorted, sort_idx] = sort(amp_data.Frequency);
        eye_movement_sorted = amp_data(sort_idx,maxspeed_mean);
        sd_sorted = amp_data(sort_idx,maxspeed_std);
        eye_movement_sorted = table2array(eye_movement_sorted);
        sd_sorted = table2array(sd_sorted);

        % Plot the line with markers
        p(d) = plot(freq_sorted+xoffset(d), eye_movement_sorted, '-o', 'Color', col{d}, 'LineWidth', 1.8, 'MarkerSize', 4);

        % Add error bars
        for i = 1:length(freq_sorted)
            x = freq_sorted(i)+xoffset(d);
            y = eye_movement_sorted(i);
            sd = sd_sorted(i);

            % Draw error bars (vertical lines)
            line([x x], [y-sd y+sd], 'Color', col{d}, 'LineStyle', '-', 'LineWidth', 1);

            % Draw caps on the error bars
            line([x-0.05 x+0.05], [y-sd y-sd], 'Color', col{d}, 'LineStyle', '-', 'LineWidth', 1);
            line([x-0.05 x+0.05], [y+sd y+sd], 'Color', col{d}, 'LineStyle', '-', 'LineWidth', 1);
        end

        
    end
    
    % Prepare legend entry
    legend_entries{d} = datestr(dates(d));
    
end
% Set consistent y-axis limits
ylim(yLimits);
ylim([0 700])
xticks(freq_sorted);xlim([0 11])

% Add labels and title
fontsize(gcf,8,'points');
xlabel('Frequency (Hz)','FontSize',10)
ylabel('Peak eye speed (\circ/s)','FontSize',10)
yticks([0 200 400 600])
yticklabels('')
xticklabels({'.5','1','2','3','8','10'})
grid on;
% legend(p,{'Day 1','Day 15','Day 20'},'Location','northwest','FontSize',10)
title(['Amplitude: ' num2str(currentamp) ' \muA'],'FontSize',10)
fontname(gcf,'Times')
hold off;
end


function [] = plotAmpSweep(data,eye,freq_plot)
col{1} = '#b3cde3';
col{2} = '#8c96c6';
col{3} = '#8856a7';
xoffset = 1*[-1 0 1];
data = data(data.Frequency==freq_plot,:);
curr_amps = unique(data.CurrentAmp(~isnan(data.CurrentAmp)));
dates = unique(data.Date);
if strcmp(eye,'R')
    maxspeed_mean = 'MaxEyeMovement_R';
    maxspeed_std = 'MaxEyeMovement_R_sd';
else
end
% Create figure
figure('Units','inches','Position', [5 5 1.96 2]);
% First pass to determine the maximum y-range needed
maxY = -Inf;
minY = Inf;

for d = 1:length(dates)
    currentDate = dates(d);
    dateData = data(data.Date == currentDate, :);

    if ~isempty(dateData)
        % Include standard deviation in y-range calculation
        maxY = max(maxY, max(table2array(dateData(:,maxspeed_mean)) + table2array(dateData(:,maxspeed_std))));
        minY = min(minY, min(table2array(dateData(:,maxspeed_mean)) - table2array(dateData(:,maxspeed_std))));
    end
end

% Define y-axis limits with a small buffer
yBuffer = 0.1 * (maxY - minY);
yLimits = [minY - yBuffer, maxY + yBuffer];

% Create the profile views
for d = 1:length(dates)
    currentDate = dates(d);
    dateData = data(data.Date == currentDate, :);

    % Create subplot for each date

    hold on;

    % Get unique amplitudes for this date
    date_freq = unique(dateData.Frequency);

    % Plot a line for each amplitude
    for a = 1:length(date_freq)
        freq = date_freq(a);
        freq_data = dateData(dateData.Frequency == freq, :);

        % Find the color index in the overall amplitude list
        color_idx = find(curr_amps == freq);

        % Sort by frequency
        [amp_sorted, sort_idx] = sort(freq_data.CurrentAmp);
        eye_movement_sorted = freq_data(sort_idx,maxspeed_mean);
        sd_sorted = freq_data(sort_idx,maxspeed_std);
        eye_movement_sorted = table2array(eye_movement_sorted);
        sd_sorted = table2array(sd_sorted);

        % Plot the line with markers
        p(d) = plot(amp_sorted+xoffset(d), eye_movement_sorted, '-o', 'Color', col{d}, 'LineWidth', 1.8, 'MarkerSize', 4);

        % Add error bars
        for i = 1:length(amp_sorted)
            x = amp_sorted(i)+xoffset(d);
            y = eye_movement_sorted(i);
            sd = sd_sorted(i);

            % Draw error bars (vertical lines)
            line([x x], [y-sd y+sd], 'Color', col{d}, 'LineStyle', '-', 'LineWidth', 1);

            % Draw caps on the error bars
            line([x-0.05 x+0.05], [y-sd y-sd], 'Color', col{d}, 'LineStyle', '-', 'LineWidth', 1);
            line([x-0.05 x+0.05], [y+sd y+sd], 'Color', col{d}, 'LineStyle', '-', 'LineWidth', 1);
        end

        
    end
    
    % Prepare legend entry
    legend_entries{d} = datestr(dates(d));
    
end
% Set consistent y-axis limits
ylim(yLimits);
ylim([0 700])
xticks(amp_sorted);xlim([0 max(amp_sorted)+5])

% Add labels and title
fontsize(gcf,8,'points');
xlabel('Amplitude (\muA)','FontSize',10)
ylabel('Peak eye speed (\circ/s)','FontSize',10)
% yticklabels('')
yticks([0 200 400 600])
grid on;
legend(p,{'Day 1','Day 15','Day 20'},'Location','northwest','FontSize',10)
title(['Frequency: ' num2str(freq_plot) ' Hz'],'FontSize',10)
fontname(gcf,'Times')
hold off;
end

function plotMisalignmentSphere(data,eye,amp,freq)
if ~isempty(amp)
    data = data(data.CurrentAmp==amp,:);
end
if ~isempty(freq)
    data = data(data.Frequency==freq,:);
end
dates = unique(data.Date);
if strcmp(eye,'R')
    maxaxis = 'MaxAxis_R';
else
    maxaxis = 'MaxAxis_L';
end
figure('Position', [100, 100, 500, 500], 'Color', 'white');
hold on
box off
% Draw a unit sphere
[X, Y, Z] = sphere(20);
s = surf(X, Y, Z, 'FaceAlpha', 1, 'EdgeAlpha', 0.15, 'FaceColor', [1 1 1]);
axis equal
axis off;
% grid on;
view(115,-5);
view(101.6875,20.6641);
% cameraproj
% Add axis lines
line_length=2;
larp = norm(line_length)*[1 0 0] / norm([1 0 0]);
ralp = norm(line_length)*[0 1 0] / norm([0 1 0]);
x = norm(line_length)*[1 1 0] / norm([1 1 0]);
y = norm(line_length)*[-1 1 0] / norm([-1 1 0]);
quiver3(0, 0, 0, larp(1), larp(2), larp(3), 1, 'Color','#07a621', 'LineWidth', 3,'ShowArrowHead','off'); % LARP
quiver3(0, 0, 0, ralp(1), ralp(2), ralp(3), 1, 'Color','#0515fa', 'LineWidth', 3,'ShowArrowHead','off'); % RALP
quiver3(0, 0, 0, x(1), x(2), x(3), 1, 'Color','k', 'LineWidth', 3,'ShowArrowHead','off'); % X
quiver3(0, 0, 0, y(1), y(2), y(3), 1, 'Color','#ebc634', 'LineWidth', 3,'ShowArrowHead','off'); % Y
quiver3(0, 0, 0, 0, 0, line_length, 1, 'Color','#f70519', 'LineWidth', 3,'ShowArrowHead','off'); % Z
for i = 1:height(data)
    if data.Date(i)==dates(1)
        col = '#b3cde3';
    elseif data.Date(i)==dates(2)
        col = '#8c96c6';
    elseif data.Date(i)==dates(3)
        col = '#8856a7';
    end
    mean_axis = mean(cell2mat(table2array(data(i,maxaxis))));
    mean_axis_norm = 1.05*mean_axis/norm(mean_axis);
    scatter3(mean_axis_norm(1), mean_axis_norm(2), mean_axis_norm(3), 30,'o', ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'none','MarkerFaceAlpha',0.6);
end
set(gcf,'Renderer','painters')
end