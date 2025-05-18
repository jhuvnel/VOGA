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
tab_LA = tab(2:2:end,{'File','Date','Frequency','CurrentAmp','Cycles','StimAxis','MaxEyeMovement_L','MaxEyeMovement_R','MaxEyeMovement_L_sd','MaxEyeMovement_R_sd','MaxVel_LL','MaxVel_LL_sd','MaxVel_RL','MaxVel_RL_sd','MaxVel_LZ','MaxVel_LZ_sd','MaxVel_RZ','MaxVel_RZ_sd','MaxVel_LR','MaxVel_LR_sd','MaxVel_RR','MaxVel_RR_sd','MaxVel','MaxVel_sd','Align_L','Align_L_sd','Align_R','Align_R_sd','Align','Align_sd'});
tab_RP = tab(1:2:end,{'File','Date','Frequency','CurrentAmp','Cycles','StimAxis','MaxEyeMovement_L','MaxEyeMovement_R','MaxEyeMovement_L_sd','MaxEyeMovement_R_sd','MaxVel_LL','MaxVel_LL_sd','MaxVel_RL','MaxVel_RL_sd','MaxVel_LZ','MaxVel_LZ_sd','MaxVel_RZ','MaxVel_RZ_sd','MaxVel_LR','MaxVel_LR_sd','MaxVel_RR','MaxVel_RR_sd','MaxVel','MaxVel_sd','Align_L','Align_L_sd','Align_R','Align_R_sd','Align','Align_sd'});

curr_amps = unique(tab_LA.CurrentAmp(~isnan(tab_LA.CurrentAmp)));

%% LA L eye - Profile View

% Extract unique dates from the table
dates = unique(tab_LA.Date);

% Create figure
figure('Position', [100, 100, 1200, 400]);

% First pass to determine the maximum y-range needed
maxY = -Inf;
minY = Inf;

for d = 1:length(dates)
    currentDate = dates(d);
    dateData = tab_LA(tab_LA.Date == currentDate, :);
    
    if ~isempty(dateData)
        % Include standard deviation in y-range calculation
        maxY = max(maxY, max(dateData.MaxEyeMovement_L + dateData.MaxEyeMovement_L_sd));
        minY = min(minY, min(dateData.MaxEyeMovement_L - dateData.MaxEyeMovement_L_sd));
    end
end

% Define y-axis limits with a small buffer
yBuffer = 0.1 * (maxY - minY);
yLimits = [minY - yBuffer, maxY + yBuffer];

% Get unique current amplitudes for color mapping
all_amplitudes = unique(tab_LA.CurrentAmp);
colormap_jet = jet(length(all_amplitudes));

% Create the profile views
for d = 1:length(dates)
    currentDate = dates(d);
    dateData = tab_LA(tab_LA.Date == currentDate, :);
    
    % Create subplot for each date
    subplot(1, length(dates), d);
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
        color_idx = find(all_amplitudes == amp);
        
        % Sort by frequency
        [freq_sorted, sort_idx] = sort(amp_data.Frequency);
        eye_movement_sorted = amp_data.MaxEyeMovement_L(sort_idx);
        sd_sorted = amp_data.MaxEyeMovement_L_sd(sort_idx);
        
        % Plot the line with markers
        p = plot(freq_sorted, eye_movement_sorted, '-o', 'Color', colormap_jet(color_idx,:), 'LineWidth', 2, 'MarkerSize', 6);
        
        % Add error bars
        for i = 1:length(freq_sorted)
            x = freq_sorted(i);
            y = eye_movement_sorted(i);
            sd = sd_sorted(i);
            
            % Draw error bars (vertical lines)
            line([x x], [y-sd y+sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-');
            
            % Draw caps on the error bars
            line([x-0.05 x+0.05], [y-sd y-sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-');
            line([x-0.05 x+0.05], [y+sd y+sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-');
        end
        
        % Prepare legend entry
        legend_entries{a} = [num2str(amp) ' uA'];
    end
    
    % Set consistent y-axis limits
    ylim(yLimits);
    
    % Add labels and title
    xlabel('Frequency (Hz)');
    ylabel('MaxEyeMovement LE (deg)');
    title(['Date: ', datestr(currentDate)]);
    grid on;
    
    % Add legend
    legend(legend_entries, 'Location', 'best');
    hold off;
end

% Add a main title
sgtitle('Profile Views of MaxEyeMovement LE by Frequency for LA Stim (with SD error bars)');

%% LA R eye - Profile View

% Create figure
figure('Position', [100, 100, 1200, 400]);

% First pass to determine the maximum y-range needed
maxY = -Inf;
minY = Inf;

for d = 1:length(dates)
    currentDate = dates(d);
    dateData = tab_LA(tab_LA.Date == currentDate, :);
    
    if ~isempty(dateData)
        % Include standard deviation in y-range calculation
        maxY = max(maxY, max(dateData.MaxEyeMovement_R + dateData.MaxEyeMovement_R_sd));
        minY = min(minY, min(dateData.MaxEyeMovement_R - dateData.MaxEyeMovement_R_sd));
    end
end

% Define y-axis limits with a small buffer
yBuffer = 0.1 * (maxY - minY);
yLimits = [minY - yBuffer, maxY + yBuffer];

% Create the profile views
for d = 1:length(dates)
    currentDate = dates(d);
    dateData = tab_LA(tab_LA.Date == currentDate, :);
    
    % Create subplot for each date
    subplot(1, length(dates), d);
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
        color_idx = find(all_amplitudes == amp);
        
        % Sort by frequency
        [freq_sorted, sort_idx] = sort(amp_data.Frequency);
        eye_movement_sorted = amp_data.MaxEyeMovement_R(sort_idx);
        sd_sorted = amp_data.MaxEyeMovement_R_sd(sort_idx);
        
        % Plot the line with markers
        p = plot(freq_sorted, eye_movement_sorted, '-o', 'Color', colormap_jet(color_idx,:), 'LineWidth', 2, 'MarkerSize', 6);
        
        % Add error bars
        for i = 1:length(freq_sorted)
            x = freq_sorted(i);
            y = eye_movement_sorted(i);
            sd = sd_sorted(i);
            
            % Draw error bars (vertical lines)
            line([x x], [y-sd y+sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-');
            
            % Draw caps on the error bars
            line([x-0.05 x+0.05], [y-sd y-sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-');
            line([x-0.05 x+0.05], [y+sd y+sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-');
        end
        
        % Prepare legend entry
        legend_entries{a} = [num2str(amp) ' uA'];
    end
    
    % Set consistent y-axis limits
    ylim(yLimits);
    
    % Add labels and title
    xlabel('Frequency (Hz)');
    ylabel('MaxEyeMovement RE (deg)');
    title(['Date: ', datestr(currentDate)]);
    grid on;
    
    % Add legend
    legend(legend_entries, 'Location', 'best');
    hold off;
end

% Add a main title
sgtitle('Profile Views of MaxEyeMovement RE by Frequency for LA Stim (with SD error bars)');

%% RP L eye - Profile View

% Extract unique dates from the table
dates = unique(tab_RP.Date);

% Create figure
figure('Position', [100, 100, 1200, 400]);

% First pass to determine the maximum y-range needed
maxY = -Inf;
minY = Inf;

for d = 1:length(dates)
    currentDate = dates(d);
    dateData = tab_RP(tab_RP.Date == currentDate, :);
    
    if ~isempty(dateData)
        % Include standard deviation in y-range calculation
        maxY = max(maxY, max(dateData.MaxEyeMovement_L + dateData.MaxEyeMovement_L_sd));
        minY = min(minY, min(dateData.MaxEyeMovement_L - dateData.MaxEyeMovement_L_sd));
    end
end

% Define y-axis limits with a small buffer
yBuffer = 0.1 * (maxY - minY);
yLimits = [minY - yBuffer, maxY + yBuffer];

% Get unique current amplitudes for color mapping
all_amplitudes = unique(tab_RP.CurrentAmp);
colormap_jet = jet(length(all_amplitudes));

% Create the profile views
for d = 1:length(dates)
    currentDate = dates(d);
    dateData = tab_RP(tab_RP.Date == currentDate, :);
    
    % Create subplot for each date
    subplot(1, length(dates), d);
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
        color_idx = find(all_amplitudes == amp);
        
        % Sort by frequency
        [freq_sorted, sort_idx] = sort(amp_data.Frequency);
        eye_movement_sorted = amp_data.MaxEyeMovement_L(sort_idx);
        sd_sorted = amp_data.MaxEyeMovement_L_sd(sort_idx);
        
        % Plot the line with markers
        p = plot(freq_sorted, eye_movement_sorted, '-o', 'Color', colormap_jet(color_idx,:), 'LineWidth', 2, 'MarkerSize', 6);
        
        % Add error bars
        for i = 1:length(freq_sorted)
            x = freq_sorted(i);
            y = eye_movement_sorted(i);
            sd = sd_sorted(i);
            
            % Draw error bars (vertical lines)
            line([x x], [y-sd y+sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-');
            
            % Draw caps on the error bars
            line([x-0.05 x+0.05], [y-sd y-sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-');
            line([x-0.05 x+0.05], [y+sd y+sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-');
        end
        
        % Prepare legend entry
        legend_entries{a} = [num2str(amp) ' uA'];
    end
    
    % Set consistent y-axis limits
    ylim(yLimits);
    
    % Add labels and title
    xlabel('Frequency (Hz)');
    ylabel('MaxEyeMovement LE (deg)');
    title(['Date: ', datestr(currentDate)]);
    grid on;
    
    % Add legend
    legend(legend_entries, 'Location', 'best');
    hold off;
end

% Add a main title
sgtitle('Profile Views of MaxEyeMovement LE by Frequency for RP Stim (with SD error bars)');

%% RP R eye - Profile View

% Create figure
figure('Position', [100, 100, 1200, 400]);

% First pass to determine the maximum y-range needed
maxY = -Inf;
minY = Inf;

for d = 1:length(dates)
    currentDate = dates(d);
    dateData = tab_RP(tab_RP.Date == currentDate, :);
    
    if ~isempty(dateData)
        % Include standard deviation in y-range calculation
        maxY = max(maxY, max(dateData.MaxEyeMovement_R + dateData.MaxEyeMovement_R_sd));
        minY = min(minY, min(dateData.MaxEyeMovement_R - dateData.MaxEyeMovement_R_sd));
    end
end

% Define y-axis limits with a small buffer
yBuffer = 0.1 * (maxY - minY);
yLimits = [minY - yBuffer, maxY + yBuffer];

% Create the profile views
for d = 1:length(dates)
    currentDate = dates(d);
    dateData = tab_RP(tab_RP.Date == currentDate, :);
    
    % Create subplot for each date
    subplot(1, length(dates), d);
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
        color_idx = find(all_amplitudes == amp);
        
        % Sort by frequency
        [freq_sorted, sort_idx] = sort(amp_data.Frequency);
        eye_movement_sorted = amp_data.MaxEyeMovement_R(sort_idx);
        sd_sorted = amp_data.MaxEyeMovement_R_sd(sort_idx);
        
        % Plot the line with markers
        p = plot(freq_sorted, eye_movement_sorted, '-o', 'Color', colormap_jet(color_idx,:), 'LineWidth', 2, 'MarkerSize', 6);
        
        % Add error bars
        for i = 1:length(freq_sorted)
            x = freq_sorted(i);
            y = eye_movement_sorted(i);
            sd = sd_sorted(i);
            
            % Draw error bars (vertical lines)
            line([x x], [y-sd y+sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-');
            
            % Draw caps on the error bars
            line([x-0.05 x+0.05], [y-sd y-sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-');
            line([x-0.05 x+0.05], [y+sd y+sd], 'Color', colormap_jet(color_idx,:), 'LineStyle', '-');
        end
        
        % Prepare legend entry
        legend_entries{a} = [num2str(amp) ' uA'];
    end
    
    % Set consistent y-axis limits
    ylim(yLimits);
    
    % Add labels and title
    xlabel('Frequency (Hz)');
    ylabel('MaxEyeMovement RE (deg)');
    title(['Date: ', datestr(currentDate)]);
    grid on;
    
    % Add legend
    legend(legend_entries, 'Location', 'best');
    hold off;
end

% Add a main title
sgtitle('Profile Views of MaxEyeMovement RE by Frequency for RP Stim (with SD error bars)');