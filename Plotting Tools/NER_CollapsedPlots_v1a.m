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
        maxY = max(maxY, max(dateData.MaxEyeMovement_L));
        minY = min(minY, min(dateData.MaxEyeMovement_L));
    end
end

% Define y-axis limits with a small buffer
yBuffer = 0.1 * (maxY - minY);
yLimits = [minY - yBuffer, maxY + yBuffer];

% Get unique frequencies for color mapping
all_frequencies = unique(tab_LA.Frequency);
colormap_jet = jet(length(all_frequencies));

% Create the profile views
for d = 1:length(dates)
    currentDate = dates(d);
    dateData = tab_LA(tab_LA.Date == currentDate, :);
    
    % Create subplot for each date
    subplot(1, length(dates), d);
    hold on;
    
    % Get unique frequencies for this date
    date_frequencies = unique(dateData.Frequency);
    
    % Create a legend entries array
    legend_entries = cell(length(date_frequencies), 1);
    
    % Plot a line for each frequency
    for f = 1:length(date_frequencies)
        freq = date_frequencies(f);
        freq_data = dateData(dateData.Frequency == freq, :);
        
        % Find the color index in the overall frequency list
        color_idx = find(all_frequencies == freq);
        
        % Sort by current amplitude
        [curr_amps_sorted, sort_idx] = sort(freq_data.CurrentAmp);
        eye_movement_sorted = freq_data.MaxEyeMovement_L(sort_idx);
        
        % Plot the line with markers
        plot(curr_amps_sorted, eye_movement_sorted, '-o', 'Color', colormap_jet(color_idx,:), 'LineWidth', 2, 'MarkerSize', 6);
        
        % Prepare legend entry
        legend_entries{f} = [num2str(freq) ' Hz'];
    end
    
    % Set consistent y-axis limits
    ylim(yLimits);
    
    % Add labels and title
    xlabel('Current Amplitude');
    ylabel('MaxEyeMovement LE');
    title(['Date: ', datestr(currentDate)]);
    grid on;
    
    % Add legend
    legend(legend_entries, 'Location', 'best');
    hold off;
end

% Add a main title
sgtitle('Profile Views of MaxEyeMovement LE by Current Amplitude for LA Stim');

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
        maxY = max(maxY, max(dateData.MaxEyeMovement_R));
        minY = min(minY, min(dateData.MaxEyeMovement_R));
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
    
    % Get unique frequencies for this date
    date_frequencies = unique(dateData.Frequency);
    
    % Create a legend entries array
    legend_entries = cell(length(date_frequencies), 1);
    
    % Plot a line for each frequency
    for f = 1:length(date_frequencies)
        freq = date_frequencies(f);
        freq_data = dateData(dateData.Frequency == freq, :);
        
        % Find the color index in the overall frequency list
        color_idx = find(all_frequencies == freq);
        
        % Sort by current amplitude
        [curr_amps_sorted, sort_idx] = sort(freq_data.CurrentAmp);
        eye_movement_sorted = freq_data.MaxEyeMovement_R(sort_idx);
        
        % Plot the line with markers
        plot(curr_amps_sorted, eye_movement_sorted, '-o', 'Color', colormap_jet(color_idx,:), 'LineWidth', 2, 'MarkerSize', 6);
        
        % Prepare legend entry
        legend_entries{f} = [num2str(freq) ' Hz'];
    end
    
    % Set consistent y-axis limits
    ylim(yLimits);
    
    % Add labels and title
    xlabel('Current Amplitude');
    ylabel('MaxEyeMovement RE');
    title(['Date: ', datestr(currentDate)]);
    grid on;
    
    % Add legend
    legend(legend_entries, 'Location', 'best');
    hold off;
end

% Add a main title
sgtitle('Profile Views of MaxEyeMovement RE by Current Amplitude for LA Stim');

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
        maxY = max(maxY, max(dateData.MaxEyeMovement_L));
        minY = min(minY, min(dateData.MaxEyeMovement_L));
    end
end

% Define y-axis limits with a small buffer
yBuffer = 0.1 * (maxY - minY);
yLimits = [minY - yBuffer, maxY + yBuffer];

% Get unique frequencies for color mapping
all_frequencies = unique(tab_RP.Frequency);
colormap_jet = jet(length(all_frequencies));

% Create the profile views
for d = 1:length(dates)
    currentDate = dates(d);
    dateData = tab_RP(tab_RP.Date == currentDate, :);
    
    % Create subplot for each date
    subplot(1, length(dates), d);
    hold on;
    
    % Get unique frequencies for this date
    date_frequencies = unique(dateData.Frequency);
    
    % Create a legend entries array
    legend_entries = cell(length(date_frequencies), 1);
    
    % Plot a line for each frequency
    for f = 1:length(date_frequencies)
        freq = date_frequencies(f);
        freq_data = dateData(dateData.Frequency == freq, :);
        
        % Find the color index in the overall frequency list
        color_idx = find(all_frequencies == freq);
        
        % Sort by current amplitude
        [curr_amps_sorted, sort_idx] = sort(freq_data.CurrentAmp);
        eye_movement_sorted = freq_data.MaxEyeMovement_L(sort_idx);
        
        % Plot the line with markers
        plot(curr_amps_sorted, eye_movement_sorted, '-o', 'Color', colormap_jet(color_idx,:), 'LineWidth', 2, 'MarkerSize', 6);
        
        % Prepare legend entry
        legend_entries{f} = [num2str(freq) ' Hz'];
    end
    
    % Set consistent y-axis limits
    ylim(yLimits);
    
    % Add labels and title
    xlabel('Current Amplitude');
    ylabel('MaxEyeMovement LE');
    title(['Date: ', datestr(currentDate)]);
    grid on;
    
    % Add legend
    legend(legend_entries, 'Location', 'best');
    hold off;
end

% Add a main title
sgtitle('Profile Views of MaxEyeMovement LE by Current Amplitude for RP Stim');

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
        maxY = max(maxY, max(dateData.MaxEyeMovement_R));
        minY = min(minY, min(dateData.MaxEyeMovement_R));
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
    
    % Get unique frequencies for this date
    date_frequencies = unique(dateData.Frequency);
    
    % Create a legend entries array
    legend_entries = cell(length(date_frequencies), 1);
    
    % Plot a line for each frequency
    for f = 1:length(date_frequencies)
        freq = date_frequencies(f);
        freq_data = dateData(dateData.Frequency == freq, :);
        
        % Find the color index in the overall frequency list
        color_idx = find(all_frequencies == freq);
        
        % Sort by current amplitude
        [curr_amps_sorted, sort_idx] = sort(freq_data.CurrentAmp);
        eye_movement_sorted = freq_data.MaxEyeMovement_R(sort_idx);
        
        % Plot the line with markers
        plot(curr_amps_sorted, eye_movement_sorted, '-o', 'Color', colormap_jet(color_idx,:), 'LineWidth', 2, 'MarkerSize', 6);
        
        % Prepare legend entry
        legend_entries{f} = [num2str(freq) ' Hz'];
    end
    
    % Set consistent y-axis limits
    ylim(yLimits);
    
    % Add labels and title
    xlabel('Current Amplitude');
    ylabel('MaxEyeMovement RE');
    title(['Date: ', datestr(currentDate)]);
    grid on;
    
    % Add legend
    legend(legend_entries, 'Location', 'best');
    hold off;
end

% Add a main title
sgtitle('Profile Views of MaxEyeMovement RE by Current Amplitude for RP Stim');