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
% tab_LA = tab(2:2:end,{'File','Date','Frequency','CurrentAmp','Cycles','StimAxis','MaxAxis_L','MaxAxis_R','MaxEyeMovement_L','MaxEyeMovement_R','MaxEyeMovement_L_sd','MaxEyeMovement_R_sd','MaxVel_LL','MaxVel_LL_sd','MaxVel_RL','MaxVel_RL_sd','MaxVel_LZ','MaxVel_LZ_sd','MaxVel_RZ','MaxVel_RZ_sd','MaxVel_LR','MaxVel_LR_sd','MaxVel_RR','MaxVel_RR_sd','MaxVel','MaxVel_sd','Align_L','Align_L_sd','Align_R','Align_R_sd','Align','Align_sd'});
% tab_RP = tab(1:2:end,{'File','Date','Frequency','CurrentAmp','Cycles','StimAxis','MaxAxis_L','MaxAxis_R','MaxEyeMovement_L','MaxEyeMovement_R','MaxEyeMovement_L_sd','MaxEyeMovement_R_sd','MaxVel_LL','MaxVel_LL_sd','MaxVel_RL','MaxVel_RL_sd','MaxVel_LZ','MaxVel_LZ_sd','MaxVel_RZ','MaxVel_RZ_sd','MaxVel_LR','MaxVel_LR_sd','MaxVel_RR','MaxVel_RR_sd','MaxVel','MaxVel_sd','Align_L','Align_L_sd','Align_R','Align_R_sd','Align','Align_sd'});
tab_LA = tab(2:2:end,{'File','Date','Frequency','CurrentAmp','Cycles','StimAxis','MaxEyeMovement_L','MaxEyeMovement_R','MaxEyeMovement_L_sd','MaxEyeMovement_R_sd','MaxVel_LL','MaxVel_LL_sd','MaxVel_RL','MaxVel_RL_sd','MaxVel_LZ','MaxVel_LZ_sd','MaxVel_RZ','MaxVel_RZ_sd','MaxVel_LR','MaxVel_LR_sd','MaxVel_RR','MaxVel_RR_sd','MaxVel','MaxVel_sd','Align_L','Align_L_sd','Align_R','Align_R_sd','Align','Align_sd'});
tab_RP = tab(1:2:end,{'File','Date','Frequency','CurrentAmp','Cycles','StimAxis','MaxEyeMovement_L','MaxEyeMovement_R','MaxEyeMovement_L_sd','MaxEyeMovement_R_sd','MaxVel_LL','MaxVel_LL_sd','MaxVel_RL','MaxVel_RL_sd','MaxVel_LZ','MaxVel_LZ_sd','MaxVel_RZ','MaxVel_RZ_sd','MaxVel_LR','MaxVel_LR_sd','MaxVel_RR','MaxVel_RR_sd','MaxVel','MaxVel_sd','Align_L','Align_L_sd','Align_R','Align_R_sd','Align','Align_sd'});

curr_amps = unique(tab_LA.CurrentAmp(~isnan(tab_LA.CurrentAmp)));
%% LA L eye

% Extract unique dates from the table
dates = unique(tab_LA.Date);

% Create figure
figure('Position', [100, 100, 1200, 400]);

% First pass to determine the maximum z-range needed
maxZ = -Inf;
minZ = Inf;

for d = 1:length(dates)
    currentDate = dates(d);
    
    % Filter data for the current date
    dateData = tab_LA(tab_LA.Date == currentDate, :);
    
    % Update global min/max for the z-axis
    if ~isempty(dateData)
        maxZ = max(maxZ, max(dateData.MaxEyeMovement_L));
        minZ = min(minZ, min(dateData.MaxEyeMovement_L));
    end
end

% Define z-axis limits with a small buffer
zBuffer = 0.05 * (maxZ - minZ);
zLimits = [minZ - zBuffer, maxZ + zBuffer];

% Second pass to create the plots with consistent z-axis scaling
for d = 1:length(dates)
    currentDate = dates(d);
    
    % Filter data for the current date
    dateData = tab_LA(tab_LA.Date == currentDate, :);
    
    % Get unique CurrentAmp and Frequency values
    currentAmps = unique(dateData.CurrentAmp);
    frequencies = unique(dateData.Frequency);
    
    % Create a grid for the surface plot
    [X, Y] = meshgrid(currentAmps, frequencies);
    Z = zeros(length(frequencies), length(currentAmps));
    
    % Fill the Z matrix with average MaxEyeMovement_L values
    for i = 1:length(currentAmps)
        for j = 1:length(frequencies)
            % Find rows matching current amplitude and frequency
            matches = dateData(dateData.CurrentAmp == currentAmps(i) & ...
                              dateData.Frequency == frequencies(j), :);
            
            % Calculate mean if there are matches
            if ~isempty(matches)
                Z(j, i) = mean(matches.MaxEyeMovement_L);
            end
        end
    end
    
    % Create subplot for each date
    subplot(1, length(dates), d);
    surf(X, Y, Z);
    
    % Set consistent z-axis limits
    zlim(zLimits);
    
    % Add labels and title
    xlabel('Current Amplitude');
    ylabel('Frequency');
    zlabel('MaxEyeMovement LE');
    title(['Date: ', datestr(currentDate)]);
    
    % Add colorbar and make the plot look nice
%     shading interp;
    colormap(jet);
    view(-10,10);
    
    % Set consistent colorbar limits across all plots
    caxis(zLimits);
end
colorbar
% Add a main title
sgtitle('Surface Plots of MaxEyeMovement LE by Date for LA Stim');

%% LA R eye

% Extract unique dates from the table
dates = unique(tab_LA.Date);

% Create figure
figure('Position', [100, 100, 1200, 400]);

% First pass to determine the maximum z-range needed
maxZ = -Inf;
minZ = Inf;

for d = 1:length(dates)
    currentDate = dates(d);
    
    % Filter data for the current date
    dateData = tab_LA(tab_LA.Date == currentDate, :);
    
    % Update global min/max for the z-axis
    if ~isempty(dateData)
        maxZ = max(maxZ, max(dateData.MaxEyeMovement_R));
        minZ = min(minZ, min(dateData.MaxEyeMovement_R));
    end
end

% Define z-axis limits with a small buffer
zBuffer = 0.05 * (maxZ - minZ);
zLimits = [minZ - zBuffer, maxZ + zBuffer];

% Second pass to create the plots with consistent z-axis scaling
for d = 1:length(dates)
    currentDate = dates(d);
    
    % Filter data for the current date
    dateData = tab_LA(tab_LA.Date == currentDate, :);
    
    % Get unique CurrentAmp and Frequency values
    currentAmps = unique(dateData.CurrentAmp);
    frequencies = unique(dateData.Frequency);
    
    % Create a grid for the surface plot
    [X, Y] = meshgrid(currentAmps, frequencies);
    Z = zeros(length(frequencies), length(currentAmps));
    
    % Fill the Z matrix with average MaxEyeMovement_L values
    for i = 1:length(currentAmps)
        for j = 1:length(frequencies)
            % Find rows matching current amplitude and frequency
            matches = dateData(dateData.CurrentAmp == currentAmps(i) & ...
                              dateData.Frequency == frequencies(j), :);
            
            % Calculate mean if there are matches
            if ~isempty(matches)
                Z(j, i) = mean(matches.MaxEyeMovement_R);
            end
        end
    end
    
    % Create subplot for each date
    subplot(1, length(dates), d);
    surf(X, Y, Z);
    
    % Set consistent z-axis limits
    zlim(zLimits);
    
    % Add labels and title
    xlabel('Current Amplitude');
    ylabel('Frequency');
    zlabel('MaxEyeMovement RE');
    title(['Date: ', datestr(currentDate)]);
%     shading interp;
    colormap(jet);
    view(-10,10);


    % Set consistent colorbar limits across all plots
    caxis(zLimits);
end
colorbar;


% Add a main title
sgtitle('Surface Plots of MaxEyeMovement RE by Date for LA Stim');


%% LA L eye

% Extract unique dates from the table
dates = unique(tab_RP.Date);

% Create figure
figure('Position', [100, 100, 1200, 400]);

% First pass to determine the maximum z-range needed
maxZ = -Inf;
minZ = Inf;

for d = 1:length(dates)
    currentDate = dates(d);
    
    % Filter data for the current date
    dateData = tab_RP(tab_RP.Date == currentDate, :);
    
    % Update global min/max for the z-axis
    if ~isempty(dateData)
        maxZ = max(maxZ, max(dateData.MaxEyeMovement_L));
        minZ = min(minZ, min(dateData.MaxEyeMovement_L));
    end
end

% Define z-axis limits with a small buffer
zBuffer = 0.05 * (maxZ - minZ);
zLimits = [minZ - zBuffer, maxZ + zBuffer];

% Second pass to create the plots with consistent z-axis scaling
for d = 1:length(dates)
    currentDate = dates(d);
    
    % Filter data for the current date
    dateData = tab_RP(tab_RP.Date == currentDate, :);
    
    % Get unique CurrentAmp and Frequency values
    currentAmps = unique(dateData.CurrentAmp);
    frequencies = unique(dateData.Frequency);
    
    % Create a grid for the surface plot
    [X, Y] = meshgrid(currentAmps, frequencies);
    Z = zeros(length(frequencies), length(currentAmps));
    
    % Fill the Z matrix with average MaxEyeMovement_L values
    for i = 1:length(currentAmps)
        for j = 1:length(frequencies)
            % Find rows matching current amplitude and frequency
            matches = dateData(dateData.CurrentAmp == currentAmps(i) & ...
                              dateData.Frequency == frequencies(j), :);
            
            % Calculate mean if there are matches
            if ~isempty(matches)
                Z(j, i) = mean(matches.MaxEyeMovement_L);
            end
        end
    end
    
    % Create subplot for each date
    subplot(1, length(dates), d);
    surf(X, Y, Z);
    
    % Set consistent z-axis limits
    zlim(zLimits);
    
    % Add labels and title
    xlabel('Current Amplitude');
    ylabel('Frequency');
    zlabel('MaxEyeMovement LE');
    title(['Date: ', datestr(currentDate)]);
    
    % Add colorbar and make the plot look nice
%     shading interp;
    colormap(jet);
    view(-10,10);
    
    % Set consistent colorbar limits across all plots
    caxis(zLimits);
end
colorbar
% Add a main title
sgtitle('Surface Plots of MaxEyeMovement LE by Date for RP Stim');

%% LA R eye

% Extract unique dates from the table
dates = unique(tab_RP.Date);

% Create figure
figure('Position', [100, 100, 1200, 400]);

% First pass to determine the maximum z-range needed
maxZ = -Inf;
minZ = Inf;

for d = 1:length(dates)
    currentDate = dates(d);
    
    % Filter data for the current date
    dateData = tab_RP(tab_RP.Date == currentDate, :);
    
    % Update global min/max for the z-axis
    if ~isempty(dateData)
        maxZ = max(maxZ, max(dateData.MaxEyeMovement_R));
        minZ = min(minZ, min(dateData.MaxEyeMovement_R));
    end
end

% Define z-axis limits with a small buffer
zBuffer = 0.05 * (maxZ - minZ);
zLimits = [minZ - zBuffer, maxZ + zBuffer];

% Second pass to create the plots with consistent z-axis scaling
for d = 1:length(dates)
    currentDate = dates(d);
    
    % Filter data for the current date
    dateData = tab_RP(tab_RP.Date == currentDate, :);
    
    % Get unique CurrentAmp and Frequency values
    currentAmps = unique(dateData.CurrentAmp);
    frequencies = unique(dateData.Frequency);
    
    % Create a grid for the surface plot
    [X, Y] = meshgrid(currentAmps, frequencies);
    Z = zeros(length(frequencies), length(currentAmps));
    
    % Fill the Z matrix with average MaxEyeMovement_L values
    for i = 1:length(currentAmps)
        for j = 1:length(frequencies)
            % Find rows matching current amplitude and frequency
            matches = dateData(dateData.CurrentAmp == currentAmps(i) & ...
                              dateData.Frequency == frequencies(j), :);
            
            % Calculate mean if there are matches
            if ~isempty(matches)
                Z(j, i) = mean(matches.MaxEyeMovement_R);
            end
        end
    end
    
    % Create subplot for each date
    subplot(1, length(dates), d);
    surf(X, Y, Z);
    
    % Set consistent z-axis limits
    zlim(zLimits);
    
    % Add labels and title
    xlabel('Current Amplitude');
    ylabel('Frequency');
    zlabel('MaxEyeMovement RE');
    title(['Date: ', datestr(currentDate)]);
%     shading interp;
    colormap(jet);
    view(-10,10);


    % Set consistent colorbar limits across all plots
    caxis(zLimits);
end
colorbar;


% Add a main title
sgtitle('Surface Plots of MaxEyeMovement RE by Date for RP Stim');

