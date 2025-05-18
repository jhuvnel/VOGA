
function start_end_idx = select_peak_roi(cyc_eyemove_r_neg,cyc_eyemove_r_pos,cyc_eyemove_l_neg,cyc_eyemove_l_pos)

fig = figure(10);
set(fig,'units','normalized','outerposition',[.05 .05 .9 .9]);

% Create the subplots and initialize variables
subplot(2,2,1)
plot(cyc_eyemove_r_neg');
title('RE Neg')
subplot(2,2,3)
plot(cyc_eyemove_r_pos')
title('RE Pos')
subplot(2,2,2)
plot(cyc_eyemove_l_neg');
title('LE Neg')
subplot(2,2,4)
plot(cyc_eyemove_l_pos')
title('LE Pos')

% Initialize variables to store start and finish points (as indices)
r_neg_x_start = NaN;
r_neg_x_finish = NaN;
l_neg_x_start = NaN;
l_neg_x_finish = NaN;
r_pos_x_start = NaN;
r_pos_x_finish = NaN;
l_pos_x_start = NaN;
l_pos_x_finish = NaN;

% Enable interactivity for all subplots and fix Y limits
for i = 1:4
    subplot(2,2,i);
    enableDefaultInteractivity(gca);
    YL = get(gca, 'YLim');
    ylim manual;
    ylim(YL);
end

% Main selection loop - allows restarting with keypress
restart = true;
while restart
    % Clear previous reference lines if they exist
    for i = 1:4
        subplot(2,2,i);
        cla reset;
        
        % Redraw the original plots
        if i == 1
            plot(cyc_eyemove_r_neg');
            title('RE Neg');
            data_length = size(cyc_eyemove_r_neg, 2);
        elseif i == 2
            plot(cyc_eyemove_l_neg');
            title('LE Neg');
            data_length = size(cyc_eyemove_l_neg, 2);
        elseif i == 3
            plot(cyc_eyemove_r_pos');
            title('RE Pos');
            data_length = size(cyc_eyemove_r_pos, 2);
        else
            plot(cyc_eyemove_l_pos');
            title('LE Pos');
            data_length = size(cyc_eyemove_l_pos, 2);
        end
        
        % Re-enable interactivity
        enableDefaultInteractivity(gca);
        YL = get(gca, 'YLim');
        ylim manual;
        ylim(YL);
    end
    
%     disp('---------------');
%     disp('Starting new selection session');
%     disp('Press any key during selection to restart');
    
    try
        % RE Neg (top-left, subplot 1) - start point
%         disp('1a. Select START point for RE Neg (top-left)');
        subplot(2,2,1);
        [x_click, ~, button] = ginput(1);
        
        % Check if key was pressed (non-mouse button)
        if button > 5
            disp('Key pressed. Restarting selection...');
            continue;
        end
        
        % Convert click x-coordinate to nearest index
        data_length = size(cyc_eyemove_r_neg, 2);
        r_neg_x_start = round(x_click);
        r_neg_x_start = max(1, min(data_length, r_neg_x_start)); % Ensure index is within bounds
        
        hold on;
        plot([r_neg_x_start r_neg_x_start], get(gca, 'YLim'), 'g-', 'LineWidth', 2);
%         disp(['Stored r_neg_x_start = ', num2str(r_neg_x_start), ' (index)']);
        
        % RE Neg (top-left, subplot 1) - finish point
%         disp('1b. Select FINISH point for RE Neg (top-left)');
        [x_click, ~, button] = ginput(1);
        
        % Check if key was pressed
        if button > 5
            disp('Key pressed. Restarting selection...');
            continue;
        end
        
        % Convert click x-coordinate to nearest index
        r_neg_x_finish = round(x_click);
        r_neg_x_finish = max(1, min(data_length, r_neg_x_finish)); % Ensure index is within bounds
        
        plot([r_neg_x_finish r_neg_x_finish], get(gca, 'YLim'), 'r-', 'LineWidth', 2);
        hold off;
%         disp(['Stored r_neg_x_finish = ', num2str(r_neg_x_finish), ' (index)']);
        
        % LE Neg (top-right, subplot 2) - start point
%         disp('2a. Select START point for LE Neg (top-right)');
        subplot(2,2,2);
        [x_click, ~, button] = ginput(1);
        
        % Check if key was pressed
        if button > 5
            disp('Key pressed. Restarting selection...');
            continue;
        end
        
        % Convert click x-coordinate to nearest index
        data_length = size(cyc_eyemove_l_neg, 2);
        l_neg_x_start = round(x_click);
        l_neg_x_start = max(1, min(data_length, l_neg_x_start)); % Ensure index is within bounds
        
        hold on;
        plot([l_neg_x_start l_neg_x_start], get(gca, 'YLim'), 'g-', 'LineWidth', 2);
%         disp(['Stored l_neg_x_start = ', num2str(l_neg_x_start), ' (index)']);
        
        % LE Neg (top-right, subplot 2) - finish point
%         disp('2b. Select FINISH point for LE Neg (top-right)');
        [x_click, ~, button] = ginput(1);
        
        % Check if key was pressed
        if button > 5
            disp('Key pressed. Restarting selection...');
            continue;
        end
        
        % Convert click x-coordinate to nearest index
        l_neg_x_finish = round(x_click);
        l_neg_x_finish = max(1, min(data_length, l_neg_x_finish)); % Ensure index is within bounds
        
        plot([l_neg_x_finish l_neg_x_finish], get(gca, 'YLim'), 'r-', 'LineWidth', 2);
        hold off;
        disp(['Stored l_neg_x_finish = ', num2str(l_neg_x_finish), ' (index)']);
        
        % RE Pos (bottom-left, subplot 3) - start point
%         disp('3a. Select START point for RE Pos (bottom-left)');
        subplot(2,2,3);
        [x_click, ~, button] = ginput(1);
        
        % Check if key was pressed
        if button > 5
            disp('Key pressed. Restarting selection...');
            continue;
        end
        
        % Convert click x-coordinate to nearest index
        data_length = size(cyc_eyemove_r_pos, 2);
        r_pos_x_start = round(x_click);
        r_pos_x_start = max(1, min(data_length, r_pos_x_start)); % Ensure index is within bounds
        
        hold on;
        plot([r_pos_x_start r_pos_x_start], get(gca, 'YLim'), 'g-', 'LineWidth', 2);
%         disp(['Stored r_pos_x_start = ', num2str(r_pos_x_start), ' (index)']);
        
        % RE Pos (bottom-left, subplot 3) - finish point
        disp('3b. Select FINISH point for RE Pos (bottom-left)');
        [x_click, ~, button] = ginput(1);
        
        % Check if key was pressed
        if button > 5
            disp('Key pressed. Restarting selection...');
            continue;
        end
        
        % Convert click x-coordinate to nearest index
        r_pos_x_finish = round(x_click);
        r_pos_x_finish = max(1, min(data_length, r_pos_x_finish)); % Ensure index is within bounds
        
        plot([r_pos_x_finish r_pos_x_finish], get(gca, 'YLim'), 'r-', 'LineWidth', 2);
        hold off;
%         disp(['Stored r_pos_x_finish = ', num2str(r_pos_x_finish), ' (index)']);
        
        % LE Pos (bottom-right, subplot 4) - start point
%         disp('4a. Select START point for LE Pos (bottom-right)');
        subplot(2,2,4);
        [x_click, ~, button] = ginput(1);
        
        % Check if key was pressed
        if button > 5
            disp('Key pressed. Restarting selection...');
            continue;
        end
        
        % Convert click x-coordinate to nearest index
        data_length = size(cyc_eyemove_l_pos, 2);
        l_pos_x_start = round(x_click);
        l_pos_x_start = max(1, min(data_length, l_pos_x_start)); % Ensure index is within bounds
        
        hold on;
        plot([l_pos_x_start l_pos_x_start], get(gca, 'YLim'), 'g-', 'LineWidth', 2);
%         disp(['Stored l_pos_x_start = ', num2str(l_pos_x_start), ' (index)']);
        
        % LE Pos (bottom-right, subplot 4) - finish point
%         disp('4b. Select FINISH point for LE Pos (bottom-right)');
        [x_click, ~, button] = ginput(1);
        
        % Check if key was pressed
        if button > 5
            disp('Key pressed. Restarting selection...');
            continue;
        end
        
        % Convert click x-coordinate to nearest index
        l_pos_x_finish = round(x_click);
        l_pos_x_finish = max(1, min(data_length, l_pos_x_finish)); % Ensure index is within bounds
        
        plot([l_pos_x_finish l_pos_x_finish], get(gca, 'YLim'), 'r-', 'LineWidth', 2);
        hold off;
%         disp(['Stored l_pos_x_finish = ', num2str(l_pos_x_finish), ' (index)']);
        
        % All selections complete, exit the loop
        restart = false;
        
    catch err
        disp('Error occurred during selection. Restarting...');
        disp(err.message);
    end
end
close(figure(10))

% Return all indices in the output matrix

start_end_idx = [r_neg_x_start r_neg_x_finish;
    l_neg_x_start l_neg_x_finish;
    r_pos_x_start r_pos_x_finish;
    l_pos_x_start l_pos_x_finish];