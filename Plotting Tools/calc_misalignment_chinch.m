function results = calc_misalignment_chinch(data,plotfromday0)
    % Get unique dates and current amplitudes
    dates = unique(data.Date);
    currents = unique(data.CurrentAmp);
    
    % Ensure the reference date (day 1) is the first date
    try
        reference_date = datetime('06-Feb-2025');
        if ~any(dates == reference_date)
            error('Reference date 06-Feb-2025 not found in data');
        end
    catch
        reference_date = datetime('22-Nov-2024');
        if ~any(dates == reference_date)
            error('Reference date 22-Nov-2024 not found in data');
        end
    end
    
    % Initialize results structure
    results = struct();
    results.AngularDeviation_L = cell(length(currents), length(dates));
    results.AngularDeviation_R = cell(length(currents), length(dates));
    results.MeanAngularDeviation_L = zeros(length(currents), length(dates));
    results.MeanAngularDeviation_R = zeros(length(currents), length(dates));
    results.StdAngularDeviation_L = zeros(length(currents), length(dates));
    results.StdAngularDeviation_R = zeros(length(currents), length(dates));
    
    % Initialize tables to store results in a more readable format
    L_misalignment_table = table();
    R_misalignment_table = table();
    
    % Track row index for tables
    row_idx = 1;
    
    % For each current amplitude
    for i = 1:length(currents)
        current = currents(i);
        
        % Get reference (day 1) data for this current
        ref_indices = data.Date == reference_date & data.CurrentAmp == current;
        
        % Skip if no reference data for this current
        if ~any(ref_indices)
            continue;
        end
        
        % Extract the axes for the reference date 
        if iscell(data.MaxAxis_L)
            ref_axes_L = cell2mat(data.MaxAxis_L(ref_indices));
            ref_axes_R = cell2mat(data.MaxAxis_R(ref_indices));
        else
            ref_axes_L = data.MaxAxis_L(ref_indices, :);
            ref_axes_R = data.MaxAxis_R(ref_indices, :);
        end
        
        if plotfromday0
            % Calculate mean reference axes
            mean_ref_axis_L = mean(ref_axes_L, 1);
            mean_ref_axis_R = mean(ref_axes_R, 1);
        else
            % NOTE: this is assuming all the data we have has the same stim
            % axis
            mean_ref_axis_L = -data.StimAxis{1};
            mean_ref_axis_R = -data.StimAxis{1};
        end
        
        % For each comparison date (excluding reference date)
        for j = 1:length(dates)
            comparison_date = dates(j);
          
            
            % Get data for this date and current
            comp_indices = data.Date == comparison_date & data.CurrentAmp == current;
            
            % Skip if no data for this date/current combination
            if ~any(comp_indices)
                continue;
            end
            
            % Extract the axes for the comparison date
            if iscell(data.MaxAxis_L)
                comp_axes_L = cell2mat(data.MaxAxis_L(comp_indices));
                comp_axes_R = cell2mat(data.MaxAxis_R(comp_indices));
            else
                comp_axes_L = data.MaxAxis_L(comp_indices, :);
                comp_axes_R = data.MaxAxis_R(comp_indices, :);
            end
            % Calculate misalignment
            [L_mean, L_std, L_all] = calc_misalignment(mean_ref_axis_L, comp_axes_L);
            [R_mean, R_std, R_all] = calc_misalignment(mean_ref_axis_R, comp_axes_R);
            R_all = rmoutliers(R_all,'percentiles',[10 80]);
            R_mean = mean(R_all);
            R_std = std(R_all);
            L_all = rmoutliers(L_all,'percentiles',[10 80]);
            L_mean = mean(L_all);
            L_std = std(L_all);
            
            
            % Figure out the date index (since we're skipping ref date)
            date_idx = find(dates == comparison_date);% - sum(dates(1:find(dates == comparison_date)) == reference_date);
            
            % Store results in the structure
            results.AngularDeviation_L{i, date_idx} = L_all;
            results.AngularDeviation_R{i, date_idx} = R_all;
            results.MeanAngularDeviation_L(i, date_idx) = L_mean;
            results.MeanAngularDeviation_R(i, date_idx) = R_mean;
            results.StdAngularDeviation_L(i, date_idx) = L_std;
            results.StdAngularDeviation_R(i, date_idx) = R_std;
            
            % Add to tables for easier viewing
            for cycle_idx = 1:length(L_all)
                L_misalignment_table(row_idx, :) = {comparison_date, current, L_all(cycle_idx)};
                R_misalignment_table(row_idx, :) = {comparison_date, current, R_all(cycle_idx)};
                row_idx = row_idx + 1;
            end
        end
    end
    
    % Set table variable names
    L_misalignment_table.Properties.VariableNames = {'Date', 'CurrentAmp', 'Misalignment'};
    R_misalignment_table.Properties.VariableNames = {'Date', 'CurrentAmp', 'Misalignment'};
    
    % Add tables to results
    results.L_misalignment_table = L_misalignment_table;
    results.R_misalignment_table = R_misalignment_table;
    
    % Add current amplitudes and dates to results for reference
    results.Currents = currents;
    results.ComparisonDates = setdiff(dates, reference_date);
    results.ReferenceDate = reference_date;
end