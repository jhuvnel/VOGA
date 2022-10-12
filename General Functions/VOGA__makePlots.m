function VOGA__makePlots
opts = {'Raw VOG','Segment','Cycle Average','Group Cycle Avg','Parameterized',...
    'Across Subjects','Sphere Plot','Edit Plot Params'};
%Default plot parameters (No annotation, No set Y-axis maximum, Show eye
%movements [not just the sync signal/gyro], No fits on the plot, LRZ axes
%[not xyz])
plot_params = {'0','','1','0','lrz'};
%The subset of options that can be run within a folder that covers one experiment within a visit (not aggregated)
fold_level_opts = {'Raw VOG','Segment','Cycle Average','Group Cycle Avg','Parameterized','Sphere Plot'}; 
%The subset of options that need to be run at the directory above the subject folders(aggregated)
exp_level_opts = {'Across Subjects'};
%Prompt the user to select an option
[ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
    'SelectionMode','single','ListSize',[150 125],'ListString',opts);
while tf % Run until the user selects cancel
    %% Set up for each figure type
    code_Path = [userpath,filesep,'VOGA'];
    params.code_Path = code_Path;
    % Get version and experimenter info from the file
    if ~any(contains(extractfield(dir(userpath),'name'),'VOGA_VerInfo.txt'))
        VOGA__setVersion;
    end
    data = readtable('VOGA_VerInfo.txt','ReadVariableNames',false);
    params.version = data{1,2}{:};
    params.Experimenter = data{2,2}{:}; 
    % Get subject info
    warning('off')
    sub_info = readtable('SubjectInfo.xlsx');
    warning('on')
    params.sub_info = sub_info;
    %Add the plot parameters
    params.annot = str2double(plot_params{1});
    params.YMax = str2double(plot_params{2});
    params.plot_eyes = str2double(plot_params{3});
    params.plot_fits = str2double(plot_params{4});
    params.lrz_xyz = plot_params{5};
    if any(contains(fold_level_opts,opts{ind}))
        new_path = uigetdir('Select the path. Should be a folder containing one experiment type for one visit.');
        if isnumeric(new_path)
            error('No folder selected.')
        end
        cd(new_path)
        if ~VOGA__checkFolders(0)
            error('Expected folder structure not present. Navigate to appropriate folder before trying again.')
        end
        params.Path = cd;
        params.Raw_Path = [cd,filesep,'Raw Files'];
        params.Seg_Path = [cd,filesep,'Segmented Files'];
        params.Cyc_Path = [cd,filesep,'Cycle Averages'];
    end
    if any(contains(exp_level_opts,opts{ind}))
        new_path = uigetdir('Select the path. Should be a folder containing one experiment type for one visit.');
        if isnumeric(new_path)
            error('No folder selected.')
        end
        cd(new_path)
    end
    if strcmp(opts{ind},'Raw VOG')
        %Select file inside this function
        plotRawVOG(params.Raw_Path,params.plot_eyes,params.lrz_xyz) 
    elseif strcmp(opts{ind},'Segment')
        %Select files first
        all_files = extractfield(dir([params.Seg_Path,filesep,'*.mat']),'name');
        [indx,tf2] = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',all_files,'SelectionMode','multiple');
        if tf2
            for i = 1:length(indx)
                load([params.Seg_Path,filesep,all_files{indx(i)}],'Data')
                plotSegment(Data);
            end
        end
    elseif strcmp(opts{ind},'Cycle Average')
        %Select files first
        all_files = extractfield(dir([params.Cyc_Path,filesep,'*Cyc*.mat']),'name');
        [indx,tf2] = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',all_files,'SelectionMode','multiple');
        if tf2
            for i = 1:length(indx)
                load([params.Cyc_Path,filesep,all_files{indx(i)}],'CycAvg')
                plotCycAvg(CycAvg,params.plot_fits,params.lrz_xyz);
            end
        end
    elseif strcmp(opts{ind},'Group Cycle Avg')
        plotGroupCycAvg(params);
    elseif strcmp(opts{ind},'Parameterized')
        plotParamResults(params);
    elseif strcmp(opts{ind},'Across Subjects')
        plotSummaryFigures(params);
    elseif strcmp(opts{ind},'Sphere Plot')
        plotSpherePlot(params);
    elseif strcmp(opts{ind},'Edit Plot Params')
        prompt = {'Descriptive annotation (0/1)','Y-Axis Limit','Show eye (0/1)','Show fits (0/1)','Plot LRZ or XYZ (lrz/xyz)'};
        temp_plot_params = inputdlg(prompt,'Set Plot Defaults',[1 50],plot_params);
        if ~isempty(temp_plot_params)
            plot_params = temp_plot_params;
        end        
    end
    [ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
        'SelectionMode','single',...
        'ListSize',[150 125],...
        'ListString',opts);
end
end