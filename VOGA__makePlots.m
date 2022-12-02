function VOGA__makePlots(plot_type,Path)
opts = {'Raw VOG','Segment','Cycle Average','Group Cycle Avg','Parameterized',...
    'Across Subjects','Sphere Plot','Edit Plot Params'};
%Default plot parameters (No annotation, No set Y-axis maximum, Show eye
%movements [not just the sync signal/gyro], No fits on the plot, LRZ axes
%[not xyz])
plot_params = {'0','','1','0','lrz'};
%The subset of options that can be run within a folder that covers one experiment within a visit (not aggregated)
fold_level_opts = {'Raw VOG','Segment','Cycle Average','Group Cycle Avg','Parameterized','Sphere Plot'};
if nargin < 1 || ~ismember(plot_type,opts)
    plot_type = '';
    %Prompt the user to select an option
    [ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
        'SelectionMode','single','ListSize',[150 125],'ListString',opts);
else
    tf = 1;
    ind = find(ismember(opts,plot_type));
end
if nargin < 2
    Path = cd;
end
while tf % Run until the user selects cancel
    %% Set up for each figure type
    code_Path = [userpath,filesep,'VOGA'];
    params.code_Path = code_Path;
    % Get version and experimenter info from the file
    if ~any(contains(extractfield(dir(userpath),'name'),'VOGA_VerInfo.txt'))
        VOGA__setVersion;
    end
    data = readtable([userpath,filesep,'VOGA_VerInfo.txt'],'ReadVariableNames',false);
    params.version = data{1,2}{:};
    params.Experimenter = data{2,2}{:}; 
    % Get subject info
    sub_info = readtable('SubjectInfo.xlsx');
    params.sub_info = sub_info;
    %Add the plot parameters
    params.annot = str2double(plot_params{1});
    params.YMax = str2double(plot_params{2});
    params.plot_eyes = str2double(plot_params{3});
    params.plot_fits = str2double(plot_params{4});
    params.lrz_xyz = plot_params{5};
    if any(contains(fold_level_opts,opts{ind}))
        if ~VOGA__makeFolders(Path,0) %Expected folders not there
                new_path = uigetdir('Select the path. Should be a folder containing one experiment type for one visit.');
            if isnumeric(new_path)
                error('No folder selected.')
            end
            cd(new_path)
            if ~VOGA__makeFolders(Path,0) %New path still not right
                error('Expected folder structure not present. Navigate to appropriate folder before trying again.')
            end
        end
        params.Path = Path;
        params.Raw_Path = [Path,filesep,'Raw Files'];
        params.Seg_Path = [Path,filesep,'Segmented Files'];
        params.Cyc_Path = [Path,filesep,'Cycle Averages'];
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
    if isempty(plot_type)
        [ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
            'SelectionMode','single','ListSize',[150 125],'ListString',opts);
    else
        tf = 0;
    end
end
end