function VOGA__makePlots(plot_type,Path)
%Plot Types
opts = {'Raw VOG','Segment','Cycle Average','Parameterized',...
    'Across Subjects','Sphere Plot'};
% Take in inputs
if nargin < 2
    Path = cd;
end
if nargin < 1 || ~ismember(plot_type,opts)
    %Prompt the user to select an option
    ind = nmlistdlg('PromptString','Select an plot to make:',...
        'SelectionMode','single','ListSize',[150 125],'ListString',opts);
    if isempty(ind)
        return;
    end
else
    ind = find(ismember(opts,plot_type));
end
%Error handle expected foldered tree not there. Not needed for summary
%figures made in "Across Subjects"
if ~contains(opts{ind},'Across Subjects')&&~VOGA__makeFolders(Path,0) 
    error('Expected folder structure not present. Navigate to appropriate folder before trying again.')
end
%% Set up for each figure type
params.code_Path = [userpath,filesep,'VOGA'];
% Get version and experimenter info from the file
VOGA_VerInfo = rows2vars(readtable([userpath,filesep,'VOGA_VerInfo.txt'],...
    'ReadVariableNames',false,'ReadRowNames',true));
params.version = VOGA_VerInfo.Version{:};
params.Experimenter = VOGA_VerInfo.Experimenter{:};
params.MVIPath = VOGA_VerInfo.Path{:};
%Check that the MVI Path is valid
if ~isfolder(params.MVIPath)
    error('Invalid MVI "Study Subject" folder path.')
end
%Get subject info
params.sub_info = readtable([params.code_Path,filesep,'MVI_Information.xlsx']);
%Add the plot parameters
params.Path = Path;
params.Raw_Path = [Path,filesep,'Raw Files'];
params.Seg_Path = [Path,filesep,'Segmented Files'];
params.Cyc_Path = [Path,filesep,'Cycle Averages'];
if strcmp(opts{ind},'Raw VOG')
    %Select file inside this function
    plotRawVOG(params.Raw_Path) 
elseif strcmp(opts{ind},'Segment')
    %Select files first
    all_files = extractfield(dir([params.Seg_Path,filesep,'*.mat']),'name');
    indx = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',all_files,'SelectionMode','multiple');
    for i = 1:length(indx)
        load([params.Seg_Path,filesep,all_files{indx(i)}],'Data')
        plotSegment(Data);
    end
elseif strcmp(opts{ind},'Cycle Average')
    %Select files first
    all_files = extractfield(dir([params.Cyc_Path,filesep,'*Cyc*.mat']),'name');
    indx = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',all_files,'SelectionMode','multiple');
    for i = 1:length(indx)
        load([params.Cyc_Path,filesep,all_files{indx(i)}],'CycAvg')
        plotCycAvg(CycAvg,params.plot_fits);
    end
elseif strcmp(opts{ind},'Parameterized')
    plotParamResults(params);
elseif strcmp(opts{ind},'Across Subjects')
    plotSummaryFigures(params);
elseif strcmp(opts{ind},'Sphere Plot')
    plotSpherePlot(params);       
end
end