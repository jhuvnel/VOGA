%% MakeCycleSummaryTable
% This version takes in all CycAvg types that have been defined and saves a
% table for each
% It can take in arguments of the path to save the table and the path where
% all the cycle averages are located. If not provided, it will allow the
% user to provide a text file with the directories and assume the tables
% should be saved in the current directory.

function [all_results,cyc_params] = MakeCycleSummaryTable(out_path,Cyc_Path,rerun,display_progress)
% Handle inputs
if nargin < 4 | ~islogical(display_progress)
    display_progress = 1;
end
if nargin < 3
    rerun = 1;
else
    rerun = ~strcmpi(rerun,'load');
end
if nargin < 1
    out_path = cd;
end
if nargin < 2 % Get the paths from a .txt file in a directory above the listed paths (avoids using differing drive names)
    [path2,path1] = uigetfile('*.txt','Select the text file with all the Cycle Average paths.');
    Cyc_Path = strcat(path1,strrep(strrep(importdata([path1,filesep,path2]),'/',filesep),'\',filesep));
end
all_results = [];
cyc_params = [];
fstruct = dir([Cyc_Path,filesep,'*.mat']);
fstruct(contains(extractfield(fstruct,'name'),'NotAnalyzeable')) = []; %remove un-analyzeable files
if isempty(fstruct)
    disp(['No Cycle Average Files found in the directories: ',Cyc_Path])    
    return;
end
files = join([extractfield(fstruct,'folder'),extractfield(fstruct,'name')],filesep);
if display_progress
    disp([num2str(length(files)),' files found.'])
end
%% Now analyze each file and get the table
tabs = cell(length(files),1);
cyc_params = cell(length(files),2);
for i = 1:length(files)
    if rerun&&display_progress
        disp([num2str(i),'/',num2str(length(files)),': ',files{i}])
    end
    a = load(files{i});
    b = fieldnames(a);
    CycAvg = a.(b{1});
    slash = strfind(files{i},filesep);
    fname = files{i}(slash(end)+1:end);
    if ~isfield(CycAvg,'name')||~strcmp(CycAvg.name,fname)||~isfield(CycAvg,'parameterized')||rerun
        CycAvg.name = fname;
        CycAvg = ParameterizeCycAvg(CycAvg);
        save(files{i},'CycAvg')
    end
    tabs{i} = CycAvg.parameterized;
    cyc_params(i,1) = files(i);
    cyc_params{i,2} = CycAvg.cycle_params;
end
all_results = vertcat(tabs{:});
%% Save
res_files = [dir([out_path,filesep,'*Results.mat']);dir([out_path,filesep,'*CycParam.mat'])];
for i = 1:length(res_files)
    delete([res_files(i).folder,filesep,res_files(i).name])
end
out_fname = [out_path,filesep,char(datetime("now"),'yyyyMMdd_HHmmss'),'_VOG'];
save([out_fname,'Results.mat'],'all_results')
save([out_fname,'CycParam.mat'],'cyc_params')
end