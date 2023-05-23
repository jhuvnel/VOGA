function VOGA__automaticAnalysis
%Made for Autoscan data
tic;
if ~VOGA__makeFolders(cd,0,0)
    disp('Folder structure not present. Generate folders, process raw data, and segment first.')
    return;
end
Path = cd;
Seg_Path = [Path,filesep,'Segmented Files'];
Cyc_Path = [Path,filesep,'Cycle Averages'];
% Cyc Avg Files
progress_tab = assessProgress(Path);
for j = 1:size(progress_tab,1)    
    fname = progress_tab{j,1}{:};
    load([Seg_Path,filesep,fname],'Data');
    [CycAvg,analyzed] = MakeCycAvg(Data,Cyc_Path,{'Save'});   
    MakeCycAvg__saveCycAvg(Cyc_Path,strrep(CycAvg.name,'CycAvg_',''),CycAvg,analyzed);
end
%% Make Table and Summary Figure
MakeCycleSummaryTable(cd,[cd,filesep,'Cycle Averages'],0);
%Make Summary Figure
plot_params = {'0','','1','0'};
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
sub_info = readtable('MVI_Information.xlsx');
params.sub_info = sub_info;
%Add the plot parameters
params.annot = str2double(plot_params{1});
params.YMax = str2double(plot_params{2});
params.plot_eyes = str2double(plot_params{3});
params.plot_fits = str2double(plot_params{4});
params.Path = cd;
params.Raw_Path = [cd,filesep,'Raw Files'];
params.Seg_Path = [cd,filesep,'Segmented Files'];
params.Cyc_Path = [cd,filesep,'Cycle Averages'];
params.which_files = 'All';
plotParamResults(params)
disp(toc)
end
