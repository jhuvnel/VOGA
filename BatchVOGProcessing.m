clear;
clc;
VOGA__makeFolders;
Path = cd;
MakeCycleSummaryTable(Path,[Path,filesep,'Cycle Averages'],1);
code_Path = [userpath,filesep,'VOGA'];
params.Path = cd;
params.Raw_Path = [cd,filesep,'Raw Files'];
params.Seg_Path = [cd,filesep,'Segmented Files'];
params.Cyc_Path = [cd,filesep,'Cycle Averages'];
params.code_Path = code_Path;
data = readtable([code_Path,filesep,'VerInfo.txt'],'ReadVariableNames',false);
params.version = data{1,2}{:};
params.Experimenter = data{2,2}{:}; 
params.annot = 1;
params.YMax = [];
plotGroupCycAvg(params);
disp('Done')