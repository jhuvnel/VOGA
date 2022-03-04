function RerunTableandFigs(Path)
params.Path = Path;
params.Raw_Path = [Path,filesep,'Raw Files'];
params.Seg_Path = [Path,filesep,'Segmented Files'];
params.Cyc_Path = [Path,filesep,'Cycle Averages'];
params.code_Path = [userpath,filesep,'VOGA'];
data = readtable([userpath,filesep,'VOGA_VerInfo.txt'],'ReadVariableNames',false);
params.version = data{1,2}{:};
params.Experimenter = data{2,2}{:}; 
params.annot = 1;
params.YMax = [];
VOGA__makeFolders(Path)
MakeCycleSummaryTable(Path,[Path,filesep,'Cycle Averages'],1);
plotGroupCycAvg(params);
end