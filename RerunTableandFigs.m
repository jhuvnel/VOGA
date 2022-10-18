if ~ispc %AIA Mac
    drive_path = '/Volumes/vnelhuman$/MVI/Study Subjects/';
else %AIA Lab Computer
    drive_path = '\\win.ad.jhu.edu\cloud\vnelhuman$\MVI\Study Subjects\';
end
folders = strrep(strcat({drive_path},table2cell(readtable('Directories.txt','Delimiter',newline,'ReadVariableNames', false))),'/',filesep); %should be in the user path
for i = 1:length(folders)
    disp([num2str(i),'/',num2str(length(folders)),': ',folders{i}])
    RerunTableandFigs(folders{i})
    Path = folders{i};
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
    close all;
end
disp('DONE!')