function VOGA__automaticAnalysis
if ~VOGA__makeFolders(cd,0,0)
    disp('Folder structure not present. Generate folders, process raw data, and segment first.')
    return;
end
%Made for Autoscan data
tic;
Path = cd;
Seg_Path = [Path,filesep,'Segmented Files'];
Cyc_Path = [Path,filesep,'Cycle Averages'];
% Cyc Avg Files
progress_tab = assessProgress(Path);
for j = 1:size(progress_tab,1)    
    fname = progress_tab{j,1}{:};
    load([Seg_Path,filesep,fname],'Data');
    try
        [CycAvg,analyzed] = MakeCycAvg(Data,Cyc_Path,{'Save'},1);   
        MakeCycAvg__saveCycAvg(Cyc_Path,strrep(CycAvg.name,'CycAvg_',''),CycAvg,analyzed,1);
    catch
        disp(fname)
    end
end
%% Make Table and Summary Figure
MakeCycleSummaryTable(cd,[cd,filesep,'Cycle Averages'],0);
%Make Summary Figure
VOGA__makePlots('Parameterized');
%AutoscanParamSelection;
disp(toc)
end
