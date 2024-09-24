function VOGA__automaticAnalysis
if ~VOGA__makeFolders(cd,0,0)
    disp('Folder structure not present. Generate folders, process raw data, and segment first.')
    return;
end
Path = cd;
Seg_Path = [Path,filesep,'Segmented Files'];
Cyc_Path = [Path,filesep,'Cycle Averages'];
%Only run when the Cycle Average file doesn't yet exist
progress_tab = assessProgress(Path);
rel_files = progress_tab.Segment(~progress_tab.CycAvgFile&~progress_tab.NotAnalyzeable);
disp('Press "Enter" for save, "m" for discard/do manually, and "l" to leave/quit.')
for j = 1:length(rel_files)
    fname = rel_files{j};
    load([Seg_Path,filesep,fname],'Data');
    try
        [CycAvg,analyzed] = AutomatedFiltering(Data,0);
        user_in = input('Keep automated analysis?: ','s');
        while ~isempty(user_in)&&~ismember(lower(user_in),'ml')
            disp('Invalid key. Only "Enter", "m" and "l" are valid entries.')
            user_in = input('Keep automated analysis?: ','s');
        end
        if strcmp(user_in,'l')
            break;
        elseif isempty(user_in)
            MakeCycAvg__saveCycAvg(Cyc_Path,strrep(CycAvg.name,'CycAvg_',''),CycAvg,analyzed,1);
        end
    catch
        disp(['Error analyzing file: ',fname])
    end
end
end