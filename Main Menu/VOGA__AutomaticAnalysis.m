%% VOGA__AutomaticAnalysis
%
% This function expects that all relevant VOG files are already Segmented.
% From there, this program applies the heuristic filtering and cycle 
% selection parameters to unanalyzed files and then allows the user to
% either save (accept) or not save (reject) the analysis. It requires
% MATLAB to be run in the directory of interest.
%
function VOGA__AutomaticAnalysis
if ~MakeFolders(cd,0,0)
    disp('Folder structure not present. Generate folders, process raw data, and segment first.')
    return;
end
%Expected folder names
Path = cd; Seg_Path = [Path,filesep,'Segmented Files']; Cyc_Path = [Path,filesep,'Cycle Averages'];
%Find all the unanalyzed segments
progress_tab = assessProgress(Path);
rel_files = progress_tab.Segment(~progress_tab.CycAvgFile&~progress_tab.NotAnalyzeable);
% For each segment, run the Automated Filtering script and allow the user
% input to determine the action
disp('Press "Enter" for save, "m" for discard/do manually, and "l" to leave/quit.')
for j = 1:length(rel_files)
    fname = rel_files{j};
    load([Seg_Path,filesep,fname],'Data');
    try %Don't crash the good files if it errors out
        [CycAvg,analyzed] = AutomatedFiltering(Data,0);
        user_in = input('Keep automated analysis?: ','s');
        while ~isempty(user_in)&&~ismember(lower(user_in),'ml') %make sure the keystrokes are valid (Empty, m or l)
            disp('Invalid key. Only "Enter", "m" and "l" are valid entries.')
            user_in = input('Keep automated analysis?: ','s');
        end
        if strcmp(user_in,'l') %End program
            break;
        elseif isempty(user_in) %Save cycle average
            MakeCycAvg__saveCycAvg(Cyc_Path,strrep(CycAvg.name,'CycAvg_',''),CycAvg,analyzed,1);
        end
    catch %Figure out why it crashed later
        disp(['Error analyzing file: ',fname])
    end
end
end