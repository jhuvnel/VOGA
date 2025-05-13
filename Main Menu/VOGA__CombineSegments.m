%% VOGA__CombineSegments
%
% This function allows a user to combine two segments that contain data
% from the same experiment paragidm. This would happen if they were
% recorded in different files or on different days. This script runs as 
% while loop combining two segments at a time until the "Cancel" button is
% selected. It requires MATLAB to be run in the directory of interest.
%
function VOGA__CombineSegments
if ~MakeFolders(cd,0,0)
    disp('Folder structure not present. Generate folders, process raw data, and segment first.')
    return;
end
%Expected folder names
Path = cd; Seg_Path = [Path,filesep,'Segmented Files'];
tf1 = 1; tf2 = 1;%Needed for the first loop to run
while(tf1&&tf2)
    % Find the segments
    progress_tab = assessProgress(Path);
    all_exp_names = progress_tab.Segment;
    % Have the user select which segments to combine
    [indx1,tf1] = nmlistdlg('PromptString','Select segments to combine. Segment 1 (combined segment name):',...
        'SelectionMode','single','ListSize',[500 500],'ListString',all_exp_names);
    if ~tf1 %Canceled on first inquiry so don't show the second
        break;
    end
    [indx2,tf2] = nmlistdlg('PromptString','Select segments to combine. Segment 2:',...
        'SelectionMode','single','ListSize',[500 500],'ListString',all_exp_names);
    if ~tf2 %Canceled on second inquiry
        break;
    elseif indx1==indx2 %selected same file twice
        uiwait(msgbox('Cannot combine a segment with itself.'))
    else %Combine selected segments and delete the second segement
        load([Seg_Path,filesep,all_exp_names{indx1}],'Data');
        Data1 = Data;
        load([Seg_Path,filesep,all_exp_names{indx2}],'Data');
        Data2 = Data;
        Data = CombineSegments(Data1,Data2);
        delete([Seg_Path,filesep,all_exp_names{indx2}])
        save([Seg_Path,filesep,all_exp_names{indx1}],'Data')
    end
end
end