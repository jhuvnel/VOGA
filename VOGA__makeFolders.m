%% VOGA__makeFolders
% This script checks for the expected folder tree and promts the user to
% creates the folders if they have not yet been made.
% If opt_make_fold is 0, no directories are created/deleted. 
% If promt_user is 0, the directories will be make automatically.
%
% Returns flag: 0 = folders not there, 1 = folders there
function flag = VOGA__makeFolders(Path,opt_make_fold,prompt_user)
if nargin < 1||isempty(Path)
    Path = cd;
end
if nargin < 2
    opt_make_fold = 1;
end
if nargin < 3
    prompt_user = 1;
end
%Rename old folder names if present
if isfolder([Path,filesep,'Raw LD VOG Files'])
    movefile([Path,filesep,'Raw LD VOG Files'],[Path,filesep,'Raw Files'])
end
if isfolder([Path,filesep,'CycAvg'])
    movefile([Path,filesep,'CycAvg'],[Path,filesep,'Cycle Averages'])
end
if isfolder([Path,filesep,'Cyc_Avg'])
    movefile([Path,filesep,'Cyc_Avg'],[Path,filesep,'Cycle Averages'])
end
%Now look for missing folders and prompt the user to create them
path_folders = extractfield(dir(Path),'name',find(extractfield(dir(Path),'isdir')));
folders = {'Raw Files','Segmented Files','Cycle Averages','Figures'};
if ~all(ismember(folders,path_folders))&&opt_make_fold %Check for folders
    if prompt_user
        make_fold = questdlg('This directory is missing one or more expected folders.','','Make Them','Exit','Make Them');
    else
        make_fold = 'Make Them';
    end
    if strcmp(make_fold,'Make Them')
        missing_fold = folders(~ismember(folders,path_folders));
        for i = 1:length(missing_fold)
            mkdir([Path,filesep,missing_fold{i}])
        end
        path_folders = extractfield(dir,'name',find(extractfield(dir,'isdir')));
    end
end
flag = all(ismember(folders,path_folders));
end