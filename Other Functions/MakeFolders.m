%% MakeFolders.m
% This script checks for the expected VOG folder tree and promts the user
% to create the folders if they have not yet been made.
% The expected folders are defined by the variable "folders" and can be
% updated in future versions if needed.
%
% % Inputs: 
% Path: the target directory, which defaults to the current directory
% opt_make_fold: 0 = no directories are created but existing ones may be renamed
%                1 (default) = directories may be created and/or renamed
% promt_user: 0 = directories will be made automatically
%             1 (default) = the user will be asked in a pop-up dialogue 
% whether to make the directories. The benefit of the dialogue is that the 
% user can cancel the process if they realize they are not in a directory 
% that needs these folders.
%
% % Outputs:
% flag: 0 = folders do not exist
%       1 = folders exist
%
function flag = MakeFolders(Path,opt_make_fold,prompt_user)
%Expected folders (change as needed)
folders = {'Raw Files','Segmented Files','Cycle Averages','Figures'};
% Input handling, set defaults
if nargin < 1||isempty(Path)
    Path = cd;
end
if nargin < 2
    opt_make_fold = 1;
end
if nargin < 3
    prompt_user = 1;
end
%Rename old folder names from the LDVOG era if present
if isfolder([Path,filesep,'Raw LD VOG Files'])
    movefile([Path,filesep,'Raw LD VOG Files'],[Path,filesep,'Raw Files'])
end
if isfolder([Path,filesep,'CycAvg'])
    movefile([Path,filesep,'CycAvg'],[Path,filesep,'Cycle Averages'])
end
if isfolder([Path,filesep,'Cyc_Avg'])
    movefile([Path,filesep,'Cyc_Avg'],[Path,filesep,'Cycle Averages'])
end
% Find the names of folders on the path
path_folders = extractfield(dir(Path),'name',find(extractfield(dir(Path),'isdir')));
if ~all(ismember(folders,path_folders))&&opt_make_fold %Check for folders and for the option to make folders
    make_fold = 'Make Them';
    if prompt_user %Have the user prompted about whether or not to make the folders
        make_fold = questdlg('This directory is missing one or more expected folders.','','Make Them','Exit','Make Them');
    end
    if strcmp(make_fold,'Make Them')
        missing_fold = folders(~ismember(folders,path_folders));
        for i = 1:length(missing_fold)
            mkdir([Path,filesep,missing_fold{i}])
        end
        path_folders = extractfield(dir,'name',find(extractfield(dir,'isdir'))); % Update the names of folders on the path
    end
end
flag = all(ismember(folders,path_folders)); %Check for folders
end