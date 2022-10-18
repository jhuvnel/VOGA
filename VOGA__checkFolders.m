%% VOGA__checkFolders
% This script checks for the expected folder tree and promts the user to
% creates the folders if they have not yet been made if opt_make_fold is 0
% Returns flag: 0 = folders not there, 1 = folders there
function flag = VOGA__checkFolders(opt_make_fold)
path_folders = extractfield(dir,'name',find(extractfield(dir,'isdir')));
expected = [any(contains(path_folders,'Raw Files')),...
    any(contains(path_folders,'Segmented Files')),...
    any(contains(path_folders,'Cycle Averages')),...
    any(contains(path_folders,'Figures')),...
    any(contains(path_folders,'CRFs'))];
if ~all(expected) %Check for folders
    if opt_make_fold
        make_fold = questdlg('This directory is missing one or more expected folders.','','Make Them','Exit','Make Them');
        if strcmp(make_fold,'Make Them')
            VOGA__makeFolders;
            flag = 1;
        else
            flag = 0;
        end
    else
        flag = 0;
    end
else
    flag = 1;
end
end