function flag = VOGA__checkFolders
    path_folders = extractfield(dir,'name',find(extractfield(dir,'isdir')));
    expected = [any(contains(path_folders,'Raw Files')),any(contains(path_folders,'Segmented Files')),any(contains(path_folders,'Cycle Averages'))];
    if ~all(expected) %Check for folders
        flag = 1;
        uiwait(msgbox('This directory is missing one or more expected folders. Run "Make Folders" to make the missing folders.'))
    else
        flag = 0;
    end
end