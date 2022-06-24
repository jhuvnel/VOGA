function flag = VOGA__checkFolders
    path_folders = extractfield(dir,'name',find(extractfield(dir,'isdir')));
    expected = [any(contains(path_folders,'Raw Files')),...
        any(contains(path_folders,'Segmented Files')),...
        any(contains(path_folders,'Cycle Averages')),...
        any(contains(path_folders,'Figures')),...
        any(contains(path_folders,'CRFs'))];
    if ~all(expected) %Check for folders
        make_fold = questdlg('This directory is missing one or more expected folders. Click "OK" to make the missing folders or "Run Anyway".','','OK','Run Anyway','OK');
        if strcmp(make_fold,'OK')
            VOGA__makeFolders
            flag = 0;
        elseif strcmp(make_fold,'Run Anyway')
            flag = 0;
        else
            flag = 1;
        end
    else
        flag = 0;
    end
end