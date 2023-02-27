%% Collapse Folder
% Takes the files in folder "folder" in the path "Path" and puts them
% directly in the Path with the name folder_file
% Third input arg decides whether to move subfolders or not
function collapseFolder(Path,folder,collapse_sub)
if ~isfolder([Path,filesep,folder]) %trust but verify
    disp(['The following directory does not exist: ',Path,filesep,folder])
    return;
end
if collapse_sub %Collapse any sub-folders
    fold_items = dir([Path,filesep,folder]);
    fold_items(ismember(extractfield(fold_items,'name'),{'.','..'})) = [];
    sub_d = extractfield(fold_items,'name',extractfield(fold_items,'isdir'));
    for i = 1:length(sub_d)
        collapseFolder([Path,filesep,folder],sub_d{i},collapse_sub) %Nifty recursion
    end    
end
%Now rename all files
fold_items = dir([Path,filesep,folder]);
fold_items(ismember(extractfield(fold_items,'name'),{'.','..'})) = [];
fnames = extractfield(fold_items,'name');
if isempty(fnames)
    disp(['The following directory was empty: ',Path,filesep,folder])
    return;
end
for i = 1:length(fnames)
    movefile([Path,filesep,folder,filesep,fnames{i}],[Path,filesep,folder,'_',fnames{i}])
end
rmdir([Path,filesep,folder])
end