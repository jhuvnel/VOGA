%% Collapse Folder
% Takes the files in folder "folder" in the path "Path" and puts them
% directly in the Path with the name folder_file
function collapseFolder(Path,folder)
if ~isfolder([Path,filesep,folder]) %trust but verify
    disp(['The following directory does not exist: ',Path,filesep,folder])
    return;
end
fold_items = dir([Path,filesep,folder]);
fold_items(ismember(extractfield(fold_items,'name'),{'.','..'})) = [];
%Collapse any sub-folders
sub_d = extractfield(fold_items,'name',extractfield(fold_items,'isdir'));
for i = 1:length(sub_d)
    collapseFolder([Path,filesep,folder],sub_d{i}) %Nifty recursion
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