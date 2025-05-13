%% VOGA__Process Raw Data
%
% This function expects that raw VOG files are in the current directory or
% specified Path. It splits larger files into smaller ones for GNO and ESC,
% gives NL descriptive names and moves all files to the Raw Files folder.
%
function VOGA__ProcessRawData(Path)
if nargin < 1
    Path = cd;
end
if ~MakeFolders(Path)
    return;
end
Raw_Path = [Path,filesep,'Raw Files'];
%Folders expected, or that should be ignored if present, edit as needed
ignored_folders = {'Raw Files','Segmented Files','Cycle Averages','Figures',...
    'Archive','trash','Summary','.','CRF'}; 
%Files that should be ignored if present, edit as needed
ignored_files = {'thumbs.db','DS_Store','Results.mat','Param.mat',...
    'Parameters.mat','Progress.mat','.ppt'};
%% Handle folders in the main path
% Some folders are expected from dumping data from the NL, ESC or GNO
% systems and are processed below. Some are the normal file tree folders
% and some are expected if there has already been some analysis.
fold_Path = extractfield(dir(Path),'name',extractfield(dir(Path),'isdir'));
fold_Path(contains(lower(fold_Path),lower(ignored_folders))) = [];
for i = 1:length(fold_Path) %Deal with each of the rest of the existing folders
    if contains(lower(fold_Path{i}),'video') %GNO videos, move the videos in Raw Files
        movefile([Path,filesep,fold_Path{i},filesep,'*'],Raw_Path)
        rmdir([Path,filesep,fold_Path{i}])
    elseif ~isempty(dir([Path,filesep,fold_Path{i},filesep,'*.dat']))||~isempty(dir([Path,filesep,fold_Path{i},filesep,'*Test.csv'])) %Detect NKI/NL folder by the presence of .dat or .csv files
        fnames = extractfield(dir([Path,filesep,fold_Path{i}]),'name',~ismember(extractfield(dir([Path,filesep,fold_Path{i}]),'name'),{'.','..'}));
        for ii = 1:length(fnames) %Add the folder name to the beginning of each file name
            movefile([Path,filesep,fold_Path{i},filesep,fnames{ii}],[Raw_Path,filesep,fold_Path{i},'_',fnames{ii}])
        end  
        rmdir([Path,filesep,fold_Path{i}])
    elseif contains(lower(fold_Path{i}),'calibration') %ESC and sometimes LDVOG, move the whole folder to Raw Files
        movefile([Path,filesep,fold_Path{i}],[Raw_Path,filesep,fold_Path{i}])
    elseif contains(Path,'ESC') %ESC3 has a lot of folders, use the recursive loop in that function to deal with it
        collapseFolder(Path,fold_Path{i},1)
    else
        disp(['Ignored folder: ',fold_Path{i}])
    end
end
%% Handle files in the main path 
% Some files are expected if there has already been some analysis. These
% should be ignored.
file_Path = extractfield(dir(Path),'name',~extractfield(dir(Path),'isdir')&...
    ~contains(extractfield(dir(Path),'name'),ignored_files));
if ~isempty(file_Path)
    if contains(Path,'ESC') %detect ESC files by this convention
        % For the 2016-2022 files, the data files don't yet have the .mat extension needed to open them
        for j = 1:length(file_Path)
            rel_file = file_Path{j};
            if ~all(ismember(lower(rel_file(find(rel_file=='.',1,'last')+1:end)),'abcdefghijklmnopqrstuvwxyz')) %not valid file extenstion
                movefile([Path,filesep,rel_file],[Raw_Path,filesep,rel_file,'.mat'])
            else
                movefile([Path,filesep,rel_file],[Raw_Path,filesep,rel_file])
            end
        end
    elseif contains(Path,'GNO') %GNO folder
        rel_file = file_Path(cellfun(@length,strfind(file_Path,'_'))==7&contains(file_Path,{'.txt','.xml','.csv'})); 
        disp('Splitting GNO files.')
        if ~isfolder([Raw_Path,filesep,'Combined Files']) %Make this folder for the processed files to go into
            mkdir([Raw_Path,filesep,'Combined Files'])
        end 
        for j = 1:length(rel_file)
            disp(['File ',num2str(j),'/',num2str(length(rel_file))])
            splitGNOfile(Path,rel_file{j},0) %Already deidentified above
            movefile([Path,filesep,rel_file{j}],[Raw_Path,filesep,'Combined Files',filesep,'Processed_',rel_file{j}]) %Add this prefix to not redo files that have already been processed
        end
        disp('Done splitting GNO files.')
    end
    file_Path = extractfield(dir(Path),'name',~extractfield(dir(Path),'isdir')&...
        ~contains(extractfield(dir(Path),'name'),ignored_files));
    %Move all non .pdf files. All ESC .pdfs have already been moved
    rel_file = file_Path(~contains(file_Path,'.pdf'));
    for j = 1:length(rel_file)
        movefile([Path,filesep,rel_file{j}],[Raw_Path,filesep,rel_file{j}])
    end    
end
%% Create Notes Files
MakeNotes(Raw_Path);
close all;
end