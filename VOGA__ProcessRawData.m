function flag = VOGA__ProcessRawData(Path)
if nargin < 1
    Path = cd;
end
if ~VOGA__makeFolders(Path)
    return;
end
Raw_Path = [Path,filesep,'Raw Files'];
%% Handle folders in the main path
fold_Path = extractfield(dir(Path),'name',extractfield(dir(Path),'isdir'));
fold_Path(contains(lower(fold_Path),{'archive','trash','.'})) = [];
fold_Path(contains(fold_Path,{'Raw Files','Segmented Files','Cycle Averages','Figures','CRFs'})) = []; %Expected folders
for i = 1:length(fold_Path)
    if contains(lower(fold_Path{i}),'video') %GNO videos, move the videos in Raw Files
        movefile([Path,filesep,fold_Path{i},filesep,'*'],Raw_Path)
        rmdir([Path,filesep,fold_Path{i}])
    elseif contains(lower(fold_Path{i}),'calibration') %ESC and sometimes LDVOG, move the whole folder
        movefile([Path,filesep,fold_Path{i}],[Raw_Path,filesep,fold_Path{i}])
    elseif ~isempty(dir([Path,filesep,fold_Path{i},filesep,'*.dat']))||~isempty(dir([Path,filesep,fold_Path{i},filesep,'*Test.csv'])) %Detect NKI/NL folder by the presence of .dat files
        fnames = extractfield(dir([Path,filesep,fold_Path{i}]),'name',~ismember(extractfield(dir([Path,filesep,fold_Path{i}]),'name'),{'.','..'}));
        for ii = 1:length(fnames)
            movefile([Path,filesep,fold_Path{i},filesep,fnames{ii}],[Path,filesep,fold_Path{i},filesep,fold_Path{i},'_',fnames{ii}])
        end  
        movefile([Path,filesep,fold_Path{i},filesep,'*'],Raw_Path)
        rmdir([Path,filesep,fold_Path{i}])
    elseif contains(Path,'ESC') %ESC3 has a lot of folders
        collapseFolder(Path,fold_Path{i},1)
    else
        disp(['Unclear how to process folder: ',fold_Path{i}])
    end
end
%% Handle files in the main path 
file_Path = extractfield(dir(Path),'name',~extractfield(dir(Path),'isdir'));
file_Path(contains(file_Path,{'thumbs.db','DS_Store'})) = [];
if ~isempty(file_Path)
    num_dash = cellfun(@length,strfind(file_Path,'_')); %GNO files have exactly 7 _ in the file names
    file_check = {'Processed','Lateral','LARP','RALP'}; %GNO pre-processed files, ready to be moved
    GNO_ext = {'.txt','.xml','.csv','.avi'};
    if contains(Path,'ESC') %detect ESC files by this convention
        % For the 2016-2022 files, the data files don't yet have the .mat
        % extension needed to open them
        for j = 1:length(file_Path)
            rel_file = file_Path{j};
            ext = rel_file(find(rel_file=='.',1,'last')+1:end);
            if ~all(ismember(lower(ext),'abcdefghijklmnopqrstuvwxyz')) %not valid file extenstion
                movefile([Path,filesep,rel_file],[Path,filesep,rel_file,'.mat'])
            end
        end
        rel_file = extractfield(dir(Path),'name',~extractfield(dir(Path),'isdir'));
        for j = 1:length(rel_file)
            movefile([Path,filesep,rel_file{j}],[Raw_Path,filesep,rel_file{j}])
        end
    elseif any((num_dash==7|(num_dash==8&contains(file_Path,file_check)))&contains(file_Path,GNO_ext))||contains(Path,'GNO')
        rel_file = file_Path(num_dash==7&contains(file_Path,GNO_ext)); %Just focus on splitting the files first
        rel_file(contains(rel_file,'.avi')) = []; %Nothing to do with videos
        disp('Splitting GNO files.')
        for j = 1:length(rel_file)
            disp(['File ',num2str(j),'/',num2str(length(rel_file))])
            splitGNOfile(Path,rel_file{j},0) %Already deidentified above
            movefile([Path,filesep,rel_file{j}],[Path,filesep,'Processed_',rel_file{j}]) %Add this prefix to not redo files that have already been processed
        end
        disp('Done splitting GNO files.')
        %Rerun this part to find new files on the path that need to be moved
        file_Path = extractfield(dir(Path),'name',~extractfield(dir(Path),'isdir'));
        if any(contains(file_Path,'Processed')) %Move processed files into a folder and move it
            if ~isfolder([Raw_Path,filesep,'Combined Files'])
                mkdir([Raw_Path,filesep,'Combined Files'])
            end        
            movefile([Path,filesep,'Processed*'],[Raw_Path,filesep,'Combined Files'])
        end
        file_Path(contains(file_Path,'Processed')) = []; 
        for j = 1:length(GNO_ext)
            if any(contains(file_Path,GNO_ext{j}))
                movefile([Path,filesep,'*',GNO_ext{j}],[Path,filesep,'Raw Files'])
            end
        end
    else %Move all non .pdf files
        rel_file = file_Path(~contains(file_Path,'.pdf'));
        for j = 1:length(rel_file)
            movefile([Path,filesep,rel_file{j}],[Raw_Path,filesep,rel_file{j}])
        end
    end
end
%% Create Notes Files
flag = MakeNotes(Raw_Path);
close all;
end