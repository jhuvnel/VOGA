function end_flag = VOGA__ProcessRawData
if ~VOGA__checkFolders(1)
    disp('VOGA instance ended.') %Only way to get here is if you hit "Exit" on the dialogue
    end_flag = true;
    return;
else
    end_flag = false;
end
%Make Notes Files, edit notes files and edit VOG trigger if needed
Path = cd;
Raw_Path = [cd,filesep,'Raw Files'];
%% Move Files into a Workable State
%% NKI/NL
%Transfer NKI Raw Files from their subfolders if they exist
NKIfold = extractfield(dir(Raw_Path),'name',find(extractfield(dir(Raw_Path),'isdir')));
NKIfold(contains(NKIfold,{'.','trash','Trash','TRASH','archive','Archive','ARCHIVE'})) = []; %ignore the invisible directories and anything labelled trash/archive
if isempty(NKIfold) %No NKI folders
    return;
end
for i = 1:length(NKIfold)
    dat_file = extractfield(dir([Raw_Path,filesep,NKIfold{i},filesep,'*Test.dat']),'name');
    avi_file = extractfield(dir([Raw_Path,filesep,NKIfold{i},filesep,'*Test.avi']),'name');
    for j = 1:length(dat_file)
        movefile([Raw_Path,filesep,NKIfold{i},filesep,dat_file{j}],[Raw_Path,filesep,NKIfold{i},'_',dat_file{j}]);
    end
    for j = 1:length(avi_file)
        movefile([Raw_Path,filesep,NKIfold{i},filesep,avi_file{j}],[Raw_Path,filesep,NKIfold{i},'_',avi_file{j}]);
    end
end 
%% Split and move GNO files into Raw Files
names = extractfield(dir(Path),'name')';
isfold = extractfield(dir(Path),'isdir')';
% See if there are any relevant files/folders at all on the path
if ~any(contains(names,{'.txt','.xml','.csv','.avi','video'}))
    disp('No GNO file types detected.')
    return;
end    
deidentify_dat = questdlg('Do the file names need to be de-identified?','','Yes','No','No');
if isempty(deidentify_dat)
    disp('Ending GNO file transfer process.')
    return;
elseif strcmp(deidentify_dat,'Yes')
    deidentify = 1;
    rm_string = inputdlg('Enter the subject name to deidentify (including the _ in the middle):');
    if isempty(rm_string)|| isempty(rm_string{:})
        disp('Ending GNO file transfer process.')
        return;
    else
        rm_string = rm_string{:};
    end
else
    deidentify = 0;
end
% Move video files
if any(contains(names,'video')&isfold)
    dir_names = names(contains(names,'video')&isfold);
    disp('Moving videos now')
    for i = 1:length(dir_names)
        cd([Path,filesep,dir_names{i}])
        movefile('*',Path) %Move everything to the root directory
        cd(Path)
        rmdir([Path,filesep,dir_names{i}])
    end
elseif any(contains(names,'.avi'))  
    disp('No video folder found but video files have been found.')
else
    disp('No video folder or files are present')
end
% Break up text/xml/csv files
file_check = {'Processed','Lateral','LARP','RALP'};
all_ext = {'.txt','.xml','.csv','.avi'};
for j = 1:3
    ext = all_ext{j};
    if any(contains(names,ext)&~contains(names,file_check))
        disp(['Parsing ',ext,' files now...'])
        fnames = names(contains(names,ext)&~contains(names,file_check))';
        for i = 1:length(fnames)
            fname = fnames{i};
            splitGNOfile(Path,fname,deidentify)
            movefile(fname,['Processed_',fname]) %Add this prefix to not redo files that have already been processed
        end
    elseif any(contains(names,ext)&contains(names,file_check))
        disp(['All ',ext,' files appear to have already been processed and are being moved to Raw Files'])    
    else
        disp(['No ',ext,' files found'])
    end
end
% Move files    
names = extractfield(dir(Path),'name')';
if any(contains(names,'Processed'))
    movefile('Processed*','Raw Files')
    cd('Raw Files')
    if exist('Combined Files','dir')~=7
        mkdir('Combined Files')
    end
    movefile('Processed*','Combined Files')
    cd ../
end      
if deidentify
    deidentify_filenames(Path,rm_string)
end
names = extractfield(dir(Path),'name')';
for j = 1:length(all_ext)
    if any(contains(names,all_ext{j}))
        movefile(['*',all_ext{j}],'Raw Files')
    end
end 
%% Move ESC files into Raw Files
names = extractfield(dir(Path),'name')';
if ~any(contains(names,{'_export'}))
    disp('No ESC file types detected.')
    return;
end    
deidentify_dat = questdlg('Do the file names need to be de-identified?','','Yes','No','No');
if isempty(deidentify_dat)
    disp('Ending GNO file transfer process.')
    return;
elseif strcmp(deidentify_dat,'Yes')
    deidentify = 1;
    rm_string = inputdlg('Enter the subject name to deidentify (including the _ in the middle):');
    if isempty(rm_string)|| isempty(rm_string{:})
        disp('Ending GNO file transfer process.')
        return;
    else
        rm_string = rm_string{:};
    end
else
    deidentify = 0;
end
if deidentify
    deidentify_filenames(Path,rm_string)
end
names = extractfield(dir(Path),'name')';
for j = 1:length(names)
    if contains(names{j},'export')
        movefile(names{j},[names{j},'.mat'])
    end
end 
movefile('*.mat','Raw Files')
movefile('*.pdf','Raw Files')
%% If GNO/ESC, open PDF for notes, close manually later
pdf_files = extractfield(dir([cd,filesep,'*.pdf']),'name');
for i = 1:length(pdf_files)
    open(pdf_files{i})
end
%% Create Notes Files
%Detect and process log files (eeVOR)
logtoNotes(Raw_Path)
%Make notes files in other ways
MakeNotes(Raw_Path)
end