function moveGNOfiles(Path)
    names = extractfield(dir(Path),'name')';
    isfold = extractfield(dir(Path),'isdir')';
    %% See if there are any relevant files/folders at all on the path
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
    %% Move video files
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
    %% Break up text/xml/csv files
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
    %% Move files
    mkdir('Combined Files')
    names = extractfield(dir(Path),'name')';
    if any(contains(names,'Processed'))
        movefile('Processed*','Combined Files')
    end
    movefile('Combined Files','Raw Files')   
    if deidentify
        deidentify_filenames(Path,rm_string)
    end
    names = extractfield(dir(Path),'name')';
    for j = 1:length(all_ext)
        if any(contains(names,all_ext{j}))
            movefile(['*',all_ext{j}],'Raw Files')
        end
    end    
end