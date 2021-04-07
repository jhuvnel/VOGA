function moveGNOfiles(Path)
    names = extractfield(dir(Path),'name')';
    isfold = extractfield(dir(Path),'isdir')';
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
        movefile('*.avi','Raw Files') 
    elseif any(contains(names,'.avi'))  
        disp('No video folder found but video files are being moved to Raw Files')
        movefile('*.avi','Raw Files') 
    else
        disp('No video folder or files are present')
    end
    %% Break up text/xml/csv files
    file_check = {'Processed'};
    all_ext = {'.txt','.xml','.csv'};
    for j = 1:3
        ext = all_ext{j};
        if any(contains(names,ext)&~contains(names,file_check))
            disp(['Parsing ',ext,' files now...'])
            fnames = names(contains(names,ext)&~contains(names,file_check));
            for i = 1:length(fnames)
                splitGNOfile(Path,fnames{i})
                movefile(fnames{i},['Processed_',fnames{i}]) %Add this prefix to not redo files that have already been processed
            end
        elseif any(contains(names,ext)&contains(names,file_check))
            disp(['All ',ext,' files appear to have already been processed and are being moved to Raw Files'])    
        else
            disp(['No ',ext,' files found'])
        end
    end
    %% Move files
    if any(contains(names,all_ext)&contains(names,file_check))
        mkdir('Combined Files')
        movefile('Processed*','Combined Files')
        movefile('Combined Files','Raw Files')
        for j = 1:length(all_ext)
            if any(contains(names,all_ext{j})&contains(names,file_check))
                movefile(['*',all_ext{j}],'Raw Files')
            end
        end
    end     
end