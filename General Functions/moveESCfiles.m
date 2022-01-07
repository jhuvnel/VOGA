function moveESCfiles(Path)
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
end