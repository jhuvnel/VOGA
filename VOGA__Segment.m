function end_flag = VOGA__Segment
if ~VOGA__checkFolders(0)
    disp('Folder structure not present. Generate folders, and process raw data" first.')
    end_flag = true;
    return;
else
    end_flag = false;
end
Raw_Path = [cd,filesep,'Raw Files'];
Seg_Path = [cd,filesep,'Segmented Files'];
file_names = extractfield(dir(Raw_Path),'name');
file_names(~contains(file_names,{'.txt','.dat','.mat'})) = [];
if isempty(file_names)
    file_names = {''};
end
Notes_ind = contains(file_names,'-Notes.txt');
VOG_ind = contains(file_names,{'SESSION','.dat','Lateral.txt','LARP.txt','RALP.txt','.mat'})&~Notes_ind;
VOG_ind_num = find(VOG_ind);
has_notes = contains(strrep(strrep(file_names(VOG_ind),'.txt',''),'.dat',''),strrep(file_names(Notes_ind),'-Notes.txt',''));
VOG_files = file_names(VOG_ind_num(has_notes));
if ~any(VOG_ind)
    uiwait(msgbox(['No LDVOG, NKI, GNO, or ESC files have been detected in ',Raw_Path]))
    return;
elseif isempty(VOG_files)
    uiwait(msgbox(['No Notes files for the LDVOG, NKI, GNO, or ESC files in ',Raw_Path]))
    return;
end
indx = nmlistdlg('PromptString','Select files to segment:','ListSize',[300 300],'ListString',VOG_files);
sel_files = VOG_files(indx);
%% Loop over each file
for i = 1:length(sel_files)
    In_Path = [Raw_Path,filesep,sel_files{i}];
    disp([num2str(i),'/',num2str(length(sel_files)),': ',sel_files{i}])
    Segment(In_Path,Seg_Path)
end
end