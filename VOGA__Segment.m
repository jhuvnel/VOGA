function VOGA__Segment(Path,seg_all)
if nargin < 1
    Path = cd;
end
if nargin < 2
    seg_all = 0;
end
if ~VOGA__makeFolders(Path)
    return;
end
Raw_Path = [Path,filesep,'Raw Files'];
Seg_Path = [Path,filesep,'Segmented Files'];
VOG_fname_pat = {'SESSION','Lateral.txt','LARP.txt','RALP.txt','.dat','.mat','ImuData'};
rel_dir = dir(Raw_Path);
rel_dir(extractfield(rel_dir,'isdir')) = [];
if isempty(rel_dir)
    disp(['No files in this path: ',Raw_Path])
    return;
end
file_names = extractfield(rel_dir,'name');
Notes_ind = contains(file_names,'-Notes.txt');
VOG_ind = contains(file_names,VOG_fname_pat)&~Notes_ind&~contains(file_names,{'Raw','LDHP','LDPC'});
VOG_ind_num = find(VOG_ind);
has_notes = contains(file_names(VOG_ind),strrep(file_names(Notes_ind),'-Notes.txt',''));
VOG_files = file_names(VOG_ind_num(has_notes));
if ~any(VOG_ind)
    disp(['No LDVOG, NL, GNO, or ESC VOG files have been detected: ',Raw_Path])
    return;
elseif isempty(VOG_files)
    disp(['No Notes files detected in: ',Raw_Path])
    return;
end
if seg_all
    sel_files = VOG_files;
else
    indx = nmlistdlg('PromptString','Select files to segment:','ListSize',[300 300],'ListString',VOG_files);
    sel_files = VOG_files(indx);
end
%% Loop over each file
for i = 1:length(sel_files)
    In_Path = [Raw_Path,filesep,sel_files{i}];
    disp([num2str(i),'/',num2str(length(sel_files)),': ',sel_files{i}])
    Segment(In_Path,Seg_Path)
end
end