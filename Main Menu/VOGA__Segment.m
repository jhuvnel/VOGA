%% VOGA__Segment
%
% This function expects that all relevant VOG files are already tagged with
% their metadata in the Raw Files folder. It looks for all VOG data types
% and the manually_segment parameter allows the user to select files to
% manually segment.
%
function VOGA__Segment(Path,manually_segment)
if nargin < 1
    Path = cd;
end
if nargin < 2
    manually_segment = 0;
end
if ~MakeFolders(Path,0)
    return;
end
Raw_Path = [Path,filesep,'Raw Files']; Seg_Path = [Path,filesep,'Segmented Files'];
%Keywords to recognize a file as a VOG file
VOG_fname_pat = {'SESSION','Lateral.txt','LARP.txt','RALP.txt','.dat','.mat','ImuData'};
file_names = extractfield(dir(Raw_Path),'name',find(~extractfield(dir(Raw_Path),'isdir')));
if isempty(file_names)
    disp(['No files in this path: ',Raw_Path])
    return;
end
Notes_ind = contains(file_names,'-Notes.txt'); %Find which .txt files have metadata
VOG_ind = contains(file_names,VOG_fname_pat)&~Notes_ind&~contains(file_names,{'Raw','LDHP','LDPC'}); %Isolate files with data
VOG_ind_num = find(VOG_ind);
has_notes = contains(file_names(VOG_ind),strrep(file_names(Notes_ind),'-Notes.txt',''));
VOG_files = file_names(VOG_ind_num(has_notes)); %Find which VOG files have associated metatdata files
if ~any(VOG_ind)
    disp(['No LDVOG, NL, GNO, or ESC VOG files have been detected: ',Raw_Path])
    return;
elseif isempty(Notes_ind)
    disp(['No Notes files detected in: ',Raw_Path])
    return;
end
indx = 1;
if length(VOG_files)>1
    indx = nmlistdlg('PromptString','Select files to segment:','ListSize',[300 300],'ListString',VOG_files);
end
% Loop over each file
for i = 1:length(indx)
    In_Path = [Raw_Path,filesep,VOG_files{indx(i)}];
    disp([num2str(i),'/',num2str(length(indx)),': ',VOG_files{indx(i)}])
    if manually_segment
        ManuallySegment(In_Path,Seg_Path)
    else
        Segment(In_Path,Seg_Path)
    end
end
end