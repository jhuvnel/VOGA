%Makes a .txt file with notes from a folder's experiments. To be used as a
%log file, especially for motion-based experiments like Rotary Chair and
%aHIT
function makeTestingLogFile(Raw_Path)
%% Check for VOG files
if nargin < 1
    Raw_Path = cd;
    if isfolder([cd,filesep,'Raw Files'])
        Raw_Path = [cd,filesep,'Raw Files']; 
    end
end
file_names = extractfield([dir([Raw_Path,filesep,'*.txt']);...
    dir([Raw_Path,filesep,'*.dat'])],'name');
if isempty(file_names)
    file_names = {''};
end
VOG_files = file_names(contains(file_names,{'SESSION','.dat'})&...
    ~contains(file_names,{'LDHP','LDPC','-Notes,txt'}));
if isempty(VOG_files)
    disp(['No VOG files found in this directory: ',Raw_Path])
    return;
end
%% Make text to save
w_notes = [{'Files'};strcat(VOG_files,': ')];
%% Save File
filePh = fopen([Raw_Path,filesep,'TestingLogFile.txt'],'w');
fprintf(filePh,'%s\n',w_notes{:});
fclose(filePh);
end