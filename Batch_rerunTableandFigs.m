if ~ispc %AIA Mac
    drive_path = '/Volumes/vnelhuman$/MVI/Study Subjects/';
else %AIA Lab Computer
    drive_path = '\\win.ad.jhu.edu\cloud\vnelhuman$\MVI\Study Subjects\';
end
folders = strrep(strcat({drive_path},table2cell(readtable('Directories.txt','Delimiter',newline,'ReadVariableNames', false))),'/',filesep); %should be in the user path
for i = 1:length(folders)
    disp([num2str(i),'/',num2str(length(folders)),': ',folders{i}])
    RerunTableandFigs(folders{i})
    close all;
end