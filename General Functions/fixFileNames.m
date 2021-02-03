%This script renames CycAvg files with a standard naming convention
files = dir('*.mat');
files = {files.name}';
disp(files)
%% Loop
order = [1,2,3,4,7,5,6,8,9];
for i = 1:length(files)
    fname = files{i};
    load(fname,'CycAvg')
    if any(contains(fieldnames(CycAvg),'name'))
        CycAvg.old_name = fname;    
    end
%     %% Remove suffix
%     if contains(fname,'_Up')
%         ind = strfind(fname,'_Up');
%         fname = fname(1:ind-1);
%     end
%     %% Change names 
%     fname = strrep(fname,'ModON','MotionMod');
%     fname = strrep(fname,'NOSTIM','NoStim');
%     fname = strrep(fname,'Sinusoid','Sine');
%     fname = strrep(fname,'0p','0.');
%     %% Reorder items
%     fparts = strsplit(fname,'-');
%     fname = strjoin(fparts(order),'-');  
%     fname = [fname,'.mat'];
    %% Save
    CycAvg.name = fname;
    delete(files{i})
    save(fname,'CycAvg')
end