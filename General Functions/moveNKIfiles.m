function moveNKIfiles(Raw_Path)
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
end