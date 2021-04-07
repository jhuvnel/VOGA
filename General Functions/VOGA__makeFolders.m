function VOGA__makeFolders
    %See if the folders already exist or need to be renamed or created
    path_folders = extractfield(dir,'name',find(extractfield(dir,'isdir')));
    if any(contains(path_folders,'Raw LD VOG Files'))
        movefile('Raw LD VOG Files','Raw Files')
    elseif ~any(contains(path_folders,'Raw Files'))
        mkdir('Raw Files')
    end
    if any(contains(path_folders,'CycAvg')) 
        movefile('CycAvg','Cycle Averages')
    elseif any(contains(path_folders,'Cyc_Avg')) 
        movefile('Cyc_Avg','Cycle Averages')
    elseif ~any(contains(path_folders,'Cycle Averages'))
        mkdir('Cycle Averages')
    end
    if ~any(contains(path_folders,'Segmented Files'))
        mkdir('Segmented Files')
    end
end

