function tab = assessProgress(path)
    sep = filesep;
    a = dir(path);
    a = {a.name};
    %Make sure the needed folders are in the directory
    if ~any(ismember(a,'Segmented Files'))
        error('No "Segmented Files" folder in this directory')
    elseif ~any(ismember(a,'Cycle Averages'))
        error('No "Cycle Averages" folder in this directory')    
    end
    b = dir([path,sep,'Segmented Files',sep,'*.mat']);
    if isempty(b)
        error('No segments present.')
    else
        b = {b.name};
    end
    % Now remove unneccesary info for ease of reading
    c = cell(length(b),1);
    for i = 1:length(b)
        r_str = b{i};
        pos = find(ismember(r_str,'-') == 1);
        c{i} = r_str(pos(2)+1:end);
    end
    %Things that have been looked at
    d1 = dir([path,sep,'Cycle Averages',sep,'*.mat']);
    if ~isempty(d1)
        d1 = {d1.name}';
        % Now remove unneccesary info for ease of reading
        e1 = cell(length(d1),1);
        for i = 1:length(d1)
            r_str = d1{i};
            pos = find(ismember(r_str,'-') == 1);
            e1{i} = strrep(r_str(pos(2)+1:end),'_CycleAvg','');
        end
    else
        e1 = {''};
    end
    %Things that couldn't be analyzed
    d2 = dir([path,sep,'Cycle Averages',sep,'Not*.mat']);
    if ~isempty(d2)
        d2 = {d2.name}';
        % Now remove unneccesary info for ease of reading
        e2 = cell(length(d2),1);
        for i = 1:length(d2)
            r_str = d2{i};
            pos = find(ismember(r_str,'-') == 1);
            e2{i} = strrep(r_str(pos(2)+1:end),'_CycleAvg','');
        end
    else
        e2 = {''};
    end
    %Just make the whole table from scratch each time based one what is in
    %segments
    f = ismember(c,e1);
    g = ismember(c,e2);
    tab = [cell2table(c),array2table([f,g])];
    tab.Properties.VariableNames = {'Segment','CycAvgFile','NotAnalyzeable'};    
    save([path,sep,'Progress.mat'],'tab')
end