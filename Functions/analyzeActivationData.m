function analyzeActivationData(Cyc_Path)
    files = extractfield(dir([Cyc_Path,filesep,'*.mat']),'name');
    rel_part = cell(length(files),3);
    for i = 1:length(files)
        temp = strsplit(files{i},'-');
        temp2 = strrep(temp{end},'.mat','');
        num = temp2(regexp(temp2,'[0-9]'));        
        rel_part(i,:) = {str2double(num),strrep(temp2,num,''),files{i}};
    end
    rel_part_tab = cell2table(rel_part);
    rel_part_tab.Properties.VariableNames = {'num','cond','fname'};
    rel_part_tab = sortrows(sortrows(rel_part_tab,'cond','ascend'),'num','ascend');
    ordered_files = rel_part_tab.fname;
    load(ordered_files{end},'CycAvg')
    last_t = CycAvg.t(end);
    Fs = CycAvg.Fs;
    full_t = 0:1/Fs:last_t;
    full_stim = NaN(1,length(full_t));
    
end