function deidentify_filenames(Path,rm_string)
names = extractfield(dir(Path),'name');
names = names(contains(names,rm_string)); 
if ~isempty(names) %Actual file names to change
    for i = 1:length(names)
        movefile(names{i},strrep(names{i},rm_string,'DEIDENTIFIED_SUBJECT'));
    end
else %Look in .mat files for raw files to rename
    names = extractfield(dir([Path,filesep,'*.mat']),'name');
    for i = 1:length(names)
        a = load(names{i});
        b = fieldnames(a);
        if strcmp(b{1},'Data') %Segment
            Data = a.('Data');
            Data.rawfile = strrep(Data.rawfile,rm_string,'DEIDENTIFIED_SUBJECT');
            Data.info.rawnotes = strrep(Data.info.rawnotes,rm_string,'DEIDENTIFIED_SUBJECT');
            Data.info.rawfile = strrep(Data.info.rawfile,rm_string,'DEIDENTIFIED_SUBJECT');
            save(names{i},'Data');
        elseif strcmp(b{1},'CycAvg') %CycAvg
            CycAvg = a.('CycAvg');
            CycAvg.info.rawnotes = strrep(CycAvg.info.rawnotes,rm_string,'DEIDENTIFIED_SUBJECT');
            CycAvg.info.rawfile = strrep(CycAvg.info.rawfile,rm_string,'DEIDENTIFIED_SUBJECT');
            save(names{i},'CycAvg');
        else
            disp(['Unregonized data type in: ',names{i}])
        end
    end
end
end