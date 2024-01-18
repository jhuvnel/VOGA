function deidentify_filenames(Path,rm_string)
names = extractfield(dir([Path,filesep,'*',rm_string,'*']),'name');
for i = 1:length(names)
    movefile([Path,filesep,names{i}],[Path,filesep,strrep(names{i},rm_string,'DEIDENTIFIED_SUBJECT')]);
end
end