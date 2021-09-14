if contains(cd,'Cycle Averages')
    Path = cd;
else
    Path = [cd,filesep,'Cycle Averages'];
end
fnames = extractfield(dir([Path,filesep,'*.mat']),'name');
for i = 1:length(fnames)
   disp(fnames{i})
   load([Path,filesep,fnames{i}],'CycAvg')
   disp(find(~CycAvg.keep_tr))
end