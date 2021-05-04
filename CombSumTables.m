% For one subject's Rotary Chair results
exp_type = 'Rotary Chair';
Path = cd; %run from that subject's folder
path_parts = split(Path,filesep);
subject = strrep(path_parts{contains(path_parts,'MVI')},'_','');
visit_folds = extractfield(dir(Path),'name',extractfield(dir(Path),'isdir')&contains(extractfield(dir(Path),'name'),'Visit'));
tabs = cell(length(visit_folds),1);
for i = 1:length(visit_folds)    
    rel_fold = extractfield(dir([Path,filesep,visit_folds{i},filesep,exp_type,filesep,'*Results.mat']),'name');
    if ~isempty(rel_fold)
        fname = [Path,filesep,visit_folds{i},filesep,exp_type,filesep,rel_fold{end}]; %get the most recent item
        load(fname,'all_results')
        tabs{i} = all_results;
    end
end
all_results = sortrows(vertcat(tabs{:}),'Date','ascend');
save([Path,filesep,datestr(now,'yyyymmdd_HHMMSS'),'_',subject,'_',strrep(exp_type,' ',''),'Results.mat'],'all_results')