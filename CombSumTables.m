%% One subject's Rotary Chair or eeVOR results
table_acq = 'rerun'; %or load
exp_type = 'Rotary Chair';
Path = cd; %run from that subject's folder
path_parts = split(Path,filesep);
subject = strrep(path_parts{contains(path_parts,'MVI')&contains(path_parts,'R')},'_','');
visit_folds = extractfield(dir(Path),'name',extractfield(dir(Path),'isdir')&contains(extractfield(dir(Path),'name'),'Visit'));
tabs = cell(length(visit_folds),1);
for i = 1:length(visit_folds)    
    rel_fold = extractfield(dir([Path,filesep,visit_folds{i},filesep,exp_type,filesep,'*Results.mat']),'name');
    if ~isempty(rel_fold) %Has a Results folder
        if strcmp(table_acq,'rerun')
            delete([Path,filesep,visit_folds{i},filesep,exp_type,filesep,'*Results.mat']) %Remove outdated versions
            disp([Path,filesep,visit_folds{i},filesep,exp_type])
            all_results = MakeCycleSummaryTable([Path,filesep,visit_folds{i},filesep,exp_type],[Path,filesep,visit_folds{i},filesep,exp_type,filesep,'Cycle Averages'],1);            
            if isempty(all_results)
                 all_results = MakeCycleSummaryTable([Path,filesep,visit_folds{i},filesep,exp_type],[Path,filesep,visit_folds{i},filesep,exp_type,filesep,'LDVOG',filesep,'Cycle Averages'],1);
            end            
        else %load
            fname = [Path,filesep,visit_folds{i},filesep,exp_type,filesep,rel_fold{end}]; %get the most recent item
            load(fname,'all_results')
        end              
        tabs{i} = all_results;
    end
end
all_results = sortrows(vertcat(tabs{:}),'Date','ascend');
%Sometimes the first few subject names are the R numbers
if ~strcmp(all_results.Subject{1},all_results.Subject{end})&&contains(all_results.Subject{end},all_results.Subject{1})
    all_results.Subject(~contains(all_results.Subject,all_results.Subject{end})) = all_results.Subject(end);
end
delete(['*',exp_type,'Results.mat']) %Remove outdated versions
save([Path,filesep,datestr(now,'yyyymmdd_HHMMSS'),'_',subject,'_',strrep(exp_type,' ',''),'Results.mat'],'all_results')
disp('Done')
%% Combine Multiple Subject's tables
exp_type = 'Rotary Chair';
if ispc %AIA lab machine
    drive_name = 'Z:';
else %AIA Mac
    drive_name = '/Volumes/vnelhuman$/MVI';
end
out_Path = [drive_name,filesep,'DATA SUMMARY',filesep,'IN PROGRESS',filesep,'VOR - ',exp_type]; %run from the folder wher you want it saved from 
Path = [drive_name,filesep,'Study Subjects'];
sub_folds = extractfield(dir(Path),'name',extractfield(dir(Path),'isdir')&contains(extractfield(dir(Path),'name'),'MVI'));
tabs = cell(length(sub_folds),1);
for i = 1:length(sub_folds)    
    rel_fold = extractfield(dir([Path,filesep,sub_folds{i},filesep,'*',strrep(exp_type,' ',''),'Results.mat']),'name');
    if ~isempty(rel_fold)
        %Load most recent version
        fname = [Path,filesep,sub_folds{i},filesep,rel_fold{end}]; %get the most recent item
        load(fname,'all_results')
        %Run again
        %
        tabs{i} = all_results;
    end
end
all_results = sortrows(vertcat(tabs{:}),'Date','ascend');
save([out_Path,filesep,datestr(now,'yyyymmdd_HHMMSS'),'_AllSubjects',strrep(exp_type,' ',''),'Results.mat'],'all_results')