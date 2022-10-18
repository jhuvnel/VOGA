function [all_results,cyc_params,fname] = CombineSummaryTables(Path,type)
%% One Visit
if contains(type,'One Visit')
    vis_loc = strfind(Path,'Visit');
    if isempty(vis_loc)
        error(['Path name does not have "Visit" in it. Cannot make a table for this visit: ',Path])
    end
    nextdash = find(Path(vis_loc:end)==filesep,1,'first');
    if ~isempty(nextdash)
        Path = Path(1:(vis_loc+nextdash-2));
    end
    path_parts = split(Path,filesep);
    subject = strrep(path_parts{contains(path_parts,'MVI')&contains(path_parts,'R')},'_','');
    vis = strrep(path_parts{contains(path_parts,'Visit')},' ','');
    exp_folds = {'Rotary Chair','eeVOR',['vHIT',filesep,'GNO']};%Rotary Chair, eeVOR, vHIT/GNO (change as needed)
    tabs = cell(length(exp_folds),1);
    params = cell(length(exp_folds),1);
    for i = 1:length(exp_folds)
        if contains(type,'New') %Rerun tables and param files if there are Cycle Averages
            cyc_files = extractfield(dir([Path,filesep,exp_folds{i},filesep,'Cycle Averages',filesep,'*.mat']),'name');
            if ~isempty(cyc_files)
                MakeCycleSummaryTable([Path,filesep,exp_folds{i}],[Path,filesep,exp_folds{i},filesep,'Cycle Averages'],1);  
            end
        end
        rel_fold = extractfield(dir([Path,filesep,exp_folds{i},filesep,'*Results.mat']),'name');
        if ~isempty(rel_fold) %Has a Results table
            fname = [Path,filesep,exp_folds{i},filesep,rel_fold{end}]; %get the most recent item
            disp(fname)
            tab_var = load(fname);
            if isfield(tab_var,'all_results')
                tabs{i} = tab_var.all_results;
                rel_fold2 = extractfield(dir([Path,filesep,exp_folds{i},filesep,'*Param.mat']),'name');
                if isempty(rel_fold2) %Does not have a param table--run again
                    try %look for cycle params and if it doesn't exist, rerun the file
                        MakeCycleSummaryTable([Path,filesep,exp_folds{i}],[Path,filesep,exp_folds{i},filesep,'Cycle Averages'],0);
                    catch 
                        MakeCycleSummaryTable([Path,filesep,exp_folds{i}],[Path,filesep,exp_folds{i},filesep,'Cycle Averages'],1);
                    end
                    rel_fold2 = extractfield(dir([Path,filesep,exp_folds{i},filesep,'*Param.mat']),'name');
                end
                fname = [Path,filesep,exp_folds{i},filesep,rel_fold2{end}]; %get the most recent item
                disp(fname)
                load(fname,'cyc_params')
                params{i} = cyc_params;
            end
        end
    end
    if ~isempty(vertcat(tabs{:}))
        [all_results,ord] = sortrows(vertcat(tabs{:}),'Date','ascend');
        cyc_params = vertcat(params{:});
        cyc_params = cyc_params(ord,:);
        %Sometimes the first few subject names are the R numbers
        if length(unique(all_results.Subject))>1
            subs = unique(all_results.Subject);
            sub = subs(contains(subs,'MVI'));
            all_results.Subject = repmat(sub(end),length(all_results.Subject),1);
        end
        fname = [Path,filesep,subject,'_',vis,'_VOG'];
        delete([Path,filesep,'*Results.mat']) %Remove outdated versions
        save([fname,'Results.mat'],'all_results')
        delete([Path,filesep,'*Param.mat']) %Remove outdated versions
        save([fname,'CycParam.mat'],'cyc_params')
    else
        all_results = [];
        cyc_params = [];
        fname = [];
    end
elseif contains(type,'All Visits')
    path_parts = split(Path,filesep);
    sub_part = contains(path_parts,'MVI')&contains(path_parts,'R');
    if isempty(sub_part)
        error('Path name does not have a valid MVI***R*** ID in it. Cannot make a table for this subject')
    end
    subject = strrep(path_parts{sub_part},'_','');
    Path = strjoin(path_parts(1:find(contains(path_parts,'MVI')&contains(path_parts,'R'))),filesep);
    visit_folds = extractfield(dir(Path),'name',extractfield(dir(Path),'isdir')&contains(extractfield(dir(Path),'name'),'Visit'));
    tabs = cell(length(visit_folds),1);
    params = cell(length(visit_folds),1);
    for i = 1:length(visit_folds) 
        if contains(type,'New')
            [all_results,cyc_params] = CombineSummaryTables([Path,filesep,visit_folds{i}],'One Visit (New)');
        else
            [all_results,cyc_params] = CombineSummaryTables([Path,filesep,visit_folds{i}],'One Visit (Existing)');
        end
        if ~isempty(all_results)             
            tabs{i} = all_results;
        end
        if ~isempty(cyc_params)             
            params{i} = cyc_params;
        end
    end
    if ~isempty(vertcat(tabs{:}))
        [all_results,ord] = sortrows(vertcat(tabs{:}),'Date','ascend');
        cyc_params = vertcat(params{:});
        cyc_params = cyc_params(ord,:);
        %Sometimes the first few subject names are the R numbers
        if length(unique(all_results.Subject))>1
            subs = unique(all_results.Subject);
            sub = subs(contains(subs,'MVI'));
            all_results.Subject = repmat(sub(1),length(all_results.Subject),1);
        end
        fname = [Path,filesep,subject,'_VOG'];
        delete([Path,filesep,'*Results.mat']) %Remove outdated versions
        save([fname,'Results.mat'],'all_results')
        delete([Path,filesep,'*Param.mat']) %Remove outdated versions
        save([fname,'CycParam.mat'],'cyc_params')
    end
elseif contains(type,'All Subjects')    
    sub_folds = extractfield(dir(Path),'name',extractfield(dir(Path),'isdir')...
        &contains(extractfield(dir(Path),'name'),'MVI')...
        &contains(extractfield(dir(Path),'name'),'_R'));
    if isempty(sub_folds)
        error('No subject folders found. Check the path (should be in "Study Subjects") and try again.')
    end
    tabs = cell(length(sub_folds),1);
    params = cell(length(sub_folds),1);
    for i = 1:length(sub_folds) 
        if contains(type,'New')
            [all_results,cyc_params] = CombineSummaryTables([Path,filesep,sub_folds{i}],'All Visits (New)');
        else
            [all_results,cyc_params] = CombineSummaryTables([Path,filesep,sub_folds{i}],'All Visits (Existing)');
        end
        if ~isempty(all_results)             
            tabs{i} = all_results;
        end
        if ~isempty(cyc_params)             
            params{i} = cyc_params;
        end
    end
    if ~isempty(vertcat(tabs{:}))
        all_results = vertcat(tabs{:});
        cyc_params = vertcat(params{:});
        fname = [Path,filesep,'ALLMVI-VOG'];
        delete([Path,filesep,'*VOGResults.mat']) %Remove outdated versions
        save([fname,'Results.mat'],'all_results')
    end
end
end