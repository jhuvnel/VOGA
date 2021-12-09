function VOGA__SummaryTable
    %% Run Once
    opts = {'Add to Existing in Folder','Make New in Folder','One Visit','All Visits','All Subject','All Experiment'};    
    [ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
                       'SelectionMode','single',...
                       'ListSize',[150 125],...
                       'ListString',opts); 
    if tf
        if strcmp(opts{ind},'Make New in Folder')
            Path = cd;
            Cyc_Path = [Path,filesep,'Cycle Averages'];
            MakeCycleSummaryTable(Path,Cyc_Path,1);        
        elseif strcmp(opts{ind},'Add to Existing in Folder')
            Path = cd;
            Cyc_Path = [Path,filesep,'Cycle Averages'];
            MakeCycleSummaryTable(Path,Cyc_Path,0);
        elseif strcmp(opts{ind},'One Visit')
            if ~contains(cd,'Visit')
                disp('Path name does not have "Visit" in it. Cannot make a table for this visit')
            else
                Path = cd;
                vis_loc = strfind(Path,'Visit');
                nextdash = find(Path(vis_loc:end)==filesep,1,'first');
                if ~isempty(nextdash)    
                    Path = Path(1:(vis_loc+nextdash-2));
                end
                path_parts = split(Path,filesep);
                subject = strrep(path_parts{contains(path_parts,'MVI')&contains(path_parts,'R')},'_','');
                vis = strrep(path_parts{contains(path_parts,'Visit')},' ','');
                exp_folds = {'Rotary Chair','eeVOR',['vHIT',filesep,'GNO']};%Rotary Chair, eeVOR, vHIT/GNO (change as needed)
                tabs = cell(length(exp_folds),1); 
                for i = 1:length(exp_folds)    
                    rel_fold = extractfield(dir([Path,filesep,exp_folds{i},filesep,'*Results.mat']),'name');
                    if ~isempty(rel_fold) %Has a Results folder
                        fname = [Path,filesep,exp_folds{i},filesep,rel_fold{end}]; %get the most recent item
                        disp(fname)
                        load(fname,'all_results')          
                        tabs{i} = all_results;
                    end
                end
                if ~isempty(vertcat(tabs{:}))
                    all_results = sortrows(vertcat(tabs{:}),'Date','ascend');
                    %Sometimes the first few subject names are the R numbers
                    if length(unique(all_results.Subject))>1
                        subs = unique(all_results.Subject);
                        sub = subs(contains(subs,'MVI'));
                        all_results.Subject = repmat(sub(1),length(all_results.Subject),1);
                    end
                    exp_types = join(unique(all_results.Experiment,'stable'),'-');
                    delete([Path,filesep,'*Results.mat']) %Remove outdated versions
                    save([Path,filesep,subject,'_',vis,'_',exp_types{:},'_Results.mat'],'all_results')
                end
            end    
        elseif strcmp(opts{ind},'All Visits')
            %Assume all visit table exists
            if ~any(contains(split(cd,filesep),'MVI')&contains(split(cd,filesep),'R'))
                disp('Path name does not have a valid MVI***R*** ID in it. Cannot make a table for this subject')
            else
                Path = cd;
                path_parts = split(Path,filesep);
                subject = strrep(path_parts{contains(path_parts,'MVI')&contains(path_parts,'R')},'_','');
                Path = strjoin(path_parts(1:find(contains(path_parts,'MVI')&contains(path_parts,'R'))),filesep);
                visit_folds = extractfield(dir(Path),'name',extractfield(dir(Path),'isdir')&contains(extractfield(dir(Path),'name'),'Visit'));
                for j = 1:length(visit_folds)  
                    vis = strrep(visit_folds{j},' ','');
                    exp_folds = {'Rotary Chair','eeVOR',['vHIT',filesep,'GNO']};%Rotary Chair, eeVOR, vHIT/GNO (change as needed)
                    tabs = cell(length(exp_folds),1); 
                    for i = 1:length(exp_folds)    
                        rel_fold = extractfield(dir([Path,filesep,visit_folds{j},filesep,exp_folds{i},filesep,'*Results.mat']),'name');
                        if ~isempty(rel_fold) %Has a Results folder
                            fname = [Path,filesep,visit_folds{j},filesep,exp_folds{i},filesep,rel_fold{end}]; %get the most recent item
                            disp(fname)
                            load(fname,'all_results') 
                            tabs{i} = all_results;
                        end
                    end
                    try
                        vertcat(tabs{:});
                    catch
                        disp(tabs)
                        error('Unequal table size')
                    end
                    if ~isempty(vertcat(tabs{:}))
                        all_results = sortrows(vertcat(tabs{:}),'Date','ascend');
                        %Sometimes the first few subject names are the R numbers
                        if length(unique(all_results.Subject))>1
                            subs = unique(all_results.Subject);
                            sub = subs(contains(subs,'MVI'));
                            all_results.Subject = repmat(sub(1),length(all_results.Subject),1);
                        end
                        exp_types = join(unique(all_results.Experiment,'stable'),'-');
                        delete([Path,filesep,visit_folds{j},filesep,'*Results.mat']) %Remove outdated versions
                        save([Path,filesep,visit_folds{j},filesep,subject,'_',vis,'_',exp_types{:},'_Results.mat'],'all_results')
                    end
                end
                Path = cd;
                path_parts = split(Path,filesep);
                subject = strrep(path_parts{contains(path_parts,'MVI')&contains(path_parts,'R')},'_','');
                Path = strjoin(path_parts(1:find(contains(path_parts,'MVI')&contains(path_parts,'R'))),filesep);
                visit_folds = extractfield(dir(Path),'name',extractfield(dir(Path),'isdir')&contains(extractfield(dir(Path),'name'),'Visit'));
                tabs = cell(length(visit_folds),1);
                for i = 1:length(visit_folds)    
                    rel_fold = extractfield(dir([Path,filesep,visit_folds{i},filesep,'*Results.mat']),'name');
                    if ~isempty(rel_fold) %Has a Results folder
                        fname = [Path,filesep,visit_folds{i},filesep,rel_fold{end}]; %get the most recent item
                        load(fname,'all_results')             
                        tabs{i} = all_results;
                    end
                end
                if ~isempty(vertcat(tabs{:}))
                    all_results = sortrows(vertcat(tabs{:}),'Date','ascend');
                    %Sometimes the first few subject names are the R numbers
                    if length(unique(all_results.Subject))>1
                        subs = unique(all_results.Subject);
                        sub = subs(contains(subs,'MVI'));
                        all_results.Subject = repmat(sub(1),length(all_results.Subject),1);
                    end
                    exp_types = join(unique(all_results.Experiment,'stable'),'-');
                    delete([Path,filesep,'*Results.mat']) %Remove outdated versions
                    save([Path,filesep,subject,'_',exp_types{:},'_Results.mat'],'all_results')
                end
            end  
        elseif strcmp(opts{ind},'All Subject')
            %Assume all visit table exists
            if ~any(contains(split(cd,filesep),'MVI')&contains(split(cd,filesep),'R'))
                disp('Path name does not have a valid MVI***R*** ID in it. Cannot make a table for this subject')
            else
                Path = cd;
                path_parts = split(Path,filesep);
                subject = strrep(path_parts{contains(path_parts,'MVI')&contains(path_parts,'R')},'_','');
                Path = strjoin(path_parts(1:find(contains(path_parts,'MVI')&contains(path_parts,'R'))),filesep);
                visit_folds = extractfield(dir(Path),'name',extractfield(dir(Path),'isdir')&contains(extractfield(dir(Path),'name'),'Visit'));
                tabs = cell(length(visit_folds),1);
                for i = 1:length(visit_folds)    
                    rel_fold = extractfield(dir([Path,filesep,visit_folds{i},filesep,'*Results.mat']),'name');
                    if ~isempty(rel_fold) %Has a Results folder
                        fname = [Path,filesep,visit_folds{i},filesep,rel_fold{end}]; %get the most recent item
                        load(fname,'all_results')             
                        tabs{i} = all_results;
                    end
                end
                if ~isempty(vertcat(tabs{:}))
                    all_results = sortrows(vertcat(tabs{:}),'Date','ascend');
                    %Sometimes the first few subject names are the R numbers
                    if length(unique(all_results.Subject))>1
                        subs = unique(all_results.Subject);
                        sub = subs(contains(subs,'MVI'));
                        all_results.Subject = repmat(sub(1),length(all_results.Subject),1);
                    end
                    exp_types = join(unique(all_results.Experiment,'stable'),'-');
                    delete([Path,filesep,'*Results.mat']) %Remove outdated versions
                    save([Path,filesep,subject,'_',exp_types{:},'_Results.mat'],'all_results')
                end
            end           
        elseif strcmp(opts{ind},'All Experiment')
            if ispc %AIA lab machine
                Path = 'Z:\Study Subjects';
            else %AIA Mac
                Path = '/Volumes/vnelhuman$/MVI/Study Subjects';
            end
            sub_folds = extractfield(dir(Path),'name',extractfield(dir(Path),'isdir')&contains(extractfield(dir(Path),'name'),'MVI'));
            tabs = cell(length(sub_folds),1);
            for i = 1:length(sub_folds)    
                rel_fold = extractfield(dir([Path,filesep,sub_folds{i},filesep,'*Results.mat']),'name');
                if ~isempty(rel_fold)
                    %Load most recent version
                    fname = [Path,filesep,sub_folds{i},filesep,rel_fold{end}]; %get the most recent item
                    load(fname,'all_results')
                    %Run again
                    tabs{i} = all_results;
                end
            end
            if ~isempty(vertcat(tabs{:}))
                all_results_allexp = vertcat(tabs{:});
                delete([Path,filesep,'*Results.mat']) %Remove outdated versions
                exp_types = unique(all_results_allexp.Experiment,'stable');
                for i = 1:length(exp_types)
                    all_results = all_results_allexp(contains(all_results_allexp.Experiment,exp_types{i}),:);
                    save([cd,filesep,datestr(now,'yyyymmdd_HHMMSS'),'_AllSubjects_',strrep(exp_types{i},' ',''),'Results.mat'],'all_results')
                end
            end
        end    
    end
end