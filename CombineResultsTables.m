function CombineResultsTables
    find_sub_fold = @(path) strcat(path,filesep,extractfield(dir(path),'name',find(extractfield(dir(path),'isdir')&~contains(extractfield(dir(path),'name'),'.'))));
    sub_fold = find_sub_fold(cd);
    has_file = false(1,length(sub_fold));
    no_fold = false(1,length(sub_fold));
    while(~all(no_fold|has_file))
        has_file = false(1,length(sub_fold));
        no_fold = false(1,length(sub_fold));
        for i = 1:length(sub_fold)
            has_file(i) = any(contains(extractfield(dir(sub_fold{i}),'name'),'Results.mat'));
            no_fold(i) = ~any(extractfield(dir(sub_fold{i}),'isdir')&~contains(extractfield(dir(sub_fold{i}),'name'),'.'));
            if ~no_fold(i)&&~has_file(i)
                sub_fold = [sub_fold;find_sub_fold(sub_fold{i})];
                no_fold(i) = 1;
            end
        end
        has_file = [has_file,false(1,length(sub_fold)-length(has_file))];
        no_fold = [no_fold,false(1,length(sub_fold)-length(no_fold))];
        sub_fold(no_fold) = [];
    end
    all_res = table();
    for i = 1:length(sub_fold)
        files = extractfield(dir(sub_fold{i}),'name');
        res_f = files(contains(files,'Results.mat'));
        load([sub_fold{i},filesep,res_f{end}],'all_results')
        all_res = [all_res;all_results];
    end
    all_results = all_res;
    exps = unique(all_results.Experiment);
    if length(exps)==1
        exp_type = exps{:};
    else
        exp_type = '';
    end
    save([cd,filesep,datestr(now,'yyyymmdd_HHMMSS'),'_',exp_type,'Results.mat'],'all_results')    
end