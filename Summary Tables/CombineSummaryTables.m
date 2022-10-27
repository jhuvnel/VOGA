%% Combine Summary Tables

%Type:
% 'Load'
% 'Rerun' (Will let you pick what to rerun and will load the rest)

function all_results = CombineSummaryTables(type,MVI_path)
%% Set defaults
if nargin < 2 || isempty(MVI_path)
    prompt = 'Select the MVI Study subject root folder.';
    MVI_path = uigetdir(prompt,prompt);
    if ~contains(MVI_path,'MVI')
        disp(['The selected path does not contain the text "MVI", so it may be wrong: ',MVI_path])
    end
end
if nargin < 1
    type = 'Load'; %Make aggregate MVI table from existing tables
end
%% Remake VOGResults.mat
%Run for each subject so that it's feasible to save the cyc_params file
%(too big to save for all subejcts at once).
% Make a new pooled MVI subject file
exps = {'Rotary Chair','eeVOR','vHIT','aHIT'};
MVI_subs = extractfield(dir([MVI_path,filesep,'MVI*_R*']),'name');
subjects = strrep(MVI_subs,'_','');
all_sub_tab = cell(length(subjects),1);
if contains(type,'Rerun') %Will want to remake figures
    data = readtable([userpath,filesep,'VOGA_VerInfo.txt'],'ReadVariableNames',false);
    params.version = data{1,2}{:};
    params.Experimenter = data{2,2}{:};
    sub_info = readtable('SubjectInfo.xlsx');
    params.sub_info = sub_info;
    [exp_indx,tf] = listdlg('PromptString','Select experiment typs to redo.',...
        'ListString',exps,'SelectionMode','multiple','ListSize',[200 300],'InitialValue',1:length(exps));
    if ~tf
        disp('No subjects selected.')
        return;
    end
    [sub_indx,tf] = listdlg('PromptString','Select subjects to redo this experiment type on.',...
        'ListString',subjects,'SelectionMode','multiple','ListSize',[200 300],'InitialValue',1:length(subjects));
    if ~tf
        disp('No subjects selected.')
        return;
    end
    rel_sub = subjects(sub_indx); 
    rel_exp = exps(exp_indx);
else
    rel_sub = {''};
end
for i = 1:length(subjects)
    disp([subjects{i},':'])
    % Make a new pooled MVI subject file with these VOG directories
    VOG_fnames = [dir([MVI_path,filesep,MVI_subs{i},filesep,'Visit*',filesep,'eeVOR',filesep,'*Results.mat']);...
        dir([MVI_path,filesep,MVI_subs{i},filesep,'Visit*',filesep,'Rotary Chair',filesep,'*Results.mat']);...
        dir([MVI_path,filesep,MVI_subs{i},filesep,'Visit*',filesep,'aHIT',filesep,'*Results.mat']);...
        dir([MVI_path,filesep,MVI_subs{i},filesep,'Visit*',filesep,'vHIT',filesep,'GNO',filesep,'*Results.mat']);...
        dir([MVI_path,filesep,MVI_subs{i},filesep,'Visit*',filesep,'vHIT',filesep,'ESC',filesep,'*Results.mat'])];
    one_sub_tab = cell(length(VOG_fnames),1);
    one_sub_params = cell(length(VOG_fnames),1);
    if ~isempty(rel_sub)&&ismember(subjects(i),rel_sub)
        redo_folder = contains({VOG_fnames.folder},rel_exp);
    else
        redo_folder = false(length(VOG_fnames),1);
    end
    for j = 1:length(VOG_fnames)
        disp([num2str(j),'/',num2str(length(VOG_fnames)),': ',VOG_fnames(j).folder])
        cyc_param_fname = extractfield(dir([VOG_fnames(j).folder,filesep,'*CycParam.mat']),'name');
        if redo_folder(j)
            VOGA__makeFolders(VOG_fnames(j).folder);
            [all_results,cyc_params] = MakeCycleSummaryTable(VOG_fnames(j).folder,[VOG_fnames(j).folder,filesep,'Cycle Averages'],1);
            params.Path = VOG_fnames(j).folder;
            plotGroupCycAvg(params);
            plotParamResults(params);
            close all;
        elseif isempty(cyc_param_fname) %Will make the cyc_params struct if missing
            [all_results,cyc_params] = MakeCycleSummaryTable(VOG_fnames(j).folder,[VOG_fnames(j).folder,filesep,'Cycle Averages'],1);
        else %Load
            load([VOG_fnames(j).folder,filesep,VOG_fnames(j).name],'all_results');
            load(strrep([VOG_fnames(j).folder,filesep,cyc_param_fname{:}],'Results','CycParam'),'cyc_params')
        end
        one_sub_tab{j} = all_results;
        one_sub_params{j} = cyc_params;
    end
    all_results = vertcat(one_sub_tab{:});
    all_results.Subject(~strcmp(all_results.Subject,subjects{i})) = subjects(i); %Make sure the subject name is consistent
    save([MVI_path,filesep,MVI_subs{i},filesep,subjects{i},'_VOGResults.mat'],'all_results')
    all_sub_tab{i} = all_results;
    cyc_params = vertcat(one_sub_params{j});
    save([MVI_path,filesep,MVI_subs{i},filesep,subjects{i},'_VOGCycParam.mat'],'cyc_params')
end
all_results = vertcat(all_sub_tab{:});
save([MVI_path,filesep,'ALLMVI-VOGResults.mat'],'all_results')
end