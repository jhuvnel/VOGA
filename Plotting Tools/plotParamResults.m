%% Plot Param Results.m
% This function makes figures for an experiment folder once it has been
% analyzed and the Results.mat table exists.
% Updated to make GroupCycAvg Figures too
% There should be only one subject in the table. Else it returns the
% function.

function plotParamResults(params)
% Initialize
close all;
load('VNELcolors.mat','colors')
all_canals = {'LA';'RP';'LP';'RA';'LH';'RH';'-X';'+X';'-Y';'+Y'};
all_canals(:,2) = reshape(repmat({'LARP','RALP','LHRH','X','Y'},2,1),[],1);
all_canals(:,3) = reshape(repmat({'l','r','z','x','y'},2,1),[],1);
all_markers = split('o d ^ p > h < s v _ * | .');
%Process arguements from params input
Path = params.Path;
Cyc_Path = params.Cyc_Path;
sub_info = params.sub_info;
Subs = sub_info.Subject;
Ears = sub_info.Ear;
% Load table in question
res_file = extractfield(dir([Path,filesep,'*Results.mat']),'name')';
if isempty(res_file)
    disp('No table with cycle parameters found on this path. Creating one now.')
    rerun = ~strcmp(questdlg('If a parameter table already exists, use that one or rerun?','','Use existing table','Rerun','Rerun'),'Use existing table');
    MakeCycleSummaryTable(Path,Cyc_Path,rerun);
    res_file = extractfield(dir([Path,filesep,'*Results.mat']),'name')';
end
load([Path,filesep,res_file{end}],'all_results')
% Confirm there is only one subject (can be mulitple visits)
if length(unique(all_results.Subject))>1
    disp('More than one subject detected from the table on the path provided. Run a different plotting script.')
    return;
end
%Find the implanted canals for this subject
imp_canals = strcat(Ears{ismember(unique(all_results.Subject),Subs)},{'A','P','H'});
all_results.ImpCanal = contains(all_results.AxisName,imp_canals);
all_results.DateStr = cellstr(datestr(all_results.Date,'yyyymmdd'));
%Set up the plot tables
exp_types = {'Sine','HIT','Exp','PulseTrain'};
plot_types = {'CycAvg','MaxVel','Param'};
all_plot_tabs = cell2table(cell(length(exp_types),length(plot_types)),'VariableNames',plot_types,'RowNames',exp_types);
%% Sinusoids
if any(contains(all_results.Type,'Sine'))
    %Isolate relevant table entries and put them in order for plotting by amplitude and frequency
    tab = sortrows(sortrows(all_results(contains(all_results.Type,'Sine')&contains(all_results.AxisName,all_canals(:,1)),:),'Amplitude','ascend'),'Frequency','ascend');
    fn = size(tab,1);
    YMax = 5*ceil(max(tab.MaxVel+tab.MaxVel_sd)/5);
    GainMax = 0.1*ceil(max(tab.Gain+tab.Gain_sd)/0.1);
    %Make one figure for each group of cycle averages (same subject, visit,
    %date, condition, goggle, axis, across frequency, and/or amplitude)
    %Make one figure per for frequency and/or amplitude sweep with MaxVel
    %of each each component for each condition and one figure with Gain and
    %Phase
    YVar = {'Gain';'Phase'};
    YLabs = {'Gain';'Phase Lead (deg)'};
    rel_labs = {'Subject','Visit','DateStr','Condition','Goggle'};
    [~,rel_col] = ismember(rel_labs,tab.Properties.VariableNames);
    file_parts = table2cell(tab(:,rel_col));
    %Parts applicable to each file
    common_cond_i = all(strcmp(file_parts,repmat(file_parts(1,:),fn,1)),1);
    common_cond = strjoin(file_parts(1,common_cond_i));
    common_cond_i(contains(rel_labs,'Condition')) = 0;
    rel_file_parts = join(file_parts(:,~common_cond_i),2);
    [conds,~,IC] = unique(rel_file_parts,'stable');
    enum = length(conds);
    [~,Ic] = ismember(tab.AxisName,all_canals(:,1));
    [freqs,~,IF] = unique(tab.Frequency);
    IF(isnan(freqs)) = 0;
    freqs(isnan(freqs)) = [];
    fnum = length(freqs);
    str_freqs = strcat(strrep(cellstr(num2str(freqs)),' ',''),'Hz');
    [amps,~,IA] = unique(tab.Amplitude);
    IA(isnan(amps)) = 0;
    amps(isnan(amps)) = [];
    anum = length(amps);
    str_amps = strcat(strrep(cellstr(num2str(amps)),' ',''),'dps');
    [freq_amps,~,IFA] = unique([tab.Frequency,tab.Amplitude],'rows','stable');
    fanum = length(freq_amps);
    str_freq_amps = join([strcat(strrep(cellstr(num2str(freq_amps(:,1))),' ',''),'Hz'),...
        strcat(strrep(cellstr(num2str(freq_amps(:,2))),' ',''),'dps')]);
    tab.IFA = IFA;
    %Make a table with all frequency and amplitude sweeps for plotting
    %Find indecies for frequency sweeps with more than 1 file
    [fa_n_mat,~,fa_n_inds] = unique([[IC,Ic,IF,0*IA];[IC,Ic,0*IF,IA]],'rows');
    GC_fa = groupcounts(fa_n_inds);
    fa_n_mat(GC_fa<2,:) = [];
    %Combo of freq and amp like on aHIT
    if isempty(fa_n_mat)
        [fa_n_mat,~,fa_n_inds] = unique([IC,Ic,0*IF,0*IA],'rows');
        GC_fa = groupcounts(fa_n_inds);
        fa_n_mat(GC_fa<2,:) = [];
    end
    if isempty(fa_n_mat)
        cycavg_plots = [];
        maxvel_plots = [];
        param_plots = [];
    else
        dat = array2table(fa_n_mat,'VariableNames',{'Cond_i','Canal_i','Freq_i','Amp_i'});
        for i = 1:size(dat,1)
            inds = find(IC==dat.Cond_i(i)&Ic==dat.Canal_i(i)&...
                (IF==dat.Freq_i(i)|0*IF==dat.Freq_i(i))&...
                (IA==dat.Amp_i(i)|0*IA==dat.Amp_i(i)));
            %Remove duplicate experiments by choosing the one with the
            %closest time stamp to the others
            [gc,gr] = groupcounts(IFA(inds));
            dup = gr(gc>1);
            for ii = 1:length(dup)
                sub_i = find(IFA(inds)==dup(ii));
                [~,k_i] = min(abs(tab.Date(inds(sub_i))-median(tab.Date(inds))));
                sub_i(k_i) = []; %Keep this one
                inds(sub_i) = []; %Delete the others
            end
            if length(inds)>1
                if dat.Freq_i(i)==0&&dat.Amp_i(i)==0 %Freq/Amp
                    freq = '';
                    amp = '';
                    name = [common_cond,' ',conds{dat.Cond_i(i)},...
                        ' ',all_canals{dat.Canal_i(i),2},' Sine'];
                    x_val = IFA(inds);
                    x_var = 'IFA';
                    x_tick = 1:fanum;
                    x_tlab = strrep(str_amps,'dps','');
                    x_lab = 'Amplitude (dps)';
                    x_scale = 'linear';
                    sub_names = str_freq_amps;
                    files = cell(fanum,1);
                elseif dat.Freq_i(i)==0 %Amp
                    freq = '';
                    amp = [str_amps{dat.Amp_i(i)},' Sine'];
                    name = [common_cond,' ',conds{dat.Cond_i(i)},...
                        ' ',all_canals{dat.Canal_i(i),2},' ',str_amps{dat.Amp_i(i)},' Sine'];
                    x_val = IF(inds);
                    x_var = 'Frequency';
                    x_tick = freqs;
                    x_tlab = strrep(str_freqs,'Hz','');
                    x_lab = 'Frequency (Hz)';
                    x_scale = 'log';
                    sub_names = str_freqs;
                    files = cell(fnum,1);
                elseif dat.Amp_i(i)==0 %Freq
                    freq = [str_freqs{dat.Freq_i(i)},' Sine'];
                    amp = '';
                    name = [common_cond,' ',conds{dat.Cond_i(i)},...
                        ' ',all_canals{dat.Canal_i(i),2},' ',str_freqs{dat.Freq_i(i)},' Sine'];
                    x_val = IA(inds);
                    x_var = 'Amplitude';
                    x_tick = amps;
                    x_tlab = strrep(str_amps,'dps','');
                    %Remove labels because of crowding
                    if ismember(200:100:500,amps)
                        x_tlab(contains(x_tlab,{'300','400'})) = {''};
                    elseif ismember(200:100:400,amps)
                        x_tlab(contains(x_tlab,{'300'})) = {''};
                    end
                    x_lab = 'Amplitude (dps)';
                    x_scale = 'log';
                    sub_names = str_amps;
                    files = cell(anum,1);
                end
                files(x_val) = tab.File(inds);
                dat.Cond{i} = conds{dat.Cond_i(i)};
                dat.Canal{i} = all_canals{dat.Canal_i(i),1};
                dat.Freq{i} = freq;
                dat.Amp{i} = amp;
                dat.Name{i} = name;
                dat.Files{i} = files;
                dat.Tab{i} = tab(inds,:);
                dat.SubNames{i} = sub_names;
                dat.YLim{i} = YMax*[-1 1];
                dat.XVar{i} = x_var;
                dat.XTick{i} = x_tick;
                dat.XTickLab{i} = x_tlab;
                dat.XLabel{i} = x_lab;
                dat.XScale{i} = x_scale;
                dat.X{i} = tab.(x_var)(inds);
                for ii = 1:length(YVar)
                    dat.(YVar{ii}){i} = tab.(YVar{ii})(inds);
                    dat.([YVar{ii},'_sd']){i} = tab.([YVar{ii},'_sd'])(inds);
                end
            end
        end
        dat(cellfun(@isempty,dat.Cond),:) = [];
        [~,ind] = unique(dat.Name,'stable');
        cycavg_plots = dat(ind,{'Name','Files','SubNames','YLim'}); %Figure name and cyc avg file names for all
        %Initialize the table for MaxVel plots
        [~,mp_i2,maxvel_plot_i] = unique(dat(:,{'Cond_i','Freq_i','Amp_i'}),'rows','stable');
        maxvel_plots = table();
        maxvel_plots.Name = strrep(strcat({[common_cond,' ']},dat.Cond(mp_i2),...
            {' '},dat.Freq(mp_i2),dat.Amp(mp_i2),{' Velocity'}),'  ',' ');
        for i=1:size(maxvel_plots,1)
            sub_tab = dat(maxvel_plot_i==i,:);
            rel_tab = cell(2,5);
            sub_t = cell(2,5);
            sub_t(sub_tab.Canal_i) = sub_tab.Canal;
            rel_tab(sub_tab.Canal_i) = sub_tab.Tab;
            rel_tab(:,all(cellfun(@isempty,sub_t),1)) = [];
            sub_t(:,all(cellfun(@isempty,sub_t),1)) = [];
            maxvel_plots.SubNames{i} = sub_t;
            maxvel_plots.Tables{i} = rel_tab;
            maxvel_plots.YLim{i} = [0 YMax];
        end
        maxvel_plots(:,{'XVar','XTick','XTickLab','XLabel','XScale'}) = dat(mp_i2,{'XVar','XTick','XTickLab','XLabel','XScale'});
        %Initialize the table for TabParam plots. For Sine, this is Gain/Phase
        %Legends true for them all
        leg_tab = cell2table([all_canals(:,1);conds],'VariableNames',{'Name'});
        leg_tab.Marker(:) = {'none'};
        leg_tab.Marker(11:end) = all_markers(1:enum);
        leg_tab.LineStyle(1:2:9) = {'-'};
        leg_tab.LineStyle(2:2:10) = {':'};
        leg_tab.LineStyle(11:end) = {'none'};
        leg_tab.Color(:) = {[0,0,0]};
        for ii = 1:10
            leg_tab.Color{ii} = colors.(['l_',all_canals{ii,3}]);
        end
        [~,mp_i2,param_plot_i] = unique(dat(:,{'Freq_i','Amp_i'}),'rows','stable');
        param_plots = table();
        param_plots.Name = strrep(strcat({[common_cond,' ']},dat.Freq(mp_i2),dat.Amp(mp_i2),{' Gain'}),'  ',' ');
        for i=1:size(param_plots,1)
            sub_tab = dat(param_plot_i==i,:);
            [sub_t,~,indC] = unique(all_canals(sub_tab.Canal_i,2),'stable');
            rel_tab = cell(length(YVar),length(sub_t));
            rel_leg = cell(length(YVar),length(sub_t));
            for ii = 1:length(sub_t)
                for j = 1:length(YVar)
                    n_tab = table();
                    n_tab.Marker = leg_tab.Marker(sub_tab.Cond_i(indC==ii)+10);
                    n_tab.LineStyle = leg_tab.LineStyle(sub_tab.Canal_i(indC==ii));
                    n_tab.Color = leg_tab.Color(sub_tab.Canal_i(indC==ii));
                    n_tab.X = sub_tab.X(indC==ii);
                    n_tab.Y = sub_tab.(YVar{j})(indC==ii);
                    n_tab.Y_sd = sub_tab.([YVar{j},'_sd'])(indC==ii);
                    rel_tab{j,ii} = n_tab;
                end
                rel_leg{1,ii} = leg_tab(unique(sub_tab.Canal_i(indC==ii)),:);
            end
            %Add legend with conditions
            if size(rel_leg,1)==1 %Only one row
                rel_leg{1,1} = [rel_leg{1,1};leg_tab(11:end,:)];
            else
                rel_leg{end,1} = leg_tab(11:end,:);
            end
            param_plots.SubNames{i} = sub_t;
            param_plots.YLim{i} = [0,GainMax;-180,180];
            param_plots.YLabel{i} = YLabs;
            param_plots.Tables{i} = rel_tab;
            param_plots.Legend{i} = rel_leg;
        end
        param_plots(:,{'XTick','XTickLab','XLabel','XScale'}) = dat(mp_i2,{'XTick','XTickLab','XLabel','XScale'});
    end
    all_plot_tabs.CycAvg{'Sine'} = cycavg_plots;
    all_plot_tabs.MaxVel{'Sine'} = maxvel_plots;
    all_plot_tabs.Param{'Sine'} = param_plots;
end
%% Impulse
if any(contains(all_results.Type,'Impulse'))
    can_lab = {'LA','LP','LH','RA','RP','RH'};
    [~,canal_ord] = ismember(can_lab,all_canals(:,1));
    canals = all_canals(canal_ord,:);
    tab = all_results(contains(all_results.Type,'Impulse'),:);
    [~,Ic] = ismember(tab.AxisName,canals(:,1));
    tab.Ic = Ic;
    tab = sortrows(tab,'Ic','ascend');
    Ic = tab.Ic;
    fn = size(tab,1);
    GainMax = 0.1*ceil(max(tab.Gain+tab.Gain_sd)/0.1);
    %Make one figure for each group of cycle averages (same subject, visit,
    %date, goggle, canal across conditions
    %Make one figure with Gain and Lantecy across canals for each condition
    YVar = {'Gain';'Latency'};
    YLabs = {'Gain';'Latency (ms)'};
    rel_labs = {'Subject','Visit','DateStr','Condition','Goggle'};
    [~,rel_col] = ismember(rel_labs,tab.Properties.VariableNames);
    file_parts = table2cell(tab(:,rel_col));
    %Parts applicable to each file
    common_cond_i = all(strcmp(file_parts,repmat(file_parts(1,:),fn,1)),1);
    common_cond = strjoin(file_parts(1,common_cond_i));
    common_cond_i(contains(rel_labs,'Condition')) = 0;
    rel_file_parts = join(file_parts(:,~common_cond_i),2);
    [conds,~,IC] = unique(rel_file_parts,'stable');
    enum = length(conds);
    [~,~,IFA] = unique([IC,Ic],'rows','stable');
    dat = array2table([1:enum,zeros(1,length(unique(Ic)));zeros(1,enum),reshape(unique(Ic),1,[])]','VariableNames',{'Cond_i','Canal_i'});
    for i = 1:size(dat,1)
        inds = find((IC==dat.Cond_i(i)|0*IC==dat.Cond_i(i))&...
            (Ic==dat.Canal_i(i)|0*Ic==dat.Canal_i(i)));
        %Remove duplicate experiments by choosing the one with the
        %closest time stamp to the others
        [gc,gr] = groupcounts(IFA(inds));
        dup = gr(gc>1);
        for ii = 1:length(dup)
            sub_i = find(IFA(inds)==dup(ii));
            [~,k_i] = min(abs(tab.Date(inds(sub_i))-median(tab.Date(inds))));
            sub_i(k_i) = []; %Keep this one
            inds(sub_i) = []; %Delete the others
        end
        if dat.Cond_i(i)==0 % All conditions in this canal
            cond = '';
            canal = canals{dat.Canal_i(i),1};
            name = [common_cond,' ',canal,' Impulse'];
        elseif dat.Canal_i(i)==0 % All canals in this condition (for param plots)
            cond = conds{dat.Cond_i(i)};
            canal = '';
            name = '';
        end
        dat.Cond{i} = cond;
        dat.Canal{i} = canal;
        dat.Name{i} = name;
        dat.Files{i} = tab.File(inds);
        dat.Tab{i} = tab(inds,:);
        dat.SubNames{i} = conds;
        dat.YLim{i} = [-50 250];
        dat.XVar{i} = 'Ic';
        dat.XTick{i} = 1:6;
        dat.XTickLab{i} = can_lab;
        dat.XLabel{i} = 'Canal';
        dat.XScale{i} = 'linear';
        dat.X{i} = Ic(inds);
        for ii = 1:length(YVar)
            dat.(YVar{ii}){i} = tab.(YVar{ii})(inds);
            dat.([YVar{ii},'_sd']){i} = tab.([YVar{ii},'_sd'])(inds);
        end
    end
    cycavg_plots = dat(~ismember(dat.Name,{''}),{'Name','Files','SubNames','YLim'}); %Figure name and cyc avg file names for all
    maxvel_plots = [];
    %Make the table for TabParam plots. For impulses, this is Gain/Latency (set above)
    %Legends with condition(s)
    leg_tab = cell2table([conds,all_markers(1:enum)],'VariableNames',{'Name','Marker'});
    leg_tab.LineStyle(:) = {'none'};
    leg_tab.Color(:) = {[0,0,0]};
    param_plots = cell2table({[common_cond,' Impulse Gain']},'VariableNames',{'Name'});
    sub_tab = dat(ismember(dat.Name,{''}),:);
    rel_tab = cell(length(YVar),1);
    rel_leg = cell(length(YVar),1);
    for j = 1:length(YVar)
        n_tab = table();
        n_tab.Marker = leg_tab.Marker(sub_tab.Cond_i);
        n_tab.LineStyle(:) = {'-'};
        n_tab.Color(:) = {[0,0,0]};
        n_tab.X = sub_tab.X;
        n_tab.Y = sub_tab.(YVar{j});
        n_tab.Y_sd = sub_tab.([YVar{j},'_sd']);
        rel_tab{j,1} = n_tab;
    end
    %Add legend with conditions
    rel_leg{end,1} = leg_tab;
    param_plots.SubNames{1} = {' '};
    param_plots.YLim{1} = [-0.1,GainMax;0,100];
    param_plots.YLabel{1} = YLabs;
    param_plots.Tables{1} = rel_tab;
    param_plots.Legend{1} = rel_leg;
    param_plots(1,{'XTick','XTickLab','XLabel','XScale'}) = sub_tab(1,{'XTick','XTickLab','XLabel','XScale'});
    all_plot_tabs.CycAvg{'HIT'} = cycavg_plots;
    all_plot_tabs.MaxVel{'HIT'} = maxvel_plots;
    all_plot_tabs.Param{'HIT'} = param_plots;
end
%% Velocity Step
if any(contains(all_results.Type,'Exponential'))
    %Isolate relevant table entries and put them in order for plotting by amplitude and frequency
    tab = sortrows(all_results(contains(all_results.Type,'Exponential'),:),'Amplitude','ascend');
    tab.AmpStr = strcat(strrep(cellstr(num2str(tab.Amplitude)),' ',''),'dps');
    fn = size(tab,1);
    YMax = 10*ceil(max(tab.MaxVel)/10);
    %Make one figure for each group of cycle averages (same subject, visit,
    %date, goggle, axis, and amplitude across condition)
    rel_labs = {'Subject','Visit','DateStr','Condition','Goggle','AmpStr'};
    [~,rel_col] = ismember(rel_labs,tab.Properties.VariableNames);
    file_parts = table2cell(tab(:,rel_col));
    %Parts applicable to each file
    common_cond_i = all(strcmp(file_parts,repmat(file_parts(1,:),fn,1)),1);
    common_cond = strjoin(file_parts(1,common_cond_i));
    common_cond_i(contains(rel_labs,'Condition')) = 0;
    rel_file_parts = join(file_parts(:,~common_cond_i),2);
    [conds,~,IC] = unique(rel_file_parts,'stable');
    [~,Ic] = ismember(tab.AxisName,all_canals(:,1));
    [~,~,IFA] = unique([IC,Ic],'rows','stable');
    dat = array2table(reshape(unique(Ic),[],1),'VariableNames',{'Canal_i'});
    for i = 1:size(dat,1)
        inds = find(Ic==dat.Canal_i(i)|0*Ic==dat.Canal_i(i));
        %Remove duplicate experiments by choosing the one with the
        %closest time stamp to the others
        [gc,gr] = groupcounts(IFA(inds));
        dup = gr(gc>1);
        for ii = 1:length(dup)
            sub_i = find(IFA(inds)==dup(ii));
            [~,k_i] = min(abs(tab.Date(inds(sub_i))-median(tab.Date(inds))));
            sub_i(k_i) = []; %Keep this one
            inds(sub_i) = []; %Delete the others
        end
        dat.Name{i} = [common_cond,' ',all_canals{dat.Canal_i(i),1},' Exponential'];
        dat.Files{i} = tab.File(inds);
        dat.SubNames{i} = conds;
        dat.YLim{i} = YMax*[-1 1];
    end
    cycavg_plots = dat(:,2:end); %Figure name and cyc avg file names for all
    maxvel_plots = [];
    param_plots = [];
    all_plot_tabs.CycAvg{'Exp'} = cycavg_plots;
    all_plot_tabs.MaxVel{'Exp'} = maxvel_plots;
    all_plot_tabs.Param{'Exp'} = param_plots;
end
%% Autoscan - STILL NEED TO FIX
if any(contains(all_results.Condition,'Autoscan'))
    %     %Sorted to already be in the correct order
    %     all_results = all_results(contains(all_results.Condition,'Autoscan'),:);
    %     temp_e = split(all_results.Electrode,'E');
    %     all_results.Enum = str2double(temp_e(:,2));
    %     all_results = sortrows(sortrows(sortrows(all_results,'CurrentAmp(uA)','ascend'),'PhaseDur(us)','ascend'),'Enum','ascend');
    %     elec_phase = strcat(all_results.Electrode,{' '},num2str(all_results.('PhaseDur(us)')),'us');
    %     elec_opts = unique(elec_phase,'stable');
    %     %Select Files to Plot
    %     canal = listdlg('PromptString','Plot which electrodes and phase durations?',...
    %         'ListString',elec_opts,'SelectionMode','multiple');
    %     if isempty(canal)
    %         return;
    %     end
    %     rel_results = all_results(ismember(elec_phase,elec_opts(canal)),:);
    %     % Make the figure name
    %     file_parts = [rel_results.Subject,rel_results.Visit,cellstr(datestr(rel_results.Date,'yyyymmdd')),...
    %         rel_results.Condition,rel_results.Goggle,strcat(strrep(cellstr(num2str(rel_results.('PulseFreq(pps)'))),' ',''),'pps')];
    %     fig_title = {strjoin([join(unique(file_parts,'stable')),strjoin(elec_opts(canal),'_')])};
    %     elec_phase = strcat(rel_results.Electrode,{' '},num2str(rel_results.('PhaseDur(us)')),'us');
    %     rel_ep = unique(elec_phase,'stable');
    %     n_row = length(rel_ep);
    %     rel_ep_inds = false(size(rel_results,1),n_row);
    %     rel_ep_s  = NaN(1,n_row);
    %     for i = 1:n_row
    %         rel_ep_inds(:,i) = contains(elec_phase,rel_ep{i});
    %         rel_ep_s(i) = find(rel_ep_inds(:,i)==1,1,'first');
    %     end
    %     n_col = sum(rel_ep_inds);
    %     f_order = rel_results.File;
    %     curr_lab = cellstr(num2str(rel_results.('CurrentAmp(uA)')));
    %     curr_lab(rel_ep_s) = strcat(curr_lab(rel_ep_s),'\muA');
    %     cyc_lab = cellstr(num2str(rel_results.Cycles));
    %     cyc_lab(rel_ep_s) = strcat('n=',cyc_lab(rel_ep_s));
    %     % Plot Current Levels
    %     fig = figure;
    %     fig.Color = [1,1,1];
    %     fig.Units = 'inches';
    %     fig.Position = [1 1 8 4];
    %     if annot
    %         annotation('textbox',[0 0 1 1],'String',['VOGA',version,...
    %             newline,Experimenter],'FontSize',5,...
    %             'EdgeColor','none','interpreter','none');
    %         annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_name],'FontSize',5,...
    %             'EdgeColor','none','interpreter','none','VerticalAlignment','bottom');
    %     end
    %     ha = gobjects(1,length(f_order));
    %     %Set params
    %     grid_on = false;
    %     if isempty(YMax)
    %         YMax = 100;
    %     end
    %     YLim = YMax*[-1 1];
    %     x_min = 0.15;
    %     x_max = 0.95;
    %     space_x = 0.01;
    %     y_min = 0.08;
    %     y_max = 0.98;
    %     space_y = 0.03;
    %     %Calculate
    %     x_wid = (x_max - x_min - space_x*(n_col-1))./n_col;
    %     y_wid = (y_max - y_min - space_y*(n_row-1))/n_row;
    %     y_pos = fliplr(y_min:(y_wid+space_y):y_max);
    %     fig_x_wid = NaN(1,length(f_order));
    %     fig_row_pos = NaN(1,length(f_order));
    %     fig_col_pos = NaN(1,length(f_order));
    %     k = 0;
    %     for i = 1:n_row
    %         x_pos = x_min:(x_wid(i)+space_x):x_max;
    %         fig_row_pos((1:n_col(i))+k) = x_pos;
    %         fig_col_pos((1:n_col(i))+k) = y_pos(i);
    %         fig_x_wid((1:n_col(i))+k) = x_wid(i);
    %         k = k+n_col(i);
    %     end
    %     for i = 1:length(f_order)
    %         ha(i) = subplot(n_row,max(n_col),i);
    %         %Load and plot
    %         b = load([Cyc_Path,filesep,f_order{i}]);
    %         a = fieldnames(b);
    %         CycAvg = b.(a{1});
    %         fields = fieldnames(CycAvg);
    %         if ~ismember('t',fields)
    %             CycAvg.t = (0:1/CycAvg.Fs:(length(CycAvg.ll_cycavg)-1)/CycAvg.Fs)';
    %         end
    %         if length(CycAvg.t) > 1000
    %             s = round(linspace(1,length(CycAvg.t),1000));
    %         else
    %             s = 1:length(CycAvg.t);
    %         end
    %         hold on
    %         %Now add the fills and standard deviations and means
    %         trac = {'l_z','r_z','l_l','r_l','l_r','r_r'};
    %         if contains(f_order{i},{'RA','LP'})%RALP
    %             curr_col = colors.l_r;
    %             rel_trac = {'l_r','r_r'};
    %         elseif contains(f_order{i},{'LH','RH'})%LHRH
    %             curr_col = colors.l_z;
    %             rel_trac = {'l_z','r_z'};
    %         elseif contains(f_order{i},{'LA','RP'})%LARP
    %             curr_col = colors.l_l;
    %             rel_trac = {'l_l','r_l'};
    %         end
    %         for j = 1:length(trac)
    %             trace = strrep(trac{j},'_','');
    %             plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s),'Color',colors.(trac{j}))
    %             plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s),'Color',colors.(trac{j}))
    %             plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s),'Color',colors.(trac{j}),'LineWidth',2);
    %         end
    %         %Plot the intended canal again so that it's in the foreground
    %         for j = 1:length(rel_trac)
    %             trace = strrep(rel_trac{j},'_','');
    %             plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s),'Color',colors.(rel_trac{j}))
    %             plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s),'Color',colors.(rel_trac{j}))
    %             plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s),'Color',colors.(rel_trac{j}),'LineWidth',2);
    %         end
    %         hold off
    %         axis([0 0.5 YLim])
    %         text(0.5,YLim(2)-0.01*diff(YLim),curr_lab{i},'Color',curr_col,'HorizontalAlignment','right','VerticalAlignment','top')
    %         text(0.5,YLim(1)+0.01*diff(YLim),cyc_lab{i},'Color','k','HorizontalAlignment','right','VerticalAlignment','bottom')
    %     end
    %     for i = 1:length(f_order)
    %         set(ha(i),'Position',[fig_row_pos(i),fig_col_pos(i),fig_x_wid(i),y_wid])
    %         if ~ismember(i,[rel_ep_s-1,length(f_order)])
    %             annotation('line',(fig_row_pos(i)+fig_x_wid(i)+0.5*space_x)*[1 1],fig_col_pos(i)+[0 y_wid],'LineWidth',1,'LineStyle','--')
    %         end
    %     end
    %     annotation('line',fig_row_pos(end)+[0 x_wid(end)],y_min-0.01*[1 1],'LineWidth',2)
    %     annotation('line',(x_max+space_x)*[1 1],y_min+[0 YMax*y_wid/(2*YLim(2))],'LineWidth',2)
    %     annotation('textbox','String','0.5s','EdgeColor','none',...
    %         'Position',[fig_row_pos(end),0,x_wid(end),y_min-0.01],'HorizontalAlignment','right','VerticalAlignment','middle')
    %     annotation('textbox','String',[num2str(YMax),newline,'\circ/s'],'EdgeColor','none',...
    %         'Position',[x_max+space_x,y_min,1-(x_max+space_x),YMax*y_wid/(2*YLim(2))],'HorizontalAlignment','center','VerticalAlignment','middle')
    %     set(ha,'XColor',[1,1,1],'YColor',[1,1,1])
    %     lab_ax2 = axes('Position',[x_min y_min x_min y_max-y_min],'XColor',[1 1 1],'YColor',[1,1,1],...
    %      'Color','none','XTick',[],'YTick',[]);
    %     ylabel(lab_ax2,unique(rel_results.Subject),'FontSize',30,'FontWeight','bold','Color','k','Position',[-0.6 0.5 0])
    %     for i = 1:length(rel_ep_s)
    %         if contains(elec_phase(i),{'RA','LP'}) %RALP
    %             curr_col = colors.l_r;
    %         elseif contains(elec_phase(i),{'LH','RH'}) %LHRH
    %             curr_col = colors.l_z;
    %         elseif contains(elec_phase(i),{'RP','LA'}) %LARP
    %             curr_col = colors.l_l;
    %         else
    %             curr_col = 'k';
    %         end
    %         ylabel(ha(rel_ep_s(i)),split(elec_phase(rel_ep_s(i)),' '),'Color',curr_col,'FontSize',20,'Position',[-0.05 0 -1]);
    %     end
    %     if(grid_on)
    %         set(ha,'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
    %     end
    %     savefig([Path,filesep,'Figures',filesep,strrep(fig_title{:},' ','-'),'.fig'])
    %
    %     %Make Figure like Boutros 2019 Figure 4 but Magnitude and Misalignment
    %     tab = all_results(contains(all_results.Condition,'Autoscan'),:);
    %     tab.Enum = str2double(extract(tab.Electrode,digitsPattern(2)|digitsPattern(1)));
    %     %Sort to be in phase dur, pulse freq, electrode, and then current order
    %     tab = sortrows(sortrows(sortrows(sortrows(tab,'CurrentAmp','ascend'),'Enum','ascend'),'PulseFreq','ascend'),'PhaseDur','ascend');
    %     fn = size(tab,1);
    %     file_parts = [tab.Subject,cellstr(datestr(tab.Date,'yyyymmdd')),...
    %         tab.Goggle,tab.Electrode,strcat(strrep(cellstr(num2str(tab.PhaseDur)),' ',''),'us'),...
    %         strcat(strrep(cellstr(num2str(tab.PulseFreq)),' ',''),'pps'),...
    %         strcat(strrep(cellstr(num2str(tab.CurrentAmp)),' ',''),'uA')];
    %     common_cond_i = all(strcmp(file_parts(:,1:end-1),repmat(file_parts(1,1:end-1),fn,1)),1); %everything but amp
    %     rel_file_parts = file_parts(:,[~common_cond_i,false]);
    %     common_cond = strjoin(file_parts(1,[common_cond_i,false]));
    %     if isempty(rel_file_parts) %only one condition
    %         conds = common_cond;
    %         IC = ones(fn,1);
    %     else
    %         [all_conds,~,IC] = unique(join(rel_file_parts,2),'stable');
    %         %Let the user pick which items to plot
    %         [indx,tf] = listdlg('ListString',all_conds,...
    %             'PromptString',['Pick the experiments to plot from ',common_cond],'ListSize',[400 300],'InitialValue',1:length(all_conds));
    %         if ~tf
    %             return;
    %         end
    %         rel_file_parts = rel_file_parts(ismember(IC,indx),:);
    %         tab = tab(ismember(IC,indx),:);
    %         fn = size(tab,1);
    %         [conds,ic,IC] = unique(join(rel_file_parts,2),'stable');
    %     end
    %     [~,can_col_i] = ismember(tab.AxisName,all_canals(:,1));
    %     can_col = all_canals(can_col_i,2);
    %     maxvel = NaN(2,fn);
    %     align = NaN(2,fn);
    %     for i = 1:fn
    %         [maxvel(1,i),eye] = max(abs([tab.(['MaxVel_L',can_col{i}])(i),tab.(['MaxVel_R',can_col{i}])(i)]));
    %         maxvel(2,i) = tab.(['MaxVel_',eyes{eye},can_col{i},'_sd'])(i);
    %         align(:,i) = [tab.(['Align_',eyes{eye}])(i);tab.(['Align_',eyes{eye},'_sd'])(i)];
    %     end
    %     %Make some bold (if you know which one was activated on)
    %     E_num = strcat('E',cellfun(@num2str,num2cell(3:11),'UniformOutput',false));
    %     indx = listdlg('ListString',E_num,...
    %         'PromptString','Pick the electrodes to bold. Press Cancel for none.','ListSize',[400 300],'SelectionMode','multiple');
    %     E_bold = E_num(indx);
    %     %Make the figure
    %     ha = gobjects(1,6);
    %     errorbarcapsize=1;
    %     figsizeinches=[7,6];
    %     XLim = [-5 105];
    %     YLim_align = [0 80];
    %     figure('Units','inch','Position',[2 2 figsizeinches],'Color',[1,1,1]);
    %     if annot
    %         annotation('textbox',[0 0 1 1],'String',[Path,newline,code_name,newline,...
    %             'VOGA',version,newline,Experimenter],'FontSize',5,...
    %             'EdgeColor','none','interpreter','none');
    %     end
    %     annotation('textbox',[0 .9 1 .1],'String',strrep(common_cond,'us','\mus'),'FontSize',14,...
    %         'HorizontalAlignment','center','EdgeColor','none');
    %     color = NaN(length(conds),3);
    %     for i = 1:length(conds)
    %         if contains(conds{i},E_bold)
    %             color(i,:) = colors.(['l_',lower(can_col{ic(i)})]);
    %         else
    %             color(i,:) = colors.(['l_',lower(can_col{ic(i)}),'_s']);
    %         end
    %     end
    %     markers = {'x','o','d','s','v','p','^','+','<'};
    %     % Posterior Plots
    %     rel_cond = find(contains(conds,{'E3','E4','E5'}));
    %     rel_pnum = [1,4];
    %     ha(rel_pnum(1)) = subplot(2,3,rel_pnum(1));
    %     hold on
    %     plot(NaN,NaN)
    %     h3 = gobjects(length(rel_cond),1);
    %     for i = 1:length(rel_cond)
    %         curr_ax = tab.CurrentAmp(IC==rel_cond(i));
    %         curr_norm = 100*(curr_ax - curr_ax(1))/(curr_ax(end)-curr_ax(1));
    %         h3(i) = errorbar(curr_norm,maxvel(1,IC==rel_cond(i)),maxvel(2,IC==rel_cond(i)),'Marker',markers{i},'Color',color(rel_cond(i),:),'LineStyle','-','LineWidth',1,'CapSize',errorbarcapsize);
    %     end
    %     hold off
    %     if ~isempty(rel_cond)
    %         leg = legend(ha(rel_pnum(1)),h3,conds(rel_cond),'location','northwest','box','off','FontSize',7);
    %         leg.ItemTokenSize(1) = 12;
    %     end
    %     ha(rel_pnum(2)) = subplot(2,3,rel_pnum(2));
    %     hold on
    %     plot(NaN,NaN)
    %     for i = 1:length(rel_cond)
    %         curr_ax = tab.CurrentAmp(IC==rel_cond(i));
    %         curr_norm = 100*(curr_ax - curr_ax(1))/(curr_ax(end)-curr_ax(1));
    %         errorbar(curr_norm,align(1,IC==rel_cond(i)),align(2,IC==rel_cond(i)),'Marker',markers{i},'Color',color(rel_cond(i),:),'LineStyle','-','LineWidth',1,'CapSize',errorbarcapsize)
    %     end
    %     hold off
    %     % Horizontal Plots
    %     rel_cond = find(contains(conds,{'E6','E7','E8'}));
    %     rel_pnum = [2,5];
    %     ha(rel_pnum(1)) = subplot(2,3,rel_pnum(1));
    %     hold on
    %     plot(NaN,NaN)
    %     h3 = gobjects(length(rel_cond),1);
    %     for i = 1:length(rel_cond)
    %         curr_ax = tab.CurrentAmp(IC==rel_cond(i));
    %         curr_norm = 100*(curr_ax - curr_ax(1))/(curr_ax(end)-curr_ax(1));
    %         h3(i) = errorbar(curr_norm,maxvel(1,IC==rel_cond(i)),maxvel(2,IC==rel_cond(i)),'Marker',markers{i},'Color',color(rel_cond(i),:),'LineStyle','-','LineWidth',1,'CapSize',errorbarcapsize);
    %     end
    %     hold off
    %     if ~isempty(rel_cond)
    %         leg = legend(ha(rel_pnum(1)),h3,conds(rel_cond),'location','northwest','box','off','FontSize',7);
    %         leg.ItemTokenSize(1) = 12;
    %     end
    %     ha(rel_pnum(2)) = subplot(2,3,rel_pnum(2));
    %     hold on
    %     plot(NaN,NaN)
    %     for i = 1:length(rel_cond)
    %         curr_ax = tab.CurrentAmp(IC==rel_cond(i));
    %         curr_norm = 100*(curr_ax - curr_ax(1))/(curr_ax(end)-curr_ax(1));
    %         errorbar(curr_norm,align(1,IC==rel_cond(i)),align(2,IC==rel_cond(i)),'Marker',markers{i},'Color',color(rel_cond(i),:),'LineStyle','-','LineWidth',1,'CapSize',errorbarcapsize)
    %     end
    %     hold off
    %     % Anterior Plots
    %     rel_cond = find(contains(conds,{'E9','E10','E11'}));
    %     rel_pnum = [3,6];
    %     ha(rel_pnum(1)) = subplot(2,3,rel_pnum(1));
    %     hold on
    %     plot(NaN,NaN)
    %     h3 = gobjects(length(rel_cond),1);
    %     for i = 1:length(rel_cond)
    %         curr_ax = tab.CurrentAmp(IC==rel_cond(i));
    %         curr_norm = 100*(curr_ax - curr_ax(1))/(curr_ax(end)-curr_ax(1));
    %         h3(i) = errorbar(curr_norm,maxvel(1,IC==rel_cond(i)),maxvel(2,IC==rel_cond(i)),'Marker',markers{i},'Color',color(rel_cond(i),:),'LineStyle','-','LineWidth',1,'CapSize',errorbarcapsize);
    %     end
    %     hold off
    %     if ~isempty(rel_cond)
    %         leg = legend(ha(rel_pnum(1)),h3,conds(rel_cond),'location','northwest','box','off','FontSize',7);
    %         leg.ItemTokenSize(1) = 12;
    %     end
    %     ha(rel_pnum(2)) = subplot(2,3,rel_pnum(2));
    %     hold on
    %     plot(NaN,NaN)
    %     for i = 1:length(rel_cond)
    %         curr_ax = tab.CurrentAmp(IC==rel_cond(i));
    %         curr_norm = 100*(curr_ax - curr_ax(1))/(curr_ax(end)-curr_ax(1));
    %         errorbar(curr_norm,align(1,IC==rel_cond(i)),align(2,IC==rel_cond(i)),'Marker',markers{i},'Color',color(rel_cond(i),:),'LineStyle','-','LineWidth',1,'CapSize',errorbarcapsize)
    %     end
    %     hold off
    %     x_min = 0.1;
    %     x_max = 0.99;
    %     x_spac = 0.01;
    %     y_min = 0.08;
    %     y_max = 0.93;
    %     y_space = 0.03;
    %     x_wid = (x_max-x_min-2*x_spac)/3;
    %     y_wid = (y_max-y_min-y_space)/2;
    %     x_pos = x_min:x_wid+x_spac:x_max;
    %     y_pos = y_min:y_wid+y_space:y_max;
    %     ha(1).Position = [x_pos(1),y_pos(2),x_wid,y_wid];
    %     ha(2).Position = [x_pos(2),y_pos(2),x_wid,y_wid];
    %     ha(3).Position = [x_pos(3),y_pos(2),x_wid,y_wid];
    %     ha(4).Position = [x_pos(1),y_pos(1),x_wid,y_wid];
    %     ha(5).Position = [x_pos(2),y_pos(1),x_wid,y_wid];
    %     ha(6).Position = [x_pos(3),y_pos(1),x_wid,y_wid];
    %     if isempty(YMax)
    %         YMax = max(reshape(cell2mat(get(ha(1:3),'Ylim')),[],1));
    %     end
    %     set(ha,'box','off','XLim',XLim)
    %     ylabel(ha(1),{'Eye Velocity';'Magnitude (\circ/s)'})
    %     ylabel(ha(4),{'Misalignment';'Angle (\circ)'})
    %     xlabel(ha(5),'% of Current Range')
    %     set(ha(1:3),'YLim',[0,YMax],'XColor','none')
    %     set(ha(4:6),'YLim',YLim_align,'XTick',0:20:100)
    %     set(ha(1),'YTick',20:20:YMax)
    %     set(ha(4),'YTick',10:10:max(YLim_align))
    %     set(ha([2,3,5,6]),'YColor','none')
    %     fig_name = inputdlg('Name this figure','',1,{[strrep(common_cond,' ','-'),'.fig']});
    %     if ~isempty(fig_name)
    %         savefig([Path,filesep,'Figures',filesep,fig_name{:}])
    %     end
    %     close;
    
    all_plot_tabs.CycAvg{'PulseTrain'} = cycavg_plots;
    all_plot_tabs.MaxVel{'PulseTrain'} = maxvel_plots;
    all_plot_tabs.Param{'PulseTrain'} = param_plots;
end
%% Make all plots
all_cycavg_plots = vertcat(all_plot_tabs.CycAvg{:});
all_maxvel_plots = vertcat(all_plot_tabs.MaxVel{:});
all_param_plots = vertcat(all_plot_tabs.Param{:});
% Group Cyc Avg
if isfolder(Cyc_Path) %Otherwise, skip these and go to the parameterized versions
    for k = 1:size(all_cycavg_plots,1)
        plot_info = all_cycavg_plots(k,:);
        fig = plotGroupCycAvg(plot_info,params);
        savefig(fig,[Path,filesep,'Figures',filesep,strrep(plot_info.Name{1},' ','-'),'.fig'])
        saveas(fig,[Path,filesep,'CRFs',filesep,strrep(plot_info.Name{1},' ','-'),'.png'])
        close;
    end
end
% MaxVel Parameterized
for k = 1:size(all_maxvel_plots,1)
    plot_info = all_maxvel_plots(k,:);
    fig = plotMaxVelAllEyeComp(plot_info,params);
    savefig(fig,[Path,filesep,'Figures',filesep,strrep(plot_info.Name{1},' ','-'),'.fig'])
    saveas(fig,[Path,filesep,'CRFs',filesep,strrep(plot_info.Name{1},' ','-'),'.png'])
    close;
end
% Parameterized Figures
for k = 1:size(all_param_plots,1)
    plot_info = all_param_plots(k,:);
    fig = plotTabParam(plot_info,params);
    savefig(fig,[Path,filesep,'Figures',filesep,strrep(plot_info.Name{1},' ','-'),'.fig'])
    saveas(fig,[Path,filesep,'CRFs',filesep,strrep(plot_info.Name{1},' ','-'),'.png'])
    close;
end
end