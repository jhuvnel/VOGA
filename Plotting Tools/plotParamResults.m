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
    YMax = 5*ceil(max(reshape(tab{:,{'MaxVel_LL','MaxVel_LR','MaxVel_LZ','MaxVel_RL','MaxVel_RR','MaxVel_RZ'}}+...
        tab{:,{'MaxVel_LL_sd','MaxVel_LR_sd','MaxVel_LZ_sd','MaxVel_RL_sd','MaxVel_RR_sd','MaxVel_RZ_sd'}},[],1))/5);
    GainMax = 0.1*ceil(max(tab.Gain+tab.Gain_sd)/0.1);
    %Make one figure for each group of cycle averages (same subject, visit,
    %date, condition, goggle, axis, across frequency, and/or amplitude)
    %Make one figure per for frequency and/or amplitude sweep with MaxVel
    %of each each component for each condition and one figure with Gain and
    %Phase
    YVar = {'Gain';'Phase'};
    YLabs = {'Gain';'Phase Lead (deg)'};
    rel_labs = {'Subject','Visit','DateStr','Experiment','Type','Condition','Goggle'};
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
                        ' ',all_canals{dat.Canal_i(i),2}];
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
                    amp = [str_amps{dat.Amp_i(i)}];
                    name = [common_cond,' ',conds{dat.Cond_i(i)},...
                        ' ',all_canals{dat.Canal_i(i),2},' ',...
                        str_amps{dat.Amp_i(i)}];
                    x_val = IF(inds);
                    x_var = 'Frequency';
                    x_tick = freqs;
                    x_tlab = strrep(str_freqs,'Hz','');
                    x_lab = 'Frequency (Hz)';
                    x_scale = 'log';
                    sub_names = str_freqs;
                    files = cell(fnum,1);
                elseif dat.Amp_i(i)==0 %Freq
                    freq = [str_freqs{dat.Freq_i(i)}];
                    amp = '';
                    name = [common_cond,' ',conds{dat.Cond_i(i)},...
                        ' ',all_canals{dat.Canal_i(i),2},' ',...
                        str_freqs{dat.Freq_i(i)}];
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
    
    fn = size(tab,1);
    GainMax = 0.1*ceil(max(tab.Gain+tab.Gain_sd)/0.1);
    %Make one figure for each group of cycle averages (same subject, visit,
    %date, goggle, canal across conditions
    %Make one figure with Gain and Lantecy across canals for each condition
    YVar = {'Gain';'Latency'};
    YLabs = {'Gain';'Latency (ms)'};
    rel_labs = {'Subject','Visit','DateStr','Experiment','Type','Condition','Goggle'};
    [~,rel_col] = ismember(rel_labs,tab.Properties.VariableNames);
    file_parts = table2cell(tab(:,rel_col));
    %Parts applicable to each file
    common_cond_i = all(strcmp(file_parts,repmat(file_parts(1,:),fn,1)),1);
    common_cond = strjoin(file_parts(1,common_cond_i));
    common_cond_i(contains(rel_labs,'Condition')) = 0;
    rel_file_parts = join(file_parts(:,~common_cond_i),2);
    [conds,~,IC] = unique(rel_file_parts,'stable');
    enum = length(conds);
    tab.IC = IC;
    tab = sortrows(sortrows(tab,'IC','ascend'),'Ic','ascend');
    IC = tab.IC;
    Ic = tab.Ic;
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
            name = [common_cond,' ',canal];
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
    param_plots = cell2table({[common_cond,' Gain']},'VariableNames',{'Name'});
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
    tab = sortrows(all_results(contains(all_results.Type,'Exponential')&~strcmp(all_results.AxisName,''),:),'Amplitude','ascend');
    tab.AmpStr = strcat(strrep(cellstr(num2str(tab.Amplitude)),' ',''),'dps');
    fn = size(tab,1);
    YMax = 10*ceil(max(tab.MaxVel)/10);
    %Make one figure for each group of cycle averages (same subject, visit,
    %date, goggle, axis, and amplitude across condition)
    rel_labs = {'Subject','Visit','DateStr','Experiment','Type','Condition','Goggle','AmpStr'};
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
        dat.Name{i} = [common_cond,' ',all_canals{dat.Canal_i(i),1}];
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
%% Autoscan
if any(contains(all_results.Condition,'Autoscan'))
    all_elec = strrep(join([reshape(repmat({'LPE','LHE','LAE','RPE','RHE','RAE'},...
        3,1),[],1),repmat(cellstr(num2str((3:11)')),2,1)]),' ','');
    %Isolate relevant table entries and put them in order for plotting by
    %pulse rate and current amplitude
    temp_e = split(all_results.Electrode,'E');
    all_results.Enum = str2double(temp_e(:,2));
    tab = sortrows(sortrows(sortrows(sortrows(all_results(contains(all_results.Condition,'Autoscan')&...
        contains(all_results.AxisName,all_canals(:,1)),:),'CurrentAmp','ascend'),...
        'PulseFreq','ascend'),'Enum','ascend'),'PhaseDur','ascend');
    tab.CurrentAmpStr = strcat(strrep(cellstr(num2str(tab.CurrentAmp)),' ',''),'uA');
    tab.PhaseDurStr = strcat(strrep(cellstr(num2str(tab.PhaseDur)),' ',''),'us');
    tab.PulseFreqStr = strcat(strrep(cellstr(num2str(tab.PulseFreq)),' ',''),'pps');
    tab.CurrentAmpNorm = tab.CurrentAmp; %Update in the for loop based on the current range available
    tab.CyclesStr = strcat('n=',strrep(cellstr(num2str(tab.Cycles)),' ',''));
    fn = size(tab,1);
    YMax = 5*ceil(max(tab.MaxVel+tab.MaxVel_sd)/5);
    W = 0.75;
    tab.Score = (W*tab.MaxVel./max(tab.MaxVel) + (1-W)*(1-tab.Align/90))*100;
    tab.Score_sd = (W*tab.MaxVel_sd./max(tab.MaxVel) + (1-W)*(tab.Align_sd/90))*100;
    %Make one figure for each group of cycle averages (same subject, visit,
    %date, condition, goggle, axis, across frequency, and/or amplitude)
    %Make one figure per for frequency and/or amplitude sweep with MaxVel
    %of each each component for each condition and one figure with Gain and
    %Phase
    YVar = {'MaxVel';'Align'};
    YLabs = {'Eye Velocity (dps)';'Misalignment (deg)';'Score (%)'};
    rel_labs = {'Subject','Visit','DateStr','Condition','Experiment','Type','Goggle',...
        'PulseFreqStr','PhaseDurStr','Electrode'};
    [~,rel_col] = ismember(rel_labs,tab.Properties.VariableNames);
    file_parts = strrep(table2cell(tab(:,rel_col)),'us','\mus');
    %Parts applicable to each file
    common_cond_i = all(strcmp(file_parts,repmat(file_parts(1,:),fn,1)),1);
    common_cond = strjoin(file_parts(1,common_cond_i));
    if all(common_cond_i)
        common_cond_i(contains(rel_labs,'Condition')) = 0;
    end
    [conds,ic,IC] = unique(join(file_parts(:,~common_cond_i),2),'stable');
    enum = length(conds);
    dat = cell2table(conds,'VariableNames',{'Cond'});
    dat.PulseFreq = tab.PulseFreqStr(ic);
    dat.PhaseDur = tab.PhaseDurStr(ic);
    dat.Electrode = tab.Electrode(ic);
    dat.Canal = tab.AxisName(ic);
    dat.Enum = tab.Enum(ic);    
    for i = 1:enum
    	inds = find(IC==i);
        tab.CurrentAmpNorm(inds) = round(100*(tab.CurrentAmp(inds)-tab.CurrentAmp(inds(1)))/diff(tab.CurrentAmp(inds([1,end]))));
        %Remove duplicate experiments by choosing the one with the
        %closest time stamp to the others
        [gc,gr] = groupcounts(tab.CurrentAmp(inds));
        dup = gr(gc>1);
        for ii = 1:length(dup)
            sub_i = find(tab.CurrentAmp(inds)==dup(ii));
            [~,k_i] = min(abs(tab.Date(inds(sub_i))-median(tab.Date(inds))));
            sub_i(k_i) = []; %Keep this one
            inds(sub_i) = []; %Delete the others
        end
        dat.Cond_noE{i} = strrep(dat.Cond{i},dat.Electrode{i},'');
        dat.Name{i} = [common_cond,' ',dat.Cond{i}];
        dat.Files{i} = tab.File(inds);
        dat.Tab{i} = tab(inds,:);
        dat.SubNames{i} = join([strrep(tab.CurrentAmpStr(inds),'uA','\muA'),tab.CyclesStr(inds)],newline)';
        dat.YLim{i} = YMax*[-1 1];
        dat.XVar{i} = 'CurrentAmpNorm';
        dat.XTick{i} = 0:25:100;
        dat.XTickLab{i} = dat.XTick{i};
        dat.XLabel{i} = '% of Current Range';
        dat.XScale{i} = 'linear';
        dat.X{i} = tab.CurrentAmpNorm(inds);
        dat.CurrentAmp{i} = tab.CurrentAmp(inds);
        for ii = 1:length(YVar)
            dat.(YVar{ii}){i} = tab.(YVar{ii})(inds);
            dat.([YVar{ii},'_sd']){i} = tab.([YVar{ii},'_sd'])(inds);
        end
        dat.Score{i} = tab.Score(inds);
        dat.MaxScore{i} = max(tab.Score(inds));
    end
    cycavg_plots = dat(:,{'Name','Files','SubNames','YLim'}); %Figure name and cyc avg file names for all
    %Initialize the table for MaxVel plots
    [mv_name,mp_i2,maxvel_plot_i] = unique(dat.Cond_noE,'stable');
    maxvel_plots = table();
    maxvel_plots.Name = strrep(strcat([common_cond,' '],mv_name,' Velocity'),'  ',' ');
    for i=1:size(maxvel_plots,1)
        sub_tab = dat(maxvel_plot_i==i,:);
        rel_tab = cell(3,3);
        rel_tab(sub_tab.Enum-2) = sub_tab.Tab;
        maxvel_plots.SubNames{i} = reshape(all_elec(contains(all_elec,sub_tab.Canal{1}(1))),3,3);
        maxvel_plots.Tables{i} = rel_tab;
        maxvel_plots.YLim{i} = [0 YMax];
    end
    maxvel_plots(:,{'XVar','XTick','XTickLab','XLabel','XScale'}) = dat(mp_i2,{'XVar','XTick','XTickLab','XLabel','XScale'});
    %Intialized legend as needed
    leg_tab = cell2table(conds,'VariableNames',{'Name'});
    leg_tab.Marker(:) = {'none'};
    leg_tab.LineStyle(:) = {'none'};
    leg_tab.Color(:) = {[0,0,0]};
    leg_tab.Marker(contains(dat.Canal,'P')) = all_markers(1:sum(contains(dat.Canal,'P')));
    leg_tab.Marker(contains(dat.Canal,'A')) = all_markers(1:sum(contains(dat.Canal,'A')));
    leg_tab.Marker(contains(dat.Canal,'H')) = all_markers(1:sum(contains(dat.Canal,'H')));
    leg_tab.Color(contains(dat.Canal,{'LP','RA'})) = {colors.l_r};
    leg_tab.Color(contains(dat.Canal,{'LA','RP'})) = {colors.l_l};
    leg_tab.Color(contains(dat.Canal,{'H'})) = {colors.l_z};
    %Initialize the table for TabParam plots
    param_plots = table();
    param_plots.Name{1} = common_cond;
    sub_tab = dat;
    sub_t = strcat(sub_tab.Canal{1}(1),{'P','H','A'});
    [~,indC] = ismember(sub_tab.Canal,sub_t);
    rel_tab = cell(length(YVar),length(sub_t));
    rel_leg = cell(length(YVar),length(sub_t));
    for ii = 1:length(sub_t) 
        rel_leg{1,ii} = leg_tab(indC==ii,:);
        for j = 1:length(YVar)
            n_tab = table();
            n_tab.Marker = leg_tab.Marker(indC==ii);
            n_tab.LineStyle(:) = {'-'};
            n_tab.Color = leg_tab.Color(indC==ii);
            n_tab.X = sub_tab.X(indC==ii);
            n_tab.Y = sub_tab.(YVar{j})(indC==ii);
            n_tab.Y_sd = sub_tab.([YVar{j},'_sd'])(indC==ii);
            rel_tab{j,ii} = n_tab;
        end        
    end
    param_plots.SubNames{1} = sub_t;
    param_plots.YLim{1} = [0,YMax;0,180;-10,110];
    param_plots.YLabel{1} = YLabs;
    param_plots.Tables{1} = rel_tab;
    param_plots.Legend{1} = rel_leg;
    param_plots(:,{'XTick','XTickLab','XLabel','XScale'}) = dat(1,{'XTick','XTickLab','XLabel','XScale'});
    all_plot_tabs.CycAvg{'PulseTrain'} = cycavg_plots;
    all_plot_tabs.MaxVel{'PulseTrain'} = maxvel_plots;
    all_plot_tabs.Param{'PulseTrain'} = param_plots;
    %Save data for later
    save('AutoscanParameters.mat','dat')
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
        fig_name = strrep(strrep(plot_info.Name{1},' ','-'),'\mu','u');
        savefig(fig,[Path,filesep,'Figures',filesep,fig_name,'.fig'])
        saveas(fig,[Path,filesep,'CRFs',filesep,fig_name,'.png'])
        close;
    end
end
% MaxVel Parameterized
for k = 1:size(all_maxvel_plots,1)
    plot_info = all_maxvel_plots(k,:);
    fig = plotMaxVelAllEyeComp(plot_info,params);
    fig_name = strrep(strrep(plot_info.Name{1},' ','-'),'\mu','u');
    savefig(fig,[Path,filesep,'Figures',filesep,fig_name,'.fig'])
    saveas(fig,[Path,filesep,'CRFs',filesep,fig_name,'.png'])
    close;
end
% Parameterized Figures
for k = 1:size(all_param_plots,1)
    plot_info = all_param_plots(k,:);
    fig = plotTabParam(plot_info,params);
    fig_name = strrep(strrep(plot_info.Name{1},' ','-'),'\mu','u');
    savefig(fig,[Path,filesep,'Figures',filesep,fig_name,'.fig'])
    saveas(fig,[Path,filesep,'CRFs',filesep,fig_name,'.png'])
    close;
end
end