%% Plot Group Cyc Avg.m
%This function makes figures with multiple cycle averages of similar
%experiments across one degree of freedom like:
% Sine w/ different frequencies and/or amplitudes
% Autoscan Current Levels

function plotGroupCycAvg(params)
Path = params.Path;
Cyc_Path = params.Cyc_Path;
code_Path = params.code_Path;
version = params.version;
Experimenter = params.Experimenter;
if isfield(params,'annot')
    annot = params.annot;
else
    annot = 1;
end
if isfield(params,'YMax')
    YMax = params.YMax;
    YMax(isnan(YMax)) = [];
else
    YMax = [];
end
close all;
load('VNELcolors.mat','colors')
code_name = ['Plotting Scripts',filesep,'plotGroupCycAvg.m'];
% Load table in question
res_file = extractfield(dir([Path,filesep,'*Results.mat']),'name')';
if isempty(res_file)
    disp('No table with cycle parameters found on this path. Creating one now.')
    rerun = ~strcmp(questdlg('If a parameter table already exists, use that one or rerun?','','Use existing table','Rerun','Rerun'),'Use existing table');
    MakeCycleSummaryTable(Path,Cyc_Path,rerun);
    res_file = extractfield(dir([Path,filesep,'*Results.mat']),'name')';
end
load(res_file{end},'all_results')
%% Sine
if any(contains(all_results.Type,'Sine'))
    all_results2 = all_results(contains(all_results.Type,'Sine'),:);
    fn = size(all_results2,1);
    % Plot Group Cycle Average Results
    cyc_files = all_results2.File;
    freqs = unique(all_results2.('Frequency(Hz)'));
    fnum = length(freqs);
    amps = unique(all_results2.('Amplitude(dps)'));
    anum = length(amps);
    exps = zeros(fnum,anum);
    for i = 1:fn
        exps(ismember(freqs,all_results2.('Frequency(Hz)')(i)),ismember(amps,all_results2.('Amplitude(dps)')(i))) = 1;
    end
    file_parts = [all_results2.Subject,all_results2.Visit,cellstr(datestr(all_results2.Date,'yyyymmdd')),all_results2.Condition,all_results2.AxisName,all_results2.Goggle,strcat(strrep(cellstr(num2str(all_results2.('Frequency(Hz)'))),' ',''),'Hz'),strcat(strrep(cellstr(num2str(all_results2.('Amplitude(dps)'))),' ',''),'dps')];
    conds1 = unique(join(file_parts(:,1:6)));
    if any(sum(exps,2)>1)
        rel_freqs = freqs(sum(exps,2)>1);
        conds_f = repmat(conds1,1,length(rel_freqs));
        for i = 1:length(rel_freqs)
            conds_f(:,i) = strcat(conds1,{[' ',num2str(rel_freqs(i)),'Hz']});
        end
        conds_f = reshape(conds_f,[],1);
    else
        conds_f = {};
    end
    if any(sum(exps,1)>1)
        rel_amps = amps(sum(exps,1)>1);
        conds_a = repmat(conds1,1,length(rel_amps));
        for i = 1:length(rel_amps)
            conds_a(:,i) = strcat(conds1,{[' ',num2str(rel_amps(i)),'dps']});
        end
        conds_a = reshape(conds_a,[],1);
    else
        conds_a = {};
    end
    conds = [conds_f;conds_a];
    if isempty(conds) %No common freq or amp
        conds = conds1;
    end
    enum = size(conds,1);
    fig_names = cell(enum,1);
    for j = 1:enum
        disp(['Making sine plot ',num2str(j),'/',num2str(enum)])
        %Set up labels
        if contains(conds(j),'Hz') %Amp Sweep
            labs = strcat(strrep(cellstr(num2str(amps)),' ',''),'dps');
        elseif contains(conds(j),'dps') %Freq Sweep
            labs = strcat(strrep(cellstr(num2str(freqs)),' ',''),'Hz');
        else %Combo
            labs = unique(join(file_parts(:,7:8)));
        end
        %Find files to plot--NEED TO FIX THIS
        rel_files = cell(length(labs),1);
        for i = 1:length(labs)
            files = cyc_files(sum(ismember(file_parts,[split(conds(j));split(labs(i))]'),2)==length([split(conds(j));split(labs(i))]));
            if ~isempty(files)
                rel_files(i) = files;
            end
        end
        %Only plot if there are more than two CycAvg for the
        %experiments
        if sum(~cellfun(@isempty,rel_files))>1
            fig_name = conds{j};
            if contains(conds{j},{'X','Y'})
                leg_text = {'Inv Stim','Left X','Right X',...
                    'Left Y','Right Y','Left Z','Right Z'};
                canals = {'lx','rx','ly','ry','lz','rz'};
            else
                leg_text = {'Inv Stim','Left LARP','Right LARP',...
                    'Left RALP','Right RALP','Left Z','Right Z'};
                canals = {'ll','rl','lr','rr','lz','rz'};
            end
            figure('Units','inches','Position',[0.2778    5.8472   17.2222    3.8333],'Color',[1,1,1])
            %Title
            annotation('textbox',[0 .9 1 .1],'String',fig_name,'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
            if annot
                annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,...
                    code_name,newline,...
                    'VOGA',version,newline,Experimenter],'FontSize',5,...
                    'EdgeColor','none','interpreter','none');
            end
            ha = gobjects(1,length(labs));
            x_space = 0.02;
            x_min = 0.04;
            x_max = 0.98;
            x_wid = (x_max-x_min-x_space*(length(labs)-1))/length(labs);
            y_height = 0.75;
            x_pos = x_min:(x_wid+x_space):x_max;
            y_pos = 0.12;
            YMax_vals = NaN(1,length(length(labs)));
            for i = 1:length(rel_files)
                ha(i) = subplot(1000,1000,1);
                ha(i).Position = [x_pos(i) y_pos x_wid y_height];
                if ~isempty(rel_files{i})
                    load([Cyc_Path,filesep,rel_files{i}],'CycAvg')
                    fields = fieldnames(CycAvg);
                    if contains(rel_files{i},{'LH','RH'})
                        rel_canals = {'lz','rz'};
                    elseif contains(rel_files{i},{'LA','RP'})
                        rel_canals = {'ll','rl'};
                    elseif contains(rel_files{i},{'RA','LP'})
                        rel_canals = {'lr','rr'};    
                    elseif contains(rel_files{i},{'X'})
                        rel_canals = {'lx','rx'}; 
                    elseif contains(rel_files{i},{'Y'})
                        rel_canals = {'ly','ry'}; 
                    else
                        rel_canals = {};
                    end
%                     if all(CycAvg.info.stim_axis == [0 0 1])||all(CycAvg.info.stim_axis == [0 0 -1])
%                         rel_canals = {'lz','rz'};
%                     elseif all(CycAvg.info.stim_axis == [1 0 0])||all(CycAvg.info.stim_axis == [-1 0 0])
%                         rel_canals = {'ll','rl'};
%                     elseif all(CycAvg.info.stim_axis == [0 1 0])||all(CycAvg.info.stim_axis == [0 -1 0])
%                         rel_canals = {'lr','rr'};    
%                     elseif all(CycAvg.info.stim_axis == [0.707 0.707 0])||all(CycAvg.info.stim_axis == [-0.707 -0.707 0])
%                         rel_canals = {'lx','rx'}; 
%                     elseif all(CycAvg.info.stim_axis == [0.707 -0.707 0])||all(CycAvg.info.stim_axis == [-0.707 0.707 0])
%                         rel_canals = {'ly','ry'}; 
%                     end
                    if ~ismember('t',fields)
                        CycAvg.t = reshape((0:1/CycAvg.Fs:(length(CycAvg.ll_cycavg)-1)/CycAvg.Fs),[],1);
                    else
                        CycAvg.t = reshape(CycAvg.t,[],1);
                    end
                    if length(CycAvg.t) > 1000
                        s = round(linspace(1,length(CycAvg.t),1000));
                    else
                        s = 1:length(CycAvg.t);
                    end
                    [aa,ab] = size(CycAvg.stim);
                    if aa == length(CycAvg.t) && ab ~=1
                        CycAvg.stim = mean(CycAvg.stim,2)';
                    elseif ab == length(CycAvg.t) && aa ~=1
                        CycAvg.stim = mean(CycAvg.stim,1);
                    end
                    h = gobjects(length(canals)+1,1);
                    h(1) = plot(CycAvg.t(s),-CycAvg.stim(s),'k');
                    hold on
                    %Now add the fills and standard deviations and means
                    max_eyevel = zeros(1,length(canals));
                    for ii = 1:length(canals)
                        max_eyevel(ii) = max(abs([reshape(CycAvg.([canals{ii},'_cycavg'])(s)-CycAvg.([canals{ii},'_cycstd'])(s),1,[]),reshape(CycAvg.([canals{ii},'_cycavg'])(s)+CycAvg.([canals{ii},'_cycstd'])(s),1,[])]));
                        fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.([canals{ii},'_cycavg'])(s)-CycAvg.([canals{ii},'_cycstd'])(s),fliplr((CycAvg.([canals{ii},'_cycavg'])(s)+CycAvg.([canals{ii},'_cycstd'])(s)))],colors.([canals{ii}(1),'_',canals{ii}(2),'_s']))
                        plot(CycAvg.t(s),CycAvg.([canals{ii},'_cycavg'])(s) + CycAvg.([canals{ii},'_cycstd'])(s),'Color',colors.([canals{ii}(1),'_',canals{ii}(2)]))
                        plot(CycAvg.t(s),CycAvg.([canals{ii},'_cycavg'])(s) - CycAvg.([canals{ii},'_cycstd'])(s),'Color',colors.([canals{ii}(1),'_',canals{ii}(2)]))
                        h(ii+1) = plot(CycAvg.t(s),CycAvg.([canals{ii},'_cycavg'])(s),'Color',colors.([canals{ii}(1),'_',canals{ii}(2)]),'LineWidth',2);
                    end
                    for ii = 1:length(rel_canals)
                        fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.([rel_canals{ii},'_cycavg'])(s)-CycAvg.([rel_canals{ii},'_cycstd'])(s),fliplr((CycAvg.([rel_canals{ii},'_cycavg'])(s)+CycAvg.([rel_canals{ii},'_cycstd'])(s)))],colors.([rel_canals{ii}(1),'_',rel_canals{ii}(2),'_s']))
                        plot(CycAvg.t(s),CycAvg.([rel_canals{ii},'_cycavg'])(s) + CycAvg.([rel_canals{ii},'_cycstd'])(s),'Color',colors.([rel_canals{ii}(1),'_',rel_canals{ii}(2)]))
                        plot(CycAvg.t(s),CycAvg.([rel_canals{ii},'_cycavg'])(s) - CycAvg.([rel_canals{ii},'_cycstd'])(s),'Color',colors.([rel_canals{ii}(1),'_',rel_canals{ii}(2)]))
                        plot(CycAvg.t(s),CycAvg.([rel_canals{ii},'_cycavg'])(s),'Color',colors.([rel_canals{ii}(1),'_',rel_canals{ii}(2)]),'LineWidth',2);
                    end
                    hold off
                else
                    h = gobjects(length(canals)+1,1);
                    h(1) = plot(NaN,NaN,'k');
                    hold on
                    for ii = 1:length(canals)
                        h(ii+1) = plot(NaN,NaN,'Color',colors.([canals{ii}(1),'_',canals{ii}(2)]),'LineWidth',2);
                    end
                    hold off
                    max_eyevel = 0;
                end
                title(labs{i})
                xlabel('Time (s)')
                YMax_vals(i) = max(max_eyevel);
                if i == 1
                    ylabel('Angular Velocity (dps)')
                    leg_reord = [2,4,6,1,3,5,7];
                    leg = legend(ha(1),h(leg_reord),leg_text(leg_reord),'NumColumns',2);
                    leg.ItemTokenSize(1) = 7;
                    leg_pos = leg.Position; %legend position
                    perc_pos = ((y_height+y_pos)-leg_pos(2))/y_height; %percentage of the axes that are just the legend
                    YMax_vals(i) = max(max_eyevel)*0.5/(0.5-perc_pos);
                end
                set(gca,'XLim',[0,CycAvg.t(end)])
            end
            if isempty(YMax)
                YLim = [-1 1]*ceil(1.1*max(YMax_vals)/10)*10; %Round up to the nearest 10;
            else
                YLim = [-YMax YMax];
            end
            set(ha,'YLim',YLim)
            fig_names{j} = [Path,filesep,'Figures',filesep,'CycleAverages_',strrep(fig_name,' ','-'),'.fig'];
            savefig(fig_names{j})
            close;
        end
    end
    fig_names(cellfun(@isempty,fig_names)) = [];
    for i = 1:length(fig_names)
        open(fig_names{i}) 
    end
end
%% Velocity Step
if any(contains(all_results.Type,'Exponential'))
    %%
    all_results2 = all_results(contains(all_results.Type,'Exponential'),:);
    % Plot Group Cycle Average Results
    cyc_files = all_results2.File;
    amps = unique(abs(all_results2.('Amplitude(dps)'))); %There should be a positive and negative version of each velocity but not needed
    anum = length(amps);
    file_parts = [all_results2.Subject,all_results2.Visit,cellstr(datestr(all_results2.Date,'yyyymmdd')),all_results2.Condition,all_results2.AxisName,all_results2.Goggle,strcat(strrep(cellstr(num2str(all_results2.('Amplitude(dps)'))),' ',''),'dps')];
    rel_file_parts = file_parts(:,1:6);
    common_cond = cell(1,size(file_parts,2));
    for i = 1:size(file_parts,2)
        if length(unique(file_parts(:,i)))==1
            common_cond(i) = unique(file_parts(:,i));
        end
    end
    rel_file_parts(:,~cellfun(@isempty,common_cond)) = [];
    common_cond(cellfun(@isempty,common_cond)) = [];
    common_cond = strjoin(common_cond);
    if size(rel_file_parts,2)>1
        conds = unique(join(rel_file_parts));
    else
        conds = unique(rel_file_parts);
    end
    if isempty(conds)
        conds = all_results2.Condition(1);
    end
    enum = length(conds);
    fig_names = cell(anum,1);
    for j = 1:anum
        disp(['Making exponential plot ',num2str(j),'/',num2str(anum)])
        rel_amp = amps(j);
        YMax = 1.1*rel_amp;
        %Set up labels
        labs = repmat(conds,1,2);
        labs(:,1) = strcat(labs(:,1),[' ',num2str(rel_amp),'dps']);
        labs(:,2) = strcat(labs(:,2),[' ',num2str(-rel_amp),'dps']);
        labs = reshape(labs,[],1);
        %Find files to plot
        rel_files = cell(length(labs),1);
        for i = 1:length(labs)
            files = cyc_files(sum(ismember(file_parts,split(labs(i))),2)==length(split(labs(i))));
            if ~isempty(files)
                rel_files(i) = files;
            end
        end
        fig_name = [common_cond,' ',num2str(rel_amp),'dps Velocity Step'];
        if all(contains(conds,{'X','Y'}))
            leg_text = {'Head Velocity','Left X','Right X',...
                'Left Y','Right Y','Left Z','Right Z'};
            canals = {'lx','rx','ly','ry','lz','rz'};
        elseif all(~contains(conds,{'X','Y'}))
            leg_text = {'Head Velocity','Left LARP','Right LARP',...
                'Left RALP','Right RALP','Left Z','Right Z'};
            canals = {'ll','rl','lr','rr','lz','rz'};
        else
            leg_text = {'Head Velocity','Left LARP','Right LARP',...
                'Left RALP','Right RALP','Left Z','Right Z',...
                'Left X','Right X','Left Y','Right Y'};
            canals = {'ll','rl','lr','rr','lz','rz','lx','rx','ly','ry'};
        end
        figure('Units','inches','Position',[0.25    0.25   17    8],'Color',[1,1,1])
        %Title
        annotation('textbox',[0 .9 1 .1],'String',fig_name,'FontSize',14,...
            'HorizontalAlignment','center','EdgeColor','none');
        if annot
            annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,...
                code_name,newline,...
                'VOGA',version,newline,Experimenter],'FontSize',5,...
                'EdgeColor','none','interpreter','none');
        end
        ha = gobjects(1,length(labs));
        x_space = 0.01;
        x_min = 0.04;
        x_max = 0.98;
        y_min = 0.08;
        y_max = 0.92;
        y_space = 0.03;
        x_wid = (x_max-x_min-x_space*(length(conds)-1))/length(conds);
        y_wid = (y_max-y_min-y_space)/2;
        x_pos = reshape(repmat(x_min:(x_wid+x_space):x_max,2,1)',[],1);
        y_pos = reshape(repmat(y_max-y_wid:-(y_wid+y_space):y_min,length(conds),1),[],1);
        for i = 1:length(labs)
            ha(i) = subplot(2,length(conds),i);
        end
        for i = 1:length(labs)
            ha(i).Position = [x_pos(i),y_pos(i),x_wid,y_wid];
        end
        for i = 1:length(labs)
            axes(ha(i))
            plot(NaN,NaN)
            if ~isempty(rel_files{i})
                load([Cyc_Path,filesep,rel_files{i}],'CycAvg')
                fields = fieldnames(CycAvg);
                if all(CycAvg.info.stim_axis == [0 0 1])||all(CycAvg.info.stim_axis == [0 0 -1])
                    rel_canals = {'lz','rz'};
                elseif all(CycAvg.info.stim_axis == [1 0 0])||all(CycAvg.info.stim_axis == [-1 0 0])
                    rel_canals = {'ll','rl'};
                elseif all(CycAvg.info.stim_axis == [0 1 0])||all(CycAvg.info.stim_axis == [0 -1 0])
                    rel_canals = {'lr','rr'};    
                elseif all(CycAvg.info.stim_axis == [0.707 0.707 0])||all(CycAvg.info.stim_axis == [-0.707 -0.707 0])
                    rel_canals = {'lx','rx'}; 
                elseif all(CycAvg.info.stim_axis == [0.707 -0.707 0])||all(CycAvg.info.stim_axis == [-0.707 0.707 0])
                    rel_canals = {'ly','ry'}; 
                end
                rel_gain = NaN(length(rel_canals),1);
                rel_tau = NaN(length(rel_canals),1);
                if ~ismember('t',fields)
                    CycAvg.t = reshape((0:1/CycAvg.Fs:(length(CycAvg.ll_cycavg)-1)/CycAvg.Fs),[],1);
                else
                    CycAvg.t = reshape(CycAvg.t,[],1);
                end
                if length(CycAvg.t) > 1000
                    s = round(linspace(1,length(CycAvg.t),1000));
                else
                    s = 1:length(CycAvg.t);
                end
                [aa,ab] = size(CycAvg.stim);
                if aa == length(CycAvg.t) && ab ~=1
                    CycAvg.stim = mean(CycAvg.stim,2)';
                elseif ab == length(CycAvg.t) && aa ~=1
                    CycAvg.stim = mean(CycAvg.stim,1);
                end
                h = gobjects(length(canals)+1,1);
                h(1) = plot(CycAvg.t(s),CycAvg.stim(s),'k');
                hold on
                %Now add the fills and standard deviations and means
                for ii = 1:length(canals)
                    h(ii+1) = plot(CycAvg.t(s),CycAvg.([canals{ii},'_cyc'])(s),'Color',colors.([canals{ii}(1),'_',canals{ii}(2)]),'LineWidth',1);
                end
                for ii = 1:length(rel_canals)
                    plot(CycAvg.t(s),CycAvg.([rel_canals{ii},'_cyc_prefilt'])(s),'.','Color',colors.([rel_canals{ii}(1),'_',rel_canals{ii}(2)]),'LineWidth',1);
                    plot(CycAvg.t(s),CycAvg.([rel_canals{ii},'_cycavg_fit'])(s),'Color',colors.([rel_canals{ii}(1),'_',rel_canals{ii}(2)]),'LineWidth',2);
                    rel_gain(ii) = abs(CycAvg.cycle_params.exp_ord1_high_params.(rel_canals{ii})(1)/rel_amp);
                    rel_tau(ii) = CycAvg.cycle_params.exp_ord1_high_params.(rel_canals{ii})(2);
                end
                hold off
                if contains(labs{i},' -')
                    text_y_pos = -YMax;
                else
                    title(strrep(labs{i},[' ',num2str(rel_amp),'dps'],''))
                    text_y_pos = -0.5*YMax;
                end
                rel_gain(isnan(rel_gain)) = 0;
                rel_tau(isnan(rel_tau)) = 0;
                text(115,text_y_pos,['Hain Product: ',num2str(round(max(rel_gain.*rel_tau),2),2)],'HorizontalAlignment','right','VerticalAlignment','bottom')
            else
                h = gobjects(length(canals)+1,1);
                h(1) = plot(NaN,NaN,'k');
                hold on
                for ii = 1:length(canals)
                    h(ii+1) = plot(NaN,NaN,'Color',colors.([canals{ii}(1),'_',canals{ii}(2)]),'LineWidth',2);
                end
                hold off
            end
        end
        leg_reord = [2,4,6,1,3,5,7];
        leg = legend(ha(1),h(leg_reord),leg_text(leg_reord),'NumColumns',2,'Location','northeast');
        leg.ItemTokenSize(1) = 7;
        xlabel(ha(end-(enum-1):end),'Time (s)')
        ylabel(ha(1:enum:end),'Angular Velocity (dps)')
        set(ha(1:end-enum),'XTickLabel',[])
        set(ha(mod(1:2*enum,enum)~=1),'YTickLabel',[])
        set(ha,'XLim',[0 120])
        set(ha(1:enum),'YLim',[-YMax/2 YMax])
        set(ha(enum+1:end),'YLim',[-YMax YMax/2])
        fig_names{j} = [Path,filesep,'Figures',filesep,'CycleAverages_',strrep(fig_name,' ','-'),'.fig'];
        savefig(fig_names{j})
    end
end
%% Autoscan
if any(contains(all_results.Condition,'Autoscan'))
    all_results = all_results(contains(all_results.Condition,'Autoscan'),:);
    electrodes = unique(all_results.Electrode);
    e_num = regexp(electrodes,'\d*','Match');
    e_num = str2double(vertcat(e_num{:}));
    [~,ia] = sort(e_num);
    electrodes = electrodes(ia);
    %Select Files to Plot
    canal = listdlg('PromptString','Which three electrodes?',...
        'ListString',electrodes,...
        'SelectionMode','multiple');
    i1 = find(contains(all_results.Electrode,electrodes(canal(1))));
    i2 = find(contains(all_results.Electrode,electrodes(canal(2))));
    i3 = find(contains(all_results.Electrode,electrodes(canal(3))));
    if ~(length(i3)==length(i1)&&length(i1)==length(i2))
        error('Unequal number of files for each electrode')
    end
    [curr1,i11] = sort(all_results.('CurrentAmp(uA)')(i1));
    i1 = i1(i11);
    [curr2,i22] = sort(all_results.('CurrentAmp(uA)')(i2));
    i2 = i2(i22);
    [curr3,i33] = sort(all_results.('CurrentAmp(uA)')(i3));
    i3 = i3(i33);
    row1 = all_results.File(i1);
    row2 = all_results.File(i2);
    row3 = all_results.File(i3);
    f_order = [row1;row2;row3];
    n_col = length(row1);  
    curr_lab = cellstr(num2str([curr1;curr2;curr3]));
    curr_lab{1} = [curr_lab{1},'\muA'];
    curr_lab{1+n_col} = [curr_lab{1+n_col},'\muA'];
    curr_lab{1+2*n_col} = [curr_lab{1+2*n_col},'\muA'];
    % Make the figure name
    file_parts = [all_results.Subject,all_results.Visit,cellstr(datestr(all_results.Date,'yyyymmdd')),all_results.Condition,all_results.Goggle,strcat(strrep(cellstr(num2str(all_results.('PulseFreq(pps)'))),' ',''),'pps'),strcat(strrep(cellstr(num2str(all_results.('PhaseDur(us)'))),' ',''),'us')];   
    fig_title = {strjoin([join(unique(file_parts,'stable')),strjoin(electrodes(canal),' ')])};
    % Plot Current Levels
    fig = figure;
    fig.Color = [1,1,1];
    fig.Units = 'inches';
    fig.Position = [1 1 8 4];
    if annot
        annotation('textbox',[0 .9 1 .1],'String',fig_title,'FontSize',14,...
            'HorizontalAlignment','center','EdgeColor','none');
        annotation('textbox',[0 0 1 1],'String',['VOGA',version,...
            newline,Experimenter],'FontSize',5,...
            'EdgeColor','none','interpreter','none');
        annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,...
            code_name],'FontSize',5,...
            'EdgeColor','none','interpreter','none','VerticalAlignment','bottom');
    end
    ha = gobjects(1,length(f_order));
    %Set params
    grid_on = false;
    if isempty(YMax)
        YMax = 100;
    end
    YLim = YMax*[-1 1];
    x_min = 0.01;
    x_max = 0.95;
    space_x = 0.01;
    y_min = 0.08;
    y_max = 0.92;
    space_y = 0.03;
    %Calculate
    x_wid = (x_max - x_min - space_x*(n_col-1))/n_col;
    fig_row_pos = repmat(x_min:(x_wid+space_x):x_max,1,3);
    y_wid = (y_max - y_min - space_y*2)/3;
    fig_col_pos = reshape(repmat(fliplr(y_min:(y_wid+space_y):y_max),n_col,1),[],1)';
    annotation('line',fig_row_pos(end)+[0 x_wid],y_min-0.01*[1 1],'LineWidth',2)
    annotation('line',(x_max+space_x)*[1 1],y_min+[0 YMax*y_wid/(2*YLim(2))],'LineWidth',2)
    annotation('textbox','String','0.5s','EdgeColor','none',...
        'Position',[fig_row_pos(end),0,x_wid,y_min-0.01],'HorizontalAlignment','right','VerticalAlignment','middle')
    annotation('textbox','String',[num2str(YMax),newline,'\circ/s'],'EdgeColor','none',...
        'Position',[x_max+space_x,y_min,1-(x_max+space_x),YMax*y_wid/(2*YLim(2))],'HorizontalAlignment','center','VerticalAlignment','middle')
    for i = 1:length(f_order)
        ha(i) = subplot(3,n_col,i);
        set(gca,'XColor','none','YColor','none')
    end
    for i = 1:length(f_order)
        axes(ha(i))
        set(ha(i),'Position',[fig_row_pos(i),fig_col_pos(i),x_wid,y_wid])
        if mod(i,n_col) > 0
            annotation('line',(fig_row_pos(i)+x_wid+0.5*space_x)*[1 1],fig_col_pos(i)+[0 y_wid],'LineWidth',1,'LineStyle','--')
        end
        %Load and plot
        b = load([Cyc_Path,filesep,f_order{i}]);
        a = fieldnames(b);
        CycAvg = b.(a{1});
        fields = fieldnames(CycAvg);
        if ~ismember('t',fields)
            CycAvg.t = (0:1/CycAvg.Fs:(length(CycAvg.ll_cycavg)-1)/CycAvg.Fs)';
        end
        if length(CycAvg.t) > 1000
            s = round(linspace(1,length(CycAvg.t),1000));
        else
            s = 1:length(CycAvg.t);
        end
        hold on
        %Now add the fills and standard deviations and means
        %Plot the intended canal again so that it's in the foreground
        if contains(f_order{i},'LP') || contains(f_order{i},'RA') %RALP
            curr_col = colors.l_r;
            %LE-LHRH
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
            %RE-LHRH
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
            %LE-LARP
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s),'Color',colors.l_l)
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s),'Color',colors.l_l)
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s),'Color',colors.l_l,'LineWidth',2);
            %RE-LARP
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s),'Color',colors.r_l)
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s),'Color',colors.r_l)
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s),'Color',colors.r_l,'LineWidth',2);
            %LE_RALP
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s),'Color',colors.l_r)
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s),'Color',colors.l_r)
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s),'Color',colors.l_r,'LineWidth',2);
            %RE-RALP
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s),'Color',colors.r_r)
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s),'Color',colors.r_r)
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s),'Color',colors.r_r,'LineWidth',2);
        elseif contains(f_order{i},'LH') || contains(f_order{i},'RH') %LHRH
            curr_col = colors.l_z;
            %LE-LARP
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s),'Color',colors.l_l)
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s),'Color',colors.l_l)
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s),'Color',colors.l_l,'LineWidth',2);
            %RE-LARP
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s),'Color',colors.r_l)
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s),'Color',colors.r_l)
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s),'Color',colors.r_l,'LineWidth',2);
            %LE_RALP
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s),'Color',colors.l_r)
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s),'Color',colors.l_r)
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s),'Color',colors.l_r,'LineWidth',2);
            %RE-RALP
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s),'Color',colors.r_r)
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s),'Color',colors.r_r)
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s),'Color',colors.r_r,'LineWidth',2);
            %LE-LHRH
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
            %RE-LHRH
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
        elseif contains(f_order{i},'RP') || contains(f_order{i},'LA') %LARP
            curr_col = colors.l_l;
            %LE-LHRH
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
            %RE-LHRH
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
            %LE_RALP
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s),'Color',colors.l_r)
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s),'Color',colors.l_r)
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s),'Color',colors.l_r,'LineWidth',2);
            %RE-RALP
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s),'Color',colors.r_r)
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s),'Color',colors.r_r)
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s),'Color',colors.r_r,'LineWidth',2);
            %LE-LARP
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s),'Color',colors.l_l)
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s),'Color',colors.l_l)
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s),'Color',colors.l_l,'LineWidth',2);
            %RE-LARP
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s),'Color',colors.r_l)
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s),'Color',colors.r_l)
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s),'Color',colors.r_l,'LineWidth',2);
        end
        hold off
        axis([0 0.5 YLim])
        text(0.5,YLim(2),curr_lab{i},'Color',curr_col,'HorizontalAlignment','right','VerticalAlignment','top')
        if(grid_on)
            set(gca,'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
        end
        if mod(i,n_col)==1
            text(0.5,YLim(1),['n=',num2str(length(CycAvg.cyclist))],'Color','k','HorizontalAlignment','right','VerticalAlignment','bottom')
        else
            text(0.5,YLim(1),num2str(length(CycAvg.cyclist)),'Color','k','HorizontalAlignment','right','VerticalAlignment','bottom')
        end
    end
    savefig([Path,filesep,'Figures',filesep,strrep(fig_title{:},' ','-'),'.fig'])
end
%% Impulses
if any(contains(all_results.Type,'Impulse')) 
    all_results = all_results(contains(all_results.Type,'Impulse'),:);
    fn = size(all_results,1);
    % Plot Group Cycle Average Results
    cyc_files = all_results.File;
    canal = cell(length(cyc_files),1);
    canal(contains(cyc_files,'LH')) = {'LH'};
    canal(contains(cyc_files,'RH')) = {'RH'};
    canal(contains(cyc_files,'LA')) = {'LA'};
    canal(contains(cyc_files,'LP')) = {'LP'};
    canal(contains(cyc_files,'RA')) = {'RA'};
    canal(contains(cyc_files,'RP')) = {'RP'};    
    file_parts = [all_results.Subject,all_results.Visit,cellstr(datestr(all_results.Date,'yyyymmdd')),all_results.Experiment,all_results.Condition,all_results.Goggle];
    joined_fparts = join(file_parts);
    conds = strcat(unique(joined_fparts),{' Impulses'});
    %conds1 = [strcat(unique(joined_fparts),{' Left Ear Impulses'});strcat(unique(joined_fparts),{' Right Ear Impulses'})];
    enum = size(conds,1); 
    rel_file_parts = file_parts;
    common_cond = cell(1,size(file_parts,2));
    for i = 1:size(file_parts,2)
        if length(unique(file_parts(:,i)))==1
            common_cond(i) = unique(file_parts(:,i));
        end
    end
    rel_file_parts(:,~cellfun(@isempty,common_cond)) = [];
    common_cond(cellfun(@isempty,common_cond)) = [];
    common_cond = strjoin(common_cond); 
    if isempty(rel_file_parts) %only one condition
        conds2 = all_results.Condition(1);
        IC = ones(fn,1);
    elseif size(rel_file_parts,2)==1
        [conds2,~,IC] = unique(rel_file_parts,'stable');
    else
        [conds2,~,IC] = unique(join(rel_file_parts),'stable');
    end
    fig_names = cell(enum+2,1);
    if isempty(YMax)
        YMax = 250;
    end
    for j = 1:enum
        fig_name = conds{j};
        disp(['Making plot ',num2str(j),'/',num2str(enum+2),': ',fig_name])
        %Set up labels
        labs = {'Left Anterior (LA)','Left Posterior (LP)','Left Horizontal (LH)','Right Posterior (RP)','Right Anterior (RA)','Right Horizontal (RH)'};
        %Find files to plot
        rel_files = cell(length(labs),1);
        for i = 1:length(labs)
            files = cyc_files(contains(joined_fparts,strrep(conds(j),' Impulses',''))&contains(canal,labs{i}(end-2:end-1))); 
            if ~isempty(files)
                rel_files(i) = files(end);
            end
        end
        % Plot all canals in both ear for each condition
        %Only plot if there is more than one file to plot
        if sum(~cellfun(@isempty,rel_files))>1            
            leg_text = {'Stim','Left LARP','Right LARP',...
                'Left RALP','Right RALP','Left X','Right X',...
                    'Left Y','Right Y','Left Z','Right Z'};
            canals = {'ll','rl','lr','rr','lx','rx','ly','ry','lz','rz'};
            figure('Units','inches','Position',[0.2778    2.8472   17.2222    7.6666],'Color',[1,1,1])
            %Title
            annotation('textbox',[0 .9 1 .1],'String',fig_name,'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
            if annot
                annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,...
                    code_name,newline,...
                    'VOGA',version,newline,Experimenter],'FontSize',5,...
                    'EdgeColor','none','interpreter','none');
            end
            ha = gobjects(1,length(labs));
            x_space = 0.01;
            x_min = 0.04;
            x_max = 0.98;
            y_min = 0.08;
            y_max = 0.92;
            y_space = 0.03;
            x_wid = (x_max-x_min-2*x_space)/3;
            y_wid = (y_max-y_min-x_space)/2;
            x_pos = repmat(x_min:(x_wid+x_space):x_max,1,2);
            y_pos = reshape(repmat(fliplr(y_min:(y_wid+y_space):y_max),3,1),[],1)';           
            for i = 1:length(rel_files)
                ha(i) = subplot(2,3,i);
            end  
            leg_ind = false(1,length(leg_text));
            leg_ind(1) = true;
            for i = 1:length(rel_files)
                axes(ha(i))
                ha(i).Position = [x_pos(i) y_pos(i) x_wid y_wid];
                if ~isempty(rel_files{i})
                    load([Cyc_Path,filesep,rel_files{i}],'CycAvg')
                    fields = fieldnames(CycAvg);
                    if ~ismember('t',fields)
                        CycAvg.t = reshape(0:1/CycAvg.Fs:(length(CycAvg.ll_cycavg)-1)/CycAvg.Fs,[],1);
                    else
                        CycAvg.t = reshape(CycAvg.t,1,[]);
                    end
                    if -min(mean(CycAvg.stim))>max(mean(CycAvg.stim))
                        invh = -1;
                        inve = 1;
                    else
                        invh = 1;
                        inve = -1;
                    end
                    if all(CycAvg.info.stim_axis == [0 0 1])||all(CycAvg.info.stim_axis == [0 0 -1])
                        rel_canals = {'lz','rz'};
                    elseif all(CycAvg.info.stim_axis == [1 0 0])||all(CycAvg.info.stim_axis == [-1 0 0])
                        rel_canals = {'ll','rl'};
                    elseif all(CycAvg.info.stim_axis == [0 1 0])||all(CycAvg.info.stim_axis == [0 -1 0])
                        rel_canals = {'lr','rr'};    
                    elseif all(CycAvg.info.stim_axis == [0.707 0.707 0])||all(CycAvg.info.stim_axis == [-0.707 -0.707 0])
                        rel_canals = {'lx','rx'}; 
                    elseif all(CycAvg.info.stim_axis == [0.707 -0.707 0])||all(CycAvg.info.stim_axis == [-0.707 0.707 0])
                        rel_canals = {'ly','ry'}; 
                    end
                    h = gobjects(length(canals)+1,1);
                    h(1) = plot(NaN,NaN,'k','LineWidth',1);
                    hold on
                    plot(CycAvg.t,invh*CycAvg.stim,'k','LineWidth',0.5);
                    for ii = 1:length(canals)
                        trace = canals{ii};                       
                         h(ii+1) = plot(NaN,NaN,'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',1);
                    end                   
%                     for ii = 1:length(canals)
%                         trace = canals{ii};                       
%                         h(ii+1) = plot(NaN,NaN,'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',1);
%                         if isfield(CycAvg,[trace,'_cyc_prefilt'])
%                             leg_ind(ii+1) = true;
%                             plot(CycAvg.t,inve*CycAvg.([trace,'_cyc_prefilt']),'Color',colors.([trace(1),'_',trace(2),'_s']),'LineWidth',0.5);
%                             plot(CycAvg.t,inve*CycAvg.([trace,'_cyc']),'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',1);                        
%                         elseif isfield(CycAvg,[trace,'_cyc'])
%                             leg_ind(ii+1) = true;
%                             plot(CycAvg.t,inve*CycAvg.([trace,'_cyc']),'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',1);                        
%                         end
%                     end 
                    gain = NaN(1,length(rel_canals));
                    lat = NaN(1,length(rel_canals));
                    for ii = 1:length(rel_canals)
                        trace = rel_canals{ii};
                        if isfield(CycAvg,[trace,'_cyc'])
                            leg_ind(find(ismember(canals,trace))+1) = true;
                            plot(CycAvg.t,inve*CycAvg.([trace,'_cyc_prefilt']),'Color',colors.([trace(1),'_',trace(2),'_s']),'LineWidth',0.5);
                            plot(CycAvg.t,inve*CycAvg.([trace,'_cyc']),'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',1);
                        end
                        gain(ii) = CycAvg.parameterized.(['Gain_',upper(trace),'_HIGH']);
                        lat(ii) = CycAvg.parameterized.(['Latency_',upper(trace)]);
                    end 
                    if isfield(CycAvg,'stim_start')
                        xline(CycAvg.stim_start,'k','LineWidth',2)
                        xline(CycAvg.stim_end,'k','LineWidth',2)
                        xline(CycAvg.trace_start,'k--','LineWidth',2)
%                         saccades = [CycAvg.cycle_params.Saccade1;CycAvg.cycle_params.Saccade2]/1000+CycAvg.stim_start;
%                         saccades(isnan(saccades)) = [];
%                         if ~isempty(saccades)
%                             plot(saccades,0.98*YMax*ones(1,length(saccades)),'k*') 
%                         end
                    end                    
                    hold off
                    text(0.99*CycAvg.t(end),-YMax/2,['Gain: ',num2str(mean(gain,'omitnan'),2),newline,'Latency (ms): ',num2str(round(mean(lat,'omitnan')),3)],'HorizontalAlignment','right','VerticalAlignment','bottom')
                else
                    h = gobjects(length(canals)+1,1);
                    h(1) = plot(NaN,NaN,'k');
                    hold on
                    for ii = 1:length(canals)
                        h(ii+1) = plot(NaN,NaN,'Color',colors.([canals{ii}(1),'_',canals{ii}(2)]),'LineWidth',1);
                    end
                    hold off
                end
                title(labs{i})
            end
            ylabel(ha(1),'Angular Velocity (dps)')
            ylabel(ha(4),'Angular Velocity (dps)')
            xlabel(ha(4),'Time (s)')
            xlabel(ha(5),'Time (s)')
            xlabel(ha(6),'Time (s)')
            set(ha([2,3,5,6]),'YTickLabel',[])
            set(ha([1,2,3]),'XTickLabel',[])
            set(ha,'YLim',[-YMax/2 YMax],'XLim',[0,CycAvg.t(end)])
            leg = legend(ha(1),h(leg_ind),leg_text(leg_ind),'NumColumns',2);
            leg.ItemTokenSize(1) = 7;
            fig_names{j} = [Path,filesep,'Figures',filesep,'CycleAverages_',strrep(fig_name,' ','-'),'.fig'];
            savefig(fig_names{j})
            close;
        end
    end
    % Plot each condition for each ear
    % ADD CODE HERE
    for j = 1:2 %Each ear
        if j == 1
            fig_name = [common_cond,' Left Ear'];
            ylabs = {'Left Horizontal (LH)','Left Anterior (LA)','Left Posterior (LP)'};            
        else
            fig_name = [common_cond,' Right Ear'];
            ylabs = {'Right Horizontal (RH)','Right Posterior (RP)','Right Anterior (RA)'};
        end
        disp(['Making plot ',num2str(enum+j),'/',num2str(enum+2),': ',fig_name])
        %Set up labels
        labs = conds2;
        long_ylabs = reshape(repmat(ylabs,length(labs),1),[],1);
        %Find files to plot
        rel_files = cell(3,length(labs));
        for i = 1:length(labs)
            for c = 1:3
                files = cyc_files(IC==i&contains(canal,ylabs{c}(end-2:end-1))); 
                if ~isempty(files)
                    rel_files(c,i) = files(end);
                end
            end
        end
        rel_files = reshape(rel_files',[],1);
        % Plot all canals in one ear for each condition
        %Only plot if there is at least one one file to plot
        if sum(~cellfun(@isempty,rel_files))>0            
            leg_text = {'Stim','Left LARP','Right LARP',...
                'Left RALP','Right RALP','Left X','Right X',...
                    'Left Y','Right Y','Left Z','Right Z'};
            leg_ind = false(1,length(leg_text));
            leg_ind(1) = true;
            canals = {'ll','rl','lr','rr','lx','rx','ly','ry','lz','rz'};
            figure('Units','inches','Position',[0.2778    2.8472   17.2222    7.6666],'Color',[1,1,1])
            %Title
            annotation('textbox',[0 .9 1 .1],'String',fig_name,'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
            if annot
                annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,...
                    code_name,newline,...
                    'VOGA',version,newline,Experimenter],'FontSize',5,...
                    'EdgeColor','none','interpreter','none');
            end
            ha = gobjects(length(rel_files),1);
            x_space = 0.02;
            y_space = 0.04;
            x_min = 0.04;
            x_max = 0.98;
            x_wid = (x_max-x_min-x_space*(length(labs)-1))/length(labs);
            x_pos = reshape(repmat(x_min:(x_wid+x_space):x_max,3,1)',[],1);
            y_min = 0.08;
            y_max = 0.90;
            y_wid = (y_max-y_min-2*y_space)/3;
            y_pos = reshape(repmat(fliplr(y_min:(y_wid+y_space):y_max),length(labs),1),[],1);
            for i = 1:length(rel_files)
                ha(i) = subplot(3,length(labs),i);
            end
            for i = 1:length(rel_files)
                axes(ha(i))
                ha(i).Position = [x_pos(i) y_pos(i) x_wid y_wid];
                if ~isempty(rel_files{i})
                    load([Cyc_Path,filesep,rel_files{i}],'CycAvg')
                    if isfield(CycAvg,'rz_cycavg')
                        CycAvg.t = reshape(0:1/CycAvg.Fs:(length(CycAvg.rz_cycavg)-1)/CycAvg.Fs,1,[]);
                    else
                        CycAvg.t = reshape(0:1/CycAvg.Fs:(length(CycAvg.lz_cycavg)-1)/CycAvg.Fs,1,[]);
                    end
                    if -min(mean(CycAvg.stim))>max(mean(CycAvg.stim))
                        invh = -1;
                        inve = 1;
                    else
                        invh = 1;
                        inve = -1;
                    end
                    if all(CycAvg.info.stim_axis == [0 0 1])||all(CycAvg.info.stim_axis == [0 0 -1])
                        rel_canals = {'lz','rz'};
                    elseif all(CycAvg.info.stim_axis == [1 0 0])||all(CycAvg.info.stim_axis == [-1 0 0])
                        rel_canals = {'ll','rl'};
                    elseif all(CycAvg.info.stim_axis == [0 1 0])||all(CycAvg.info.stim_axis == [0 -1 0])
                        rel_canals = {'lr','rr'};    
                    elseif all(CycAvg.info.stim_axis == [0.707 0.707 0])||all(CycAvg.info.stim_axis == [-0.707 -0.707 0])
                        rel_canals = {'lx','rx'}; 
                    elseif all(CycAvg.info.stim_axis == [0.707 -0.707 0])||all(CycAvg.info.stim_axis == [-0.707 0.707 0])
                        rel_canals = {'ly','ry'}; 
                    end
                    h = gobjects(length(canals)+1,1);
                    h(1) = plot(NaN,NaN,'k','LineWidth',1);
                    hold on
                    plot(CycAvg.t,invh*CycAvg.stim,'k','LineWidth',0.5);
%                     for ii = 1:length(canals)
%                         trace = canals{ii};
%                         h(ii+1) = plot(NaN,NaN,'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',1);
%                         if isfield(CycAvg,[trace,'_cyc_prefilt'])
%                             leg_ind(ii+1) = true;
%                             plot(CycAvg.t,inve*CycAvg.([trace,'_cyc_prefilt']),'Color',colors.([trace(1),'_',trace(2),'_s']),'LineWidth',0.5);
%                             plot(CycAvg.t,inve*CycAvg.([trace,'_cyc']),'Color',colors.([trace(1),'_',trace(2),'']),'LineWidth',1);                  
%                         elseif isfield(CycAvg,[trace,'_cyc'])
%                             leg_ind(ii+1) = true;
%                             plot(CycAvg.t,inve*CycAvg.([trace,'_cyc']),'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',1);
%                         end
%                     end  
                    gain = NaN(1,length(rel_canals));
                    lat = NaN(1,length(rel_canals));
                    %Plot the most salient canal again
                    for ii = 1:length(rel_canals)
                        trace = rel_canals{ii};
                        if isfield(CycAvg,[trace,'_cyc_prefilt'])
                            plot(CycAvg.t,inve*CycAvg.([trace,'_cyc_prefilt']),'Color',colors.([trace(1),'_',trace(2),'_s']),'LineWidth',0.5);
                            plot(CycAvg.t,inve*CycAvg.([trace,'_cyc']),'Color',colors.([trace(1),'_',trace(2),'']),'LineWidth',1);                  
                        elseif isfield(CycAvg,[trace,'_cyc'])
                            plot(CycAvg.t,inve*CycAvg.([trace,'_cyc']),'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',1);
                        end
                        gain(ii) = CycAvg.parameterized.(['Gain_',upper(trace),'_HIGH']);
                        lat(ii) = CycAvg.parameterized.(['Latency_',upper(trace)]);
                    end 
                    if isfield(CycAvg,'stim_start')
                        xline(CycAvg.stim_start,'k','LineWidth',2)
                        xline(CycAvg.stim_end,'k','LineWidth',2)
                        xline(CycAvg.trace_start,'k--','LineWidth',2)
                    end
                    hold off
                    text(0.99*CycAvg.t(end),-YMax/2,['Gain: ',num2str(mean(gain,'omitnan'),2),newline,'Latency (ms): ',num2str(round(mean(lat,'omitnan')),3)],'HorizontalAlignment','right','VerticalAlignment','bottom')
                else
                    h = gobjects(length(canals)+1,1);
                    h(1) = plot(NaN,NaN,'k');
                    hold on
                    for ii = 1:length(canals)
                        h(ii+1) = plot(NaN,NaN,'Color',colors.([canals{ii}(1),'_',canals{ii}(2)]),'LineWidth',1);
                    end
                    hold off
                end 
                if i<=length(labs)
                    title(ha(i),[labs(i);long_ylabs(i)])
                else
                    title(ha(i),long_ylabs(i))
                end
            end
            set(ha,'XLim',[0,CycAvg.t(end)],'YLim',[-YMax/2 YMax])
            plot_num = reshape(1:length(rel_files),[],3)';
            if size(plot_num,2)>1
                set(ha(plot_num(:,2:end)),'YTickLabel',[])  
            end
            set(ha(plot_num(1:end-1,:)),'XTickLabel',[])
            ylabel(ha(plot_num(1,1)),'Angular Velocity (dps)')
            ylabel(ha(plot_num(2,1)),'Angular Velocity (dps)')
            ylabel(ha(plot_num(3,1)),'Angular Velocity (dps)')
            for i = 1:size(plot_num,2)
                xlabel(ha(plot_num(end,i)),'Time (s)')
            end
            leg = legend(ha(1),h(leg_ind),leg_text(leg_ind),'NumColumns',2);
            leg.ItemTokenSize(1) = 7;
            fig_names{j+enum} = [Path,filesep,'Figures',filesep,'CycleAverages_',strrep(fig_name,' ','-'),'.fig'];
            savefig(fig_names{j+enum})
            close;
        end
    end
    fig_names(cellfun(@isempty,fig_names)) = [];
    for i = 1:length(fig_names)
        open(fig_names{i}) 
    end
end
end