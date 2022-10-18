%% Plot Param Results.m
% This function makes figures from the Results.mat tables.

function plotParamResults(params)
Path = params.Path;
Cyc_Path = params.Cyc_Path;
code_Path = params.code_Path;
version = params.version;
Experimenter = params.Experimenter;
sub_info = params.sub_info;
Subs = sub_info{:,1};
Ears = sub_info{:,2};
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
% Initialize
close all;
load('VNELcolors.mat','colors')
code_name = ['Plotting Scripts',filesep,'plotParamResults.m'];
% Load table in question
res_file = extractfield(dir([Path,filesep,'*Results.mat']),'name')';
if isempty(res_file)
    disp('No table with cycle parameters found on this path. Creating one now.')
    rerun = ~strcmp(questdlg('If a parameter table already exists, use that one or rerun?','','Use existing table','Rerun','Rerun'),'Use existing table');
    MakeCycleSummaryTable(Path,Cyc_Path,rerun);
    res_file = extractfield(dir([Path,filesep,'*Results.mat']),'name')';
end
load(res_file{end},'all_results')
if any(contains(all_results.Type,'Sine'))
    %% Sinusoids
    all_canals = {'LARP','RALP','LHRH','Y','X'}; %Preferred order
    ear_canals = {'LARP','RALP','LHRH'};
    all_results = all_results(contains(all_results.Type,'Sine')&contains(all_results.AxisName,all_canals),:);
%     if any(contains(all_results.Experiment,'RotaryChair'))
%         load('RotaryChairNormativeData.mat','norm_dat');
%         freq = norm_dat.freq;
%         norm_gain_m = norm_dat.gain;
%         norm_gain_std = norm_dat.gain_std;
%         norm_phase_m = norm_dat.phase;
%         norm_phase_std = norm_dat.phase_std;
%     end
    fn = size(all_results,1);
    % Plot Group Cycle Average Results
    canals = all_canals(ismember(all_canals,unique(all_results.AxisName)));
    cnum = length(canals);
    canal_let = lower(strrep(strrep(strrep(canals,'LARP','L'),'RALP','R'),'LHRH','Z'));
    freqs = unique(all_results.('Frequency(Hz)'));
    fnum = length(freqs);
    amps = unique(all_results.('Amplitude(dps)'));
    anum = length(amps);
    exps = zeros(fnum,anum);
    maxvel_poscyc = NaN(1,fn);
    maxvel_poscyc_sd = NaN(1,fn);
    maxvel_negcyc = NaN(1,fn);
    maxvel_negcyc_sd = NaN(1,fn);
    phase = NaN(1,fn);
    phase_sd = NaN(1,fn);
    align_poscyc = NaN(1,fn);
    align_poscyc_sd = NaN(1,fn);
    align_negcyc = NaN(1,fn);
    align_negcyc_sd = NaN(1,fn);
    for i = 1:fn
        exps(ismember(freqs,all_results.('Frequency(Hz)')(i)),ismember(amps,all_results.('Amplitude(dps)')(i))) = 1;
        rel_col = strrep(strrep(strrep(all_results.AxisName{i},'LARP','L'),'RALP','R'),'LHRH','Z');
        [~,eye] = max(abs([all_results.(['MaxVel_L',rel_col,'_HIGH'])(i),all_results.(['MaxVel_R',rel_col,'_HIGH'])(i),all_results.(['MaxVel_L',rel_col,'_LOW'])(i),all_results.(['MaxVel_R',rel_col,'_LOW'])(i)]));
        if mod(eye,2)==1 %Left
            eye_s = 'L';
        else %Right
            eye_s = 'R';
        end
        maxvel_poscyc(i) = all_results.(['MaxVel_',eye_s,rel_col,'_HIGH'])(i);
        maxvel_poscyc_sd(i) = all_results.(['MaxVel_',eye_s,rel_col,'_HIGH_sd'])(i);
        maxvel_negcyc(i) =  all_results.(['MaxVel_',eye_s,rel_col,'_LOW'])(i);
        maxvel_negcyc_sd(i) = all_results.(['MaxVel_',eye_s,rel_col,'_LOW_sd'])(i);
        phase(i) = all_results.(['Phase_',eye_s])(i);
        phase_sd(i) = all_results.(['Phase_',eye_s,'_sd'])(i);
        align_poscyc(i) = all_results.(['Align_',eye_s,'_HIGH'])(i);
        align_poscyc_sd(i) = all_results.(['Align_',eye_s,'_HIGH_sd'])(i);
        align_negcyc(i) =  all_results.(['Align_',eye_s,'_LOW'])(i);
        align_negcyc_sd(i) = all_results.(['Align_',eye_s,'_LOW_sd'])(i);
    end
    file_parts = [all_results.Subject,all_results.Visit,cellstr(datestr(all_results.Date,'yyyymmdd')),...
        all_results.Condition,all_results.Goggle,strcat(strrep(cellstr(num2str(all_results.('Frequency(Hz)'))),' ',''),'Hz'),...
        strcat(strrep(cellstr(num2str(all_results.('Amplitude(dps)'))),' ',''),'dps')];
    rel_file_parts = file_parts(:,1:5); %everything but freq/amp
    common_cond = cell(1,5);
    for i = 1:5
        if length(unique(file_parts(:,i)))==1
            common_cond(i) = unique(file_parts(:,i));
        end
    end
    rel_file_parts(:,~cellfun(@isempty,common_cond)) = [];
    common_cond(cellfun(@isempty,common_cond)) = [];
    common_cond = strjoin(common_cond); 
    if isempty(rel_file_parts) %only one condition
        conds = all_results.Condition(1);
        IC = ones(fn,1);
    elseif size(rel_file_parts,2)==1
        [conds,~,IC] = unique(rel_file_parts,'stable');
    else
        [conds,~,IC] = unique(join(rel_file_parts),'stable');
    end
    enum = length(conds);
    %Position information for the graph types
    max_y = 0.80;
    min_y = 0.15;
    %spac_y = 0.01;
    min_x = 0.08;
    max_x = 0.90;
    spac_x = 0.01;
    wid_x = (max_x-min_x-2*spac_x)/3;
    x_pos = min_x:(wid_x+spac_x):max_x-wid_x; 
    wid_y = max_y-min_y;
    y_pos = min_y;    
    if isempty(YMax)
        YMax = 5*ceil(max([maxvel_poscyc+maxvel_poscyc_sd,maxvel_negcyc+maxvel_negcyc_sd])/5);
    end
    all_markers = 'o*.sd^v><ph_|';
    if any(sum(exps,1)>1) %Frequency Sweep
        rel_amps = amps(sum(exps,1)>1);
        for a = 1:length(rel_amps)            
            fig_name = [common_cond,' ',num2str(rel_amps(a)),'dps Sine Frequency Sweep'];
            figure('Units','inch','Position',[2 2 7 3],'Color',[1,1,1]);
            annotation('textbox',[0 0 1 1],'String',[Path,newline,code_Path,filesep,...
            'plotParamResults.m',newline,...
            'VOGA',version,newline,Experimenter],'FontSize',5,...
            'EdgeColor','none','interpreter','none','Color',[0.5,0.5,0.5]);
            annotation('textbox',[0 .9 1 .1],'String',fig_name,'FontSize',14,...
                    'HorizontalAlignment','center','EdgeColor','none');
            ha = gobjects(1,3);
            h1 = gobjects(1,cnum);
            h2 = gobjects(1,enum);
            ha(1) = subplot(1,3,1);
            ha(1).Position = [x_pos(1),y_pos,wid_x,wid_y];
            ha(2) = subplot(1,3,2);
            ha(2).Position = [x_pos(2),y_pos,wid_x,wid_y];
            ha(3) = subplot(1,3,3);
            ha(3).Position = [x_pos(3),y_pos,wid_x,wid_y];
            c_bool = false(1,cnum);
            e_bool = false(1,enum);
            for c = 1:cnum
                for i = 1:enum
                    rel_i = find(all_results.('Amplitude(dps)')==rel_amps(a)&contains(all_results.AxisName,canals{c})&IC==i);
                    if length(rel_i)>1
                        [f_ax,f_i] = sort(all_results.('Frequency(Hz)')(rel_i)); 
                        rel_i = rel_i(f_i);
                        axes(ha(1))
                        hold on
                        errorbar(f_ax,maxvel_poscyc(rel_i),maxvel_poscyc_sd(rel_i),[all_markers(i),'-'],'Color',colors.(['l_',canal_let{c}]));
                        hold off
                        axes(ha(2))
                        hold on
                        errorbar(f_ax,maxvel_negcyc(rel_i),maxvel_negcyc_sd(rel_i),[all_markers(i),'-'],'Color',colors.(['l_',canal_let{c}]));
                        h1(c) = plot(NaN,NaN,'LineWidth',2,'Color',colors.(['l_',canal_let{c}]));
                        hold off
                        axes(ha(3))
                        hold on
                        errorbar(f_ax,phase(rel_i),phase_sd(rel_i),[all_markers(i),'-'],'Color',colors.(['l_',canal_let{c}]));
                        h2(i) = plot(NaN,NaN,['k',all_markers(i)]);
                        hold off
                        c_bool(c) = true;
                        e_bool(i) = true;
                    end
                end
            end            
            linkaxes(ha,'x')
            linkaxes(ha(1:2),'y')
            set(ha,'xscale','log','XTick',freqs,'XTickLabelRotation',0,'XMinorTick','off','XLim',[min(freqs)*0.8 max(freqs)*1.2])
            set(ha(1:2),'YLim',[0 YMax])
            set(ha(2),'YTickLabel',[])
            set(ha(3),'YAxisLocation','right','YLim',[-180 360])
            title(ha(1),'Positive Half-Cycle')
            title(ha(2),'Negative Half-Cycle')
            title(ha(3),'Phase Lead')
            ylabel(ha(1),'Maxmimum Eye Velocity (dps)')
            ylabel(ha(3),'Phase (deg)')
            xlabel(ha(2),'Frequency (Hz)')
            %Make legends
            axes(ha(2))
            hold on
            for c = 1:cnum
                h1(c) = plot(NaN,NaN,'LineWidth',2,'Color',colors.(['l_',canal_let{c}]));
            end
            hold off
            leg1 = legend(h1(c_bool),upper(canal_let(c_bool)),'NumColumns',sum(c_bool),'Location','northwest');
            leg1.ItemTokenSize(1) = 5;
            title(leg1,'Axis')
            axes(ha(3))
            hold on
            for i = 1:enum
                h2(i) = plot(NaN,NaN,['k',all_markers(i)]);
            end
            hold off
            leg2 = legend(h2(e_bool),conds(e_bool),'NumColumns',1,'Location','northwest');
            leg2.ItemTokenSize(1) = 5;
            title(leg2,'Conditions')
            savefig([Path,filesep,'Figures',filesep,strrep(fig_name,' ','-'),'.fig'])
            close;
        end
    end  
    if any(sum(exps,2)>1) %Amplitude Sweep
        rel_freqs = freqs(sum(exps,2)>1);      
        for f = 1:length(rel_freqs)
            fig_name = [common_cond,' ',num2str(rel_freqs(f)),'Hz Sine Amplitude Sweep'];
            figure('Units','inch','Position',[2 2 7 3],'Color',[1,1,1]);
            annotation('textbox',[0 0 1 1],'String',[Path,newline,code_Path,filesep,...
            'plotParamResults.m',newline,...
            'VOGA',version,newline,Experimenter],'FontSize',5,...
            'EdgeColor','none','interpreter','none','Color',[0.5,0.5,0.5]);
            annotation('textbox',[0 .9 1 .1],'String',fig_name,'FontSize',14,...
                    'HorizontalAlignment','center','EdgeColor','none');
            ha = gobjects(1,3);
            h1 = gobjects(1,cnum);
            h2 = gobjects(1,enum);
            ha(1) = subplot(1,3,1);
            ha(1).Position = [x_pos(1),y_pos,wid_x,wid_y];
            ha(2) = subplot(1,3,2);
            ha(2).Position = [x_pos(2),y_pos,wid_x,wid_y];
            ha(3) = subplot(1,3,3);
            ha(3).Position = [x_pos(3),y_pos,wid_x,wid_y];
            c_bool = false(1,cnum);
            e_bool = false(1,enum);
            for c = 1:cnum
                for i = 1:enum
                    rel_i = find(all_results.('Frequency(Hz)')==rel_freqs(f)&contains(all_results.AxisName,canals{c})&IC==i);
                    if length(rel_i)>1
                        [a_ax,a_i] = sort(all_results.('Amplitude(dps)')(rel_i));
                        rel_i = rel_i(a_i);
                        axes(ha(1))
                        hold on
                        errorbar(a_ax,maxvel_poscyc(rel_i),maxvel_poscyc_sd(rel_i),[all_markers(i),'-'],'Color',colors.(['l_',canal_let{c}]));
                        hold off
                        axes(ha(2))
                        hold on
                        errorbar(a_ax,maxvel_negcyc(rel_i),maxvel_negcyc_sd(rel_i),[all_markers(i),'-'],'Color',colors.(['l_',canal_let{c}]));
                        hold off
                        axes(ha(3))
                        hold on
                        errorbar(a_ax,phase(rel_i),phase_sd(rel_i),[all_markers(i),'-'],'Color',colors.(['l_',canal_let{c}]));
                        hold off
                        c_bool(c) = true;
                        e_bool(i) = true;
                    end
                end
            end            
            linkaxes(ha,'x')
            linkaxes(ha(1:2),'y')
            set(ha,'XTick',amps,'XTickLabelRotation',0,'XMinorTick','off','XLim',[0 max(amps)*1.05])
            set(ha(1:2),'YLim',[0 YMax])
            set(ha(2),'YTickLabel',[])
            set(ha(3),'YAxisLocation','right','YLim',[-180 360])
            title(ha(1),'Positive Half-Cycle')
            title(ha(2),'Negative Half-Cycle')
            title(ha(3),'Phase Lead')
            ylabel(ha(1),'Maxmimum Eye Velocity (dps)')
            ylabel(ha(3),'Phase (deg)')
            xlabel(ha(2),'Head Velocity (dps)')
            XTickLab = get(ha(1),'XTickLabel');
            XTickLab(amps<100) = {''};
            set(ha,'XTickLabel',XTickLab)
            %Make legends
            axes(ha(2))
            hold on
            for c = 1:cnum
                h1(c) = plot(NaN,NaN,'LineWidth',2,'Color',colors.(['l_',canal_let{c}]));
            end
            hold off
            leg1 = legend(h1(c_bool),upper(canal_let(c_bool)),'NumColumns',sum(c_bool),'Location','northwest');
            leg1.ItemTokenSize(1) = 5;
            title(leg1,'Axis')
            axes(ha(3))
            hold on
            for i = 1:enum
                h2(i) = plot(NaN,NaN,['k',all_markers(i)]);
            end
            hold off
            leg2 = legend(h2(e_bool),conds(e_bool),'NumColumns',1,'Location','northwest');
            leg2.ItemTokenSize(1) = 5;
            title(leg2,'Conditions')
            savefig([Path,filesep,'Figures',filesep,strrep(fig_name,' ','-'),'.fig'])
            close;
            %% Just the canals
            figure('Units','inch','Position',[2 2 7 3],'Color',[1,1,1]);
            annotation('textbox',[0 0 1 1],'String',[Path,newline,code_Path,filesep,...
            'plotParamResults.m',newline,...
            'VOGA',version,newline,Experimenter],'FontSize',5,...
            'EdgeColor','none','interpreter','none','Color',[0.5,0.5,0.5]);
            annotation('textbox',[0 .9 1 .1],'String',fig_name,'FontSize',14,...
                    'HorizontalAlignment','center','EdgeColor','none');
            ha = gobjects(1,3);
            h1 = gobjects(1,cnum);
            h2 = gobjects(1,enum);
            ha(1) = subplot(1,3,1);
            ha(1).Position = [x_pos(1),y_pos,wid_x,wid_y];
            ha(2) = subplot(1,3,2);
            ha(2).Position = [x_pos(2),y_pos,wid_x,wid_y];
            ha(3) = subplot(1,3,3);
            ha(3).Position = [x_pos(3),y_pos,wid_x,wid_y];
            e_bool = false(1,enum);
            for i = 1:enum
                rel_i = find(all_results.('Frequency(Hz)')==rel_freqs(f)&contains(all_results.AxisName,{'LHRH','LARP','RALP'})&IC==i);
                if length(rel_i)>1
                    [a_ax,a_i] = sort(all_results.('Amplitude(dps)')(rel_i));
                    rel_i = rel_i(a_i);
                    ind_L = contains(all_results.AxisName(rel_i),'LARP');
                    ind_R = contains(all_results.AxisName(rel_i),'RALP');
                    ind_Z = contains(all_results.AxisName(rel_i),'LHRH');
                    sub_ind = false(1,length(Subs));
                    for s = 1:length(Subs)
                        sub_ind(s) = contains([fig_name,' ',conds{i}],Subs{s});
                    end
                    ear = Ears{sub_ind};
                    axes(ha(2))
                    hold on
                    errorbar(a_ax(ind_L),phase(rel_i(ind_L)),phase_sd(rel_i(ind_L)),[all_markers(i),'-'],'Color',colors.l_l);
                    errorbar(a_ax(ind_R),phase(rel_i(ind_R)),phase_sd(rel_i(ind_R)),[all_markers(i),'-'],'Color',colors.l_r);
                    errorbar(a_ax(ind_Z),phase(rel_i(ind_Z)),phase_sd(rel_i(ind_Z)),[all_markers(i),'-'],'Color',colors.l_z);
                    hold off
                    if strcmp(ear,'L') % +Z, -L and -R
                        axes(ha(1))
                        hold on
                        errorbar(a_ax(ind_L),maxvel_negcyc(rel_i(ind_L)),maxvel_negcyc_sd(rel_i(ind_L)),[all_markers(i),'-'],'Color',colors.l_l);
                        errorbar(a_ax(ind_R),maxvel_negcyc(rel_i(ind_R)),maxvel_negcyc_sd(rel_i(ind_R)),[all_markers(i),'-'],'Color',colors.l_r);
                        errorbar(a_ax(ind_Z),maxvel_poscyc(rel_i(ind_Z)),maxvel_poscyc_sd(rel_i(ind_Z)),[all_markers(i),'-'],'Color',colors.l_z);
                        hold off
                        axes(ha(3))
                        hold on
                        errorbar(a_ax(ind_L),align_negcyc(rel_i(ind_L)),align_negcyc_sd(rel_i(ind_L)),[all_markers(i),'-'],'Color',colors.l_l);
                        errorbar(a_ax(ind_R),align_negcyc(rel_i(ind_R)),align_negcyc_sd(rel_i(ind_R)),[all_markers(i),'-'],'Color',colors.l_r);
                        errorbar(a_ax(ind_Z),align_poscyc(rel_i(ind_Z)),align_poscyc_sd(rel_i(ind_Z)),[all_markers(i),'-'],'Color',colors.l_z);
                        hold off
                    elseif strcmp(ear,'R') % -Z, +L and +R
                        axes(ha(1))
                        hold on
                        errorbar(a_ax(ind_L),maxvel_poscyc(rel_i(ind_L)),maxvel_poscyc_sd(rel_i(ind_L)),[all_markers(i),'-'],'Color',colors.l_l);
                        errorbar(a_ax(ind_R),maxvel_poscyc(rel_i(ind_R)),maxvel_poscyc_sd(rel_i(ind_R)),[all_markers(i),'-'],'Color',colors.l_r);
                        errorbar(a_ax(ind_Z),maxvel_negcyc(rel_i(ind_Z)),maxvel_negcyc_sd(rel_i(ind_Z)),[all_markers(i),'-'],'Color',colors.l_z);
                        hold off
                        axes(ha(3))
                        hold on
                        errorbar(a_ax(ind_L),align_poscyc(rel_i(ind_L)),align_poscyc_sd(rel_i(ind_L)),[all_markers(i),'-'],'Color',colors.l_l);
                        errorbar(a_ax(ind_R),align_poscyc(rel_i(ind_R)),align_poscyc_sd(rel_i(ind_R)),[all_markers(i),'-'],'Color',colors.l_r);
                        errorbar(a_ax(ind_Z),align_negcyc(rel_i(ind_Z)),align_negcyc_sd(rel_i(ind_Z)),[all_markers(i),'-'],'Color',colors.l_z);
                        hold off
                    end
                    e_bool(i) = true;
                end
            end           
            linkaxes(ha,'x')
            linkaxes(ha(2:3),'y')
            set(ha,'XTick',amps,'XTickLabelRotation',0,'XMinorTick','off','XLim',[0 max(amps)*1.05])
            set(ha(1),'YLim',[0 YMax])
            set(ha(2:3),'YLim',[-180 360])
            set(ha(2),'YTickLabel',[])
            set(ha(2:3),'YAxisLocation','right')
            title(ha(1),'Maximum Eye Velocity')
            title(ha(2),'Phase Lead')
            title(ha(3),'Misalignment')
            ylabel(ha(1),'Velocity (dps)')
            ylabel(ha(3),'Degrees')
            xlabel(ha(2),'Head Velocity (dps)')
            XTickLab = get(ha(1),'XTickLabel');
            XTickLab(amps<100) = {''};
            set(ha,'XTickLabel',XTickLab)
            %Make legends
            h1 = [];
            axes(ha(2))
            hold on
            h1(1) = plot(NaN,NaN,'LineWidth',2,'Color',colors.l_l);
            h1(2) = plot(NaN,NaN,'LineWidth',2,'Color',colors.l_r);
            h1(3) = plot(NaN,NaN,'LineWidth',2,'Color',colors.l_z);
            hold off
            leg1 = legend(h1,{'LARP','RALP','LHRH'},'NumColumns',3,'Location','northwest');
            leg1.ItemTokenSize(1) = 5;
            title(leg1,'Axis')
            axes(ha(3))
            hold on
            for i = 1:enum
                h2(i) = plot(NaN,NaN,['k',all_markers(i)]);
            end
            hold off
            leg2 = legend(h2(e_bool),conds(e_bool),'NumColumns',1,'Location','northwest');
            leg2.ItemTokenSize(1) = 5;
            title(leg2,'Conditions')
            savefig([Path,filesep,'Figures',filesep,strrep(fig_name,' ','-'),'-CanalElectrodesOnly.fig'])
            close;
        end
    end                
elseif(any(contains(all_results.Condition,'Autoscan')))
    %% Make Figure like Boutros 2019 Figure 4 but Magnitude and Misalignment
    %Now figure out which files to plot
    if ~ismember(all_results.Properties.VariableNames,'PulseFreq(pps)') 
        all_exps = all_results.Condition;
        parts = split(all_exps{1},' ');
        exp_name = [all_results.Subject{1},' ',all_results.Visit{1},' ',datestr(all_results.Date(1),'yyyymmdd'),' Autoscan ',parts{contains(parts,'us')},' ',parts{contains(parts,'pps')},' ',all_results.Goggle{1}];
    else
        all_exps = strcat(cellstr(datestr(all_results.Date,'yyyymmdd')),{' '},all_results.Electrode,{' '},cellstr(num2str(all_results.('PhaseDur(us)'))),{'us '},cellstr(num2str(all_results.('PulseFreq(pps)'))),{'pps '},cellstr(num2str(all_results.('CurrentAmp(uA)'))),{'uA '},all_results.Goggle);
        exp_name = [all_results.Subject{1},' ',all_results.Visit{1},' ',datestr(all_results.Date(1),'yyyymmdd'),' Autoscan ',num2str(all_results.('PhaseDur(us)')(1)),'us ',num2str(all_results.('PulseFreq(pps)')(1)),'pps ',all_results.Goggle{1}];
    end
    E_inds = false(length(all_exps),9);
    if isfield(params,'which_files')
        which_files = params.which_files;
    else    
        which_files = questdlg('Plot all the files in this directory or manually select them?','','All','Select','Select');
    end
    for i = 1:9
        E_sub_i = find(contains(all_exps,['E',num2str(i+2)]));
        if ~isempty(E_sub_i)
            if strcmp(which_files,'Select')
                indx = listdlg('ListString',all_exps(E_sub_i),'PromptString',['Pick the files for E',num2str(i+2)],'ListSize',[400 300]);
                E_inds(E_sub_i(indx),i) = true;
            else
                E_inds(E_sub_i,i) = true;
            end
        else
            disp(['No files found for E',num2str(i+2)])
        end
    end
    %Find the cycle numbers for each column
    N = [min(all_results.Cycles(any(E_inds(:,1:3),2))),...
        min(all_results.Cycles(any(E_inds(:,4:6),2))),...
        min(all_results.Cycles(any(E_inds(:,7:9),2)))];
    %Make some bold (if you know which one was activated on)
    E_bold = false(1,9);
    indx = listdlg('ListString',strcat('E',cellfun(@num2str,num2cell(3:11),'UniformOutput',false)),...
        'PromptString','Pick the electrodes to bold. Press Cancel for none.','ListSize',[400 300],'SelectionMode','multiple');
    E_bold(indx) = true;
    %Make the figure
    ha = gobjects(1,6);
    %markerbig=5;
    %markersmall=4;
    %linethick=2;
    %linethin=1;
    errorbarcapsize=1;
    figsizeinches=[7,6];
    XLim = [-5 105];
    YLim_align = [0 80];
    %figsizeinchesBoxplot=[2.3,4];
    figure('Units','inch','Position',[2 2 figsizeinches],'Color',[1,1,1]);%CDS083119a
    if annot
        annotation('textbox',[0 0 1 1],'String',[Path,newline,code_Path,filesep,...
            code_name,newline,...
            'VOGA',version,newline,Experimenter],'FontSize',5,...
            'EdgeColor','none','interpreter','none');
    end
    annotation('textbox',[0 .9 1 .1],'String',strrep(exp_name,'us','\mus'),'FontSize',14,...
        'HorizontalAlignment','center','EdgeColor','none');
    ha(1) = subplot(2,3,1);
    ha(2) = subplot(2,3,2);
    ha(3) = subplot(2,3,3);
    ha(4) = subplot(2,3,4);
    ha(5) = subplot(2,3,5);
    ha(6) = subplot(2,3,6);
    x_min = 0.1;
    x_max = 0.99;
    x_space = 0.01;
    y_min = 0.08;
    y_max = 0.93;
    y_space = 0.03;
    x_wid = (x_max-x_min-2*x_space)/3;
    y_wid = (y_max-y_min-y_space)/2;
    x_pos = x_min:x_wid+x_space:x_max;
    y_pos = y_min:y_wid+y_space:y_max;
    ha(1).Position = [x_pos(1),y_pos(2),x_wid,y_wid];
    ha(2).Position = [x_pos(2),y_pos(2),x_wid,y_wid];
    ha(3).Position = [x_pos(3),y_pos(2),x_wid,y_wid];
    ha(4).Position = [x_pos(1),y_pos(1),x_wid,y_wid];
    ha(5).Position = [x_pos(2),y_pos(1),x_wid,y_wid];
    ha(6).Position = [x_pos(3),y_pos(1),x_wid,y_wid];
    markers = {'x','o','d'};
    %Set colors (faded vs normal)
    if any(contains(all_exps,{'LP';'LH';'LA'})) %Left sided
        color_s = [repmat(colors.l_r_s,3,1);repmat(colors.l_z_s,3,1);repmat(colors.l_l_s,3,1)];
        color = [repmat(colors.l_r,3,1);repmat(colors.l_z,3,1);repmat(colors.l_l,3,1)];
        color(~E_bold,:) = color_s(~E_bold,:);
    else %Right Sided
        color_s = [repmat(colors.l_l_s,3,1);repmat(colors.l_z_s,3,1);repmat(colors.l_r_s,3,1)];
        color = [repmat(colors.l_l,3,1);repmat(colors.l_z,3,1);repmat(colors.l_r,3,1)];
        color(~E_bold,:) = color_s(~E_bold,:);
    end
    %Plot each canal
    for i = 1:3
        rel_tab.E1 = all_results(E_inds(:,3*i-2),:);
        rel_tab.E2 = all_results(E_inds(:,3*i-1),:);
        rel_tab.E3 = all_results(E_inds(:,3*i),:);
        rel_bold = E_bold(3*i-2:3*i);
        i_ord = [find(~rel_bold),find(rel_bold)]; %order of plotting so bold in front
        h = gobjects(1,length(markers));
        labs = cell(1,3);
        %Plot legend
        axes(ha(i))
        hold on
        for j = 1:3
            %Make the fake plots for the labels first
            h(j) = plot(NaN,NaN,'Marker',markers{j},'Color',color(3*i-3+j,:),'LineWidth',1);
            %Make the labels
            if ~isempty(rel_tab.(['E',num2str(j)])) %Something in this one
                %Old way
                %exp = strsplit(rel_tab.(['E',num2str(j)]).Condition{1});
                %labs{1,j} = strrep([exp{contains(exp,'E')},', ',exp{contains(exp,'us')},'/phase'],'u','\mu');
                % New way w/ table
                labs{1,j} = [rel_tab.(['E',num2str(j)]).Electrode{1},', ',num2str(rel_tab.(['E',num2str(j)]).('PhaseDur(us)')(1)),'\mus/phase'];
            else %Remove from the plotting order
                i_ord(i_ord==j) = [];
            end
        end
        h(cellfun(@isempty,labs)) = []; %Take out an electrode with no files
        labs(cellfun(@isempty,labs)) = [];
        hold off
        %Actually plot when values are present
        for j = 1:length(i_ord)
% OLD WAY            
%             exps = rel_tab.(['E',num2str(i_ord(j))]).Condition;
%             curr = NaN(1,length(exps));
%             for q = 1:length(curr)
%                 exp = strsplit(exps{q});
%                 curr(q) = str2double(strrep(exp{contains(exp,'uA')},'uA',''));
%             end
            % New way with table
            exp = rel_tab.(['E',num2str(i_ord(j))]).Electrode(1);
            curr = rel_tab.(['E',num2str(i_ord(j))]).('CurrentAmp(uA)');
            
            [curr_norm,curr_i] = sort(100*(curr-min(curr))/(max(curr)-min(curr)));
            %Determine canal
            if any(contains(exp,{'LP','RA'})) %RALP
                canal = 'R';
            elseif any(contains(exp,{'LH','RH'})) %LHRH
                canal = 'Z';
            else %LARP
                canal = 'L';
            end
            %Extract the relevant vectors
            L_Vel = abs(rel_tab.(['E',num2str(i_ord(j))]).(['MaxVel_L',canal,'_HIGH'])(curr_i));
            L_Vel_sd = rel_tab.(['E',num2str(i_ord(j))]).(['MaxVel_L',canal,'_HIGH_sd'])(curr_i);
            R_Vel = abs(rel_tab.(['E',num2str(i_ord(j))]).(['MaxVel_R',canal,'_HIGH'])(curr_i));
            R_Vel_sd = rel_tab.(['E',num2str(i_ord(j))]).(['MaxVel_R',canal,'_HIGH_sd'])(curr_i);
            L_Align = rel_tab.(['E',num2str(i_ord(j))]).Align_L_HIGH(curr_i);
            L_Align_sd = rel_tab.(['E',num2str(i_ord(j))]).Align_L_HIGH_sd(curr_i);
            R_Align = rel_tab.(['E',num2str(i_ord(j))]).Align_R_HIGH(curr_i);
            R_Align_sd = rel_tab.(['E',num2str(i_ord(j))]).Align_R_HIGH_sd(curr_i);
            %Choose the larger eye for each current
            [~,eye] = max([L_Vel,R_Vel],[],2);
            Vel = L_Vel;
            Vel(eye==2) = R_Vel(eye==2);
            Vel_sd = L_Vel_sd;
            Vel_sd(eye==2) = R_Vel_sd(eye==2);
            Align = L_Align;
            Align(eye==2) = R_Align(eye==2);
            Align_sd = L_Align_sd;
            Align_sd(eye==2) = R_Align_sd(eye==2);
            %Plot velocity
            axes(ha(i))
            hold on
            errorbar(curr_norm,Vel,Vel_sd,'Color',color(3*i-3+i_ord(j),:),'LineStyle','none','LineWidth',1,'CapSize',errorbarcapsize)
            plot(curr_norm,Vel,'Marker',markers{i_ord(j)},'Color',color(3*i-3+i_ord(j),:),'LineWidth',1)
            hold off
            %Plot alignment
            axes(ha(i+3))
            hold on
            errorbar(curr_norm,Align,Align_sd,'Color',color(3*i-3+i_ord(j),:),'LineStyle','none','LineWidth',1,'CapSize',errorbarcapsize)
            plot(curr_norm,Align,'Marker',markers{i_ord(j)},'Color',color(3*i-3+i_ord(j),:),'LineWidth',1)
            hold off
        end
        axes(ha(i))
        leg = legend(h,labs,'box','off','Location','northwest','FontSize',7);
        leg.ItemTokenSize(1) = 12;
        set(gca,'box','off')
        if i == 1
            ylabel({'Eye Velocity';'Magnitude (\circ/s)'})
        else
            set(gca,'YColor','none')
        end
        set(gca,'XLim',XLim,'XColor','none')
        axes(ha(i+3))
        set(gca,'YLim',YLim_align)
        if i == 1
            set(gca,'YTick',10:10:max(YLim_align))
            ylabel({'Misalignment';'Angle (\circ)'})
        else
            set(gca,'YColor','none')
        end
        set(gca,'XLim',XLim,'XTick',0:20:100)
        if i==2
            xlabel('% of Current Range')
        end
    end
    if isempty(YMax)
       YMax = max(reshape(cell2mat(get(ha(1:3),'Ylim')),[],1));
    end    
    set(ha(1:3),'YLim',[0,YMax])
    set(ha(1),'YTick',20:20:YMax)    
    fig_name = inputdlg('Name this figure','',1,{[strrep(exp_name,' ','-'),'.fig']});
    if ~isempty(fig_name)
        savefig([Path,filesep,'Figures',filesep,fig_name{:}])
    end
    close;
elseif(any(contains(all_results.Type,'Impulse')))
    %% Impulse
    all_canals = {'LA','LP','LH','RP','RA','RH'}; %Preferred order
    all_results = all_results(contains(all_results.Type,'Impulse'),:);
    fn = size(all_results,1);
    file_parts = [all_results.Subject,all_results.Visit,cellstr(datestr(all_results.Date,'yyyymmdd')),...
        all_results.Experiment,all_results.Condition,all_results.Goggle];
    rel_file_parts = file_parts; %everything but freq/amp
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
        conds = all_results.Condition(1);
        IC = ones(fn,1);
    elseif size(rel_file_parts,2)==1
        [conds,~,IC] = unique(rel_file_parts,'stable');
    else
        [conds,~,IC] = unique(join(rel_file_parts),'stable');
    end
    enum = length(conds);
    %Organize the gains and latencies
    gain = NaN(enum,length(all_canals));
    gain_sd = NaN(enum,length(all_canals));
    lat = NaN(enum,length(all_canals));
    lat_sd = NaN(enum,length(all_canals));
    for i = 1:fn
        if all(all_results.StimAxis{i}==[0,0,1])
             rel_col = 'Z';
             c = 3;
        elseif all(all_results.StimAxis{i}==[0,0,-1])
             rel_col = 'Z';
             c = 6;
        elseif all(all_results.StimAxis{i}==[1,0,0])
             rel_col = 'L';
             c = 4;
        elseif all(all_results.StimAxis{i}==[-1,0,0])
             rel_col = 'L';
             c = 1;
        elseif all(all_results.StimAxis{i}==[0,-1,0])
             rel_col = 'R';
             c = 2;
        elseif all(all_results.StimAxis{i}==[0,1,0])
             rel_col = 'R';
             c = 5;
        end
        [~,eye] = max(abs([all_results.(['Gain_L',rel_col,'_HIGH'])(i),all_results.(['Gain_R',rel_col,'_HIGH'])(i)]));
        if mod(eye,2)==1 %Left
            eye_s = 'L';
        else %Right
            eye_s = 'R';
        end
        gain(IC(i),c) = all_results.(['Gain_',eye_s,rel_col,'_HIGH'])(i);
        gain_sd(IC(i),c) = all_results.(['Gain_',eye_s,rel_col,'_HIGH_sd'])(i);
        lat(IC(i),c) = all_results.(['Latency_',eye_s,rel_col,])(i);
        lat_sd(IC(i),c) = all_results.(['Latency_',eye_s,rel_col,'_sd'])(i);
    end            
    %Position information for the graph types
    max_y = 0.80;
    min_y = 0.15;
    %spac_y = 0.01;
    min_x = 0.08;
    max_x = 0.90;
    spac_x = 0.01;
    wid_x = (max_x-min_x-spac_x)/2;
    x_pos = min_x:(wid_x+spac_x):max_x-wid_x; 
    wid_y = max_y-min_y;
    y_pos = min_y;    
    all_markers = 'o*.sd^v><ph_|';
    fig_name = [common_cond,' Head Impulse Test'];
    figure('Units','inch','Position',[2 2 7 3],'Color',[1,1,1]);
    annotation('textbox',[0 0 1 1],'String',[Path,newline,code_Path,filesep,...
    'plotParamResults.m',newline,...
    'VOGA',version,newline,Experimenter],'FontSize',5,...
    'EdgeColor','none','interpreter','none','Color',[0.5,0.5,0.5]);
    annotation('textbox',[0 .9 1 .1],'String',fig_name,'FontSize',14,...
            'HorizontalAlignment','center','EdgeColor','none');
    ha = gobjects(1,2);
    h1 = gobjects(1,enum);
    ha(1) = subplot(1,2,1);
    ha(1).Position = [x_pos(1),y_pos,wid_x,wid_y];
    ha(2) = subplot(1,2,2);
    ha(2).Position = [x_pos(2),y_pos,wid_x,wid_y];
    for i = 1:enum
        val = ~isnan(gain(i,:));
        x_ax = 1:6;
        x_ax(~val) = []; 
        axes(ha(1))
        hold on
        errorbar(x_ax,gain(i,val),gain_sd(i,val),[all_markers(i),'-'],'Color','k');
        hold off
        axes(ha(2))
        hold on
        errorbar(x_ax,lat(i,val),lat_sd(i,val),[all_markers(i),'-'],'Color','k');
        hold off
    end           
    linkaxes(ha,'x')
    set(ha(2),'YAxisLocation','right')
    title(ha(1),'Gain')
    title(ha(2),'Latency (ms)')
    ylabel(ha(1),'Eye Position/Head Position')
    ylabel(ha(2),'Time Between Head and Eye Motion')
    xlabel(ha(1),'Canal')
    xlabel(ha(2),'Canal')
    set(ha,'XTick',1:6,'XTickLabel',all_canals,'XTickLabelRotation',0,'XMinorTick','off')
    set(ha,'XLim',[0.5 6.5])
    ymax = get(ha(1),'YLim');
    set(ha(1),'YLim',[0 max(ymax(2),1.1)])
    set(ha(2),'YLim',[-7 200])
    %Make legends
    axes(ha(2))
    hold on
    for i = 1:enum
        h1(i) = plot(NaN,NaN,['k',all_markers(i)]);
    end
    hold off
    leg = legend(h1,conds,'NumColumns',1,'Location','northwest');
    leg.ItemTokenSize(1) = 5;
    title(leg,'Conditions')
    savefig([Path,filesep,'Figures',filesep,strrep(fig_name,' ','-'),'.fig'])   
end
end