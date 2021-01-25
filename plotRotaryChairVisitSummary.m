function plotRotaryChairVisitSummary(path,Cyc_Path,code_Path)
if nargin == 0
    code_Path = '/Volumes/MVI/DATA SUMMARY/IN PROGRESS/VOG Analysis Scripts/LDVOG_Neurolign';
    path = cd;
    Cyc_Path = [path,filesep,'Cycle Averages'];
end
%% Colors
% Normal colors
colors.l_x = [237,150,33]/255;
colors.l_y = [125,46,143]/255;
colors.l_z = [1 0 0];
colors.l_l = [0,128,0]/255;
colors.l_r = [0 0 1];
colors.r_x = [237,204,33]/255;
colors.r_y = [125,46,230]/255;
colors.r_z = [1,0,1];
colors.r_l = [0 1 0];
colors.r_r = [64,224,208]/255;
% Faded colors
colors.l_x_s = colors.l_x + 0.5*(1-colors.l_x);
colors.l_y_s = colors.l_y + 0.5*(1-colors.l_y);
colors.l_z_s = colors.l_z + 0.5*(1-colors.l_z);
colors.l_l_s = colors.l_l + 0.5*(1-colors.l_l);
colors.l_r_s = colors.l_r + 0.5*(1-colors.l_r);
colors.r_x_s = colors.r_x + 0.5*(1-colors.r_x);
colors.r_y_s = colors.r_y + 0.5*(1-colors.r_y);
colors.r_z_s = colors.r_z + 0.5*(1-colors.r_z);
colors.r_l_s = colors.r_l + 0.5*(1-colors.r_l);
colors.r_r_s = colors.r_r + 0.5*(1-colors.r_r);
%% Plot all Cycle Averages for an Experiment Type
temp = dir([Cyc_Path,filesep,'Sine*.mat']);
cyc_files = {temp.name}';
cyc_files(contains(cyc_files,'NotAnalyzeable')) = [];
file_parts = cell(length(cyc_files),2);
for i = 1:length(cyc_files)
    fname = strrep(strrep(cyc_files{i},'CycAvg_',''),'.mat','');
    fparts = split(fname,'-');
    file_parts(i,2) = {strrep(fparts{contains(fparts,'Hz')},'p','.')};
    fparts(contains(fparts,'Hz')) = [];
    fparts(contains(fparts,'dps')) = [];
    file_parts(i,1) = {strjoin(fparts,' ')};
end
exp_name = unique(file_parts(:,1));
enum = length(exp_name);
[~,indf] = sort(cellfun(@str2double,strrep(unique(file_parts(:,2)),'Hz',''))); %make sure it's in numerical order
freqs = unique(file_parts(:,2));
freqs = freqs(indf);
fnum = length(freqs);
for j = 1:enum
    YLim = [-100 100];
    figure('Units','inches','Position',[0.2778    5.8472   17.2222    3.8333],'Color',[1,1,1])
    %Title
    annotation('textbox',[0 .9 1 .1],'String',exp_name{j},'FontSize',14,...
    'HorizontalAlignment','center','EdgeColor','none');
    annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,'plotRotaryChairVisitSummary.m'],'FontSize',5,...
    'EdgeColor','none','interpreter','none');
    ha = gobjects(1,fnum);
    x_space = 0.01;
    x_min = 0.04;
    x_max = 0.98;
    x_wid = (x_max-x_min-x_space*(fnum-1))/fnum;
    y_height = 0.75;
    x_pos = x_min:(x_wid+x_space):x_max;
    y_pos = 0.12;
    for i = 1:fnum
        ha(i) = subplot(1000,1000,1);
        ha(i).Position = [x_pos(i) y_pos x_wid y_height];
        if sum(contains(cyc_files,strrep(exp_name{j},' ','-'))&contains(cyc_files,['-',freqs{i}]))==1
            load([Cyc_Path,filesep,cyc_files{contains(cyc_files,strrep(exp_name{j},' ','-'))&contains(cyc_files,['-',freqs{i}])}],'CycAvg')
            fields = fieldnames(CycAvg);       
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
            h(1) = plot(CycAvg.t(s),CycAvg.stim(s),'k');
            hold on
            %Now add the fills and standard deviations and means
            %LE-LHRH
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s)))],colors.l_z_s)
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s)))],colors.l_z_s)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
            h(2) = plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
            %RE-LHRH
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s)))],colors.r_z_s)
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s)))],colors.r_z_s)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
            h(3) = plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
            %LE-LARP
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ll_cycavg(s),fliplr((CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s)))],colors.l_l_s)
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ll_cycavg(s),fliplr((CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s)))],colors.l_l_s)
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s),'Color',colors.l_l)
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s),'Color',colors.l_l)
            h(4) = plot(CycAvg.t(s),CycAvg.ll_cycavg(s),'Color',colors.l_l,'LineWidth',2);
            %RE-LARP
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rl_cycavg(s),fliplr((CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s)))],colors.r_l_s)
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rl_cycavg(s),fliplr((CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s)))],colors.r_l_s)
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s),'Color',colors.r_l)
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s),'Color',colors.r_l)
            h(5) = plot(CycAvg.t(s),CycAvg.rl_cycavg(s),'Color',colors.r_l,'LineWidth',2);
            %LE_RALP
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lr_cycavg(s),fliplr((CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s)))],colors.l_r_s)
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lr_cycavg(s),fliplr((CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s)))],colors.l_r_s)
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s),'Color',colors.l_r)
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s),'Color',colors.l_r)
            h(6) = plot(CycAvg.t(s),CycAvg.lr_cycavg(s),'Color',colors.l_r,'LineWidth',2);
            %RE-RALP
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rr_cycavg(s),fliplr((CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s)))],colors.r_r_s)
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rr_cycavg(s),fliplr((CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s)))],colors.r_r_s)
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s),'Color',colors.r_r)
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s),'Color',colors.r_r)
            h(7) = plot(CycAvg.t(s),CycAvg.rr_cycavg(s),'Color',colors.r_r,'LineWidth',2);
            hold off
        else 
            plot(NaN,NaN)
        end
        set(gca,'XLim',[0 1/str2double(strrep(freqs{i},'Hz',''))])
        set(gca,'YLim',YLim)
        title(freqs{i})
        xlabel('Time (s)')
        if i == 1
            ylabel('Angular Velocity (dps)')            
        else
            set(gca,'YTickLabel',[])
        end
    end  
    legend(ha(1),h,'Head Z','Left Z','Right Z','Left LARP','Right LARP','Left RALP','Right RALP')
    savefig([path,filesep,strrep(exp_name{j},' ','-'),'.fig'])
end
%% Plot stimulation half-cycle gain over frequency for each of the conditions
% Normative Rotary Chair Data
freq = [0.005 0.010 0.020 0.050 0.100 0.200 0.500 1.000];
%Normative Data from Wall et. al. 1984 on patients 50-69
norm_gain_m = [0.2175	0.3576	0.4803	0.5634	0.56	0.5445	0.5678	0.698];
norm_gain_std = [0.0757	0.1235	0.1474	0.2243	0.2363	0.2509	0.2121	0.2777];
norm_phase_m = [70.9478	43.6413	25.9262	10.7906	4.084	-4.8461	-11.4084 -11.3412];
norm_phase_std = [7.5044	6.526	6.714	8.0519	4.5136	7.6279	6.2448	11.3167];
%This can only be done if there is a RotaryChairResults.mat file
%See if one already exists
temp = dir([path,filesep,'*RotaryChairResults.mat']);
res_file = {temp.name}';
if isempty(res_file)
    disp('First, make the table with the cycle sine fits...')
    MakeCycleSineFitTable(path,Cyc_Path);
    temp = dir([path,filesep,'*RotaryChairResults.mat']);
    res_file = {temp.name}';
end
load(res_file{end},'all_results')
%Get the stim ear
load(all_results{1,1}{:},'CycAvg')
ear = CycAvg.info.ear;
%Figure out which frequenicies you need
[freq_ax,indf] = sort(cellfun(@str2double,strrep(unique(all_results.Frequency)','Hz','')));
freqs = unique(all_results.Frequency);
freqs = freqs(indf);
fnum = length(freqs);
%Figure out which experiments to plot
[dateconds,id] = unique(strcat(cellstr(datestr(all_results.Date)),{' '},all_results.Condition),'stable');
conds = unique(all_results.Condition,'stable');
enum = max([length(dateconds),length(conds)]);
if length(dateconds)==length(conds) %Don't need the date
    exp_name = conds';
else
    exp_name = dateconds';
end
%Create arrays for graphing
excludedatawithgainbelow=0.015;
%Make the items to plot
gain = NaN(enum,fnum);
gain_sd = NaN(enum,fnum);
phase = NaN(enum,fnum);
phase_sd = NaN(enum,fnum);
for i = 1:enum
    inds = find(all_results.Date==all_results.Date(id(i))&contains(all_results.Condition,all_results.Condition(id(i))));
    [~,f_ind] = ismember(all_results.Frequency(inds),freqs);
    if strcmp(ear,'L')
        [sub_gain,eye] = max([all_results.LHRH_L_POS_Gain(inds)';all_results.LHRH_R_POS_Gain(inds)']);
        sub_gain_sd = all_results.LHRH_L_POS_Gain_sd(inds)';
        sub_gain_sd(eye==2) = all_results.LHRH_R_POS_Gain_sd(inds(eye==2));
    elseif strcmp(ear,'R')
        [sub_gain,eye] = max([all_results.LHRH_L_NEG_Gain(inds)';all_results.LHRH_R_NEG_Gain(inds)']);
        sub_gain_sd = all_results.LHRH_L_NEG_Gain_sd(inds)';
        sub_gain_sd(eye==2) = all_results.LHRH_R_NEG_Gain_sd(inds(eye==2));
    end
    sub_phase = all_results.L_Phase(inds);
    sub_phase(eye==2) = all_results.R_Phase(inds(eye==2));
    sub_phase_sd = all_results.L_Phase_sd(inds);
    sub_phase_sd(eye==2) = all_results.R_Phase_sd(inds(eye==2));
    sub_phase(sub_gain < excludedatawithgainbelow) = NaN;
    sub_phase_sd(sub_gain < excludedatawithgainbelow) = NaN;
    %remove SD for n=1
    n = all_results.n_cycles(inds);
    sub_gain_sd(n==1) = NaN;
    sub_phase_sd(n==1) = NaN;
    %Put into mats
    gain(i,f_ind) = sub_gain;
    gain_sd(i,f_ind) = sub_gain_sd;
    phase(i,f_ind) = sub_phase;
    phase_sd(i,f_ind) = sub_phase_sd;
end
% Plot
h=gobjects(1,enum+3);
ha = gobjects(1,2);
logxshift=1.03; %how much to multiply the x value to offset its marker rightward (divide to move left)
%markerbig=5;
%markersmall=4;
linethick=2;
linethin=1;
errorbarcapsize=1;
graydark=0.85;
graylight=0.95;
figsizeinches=[7,6];
%figsizeinchesBoxplot=[2.3,4];
plot_offset = [logxshift^-1, logxshift^-2, 1, logxshift, logxshift^2];
figure('Units','inch','Position',[2 2 figsizeinches],'Color',[1,1,1]);%CDS083119a
annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,'plotRotaryChairVisitSummary.m'],'FontSize',5,...
    'EdgeColor','none','interpreter','none','VerticalAlignment','bottom');
ha(1) = subplot(2,1,1);
ha(1).Position = [0.1,0.55,0.85,0.4];
ha(2) = subplot(2,1,2);
ha(2).Position = [0.1,0.10,0.85,0.4];
axes(ha(1))
hold on
h(enum+3) = fill([freq,fliplr(freq)],[norm_gain_m+2*norm_gain_std,fliplr(norm_gain_m-2*norm_gain_std)],graylight*[1 1 1],'LineStyle','none');
h(enum+2) = fill([freq,fliplr(freq)],[norm_gain_m+1*norm_gain_std,fliplr(norm_gain_m-1*norm_gain_std)],graydark*[1 1 1],'LineStyle','none');
h(enum+1) = plot(freq,norm_gain_m,'k--','LineWidth',linethick);
for i = 1:enum
    h(i) = plot(plot_offset(i)*freq_ax,gain(i,:),'LineWidth',linethick);
    errorbar(plot_offset(i)*freq_ax,gain(i,:),gain_sd(i,:),'Color',h(i).Color,'LineStyle','none','LineWidth',linethin,'CapSize',errorbarcapsize)
end
hold off
legend(h,[exp_name,{'Normal mean','Normal±1SD','Normal±2SD'}],'Location','NorthWest','NumColumns',2)
title([all_results.Subject{1},' ',all_results.Visit{1},' Rotary Chair Sinusoids'])
ylabel('Horizontal VOR Gain')
set(ha(1),'XTick',freq_ax,'XTickLabel',[])
axes(ha(2))
hold on
fill([freq,fliplr(freq)],[norm_phase_m+2*norm_phase_std,fliplr(norm_phase_m-2*norm_phase_std)],graylight*[1 1 1],'LineStyle','none')
fill([freq,fliplr(freq)],[norm_phase_m+1*norm_phase_std,fliplr(norm_phase_m-1*norm_phase_std)],graydark*[1 1 1],'LineStyle','none')
plot(freq,norm_phase_m,'k--','LineWidth',linethick)
for i = 1:enum
    plot(plot_offset(i)*freq_ax,phase(i,:),'Color',h(i).Color,'LineWidth',linethick);
    errorbar(plot_offset(i)*freq_ax,phase(i,:),phase_sd(i,:),'Color',h(i).Color,'LineStyle','none','LineWidth',linethin,'CapSize',errorbarcapsize)
end
%plot([0.01,3],[0,0],'k:')
hold off
ylabel('Phase Lead (deg)')
set(ha(2),'XTick',freq_ax,'XTickLabel',freq_ax)
xlabel('Frequency [Hz]')
set(ha,'XLim',[0.9*freq_ax(1),2.2],'YTickLabelMode','auto','XScale','log')
set(ha(1),'YLim',[0 0.85],'YTick',0.1:0.1:0.8)
set(ha(2),'YLim',[-30,100],'YTick',-20:20:100)
savefig([path,filesep,all_results.Subject{1},'-',all_results.Visit{1},'-RotaryChair-Sine-LHRH-OverFreq','.fig'])
end