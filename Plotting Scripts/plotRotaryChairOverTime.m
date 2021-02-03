%% Plot Rotary Chair Over Time
%Plots Gain (pos and neg half-cycle), phase, and alignment (post and neg
%half-cycle) for one subject and multiple conditions over time by frequency
clc;
clear;
close all;
%% Normative Rotary Chair Data
freq = [0.005 0.010 0.020 0.050 0.100 0.200 0.500 1.000];
%Normative Data from Wall et. al. 1984 on patients 50-69
norm_gain_m = [0.2175	0.3576	0.4803	0.5634	0.56	0.5445	0.5678	0.698];
norm_gain_std = [0.0757	0.1235	0.1474	0.2243	0.2363	0.2509	0.2121	0.2777];
norm_phase_m = [70.9478	43.6413	25.9262	10.7906	4.084	-4.8461	-11.4084 -11.3412];
norm_phase_std = [7.5044	6.526	6.714	8.0519	4.5136	7.6279	6.2448	11.3167];
%% Subject Information
% Add to as needed - current as of subject 8 and all visits pre-COVID (May
% 2020). Any devices that have been mailed to patients are listed as "next"
% because they don't need to be in these figures.
Ears = {'L','L','L','L','R','R','L','R','L'};
%All parameter changes 
Param_Changes(1) = {datetime([2016,09,07;...
                              2016,09,12;...
                              2016,10,19;...
                              2017,01,18])};     %next 2020,05,09  
Param_Changes(2) = {datetime([2016,11,30;...
                              2016,12,02;...
                              2017,1,11])};     
Param_Changes(3) = {datetime([2017,2,23;...
                              2017,3,08;...
                              2019,3,27])};                           
Param_Changes(4) = {datetime([2018,1,5;...
                              2018,1,9;...
                              2018,16,2])};    %next 2020,04,04
Param_Changes(5) = {datetime([2018,9,12;...
                              2020,3,31])};     %next 2020,05,02                           
Param_Changes(6) = {datetime([2018,9,26;...
                              2020,02,28])}; 
Param_Changes(7) = {datetime([2019,2,6;...
                              2019,7,2])};     %next 2020,05,03                      
Param_Changes(8) = {datetime([2019,10,09;...
                              2019,11,01;...
                              2019,11,27;...
                              2020,1,16;...
                              2020,3,6])};   
Param_Changes(9) = {datetime([2020,10,15])};                          
%% Get information of what to plot
[path2,path1] = uigetfile('*.mat','Select Parameter Fit File');
load([path1,filesep,path2],'all_results')
%Subject Information
subjects = unique(all_results.Subject);
sub_num = nmlistdlg('PromptString','Subject:','SelectionMode','single','ListString',subjects);
subject = subjects(sub_num);
tab = sortrows(all_results(contains(all_results.Subject,subject),:),'Date','ascend');
sub = str2double(subject{:}(6));

%Candidate
% subject = {'R908'};
% sub = 9;
% tab = sortrows(all_results,'Date','ascend');
ear = Ears{sub};
param_dates = Param_Changes{sub};
visits = unique(tab.Visit,'stable');
freqs = {'0.01Hz','0.02Hz','0.05Hz','0.1Hz','0.2Hz','0.5Hz','1Hz','2Hz'};
%% Plot Gain (one subject, ipsi half-cycle)
%Plot values
XLim = [4 1000];
xtick = [5,28*[0,1,6,24]+10];
xticklab = [-1,0,1,6,24];
YLim_lowf = [-0.01,0.6];
YLim_highf = [-0.025,1.0];
%Stimulation Half-Cycle
figure(1)
colors = get(gca,'ColorOrder');
colors = [colors(3:end,:);colors(1:2,:)];
ha = gobjects(length(freqs),1);
h = [];
subplot(1,1,1)
delete(findall(gcf,'type','annotation'))
annotation('textbox',[0 .9 1 .1],'String',[subject{:},' Horizontal Rotary Chair: Implant Stimulation Half-Cycle'],'FontSize',12,...
    'HorizontalAlignment','center','EdgeColor','none');
for i = 1:length(freqs)
    %Isolate values
    sub_tab = tab(strcmp(tab.Frequency,freqs{i}),:);
    act_days = days(sub_tab.Date-param_dates(1))+10;
    act_days(act_days < 5) = 5; %Move Pre-Op Value 
    %Make axes
    ha(i) = subplot(2,4,i);
    hold on
    h(1) = plot(NaN,NaN,'rx');
    h(2) = plot(NaN,NaN,'mx');
    h(3) = plot(NaN,NaN,'r*-');
    h(4) = plot(NaN,NaN,'m*-');
    h(5) = plot(NaN,NaN,'ro:');
    h(6) = plot(NaN,NaN,'mo:');
    leg_lab = {'Left Eye - Pre-Activation','Right Eye - Pre-Activation','Left Eye - Modulation ON','Right Eye - Modution ON','Left Eye - Modulation OFF','Right Eye - Modulation OFF','Pre-Operative','Post-Operative/Pre-Activation'};
    set(gca,'XLim',XLim,'XTick',xtick,'xscale','log','FontSize',8,'box','on')
    set(gca,'XGrid',1,'XMinorGrid',0,'YGrid',1,'YMinorGrid',0)
    if i > 4 %bottom row of graphs
        YLim = YLim_highf;
        set(gca,'Ylim',YLim,'XTickLabels',xticklab)
        xlabel('Months Since Activation','FontSize',8)
        set(gca,'Position',[0.06+(0.225*(i-5)),0.05,0.215,0.4])
        bh = 0.1;
    else 
        YLim = YLim_lowf;
        set(gca,'Ylim',YLim,'XTickLabels',[])
        set(gca,'Position',[0.06+0.225*(i-1),0.5,0.215,0.4])
        bh = 0.05;
    end
    set(gca,'YTick',0:0.2:YLim(2),'YTickLabel',[])
    if i == 1 || i ==5
        set(gca,'YTickLabel',0:0.2:YLim(2))
        ylabel('Gain w/re to Head Velocity')
    elseif i == 4 || i == 8
        set(gca,'YTickLabel',100*(0:0.2:YLim(2)))
        set(gca,'YAxisLocation','right')
        ylabel('Eye Velocity (dps)')
    end    
    set(gca,'XMinorTick','off')
    title(freqs{i})
    
    %Plot
    h(7) = fill([4 6 6 4],YLim(2)+[-bh,-bh,0,0],0.75*[1,1,1]);
    xline(6);
    h(8) = fill([6 10 10 6],YLim(2)+[-bh,-bh,0,0],0.85*[1,1,1]);
    xline(10);
    for j = 1:length(param_dates)
        n = mod(j,size(colors,1));
        n(n==0) = size(colors,1);
        if j == length(param_dates)
            time = XLim(2)*ones(1,4);
            time([1,4]) = days(param_dates(j)-param_dates(1))+10;
        else
            time = days([param_dates(j) param_dates(j+1) param_dates(j+1) param_dates(j)]-param_dates(1))+10;
        end
        h(8+j) = fill(time,YLim(2)+[-bh,-bh,0,0],colors(n,:));
        leg_lab{8+j} = ['Mapping #',num2str(j)];
        xline(time(2));
    end    
    if ~isempty(sub_tab)
        i1 = contains(sub_tab.Condition,'Pre-Op'); %Visit 0
        i2 = contains(sub_tab.Condition,'NoStim'); %Visit 3 and anothers
        i3 = contains(sub_tab.Condition,'MotionMod'); 
        i4 = contains(sub_tab.Condition,'ConstantRate');
        switch ear
            case 'L'
                errorbar(act_days(i1),sub_tab.LHRH_L_POS_Gain(i1),sub_tab.LHRH_L_POS_Gain_sd(i1),'rx')
                errorbar(act_days(i1),sub_tab.LHRH_R_POS_Gain(i1),sub_tab.LHRH_R_POS_Gain_sd(i1),'mx')
                errorbar(act_days(i2),sub_tab.LHRH_L_POS_Gain(i2),sub_tab.LHRH_L_POS_Gain_sd(i2),'rx')
                errorbar(act_days(i2),sub_tab.LHRH_R_POS_Gain(i2),sub_tab.LHRH_R_POS_Gain_sd(i2),'mx')
                errorbar(act_days(i3),sub_tab.LHRH_L_POS_Gain(i3),sub_tab.LHRH_L_POS_Gain_sd(i3),'r*-')
                errorbar(act_days(i3),sub_tab.LHRH_R_POS_Gain(i3),sub_tab.LHRH_R_POS_Gain_sd(i3),'m*-')
                errorbar(act_days(i4),sub_tab.LHRH_L_POS_Gain(i4),sub_tab.LHRH_L_POS_Gain_sd(i4),'ro:')
                errorbar(act_days(i4),sub_tab.LHRH_R_POS_Gain(i4),sub_tab.LHRH_R_POS_Gain_sd(i4),'mo:')
            case 'R'
                errorbar(act_days(i1),sub_tab.LHRH_L_NEG_Gain(i1),sub_tab.LHRH_L_NEG_Gain_sd(i1),'rx')
                errorbar(act_days(i1),sub_tab.LHRH_R_NEG_Gain(i1),sub_tab.LHRH_R_NEG_Gain_sd(i1),'mx')
                errorbar(act_days(i2),sub_tab.LHRH_L_NEG_Gain(i2),sub_tab.LHRH_L_NEG_Gain_sd(i2),'rx')
                errorbar(act_days(i2),sub_tab.LHRH_R_NEG_Gain(i2),sub_tab.LHRH_R_NEG_Gain_sd(i2),'mx')
                errorbar(act_days(i3),sub_tab.LHRH_L_NEG_Gain(i3),sub_tab.LHRH_L_NEG_Gain_sd(i3),'r*-')
                errorbar(act_days(i3),sub_tab.LHRH_R_NEG_Gain(i3),sub_tab.LHRH_R_NEG_Gain_sd(i3),'m*-')
                errorbar(act_days(i4),sub_tab.LHRH_L_NEG_Gain(i4),sub_tab.LHRH_L_NEG_Gain_sd(i4),'ro:')
                errorbar(act_days(i4),sub_tab.LHRH_R_NEG_Gain(i4),sub_tab.LHRH_R_NEG_Gain_sd(i4),'mo:')
        end
    end
    hold off
    if i == 1
        legend(h,leg_lab,'Location','northwest','NumColumns',1)
    end
end
ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]))
%savefig([subject{:},'_RotaryChairGain_AcrossTime.fig'])
%% Plot Gain (one subject, contra half-cycle)
%Plot values
XLim = [4 1000];
xtick = [5,28*[0,1,6,24]+10];
xticklab = [-1,0,1,6,24];
YLim_lowf = [-0.01,0.6];
YLim_highf = [-0.025,1.0];
%Inhibition Half-Cycle
figure(2)
colors = get(gca,'ColorOrder');
colors = [colors(3:end,:);colors(1:2,:)];
ha = gobjects(length(freqs),1);
subplot(1,1,1)
delete(findall(gcf,'type','annotation'))
annotation('textbox',[0 .9 1 .1],'String',[subject{:},' Horizontal Rotary Chair: Implant Inhibition Half-Cycle'],'FontSize',12,...
    'HorizontalAlignment','center','EdgeColor','none');
for i = 1:length(freqs)
    %Isolate values
    sub_tab = tab(strcmp(tab.Frequency,freqs{i}),:);
    act_days = days(sub_tab.Date-param_dates(1))+10;
    act_days(act_days < 5) = 5; %Move Pre-Op Value 
    %Make axes
    ha(i) = subplot(2,4,i);
    hold on
    set(gca,'XLim',XLim,'XTick',xtick,'xscale','log','FontSize',8,'box','on')
    set(gca,'XGrid',1,'XMinorGrid',0,'YGrid',1,'YMinorGrid',0)
    if i > 4 %bottom row of graphs
        YLim = YLim_highf;
        set(gca,'Ylim',YLim,'XTickLabels',xticklab)
        xlabel('Months Since Activation','FontSize',8)
        set(gca,'Position',[0.06+(0.225*(i-5)),0.05,0.215,0.4])
        bh = 0.1;
    else 
        YLim = YLim_lowf;
        set(gca,'Ylim',YLim,'XTickLabels',[])
        set(gca,'Position',[0.06+0.225*(i-1),0.5,0.215,0.4])
        bh = 0.05;
    end
    set(gca,'YTick',0:0.2:YLim(2),'YTickLabel',[])
    if i == 1 || i ==5
        set(gca,'YTickLabel',0:0.2:YLim(2))
        ylabel('Gain w/re to Head Velocity')
    elseif i == 4 || i == 8
        set(gca,'YTickLabel',100*(0:0.2:YLim(2)))
        set(gca,'YAxisLocation','right')
        ylabel('Eye Velocity (dps)')
    end    
    set(gca,'XMinorTick','off')
    title(freqs{i})
       
    %Plot
    fill([4 6 6 4],YLim(2)+[-bh,-bh,0,0],0.75*[1,1,1]);
    xline(6);
    fill([6 10 10 6],YLim(2)+[-bh,-bh,0,0],0.85*[1,1,1]);
    xline(10);
    for j = 1:length(param_dates)
        n = mod(j,size(colors,1));
        n(n==0) = size(colors,1);
        if j == length(param_dates)
            time = XLim(2)*ones(1,4);
            time([1,4]) = days(param_dates(j)-param_dates(1))+10;
        else
            time = days([param_dates(j) param_dates(j+1) param_dates(j+1) param_dates(j)]-param_dates(1))+10;
        end
        fill(time,YLim(2)+[-bh,-bh,0,0],colors(n,:));
        xline(time(2));
    end    
    if ~isempty(sub_tab)
        i1 = contains(sub_tab.Condition,'Pre-Op'); %Visit 0
        i2 = contains(sub_tab.Condition,'NoStim'); %Visit 3 and anothers
        i3 = contains(sub_tab.Condition,'MotionMod'); 
        i4 = contains(sub_tab.Condition,'ConstantRate');
        switch ear
            case 'R'
                errorbar(act_days(i1),sub_tab.LHRH_L_POS_Gain(i1),sub_tab.LHRH_L_POS_Gain_sd(i1),'rx')
                errorbar(act_days(i1),sub_tab.LHRH_R_POS_Gain(i1),sub_tab.LHRH_R_POS_Gain_sd(i1),'mx')
                errorbar(act_days(i2),sub_tab.LHRH_L_POS_Gain(i2),sub_tab.LHRH_L_POS_Gain_sd(i2),'rx')
                errorbar(act_days(i2),sub_tab.LHRH_R_POS_Gain(i2),sub_tab.LHRH_R_POS_Gain_sd(i2),'mx')
                errorbar(act_days(i3),sub_tab.LHRH_L_POS_Gain(i3),sub_tab.LHRH_L_POS_Gain_sd(i3),'r*-')
                errorbar(act_days(i3),sub_tab.LHRH_R_POS_Gain(i3),sub_tab.LHRH_R_POS_Gain_sd(i3),'m*-')
                errorbar(act_days(i4),sub_tab.LHRH_L_POS_Gain(i4),sub_tab.LHRH_L_POS_Gain_sd(i4),'ro:')
                errorbar(act_days(i4),sub_tab.LHRH_R_POS_Gain(i4),sub_tab.LHRH_R_POS_Gain_sd(i4),'mo:')
            case 'L'
                errorbar(act_days(i1),sub_tab.LHRH_L_NEG_Gain(i1),sub_tab.LHRH_L_NEG_Gain_sd(i1),'rx')
                errorbar(act_days(i1),sub_tab.LHRH_R_NEG_Gain(i1),sub_tab.LHRH_R_NEG_Gain_sd(i1),'mx')
                errorbar(act_days(i2),sub_tab.LHRH_L_NEG_Gain(i2),sub_tab.LHRH_L_NEG_Gain_sd(i2),'rx')
                errorbar(act_days(i2),sub_tab.LHRH_R_NEG_Gain(i2),sub_tab.LHRH_R_NEG_Gain_sd(i2),'mx')
                errorbar(act_days(i3),sub_tab.LHRH_L_NEG_Gain(i3),sub_tab.LHRH_L_NEG_Gain_sd(i3),'r*-')
                errorbar(act_days(i3),sub_tab.LHRH_R_NEG_Gain(i3),sub_tab.LHRH_R_NEG_Gain_sd(i3),'m*-')
                errorbar(act_days(i4),sub_tab.LHRH_L_NEG_Gain(i4),sub_tab.LHRH_L_NEG_Gain_sd(i4),'ro:')
                errorbar(act_days(i4),sub_tab.LHRH_R_NEG_Gain(i4),sub_tab.LHRH_R_NEG_Gain_sd(i4),'mo:')
        end
    end
    hold off
end
%savefig([subject{:},'_RotaryChairGain_ContraHalfCycle.fig'])
%%  Plot phase (one subject)
excludedatawithgainbelow=0.015;
%Plot values
XLim = [4 1000];
xtick = [5,28*[0,1,6,24]+10];
xticklab = [-1,0,1,6,24];
YLim = [-180 180];
ytick = YLim(1):45:YLim(2);
%Inhibition Half-Cycle
figure(3)
colors = get(gca,'ColorOrder');
colors = [colors(3:end,:);colors(1:2,:)];
ha = gobjects(length(freqs),1);
subplot(1,1,1)
delete(findall(gcf,'type','annotation'))
annotation('textbox',[0 .9 1 .1],'String',[subject{:},' Horizontal Rotary Chair'],'FontSize',12,...
    'HorizontalAlignment','center','EdgeColor','none');
for i = 1:length(freqs)
    %Isolate values
    sub_tab = tab(strcmp(tab.Frequency,freqs{i}),:);
    act_days = days(sub_tab.Date-param_dates(1))+10;
    act_days(act_days < 5) = 5; %Move Pre-Op Value 
    %Make axes
    ha(i) = subplot(2,4,i);
    hold on
    set(gca,'XLim',XLim,'XTick',xtick,'xscale','log','FontSize',8,'box','on')
    set(gca,'XGrid',1,'XMinorGrid',0,'YGrid',1,'YMinorGrid',0)
    if i > 4 %bottom row of graphs
        set(gca,'Ylim',YLim,'XTickLabels',xticklab)
        xlabel('Months Since Activation','FontSize',8)
        set(gca,'Position',[0.06+(0.225*(i-5)),0.05,0.215,0.4])
        bh = 10;
    else 
        set(gca,'Ylim',YLim,'XTickLabels',[])
        set(gca,'Position',[0.06+0.225*(i-1),0.5,0.215,0.4])
        bh = 10;
    end
    set(gca,'YTick',ytick,'YTickLabel',[])
    if i == 1 || i ==5
        set(gca,'YTickLabel',ytick)
        ylabel('Phase w/re to Head Velocity')
    end    
    set(gca,'XMinorTick','off')
    title(freqs{i})
       
    %Plot
    fill([4 6 6 4],YLim(2)+[-bh,-bh,0,0],0.75*[1,1,1]);
    xline(6);
    fill([6 10 10 6],YLim(2)+[-bh,-bh,0,0],0.85*[1,1,1]);
    xline(10);
    for j = 1:length(param_dates)
        n = mod(j,size(colors,1));
        n(n==0) = size(colors,1);
        if j == length(param_dates)
            time = XLim(2)*ones(1,4);
            time([1,4]) = days(param_dates(j)-param_dates(1))+10;
        else
            time = days([param_dates(j) param_dates(j+1) param_dates(j+1) param_dates(j)]-param_dates(1))+10;
        end
        fill(time,YLim(2)+[-bh,-bh,0,0],colors(n,:));
        xline(time(2));
    end 
    plot(XLim,[0 0],'k-')
    if ~isempty(sub_tab)
        i1 = contains(sub_tab.Condition,'Pre-Op'); %Visit 0
        i2 = contains(sub_tab.Condition,'NoStim'); %Visit 3 and anothers
        i3 = contains(sub_tab.Condition,'MotionMod'); 
        i4 = contains(sub_tab.Condition,'ConstantRate');
        switch ear
            case 'L'
                keep_L = sub_tab.LHRH_L_POS_Gain > excludedatawithgainbelow;
                keep_R = sub_tab.LHRH_R_POS_Gain > excludedatawithgainbelow;
            case 'R'
                keep_L = sub_tab.LHRH_L_NEG_Gain > excludedatawithgainbelow;
                keep_R = sub_tab.LHRH_R_NEG_Gain > excludedatawithgainbelow;
        end
        errorbar(act_days(i1&keep_L),sub_tab.L_Phase(i1&keep_L),sub_tab.L_Phase_sd(i1&keep_L),'rx')
        errorbar(act_days(i1&keep_R),sub_tab.R_Phase(i1&keep_R),sub_tab.R_Phase_sd(i1&keep_R),'mx')
        errorbar(act_days(i2&keep_L),sub_tab.L_Phase(i2&keep_L),sub_tab.L_Phase_sd(i2&keep_L),'rx')
        errorbar(act_days(i2&keep_R),sub_tab.R_Phase(i2&keep_R),sub_tab.R_Phase_sd(i2&keep_R),'mx')
        errorbar(act_days(i3&keep_L),sub_tab.L_Phase(i3&keep_L),sub_tab.L_Phase_sd(i3&keep_L),'r*-')
        errorbar(act_days(i3&keep_R),sub_tab.R_Phase(i3&keep_R),sub_tab.R_Phase_sd(i3&keep_R),'m*-')
        errorbar(act_days(i4&keep_L),sub_tab.L_Phase(i4&keep_L),sub_tab.L_Phase_sd(i4&keep_L),'ro:')
        errorbar(act_days(i4&keep_R),sub_tab.R_Phase(i4&keep_R),sub_tab.R_Phase_sd(i4&keep_R),'mo:')        
    end
    hold off
end
%savefig([subject{:},'_RotaryChairPhase.fig'])
%% Plot Alignment (one subject, ipsi half-cycle)
%Plot values
XLim = [4 1000];
xtick = [5,28*[0,1,6,24]+10];
xticklab = [-1,0,1,6,24];
YLim = [-5 180];
ytick = 0:20:YLim(2);
bh = 10;
%Stimulation Half-Cycle
figure(4)
colors = get(gca,'ColorOrder');
colors = [colors(3:end,:);colors(1:2,:)];
ha = gobjects(length(freqs),1);
subplot(1,1,1)
delete(findall(gcf,'type','annotation'))
annotation('textbox',[0 .9 1 .1],'String',[subject{:},' Horizontal Rotary Chair: Implant Stimulation Half-Cycle'],'FontSize',12,...
    'HorizontalAlignment','center','EdgeColor','none');
for i = 1:length(freqs)
    %Isolate values
    sub_tab = tab(strcmp(tab.Frequency,freqs{i}),:);
    act_days = days(sub_tab.Date-param_dates(1))+10;
    act_days(act_days < 5) = 5; %Move Pre-Op Value 
    %Make axes
    ha(i) = subplot(2,4,i);
    hold on
    set(gca,'XLim',XLim,'XTick',xtick,'xscale','log','FontSize',8,'box','on')
    set(gca,'XGrid',1,'XMinorGrid',0,'YGrid',1,'YMinorGrid',0)
    if i > 4 %bottom row of graphs
        set(gca,'Ylim',YLim,'XTickLabels',xticklab)
        xlabel('Months Since Activation','FontSize',8)
        set(gca,'Position',[0.06+(0.225*(i-5)),0.05,0.215,0.4])
    else 
        set(gca,'Ylim',YLim,'XTickLabels',[])
        set(gca,'Position',[0.06+0.225*(i-1),0.5,0.215,0.4])
    end
    set(gca,'YTick',ytick,'YTickLabel',[])
    if i == 1 || i ==5
        set(gca,'YTickLabel',ytick)
        ylabel('Alignment w/re to Target Canal')
    end    
    set(gca,'XMinorTick','off')
    title(freqs{i})
    
    %Plot
    fill([4 6 6 4],YLim(2)+[-bh,-bh,0,0],0.75*[1,1,1]);
    xline(6);
    fill([6 10 10 6],YLim(2)+[-bh,-bh,0,0],0.85*[1,1,1]);
    xline(10);
    for j = 1:length(param_dates)
        n = mod(j,size(colors,1));
        n(n==0) = size(colors,1);
        if j == length(param_dates)
            time = XLim(2)*ones(1,4);
            time([1,4]) = days(param_dates(j)-param_dates(1))+10;
        else
            time = days([param_dates(j) param_dates(j+1) param_dates(j+1) param_dates(j)]-param_dates(1))+10;
        end
        fill(time,YLim(2)+[-bh,-bh,0,0],colors(n,:));
        xline(time(2));
    end    
    plot(XLim,[0 0],'k-')
    if ~isempty(sub_tab)
        i1 = contains(sub_tab.Condition,'Pre-Op'); %Visit 0
        i2 = contains(sub_tab.Condition,'NoStim'); %Visit 3 and anothers
        i3 = contains(sub_tab.Condition,'MotionMod'); 
        i4 = contains(sub_tab.Condition,'ConstantRate');
        switch ear
            case 'L'
                errorbar(act_days(i1),sub_tab.L_POS_Align(i1),sub_tab.L_POS_Align_sd(i1),'rx')
                errorbar(act_days(i1),sub_tab.R_POS_Align(i1),sub_tab.R_POS_Align_sd(i1),'mx')
                errorbar(act_days(i2),sub_tab.L_POS_Align(i2),sub_tab.L_POS_Align_sd(i2),'rx')
                errorbar(act_days(i2),sub_tab.R_POS_Align(i2),sub_tab.R_POS_Align_sd(i2),'mx')
                errorbar(act_days(i3),sub_tab.L_POS_Align(i3),sub_tab.L_POS_Align_sd(i3),'r*-')
                errorbar(act_days(i3),sub_tab.R_POS_Align(i3),sub_tab.R_POS_Align_sd(i3),'m*-')
                errorbar(act_days(i4),sub_tab.L_POS_Align(i4),sub_tab.L_POS_Align_sd(i4),'ro:')
                errorbar(act_days(i4),sub_tab.R_POS_Align(i4),sub_tab.R_POS_Align_sd(i4),'mo:')
            case 'R'
                errorbar(act_days(i1),sub_tab.L_NEG_Align(i1),sub_tab.L_NEG_Align_sd(i1),'rx')
                errorbar(act_days(i1),sub_tab.R_NEG_Align(i1),sub_tab.R_NEG_Align_sd(i1),'mx')
                errorbar(act_days(i2),sub_tab.L_NEG_Align(i2),sub_tab.L_NEG_Align_sd(i2),'rx')
                errorbar(act_days(i2),sub_tab.R_NEG_Align(i2),sub_tab.R_NEG_Align_sd(i2),'mx')
                errorbar(act_days(i3),sub_tab.L_NEG_Align(i3),sub_tab.L_NEG_Align_sd(i3),'r*-')
                errorbar(act_days(i3),sub_tab.R_NEG_Align(i3),sub_tab.R_NEG_Align_sd(i3),'m*-')
                errorbar(act_days(i4),sub_tab.L_NEG_Align(i4),sub_tab.L_NEG_Align_sd(i4),'ro:')
                errorbar(act_days(i4),sub_tab.R_NEG_Align(i4),sub_tab.R_NEG_Align_sd(i4),'mo:')
        end
    end
    hold off
end
%savefig([subject{:},'_RotaryChairAlignment_IpsiHalfCycle.fig'])
%% Plot Alignment (one subject, contra half-cycle)
%Plot values
XLim = [4 1000];
xtick = [5,28*[0,1,6,24]+10];
xticklab = [-1,0,1,6,24];
YLim = [-5 180];
ytick = 0:20:YLim(2);
bh = 10;
%Inhibition Half-Cycle
figure(5)
colors = get(gca,'ColorOrder');
colors = [colors(3:end,:);colors(1:2,:)];
ha = gobjects(length(freqs),1);
subplot(1,1,1)
delete(findall(gcf,'type','annotation'))
annotation('textbox',[0 .9 1 .1],'String',[subject{:},' Horizontal Rotary Chair: Implant Inhibition Half-Cycle'],'FontSize',12,...
    'HorizontalAlignment','center','EdgeColor','none');
for i = 1:length(freqs)
    %Isolate values
    sub_tab = tab(strcmp(tab.Frequency,freqs{i}),:);
    act_days = days(sub_tab.Date-param_dates(1))+10;
    act_days(act_days < 5) = 5; %Move Pre-Op Value 
    %Make axes
    ha(i) = subplot(2,4,i);
    hold on
    set(gca,'XLim',XLim,'XTick',xtick,'xscale','log','FontSize',8,'box','on')
    set(gca,'XGrid',1,'XMinorGrid',0,'YGrid',1,'YMinorGrid',0)
    if i > 4 %bottom row of graphs
        set(gca,'Ylim',YLim,'XTickLabels',xticklab)
        xlabel('Months Since Activation','FontSize',8)
        set(gca,'Position',[0.06+(0.225*(i-5)),0.05,0.215,0.4])
    else 
        set(gca,'Ylim',YLim,'XTickLabels',[])
        set(gca,'Position',[0.06+0.225*(i-1),0.5,0.215,0.4])
    end
    set(gca,'YTick',ytick,'YTickLabel',[])
    if i == 1 || i ==5
        set(gca,'YTickLabel',ytick)
        ylabel('Alignment w/re to Target Canal')
    end    
    set(gca,'XMinorTick','off')
    title(freqs{i})
    
    %Plot
    fill([4 6 6 4],YLim(2)+[-bh,-bh,0,0],0.75*[1,1,1]);
    xline(6);
    fill([6 10 10 6],YLim(2)+[-bh,-bh,0,0],0.85*[1,1,1]);
    xline(10);
    for j = 1:length(param_dates)
        n = mod(j,size(colors,1));
        n(n==0) = size(colors,1);
        if j == length(param_dates)
            time = XLim(2)*ones(1,4);
            time([1,4]) = days(param_dates(j)-param_dates(1))+10;
        else
            time = days([param_dates(j) param_dates(j+1) param_dates(j+1) param_dates(j)]-param_dates(1))+10;
        end
        fill(time,YLim(2)+[-bh,-bh,0,0],colors(n,:));
        xline(time(2));
    end    
    plot(XLim,[0 0],'k-')
    if ~isempty(sub_tab)
        i1 = contains(sub_tab.Condition,'Pre-Op'); %Visit 0
        i2 = contains(sub_tab.Condition,'NoStim'); %Visit 3 and anothers
        i3 = contains(sub_tab.Condition,'MotionMod'); 
        i4 = contains(sub_tab.Condition,'ConstantRate');
        switch ear
            case 'R'
                errorbar(act_days(i1),sub_tab.L_POS_Align(i1),sub_tab.L_POS_Align_sd(i1),'rx')
                errorbar(act_days(i1),sub_tab.R_POS_Align(i1),sub_tab.R_POS_Align_sd(i1),'mx')
                errorbar(act_days(i2),sub_tab.L_POS_Align(i2),sub_tab.L_POS_Align_sd(i2),'rx')
                errorbar(act_days(i2),sub_tab.R_POS_Align(i2),sub_tab.R_POS_Align_sd(i2),'mx')
                errorbar(act_days(i3),sub_tab.L_POS_Align(i3),sub_tab.L_POS_Align_sd(i3),'r*-')
                errorbar(act_days(i3),sub_tab.R_POS_Align(i3),sub_tab.R_POS_Align_sd(i3),'m*-')
                errorbar(act_days(i4),sub_tab.L_POS_Align(i4),sub_tab.L_POS_Align_sd(i4),'ro:')
                errorbar(act_days(i4),sub_tab.R_POS_Align(i4),sub_tab.R_POS_Align_sd(i4),'mo:')
            case 'L'
                errorbar(act_days(i1),sub_tab.L_NEG_Align(i1),sub_tab.L_NEG_Align_sd(i1),'rx')
                errorbar(act_days(i1),sub_tab.R_NEG_Align(i1),sub_tab.R_NEG_Align_sd(i1),'mx')
                errorbar(act_days(i2),sub_tab.L_NEG_Align(i2),sub_tab.L_NEG_Align_sd(i2),'rx')
                errorbar(act_days(i2),sub_tab.R_NEG_Align(i2),sub_tab.R_NEG_Align_sd(i2),'mx')
                errorbar(act_days(i3),sub_tab.L_NEG_Align(i3),sub_tab.L_NEG_Align_sd(i3),'r*-')
                errorbar(act_days(i3),sub_tab.R_NEG_Align(i3),sub_tab.R_NEG_Align_sd(i3),'m*-')
                errorbar(act_days(i4),sub_tab.L_NEG_Align(i4),sub_tab.L_NEG_Align_sd(i4),'ro:')
                errorbar(act_days(i4),sub_tab.R_NEG_Align(i4),sub_tab.R_NEG_Align_sd(i4),'mo:')
        end
    end
    hold off
end
%savefig([subject{:},'_RotaryChairAlignment_ContraHalfCycle.fig'])