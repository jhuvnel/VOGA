%% Plot Rotary Chair Combined Subjects
%Plots Gain (stim half-cycle), phase, and alignment (stim
%half-cycle) for multiple subject and multiple conditions (most recent motion mod/constant rate)
%over frequency
clc;
clear;
close all;
%% Normative Rotary Chair Data
notes = 'Normative Data from Wall et. al. 1984 on patients 50-69';
freq = [0.005 0.010 0.020 0.050 0.100 0.200 0.500 1.000];
gain = [0.2175	0.3576	0.4803	0.5634	0.56	0.5445	0.5678	0.698];
gain_std = [0.0757	0.1235	0.1474	0.2243	0.2363	0.2509	0.2121	0.2777];
phase = [70.9478	43.6413	25.9262	10.7906	4.084	-4.8461	-11.4084 -11.3412];
phase_std = [7.5044	6.526	6.714	8.0519	4.5136	7.6279	6.2448	11.3167];
norm_dat.notes = notes;
norm_dat.freq = freq;
norm_dat.gain = gain;
norm_dat.gain_std = gain_std;
norm_dat.phase = phase;
norm_dat.phase_std = phase_std;
save('RotaryChairNormativeData.mat','norm_dat')
%%
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
Ears = {'L','L','L','L','R','R','L','R'};
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
%% Get information of what to plot
[path2,path1] = uigetfile('*.mat','Select Parameter Fit File');
load([path1,filesep,path2],'all_results')
freqs = {'0.01Hz','0.02Hz','0.05Hz','0.1Hz','0.2Hz','0.5Hz','1Hz','2Hz'};
subjects = unique(all_results.Subject);
excludedatawithgainbelow=0.01;
markers = ['x','d','o','^','p','s','+','h'];
%% Run for each subject
%Make the items to plot
gain = NaN(4,length(freqs),length(subjects));
gain_sd = NaN(4,length(freqs),length(subjects));
phase = NaN(4,length(freqs),length(subjects));
phase_sd = NaN(4,length(freqs),length(subjects));
sub_num = cell(1,length(subjects));
for j = 1:length(subjects)
    subject = subjects(j);
    tab = sortrows(all_results(contains(all_results.Subject,subject),:),'Date','ascend');
    visits = unique(tab.Visit,'stable');
    %Subject Information
    sub = str2double(subject{:}(6));
    sub_num{j} = num2str(sub);
    ear = Ears{sub};
    param_dates = Param_Changes{sub};
    if ismember('Visit0',visits) %Don't use subjects w/o pre-op data        
        % Make arrays for plotting
        %Find the relevant index
        ind = [];
        ind(:,1) = contains(tab.Condition,'Pre-Op');
        ind(:,2) = contains(tab.Condition,'NoStim')&tab.Date==tab.Date(find(contains(tab.Condition,'NoStim'),1,'first'));
        m_dat = unique(tab.Date(contains(tab.Condition,'MotionMod')),'stable');
        ind(:,3) = contains(tab.Condition,'MotionMod')&tab.Date==m_dat(end);
        c_dat = unique(tab.Date(contains(tab.Condition,'ConstantRate')),'stable');
        ind(:,4) = contains(tab.Condition,'ConstantRate')&tab.Date==c_dat(end);
        %Sometimes you can't record a Motion Mod or Constant Rate value for a visit
        %in which case you want that value from the previous visit
        if length(find(ind(:,3)))<length(freqs)
            miss_freqs = freqs(~ismember(freqs,tab.Frequency(logical(ind(:,3)))));
            ind(:,3) = ind(:,3)|contains(tab.Condition,'MotionMod')&tab.Date==c_dat(end-1)&ismember(tab.Frequency,miss_freqs);
        end
        if length(find(ind(:,4)))<length(freqs)
            miss_freqs = freqs(~ismember(freqs,tab.Frequency(logical(ind(:,4)))));
            ind(:,4) = ind(:,4)|contains(tab.Condition,'ConstantRate')&tab.Date==c_dat(end-1)&ismember(tab.Frequency,miss_freqs);
        end
        for i = 1:4
            inds = find(ind(:,i));
            [~,f_ind] = ismember(tab.Frequency(inds),freqs);
            if strcmp(ear,'L')
                [sub_gain,eye] = max([tab.LHRH_L_POS_Gain(inds)';tab.LHRH_R_POS_Gain(inds)']);
                sub_gain_sd = tab.LHRH_L_POS_Gain_sd(inds)';
                sub_gain_sd(eye==2) = tab.LHRH_R_POS_Gain_sd(inds(eye==2));
            elseif strcmp(ear,'R')
                [sub_gain,eye] = max([tab.LHRH_L_NEG_Gain(inds)';tab.LHRH_R_NEG_Gain(inds)']);
                sub_gain_sd = tab.LHRH_L_NEG_Gain_sd(inds)';
                sub_gain_sd(eye==2) = tab.LHRH_R_NEG_Gain_sd(inds(eye==2));
            end
            sub_phase = tab.L_Phase(inds);
            sub_phase(eye==2) = tab.R_Phase(inds(eye==2));
            sub_phase_sd = tab.L_Phase_sd(inds);
            sub_phase_sd(eye==2) = tab.R_Phase_sd(inds(eye==2));
            sub_phase(sub_gain < excludedatawithgainbelow) = NaN;
            sub_phase_sd(sub_gain < excludedatawithgainbelow) = NaN;
            %Set the values
            gain(i,f_ind,j) = sub_gain;
            gain_sd(i,f_ind,j) = sub_gain_sd;
            phase(i,f_ind,j) = sub_phase;
            phase_sd(i,f_ind,j) = sub_phase_sd;   
        end
    end
end

%Aggregate values
m_gain = mean(gain,3,'omitnan');
m_gain_sd = std(gain,[],3,'omitnan');
m_phase = circ_mean(phase*pi/180,[],3)*180/pi;
m_phase_sd = circ_std(phase*pi/180,[],[],3)*180/pi;
m_phase(isnan(m_gain)|m_gain<excludedatawithgainbelow) = NaN;
m_phase_sd(isnan(m_gain)|m_gain<excludedatawithgainbelow) = NaN;
%% Plot
fplot = 1:8; %Which frequencies to plot
logxshift=1.03; %how much to multiply the x value to offset its marker rightward (divide to move left)
markerbig=8;
markersmall=6;
linethick=2;
linethin=1;
errorbarcapsize=1;
graydark=0.85;
graylight=0.95;
figsizeinches=[7,6];
figsizeinchesBoxplot=[2.3,4];
freq_ax = cellfun(@str2num,strrep(freqs(fplot),'Hz',''));
plot_colors = [0,0,0; ...
              [1,1,1]*110/255; ...
              1,0,0;  ...
              1,0,0];
plot_linestyle = {'-','-.','-',':'};
plot_offset = [logxshift^-1, logxshift^-2, 1, logxshift, logxshift^2];
fig = figure;
set(fig,'Units','inch','Position',[2 2 figsizeinches]);%CDS083119a
ha(1) = subplot(2,1,1);
ha(1).Position = [0.1,0.55,0.85,0.4];
ha(2) = subplot(2,1,2);
ha(2).Position = [0.1,0.10,0.85,0.4];
set(ha,'XLim',[0.9*freq_ax(1),2.2],'YTickLabelMode','auto','XScale','log')
set(ha(1),'YLim',[0 0.85],'YTick',0.1:0.1:0.8)
set(ha(2),'YLim',[-30,100],'YTick',-20:20:100)
axes(ha(1))
hold on
h(7) = fill([freq,fliplr(freq)],[norm_gain_m+2*norm_gain_std,fliplr(norm_gain_m-2*norm_gain_std)],graylight*[1 1 1],'LineStyle','none');
h(6) = fill([freq,fliplr(freq)],[norm_gain_m+1*norm_gain_std,fliplr(norm_gain_m-1*norm_gain_std)],graydark*[1 1 1],'LineStyle','none');
h(5) = plot(freq,norm_gain_m,'k--','LineWidth',linethick);
for i = 1:4
    h(i) = plot(plot_offset(i)*freq_ax,m_gain(i,fplot),'Color',plot_colors(i,:),'LineStyle',plot_linestyle{i},'LineWidth',linethick);
    errorbar(plot_offset(i)*freq_ax,m_gain(i,fplot),m_gain_sd(i,fplot),'Color',plot_colors(i,:),'LineStyle','none','LineWidth',linethin,'CapSize',errorbarcapsize)
end
if sub == 1
    freqs2 = [0.01;0.02;0.04;0.08;0.16;0.32;0.64];
    gains2 = [0.0083;0.0250;0.0292;0.0458;0.0500;0.0333;0.0625];
    plot(freqs2,gains2,'kx','LineWidth',linethick)
end
for i = 1:length(subjects)
    plot(plot_offset(1)*freq_ax,gain(1,fplot,i),'k.','Marker',markers(i),'MarkerFaceColor','k','MarkerSize',markersmall,'LineWidth',linethin);
    plot(plot_offset(2)*freq_ax,gain(2,fplot,i),'.','Color',[110,110,110]/255,'Marker',markers(i),'MarkerSize',markersmall,'LineWidth',linethin);
    plot(plot_offset(3)*freq_ax,gain(3,fplot,i),'r.','Marker',markers(i),'MarkerFaceColor','r','MarkerSize',markersmall,'LineWidth',linethin);
    plot(plot_offset(4)*freq_ax,gain(4,fplot,i),'r.','Marker',markers(i),'MarkerSize',markersmall,'LineWidth',linethin);
end
hold off
legend(h,{'Pre-Op','Post-Op','MVI ON','MVI OFF','Normal mean','Normal±1SD','Normal±2SD'},'Location','NorthWest','NumColumns',2)
ylabel('Horizontal VOR Gain')
set(ha(1),'XTick',freq_ax,'XTickLabel',[])
axes(ha(2))
hold on
fill([freq,fliplr(freq)],[norm_phase_m+2*norm_phase_std,fliplr(norm_phase_m-2*norm_phase_std)],graylight*[1 1 1],'LineStyle','none')
fill([freq,fliplr(freq)],[norm_phase_m+1*norm_phase_std,fliplr(norm_phase_m-1*norm_phase_std)],graydark*[1 1 1],'LineStyle','none')
plot(freq,norm_phase_m,'k--','LineWidth',linethick)
for i = 1:4
    plot(plot_offset(i)*freq_ax,m_phase(i,fplot),'Color',plot_colors(i,:),'LineStyle',plot_linestyle{i},'LineWidth',linethick);
    errorbar(plot_offset(i)*freq_ax,m_phase(i,fplot),m_phase_sd(i,fplot),'Color',plot_colors(i,:),'LineStyle','none','LineWidth',linethin,'CapSize',errorbarcapsize)
end
for i = 1:length(subjects)
    plot(plot_offset(1)*freq_ax,phase(1,fplot,i),'k.','Marker',markers(i),'MarkerFaceColor','k','MarkerSize',markersmall,'LineWidth',linethin);
    plot(plot_offset(2)*freq_ax,phase(2,fplot,i),'.','Color',[110,110,110]/255,'Marker',markers(i),'MarkerSize',markersmall,'LineWidth',linethin);
    plot(plot_offset(3)*freq_ax,phase(3,fplot,i),'r.','Marker',markers(i),'MarkerFaceColor','r','MarkerSize',markersmall,'LineWidth',linethin);
    plot(plot_offset(4)*freq_ax,phase(4,fplot,i),'r.','Marker',markers(i),'MarkerSize',markersmall,'LineWidth',linethin);
    h2(i) = plot(NaN,NaN,'k.','Marker',markers(i),'MarkerSize',markersmall,'LineWidth',linethin);
end
%plot([0.01,3],[0,0],'k:')
hold off
leg = legend(h2,sub_num,'Location','northwest','NumColumns',2);
leg.ItemTokenSize(1) = 6;
title(leg,'Subjects')
ylabel('Phase Lead (deg)')
set(ha(2),'XTick',freq_ax,'XTickLabel',strrep(freqs(fplot),'Hz',''))
xlabel('Frequency [Hz]')
savefig('Aggregated_RotaryChair_SummaryAcrossFreq.fig')