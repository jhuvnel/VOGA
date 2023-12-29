%% Run in a Cycle Average Folder to test agreement with CycAvg MaxVel 
tic;
Cyc_Path = [cd,filesep,'Cycle Averages',filesep];
files = extractfield(dir([Cyc_Path,'CycAvg*Sine*.mat']),'name');
txt_var = {'File','Type','AxisName'};
num_var = {'Frequency','Amplitude','MaxVel','MaxVel_sd','MedMaxVel'};
med_results = [cell2table(cell(length(files)*2,length(txt_var)),'VariableNames',txt_var),...
    array2table(NaN(length(files)*2,length(num_var)),'VariableNames',num_var)];
ax_name = {'LHRH','LARP','RALP','X','Y'};
for i = 1:length(files)
    load([Cyc_Path,files{i}])
    ax = strrep(ax_name{cellfun(@(x) contains(files{i},x),ax_name)},'LHRH','Z');
    Data = angpos2angvel(CycAvg.Data);
    cyc_inds = CycAvg.Data_allcyc.keep_inds;    
    traces = tvd1d([median(Data.(['LE_Vel_',ax])(cyc_inds),2,'omitnan');...
        median(Data.(['RE_Vel_',ax])(cyc_inds),2,'omitnan')],2);
    med_results([2*i-1,2*i],[txt_var,num_var(1:end-1)]) = ...
        CycAvg.parameterized(:,[txt_var,num_var(1:end-1)]);
    med_results.MedMaxVel([2*i-1,2*i]) = [abs(min(traces));abs(max(traces))];
end
med_results.VelDiff = med_results.MedMaxVel-med_results.MaxVel;
med_results.PercDiff = 100*abs(med_results.MedMaxVel-med_results.MaxVel)./med_results.MaxVel;
med_results.Zscore = abs(med_results.MedMaxVel-med_results.MaxVel)./med_results.MaxVel_sd;
med_results = sortrows(sortrows(med_results,"Amplitude","ascend"),"AxisName","descend");
disp(toc)
writetable(med_results,'MedianCycResults.xlsx')
histogram(med_results.VelDiff,-20:5:140)
xlabel('Error between maximum cycle median and average')
ylabel('Occurances')
%% LHRH 300dps plot
load([Cyc_Path,files{contains(files,'LHRH')&contains(files,'300dps')}])
Data = angpos2angvel(CycAvg.Data);
cyc_inds = CycAvg.Data_allcyc.keep_inds;    
LE_med = tvd1d(median(Data.LE_Vel_Z(cyc_inds),2,'omitnan'),2);
RE_med = tvd1d(median(Data.RE_Vel_Z(cyc_inds),2,'omitnan'),2);
figure;
h1 = gobjects(6,1);
h1(2) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.lz_cycavg-CycAvg.lz_cycstd,fliplr(CycAvg.lz_cycavg+CycAvg.lz_cycstd)],[1,0,0],'FaceAlpha',0.3,'EdgeColor','r');
hold on
h1(5) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.rz_cycavg-CycAvg.rz_cycstd,fliplr(CycAvg.rz_cycavg+CycAvg.rz_cycstd)],[1,0,1],'FaceAlpha',0.3,'EdgeColor','m');
h1(1) = plot(CycAvg.t,CycAvg.lz_cycavg,'r','LineWidth',2);
h1(4) = plot(CycAvg.t,CycAvg.rz_cycavg,'m','LineWidth',2);
h1(3) = plot(CycAvg.t,LE_med,'r:','LineWidth',2);
h1(6) = plot(CycAvg.t,RE_med,'m:','LineWidth',2);
hold off
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
leg = legend(h1,{'LEye Mean','LEye SD','LEye Median','REye Mean','REye SD','REye Median'},'NumColumns',2,'Location','northwest');
title(leg,'Horizontal Response')
leg.ItemTokenSize(1) = 20;
%% LHRH 0.1Hz plot
load([Cyc_Path,files{contains(files,'LHRH')&contains(files,'0.1Hz')}])
Data = angpos2angvel(CycAvg.Data);
cyc_inds = CycAvg.Data_allcyc.keep_inds;    
LE_med = tvd1d(median(Data.LE_Vel_Z(cyc_inds),2,'omitnan'),2);
RE_med = tvd1d(median(Data.RE_Vel_Z(cyc_inds),2,'omitnan'),2);
figure;
h1 = gobjects(6,1);
hold on
h1(2) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.lz_cycavg-CycAvg.lz_cycstd,fliplr(CycAvg.lz_cycavg+CycAvg.lz_cycstd)],[1,0,0],'FaceAlpha',0.3,'EdgeColor','r');
h1(5) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.rz_cycavg-CycAvg.rz_cycstd,fliplr(CycAvg.rz_cycavg+CycAvg.rz_cycstd)],[1,0,1],'FaceAlpha',0.3,'EdgeColor','m');
h1(1) = plot(CycAvg.t,CycAvg.lz_cycavg,'r','LineWidth',2);
h1(4) = plot(CycAvg.t,CycAvg.rz_cycavg,'m','LineWidth',2);
h1(3) = plot(CycAvg.t,LE_med,'r:','LineWidth',1.5);
h1(6) = plot(CycAvg.t,RE_med,'m:','LineWidth',1.5);
hold off
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
set(gca,'YLim',[-135 85])
leg = legend(h1,{'LEye Mean','LEye SD','LEye Median','REye Mean','REye SD','REye Median'},'NumColumns',2,'Location','southwest');
title(leg,'Horizontal Response')
leg.ItemTokenSize(1) = 20;
%% LARP 200dps
load([Cyc_Path,files{contains(files,'LARP')&contains(files,'200dps')}])
Data = angpos2angvel(CycAvg.Data);
cyc_inds = CycAvg.Data_allcyc.keep_inds;    
LE_med = tvd1d(median(Data.LE_Vel_LARP(cyc_inds),2,'omitnan'),2);
RE_med = tvd1d(median(Data.RE_Vel_LARP(cyc_inds),2,'omitnan'),2);
figure;
h1 = gobjects(8,1);
hold on
h1(4) = plot(NaN,NaN,'Color',[0 0.5 0],'LineWidth',0.5);
h1(8) = plot(NaN,NaN,'g','LineWidth',0.5);
%plot(CycAvg.t,Data.LE_Vel_LARP(cyc_inds),'Color',[0 0.5 0],'LineWidth',0.5)
%plot(CycAvg.t,Data.RE_Vel_LARP(cyc_inds),'g','LineWidth',0.5)
plot(CycAvg.t,CycAvg.Data_allcyc.LE_Vel_LARP,'Color',[0 0.5 0],'LineWidth',0.5)
plot(CycAvg.t,CycAvg.Data_allcyc.RE_Vel_LARP,'g','LineWidth',0.5)
h1(2) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.ll_cycavg-CycAvg.ll_cycstd,fliplr(CycAvg.ll_cycavg+CycAvg.ll_cycstd)],[0,0.5,0],'FaceAlpha',0.3,'EdgeColor',[0,0.5,0]);
h1(6) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.rl_cycavg-CycAvg.rl_cycstd,fliplr(CycAvg.rl_cycavg+CycAvg.rl_cycstd)],[0,1,0],'FaceAlpha',0.3,'EdgeColor',[0,1,0]);
h1(1) = plot(CycAvg.t,CycAvg.ll_cycavg,'Color',[0,0.5,0],'LineWidth',2);
h1(5) = plot(CycAvg.t,CycAvg.rl_cycavg,'g','LineWidth',2);
h1(3) = plot(CycAvg.t,LE_med,':','Color',[0,0.5,0],'LineWidth',2);
h1(7) = plot(CycAvg.t,RE_med,'g:','LineWidth',2);
hold off
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
set(gca,'YLim',[-120 180])
leg = legend(h1,{'LEye Mean','LEye SD','LEye Median','LEye FiltCycles','REye Mean','REye SD','REye Median','REye FiltCycles'},'NumColumns',2,'Location','northwest');
title(leg,'LARP Response')
leg.ItemTokenSize(1) = 20;
%% 
