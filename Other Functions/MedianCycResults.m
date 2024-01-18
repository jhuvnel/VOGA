%% Run in a Cycle Average Folder to test timing and agreement with CycAvg MaxVel 
Cyc_Path = [cd,filesep,'Cycle Averages',filesep];
files = extractfield(dir([Cyc_Path,'CycAvg*Sine*.mat']),'name');
txt_var = {'File','Type','AxisName'};
num_var = {'Frequency','Amplitude','MaxVel','MaxVel_sd','MedMaxVel'};
ax_name = {'LHRH','LARP','RALP','X','Y'};
med_results = [cell2table(cell(length(files)*2,length(txt_var)),'VariableNames',txt_var),...
    array2table(NaN(length(files)*2,length(num_var)),'VariableNames',num_var)];
tic;
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
%writetable(med_results,'MedianCycResults.xlsx')
histogram(med_results.VelDiff,-20:5:140)
xlabel('Error between maximum cycle median and average')
ylabel('Occurances')
%% LHRH 300dps plot
load([Cyc_Path,files{contains(files,'LHRH')&contains(files,'300dps')}])
figure;
h1 = gobjects(6,1);
h1(2) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.lz_cycavg-CycAvg.lz_cycstd,fliplr(CycAvg.lz_cycavg+CycAvg.lz_cycstd)],[1,0,0],'FaceAlpha',0.3,'EdgeColor','r');
hold on
h1(5) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.rz_cycavg-CycAvg.rz_cycstd,fliplr(CycAvg.rz_cycavg+CycAvg.rz_cycstd)],[1,0,1],'FaceAlpha',0.3,'EdgeColor','m');
h1(1) = plot(CycAvg.t,CycAvg.lz_cycavg,'r','LineWidth',2);
h1(4) = plot(CycAvg.t,CycAvg.rz_cycavg,'m','LineWidth',2);
h1(3) = plot(CycAvg.t,CycAvg.cycle_params.lz_cycmed,'r:','LineWidth',2);
h1(6) = plot(CycAvg.t,CycAvg.cycle_params.rz_cycmed,'m:','LineWidth',2);
hold off
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
leg = legend(h1,{'LEye Mean','LEye SD','LEye Median','REye Mean','REye SD','REye Median'},'NumColumns',2,'Location','northwest');
title(leg,'Horizontal Response')
leg.ItemTokenSize(1) = 20;
%% LHRH 0.1Hz plot
load([Cyc_Path,files{contains(files,'LHRH')&contains(files,'0.1Hz')}])
figure;
h1 = gobjects(6,1);
hold on
h1(2) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.lz_cycavg-CycAvg.lz_cycstd,fliplr(CycAvg.lz_cycavg+CycAvg.lz_cycstd)],[1,0,0],'FaceAlpha',0.3,'EdgeColor','r');
h1(5) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.rz_cycavg-CycAvg.rz_cycstd,fliplr(CycAvg.rz_cycavg+CycAvg.rz_cycstd)],[1,0,1],'FaceAlpha',0.3,'EdgeColor','m');
h1(1) = plot(CycAvg.t,CycAvg.lz_cycavg,'r','LineWidth',2);
h1(4) = plot(CycAvg.t,CycAvg.rz_cycavg,'m','LineWidth',2);
h1(3) = plot(CycAvg.t,CycAvg.cycle_params.lz_cycmed,'r:','LineWidth',2);
h1(6) = plot(CycAvg.t,CycAvg.cycle_params.rz_cycmed,'m:','LineWidth',2);
hold off
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
set(gca,'YLim',[-135 85])
leg = legend(h1,{'LEye Mean','LEye SD','LEye Median','REye Mean','REye SD','REye Median'},'NumColumns',2,'Location','southwest');
title(leg,'Horizontal Response')
leg.ItemTokenSize(1) = 20;
%% LARP 200dps
load([Cyc_Path,files{contains(files,'LARP')&contains(files,'200dps')}])
figure;
h1 = gobjects(8,1);
hold on
h1(4) = plot(NaN,NaN,'Color',[0 0.5 0],'LineWidth',0.5);
h1(8) = plot(NaN,NaN,'g','LineWidth',0.5);
plot(CycAvg.t,CycAvg.Data_allcyc.LE_Vel_LARP,'Color',[0 0.5 0],'LineWidth',0.5)
plot(CycAvg.t,CycAvg.Data_allcyc.RE_Vel_LARP,'g','LineWidth',0.5)
h1(2) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.ll_cycavg-CycAvg.ll_cycstd,fliplr(CycAvg.ll_cycavg+CycAvg.ll_cycstd)],[0,0.5,0],'FaceAlpha',0.3,'EdgeColor',[0,0.5,0]);
h1(6) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.rl_cycavg-CycAvg.rl_cycstd,fliplr(CycAvg.rl_cycavg+CycAvg.rl_cycstd)],[0,1,0],'FaceAlpha',0.3,'EdgeColor',[0,1,0]);
h1(1) = plot(CycAvg.t,CycAvg.ll_cycavg,'Color',[0,0.5,0],'LineWidth',2);
h1(5) = plot(CycAvg.t,CycAvg.rl_cycavg,'g','LineWidth',2);
h1(3) = plot(CycAvg.t,CycAvg.cycle_params.ll_cycmed,':','Color',[0,0.5,0],'LineWidth',2);
h1(7) = plot(CycAvg.t,CycAvg.cycle_params.rl_cycmed,'g:','LineWidth',2);
hold off
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
set(gca,'YLim',[-120 180])
leg = legend(h1,{'LEye Mean','LEye SD','LEye Median','LEye FiltCycles','REye Mean','REye SD','REye Median','REye FiltCycles'},'NumColumns',2,'Location','northwest');
title(leg,'LARP Response')
leg.ItemTokenSize(1) = 20;
%% All Cycle Averages
load('20240103_113902_eeVORCycParam.mat')
load('VNELcolors.mat','colors')
ax_name = {'LHRH','LARP','RALP','X','Y'};
cyc_params(~contains(cyc_params(:,1),'Sine'),:) = [];
%Hardcoded
fig_ord = [3,5,1,2,4,24,25,26,27,28;8,10,6,7,9,17,19,15,16,18;13,14,11,12,NaN,22,23,20,21,NaN];
fig = figure('Units','inches','Position',[0.5 0.5 14 6],'Color','w');
%%
k=1;
for i = 1:3
    for j = 1:10
        subplot(3,10,k)
        if ~isnan(fig_ord(i,j))
            ax1 = ax_name{cellfun(@(x) contains(cyc_params{fig_ord(i,j),1},x),ax_name)};
            ax = lower(strrep(strrep(strrep(ax1,'LHRH','Z'),'LARP','L'),'RALP','R'));
            amp = extract(cyc_params{fig_ord(i,j),1},digitsPattern+"dps");
            amp = amp{:};
            freq = extract(cyc_params{fig_ord(i,j),1},(("0."+digitsPattern)|digitsPattern)+"Hz");
            freq = freq{:};
            colorL = colors.(['l_',ax]);
            colorR = colors.(['r_',ax]);
            CycAvg = cyc_params{fig_ord(i,j),2};
            hold on
            h1(2) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.(['l',ax,'_cycavg'])-CycAvg.(['l',ax,'_cycstd']),fliplr(CycAvg.(['l',ax,'_cycavg'])+CycAvg.(['l',ax,'_cycstd']))],colorL,'FaceAlpha',0.3,'EdgeColor','none');
            h1(5) = fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.(['r',ax,'_cycavg'])-CycAvg.(['r',ax,'_cycstd']),fliplr(CycAvg.(['r',ax,'_cycavg'])+CycAvg.(['r',ax,'_cycstd']))],colorR,'FaceAlpha',0.3,'EdgeColor','none');
            h1(1) = plot(CycAvg.t,CycAvg.(['l',ax,'_cycavg']),'Color',colorL,'LineWidth',2);
            h1(4) = plot(CycAvg.t,CycAvg.(['r',ax,'_cycavg']),'Color',colorR,'LineWidth',2);
            h1(3) = plot(CycAvg.t,CycAvg.(['l',ax,'_cycmed']),':','Color',colorL,'LineWidth',2);
            h1(6) = plot(CycAvg.t,CycAvg.(['r',ax,'_cycmed']),':','Color',colorR,'LineWidth',2);
            hold off
            title([ax1,' ',freq,' ',amp])
        end
        xlabel('Time (s)')
        ylabel('Angular Velocity (deg/s)')
        set(gca,'FontSize',8)
        k = k+1;
    end
end
%% Seperate figure for the legends
%Standard VNEL color figures + legend with the mean, SD and median
figure;
trac = {'z','l','r','x','y'};
titles = {'LHRH','LARP','RALP','X','Y'};
hold on
for i = 1:length(trac)
    plot([2*i-2,2*i-0.5],[0.5 0.5],'Color',colors.(['l_',trac{i}]),'LineWidth',2);
    plot([2*i-2,2*i-0.5],[-0.5 -0.5],'Color',colors.(['r_',trac{i}]),'LineWidth',2);
    text(2*i-1.25,1,titles{i},'HorizontalAlignment','center','VerticalAlignment','middle')
end
h(1) = plot(NaN,NaN,'k','LineWidth',2);
h(2) = fill([-100,-99,-99,-100],[0 0 1 1],'k','FaceAlpha',0.3,'EdgeColor','none');
h(3) = plot(NaN,NaN,'k:','LineWidth',2);
hold off
set(gca,'XLim',[-5 15],'Ylim',[-10 10])
text(-0.2,-0.5,'Right Eye','HorizontalAlignment','right','VerticalAlignment','middle')
text(-0.2,0.5,'Left Eye','HorizontalAlignment','right','VerticalAlignment','middle')
legend(h,{'Mean','±SD','Median'},'Location','east')