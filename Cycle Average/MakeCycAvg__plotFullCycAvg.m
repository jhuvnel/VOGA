%% Make Full Plot Cyc Avg
function ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim_Pos,YLim_Vel,traces_pos,traces_vel,CycAvg)
%% Initialize Figure
if isempty(ha) %first time running
    ha = gobjects(5,1);
    if isfield(CycAvg.Data_rawvel,'t_cond')
        XLim_Long = CycAvg.Data_rawvel.t_cond([1,end]);
    else
        XLim_Long = CycAvg.Data_rawpos.te([1,end]);
    end
    XLim_Short = CycAvg.t([1,end]);
    ha(1) = subplot(3,3,[1 2 3]);
    ha(2) = subplot(3,3,[4 5 6]);
    ha(3) = subplot(3,3,7);
    ha(4) = subplot(3,3,8);
    ha(5) = subplot(3,3,9);
    %Link Axes
    linkaxes(ha([1,2]),'x')
    linkaxes(ha([2,3]),'y')
    linkaxes(ha([3,4,5]),'xy')
    %Set Labels
    set(ha(1:2),'XLim',XLim_Long)
    set(ha(3:5),'XLim',XLim_Short)
    set(ha(4:5),'YTickLabel',[])
    ylabel(ha(1),'Angular Position (deg)','FontWeight','bold')
    ylabel(ha([2,3]),'Angular Velocity (dps)','FontWeight','bold')
    xlabel(ha([2,4]),'Time (s)')
    title(ha(3),'All Filtered Cycles')
    title(ha(5),'Cycle Averages')
    if ismember(type,[1,4]) %Position, Velocity, and Cycle Averaging
        ha(1).Position = [0.06 0.67 0.92 0.27];
        ha(2).Position = [0.06 0.39 0.92 0.27];
        ha(3).Position = [0.06 0.05 0.30 0.27];
        ha(4).Position = [0.37 0.05 0.30 0.27];
        ha(5).Position = [0.68 0.05 0.30 0.27];
    elseif ismember(type,2) %Position and Velocity, no Cycle Averaging
        set(ha(3:5),'Visible','off')
        ha(1).Position = [0.06 0.52 0.92 0.43];
        ha(2).Position = [0.06 0.05 0.92 0.43];
    elseif ismember(type,3) %Velocity and Cycle Averaging, no Position
        ha(1).Position = [0.06 0.55 0.92 0.40];
        set(ha(1),'YTick',[])
        ylabel(ha(1),'')
        ha(2).Position = [0.06 0.55 0.92 0.40];
        ha(3).Position = [0.06 0.05 0.30 0.42];
        ha(4).Position = [0.37 0.05 0.30 0.42];
        ha(5).Position = [0.68 0.05 0.30 0.42];
    end
end
%% Update each time
if ismember(type,[3,4])
    te = CycAvg.Data_rawvel.t_cond;
    ts = CycAvg.Data_rawvel.t_cond;
    stim = CycAvg.Data_rawvel.stim_cond;
    keep_inds= reshape(1:size(CycAvg.Data_allcyc.keep_inds,1)*size(CycAvg.Data_allcyc.keep_inds,2),size(CycAvg.Data_allcyc.keep_inds,1),[]);
else
    te = CycAvg.Data_rawpos.te;
    ts = CycAvg.Data_rawpos.ts;
    stim = CycAvg.Data_rawvel.stim;
    keep_inds = CycAvg.Data_allcyc.keep_inds;
end
t_s = CycAvg.t;
stims = CycAvg.stim;
keep_tr = CycAvg.keep_tr;
detec_tr = CycAvg.detec_tr;
fill_color = repmat({colors.cyc_rm},size(keep_inds,2),1);
fill_color(keep_tr) = {colors.cyc_keep};
set(ha(1),'YLim',YLim_Pos)
set(ha(2:5),'YLim',YLim_Vel)
title(ha(4),['Accepted Cycles: ',num2str(sum(keep_tr)),'/',num2str(length(keep_tr))])
if get(ha(1),'Visible')==1
    set(ha(1),'XTick',te(round(floor(mean(keep_inds)),0)),'XTickLabel',1:length(keep_tr),'Xaxislocation','top')
end
%Subfigure 1: Plot Raw/Filtered Position Data
axes(ha(1))
cla;
h = gobjects(1,length(traces_pos)+1);
hold on
if ismember(type,[1,3,4]) %Has cycles
    for j = 1:size(keep_inds,2)
        fill([te(keep_inds(1,j)),te(keep_inds(end,j)),te(keep_inds(end,j)),te(keep_inds(1,j))]',[500,500,-500,-500]',fill_color{j});
    end
    if ~isempty(detec_tr)
        for i = 1:length(detec_tr)
            plot(te(keep_inds([1,end],detec_tr(i))),[YLim_Pos(2),YLim_Pos(2)],'-b','LineWidth',2)
        end
    end
end
h(1) = plot(ts,stim,'k','LineWidth',line_wid.norm);
% Raw Traces
for i = 1:length(traces_pos)
    if isfield(CycAvg.Data_rawpos,[traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end),'_cond'])
        plot(te,CycAvg.Data_rawpos.([traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end),'_cond']),'Color',colors.([lower(traces_pos{i}(1)),'_',lower(traces_pos{i}(2)),'_s']),'LineWidth',line_wid.norm)
    elseif isfield(CycAvg.Data_rawpos,[traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end)])
        plot(te,CycAvg.Data_rawpos.([traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end)]),'Color',colors.([lower(traces_pos{i}(1)),'_',lower(traces_pos{i}(2)),'_s']),'LineWidth',line_wid.norm)
    end
end
% Filtered Traces
for i = 1:length(traces_pos)
    if isfield(CycAvg.Data_filtpos,[traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end),'_cond'])
        plot(ts,CycAvg.Data_filtpos.([traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end),'_cond']),'Color',colors.([lower(traces_pos{i}(1)),'_',lower(traces_pos{i}(2))]),'LineWidth',line_wid.norm);
    elseif isfield(CycAvg.Data_filtpos,[traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end)])
        plot(ts,CycAvg.Data_filtpos.([traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end)]),'Color',colors.([lower(traces_pos{i}(1)),'_',lower(traces_pos{i}(2))]),'LineWidth',line_wid.norm);
    end
    h(i+1) = plot(NaN,NaN,'Color',colors.([lower(traces_pos{i}(1)),'_',lower(traces_pos{i}(2))]));
end
hold off
leg1 = legend(h,[{'Stim'};reshape(traces_pos,[],1)],'NumColumns',length(traces_pos)+1);
leg1.ItemTokenSize(1) = 7;

% Subfigure 2: Plot Raw/Filtered Velocity Data
axes(ha(2))
cla;
h1 = gobjects(1,length(traces_vel)+1);
hold on
if ismember(type,[1,3,4]) %Has cycles
    for j = 1:size(keep_inds,2)
        fill([te(keep_inds(1,j)),te(keep_inds(end,j)),te(keep_inds(end,j)),te(keep_inds(1,j))]',[500,500,-500,-500]',fill_color{j});
    end
    linemark_type = '-';
elseif ismember(type,2) %Long files
    linemark_type = '.';
end
h1(1) = plot(ts,stim,'k','LineWidth',line_wid.norm);
% Raw Traces
for i = 1:length(traces_vel)
    if isfield(CycAvg.Data_rawvel,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end),'_cond'])
        plot(ts,CycAvg.Data_rawvel.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end),'_cond']),linemark_type,'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2)),'_s']),'LineWidth',line_wid.norm)
    elseif isfield(CycAvg.Data_rawvel,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
        plot(ts,CycAvg.Data_rawvel.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)]),linemark_type,'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2)),'_s']),'LineWidth',line_wid.norm)
    end
end
% Filtered Traces
for i = 1:length(traces_vel)
    if isfield(CycAvg.Data_filtvel,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end),'_cond'])
        plot(ts,CycAvg.Data_filtvel.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end),'_cond']),linemark_type,'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm);
    elseif isfield(CycAvg.Data_filtvel,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
        plot(ts,CycAvg.Data_filtvel.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)]),linemark_type,'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm);
    end
    h1(i+1) = plot(NaN,NaN,linemark_type,'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]));
end
hold off
leg2 = legend(h1,[{'Stim'};reshape(traces_vel,[],1)],'NumColumns',length(traces_vel)+1);
leg2.ItemTokenSize(1) = 7;

% Subfigure 3: All Cycle Aligned Traces
axes(ha(3))
cla;
if get(ha(3),'Visible')
    hold on
    plot(t_s,CycAvg.Data_allcyc.stim,'k','LineWidth',line_wid.norm)
    %All Filtered Velocity Data
    for i = 1:length(traces_vel)
        if isfield(CycAvg.Data_allcyc,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
            plot(t_s,CycAvg.Data_rawvel.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])(CycAvg.Data_allcyc.keep_inds),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2)),'_s']),'LineWidth',line_wid.norm)
        end
    end
    for i = 1:length(traces_vel)
        if isfield(CycAvg.Data_allcyc,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
            plot(t_s,CycAvg.Data_allcyc.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)]),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm)
        end
    end
    hold off
end
% Subfigure 4: Accepted Cycle Aligned Traces
axes(ha(4))
cla;
if get(ha(4),'Visible')
    hold on
    if isfield(CycAvg,'stim_cyc') %Mutiple head traces to show
        plot(t_s,CycAvg.stim_cyc,'k','LineWidth',line_wid.norm)
    else
        plot(t_s,stims,'k','LineWidth',line_wid.norm)
    end
    %Only Accepted Trace Filtered Velocity Data
    for i = 1:length(traces_vel)
        trac = [traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)];
        if isfield(CycAvg.Data_allcyc,[trac,'_smooth'])
            plot(t_s,CycAvg.Data_allcyc.([trac,'_smooth'])(:,keep_tr),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2)),'_s']),'LineWidth',line_wid.norm)
            plot(t_s,CycAvg.Data_allcyc.(trac)(:,keep_tr),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm)
        elseif isfield(CycAvg.Data_allcyc,trac)
            plot(t_s,CycAvg.Data_allcyc.(trac)(:,keep_tr),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm)
        end
        if isfield(CycAvg,[lower(traces_vel{i}(1:2)),'_saccade_time'])
            plot(CycAvg.([lower(traces_vel{i}(1:2)),'_saccade_time']),(0.99*diff(YLim_Vel)+YLim_Vel(1))*ones(length(CycAvg.([lower(traces_vel{i}(1:2)),'_saccade_time'])),1),'*','MarkerSize',5,'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2)),'_s']))
        end
    end
    hold off
end
% Subfigure 5: Cycle Averages
axes(ha(5))
cla;
if get(ha(5),'Visible')
    s = unique(round(linspace(1,length(CycAvg.t),1000)));
    hold on
    if isfield(CycAvg,'stim_cyc') %Mutiple head traces to show
        trace = 'stim';
        fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[(CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s)),fliplr((CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s)))],0.5*[1,1,1])
        plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s),'Color','k')
        plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s),'Color','k')
        plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s),'Color','k','LineWidth',line_wid.bold);
    else %Only 1
        plot(CycAvg.t(s),CycAvg.stim(s),'k','LineWidth',line_wid.bold);
    end
    %Only Selected Filtered Velocity Data
    for i = 1:length(traces_vel)
        trace = lower(traces_vel{i}(1:2));
        if isfield(CycAvg,[trace,'_cyc_fit'])
            plot(CycAvg.t(s),CycAvg.([trace,'_cyc'])(:,s),'Color',colors.([trace(1),'_',trace(2),'_s']),'LineWidth',0.5);
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[(CycAvg.([trace,'_cycavg_fit'])(s) - CycAvg.([trace,'_cycstd_fit'])(s)),fliplr((CycAvg.([trace,'_cycavg_fit'])(s) + CycAvg.([trace,'_cycstd_fit'])(s)))],colors.([trace(1),'_',trace(2),'_s']))
            plot(CycAvg.t(s),CycAvg.([trace,'_cycavg_fit'])(s) + CycAvg.([trace,'_cycstd_fit'])(s),'Color',colors.([trace(1),'_',trace(2)]))
            plot(CycAvg.t(s),CycAvg.([trace,'_cycavg_fit'])(s) - CycAvg.([trace,'_cycstd_fit'])(s),'Color',colors.([trace(1),'_',trace(2)]))
            plot(CycAvg.t(s),CycAvg.([trace,'_cycavg_fit'])(s),'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',line_wid.bold);
        elseif isfield(CycAvg,[trace,'_cycavg'])
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[(CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s)),fliplr((CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s)))],colors.([trace(1),'_',trace(2),'_s']))
            plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s),'Color',colors.([trace(1),'_',trace(2)]))
            plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s),'Color',colors.([trace(1),'_',trace(2)]))
            plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s),'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',line_wid.bold);
        end
    end
    hold off
end
end