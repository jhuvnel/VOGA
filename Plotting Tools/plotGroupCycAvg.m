function fig = plotGroupCycAvg(plot_info,params)
%Load inputs
load('VNELcolors.mat','colors')
fig_name = plot_info.Name{1};
rel_files = plot_info.Files{1};
sub_t = plot_info.SubNames{1};
YLim = plot_info.YLim{1};
Path = params.Path;
Cyc_Path = params.Cyc_Path;
fnum = length(rel_files);
code_name = ['Plotting Scripts',filesep,'plotGroupCycAvg.m'];
version = params.version;
Experimenter = params.Experimenter;
%Plot either LRX or XYZ and put the relevant axis last so it's on top
canal_ax = {{'LA','RP'},{'RA','LP'},{'LH','RH'},'X','Y'};
canal_tr = {{'lr';'rr';'lz';'rz';'ll';'rl'};...
            {'ll';'rl';'lz';'rz';'lr';'rr'};...
            {'lr';'rr';'ll';'rl';'lz';'rz'};...
            {'ly';'ry';'lz';'rz';'lx';'rx'};...
            {'lx';'rx';'lz';'rz';'ly';'ry'}};
p_tr = canal_tr{find(cellfun(@(x) contains(fig_name,x),canal_ax),1,'first')};
if contains(fig_name,'Impulse') %Just want the relevant canal for Impulses
    p_tr = p_tr(5:6);
end
%Set figure and axes sizing
fig_len = min([16,4*fnum]); %Make the figure smaller so the aspect ratio never gets too weird
x_spac = 0.01;
x_min = 0.07*8/fig_len;
x_max = 0.99;
y_pos = 0.12;
y_max = 0.82;
y_wid = y_max-y_pos;
x_wid = (x_max-x_min-x_spac*(fnum-1))/fnum;
x_pos = x_min:(x_wid+x_spac):x_max;
fig = figure('Units','inches','Position',[0.25,5.5,fig_len,4],'Color',[1,1,1]);
annotation('textbox',[0 .9 1 .1],'String',fig_name,'FontSize',14,...
    'HorizontalAlignment','center','EdgeColor','none');
annotation('textbox',[0 0 1 1],'String',[Path,newline,code_name,newline,...
    'VOGA',version,newline,Experimenter],'FontSize',5,...
    'EdgeColor','none','interpreter','none');
ha = gobjects(1,fnum);
p_tr_bool = false(1,length(p_tr));
for i = 1:fnum
    ha(i) = subplot(1,fnum,i);
    hold on
    if ~isempty(rel_files{i})
        load([Cyc_Path,filesep,rel_files{i}],'CycAvg')
        if ~isfield(CycAvg,'t')
            CycAvg.t = (0:1/CycAvg.Fs:(length(CycAvg.ll_cycavg)-1)/CycAvg.Fs);
        end
        CycAvg.t = reshape(CycAvg.t,[],1);
        s = unique(round(linspace(1,length(CycAvg.t),1000))); %plot a subset of points (1000) if there are more than 1000 points
        [aa,ab] = size(CycAvg.stim);
        if aa == length(CycAvg.t) && ab ~=1
            CycAvg.stim = mean(CycAvg.stim,2)';
        elseif ab == length(CycAvg.t) && aa ~=1
            CycAvg.stim = mean(CycAvg.stim,1);
        end
        if contains(fig_name,'Sine')
            plot(CycAvg.t(s),-CycAvg.stim(s),'k');
            for ii = 1:length(p_tr)
                tr = p_tr{ii};
                if isfield(CycAvg,[tr,'_cycavg'])
                    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.([tr,'_cycavg'])(s)-CycAvg.([tr,'_cycstd'])(s),fliplr((CycAvg.([tr,'_cycavg'])(s)+CycAvg.([tr,'_cycstd'])(s)))],colors.([tr(1),'_',tr(2),'_s']))
                    plot(CycAvg.t(s),CycAvg.([tr,'_cycavg'])(s) + CycAvg.([tr,'_cycstd'])(s),'Color',colors.([tr(1),'_',tr(2)]))
                    plot(CycAvg.t(s),CycAvg.([tr,'_cycavg'])(s) - CycAvg.([tr,'_cycstd'])(s),'Color',colors.([tr(1),'_',tr(2)]))
                    plot(CycAvg.t(s),CycAvg.([tr,'_cycavg'])(s),'Color',colors.([tr(1),'_',tr(2)]),'LineWidth',2);
                else
                    p_tr_bool(ii) = true;
                end
            end
            text_num = [CycAvg.parameterized.AxisName,cellstr(num2str(round(CycAvg.parameterized.MaxVel,2),2)),cellstr(num2str(round(CycAvg.parameterized.Phase,0),3))];
            disp_text = ['MaxVel ',text_num{1,1},'(dps): ',text_num{1,2},newline,...
                'MaxVel ',text_num{2,1},'(dps): ',text_num{2,2},newline,...
                'Phase(deg): ',text_num{1,3}];
            box_bool = 'on';
            grid_bool = 'on';
        elseif contains(fig_name,'Impulse')
            n = -1;
            if -min(mean(CycAvg.stim))>max(mean(CycAvg.stim))
                n = 1;
            end
            plot(CycAvg.t,-n*CycAvg.stim_cyc,'k');
            for ii = 1:length(p_tr)
                tr = p_tr{ii};
                if isfield(CycAvg,[tr,'_cyc'])
                    plot(CycAvg.t,n*CycAvg.([tr,'_cyc_prefilt']),'Color',[colors.([tr(1),'_',tr(2),'_s']),0.5])
                    plot(CycAvg.t,n*CycAvg.([tr,'_cyc']),'Color',colors.([tr(1),'_',tr(2)]))
                else
                    p_tr_bool(ii) = true;
                end
            end
            disp_text = ['Gain: ',num2str(round(CycAvg.parameterized.Gain,2)),newline,...
                'Lat(ms): ',num2str(round(CycAvg.parameterized.Latency,0))];
            box_bool = 'on';
            grid_bool = 'on';
        elseif contains(fig_name,'Exponential')
            plot(CycAvg.t(s),-CycAvg.stim(s),'k');
            for ii = 1:length(p_tr)
                tr = p_tr{ii};
                if isfield(CycAvg,[tr,'_cyc'])
                    plot(CycAvg.t(s),CycAvg.([tr,'_cyc_prefilt'])(s),'.','Color',colors.([tr(1),'_',tr(2),'_s']))
                    plot(CycAvg.t(s),CycAvg.([tr,'_cycavg_fit'])(s),'Color',colors.([tr(1),'_',tr(2)]),'LineWidth',2)
                else
                    p_tr_bool(ii) = true;
                end
            end
            disp_text = ['MaxVel(dps): ',num2str(round(CycAvg.parameterized.MaxVel(1),2),2),newline,...
                'Tau(s): ',num2str(round(CycAvg.parameterized.Tau(1),2),2)];
            box_bool = 'on';
            grid_bool = 'on';
        elseif contains(fig_name,'Autoscan')
            for ii = 1:length(p_tr)
                tr = p_tr{ii};
                if isfield(CycAvg,[tr,'_cyc'])
                    plot(CycAvg.t,CycAvg.([tr,'_cyc']),'Color',colors.([tr(1),'_',tr(2)]))
                else
                    p_tr_bool(ii) = true;
                end
            end
            disp_text = sub_t{i};
            sub_t{i} = '';      
            box_bool = 'on';
            grid_bool = 'on';
        end        
        XLim = CycAvg.t([1,end]);
        set(gca,'XLim',XLim)
        text(0.99*diff(XLim)+XLim(1),YLim(1),disp_text,...
            'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',7)
    end
    hold off
    title(sub_t{i})
end
p_tr(p_tr_bool) = [];
hold on
h1 = gobjects(1,length(p_tr)+1);
h1(1) = plot(NaN,NaN,'k');
for ii = 1:length(p_tr)
    h1(1+ii) = plot(NaN,NaN,'Color',colors.([p_tr{ii}(1),'_',p_tr{ii}(2)]),'LineWidth',2);
end
hold off
for i = 1:fnum
    ha(i).Position = [x_pos(i),y_pos,x_wid,y_wid];
end
ylabel(ha(1),'Angular Velocity (dps)','FontSize',12)
xlabel(ha,'Time (s)')
set(ha(2:fnum),'YTickLabel',[])
set(ha,'box',box_bool,'YLim',YLim,'XGrid',grid_bool,'YGrid',grid_bool)
leg = legend(ha(1),h1,[{'Stim'};upper(p_tr)],'NumColumns',length(h1),'box','off');
leg.ItemTokenSize(1) = 7;
leg.Position = [0.5-(0.1635/2),y_max+0.05,0.1635,0.0516];
end