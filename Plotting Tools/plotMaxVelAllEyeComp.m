function fig = plotMaxVelAllEyeComp(plot_info,params)
%Load inputs
load('VNELcolors.mat','colors')
fig_name = plot_info.Name{1};
sub_t = plot_info.SubNames{1};
rel_tab = plot_info.Tables{1};
YLim = plot_info.YLim{1};
XTick = plot_info.XTick{1};
XTickLab = plot_info.XTickLab{1};
x_ax = plot_info.XVar{1};
XLab = plot_info.XLabel{1};
XScale = plot_info.XScale{1};
Path = params.Path;
fnum = size(sub_t,2);
code_name = ['Plotting Scripts',filesep,'plotMaxVelAllEyeComp.m'];
version = params.version;
Experimenter = params.Experimenter;
%Plot either LRX or XYZ and put the relevant axis last so it's on top
canal_ax = {{'LA','RP'},{'RA','LP'},{'LH','RH'},'X','Y'};
canal_tr = {{'ll';'rl'};{'lr';'rr'};{'lz';'rz'};{'lx';'rx'};{'ly';'ry'}};
flat_can_tr = vertcat(canal_tr{:});
xyz_lrz_tr = {{'LA','RP','RA','LP','LH','RH'},{'X','Y'}};
tr_ord = {{'ll';'rl';'lr';'rr';'lz';'rz'};{'lx';'rx';'ly';'ry';'lz';'rz'}};
%Set figure and axes sizing
x_spac = 0.01;
x_min = 0.06;
x_max = 1-x_spac;
y_min = 0.07;
y_max = 0.87;
y_spac = 0.05;
y_wid = (y_max-y_min-y_spac)/2;
x_wid = (x_max-x_min-x_spac*(fnum-1))/fnum;
x_pos = x_min:(x_wid+x_spac):x_max;
y_pos = fliplr(y_min:(y_wid+y_spac):y_max);
fig = figure('Units','inches','Position',[0.5,0.5,11,7],'Color',[1,1,1]);
annotation('textbox',[0 .9 1-0.09 .1],'String',fig_name,'FontSize',14,...
    'HorizontalAlignment','center','EdgeColor','none');
annotation('textbox',[0 0 1 1],'String',[Path,newline,code_name,newline,...
    'VOGA',version,newline,Experimenter],'FontSize',5,...
    'EdgeColor','none','interpreter','none');
ha = gobjects(2,fnum);
h1 = gobjects(1,length(flat_can_tr));
for j = 1:fnum
    for i = 1:2
        ind = sub2ind([fnum,2],j,i);        
        ha(i,j) = subplot(2,fnum,ind);
        hold on
        if ~isempty(rel_tab{i,j})
            p_tr = tr_ord{find(cellfun(@(x) contains(sub_t{i,j},x),xyz_lrz_tr),1,'first')};
            rel_tr = canal_tr{find(cellfun(@(x) contains(sub_t{i,j},x),canal_ax),1,'first')};
            sub_tab = rel_tab{i,j};
            for ii = 1:length(p_tr)
                tr = p_tr{ii};
                errorbar(sub_tab.(x_ax),abs(sub_tab.(['MaxVel_',upper(tr)])),sub_tab.(['MaxVel_',upper(tr),'_sd']),'Color',colors.([tr(1),'_',tr(2)]));
                plot(sub_tab.(x_ax),abs(sub_tab.(['MaxVel_',upper(tr)])),'o','Color',colors.([tr(1),'_',tr(2)]));
            end
            for ii = 1:length(rel_tr) %Plot the relevant traces front and center
                tr = rel_tr{ii};
                errorbar(sub_tab.(x_ax),abs(sub_tab.(['MaxVel_',upper(tr)])),sub_tab.(['MaxVel_',upper(tr),'_sd']),'Color',colors.([tr(1),'_',tr(2)]),'LineWidth',1.5);
                plot(sub_tab.(x_ax),abs(sub_tab.(['MaxVel_',upper(tr)])),'o','Color',colors.([tr(1),'_',tr(2)]),'LineWidth',1.5);
            end
        end
        hold off
        title(sub_t{i,j},'FontSize',14,'Color',colors.(['l_',tr(2)]))
    end
end
for j = 1:fnum
    for i = 1:2
        ha(i,j).Position = [x_pos(j),y_pos(i),x_wid,y_wid];
    end
end
hold on
for ii = 1:length(flat_can_tr)
    h1(ii) = plot(NaN,NaN,'Color',colors.([flat_can_tr{ii}(1),'_',flat_can_tr{ii}(2)]),'LineWidth',2);
end
hold off
ylabel(ha(:,1),'Maximum Velocity (dps)','FontSize',12)
xlabel(ha(2,:),XLab,'FontSize',12)
set(ha,'XTick',XTick,'XLim',[min(XTick)*0.8 max(XTick)*1.2],...
    'XTickLabel',XTickLab,'XTickLabelRotation',0)
set(ha,'YLim',YLim,'xscale',XScale,'box','on','xminortick','off')
set(ha(:,2:fnum),'YTickLabel',[])
set(ha(1,:),'XTickLabel',[])
leg_reord = [1:2:9,2:2:10];
leg = legend(ha(1),h1(leg_reord),upper(flat_can_tr(leg_reord)),'NumColumns',2,'box','off');
leg.ItemTokenSize(1) = 7;
leg.Position =  [1-0.08,0.88,0.07,0.11];
end