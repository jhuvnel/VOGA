function fig = plotTabParam(plot_info,params)
%Load inputs
fig_name = plot_info.Name{1};
sub_t = plot_info.SubNames{1};
rel_tab = plot_info.Tables{1};
YLim = plot_info.YLim{1};
YLabel = plot_info.YLabel{1};
XTick = plot_info.XTick{1};
XTickLab = plot_info.XTickLab{1};
XLab = plot_info.XLabel{1};
XScale = plot_info.XScale{1};
rel_leg = plot_info.Legend{1};
Path = params.Path;
n_row = size(rel_tab,1);
n_col = size(rel_tab,2);
code_name = ['Plotting Scripts',filesep,'plotTabParam.m'];
version = params.version;
Experimenter = params.Experimenter;
%Set figure and axes sizing
x_spac = 0.01;
x_min = 0.06;
x_max = 1-x_spac;
y_min = 0.07;
y_max = 0.92;
y_spac = 0.02;
y_wid = (y_max-y_min-y_spac*(n_row-1))/n_row;
x_wid = (x_max-x_min-x_spac*(n_col-1))/n_col;
x_pos = x_min:(x_wid+x_spac):x_max;
y_pos = fliplr(y_min:(y_wid+y_spac):y_max);
fig = figure('Units','inches','Position',[0.5,0.5,11,7],'Color',[1,1,1]);
annotation('textbox',[0 .9 1-0.09 .1],'String',fig_name,'FontSize',14,...
    'HorizontalAlignment','center','EdgeColor','none');
annotation('textbox',[0 0 1 1],'String',[Path,newline,code_name,newline,...
    'VOGA',version,newline,Experimenter],'FontSize',5,...
    'EdgeColor','none','interpreter','none');
ha = gobjects(n_row,n_col);
leg_loc = 'northeast';
if contains(fig_name,'Autoscan')
    leg_loc = 'northwest';
end
for j = 1:n_col
    for i = 1:n_row
        ind = sub2ind([n_col,2],j,i);        
        ha(i,j) = subplot(n_row,n_col,ind);
        hold on
        sub_tab = rel_tab{i,j}; 
        sub_leg = rel_leg{i,j};
        for ii = 1:size(sub_tab,1)
            errorbar(sub_tab.X{ii},sub_tab.Y{ii},sub_tab.Y_sd{ii},...
                'Color',sub_tab.Color{ii},'LineStyle',sub_tab.LineStyle{ii},'Marker',sub_tab.Marker{ii},'MarkerFaceColor','w');
        end
        if ~isempty(sub_leg)
            h = gobjects(size(sub_leg,1),1);
            for jj = 1:size(sub_leg,1)
                h(jj) = plot(NaN,NaN,'Color',sub_leg.Color{jj},'LineStyle',sub_leg.LineStyle{jj},'Marker',sub_leg.Marker{jj},'LineWidth',1.5);
            end
            leg = legend(ha(i,j),h,sub_leg.Name,'Location',leg_loc,'NumColumns',1,'box','off');
            leg.ItemTokenSize(1) = 7;
        end
        hold off
    end
    title(ha(1,j),sub_t{j},'FontSize',14,'Color',sub_tab.Color{ii})
end
for j = 1:n_col
    for i = 1:n_row
        ha(i,j).Position = [x_pos(j),y_pos(i),x_wid,y_wid];
    end
end
for i = 1:n_row
    ylabel(ha(i,1),YLabel{i},'FontSize',12)
    set(ha(i,:),'YLim',YLim(i,:))
end
xlabel(ha(2,:),XLab,'FontSize',12)
set(ha,'XTick',XTick,'XTickLabel',XTickLab,'xscale',XScale,'box','on',...
    'XTickLabelRotation',0,'xminortick','off')
if contains(XScale,'log')
    set(ha,'XLim',[0.8*XTick(1),1.2*XTick(end)])
else
    buff = diff(XTick([1,end]))/length(XTick)/2;
    set(ha,'XLim',[XTick(1)-buff,XTick(end)+buff])
end
set(ha(:,2:n_col),'YTickLabel',[])
set(ha(1,:),'XTickLabel',[])
end