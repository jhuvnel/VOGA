function plotActivationElectrodeResponses(params)
sphere_plot_label_offset_from_sides = 0.1; %May need to be adjusted so that labels are visible
if nargin < 1
    Path = cd;
    Cyc_Path = [cd,filesep,'Cycle Averages'];
else
    Path = params.Path;
    Cyc_Path = params.Cyc_Path;
end
load('VNELcolors.mat','colors')
tab_file = dir([cd,filesep,'*Results.mat']);
if isempty(tab_file)
    error('No Results.mat file was found in this directory.')
end
load([cd,filesep,tab_file(end).name],'all_results')
%Sorted to already be in the correct order
all_results = all_results(contains(all_results.Condition,'Autoscan'),:);
temp_e = split(all_results.Electrode,'E');
all_results.Enum = str2double(temp_e(:,2));
all_results = sortrows(sortrows(sortrows(all_results,'CurrentAmp','ascend'),'PhaseDur','ascend'),'Enum','ascend');
elec_phase = strcat(all_results.Electrode,{' '},strrep(cellstr(num2str(all_results.PhaseDur)),' ',''),'us');
elec_opts = unique(elec_phase,'stable');
%Select Files to Plot
canal = listdlg('PromptString','Plot which electrodes and phase durations?',...
    'ListString',elec_opts,'SelectionMode','multiple');
if isempty(canal)
    return;
end
sel_el = elec_opts(canal);
rel_results = all_results(ismember(elec_phase,sel_el) & all_results.PulseFreq == 200,:); % only plot 200 pps scans
elec_phase = strcat(rel_results.Electrode,{' '},num2str(rel_results.('PhaseDur')),'us');
rel_ep = unique(elec_phase,'stable');
% remove duplicate autoscan levels if scan was repeated
for i = 1:length(rel_ep)
    % check for repeat scans and keep the most recent one only
    % also consider if only some current levels were repeated
    currents = rel_results.CurrentAmp(ismember(elec_phase, rel_ep{i}));
    rel_currents = unique(currents);
    if length(unique(rel_results.Date(ismember(elec_phase, rel_ep{i})))) ~= 1
        for j = 1:length(rel_currents)
            if sum(currents == rel_currents(j)) > 1 % only throw out this current amplitude level if it was repeated
                rel_results(ismember(elec_phase, rel_ep{i}) & ...
                    rel_results.Date ~= max(rel_results.Date(ismember(elec_phase, rel_ep{i}))) ...
                    & rel_results.CurrentAmp == rel_currents(j), :) = [];
                elec_phase = strcat(rel_results.Electrode,{' '},num2str(rel_results.('PhaseDur')),'us');
            end
        end
    end
end

% Make the figure name
file_parts = [rel_results.Subject,rel_results.Visit,cellstr(datetime(rel_results.Date,'Format','yyyyMMdd')),...
    rel_results.Condition,rel_results.Goggle,strcat(strrep(cellstr(num2str(rel_results.('PulseFreq'))),' ',''),'pps')];
fig_title = {strjoin([join(unique(file_parts,'stable')),strjoin(elec_opts(canal),'_')])};

n_row = length(rel_ep);
rel_ep_inds = false(size(rel_results,1),n_row);
rel_ep_s  = NaN(1,n_row);
for i = 1:n_row
    rel_ep_inds(:,i) = contains(elec_phase,rel_ep{i});
    rel_ep_s(i) = find(rel_ep_inds(:,i)==1,1,'first');
end
n_col = sum(rel_ep_inds);
f_order = rel_results.File;
curr_lab = cellstr(num2str(rel_results.('CurrentAmp')));
curr_lab(rel_ep_s) = strcat(curr_lab(rel_ep_s),'\muA');
cyc_lab = cellstr(num2str(rel_results.Cycles));
cyc_lab(rel_ep_s) = strcat('n=',cyc_lab(rel_ep_s));
LRZtrac = {'MaxVel_LL','MaxVel_LR','MaxVel_LZ','MaxVel_RL','MaxVel_RR','MaxVel_RZ'};
YMax = 10*ceil(max(abs(reshape(rel_results{:,LRZtrac}+rel_results{:,strcat(LRZtrac,'_sd')},[],1)))/10);
%% Plot Current Levels
fig = figure;
set(fig,'Color',[1,1,1],'Units','inches','Position',[1 1 5 3])
ha = gobjects(1,length(f_order));
%Set params
YLim = YMax*[-1 1];
x_min = 0.05;
x_max = 0.92;
space_x = 0.01;
y_min = 0.08;
y_max = 0.98;
space_y = 0.03;
%Calculate
x_wid = (x_max - x_min - space_x*(n_col-1))./n_col;
y_wid = (y_max - y_min - space_y*(n_row-1))/n_row;
y_pos = fliplr(y_min:(y_wid+space_y):y_max);
fig_x_wid = NaN(1,length(f_order));
fig_row_pos = NaN(1,length(f_order));
fig_col_pos = NaN(1,length(f_order));
k = 0;
for i = 1:n_row
    x_pos = x_min:(x_wid(i)+space_x):x_max;
    fig_row_pos((1:n_col(i))+k) = x_pos;
    fig_col_pos((1:n_col(i))+k) = y_pos(i);
    fig_x_wid((1:n_col(i))+k) = x_wid(i);
    k = k+n_col(i);
end
for i = 1:length(f_order)
    ha(i) = subplot(n_row,max(n_col),i);
    %Load and plot
    b = load([Cyc_Path,filesep,f_order{i}]);
    a = fieldnames(b);
    CycAvg = b.(a{1});
    fields = fieldnames(CycAvg);
    if ~ismember('t',fields)
        CycAvg.t = (0:1/CycAvg.Fs:(length(CycAvg.ll_cycavg)-1)/CycAvg.Fs)';
    end
    if length(CycAvg.t) > 1000
        s = round(linspace(1,length(CycAvg.t),1000));
    else
        s = 1:length(CycAvg.t);
    end
    hold on
    %Now add the fills and standard deviations and means
    trac = {'l_z','r_z','l_l','r_l','l_r','r_r'};
    if contains(f_order{i},{'RA','LP'})%RALP
        curr_col = colors.l_r;
        rel_trac = {'l_r','r_r'};
    elseif contains(f_order{i},{'LH','RH'})%LHRH
        curr_col = colors.l_z;
        rel_trac = {'l_z','r_z'};
    elseif contains(f_order{i},{'LA','RP'})%LARP
        curr_col = colors.l_l;
        rel_trac = {'l_l','r_l'};
    end
    for j = 1:length(trac)
        trace = strrep(trac{j},'_','');
        plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s),'Color',colors.(trac{j}))
        plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s),'Color',colors.(trac{j}))
        plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s),'Color',colors.(trac{j}),'LineWidth',2);
    end
    %Plot the intended canal again so that it's in the foreground
    for j = 1:length(rel_trac)
        trace = strrep(rel_trac{j},'_','');
        plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s),'Color',colors.(rel_trac{j}))
        plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s),'Color',colors.(rel_trac{j}))
        plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s),'Color',colors.(rel_trac{j}),'LineWidth',2);
    end
    hold off
    axis([0 0.5 YLim])
    text(0.5,YLim(2)-0.01*diff(YLim),curr_lab{i},'Color',curr_col,'HorizontalAlignment','right','VerticalAlignment','top')
    text(0.5,YLim(1)+0.01*diff(YLim),cyc_lab{i},'Color','k','HorizontalAlignment','right','VerticalAlignment','bottom')
end
for i = 1:length(f_order)
    set(ha(i),'Position',[fig_row_pos(i),fig_col_pos(i),fig_x_wid(i),y_wid])
    if ~ismember(i,[rel_ep_s-1,length(f_order)])
        annotation('line',(fig_row_pos(i)+fig_x_wid(i)+0.5*space_x)*[1 1],fig_col_pos(i)+[0 y_wid],'LineWidth',1,'LineStyle','--')
    end
end
annotation('line',fig_row_pos(end)+[0 x_wid(end)],y_min-0.01*[1 1],'LineWidth',2)
annotation('line',(x_max+space_x)*[1 1],y_min+[0 YMax*y_wid/(2*YLim(2))],'LineWidth',2)
annotation('textbox','String','0.5s','EdgeColor','none',...
    'Position',[fig_row_pos(end),0,x_wid(end),y_min-0.01],'HorizontalAlignment','right','VerticalAlignment','middle')
annotation('textbox','String',[num2str(YMax),newline,'\circ/s'],'EdgeColor','none',...
    'Position',[x_max+space_x,y_min,1-(x_max+space_x),YMax*y_wid/(2*YLim(2))],'HorizontalAlignment','left','VerticalAlignment','middle')
set(ha,'XColor','none','YColor','none')
sub = unique(rel_results.Subject);
annotation('textarrow',[0.02 0.02],[0.5 0.5],...
    'String',sub{1}(1:6), 'HeadStyle', 'none', 'LineStyle', 'none',...
    'FontSize',20,'color','k','FontWeight','bold','TextRotation',90,...
    'HorizontalAlignment','center');
savefig(fig,[Path,filesep,'Figures',filesep,strrep(fig_title{:},' ','-'),'.fig'])
saveas(fig,[Path,filesep,'Figures',filesep,strrep(fig_title{:},' ','-'),'.svg'])
%% Plot sphere plot
stim_ear = char(extract(all_results.AxisName(1),"L"|"R"));
last_col = cumsum(n_col);
e_name = cell(length(last_col),1);
fig2 = figure;
set(fig2,'Color',[1,1,1],'Units','inches','Position',[6 1 3 3])
for i = 1:length(last_col)
    load([Cyc_Path,filesep,f_order{last_col(i)}],'CycAvg')
    if contains(f_order{last_col(i)},{'RA','LP'})%RALP
        curr_col = [colors.l_r;colors.r_r];
    elseif contains(f_order{last_col(i)},{'LH','RH'})%LHRH
        curr_col = [colors.l_z;colors.r_z];
    elseif contains(f_order{last_col(i)},{'LA','RP'})%LARP
        curr_col = [colors.l_l;colors.r_l];
    end
    fig2 = MakeSpherePlot(CycAvg,fig2,2,0,1,1,curr_col,stim_ear);
    %for i = 1:length(last_col)
    e_name(i) = rel_results.Electrode(last_col(i)); 
end
set(gca,'Position',[-0.25 -0.25 1.5 1.5])
%Annotations
e_name = split(e_name,'E');
ofs = sphere_plot_label_offset_from_sides; 
if strcmp(stim_ear,'L')
    %Horizonal Annotation on the bottom
    annotation('textarrow',[0.65,0.55],0.1*[1 1],'HeadStyle','none','Color',colors.l_z,'LineWidth',1.5,'String',' LE');
    annotation('textarrow',[0.65,0.55],0.06*[1 1],'HeadStyle','none','Color',colors.r_z,'LineWidth',1.5,'String',' RE');
    annotation('textbox',[0.52 0.11 0.45 0.1],'String',['E',e_name{contains(e_name(:,1),'H'),2},' (LH)'],...
        'EdgeColor','none','Color',colors.l_z,'HorizontalAlignment','left');
    %Anterior Annotation on the top left
    annotation('textarrow',ofs+[0.1,0],0.92*[1 1],'HeadStyle','none','Color',colors.l_l,'LineWidth',1.5,'String',' LE');
    annotation('textarrow',ofs+[0.1,0],0.88*[1 1],'HeadStyle','none','Color',colors.r_l,'LineWidth',1.5,'String',' RE');
    annotation('textbox',[ofs-0.02 0.92 0.45 0.1],'String',['E',e_name{contains(e_name(:,1),'A'),2},' (LA)'],...
        'EdgeColor','none','Color',colors.l_l,'HorizontalAlignment','left');
    %Posterior Annotation on the top left
    annotation('textarrow',1-([0,0.1]+ofs),0.92*[1 1],'HeadStyle','none','Color',colors.l_r,'LineWidth',1.5,'String',' LE');
    annotation('textarrow',1-([0,0.1]+ofs),0.88*[1 1],'HeadStyle','none','Color',colors.r_r,'LineWidth',1.5,'String',' RE');
    annotation('textbox',[0.88-ofs 0.92 0.45 0.1],'String',['E',e_name{contains(e_name(:,1),'P'),2},' (LP)'],...
        'EdgeColor','none','Color',colors.l_r,'HorizontalAlignment','left');
else
    %Horizonal Annotation on the top
    annotation('textarrow',[0.65,0.55],0.92*[1 1],'HeadStyle','none','Color',colors.l_z,'LineWidth',1.5,'String',' LE');
    annotation('textarrow',[0.65,0.55],0.88*[1 1],'HeadStyle','none','Color',colors.r_z,'LineWidth',1.5,'String',' RE');
    annotation('textbox',[0.52 0.92 0.45 0.1],'String',['E',e_name{contains(e_name(:,1),'H'),2},' (RH)'],...
        'EdgeColor','none','Color',colors.l_z,'HorizontalAlignment','left');
    %Posterior Annotation on the bottom left
    annotation('textarrow',[0.15,0.05],0.1*[1 1],'HeadStyle','none','Color',colors.l_l,'LineWidth',1.5,'String',' LE');
    annotation('textarrow',[0.15,0.05],0.06*[1 1],'HeadStyle','none','Color',colors.r_l,'LineWidth',1.5,'String',' RE');
    annotation('textbox',[0.02 0.11 0.45 0.1],'String',['E',e_name{contains(e_name(:,1),'P'),2},' (RP)'],...
        'EdgeColor','none','Color',colors.l_l,'HorizontalAlignment','left');
    %Anterior Annotation on the bottom right
    annotation('textarrow',[0.90,0.80],0.1*[1 1],'HeadStyle','none','Color',colors.l_r,'LineWidth',1.5,'String',' LE');
    annotation('textarrow',[0.90,0.80],0.06*[1 1],'HeadStyle','none','Color',colors.r_r,'LineWidth',1.5,'String',' RE');
    annotation('textbox',[0.78 0.11 0.45 0.1],'String',['E',e_name{contains(e_name(:,1),'A'),2},' (RA)'],...
        'EdgeColor','none','Color',colors.l_r,'HorizontalAlignment','left');
end
savefig(fig2,[Path,filesep,'Figures',filesep,strrep(fig_title{:},' ','-'),'SpherePlot.fig'])
saveas(fig2,[Path,filesep,'Figures',filesep,strrep(fig_title{:},' ','-'),'SpherePlot.svg'])