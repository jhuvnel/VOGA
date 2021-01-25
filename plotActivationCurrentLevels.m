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
% Load File
code_Path = '/Volumes/MVI/DATA SUMMARY/IN PROGRESS/VOG Analysis Scripts/LDVOG_Neurolign';
path = cd;
Cyc_Path = [path,filesep,'Cycle Averages'];
files = [dir([Cyc_Path,filesep,'*CurrentFitting*.mat']);dir([Cyc_Path,filesep,'*Autoscan*.mat'])];
fnames = unique({files.name}');
%% Plot 3 canals at a time
%Select Files to Plot

canal = listdlg('PromptString','Which canal?','ListString',{'LARP','RALP','LHRH'},'SelectionMode','single');
if isempty(canal) % No constraints on file type   
    indx = listdlg('ListString',fnames,'PromptString','Pick the files to plot in row 1','ListSize',[400 300]);
    row1 = fnames(indx);
    indx = listdlg('ListString',fnames,'PromptString','Pick the files to plot in row 2','ListSize',[400 300]);
    row2 = fnames(indx);
    indx = listdlg('ListString',fnames,'PromptString','Pick the files to plot in row 3','ListSize',[400 300]);
    row3 = fnames(indx);
else %One per canal
    switch canal 
        case 1 %LARP
            row1 = fnames(contains(fnames,'RPE3')|contains(fnames,'LAE9'));
            row2 = fnames(contains(fnames,'RPE4')|contains(fnames,'LAE10'));
            row3 = fnames(contains(fnames,'RPE5')|contains(fnames,'LAE11'));
        case 2 %RALP
            row1 = fnames(contains(fnames,'RAE9')|contains(fnames,'LPE3'));
            row2 = fnames(contains(fnames,'RAE10')|contains(fnames,'LPE4'));
            row3 = fnames(contains(fnames,'RAE11')|contains(fnames,'LPE5'));
        case 3 %LHRH
            row1 = fnames(contains(fnames,'RHE6')|contains(fnames,'LHE6'));
            row2 = fnames(contains(fnames,'RHE7')|contains(fnames,'LHE7'));
            row3 = fnames(contains(fnames,'RHE8')|contains(fnames,'LHE8'));
    end
end

if ~(length(row3)==length(row1)&&length(row1)==length(row2))
    error('Unequal number of files selected for each canal')
end
all_canals = [row1,row2,row3];
n_col = length(row1);
%Determine the order
all_curr = zeros(size(all_canals));
for i = 1:n_col
    for j = 1:3
        file = all_canals{i,j};
        curr = file(strfind(file,'uA')-3:strfind(file,'uA')-1);
        if ~ismember(curr(1),'123456789')
            curr = curr(2:3);
        end
        all_curr(i,j) = str2double(curr);
    end
end
[~,i1] = sort(all_curr(:,1));
[~,i2] = sort(all_curr(:,2));
[~,i3] = sort(all_curr(:,3));
f_order = [row1(i1);row2(i2);row3(i3)];
curr_lab = cellstr(num2str([sort(all_curr(:,1));sort(all_curr(:,2));sort(all_curr(:,3))]));
curr_lab{1} = [curr_lab{1},'\muA'];
curr_lab{1+n_col} = [curr_lab{1+n_col},'\muA'];
curr_lab{1+2*n_col} = [curr_lab{1+2*n_col},'\muA'];
% Plot Current Levels
fig = figure;
fig.Color = [1,1,1];
fig.Units = 'inches';
fig.Position = [1 1 8 4];
ha = gobjects(1,length(f_order));
%Set params
grid_on = true;
yscale = 75;
YLim = yscale*[-1 1];
x_min = 0.01;
x_max = 0.93;
space_x = 0.01;
y_min = 0.08;
y_max = 0.99;
space_y = 0.05;
%Calculate
x_wid = (x_max - x_min - space_x*(n_col-1))/n_col;
fig_row_pos = repmat(x_min:(x_wid+space_x):x_max,1,3);
y_wid = (y_max - y_min - space_y*2)/3;
fig_col_pos = reshape(repmat(fliplr(y_min:(y_wid+space_y):y_max),n_col,1),[],1)';
for i = 1:length(f_order)
    ha(i) = subplot(3,n_col,i);
    set(gca,'XColor','none','YColor','none')
end
for i = 1:length(f_order)
    axes(ha(i))  
    set(ha(i),'Position',[fig_row_pos(i),fig_col_pos(i),x_wid,y_wid]) 
    if mod(i,n_col) > 0
        annotation('line',(fig_row_pos(i)+x_wid+0.5*space_x)*[1 1],fig_col_pos(i)+[0 y_wid],'LineWidth',1,'LineStyle','--') 
    end
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
    %Plot the intended canal again so that it's in the foreground
    if contains(f_order{i},'LP') || contains(f_order{i},'RA') %RALP
        curr_col = colors.l_r;
        %LE-LHRH
        plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
        plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
        plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
        %RE-LHRH
        plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
        plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
        plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
        %LE-LARP
        plot(CycAvg.t(s),CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s),'Color',colors.l_l)
        plot(CycAvg.t(s),CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s),'Color',colors.l_l)
        plot(CycAvg.t(s),CycAvg.ll_cycavg(s),'Color',colors.l_l,'LineWidth',2);
        %RE-LARP
        plot(CycAvg.t(s),CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s),'Color',colors.r_l)
        plot(CycAvg.t(s),CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s),'Color',colors.r_l)
        plot(CycAvg.t(s),CycAvg.rl_cycavg(s),'Color',colors.r_l,'LineWidth',2);
        %LE_RALP
        plot(CycAvg.t(s),CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s),'Color',colors.l_r)
        plot(CycAvg.t(s),CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s),'Color',colors.l_r)
        plot(CycAvg.t(s),CycAvg.lr_cycavg(s),'Color',colors.l_r,'LineWidth',2);
        %RE-RALP
        plot(CycAvg.t(s),CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s),'Color',colors.r_r)
        plot(CycAvg.t(s),CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s),'Color',colors.r_r)
        plot(CycAvg.t(s),CycAvg.rr_cycavg(s),'Color',colors.r_r,'LineWidth',2);
    elseif contains(f_order{i},'LH') || contains(f_order{i},'RH') %LHRH
        curr_col = colors.l_z;
        %LE-LARP
        plot(CycAvg.t(s),CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s),'Color',colors.l_l)
        plot(CycAvg.t(s),CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s),'Color',colors.l_l)
        plot(CycAvg.t(s),CycAvg.ll_cycavg(s),'Color',colors.l_l,'LineWidth',2);
        %RE-LARP
        plot(CycAvg.t(s),CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s),'Color',colors.r_l)
        plot(CycAvg.t(s),CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s),'Color',colors.r_l)
        plot(CycAvg.t(s),CycAvg.rl_cycavg(s),'Color',colors.r_l,'LineWidth',2);
        %LE_RALP
        plot(CycAvg.t(s),CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s),'Color',colors.l_r)
        plot(CycAvg.t(s),CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s),'Color',colors.l_r)
        plot(CycAvg.t(s),CycAvg.lr_cycavg(s),'Color',colors.l_r,'LineWidth',2);
        %RE-RALP
        plot(CycAvg.t(s),CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s),'Color',colors.r_r)
        plot(CycAvg.t(s),CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s),'Color',colors.r_r)
        plot(CycAvg.t(s),CycAvg.rr_cycavg(s),'Color',colors.r_r,'LineWidth',2);
        %LE-LHRH
        plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
        plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
        plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
        %RE-LHRH
        plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
        plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
        plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
    elseif contains(f_order{i},'RP') || contains(f_order{i},'LA') %LARP
        curr_col = colors.l_l;
        %LE-LHRH
        plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
        plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
        plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
        %RE-LHRH
        plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
        plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
        plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
        %LE_RALP
        plot(CycAvg.t(s),CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s),'Color',colors.l_r)
        plot(CycAvg.t(s),CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s),'Color',colors.l_r)
        plot(CycAvg.t(s),CycAvg.lr_cycavg(s),'Color',colors.l_r,'LineWidth',2);
        %RE-RALP
        plot(CycAvg.t(s),CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s),'Color',colors.r_r)
        plot(CycAvg.t(s),CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s),'Color',colors.r_r)
        plot(CycAvg.t(s),CycAvg.rr_cycavg(s),'Color',colors.r_r,'LineWidth',2);
        %LE-LARP
        plot(CycAvg.t(s),CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s),'Color',colors.l_l)
        plot(CycAvg.t(s),CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s),'Color',colors.l_l)
        plot(CycAvg.t(s),CycAvg.ll_cycavg(s),'Color',colors.l_l,'LineWidth',2);
        %RE-LARP
        plot(CycAvg.t(s),CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s),'Color',colors.r_l)
        plot(CycAvg.t(s),CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s),'Color',colors.r_l)
        plot(CycAvg.t(s),CycAvg.rl_cycavg(s),'Color',colors.r_l,'LineWidth',2);
    end
    hold off  
    axis([0 0.5 YLim])           
    text(0.5,YLim(2),curr_lab{i},'Color',curr_col,'HorizontalAlignment','right','VerticalAlignment','top')
    if(grid_on)
        set(gca,'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
    end    
    if mod(i,n_col)==1
        text(0.5,YLim(1),['n=',num2str(length(CycAvg.cyclist))],'Color','k','HorizontalAlignment','right','VerticalAlignment','bottom')
    else
        text(0.5,YLim(1),num2str(length(CycAvg.cyclist)),'Color','k','HorizontalAlignment','right','VerticalAlignment','bottom')
    end
end
annotation('line',fig_row_pos(end)+[0 x_wid],y_min-0.01*[1 1],'LineWidth',2) 
annotation('line',(x_max+space_x)*[1 1],y_min+[0 yscale*y_wid/(2*YLim(2))],'LineWidth',2) 
annotation('textbox','String','0.5s','EdgeColor','none',...
    'Position',[fig_row_pos(end),0,x_wid,y_min-0.01],'HorizontalAlignment','right','VerticalAlignment','middle')
annotation('textbox','String',[num2str(yscale),'\circ/s'],'EdgeColor','none',...
    'Position',[x_max+space_x,y_min,1-(x_max+space_x),yscale*y_wid/(2*YLim(2))],'HorizontalAlignment','center','VerticalAlignment','middle')
%% Plot 3D Alignment
hg = figure;
colors_H = [colors.l_z;colors.r_z];
if all(contains(fnames,{'LH','LA','LP'})) %Left side
    colors_P = [colors.l_r;colors.r_r];
    colors_A = [colors.l_l;colors.r_l];
    ear = 'L';
else   
    colors_P = [colors.l_l;colors.r_l];
    colors_A = [colors.l_r;colors.r_r];
    ear = 'R';
end
%Posterior
b = load(row1{i1(end)});
a = fieldnames(b);
CycAvg = b.(a{1});
hg = MakeSpherePlot(CycAvg,hg,2,0,0,1,colors_P,ear);
%Anterior
b = load(row3{i3(end)});
a = fieldnames(b);
CycAvg = b.(a{1});
hg = MakeSpherePlot(CycAvg,hg,2,0,0,1,colors_A,ear);
%Horizontal
b = load(row2{i2(end)});
a = fieldnames(b);
CycAvg = b.(a{1});
hg = MakeSpherePlot(CycAvg,hg,2,0,1,1,colors_H,ear);
%% Plot All Electrodes at one Phase Duration
%This can only be done if there is a eeVORResults.mat file
%See if one already exists
temp = dir([path,filesep,'*Results.mat']);
res_file = {temp.name}';
if isempty(res_file)
    disp('First, make the table with the pulse stim parameters...')
    MakeCycleSummaryTable(path,Cyc_Path);
    temp = dir([path,filesep,'*Results.mat']);
    res_file = {temp.name}';
end
load(res_file{end},'all_results')
%Get the stim ear
if contains(all_results{1,1}{:},filesep)
    load(all_results{1,1}{:},'CycAvg')
else
    load([Cyc_Path,filesep,all_results{1,1}{:}],'CycAvg')
end
ear = CycAvg.info.ear;
%Now figure out which files to plot
all_exps = all_results.Condition;
E_inds = false(length(all_exps),9); 
which_files = questdlg('Plot all the files in this directory or manually select them?','','All','Select','Select');
for i = 1:9
    E_sub_i = find(contains(all_exps,['E',num2str(i+2)]));
    if ~isempty(E_sub_i)
        if strcmp(which_files,'Select')
            indx = listdlg('ListString',all_exps(E_sub_i),'PromptString',['Pick the files for E',num2str(i+2)],'ListSize',[400 300]);
            E_inds(E_sub_i(indx),i) = true;
        else
            E_inds(E_sub_i,i) = true;
        end
    else
        disp(['No files found for E',num2str(i+2)])
    end
end
%Make some bold (best if you know which one was activated on)
E_bold = false(1,9);
indx = listdlg('ListString',strcat('E',cellfun(@num2str,num2cell(3:11),'UniformOutput',false)),...
    'PromptString','Pick the electrodes to bold. Press Cancel for none.','ListSize',[400 300],'SelectionMode','multiple');
E_bold(indx) = true;
%Make the figure
ha = gobjects(1,6);
%markerbig=5;
%markersmall=4;
linethick=2;
linethin=1;
errorbarcapsize=1;
figsizeinches=[7,6];
XLim = [-5 105];
YLim_vel = [0,100];
YLim_align = [0 80];

%figsizeinchesBoxplot=[2.3,4];
figure('Units','inch','Position',[2 2 figsizeinches],'Color',[1,1,1]);%CDS083119a
annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,'plotActivationCurrentLevels.m'],'FontSize',5,...
    'EdgeColor','none','interpreter','none','VerticalAlignment','bottom');
ha(1) = subplot(2,3,1);
ha(1).Position = [0.1,0.57,0.29,0.4];
ha(2) = subplot(2,3,2);
ha(2).Position = [0.4,0.57,0.29,0.4];
ha(3) = subplot(2,3,3);
ha(3).Position = [0.7,0.57,0.29,0.4];
ha(4) = subplot(2,3,4);
ha(4).Position = [0.1,0.12,0.29,0.4];
ha(5) = subplot(2,3,5);
ha(5).Position = [0.4,0.12,0.29,0.4];
ha(6) = subplot(2,3,6);
ha(6).Position = [0.7,0.12,0.29,0.4];
markers = {'x','o','d'};
%Set colors (faded vs normal)
if any(contains(all_exps,{'LP';'LH';'LA'})) %Left sided
    color_s = [repmat(colors.l_r_s,3,1);repmat(colors.l_z_s,3,1);repmat(colors.l_l_s,3,1)];
    color = [repmat(colors.l_r,3,1);repmat(colors.l_z,3,1);repmat(colors.l_l,3,1)];
    color(~E_bold,:) = color_s(~E_bold,:);
else %Right Sided
    color_s = [repmat(colors.l_l_s,3,1);repmat(colors.l_z_s,3,1);repmat(colors.l_r_s,3,1)];
    color = [repmat(colors.l_l,3,1);repmat(colors.l_z,3,1);repmat(colors.l_r,3,1)];
    color(~E_bold,:) = color_s(~E_bold,:);
end
%Plot each canal
for i = 1:3
    rel_tab.E1 = all_results(E_inds(:,3*i-2),:);
    rel_tab.E2 = all_results(E_inds(:,3*i-1),:);
    rel_tab.E3 = all_results(E_inds(:,3*i),:);
    rel_bold = E_bold(3*i-2:3*i);
    i_ord = [find(~rel_bold),find(rel_bold)]; %order of plotting so bold in front
    h = gobjects(1,length(markers));
    labs = cell(1,3);  
    %Plot legend
    axes(ha(i))
    hold on
    for j = 1:3
        %Make the fake plots for the labels first
        h(j) = plot(NaN,NaN,'Marker',markers{j},'Color',color(3*i-3+j,:),'LineWidth',1);
        %Make the labels
        if ~isempty(rel_tab.(['E',num2str(j)])) %Something in this one
            exp = strsplit(rel_tab.(['E',num2str(j)]).Condition{1});
            labs{1,j} = strrep([exp{contains(exp,'E')},', ',exp{contains(exp,'us')},'/phase'],'u','\mu');
        else %Remove from the plotting order
            i_ord(i_ord==j) = [];
        end
    end
    h(cellfun(@isempty,labs)) = []; %Take out an electrode with no files
    labs(cellfun(@isempty,labs)) = [];  
    hold off
    %Actually plot when values are present
    for j = 1:length(i_ord)
        exps = rel_tab.(['E',num2str(i_ord(j))]).Condition;
        curr = NaN(1,length(exps));
        for q = 1:length(curr)
            exp = strsplit(exps{q});
            curr(q) = str2double(strrep(exp{contains(exp,'uA')},'uA',''));
        end
        [curr_norm,curr_i] = sort(100*(curr-min(curr))/(max(curr)-min(curr))); 
        %Determine canal
        if any(contains(exp,{'LP','RA'})) %RALP
            canal = 'RALP';
        elseif any(contains(exp,{'LH','RH'})) %LHRH
            canal = 'LHRH';
        else %LARP
            canal = 'LARP';
        end        
        %Extract the relevant vectors
        L_Vel = abs(rel_tab.(['E',num2str(i_ord(j))]).(['L_',canal,'_MaxVel'])(curr_i));
        L_Vel_sd = rel_tab.(['E',num2str(i_ord(j))]).(['L_',canal,'_MaxVel_sd'])(curr_i);
        R_Vel = abs(rel_tab.(['E',num2str(i_ord(j))]).(['R_',canal,'_MaxVel'])(curr_i));
        R_Vel_sd = rel_tab.(['E',num2str(i_ord(j))]).(['R_',canal,'_MaxVel_sd'])(curr_i);        
        L_Align = rel_tab.(['E',num2str(i_ord(j))]).L_Align(curr_i);
        L_Align_sd = rel_tab.(['E',num2str(i_ord(j))]).L_Align_sd(curr_i);
        R_Align = rel_tab.(['E',num2str(i_ord(j))]).R_Align(curr_i);
        R_Align_sd = rel_tab.(['E',num2str(i_ord(j))]).R_Align_sd(curr_i);
        %Choose the larger eye for each current
        [~,eye] = max([L_Vel,R_Vel],[],2);
        Vel = L_Vel;
        Vel(eye==2) = R_Vel(eye==2);
        Vel_sd = L_Vel_sd;
        Vel_sd(eye==2) = R_Vel_sd(eye==2);
        Align = L_Align;
        Align(eye==2) = R_Align(eye==2);
        Align_sd = L_Align_sd;
        Align_sd(eye==2) = R_Align_sd(eye==2);
        %Plot velocity
        axes(ha(i))
        hold on
        errorbar(curr_norm,Vel,Vel_sd,'Color',color(3*i-3+i_ord(j),:),'LineStyle','none','LineWidth',1,'CapSize',errorbarcapsize) 
        plot(curr_norm,Vel,'Marker',markers{i_ord(j)},'Color',color(3*i-3+i_ord(j),:),'LineWidth',1)
        hold off
        %Plot alignment
        axes(ha(i+3))
        hold on
        errorbar(curr_norm,Align,Align_sd,'Color',color(3*i-3+i_ord(j),:),'LineStyle','none','LineWidth',1,'CapSize',errorbarcapsize) 
        plot(curr_norm,Align,'Marker',markers{i_ord(j)},'Color',color(3*i-3+i_ord(j),:),'LineWidth',1)
        hold off        
    end    
    axes(ha(i))
    leg = legend(h,labs,'box','off','Location','northwest','FontSize',7);
    leg.ItemTokenSize(1) = 12;
    set(gca,'box','off')   
    set(gca,'YLim',YLim_vel)
    if i == 1
        set(gca,'YTick',20:20:max(YLim_vel))
        ylabel({'Eye Velocity';'Magnitude (\circ/s)'})
    else
        set(gca,'YColor','none')
    end
    set(gca,'XLim',XLim,'XColor','none')
    axes(ha(i+3))
    set(gca,'YLim',YLim_align)
    if i == 1
        set(gca,'YTick',10:10:max(YLim_align))
        ylabel({'Misalignemnt';'Angle (\circ)'})
    else
        set(gca,'YColor','none')
    end
    set(gca,'XLim',XLim,'XTick',0:20:100)
    if i==2
        xlabel('% of Current Range')
    end    
end