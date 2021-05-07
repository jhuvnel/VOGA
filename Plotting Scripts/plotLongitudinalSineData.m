%% Load in files
exp_type = 'eeVOR'; %eeVOR/Rotary Chair
load('VNELColors.mat','colors')
warning('off')
sub_info = readtable('SubjectInfo.xlsx');
warning('on')
Subs = sub_info{:,1};
Ears = sub_info{:,2};
activ_dates = sub_info{:,3};
if ispc
    drive_name = 'Z:';
else
    drive_name = '/Volumes/vnelhuman$/MVI';
end
out_path = [drive_name,filesep,'DATA SUMMARY',filesep,'IN PROGRESS',filesep,'VOR - ',exp_type];
res_files = extractfield(dir([out_path,filesep,'*Results.mat']),'name');
load([out_path,filesep,res_files{end}],'all_results')
disp(['File: ',out_path,filesep,res_files{end}])
% See what items are in the file
all_results = all_results(contains(all_results.Type,'Sine')&~contains(all_results.Condition,'Light'),:); %Constrain to sinusoids in the dark
tab_labs = all_results.Properties.VariableNames;
desc_labs = tab_labs(find(contains(tab_labs,'File'))+1:find(contains(tab_labs,'Cycles'))-1);
desc_labs(contains(desc_labs,'Date')) = [];
labs = strrep(strrep(desc_labs,'(','_'),')','');
for i = 1:length(desc_labs)
    column = all_results.(desc_labs{i});
    if isnumeric(column)
        column(isnan(column)) = [];
        contents.(labs{i}) = unique(column);
    elseif isnumeric(column{1})
        column = cell2mat(column); 
        contents.(labs{i}) = unique(column,'rows');
    else
        column(cellfun(@isempty,column)) = [];
        if all(isempty(column))
            contents.(labs{i}) = [];
        else
            contents.(labs{i}) = unique(column);
        end
    end    
    disp([desc_labs{i},': '])
    disp(contents.(labs{i}))
end
%% Plot Magnitude and Phase over Time for a Frequency Sweep
fanum = length(contents.Frequency_Hz);
amp = 100; %Hardcoded right now
max_y = 0.94;
min_y = 0.05;
spac_y = 0.01;
min_x = 0.07;
max_x = 0.95;  
spac_x = 0.01;
wid_y = (max_y-min_y-spac_y*(fanum-1))/(fanum);
y_pos = max_y-wid_y:-(wid_y+spac_y):min_y;
wid_x = (max_x-min_x-2*spac_x)/3;
x_pos = min_x:(wid_x+spac_x):max_x-wid_x;
%For each select subject
[idx,tf] = listdlg('PromptString','Select Subjects to plot','SelectionMode','multiple','ListString',contents.Subject,'InitialValue',1:length(contents.Subject));
if ~tf
    error('No subjects selected')
end
for s = 1:length(idx)
    %For each frequency + amplitude combo (subplots)
    figure('Units','inches','Position',[1,0.5,7.5,9]);
    ha = gobjects(fanum,3);
    for fa = 1:fanum
        ha(fa,1) = subplot(fanum,3,3*fa-2); 
        ha(fa,2) = subplot(fanum,3,3*fa-1); 
        ha(fa,3) = subplot(fanum,3,3*fa);
    end
    for fa = 1:fanum
        ha(fa,1).Position = [x_pos(1),y_pos(fa),wid_x,wid_y];
        ha(fa,2).Position = [x_pos(2),y_pos(fa),wid_x,wid_y];
        ha(fa,3).Position = [x_pos(3),y_pos(fa),wid_x,wid_y];
    end
    for fa = 1:fanum     
        for c = 1:length(contents.Condition)
            sub_tab = all_results(contains(all_results.Subject,contents.Subject{idx(s)})&all_results.('Amplitude(dps)')==amp&all_results.('Frequency(Hz)')==contents.Frequency_Hz(fa)&contains(all_results.Condition,contents.Condition{c}),:);
            switch contents.Condition{c}
                case 'ConstantRate'
                    line_prop = 'o:';
                case 'NoStim'
                    line_prop = 'x--';
                case 'MotionMod'
                    line_prop = '.-';
                case 'eeVOR'
                    line_prop = 'd-.';
                otherwise
                    line_prop = '*-';
            end            
            if ~isempty(sub_tab)
                if contains(sub_tab.AxisName(1),{'X','Y'}) %XYZ
                    p_traces = {'l_x','r_x','l_y','r_y','l_z','r_z'};
                    traces = upper(strrep(p_traces,'_',''));
                else
                    p_traces = {'l_l','r_l','l_r','r_r','l_z','r_z'};
                    traces = upper(strrep(p_traces,'_',''));
                end
                if contains(sub_tab.AxisName(1),'LHRH')
                    let = 'z';
                else
                    let = lower(sub_tab.AxisName(1));
                end
                t = days(sub_tab.Date - activ_dates(contains(Subs,contents.Subject{s})))/365.25;
                axes(ha(fa,1))%Left side of plots is Pos half-cycle vel
                hold on
                for i = 1:length(traces)
                    errorbar(t,sub_tab.(['MaxVel_',upper(strrep(traces{i},'_','')),'_HIGH']),sub_tab.(['MaxVel_',traces{i},'_HIGH_sd']),line_prop,'Color',colors.(p_traces{i}))
                end
                hold off                
                axes(ha(fa,2)) %Middle side of plots is NEG half-cycle vel
                hold on
                for i = 1:length(traces)
                    errorbar(t,sub_tab.(['MaxVel_',traces{i},'_LOW']),sub_tab.(['MaxVel_',traces{i},'_LOW_sd']),line_prop,'Color',colors.(p_traces{i}))
                end
                hold off
                axes(ha(fa,3))
                hold on
                errorbar(t,sub_tab.Phase_L,sub_tab.Phase_L_sd,line_prop,'Color',colors.(['l_',let]))
                errorbar(t,sub_tab.Phase_R,sub_tab.Phase_R_sd,line_prop,'Color',colors.(['r_',let]))
                hold off
            end
        end
        ylabel(ha(fa,1),[num2str(contents.Frequency_Hz(fa)),'Hz'],'FontWeight','bold')
        set(ha(fa,3),'YAxisLocation','right')
        if fa~=fanum 
            set(ha(fa,1),'XTickLabel',[])
            set(ha(fa,2),'XTickLabel',[],'YTickLabel',[])
            set(ha(fa,3),'XTickLabel',[])
        else
            set(ha(fa,2),'YTickLabel',[])
            xlabel(ha(fa,2),'Time Since Activation (yrs)')
        end
        linkaxes([ha(fa,1),ha(fa,2)],'y')
    end
    title(ha(1,1),'Positive Half-Cycle Velocity')
    title(ha(1,2),'Negative Half-Cycle Velocity')
    title(ha(1,3),'Phase Lead')
    linkaxes(ha,'x')
    linkaxes(ha(:,3),'y')
    %Make the legend
    axes(ha(1,2))
    h1 = gobjects(6,1);
    hold on
    for i = 1:length(traces)
        h1(i) = plot(NaN,NaN,'Color',colors.(p_traces{i}),'LineWidth',2);
    end
    hold off
    leg1 = legend(h1,traces,'NumColumns',3,'Location','northwest');
    leg1.ItemTokenSize(1) = 5;
    title(leg1,'Eye Velocity')    
    axes(ha(1,3))
    h2 = gobjects(5,1);
    conds = {'Modulated','Constant','None','Electrical','Misc'};
    line_props = {'.-','o:','x--','d-.','*-'};
    hold on
    for i = 1:length(line_props)
        h2(i) = plot(NaN,NaN,line_props{i},'Color','k');
    end
    hold off
    leg2 = legend(h2,conds,'NumColumns',1,'Location','northwest');
    leg2.ItemTokenSize(1) = 20;
    title(leg2,'Conditions')    
    fig_name = [contents.Subject{idx(s)},' ',strjoin(contents.Experiment,'/'),' Maximum Eye Velocity and Phase'];
    annotation('textbox',[0 .9 1 .1],'String',fig_name,'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
    savefig([out_path,filesep,datestr(now,'yyyymmdd_HHMMSS'),'-',strrep(fig_name,' ','-'),'-FreqSweepOverTime.fig']);
    close;
end
%% Plot Magnitude and Phase over Time for a Frequency Sweep
fanum = length(contents.Amplitude_dps);
freq = 2; %Hardcoded right now
max_y = 0.94;
min_y = 0.05;
spac_y = 0.01;
min_x = 0.07;
max_x = 0.95;  
spac_x = 0.01;
wid_y = (max_y-min_y-spac_y*(fanum-1))/(fanum);
y_pos = max_y-wid_y:-(wid_y+spac_y):min_y;
wid_x = (max_x-min_x-2*spac_x)/3;
x_pos = min_x:(wid_x+spac_x):max_x-wid_x;
%For each select subject
[idx,tf] = listdlg('PromptString','Select Subjects to plot','SelectionMode','multiple','ListString',contents.Subject,'InitialValue',1:length(contents.Subject));
if ~tf
    error('No subjects selected')
end
for s = 1:length(idx)
    %For each frequency + amplitude combo (subplots)
    figure('Units','inches','Position',[1,0.5,7.5,9]);
    ha = gobjects(fanum,3);
    for fa = 1:fanum
        ha(fa,1) = subplot(fanum,3,3*fa-2); 
        ha(fa,2) = subplot(fanum,3,3*fa-1); 
        ha(fa,3) = subplot(fanum,3,3*fa);
    end
    for fa = 1:fanum
        ha(fa,1).Position = [x_pos(1),y_pos(fa),wid_x,wid_y];
        ha(fa,2).Position = [x_pos(2),y_pos(fa),wid_x,wid_y];
        ha(fa,3).Position = [x_pos(3),y_pos(fa),wid_x,wid_y];
    end
    for fa = 1:fanum     
        for c = 1:length(contents.Condition)
            sub_tab = all_results(contains(all_results.Subject,contents.Subject{idx(s)})&all_results.('Amplitude(dps)')==contents.Amplitude_dps(fa)&all_results.('Frequency(Hz)')==freq&contains(all_results.Condition,contents.Condition{c}),:);
            switch contents.Condition{c}
                case 'ConstantRate'
                    line_prop = 'o:';
                case 'NoStim'
                    line_prop = 'x--';
                case 'MotionMod'
                    line_prop = '.-';
                case 'eeVOR'
                    line_prop = 'd-.';
                otherwise
                    line_prop = '*-';
            end            
            if ~isempty(sub_tab)
                if contains(sub_tab.AxisName(1),{'X','Y'}) %XYZ
                    p_traces = {'l_x','r_x','l_y','r_y','l_z','r_z'};
                    traces = upper(strrep(p_traces,'_',''));
                else
                    p_traces = {'l_l','r_l','l_r','r_r','l_z','r_z'};
                    traces = upper(strrep(p_traces,'_',''));
                end
                if contains(sub_tab.AxisName(1),'LHRH')
                    let = 'z';
                else
                    let = lower(sub_tab.AxisName{1});
                end
                t = days(sub_tab.Date - activ_dates(contains(Subs,contents.Subject{s})))/365.25;
                axes(ha(fa,1))%Left side of plots is Pos half-cycle vel
                hold on
                for i = 1:length(traces)
                    errorbar(t,sub_tab.(['MaxVel_',upper(strrep(traces{i},'_','')),'_HIGH']),sub_tab.(['MaxVel_',traces{i},'_HIGH_sd']),line_prop,'Color',colors.(p_traces{i}))
                end
                hold off                
                axes(ha(fa,2)) %Middle side of plots is NEG half-cycle vel
                hold on
                for i = 1:length(traces)
                    errorbar(t,sub_tab.(['MaxVel_',traces{i},'_LOW']),sub_tab.(['MaxVel_',traces{i},'_LOW_sd']),line_prop,'Color',colors.(p_traces{i}))
                end
                hold off
                axes(ha(fa,3))
                hold on
                errorbar(t,sub_tab.Phase_L,sub_tab.Phase_L_sd,line_prop,'Color',colors.(['l_',let]))
                errorbar(t,sub_tab.Phase_R,sub_tab.Phase_R_sd,line_prop,'Color',colors.(['r_',let]))
                hold off
            end
        end
        ylabel(ha(fa,1),[num2str(contents.Amplitude_dps(fa)),'dps'],'FontWeight','bold')
        set(ha(fa,3),'YAxisLocation','right')
        if fa~=fanum 
            set(ha(fa,1),'XTickLabel',[])
            set(ha(fa,2),'XTickLabel',[],'YTickLabel',[])
            set(ha(fa,3),'XTickLabel',[])
        else
            set(ha(fa,2),'YTickLabel',[])
            xlabel(ha(fa,2),'Time Since Activation (yrs)')
        end
        linkaxes([ha(fa,1),ha(fa,2)],'y')
    end
    title(ha(1,1),'Positive Half-Cycle Velocity')
    title(ha(1,2),'Negative Half-Cycle Velocity')
    title(ha(1,3),'Phase Lead')
    linkaxes(ha,'x')
    linkaxes(ha(:,3),'y')
    %Make the legend
    axes(ha(1,2))
    h1 = gobjects(6,1);
    hold on
    for i = 1:length(traces)
        h1(i) = plot(NaN,NaN,'Color',colors.(p_traces{i}),'LineWidth',2);
    end
    hold off
    leg1 = legend(h1,traces,'NumColumns',3,'Location','northwest');
    leg1.ItemTokenSize(1) = 5;
    title(leg1,'Eye Velocity')    
    axes(ha(1,3))
    h2 = gobjects(5,1);
    conds = {'Modulated','Constant','None','Electrical','Misc'};
    line_props = {'.-','o:','x--','d-.','*-'};
    hold on
    for i = 1:length(line_props)
        h2(i) = plot(NaN,NaN,line_props{i},'Color','k');
    end
    hold off
    leg2 = legend(h2,conds,'NumColumns',1,'Location','northwest');
    leg2.ItemTokenSize(1) = 20;
    title(leg2,'Conditions')    
    fig_name = [contents.Subject{idx(s)},' ',strjoin(contents.Experiment,'/'),' Maximum Eye Velocity and Phase'];
    annotation('textbox',[0 .9 1 .1],'String',fig_name,'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
    savefig([out_path,filesep,datestr(now,'yyyymmdd_HHMMSS'),'-',strrep(fig_name,' ','-'),'-AmpSweepOverTime.fig']);
    close;
end