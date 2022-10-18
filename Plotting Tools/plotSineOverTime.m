load('MVI009R908_aHIT-Manual-RotaryChair-eeVOR_Results.mat')
warning('off')
sub_info = readtable('SubjectInfo.xlsx');
warning('on')
Subs = sub_info{:,1};
Ears = sub_info{:,2};
activ_dates = sub_info{:,3};
%% Display the items in the table
all_results = all_results(contains(all_results.Type,'Sine')&~contains(all_results.Condition,'Light')&contains(all_results.AxisName,{'X','Y','LHRH','LARP','RALP'}),:); %Constrain to sinusoids in the dark
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
%% For Each Subject and Stimulation Axis, plot maximum velocity and phase over time
%Select subject(s)
if length(contents.Subject) == 1
    idx = 1;
else
    [idx,tf] = listdlg('PromptString','Select Subjects to plot','SelectionMode','multiple','ListString',contents.Subject,'InitialValue',1:length(contents.Subject));
    if ~tf
        error('No subjects selected')
    end
end
colors = [colororder;colororder-0.25*(colororder)];
freqs = contents.Frequency_Hz;
amps = contents.Amplitude_dps;
conds = {'MotionMod','ConstantRate','NoStim','eeVOR'};
cond_names = {'Modulated','Constant','None','Electrical'};
cond_mark = {'.','o','x','*'};
max_y = 0.94;
min_y = 0.05;
spac_y = 0.01;
min_x = 0.07;
max_x = 0.95;
spac_x = 0.01;
wid_x = (max_x-min_x);
x_pos = min_x;
wid_y = (max_x-min_x-2*spac_x)/3;
y_pos = min_y:(wid_y+spac_y):max_y-wid_y;
table_labs = all_results.Properties.VariableNames;
for s = 1:length(idx)
    for ax = 1:length(contents.AxisName)
        if strcmp(contents.AxisName{ax},'LHRH')
            let = 'Z';
        else
            let = contents.AxisName{ax}(1);
        end
        %Create just a table of what's needed
        tab = all_results(contains(all_results.Subject,contents.Subject{idx(s)})&...
            contains(all_results.AxisName,contents.AxisName{ax}),:);
        %Only use the eye with the biggest response
        [~,eyeinds] = max([tab.(['MaxVel_L',let,'_HIGH']),tab.(['MaxVel_R',let,'_HIGH']),tab.(['MaxVel_L',let,'_LOW']),tab.(['MaxVel_R',let,'_LOW'])],[],2);
        eye = logical(mod(eyeinds,2));
        tab.MaxVel_HIGH = tab.(['MaxVel_R',let,'_HIGH']);
        tab.MaxVel_HIGH(eye) = tab.(['MaxVel_L',let,'_HIGH'])(eye);
        tab.MaxVel_HIGH_sd = tab.(['MaxVel_R',let,'_HIGH_sd']);
        tab.MaxVel_HIGH_sd(eye) = tab.(['MaxVel_L',let,'_HIGH_sd'])(eye);
        tab.MaxVel_LOW = tab.(['MaxVel_R',let,'_LOW']);
        tab.MaxVel_LOW(eye) = tab.(['MaxVel_L',let,'_LOW'])(eye);
        tab.MaxVel_LOW_sd = tab.(['MaxVel_R',let,'_LOW_sd']);
        tab.MaxVel_LOW_sd(eye) = tab.(['MaxVel_L',let,'_LOW_sd'])(eye);
        tab.Align_HIGH = tab.Align_R_HIGH;
        tab.Align_HIGH(eye) = tab.Align_L_HIGH(eye);
        tab.Align_HIGH_sd = tab.Align_R_HIGH_sd;
        tab.Align_HIGH_sd(eye) = tab.Align_L_HIGH_sd(eye);
        tab.Align_LOW = tab.Align_R_LOW;
        tab.Align_LOW(eye) = tab.Align_L_LOW(eye);
        tab.Align_LOW_sd = tab.Align_R_LOW_sd;
        tab.Align_LOW_sd(eye) = tab.Align_L_LOW_sd(eye);
        tab.Phase = tab.Phase_R;
        tab.Phase(eye) = tab.Phase_L(eye);
        tab.Phase_sd = tab.Phase_R_sd;
        tab.Phase_sd(eye) = tab.Phase_L_sd(eye);   
        if length(unique(tab.('Frequency(Hz)')))>1
            %% Frequency sweep
            figure('Units','inches','Position',[1,0.5,8,4]);
            fig_name = [contents.Subject{idx(s)},' ',strjoin(contents.Experiment,'/'),' ',contents.AxisName{ax}];
            annotation('textbox',[0 .9 1 .1],'String',[fig_name,' Maximum Eye Velocity and Phase'],'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
            ha = gobjects(1,3);
            ha(1) = subplot(3,1,1);
            ha(2) = subplot(3,1,2);
            ha(3) = subplot(3,1,3);
            f_rel = false(1,length(freqs));
            c_rel = false(1,length(conds));
            for c = 1:length(conds)
                for fa = 1:length(freqs)
                    sub_tab = tab(tab.('Frequency(Hz)')==freqs(fa)&...
                        tab.('Amplitude(dps)')==100&...
                        contains(tab.Condition,conds{c}),:);
                    if ~isempty(sub_tab)
                        f_rel(fa) = true;
                        c_rel(c) = true;
                        t = days(sub_tab.Date - activ_dates(contains(Subs,contents.Subject{s})))/365.25;
                        axes(ha(1))%High half-cycle vel
                        hold on
                        errorbar(t,sub_tab.MaxVel_HIGH,sub_tab.MaxVel_HIGH_sd,['--',cond_mark{c}],'Color',colors(fa,:))
                        errorbar(t,sub_tab.MaxVel_LOW,sub_tab.MaxVel_LOW_sd,[':',cond_mark{c}],'Color',colors(fa,:))
                        hold off
                        axes(ha(2))
                        hold on
                        errorbar(t,sub_tab.Phase,sub_tab.Phase_sd,['-',cond_mark{c}],'Color',colors(fa,:))
                        hold off
                        axes(ha(3))
                        hold on
                        errorbar(t,sub_tab.Align_HIGH,sub_tab.Align_HIGH_sd,['--',cond_mark{c}],'Color',colors(fa,:))
                        errorbar(t,sub_tab.Align_LOW,sub_tab.Align_LOW_sd,[':',cond_mark{c}],'Color',colors(fa,:))
                        hold off
                    end
                end
            end   
            %%
            title(ha(1),'Maximum Velocity')
            title(ha(1,2),'Phase Lead')
            title(ha(1,3),'Misalignment')
            linkaxes(ha,'x')
            %Make the legends
            axes(ha(2))
            h1 = gobjects(length(freqs),1);
            h1_lab = strcat(strrep(cellstr(num2str(freqs(f_rel))),' ',''),'Hz');
            hold on
            for i = 1:length(freqs)
            	h1(i) = plot(NaN,NaN,'Color',colors(i,:),'LineWidth',1);
            end
            hold off
            leg1 = legend(h1(f_rel),h1_lab,'NumColumns',3,'Location','northwest');
            leg1.ItemTokenSize(1) = 5;
            title(leg1,'Frequencies')
            
            axes(ha(3))
            h2 = gobjects(6,1);
            h2_labs = [cond_names,{['+',let],['-',let]}];
            hold on
            for i = 1:length(conds)
                h2(i) = plot(NaN,NaN,cond_mark{i},'Color','k');
            end
            h2(5) = plot(NaN,NaN,'--','Color','k');
            h2(6) = plot(NaN,NaN,':','Color','k');
            hold off
            leg2 = legend(h2(c_rel),h2_labs(c_rel),'NumColumns',1,'Location','northwest');
            leg2.ItemTokenSize(1) = 20;
            title(leg2,'Conditions')
            %savefig([out_path,filesep,datestr(now,'yyyymmdd_HHMMSS'),'-',strrep(fig_name,' ','-'),'-MaxEyeVelandPhase',plot_type,'OverTime.fig']);
            %close;
        end
    end
end