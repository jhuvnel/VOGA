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
%Constrain to sinusoids in the dark and in a X/Y/Z/L/R planes
all_results = all_results(contains(all_results.Type,'Sine')&~contains(all_results.Condition,'Light')&contains(all_results.AxisName,{'X','Y','LHRH','LARP','RALP'}),:); 
% See what items are in the file
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
%% Plot Gain in each Half-Cycle and Phase over Frequency/Amplitude Sweep for Visit 0, Visit 3, and most recent visit
%Select subject(s)
[idx,tf] = listdlg('PromptString','Select Subjects to plot','SelectionMode','multiple','ListString',contents.Subject,'InitialValue',1:length(contents.Subject));
if ~tf
    error('No subjects selected')
end
for plots = 1:2
    if plots==1
        freqs = contents.Frequency_Hz;
        amps = 100; %Hardcoded for 100dps
        fanum = length(freqs);
        x_val = 'Frequency(Hz)';
        xlab = 'Frequency (Hz)';
        plot_type = 'FreqSweep';
    elseif plots == 2
        amps = contents.Amplitude_dps;
        fanum = length(freqs);
        freqs = 2*ones(1,fanum); %Hardcoded for 2Hz
        x_val = 'Amplitude(dps)';
        xlab = 'Amplitude (dps)';
        plot_type = 'AmpSweep';
    end
    if fanum > 1
        axnum = length(contents.AxisName);
        ncol = 3;
        max_y = 0.94;
        min_y = 0.05;
        spac_y = 0.01;
        min_x = 0.07;
        max_x = 0.95;
        spac_x = 0.01;
        wid_y = (max_y-min_y-spac_y*(axnum-1))/axnum;
        y_pos = max_y-wid_y:-(wid_y+spac_y):min_y;
        wid_x = (max_x-min_x-spac_x*(ncol-1))/ncol;
        x_pos = min_x:(wid_x+spac_x):max_x-wid_x;
        for s = 1:length(idx)
            ear = Ears{contains(Subs,contents.Subject{idx(s)})};
            figure('Units','inches','Position',[1,0.5,7.5,9]);
            ha = gobjects(axnum,ncol);
            for ax = 1:axnum
                for col = 1:ncol
                    ha(ax,col) = subplot(axnum,ncol,ncol*ax+col-ncol);
                end
            end
            for ax=1:axnum
                for col = 1:ncol
                    ha(ax,col).Position = [x_pos(col),y_pos(ax),wid_x,wid_y];
                end
            end
            fig_name = [contents.Subject{idx(s)},' ',strjoin(contents.Experiment,'/')];
            annotation('textbox',[0 .9 1 .1],'String',[fig_name,' Half-Cycle Gains and Phase'],'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
            sub_tab = all_results(contains(all_results.Subject,contents.Subject{idx(s)})&...
                    ismember(all_results.('Amplitude(dps)'),amps)&...
                    ismember(all_results.('Frequency(Hz)'),freqs),:);
            sub_tab = sub_tab(contains(sub_tab.Visit,{'Visit0','Visit3',sub_tab.Visit{end}}),:); %Only Visits 0, 3, and most recent
            conds = unique(strcat(sub_tab.Visit,{' '},datestr(sub_tab.Date,'yyyymmdd'),{' '},sub_tab.Condition));
            conds2 = cell(length(conds),3);
            for i = 1:length(conds)
                conds2(i,:) = split(conds(i),' ');
                conds2{i,2} = datetime(conds2{i,2},'InputFormat','yyyyMMdd');
            end
            markers = 'ods^v+><ph.*x_|'; %all the MATLAB built-in markers
            line_prop = cell(length(conds),1);
            for i = 1:length(conds)
                if contains(conds{i},'ConstantRate')
                    lin = ':';
                elseif contains(conds{i},'NoStim')
                    lin = '--';
                elseif contains(conds{i},'MotionMod')
                    lin = '-';  
                elseif contains(conds{i},'eeVOR')
                    lin = '-.';
                else
                    lin = '-';
                end                
                line_prop(i,1) = {[markers(i),lin]};
            end 
            for ax = 1:axnum
                if contains(contents.AxisName{ax},'LHRH')
                    let = 'Z';
                else
                    let = contents.AxisName{ax}(1);                    
                end
                for c = 1:length(conds)
                    sub_tab2 = sub_tab(contains(sub_tab.AxisName,contents.AxisName{ax})&...
                        contains(sub_tab.Visit,conds2{c,1})&...
                        contains(sub_tab.Condition,conds2{c,3})&...
                        sub_tab.Date==conds2{c,2},:);
                    if ~isempty(sub_tab2)
                        color_l = 'k';
                        color_r = 'k';
                        t = sub_tab2.(x_val);
                        if contains(ear,{'L','l'}) %Left Sided Implant: Stim = Pos, Inhib = Neg
                            axes(ha(ax,1))
                            hold on
                            errorbar(t,sub_tab2.(['Gain_L',let,'_HIGH']),sub_tab2.(['Gain_L',let,'_HIGH_sd']),line_prop{c},'Color',color_l)
                            errorbar(t,sub_tab2.(['Gain_R',let,'_HIGH']),sub_tab2.(['Gain_R',let,'_HIGH_sd']),line_prop{c},'Color',color_r)
                            hold off
                            axes(ha(ax,2))
                            hold on
                            errorbar(t,sub_tab2.(['Gain_L',let,'_LOW']),sub_tab2.(['Gain_L',let,'_LOW_sd']),line_prop{c},'Color',color_l)
                            errorbar(t,sub_tab2.(['Gain_R',let,'_LOW']),sub_tab2.(['Gain_R',let,'_LOW_sd']),line_prop{c},'Color',color_r)
                            hold off
                        else %Right Sided Implant: Stim = Neg, Inhib = Pos
                            axes(ha(ax,1))
                            hold on
                            errorbar(t,sub_tab2.(['Gain_L',let,'_LOW']),sub_tab2.(['Gain_L',let,'_LOW_sd']),line_prop{c},'Color',color_l)
                            errorbar(t,sub_tab2.(['Gain_R',let,'_LOW']),sub_tab2.(['Gain_R',let,'_LOW_sd']),line_prop{c},'Color',color_r)
                            hold off
                            axes(ha(ax,2))
                            hold on
                            errorbar(t,sub_tab2.(['Gain_L',let,'_HIGH']),sub_tab2.(['Gain_L',let,'_HIGH_sd']),line_prop{c},'Color',color_l)
                            errorbar(t,sub_tab2.(['Gain_R',let,'_HIGH']),sub_tab2.(['Gain_R',let,'_HIGH_sd']),line_prop{c},'Color',color_r)
                            hold off
                        end
                        axes(ha(ax,3))
                        hold on
                        errorbar(t,sub_tab2.Phase_L,sub_tab2.Phase_L_sd,line_prop{c},'Color',color_l)
                        errorbar(t,sub_tab2.Phase_R,sub_tab2.Phase_R_sd,line_prop{c},'Color',color_r)
                        hold off
                    end
                end
                ylabel(ha(ax,1),contents.AxisName{ax},'FontWeight','bold')
                set(ha(ax,3),'YAxisLocation','right')
                if ax~=axnum
                    set(ha(ax,1),'XTickLabel',[])
                    set(ha(ax,2),'XTickLabel',[],'YTickLabel',[])
                    set(ha(ax,3),'XTickLabel',[])
                else
                    set(ha(ax,2),'YTickLabel',[])
                    xlabel(ha(ax,2),xlab)
                end
                linkaxes([ha(ax,1),ha(ax,2)],'y')
            end
            title(ha(1,1),'Stimulation Half-Cycle Gain')
            title(ha(1,2),'Inhibition Half-Cycle Gain')
            title(ha(1,3),'Phase Lead')
            linkaxes(ha,'x')
            linkaxes(ha(:,3),'y')
            if plots==1
                set(ha,'xscale','log')
                set(ha,'XTick',freqs,'XMinorTick','off')
                set(ha,'XLim',[0.9*freqs(1),1.1*freqs(end)])
            end
            %Make the legend
            axes(ha(1,1))
            h1 = gobjects(length(conds),1);
            hold on
            for i = 1:length(conds)
                h1(i) = plot(NaN,NaN,line_prop{i},'Color','k');
            end
            hold off
            leg1 = legend(h1,conds,'NumColumns',length(conds),'Location','northwest');
            leg1.ItemTokenSize(1) = 10;
            title(leg1,'Conditions')            
            savefig([out_path,filesep,datestr(now,'yyyymmdd_HHMMSS'),'-',strrep(fig_name,' ','-'),plot_type,'.fig']);
            close;
        end
    end
end