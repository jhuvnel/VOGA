%% Plot Group Cyc Avg.m
%This function makes figures with multiple cycle averages of similar
%experiments across one degree of freedom like:
% Sine w/ different frequency (SineFreq)
% Sine w/ different amplitude (SineAmp)
% Sine combinatations of freq and amp (SineManual)
% Autoscan Current Levels (Autoscan)

function plotGroupCycAvg(params)
    Path = params.Path;
    Cyc_Path = params.Cyc_Path;
    code_Path = params.code_Path;
    version = params.version;
    Experimenter = params.Experimenter;
    if isfield(params,'annot')
        annot = params.annot;
    else
        annot = 1;
    end
    if isfield(params,'YMax')
        YMax = params.YMax;
        YMax(isnan(YMax)) = [];
    else
        YMax = [];
    end
    close all;    
    load('VNELcolors.mat','colors')
    code_name = ['Plotting Scripts',filesep,'plotGroupCycAvg.m'];
    warning('off')
    sub_info = readtable('SubjectInfo.xlsx');
    warning('on')
    Subs = sub_info{:,1};
    Ears = sub_info{:,2};
    % Load table in question
    res_file = extractfield(dir([Path,filesep,'*Results.mat']),'name')';
    if isempty(res_file)
        disp('No table with cycle parameters found on this path. Creating one now.')
        rerun = ~strcmp(questdlg('If a parameter table already exists, use that one or rerun?','','Use existing table','Rerun','Rerun'),'Use existing table');
        MakeCycleSummaryTable(Path,Cyc_Path,rerun);
        res_file = extractfield(dir([Path,filesep,'*Results.mat']),'name')';
    end
    load(res_file{end},'all_results')
    if any(contains(all_results.Type,'Sine'))
        %% Sinusoids
        % Plot Group Cycle Average Results
        cyc_files = all_results.File(contains(all_results.Type,{'Sine','Sinusoid'}));
        freqs = unique(all_results.('Frequency(Hz)'));
        fnum = length(freqs);
        amps = unique(all_results.('Amplitude(dps)'));
        anum = length(amps);
        exps = zeros(fnum,anum);
        for i = 1:length(cyc_files)
            exps(ismember(freqs,all_results.('Frequency(Hz)')(i)),ismember(amps,all_results.('Amplitude(dps)')(i))) = 1;
        end
        file_parts = [all_results.Subject,all_results.Visit,cellstr(datestr(all_results.Date,'yyyymmdd')),all_results.Condition,all_results.AxisName,all_results.Goggle,strcat(strrep(cellstr(num2str(all_results.('Frequency(Hz)'))),' ',''),'Hz'),strcat(strrep(cellstr(num2str(all_results.('Amplitude(dps)'))),' ',''),'dps')];
        conds1 = unique(join(file_parts(:,1:6)));
        if any(sum(exps,2)>1)
            rel_freqs = freqs(sum(exps,2)>1);
            conds_f = repmat(conds1,1,length(rel_freqs));
            for i = 1:length(rel_freqs)
                conds_f(:,i) = strcat(conds1,{[' ',num2str(rel_freqs(i)),'Hz']});
            end
            conds_f = reshape(conds_f,[],1);
        else
            conds_f = {};
        end
        if any(sum(exps,1)>1)
            rel_amps = amps(sum(exps,1)>1);
            conds_a = repmat(conds1,1,length(rel_amps));
            for i = 1:length(rel_amps)
                conds_a(:,i) = strcat(conds1,{[' ',num2str(rel_amps(i)),'dps']});
            end
            conds_a = reshape(conds_a,[],1);
        else
            conds_a = {};
        end            
        conds = [conds_f;conds_a];
        if isempty(conds) %No common freq or amp
            conds = conds1;
        end
        enum = size(conds,1);
        for j = 1:enum
            %Set up labels
            if contains(conds(j),'Hz') %Amp Sweep
                labs = strcat(strrep(cellstr(num2str(amps)),' ',''),'dps');                    
            elseif contains(conds(j),'dps') %Freq Sweep
                labs = strcat(strrep(cellstr(num2str(freqs)),' ',''),'Hz');
            else %Combo
                labs = unique(join(file_parts(:,7:8)));
            end
            %Find files to plot--NEED TO FIX THIS
            rel_files = cell(length(labs),1);
            for i = 1:length(labs)                
                files = cyc_files(sum(ismember(file_parts,[split(conds(j));split(labs(i))]'),2)==length([split(conds(j));split(labs(i))]));
                if ~isempty(files)
                    rel_files(i) = files;
                end
            end
            %Only plot if there are more than two CycAvg for the
            %experiments
            if sum(~cellfun(@isempty,rel_files))>1
                fig_name = conds{j};
                if contains(conds{j},{'X','Y'})
                    leg_text = {'Inv Stim','Left X','Right X',...
                        'Left Y','Right Y','Left Z','Right Z'};
                    canals = {'lx','rx','ly','ry','lz','rz'};
                else
                    leg_text = {'Inv Stim','Left LARP','Right LARP',...
                        'Left RALP','Right RALP','Left Z','Right Z'};
                    canals = {'ll','rl','lr','rr','lz','rz'};
                end
                figure('Units','inches','Position',[0.2778    5.8472   17.2222    3.8333],'Color',[1,1,1])
                %Title
                annotation('textbox',[0 .9 1 .1],'String',fig_name,'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
                if annot
                    annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,...
                        code_name,newline,...
                        'VOGA',version,newline,Experimenter],'FontSize',5,...
                    'EdgeColor','none','interpreter','none');
                end
                ha = gobjects(1,length(labs));
                x_space = 0.02;
                x_min = 0.04;
                x_max = 0.98;
                x_wid = (x_max-x_min-x_space*(length(labs)-1))/length(labs);
                y_height = 0.75;
                x_pos = x_min:(x_wid+x_space):x_max;
                y_pos = 0.12;
                YMax_vals = NaN(1,length(length(labs)));
                for i = 1:length(rel_files)                  
                    ha(i) = subplot(1000,1000,1);
                    ha(i).Position = [x_pos(i) y_pos x_wid y_height];
                    if ~isempty(rel_files{i})
                        load([Cyc_Path,filesep,rel_files{i}],'CycAvg')
                        fields = fieldnames(CycAvg);       
                        if ~ismember('t',fields)
                            CycAvg.t = reshape((0:1/CycAvg.Fs:(length(CycAvg.ll_cycavg)-1)/CycAvg.Fs),[],1);
                        else
                            CycAvg.t = reshape(CycAvg.t,[],1);
                        end
                        if length(CycAvg.t) > 1000
                            s = round(linspace(1,length(CycAvg.t),1000));
                        else
                            s = 1:length(CycAvg.t);
                        end
                        [aa,ab] = size(CycAvg.stim);
                        if aa == length(CycAvg.t) && ab ~=1
                            CycAvg.stim = mean(CycAvg.stim,2)';
                        elseif ab == length(CycAvg.t) && aa ~=1
                            CycAvg.stim = mean(CycAvg.stim,1);
                        end
                        h = gobjects(length(canals)+1,1);
                        h(1) = plot(CycAvg.t(s),-CycAvg.stim(s),'k');
                        hold on
                        %Now add the fills and standard deviations and means
                        max_eyevel = zeros(1,length(canals));
                        for ii = 1:length(canals)
                            max_eyevel(ii) = max(abs([reshape(CycAvg.([canals{ii},'_cycavg'])(s)-CycAvg.([canals{ii},'_cycstd'])(s),1,[]),reshape(CycAvg.([canals{ii},'_cycavg'])(s)+CycAvg.([canals{ii},'_cycstd'])(s),1,[])]));
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.([canals{ii},'_cycavg'])(s)-CycAvg.([canals{ii},'_cycstd'])(s),fliplr((CycAvg.([canals{ii},'_cycavg'])(s)+CycAvg.([canals{ii},'_cycstd'])(s)))],colors.([canals{ii}(1),'_',canals{ii}(2),'_s']))
                            plot(CycAvg.t(s),CycAvg.([canals{ii},'_cycavg'])(s) + CycAvg.([canals{ii},'_cycstd'])(s),'Color',colors.([canals{ii}(1),'_',canals{ii}(2)]))
                            plot(CycAvg.t(s),CycAvg.([canals{ii},'_cycavg'])(s) - CycAvg.([canals{ii},'_cycstd'])(s),'Color',colors.([canals{ii}(1),'_',canals{ii}(2)]))
                            h(ii+1) = plot(CycAvg.t(s),CycAvg.([canals{ii},'_cycavg'])(s),'Color',colors.([canals{ii}(1),'_',canals{ii}(2)]),'LineWidth',2);
                        end                        
                        hold off
                    else 
                        h = gobjects(length(canals)+1,1);
                        h(1) = plot(NaN,NaN,'k');
                        hold on
                        for ii = 1:length(canals)
                             h(ii+1) = plot(NaN,NaN,'Color',colors.([canals{ii}(1),'_',canals{ii}(2)]),'LineWidth',2);
                        end   
                        hold off
                        max_eyevel = 0;
                    end
                    title(labs{i})
                    xlabel('Time (s)')
                    YMax_vals(i) = max(max_eyevel);
                    if i == 1
                        ylabel('Angular Velocity (dps)') 
                        leg_reord = [2,4,6,1,3,5,7];
                        leg = legend(ha(1),h(leg_reord),leg_text(leg_reord),'NumColumns',2);
                        leg.ItemTokenSize(1) = 7;
                        leg_pos = leg.Position; %legend position
                        perc_pos = ((y_height+y_pos)-leg_pos(2))/y_height; %percentage of the axes that are just the legend
                        YMax_vals(i) = max(max_eyevel)*0.5/(0.5-perc_pos);
                    end 
                    set(gca,'XLim',[0,CycAvg.t(end)])
                end 
                if isempty(YMax)
                    YLim = [-1 1]*ceil(1.1*max(YMax_vals)/10)*10; %Round up to the nearest 10; 
                else
                    YLim = [-YMax YMax]; 
                end
                set(ha,'YLim',YLim)               
                savefig([Path,filesep,'CycleAverages_',strrep(fig_name,' ','-'),'.fig'])
                close;
            end
        end
        % Parameterized Graphs
        
        
        
    elseif any(contains(all_results.Condition,'Autoscan'))
        %% Autoscan
            fnames = unique(extractfield([dir([Cyc_Path,filesep,'*CurrentFitting*.mat']);...
            dir([Cyc_Path,filesep,'*Autoscan*.mat'])],'name'));
            %Select Files to Plot
            canal = listdlg('PromptString','Which canal?',...
                'ListString',{'LARP','RALP','LHRH','Choose manually'},...
                'SelectionMode','single');
            if ~isempty(canal) % No constraints on file type   
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
                    case 4 %Choose manually
                        indx = listdlg('ListString',fnames,'PromptString','Pick the files to plot in row 1','ListSize',[400 300]);
                        row1 = fnames(indx);
                        indx = listdlg('ListString',fnames,'PromptString','Pick the files to plot in row 2','ListSize',[400 300]);
                        row2 = fnames(indx);
                        indx = listdlg('ListString',fnames,'PromptString','Pick the files to plot in row 3','ListSize',[400 300]);
                        row3 = fnames(indx);
                end
                if ~(length(row3)==length(row1)&&length(row1)==length(row2))
                    error('Unequal number of files selected for each canal')
                end
                fig_title = inputdlg('Name this figure');
                all_canals = [row1,row2,row3];
                n_col = length(row1);
                %Determine the order
                all_curr = zeros(size(all_canals));
                for i = 1:n_col
                    for j = 1:3
                        fparts = split(all_canals{i,j},'-');
                        all_curr(i,j) = str2double(strrep(strrep(fparts{contains(fparts,'uA')},'uA',''),'.mat',''));      
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
                annotation('textbox',[0 .9 1 .1],'String',fig_title,'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
                annotation('textbox',[0 0 1 1],'String',['VOGA',version,...
                            newline,Experimenter],'FontSize',5,...
                            'EdgeColor','none','interpreter','none');
                annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,...
                                code_name],'FontSize',5,...
                            'EdgeColor','none','interpreter','none','VerticalAlignment','bottom');       
                ha = gobjects(1,length(f_order));
                %Set params
                grid_on = true;
                if isempty(YMax)
                    YMax = 100;
                end
                YLim = YMax*[-1 1];
                x_min = 0.01;
                x_max = 0.95;
                space_x = 0.01;
                y_min = 0.08;
                y_max = 0.92;
                space_y = 0.03;
                %Calculate
                x_wid = (x_max - x_min - space_x*(n_col-1))/n_col;
                fig_row_pos = repmat(x_min:(x_wid+space_x):x_max,1,3);
                y_wid = (y_max - y_min - space_y*2)/3;
                fig_col_pos = reshape(repmat(fliplr(y_min:(y_wid+space_y):y_max),n_col,1),[],1)';
                annotation('line',fig_row_pos(end)+[0 x_wid],y_min-0.01*[1 1],'LineWidth',2) 
                annotation('line',(x_max+space_x)*[1 1],y_min+[0 YMax*y_wid/(2*YLim(2))],'LineWidth',2) 
                annotation('textbox','String','0.5s','EdgeColor','none',...
                    'Position',[fig_row_pos(end),0,x_wid,y_min-0.01],'HorizontalAlignment','right','VerticalAlignment','middle')
                annotation('textbox','String',[num2str(YMax),newline,'\circ/s'],'EdgeColor','none',...
                    'Position',[x_max+space_x,y_min,1-(x_max+space_x),YMax*y_wid/(2*YLim(2))],'HorizontalAlignment','center','VerticalAlignment','middle')
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
                savefig([Path,filesep,strrep(fig_title{:},' ','-'),'.fig'])
            end
    end
end