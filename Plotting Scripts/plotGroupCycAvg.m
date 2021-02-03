%% Plot Group Cyc Avg.m
%This function makes figures with multiple cycle averages of similar
%experiments across one degree of freedom like:
% Sine w/ different frequency (SineFreq)
% Sine w/ different amplitude (SineAmp)
% Pulse rate w/ different frequency (PFM)
% Pulse rate w/ different amplitude (PAM)
% Autoscan Current Levels (Autoscan)

function plotGroupCycAvg(type,path,Cyc_Path,code_Path,version,Experimenter)
    load('VNELcolors.mat','colors')
    code_name = 'plotGroupCycAvg.m';
    switch type
        case 'SineFreq'
            cyc_files = extractfield(dir([Cyc_Path,filesep,'*Sine*.mat']),'name');
            cyc_files(contains(cyc_files,'NotAnalyzeable')) = [];
            file_parts = cell(length(cyc_files),2);
            for i = 1:length(cyc_files)
                fname = strrep(strrep(cyc_files{i},'CycAvg_',''),'.mat','');
                fparts = split(fname,'-');
                file_parts(i,2) = {strrep(fparts{contains(fparts,'Hz')},'p','.')};
                fparts(contains(fparts,'Hz')) = [];
                file_parts(i,1) = {strjoin(fparts,' ')};
            end
            exp_name = unique(file_parts(:,1));
            enum = length(exp_name);
            [~,indf] = sort(cellfun(@str2double,strrep(unique(file_parts(:,2)),'Hz',''))); %make sure it's in numerical order
            freqs = unique(file_parts(:,2));
            freqs = freqs(indf);
            fnum = length(freqs);
            for j = 1:enum
                if contains(exp_name{j},{' X ',' Y '})
                    canal = 'XY';
                else
                    canal = 'LRZ';
                end
                figure('Units','inches','Position',[0.2778    5.8472   17.2222    3.8333],'Color',[1,1,1])
                %Title
                annotation('textbox',[0 .9 1 .1],'String',exp_name{j},'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
                annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,...
                    code_name,newline,...
                    'VOGA',version,newline,Experimenter],'FontSize',5,...
                'EdgeColor','none','interpreter','none');
                ha = gobjects(1,fnum);
                x_space = 0.01;
                x_min = 0.04;
                x_max = 0.98;
                x_wid = (x_max-x_min-x_space*(fnum-1))/fnum;
                y_height = 0.75;
                x_pos = x_min:(x_wid+x_space):x_max;
                y_pos = 0.12;
                for i = 1:fnum
                    freq = str2double(strrep(freqs{i},'Hz',''));
                    ha(i) = subplot(1000,1000,1);
                    ha(i).Position = [x_pos(i) y_pos x_wid y_height];
                    if sum(contains(cyc_files,strrep(exp_name{j},' ','-'))&contains(cyc_files,['-',freqs{i}]))==1
                        load([Cyc_Path,filesep,cyc_files{contains(cyc_files,strrep(exp_name{j},' ','-'))&contains(cyc_files,['-',freqs{i}])}],'CycAvg')
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
                        h(1) = plot(CycAvg.t(s),CycAvg.stim(s),'k');
                        hold on
                        %Now add the fills and standard deviations and means
                        if strcmp(canal,'XY')
                            leg_text = {'Stimulus','Left X','Right X',...
                                'Left Y','Right Y','Left Z','Right Z'};
                            %LE-X
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lx_cycavg(s),fliplr((CycAvg.lx_cycavg(s) + CycAvg.lx_cycstd(s)))],colors.l_x_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lx_cycavg(s),fliplr((CycAvg.lx_cycavg(s) - CycAvg.lx_cycstd(s)))],colors.l_x_s)
                            plot(CycAvg.t(s),CycAvg.lx_cycavg(s) + CycAvg.lx_cycstd(s),'Color',colors.l_x)
                            plot(CycAvg.t(s),CycAvg.lx_cycavg(s) - CycAvg.lx_cycstd(s),'Color',colors.l_x)
                            h(2) = plot(CycAvg.t(s),CycAvg.lx_cycavg(s),'Color',colors.l_x,'LineWidth',2);
                            %RE-X
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rx_cycavg(s),fliplr((CycAvg.rx_cycavg(s) + CycAvg.rx_cycstd(s)))],colors.r_x_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rx_cycavg(s),fliplr((CycAvg.rx_cycavg(s) - CycAvg.rx_cycstd(s)))],colors.r_x_s)
                            plot(CycAvg.t(s),CycAvg.rx_cycavg(s) + CycAvg.rx_cycstd(s),'Color',colors.r_x)
                            plot(CycAvg.t(s),CycAvg.rx_cycavg(s) - CycAvg.rx_cycstd(s),'Color',colors.r_x)
                            h(3) = plot(CycAvg.t(s),CycAvg.rx_cycavg(s),'Color',colors.r_x,'LineWidth',2);
                            %LE-Y
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ly_cycavg(s),fliplr((CycAvg.ly_cycavg(s) + CycAvg.ly_cycstd(s)))],colors.l_y_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ly_cycavg(s),fliplr((CycAvg.ly_cycavg(s) - CycAvg.ly_cycstd(s)))],colors.l_y_s)
                            plot(CycAvg.t(s),CycAvg.ly_cycavg(s) + CycAvg.ly_cycstd(s),'Color',colors.l_y)
                            plot(CycAvg.t(s),CycAvg.ly_cycavg(s) - CycAvg.ly_cycstd(s),'Color',colors.l_y)
                            h(4) = plot(CycAvg.t(s),CycAvg.ly_cycavg(s),'Color',colors.l_y,'LineWidth',2);
                            %RE-Y
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ry_cycavg(s),fliplr((CycAvg.ry_cycavg(s) + CycAvg.ry_cycstd(s)))],colors.r_y_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ry_cycavg(s),fliplr((CycAvg.ry_cycavg(s) - CycAvg.ry_cycstd(s)))],colors.r_y_s)
                            plot(CycAvg.t(s),CycAvg.ry_cycavg(s) + CycAvg.ry_cycstd(s),'Color',colors.r_y)
                            plot(CycAvg.t(s),CycAvg.ry_cycavg(s) - CycAvg.ry_cycstd(s),'Color',colors.r_y)
                            h(5) = plot(CycAvg.t(s),CycAvg.ry_cycavg(s),'Color',colors.r_y,'LineWidth',2);
                            %LE-LHRH
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s)))],colors.l_z_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s)))],colors.l_z_s)
                            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
                            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
                            h(6) = plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
                            %RE-LHRH
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s)))],colors.r_z_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s)))],colors.r_z_s)
                            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
                            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
                            h(7) = plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
                        else
                            leg_text = {'Stimulus','Left LARP','Right LARP',...
                                'Left RALP','Right RALP','Left Z','Right Z'};
                            %LE-LARP
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ll_cycavg(s),fliplr((CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s)))],colors.l_l_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ll_cycavg(s),fliplr((CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s)))],colors.l_l_s)
                            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s),'Color',colors.l_l)
                            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s),'Color',colors.l_l)
                            h(2) = plot(CycAvg.t(s),CycAvg.ll_cycavg(s),'Color',colors.l_l,'LineWidth',2);
                            %RE-LARP
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rl_cycavg(s),fliplr((CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s)))],colors.r_l_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rl_cycavg(s),fliplr((CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s)))],colors.r_l_s)
                            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s),'Color',colors.r_l)
                            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s),'Color',colors.r_l)
                            h(3) = plot(CycAvg.t(s),CycAvg.rl_cycavg(s),'Color',colors.r_l,'LineWidth',2);
                            %LE_RALP
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lr_cycavg(s),fliplr((CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s)))],colors.l_r_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lr_cycavg(s),fliplr((CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s)))],colors.l_r_s)
                            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s),'Color',colors.l_r)
                            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s),'Color',colors.l_r)
                            h(4) = plot(CycAvg.t(s),CycAvg.lr_cycavg(s),'Color',colors.l_r,'LineWidth',2);
                            %RE-RALP
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rr_cycavg(s),fliplr((CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s)))],colors.r_r_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rr_cycavg(s),fliplr((CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s)))],colors.r_r_s)
                            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s),'Color',colors.r_r)
                            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s),'Color',colors.r_r)
                            h(5) = plot(CycAvg.t(s),CycAvg.rr_cycavg(s),'Color',colors.r_r,'LineWidth',2);
                            %LE-LHRH
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s)))],colors.l_z_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s)))],colors.l_z_s)
                            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
                            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
                            h(6) = plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
                            %RE-LHRH
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s)))],colors.r_z_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s)))],colors.r_z_s)
                            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
                            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
                            h(7) = plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
                        end
                        hold off
                        XLim = [0 1/freq];
                        YLim = [min(CycAvg.stim) max(CycAvg.stim)];
                    else 
                        plot(NaN,NaN)
                        XLim = [0 1/freq];
                        YLim = [-100 100];
                    end
                    set(gca,'XLim',XLim)
                    set(gca,'YLim',YLim)
                    title(freqs{i})
                    xlabel('Time (s)')
                    if i == 1
                        ylabel('Angular Velocity (dps)')            
                    else
                        set(gca,'YTickLabel',[])
                    end
                end  
                legend(ha(1),h,leg_text)
                fig_name = inputdlg('Name this figure','',1,{[strrep(exp_name{j},' ','-'),'.fig']});
                savefig([path,filesep,fig_name{:}])
            end
        case 'SineAmp'
            cyc_files = extractfield(dir([Cyc_Path,filesep,'*Sine*.mat']),'name');
            cyc_files(contains(cyc_files,'NotAnalyzeable')) = [];
            file_parts = cell(length(cyc_files),2);
            for i = 1:length(cyc_files)
                fname = strrep(strrep(cyc_files{i},'CycAvg_',''),'.mat','');
                fparts = split(fname,'-');
                file_parts(i,2) = fparts(contains(fparts,'dps'));
                fparts(contains(fparts,'dps')) = [];
                file_parts(i,1) = {strjoin(fparts,' ')};
            end
            exp_name = unique(file_parts(:,1));
            enum = length(exp_name);
            [~,indf] = sort(cellfun(@str2double,strrep(unique(file_parts(:,2)),'dps',''))); %make sure it's in numerical order
            amps = unique(file_parts(:,2));
            amps = amps(indf);
            anum = length(amps);
            for j = 1:enum
                if contains(exp_name{j},{' X ',' Y '})
                    canal = 'XY';
                else
                    canal = 'LRZ';
                end
                figure('Units','inches','Position',[0.2778    5.8472   17.2222    3.8333],'Color',[1,1,1])
                %Title
                annotation('textbox',[0 .9 1 .1],'String',exp_name{j},'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
                annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,...
                    code_name,newline,...
                    'VOGA',version,newline,Experimenter],'FontSize',5,...
                'EdgeColor','none','interpreter','none');
                ha = gobjects(1,anum);                
                x_space = 0.025;
                x_min = 0.03;
                x_max = 0.99;
                x_wid = (x_max-x_min-x_space*(anum-1))/anum;
                y_height = 0.75;
                x_pos = x_min:(x_wid+x_space):x_max;
                y_pos = 0.12;
                for i = 1:anum
                    amp = str2double(strrep(amps{i},'dps',''));
                    ha(i) = subplot(1000,1000,1);
                    ha(i).Position = [x_pos(i) y_pos x_wid y_height];
                    if sum(contains(cyc_files,strrep(exp_name{j},' ','-'))&contains(cyc_files,['-',amps{i}]))==1
                        load([Cyc_Path,filesep,cyc_files{contains(cyc_files,strrep(exp_name{j},' ','-'))&contains(cyc_files,['-',amps{i}])}],'CycAvg')
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
                        h(1) = plot(CycAvg.t(s),CycAvg.stim(s),'k');
                        hold on
                        %Now add the fills and standard deviations and means
                        if strcmp(canal,'XY')
                            leg_text = {'Stimulus','Left X','Right X',...
                                'Left Y','Right Y','Left Z','Right Z'};
                            %LE-X
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lx_cycavg(s),fliplr((CycAvg.lx_cycavg(s) + CycAvg.lx_cycstd(s)))],colors.l_x_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lx_cycavg(s),fliplr((CycAvg.lx_cycavg(s) - CycAvg.lx_cycstd(s)))],colors.l_x_s)
                            plot(CycAvg.t(s),CycAvg.lx_cycavg(s) + CycAvg.lx_cycstd(s),'Color',colors.l_x)
                            plot(CycAvg.t(s),CycAvg.lx_cycavg(s) - CycAvg.lx_cycstd(s),'Color',colors.l_x)
                            h(2) = plot(CycAvg.t(s),CycAvg.lx_cycavg(s),'Color',colors.l_x,'LineWidth',2);
                            %RE-X
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rx_cycavg(s),fliplr((CycAvg.rx_cycavg(s) + CycAvg.rx_cycstd(s)))],colors.r_x_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rx_cycavg(s),fliplr((CycAvg.rx_cycavg(s) - CycAvg.rx_cycstd(s)))],colors.r_x_s)
                            plot(CycAvg.t(s),CycAvg.rx_cycavg(s) + CycAvg.rx_cycstd(s),'Color',colors.r_x)
                            plot(CycAvg.t(s),CycAvg.rx_cycavg(s) - CycAvg.rx_cycstd(s),'Color',colors.r_x)
                            h(3) = plot(CycAvg.t(s),CycAvg.rx_cycavg(s),'Color',colors.r_x,'LineWidth',2);
                            %LE-Y
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ly_cycavg(s),fliplr((CycAvg.ly_cycavg(s) + CycAvg.ly_cycstd(s)))],colors.l_y_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ly_cycavg(s),fliplr((CycAvg.ly_cycavg(s) - CycAvg.ly_cycstd(s)))],colors.l_y_s)
                            plot(CycAvg.t(s),CycAvg.ly_cycavg(s) + CycAvg.ly_cycstd(s),'Color',colors.l_y)
                            plot(CycAvg.t(s),CycAvg.ly_cycavg(s) - CycAvg.ly_cycstd(s),'Color',colors.l_y)
                            h(4) = plot(CycAvg.t(s),CycAvg.ly_cycavg(s),'Color',colors.l_y,'LineWidth',2);
                            %RE-Y
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ry_cycavg(s),fliplr((CycAvg.ry_cycavg(s) + CycAvg.ry_cycstd(s)))],colors.r_y_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ry_cycavg(s),fliplr((CycAvg.ry_cycavg(s) - CycAvg.ry_cycstd(s)))],colors.r_y_s)
                            plot(CycAvg.t(s),CycAvg.ry_cycavg(s) + CycAvg.ry_cycstd(s),'Color',colors.r_y)
                            plot(CycAvg.t(s),CycAvg.ry_cycavg(s) - CycAvg.ry_cycstd(s),'Color',colors.r_y)
                            h(5) = plot(CycAvg.t(s),CycAvg.ry_cycavg(s),'Color',colors.r_y,'LineWidth',2);
                            %LE-LHRH
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s)))],colors.l_z_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s)))],colors.l_z_s)
                            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
                            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
                            h(6) = plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
                            %RE-LHRH
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s)))],colors.r_z_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s)))],colors.r_z_s)
                            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
                            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
                            h(7) = plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
                        else
                            leg_text = {'Stimulus','Left LARP','Right LARP',...
                                'Left RALP','Right RALP','Left Z','Right Z'};
                            %LE-LARP
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ll_cycavg(s),fliplr((CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s)))],colors.l_l_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ll_cycavg(s),fliplr((CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s)))],colors.l_l_s)
                            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s),'Color',colors.l_l)
                            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s),'Color',colors.l_l)
                            h(2) = plot(CycAvg.t(s),CycAvg.ll_cycavg(s),'Color',colors.l_l,'LineWidth',2);
                            %RE-LARP
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rl_cycavg(s),fliplr((CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s)))],colors.r_l_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rl_cycavg(s),fliplr((CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s)))],colors.r_l_s)
                            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s),'Color',colors.r_l)
                            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s),'Color',colors.r_l)
                            h(3) = plot(CycAvg.t(s),CycAvg.rl_cycavg(s),'Color',colors.r_l,'LineWidth',2);
                            %LE_RALP
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lr_cycavg(s),fliplr((CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s)))],colors.l_r_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lr_cycavg(s),fliplr((CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s)))],colors.l_r_s)
                            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s),'Color',colors.l_r)
                            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s),'Color',colors.l_r)
                            h(4) = plot(CycAvg.t(s),CycAvg.lr_cycavg(s),'Color',colors.l_r,'LineWidth',2);
                            %RE-RALP
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rr_cycavg(s),fliplr((CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s)))],colors.r_r_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rr_cycavg(s),fliplr((CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s)))],colors.r_r_s)
                            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s),'Color',colors.r_r)
                            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s),'Color',colors.r_r)
                            h(5) = plot(CycAvg.t(s),CycAvg.rr_cycavg(s),'Color',colors.r_r,'LineWidth',2);
                            %LE-LHRH
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s)))],colors.l_z_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s)))],colors.l_z_s)
                            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
                            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
                            h(6) = plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
                            %RE-LHRH
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s)))],colors.r_z_s)
                            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s)))],colors.r_z_s)
                            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
                            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
                            h(7) = plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
                        end
                        hold off
                        XLim = [0 CycAvg.t(end)];
                        YLim = [-amp amp];
                    else 
                        plot(NaN,NaN)
                        XLim = [0 0.5];%2Hz
                        YLim = [-amp amp];
                    end
                    set(gca,'XLim',XLim)
                    set(gca,'YLim',YLim,'YTick',[-amp -amp/2 0 amp/2 amp])
                    title(amps{i})
                    xlabel('Time (s)')
                    if i == 1
                        ylabel('Angular Velocity (dps)')            
                    end
                end  
                legend(ha(1),h,leg_text)
                fig_name = inputdlg('Name this figure','',1,{[strrep(exp_name{j},' ','-'),'.fig']});
                savefig([path,filesep,fig_name{:}])
            end
        case 'PFM'
            disp('No code written for PFM yet')
        case 'PAM'
            disp('No code written for PAM yet')
        case 'Autoscan'
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
                yscale = 75;
                YLim = yscale*[-1 1];
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
                annotation('line',(x_max+space_x)*[1 1],y_min+[0 yscale*y_wid/(2*YLim(2))],'LineWidth',2) 
                annotation('textbox','String','0.5s','EdgeColor','none',...
                    'Position',[fig_row_pos(end),0,x_wid,y_min-0.01],'HorizontalAlignment','right','VerticalAlignment','middle')
                annotation('textbox','String',[num2str(yscale),newline,'\circ/s'],'EdgeColor','none',...
                    'Position',[x_max+space_x,y_min,1-(x_max+space_x),yscale*y_wid/(2*YLim(2))],'HorizontalAlignment','center','VerticalAlignment','middle')
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
                savefig([path,filesep,strrep(fig_title{:},' ','-'),'.fig'])
            end
    end
end