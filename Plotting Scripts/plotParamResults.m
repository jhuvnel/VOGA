%% Plot Param Results.m
% This function makes figures from the Results.mat tables 
% It detects the different possible figures based on the table given.

%type = 'SineAmpVelLRZ';
%type = 'SineAmpVelXYZ';
%type = 'Autoscan';

function plotParamResults(type,path,code_Path,version,Experimenter,annot,YMax)
    if nargin < 7
        YMax = 100;
    end
    switch type
        case 'SineAmpVelLRZ'   
            %% Make Figure like Boutros 2019 Figure 6
            % Initialize
            close all;
            load('VNELcolors.mat','colors')
            code_name = ['Plotting Scripts',filesep,'plotParamResults.m'];
            warning('off')
            sub_info = readtable('SubjectInfo.xlsx');
            warning('on')
            Subs = sub_info{:,1};
            Ears = sub_info{:,2};
            % Load table in question
            res_file = extractfield(dir([path,filesep,'*Results.mat']),'name')';
            if isempty(res_file)
                disp('No table with cycle parameters found on this path.')
                return;
            end
            load(res_file{end},'all_results') 
            %Pick files to run
            [indx,tf] = listdlg('ListString',all_results.File,...
                'PromptString','Pick the files to plot','ListSize',[400 300],...
                'SelectionMode','multiple');
            if tf == 0
                return;
            end
            all_results2 = all_results(indx,:);
            ear = Ears{ismember(Subs,all_results2.Subject{end})}; %only one subject expected
            canals = {'RALP','LHRH','LARP'};
            amps = {'20dps','50dps','100dps','200dps','300dps','400dps'};
            freq = '2Hz';
            exp_ind = zeros(length(amps),length(canals));
            for i = 1:length(amps)
                for j = 1:length(canals)
                    temp = find(contains(all_results2.File,amps{i})&contains(all_results2.File,canals{j})&contains(all_results2.File,['-',freq]));
                    if length(temp)==1
                        exp_ind(i,j) = temp;
                    elseif length(temp)>1 
                        exp_ind(i,j) = NaN;
                    end
                end
            end
            i1 = exp_ind(find(exp_ind,1,'first'));
            exp_parts = split(all_results2.Condition{i1},' ');
            exp_parts(contains(exp_parts,{'dps','LHRH','LARP','RALP'})) = [];
            exp_parts = strjoin(exp_parts,' ');
            exp_name = [all_results2.Subject{i1},' ',all_results2.Visit{i1},' ',datestr(all_results2.Date(i1),'yyyymmdd'),' ',all_results2.Experiment{i1},' ',exp_parts,' ',all_results2.Goggle{i1}];
            %Look at the stimulation half-cycle
            if strcmp(ear,'L')
                L_RZL = [all_results2.MaxVel_LR_LOW,all_results2.MaxVel_LZ_HIGH,all_results2.MaxVel_LL_LOW];
                L_RZL_sd = [all_results2.MaxVel_LR_LOW_sd,all_results2.MaxVel_LZ_HIGH_sd,all_results2.MaxVel_LL_LOW_sd];
                R_RZL = [all_results2.MaxVel_RR_LOW,all_results2.MaxVel_RZ_HIGH,all_results2.MaxVel_RL_LOW];
                R_RZL_sd = [all_results2.MaxVel_RR_LOW_sd,all_results2.MaxVel_RZ_HIGH_sd,all_results2.MaxVel_RL_LOW_sd];
            elseif strcmp(ear,'R')
                L_RZL = [all_results2.MaxVel_LR_HIGH,all_results2.MaxVel_LZ_LOW,all_results2.MaxVel_LL_HIGH];
                L_RZL_sd = [all_results2.MaxVel_LR_HIGH_sd,all_results2.MaxVel_LZ_LOW_sd,all_results2.MaxVel_LL_HIGH_sd];
                R_RZL = [all_results2.MaxVel_RR_HIGH,all_results2.MaxVel_RZ_LOW,all_results2.MaxVel_RL_HIGH];
                R_RZL_sd = [all_results2.MaxVel_RR_HIGH_sd,all_results2.MaxVel_RZ_LOW_sd,all_results2.MaxVel_RL_HIGH_sd];
            else
               error('Subject implant ear was set to niether left nor right.') 
            end
            %Make the figure
            ha = gobjects(1,6);
            linethick=3;
            linethin=1;
            errorbarcapsize=1;
            figsizeinches=[7,6];
            XLim = [0 110];
            YLim = [0,YMax];
            pdom = [5,10,25,50,75,100];
            XTick = [5,25,50,75,100];
            %figsizeinchesBoxplot=[2.3,4];
            figure('Units','inch','Position',[2 2 figsizeinches],'Color',[1,1,1]);%CDS083119a
            if annot
                annotation('textbox',[0 0 1 1],'String',[path,newline,code_Path,filesep,...
                        code_name,newline,...
                        'VOGA',version,newline,Experimenter],'FontSize',5,...
                    'EdgeColor','none','interpreter','none');
            end
            annotation('textbox',[0 .9 1 .1],'String',exp_name,'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
            ha(1) = subplot(2,3,1);
            ha(2) = subplot(2,3,2);
            ha(3) = subplot(2,3,3);
            ha(4) = subplot(2,3,4);
            ha(5) = subplot(2,3,5);
            ha(6) = subplot(2,3,6);
            x_min = 0.12;
            x_max = 0.99;
            x_space = 0.01;
            y_min = 0.09;
            y_max = 0.90;
            y_space = 0.03;
            x_wid = (x_max-x_min-x_space*2)/3;
            y_wid = (y_max-y_min-y_space*1)/2;
            x_pos = x_min:(x_wid+x_space):x_max;
            y_pos = y_min:(y_wid+y_space):y_max;
            ha(1).Position = [x_pos(1),y_pos(2),x_wid,y_wid];
            ha(2).Position = [x_pos(2),y_pos(2),x_wid,y_wid];
            ha(3).Position = [x_pos(3),y_pos(2),x_wid,y_wid];
            ha(4).Position = [x_pos(1),y_pos(1),x_wid,y_wid];
            ha(5).Position = [x_pos(2),y_pos(1),x_wid,y_wid];
            ha(6).Position = [x_pos(3),y_pos(1),x_wid,y_wid];
            set(ha,'YLim',YLim,'XLim',XLim,'FontSize',12)
            set(ha,'XTick',XTick,'YTick',20:20:YLim(2))
            set(ha(1:3),'XTickLabel',[])
            set(ha([2,3,5,6]),'YTickLabel',[],'YColor','none')
            title(ha(1),'RALP','FontSize',16,'Color',colors.l_r)
            title(ha(2),'LHRH','FontSize',16,'Color',colors.l_z)
            title(ha(3),'LARP','FontSize',16,'Color',colors.l_l)
            xlabel(ha(5),'% Modulation Depth','FontSize',14)
            ylabel(ha(1),'Left Eye','FontSize',14)
            ylabel(ha(4),'Right Eye','FontSize',14)
            text(ha(4),-40,30,'VOR Component Magnitude (\circ/s)','Rotation',90,'FontSize',16,'FontWeight','bold')
            colorsL = [colors.l_r;colors.l_z;colors.l_l];
            colorsR = [colors.r_r;colors.r_z;colors.r_l];
            hL = gobjects(1,3);
            hR = gobjects(1,3);
            for i = 1:3 %tested canal
                %Left Eye
                axes(ha(i))
                hold on
                for j = 1:3
                    errorbar(pdom,L_RZL(exp_ind(:,i),j),L_RZL_sd(exp_ind(:,i),j),'Color',colorsL(j,:),'LineStyle','none','LineWidth',1,'CapSize',errorbarcapsize) 
                    plot(pdom,L_RZL(exp_ind(:,i),j),'Color',colorsL(j,:),'LineWidth',linethin)
                end
                hL(i) = plot(pdom,L_RZL(exp_ind(:,i),i),'Color',colorsL(i,:),'LineWidth',linethick);
                hold off                    
                %Right Eye
                axes(ha(i+3))
                hold on
                for j = 1:3
                    errorbar(pdom,R_RZL(exp_ind(:,i),j),R_RZL_sd(exp_ind(:,i),j),'Color',colorsR(j,:),'LineStyle','none','LineWidth',1,'CapSize',errorbarcapsize) 
                    plot(pdom,R_RZL(exp_ind(:,i),j),'Color',colorsR(j,:),'LineWidth',linethin)
                end
                hR(i) = plot(pdom,R_RZL(exp_ind(:,i),i),'Color',colorsR(i,:),'LineWidth',linethick);
                hold off
            end
            leg1 = legend(ha(1),hL([2,3,1]),canals([2,3,1]),'box','off','FontSize',12,'Location','northeast');
            title(leg1,'Left Eye Components','FontWeight','normal')
            leg1.ItemTokenSize(1) = 15;
            leg1.Position = [0.0644   0.7611    0.2470    0.1551];
            leg1.Title.NodeChildren.Position = [0.7500    0.8619         0];
            leg2 = legend(ha(4),hR([2,3,1]),canals([2,3,1]),'box','off','FontSize',12,'Location','northeast');
            title(leg2,'Right Eye Components','FontWeight','normal');
            leg2.ItemTokenSize(1) = 15;
            leg2.Position = [0.0544    0.3410    0.2629    0.1551];
            leg2.Title.NodeChildren.Position = [0.7700    0.8619         0];
            fig_name = inputdlg('Name this figure','',1,{['SineAmp-',strrep(exp_name,' ','-'),'.fig']});
            if ~isempty(fig_name)
                savefig([path,filesep,fig_name{:}])
            end
            close;
        case 'SineAmpVelXYZ'   
            %% Make Figure like Boutros 2019 Figure 6 but X, Y and Z
            % Initialize
            close all;
            load('VNELcolors.mat','colors')
            code_name = ['Plotting Scripts',filesep,'plotParamResults.m'];
            warning('off')
            sub_info = readtable('SubjectInfo.xlsx');
            warning('on')
            Subs = sub_info{:,1};
            Ears = sub_info{:,2};
            % Load table in question
            res_file = extractfield(dir([path,filesep,'*Results.mat']),'name')';
            if isempty(res_file)
                disp('No table with cycle parameters found on this path.')
                return;
            end
            load(res_file{end},'all_results') 
            %Pick files to graph
            [indx,tf] = listdlg('ListString',all_results.File,...
                'PromptString','Pick the files to plot','ListSize',[400 300],...
                'SelectionMode','multiple');
            if tf == 0
                return;
            end
            all_results2 = all_results(indx,:);
            ear = Ears{ismember(Subs,all_results2.Subject{end})}; %only one subject expected
            canals = {'X','Y','LHRH'};
            amps = {'20dps','50dps','100dps','200dps','300dps','400dps'};
            freq = '2Hz';
            exp_ind = zeros(length(amps),length(canals));
            for i = 1:length(amps)
                for j = 1:length(canals)
                    temp = find(contains(all_results2.File,amps{i})&contains(all_results2.File,canals{j})&contains(all_results2.File,['-',freq]));
                    if length(temp)==1
                        exp_ind(i,j) = temp;
                    elseif length(temp)>1 
                        exp_ind(i,j) = NaN;
                    end
                end
            end
            i1 = exp_ind(find(exp_ind,1,'first'));
            exp_parts = split(all_results2.Condition{i1},' ');
            exp_parts(contains(exp_parts,{'dps','X','Y','LHRH'})) = [];
            exp_parts = strjoin(exp_parts,' ');
            exp_name = [all_results2.Subject{i1},' ',all_results2.Visit{i1},' ',datestr(all_results2.Date(i1),'yyyymmdd'),' ',all_results2.Experiment{i1},' ',exp_parts,' ',all_results2.Goggle{i1}];
            %Look at the stimulation half-cycle
            if strcmp(ear,'L')
                L_XYZ = [all_results2.MaxVel_LX_LOW,all_results2.MaxVel_LY_HIGH,all_results2.MaxVel_LZ_LOW];
                L_XYZ_sd = [all_results2.MaxVel_LX_LOW_sd,all_results2.MaxVel_LY_HIGH_sd,all_results2.MaxVel_LZ_LOW_sd];
                R_XYZ = [all_results2.MaxVel_RX_LOW,all_results2.MaxVel_RY_HIGH,all_results2.MaxVel_RZ_LOW];
                R_XYZ_sd = [all_results2.MaxVel_RX_LOW_sd,all_results2.MaxVel_RY_HIGH_sd,all_results2.MaxVel_RZ_LOW_sd];
            elseif strcmp(ear,'R')
                L_XYZ = [all_results2.MaxVel_LX_HIGH,all_results2.MaxVel_LY_LOW,all_results2.MaxVel_LZ_HIGH];
                L_XYZ_sd = [all_results2.MaxVel_LX_HIGH_sd,all_results2.MaxVel_LY_LOW_sd,all_results2.MaxVel_LZ_HIGH_sd];
                R_XYZ = [all_results2.MaxVel_RX_HIGH,all_results2.MaxVel_RY_LOW,all_results2.MaxVel_RZ_HIGH];
                R_XYZ_sd = [all_results2.MaxVel_RX_HIGH_sd,all_results2.MaxVel_RY_LOW_sd,all_results2.MaxVel_RZ_HIGH_sd];
            else
               error('Subject implant ear was set to neither left nor right.') 
            end
            %Make the figure
            ha = gobjects(1,6);
            linethick=3;
            linethin=1;
            errorbarcapsize=1;
            figsizeinches=[7,6];
            XLim = [0 110];
            YLim = [0,YMax];
            pdom = [5,10,25,50,75,100];
            XTick = [5,25,50,75,100];
            %figsizeinchesBoxplot=[2.3,4];
            figure('Units','inch','Position',[2 2 figsizeinches],'Color',[1,1,1]);%CDS083119a
            if annot
                annotation('textbox',[0 0 1 1],'String',[path,newline,code_Path,filesep,...
                        code_name,newline,...
                        'VOGA',version,newline,Experimenter],'FontSize',5,...
                    'EdgeColor','none','interpreter','none');
            end
            annotation('textbox',[0 .9 1 .1],'String',exp_name,'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
            ha(1) = subplot(2,3,1);
            ha(2) = subplot(2,3,2);
            ha(3) = subplot(2,3,3);
            ha(4) = subplot(2,3,4);
            ha(5) = subplot(2,3,5);
            ha(6) = subplot(2,3,6);
            x_min = 0.12;
            x_max = 0.99;
            x_space = 0.01;
            y_min = 0.09;
            y_max = 0.90;
            y_space = 0.03;
            x_wid = (x_max-x_min-x_space*2)/3;
            y_wid = (y_max-y_min-y_space*1)/2;
            x_pos = x_min:(x_wid+x_space):x_max;
            y_pos = y_min:(y_wid+y_space):y_max;
            ha(1).Position = [x_pos(1),y_pos(2),x_wid,y_wid];
            ha(2).Position = [x_pos(2),y_pos(2),x_wid,y_wid];
            ha(3).Position = [x_pos(3),y_pos(2),x_wid,y_wid];
            ha(4).Position = [x_pos(1),y_pos(1),x_wid,y_wid];
            ha(5).Position = [x_pos(2),y_pos(1),x_wid,y_wid];
            ha(6).Position = [x_pos(3),y_pos(1),x_wid,y_wid];
            set(ha,'YLim',YLim,'XLim',XLim,'FontSize',12)
            set(ha,'XTick',XTick,'YTick',20:20:YLim(2))
            set(ha(1:3),'XTickLabel',[])
            set(ha([2,3,5,6]),'YTickLabel',[],'YColor','none')
            title(ha(1),'X','FontSize',16,'Color',colors.l_x)
            title(ha(2),'Y','FontSize',16,'Color',colors.l_y)
            title(ha(3),'Z','FontSize',16,'Color',colors.l_z)
            xlabel(ha(5),'% Modulation Depth','FontSize',14)
            ylabel(ha(1),'Left Eye','FontSize',14)
            ylabel(ha(4),'Right Eye','FontSize',14)
            text(ha(4),-40,YMax/2,'VOR Component Magnitude (\circ/s)','Rotation',90,'FontSize',16,'FontWeight','bold')
            colorsL = [colors.l_x;colors.l_y;colors.l_z];
            colorsR = [colors.r_x;colors.r_y;colors.r_z];
            hL = gobjects(1,3);
            hR = gobjects(1,3);
            for i = 1:3 %tested canal
                %Left Eye
                axes(ha(i))
                hold on
                for j = 1:3
                    errorbar(pdom,L_XYZ(exp_ind(:,i),j),L_XYZ_sd(exp_ind(:,i),j),'Color',colorsL(j,:),'LineStyle','none','LineWidth',1,'CapSize',errorbarcapsize) 
                    plot(pdom,L_XYZ(exp_ind(:,i),j),'Color',colorsL(j,:),'LineWidth',linethin)
                end
                hL(i) = plot(pdom,L_XYZ(exp_ind(:,i),i),'Color',colorsL(i,:),'LineWidth',linethick);
                hold off                    
                %Right Eye
                axes(ha(i+3))
                hold on
                for j = 1:3
                    errorbar(pdom,R_XYZ(exp_ind(:,i),j),R_XYZ_sd(exp_ind(:,i),j),'Color',colorsR(j,:),'LineStyle','none','LineWidth',1,'CapSize',errorbarcapsize) 
                    plot(pdom,R_XYZ(exp_ind(:,i),j),'Color',colorsR(j,:),'LineWidth',linethin)
                end
                hR(i) = plot(pdom,R_XYZ(exp_ind(:,i),i),'Color',colorsR(i,:),'LineWidth',linethick);
                hold off
            end
            leg1 = legend(ha(1),hL,canals,'box','off','FontSize',12,'Location','northeast');
            title(leg1,'Left Eye Components','FontWeight','normal')
            leg1.ItemTokenSize(1) = 15;
            leg1.Position = [0.0644   0.7611    0.2470    0.1551];
            leg1.Title.NodeChildren.Position = [0.7500    0.8619         0];
            leg2 = legend(ha(4),hR,canals,'box','off','FontSize',12,'Location','northeast');
            title(leg2,'Right Eye Components','FontWeight','normal');
            leg2.ItemTokenSize(1) = 15;
            leg2.Position = [0.0544    0.3410    0.2629    0.1551];
            leg2.Title.NodeChildren.Position = [0.7700    0.8619         0];
            fig_name = inputdlg('Name this figure','',1,{['SineAmp-',strrep(exp_name,' ','-'),'.fig']});
            if ~isempty(fig_name)
                savefig([path,filesep,fig_name{:}])
            end
            close;
        case 'Autoscan'
            %% Make Figure like Boutros 2019 Figure 4 but Magnitude and Misalignment
            % Initialize
            close all;
            load('VNELcolors.mat','colors')
            code_name = ['Plotting Scripts',filesep,'plotParamResults.m'];
%             warning('off')
%             sub_info = readtable('SubjectInfo.xlsx');
%             warning('on')
%             Subs = sub_info{:,1};
%             Ears = sub_info{:,2};
            % Load table in question
            res_file = extractfield(dir([path,filesep,'*Results.mat']),'name')';
            if isempty(res_file)
                disp('No table with cycle parameters found on this path.')
                return;
            end
            load(res_file{end},'all_results') 
            %Now figure out which files to plot
            all_exps = all_results.Condition;
            parts = split(all_exps{1},' ');
            exp_name = [all_results.Subject{1},' ',all_results.Visit{1},' ',datestr(all_results.Date(1),'yyyymmdd'),' Autoscan ',parts{contains(parts,'us')},' ',parts{contains(parts,'pps')},' ',all_results.Goggle{1}];
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
            %Find the cycle numbers for each column
            N = [min(all_results.Cycles(any(E_inds(:,1:3),2))),...
                min(all_results.Cycles(any(E_inds(:,4:6),2))),...
                min(all_results.Cycles(any(E_inds(:,7:9),2)))];
            %Make some bold ( if you know which one was activated on)
            E_bold = false(1,9);
            indx = listdlg('ListString',strcat('E',cellfun(@num2str,num2cell(3:11),'UniformOutput',false)),...
                'PromptString','Pick the electrodes to bold. Press Cancel for none.','ListSize',[400 300],'SelectionMode','multiple');
            E_bold(indx) = true;
            %Make the figure
            ha = gobjects(1,6);
            %markerbig=5;
            %markersmall=4;
            %linethick=2;
            %linethin=1;
            errorbarcapsize=1;
            figsizeinches=[7,6];
            XLim = [-5 105];
            YLim_vel = [0,YMax];
            YLim_align = [0 80];
            %figsizeinchesBoxplot=[2.3,4];
            figure('Units','inch','Position',[2 2 figsizeinches],'Color',[1,1,1]);%CDS083119a
            if annot
                annotation('textbox',[0 0 1 1],'String',[path,newline,code_Path,filesep,...
                        code_name,newline,...
                        'VOGA',version,newline,Experimenter],'FontSize',5,...
                    'EdgeColor','none','interpreter','none');
            end
            annotation('textbox',[0 .9 1 .1],'String',strrep(exp_name,'us','\mus'),'FontSize',14,...
                'HorizontalAlignment','center','EdgeColor','none');
            ha(1) = subplot(2,3,1);
            ha(2) = subplot(2,3,2);
            ha(3) = subplot(2,3,3);
            ha(4) = subplot(2,3,4);
            ha(5) = subplot(2,3,5);
            ha(6) = subplot(2,3,6);
            x_min = 0.1;
            x_max = 0.99;
            x_space = 0.01;
            y_min = 0.08;
            y_max = 0.93;
            y_space = 0.03;
            x_wid = (x_max-x_min-2*x_space)/3;
            y_wid = (y_max-y_min-y_space)/2;
            x_pos = x_min:x_wid+x_space:x_max;
            y_pos = y_min:y_wid+y_space:y_max;
            ha(1).Position = [x_pos(1),y_pos(2),x_wid,y_wid];
            ha(2).Position = [x_pos(2),y_pos(2),x_wid,y_wid];
            ha(3).Position = [x_pos(3),y_pos(2),x_wid,y_wid];
            ha(4).Position = [x_pos(1),y_pos(1),x_wid,y_wid];
            ha(5).Position = [x_pos(2),y_pos(1),x_wid,y_wid];
            ha(6).Position = [x_pos(3),y_pos(1),x_wid,y_wid];
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
                        canal = 'R';
                    elseif any(contains(exp,{'LH','RH'})) %LHRH
                        canal = 'Z';
                    else %LARP
                        canal = 'L';
                    end        
                    %Extract the relevant vectors
                    L_Vel = abs(rel_tab.(['E',num2str(i_ord(j))]).(['MaxVel_L',canal,'_HIGH'])(curr_i));
                    L_Vel_sd = rel_tab.(['E',num2str(i_ord(j))]).(['MaxVel_L',canal,'_HIGH_sd'])(curr_i);
                    R_Vel = abs(rel_tab.(['E',num2str(i_ord(j))]).(['MaxVel_R',canal,'_HIGH'])(curr_i));
                    R_Vel_sd = rel_tab.(['E',num2str(i_ord(j))]).(['MaxVel_R',canal,'_HIGH_sd'])(curr_i);        
                    L_Align = rel_tab.(['E',num2str(i_ord(j))]).Align_L_HIGH(curr_i);
                    L_Align_sd = rel_tab.(['E',num2str(i_ord(j))]).Align_L_HIGH_sd(curr_i);
                    R_Align = rel_tab.(['E',num2str(i_ord(j))]).Align_R_HIGH(curr_i);
                    R_Align_sd = rel_tab.(['E',num2str(i_ord(j))]).Align_R_HIGH_sd(curr_i);
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
                    ylabel({'Misalignment';'Angle (\circ)'})
                else
                    set(gca,'YColor','none')
                end
                set(gca,'XLim',XLim,'XTick',0:20:100)
                if i==2
                    xlabel('% of Current Range')
                end    
            end
            fig_name = inputdlg('Name this figure','',1,{[strrep(exp_name,' ','-'),'.fig']});
            if ~isempty(fig_name)
                savefig([path,filesep,fig_name{:}])
            end
            close;
        case 'SpherePlot'
            %% Disco Ball plot
            close all;
            load('VNELcolors.mat','colors')
            code_name = ['Plotting Scripts',filesep,'plotParamResults.m'];
            warning('off')
            sub_info = readtable('SubjectInfo.xlsx');
            warning('on')
            Subs = sub_info{:,1};
            Ears = sub_info{:,2};
            % Load table in question
            res_file = extractfield(dir([path,filesep,'*Results.mat']),'name')';
            if isempty(res_file)
                disp('No table with cycle parameters found on this path.')
                return;
            end
            load(res_file{end},'all_results') 
            %Pick files to graph
            [indx,tf] = listdlg('ListString',all_results.File,...
                'PromptString','Pick the files to plot','ListSize',[400 300],...
                'SelectionMode','multiple');
            if tf == 0
                return;
            end
            all_results2 = all_results(indx,:);
            subject = all_results2.Subject{1};
            hg = figure;
            Function = 2;
            plotstimaxis = 0;
            plotelecaxis = 1;
            normlen = 1;
            stim_ear = Ears{ismember(Subs,subject)};
            for i = 1:size(all_results2,1)
                load([Cyc_Path,filesep,all_results2.File{i}],'CycAvg')
                %Determine canal
                if any(contains(all_results2.Condition{i},{'LP','RA'})) %RALP
                    plot_colors = [colors.l_r;colors.r_r]; 
                elseif any(contains(all_results2.Condition{i},{'LH','RH'})) %LHRH
                    plot_colors = [colors.l_z;colors.r_z]; 
                elseif any(contains(all_results2.Condition{i},{'LA','RP'})) %LARP
                    plot_colors = [colors.l_l;colors.r_l]; 
                else
                    plot_colors = [0,0,0;0.5,0.5,0.5]; %black and gray
                end
                hg = MakeSpherePlot(CycAvg,hg,Function,plotstimaxis,plotelecaxis,normlen,plot_colors,stim_ear);
            end
    end
end