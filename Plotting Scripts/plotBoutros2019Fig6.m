%% Code for making figures likes Boutros 2019 Fig 6
function plotBoutros2019Fig6(params)
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
    res_file = extractfield(dir([Path,filesep,'*Results.mat']),'name')';
    if isempty(res_file)
        disp('No table with cycle parameters found on this path. Creating one now.')
        rerun = ~strcmp(questdlg('If a parameter table already exists, use that one or rerun?','','Use existing table','Rerun','Rerun'),'Use existing table');
        MakeCycleSummaryTable(Path,Cyc_Path,rerun);
        res_file = extractfield(dir([Path,filesep,'*Results.mat']),'name')';
    end
    load(res_file{end},'all_results')
     %% Make Figure like Boutros 2019 Figure 6 
            %Pick files to run
%             [indx,tf] = listdlg('ListString',all_results.File,...
%                 'PromptString','Pick the files to plot','ListSize',[400 300],...
%                 'SelectionMode','multiple');
%             if tf == 0
%                 return;
%             end
%            all_results2 = all_results(indx,:);
            %all_results2 = all_results(contains(all_results.Goggle,'NKI'),:);
            all_results2 = all_results;
            ear = Ears{ismember(Subs,all_results2.Subject{end})}; %only one subject expected
            canals = {'RALP','LHRH','LARP'};
            amps = {'20dps','50dps','100dps','200dps','300dps','500dps'};
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
                annotation('textbox',[0 0 1 1],'String',[Path,newline,code_Path,filesep,...
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
            text(ha(4),-40,15,'VOR Component Magnitude (\circ/s)','Rotation',90,'FontSize',16,'FontWeight','bold')
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
            fig_name ={['SineAmpLRZ-',strrep(exp_name,' ','-'),'.fig']};
            %fig_name = inputdlg('Name this figure','',1,fig_name);
            if ~isempty(fig_name)
                savefig([Path,filesep,fig_name{:}])
            end
            close;
end
    