%% Plot Param Results.m
% This function makes figures from the Results.mat tables 
% It detects the different possible figures based on the table given.

function plotParamResults(type,path,Cyc_Path,code_Path,version,Experimenter)
    load('VNELcolors.mat','colors')
    code_name = 'plotParamResults.m';
    Ears = {'L','L','L','L','R','R','L','R','L'};
    switch type
        case 'Autoscan'
            res_file = extractfield(dir([path,filesep,'*Results.mat']),'name')';
            if isempty(res_file)
                disp('No table with the cycle parameters found on this path.')
                return;
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
            %Make some bold ( if you know which one was activated on)
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
    end
end