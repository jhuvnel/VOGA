function plotSummaryFigures(params)
opts = {'Rotary Chair/Sine','vHIT/GNO',...
    'Candidate: Rotary Chair/Sine','Candidate: vHIT/GNO'};
[ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
    'SelectionMode','multiple',...
    'ListSize',[200 125],...
    'ListString',opts);
if tf
    sub_mark = 'xdo^ps+hv<>|_';
    sub_info = params.sub_info;
    all_subs = sub_info.Subject;
    sub_ears = sub_info.Ear;
    results_mat = dir('*VOGResults.mat');
    if isempty(results_mat)
        error('No Results files found in this directory')
    end
    load(results_mat(end).name,'all_results')
    %Make sure it's sorted by Subject first and then by date
    all_results = sortrows(sortrows(all_results,'Date','ascend'),'Subject','ascend');
    all_results_static = all_results; %Make a static copy for the loop
    %Select one candidate
    if any(contains(opts{ind},'Candidate'))
        all_candidate_fold = extractfield(dir('R*'),'name',extractfield(dir('R*'),'isdir'));
        if isempty(all_candidate_fold)
           disp('Found no folders in the format R### as expected.')
           return;
        end
        [ind2,tf2] = nmlistdlg('PromptString','Select a candidate to add to the plot:',...
        'SelectionMode','single','ListSize',[150 125],'ListString',all_candidate_fold);
        if ~tf2
            return;
        end
        candidate = all_candidate_fold{ind2};
    end
    for ii = 1:length(ind)
        %% Setup
        if contains(opts{ind(ii)},'Rotary')
            all_results = all_results_static(contains(all_results_static.Experiment,'RotaryChair')&contains(all_results_static.Type,'Sine'),:);
            if isempty(all_results)
                error('No Rotary Chair Results files found in this directory')
            end
            % Extract relevant normative values
            load('RotaryChairNormativeData.mat','norm_dat')
            %all_freqs = unique(all_results.('Frequency(Hz)'));
            % Limit summary stats to 0.05-1Hz (0.05, 0.1, 0.2, 0.5, 1)
            freqs = [0.05,0.1,0.2,0.5,1];
            remove_phase_with_gain_below = 0.025;
            %Make structs for each condition
            conds = {'pre','post','mmo','cro'};
            cond_logic = [contains(all_results.Visit,'Visit0')&~contains(all_results.Condition,'Light'),...
                contains(all_results.Visit,'Visit3')&contains(all_results.Condition,'NoStim'),...
                contains(all_results.Condition,'Motion'),...
                contains(all_results.Condition,'Constant')];
            %Make the arrays
            %Summary stats
            for c = 1:length(conds)
                %Initialize arrays
                gain.(conds{c}) = NaN(length(all_subs),length(freqs));
                gain_sd.(conds{c}) = NaN(length(all_subs),length(freqs));
                phase.(conds{c}) = NaN(length(all_subs),length(freqs));
                phase_sd.(conds{c}) = NaN(length(all_subs),length(freqs));  
                %Populate them with the available numbers
                for i = 1:length(all_subs)
                    for j = 1:length(freqs)
                        subtab = all_results(find(contains(all_results.Subject,all_subs(i))...
                            &all_results.Frequency==freqs(j)&cond_logic(:,c)&...
                            contains(all_results.AxisName,sub_ears(i)),1,'last'),:);
                        if ~isempty(subtab)
                            gain.(conds{c})(i,j) = subtab.Gain;
                            gain_sd.(conds{c})(i,j) = subtab.Gain_sd;
                            phase.(conds{c})(i,j) = subtab.Phase;
                            phase_sd.(conds{c})(i,j) = subtab.Phase_sd;
                        end
                    end
                end
                %Make summary lines
                phase.(conds{c})(gain.(conds{c})<remove_phase_with_gain_below) = NaN;
                gain_mean.(conds{c}) = mean(gain.(conds{c}),1,'omitnan');
                gain_std.(conds{c}) = std(gain.(conds{c}),[],1,'omitnan');
                phase_mean.(conds{c}) = mean(phase.(conds{c}),1,'omitnan');
                phase_std.(conds{c}) = std(phase.(conds{c}),[],1,'omitnan');
            end
        elseif contains(opts{ind(ii)},'vHIT')
            
        end
        %% Main menu
        switch opts{ind(ii)}
            case 'Rotary Chair/Sine'
                %% Over Freq
                spread = 0.05;
                mark_size = 7;
                h1 = gobjects(6,1);
                h2 = gobjects(length(all_subs),1);
                fig1 = figure;
                set(fig1,'Units','inches','Position',[1 1 6 4],'Color',[1 1 1])
                ha(1) = subplot(2,1,1);
                plot(NaN,NaN)
                hold on
                h1(6) = fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.gain-norm_dat.gain_std,fliplr(norm_dat.gain+norm_dat.gain_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
                h1(5) = plot(norm_dat.freq,norm_dat.gain,'k--','LineWidth',2);
                h1(1) = errorbar((1-1.5*spread)*freqs,gain_mean.pre,gain_std.pre,'k-','LineWidth',1.5);
                h1(2) = errorbar((1-0.5*spread)*freqs,gain_mean.post,gain_std.post,'-.','Color',0.5*[1,1,1],'LineWidth',1);
                h1(3) = errorbar((1+0.5*spread)*freqs,gain_mean.cro,gain_std.cro,'r:','LineWidth',1);
                h1(4) = errorbar((1+1.5*spread)*freqs,gain_mean.mmo,gain_std.mmo,'r-','LineWidth',1.5);
                for i = 1:length(all_subs)
                    plot((1-0.5*spread)*freqs,gain.post(i,:),'.','Color',0.5*[1,1,1],'Marker',sub_mark(i),'MarkerSize',mark_size)
                    plot((1-1.5*spread)*freqs,gain.pre(i,:),'k.','Marker',sub_mark(i),'MarkerSize',mark_size,'MarkerFaceColor','k')
                    plot((1+1.5*spread)*freqs,gain.mmo(i,:),'r.','Marker',sub_mark(i),'MarkerSize',mark_size,'MarkerFaceColor','r')
                    plot((1+0.5*spread)*freqs,gain.cro(i,:),'r.','Marker',sub_mark(i),'MarkerSize',mark_size)
                end
                hold off
                set(gca,'xscale','log','xminortick','off','XTick',freqs,'XTickLabel',[],'box','on','Layer','top')
                axis([0.041 1.2 -0.025 0.79])
                ylabel('Horizontal VOR Gain')
                leg1 = legend(h1,{'Pre-op','Post-op','MVI OFF','MVI ON','Norm Mean','Norm±SD'},'Location','north','NumColumns',3);
                title(leg1,'Condition')
                leg1.ItemTokenSize(1) = 15;
                ha(2) = subplot(2,1,2);
                plot(NaN,NaN)
                hold on
                fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.phase-norm_dat.phase_std,fliplr(norm_dat.phase+norm_dat.phase_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
                plot(norm_dat.freq,norm_dat.phase,'k--','LineWidth',2);
                errorbar((1-0.5*spread)*freqs,phase_mean.post,phase_std.post,'-.','Color',0.5*[1,1,1],'LineWidth',1);
                errorbar((1-1.5*spread)*freqs,phase_mean.pre,phase_std.pre,'k-','LineWidth',1.5);
                errorbar((1+1.5*spread)*freqs,phase_mean.mmo,phase_std.mmo,'r-','LineWidth',1.5);
                errorbar((1+0.5*spread)*freqs,phase_mean.cro,phase_std.cro,'r:','LineWidth',1);
                for i = 1:length(all_subs)
                    plot((1-0.5*spread)*freqs,phase.post(i,:),'.','Color',0.5*[1,1,1],'Marker',sub_mark(i),'MarkerSize',mark_size)
                    plot((1-1.5*spread)*freqs,phase.pre(i,:),'k.','Marker',sub_mark(i),'MarkerSize',mark_size,'MarkerFaceColor','k')
                    plot((1+1.5*spread)*freqs,phase.mmo(i,:),'r.','Marker',sub_mark(i),'MarkerSize',mark_size,'MarkerFaceColor','r')
                    plot((1+0.5*spread)*freqs,phase.cro(i,:),'r.','Marker',sub_mark(i),'MarkerSize',mark_size)
                    h2(i) = plot(NaN,NaN,'k.','Marker',sub_mark(i),'MarkerSize',mark_size,'LineWidth',1);
                end
                hold off
                set(gca,'xscale','log','xminortick','off','XTick',freqs,'box','on','Layer','top')
                axis([0.041 1.2 -25 135])
                xlabel('Frequency (Hz)')
                ylabel('Phase Lead (deg)')
                leg2 = legend(h2,split(cellstr(num2str(1:length(all_subs)))),'Location','north','NumColumns',ceil(length(all_subs)/1));
                leg2.ItemTokenSize(1) = 8;
                title(leg2,'Subjects')
                ha(1).Position = [0.09,0.56,0.98-0.09,0.425];
                ha(2).Position = [0.09,0.12,0.98-0.09,0.425];
                %Figure letter labels
                annot_wid_x = 0.04;
                annot_wid_y = 0.06;
                annot_pos_x = (0.09+0.009)*ones(1,2);
                annot_pos_y = [0.985, 0.545]-0.075;
                annot_string = {'A','B'};
                for i = 1:length(annot_string)
                    annotation('textbox',[annot_pos_x(i),annot_pos_y(i),annot_wid_x,annot_wid_y],...
                        'String',annot_string{i},'HorizontalAlignment','center',...
                        'VerticalAlignment','middle','FontWeight','bold','FontSize',20,...
                        'Fitboxtotext','off','BackgroundColor',[1,1,1]);
                end 
                fname1 = ['Summary Figures',filesep,datestr(now,'yyyymmdd'),'_SummaryRotaryChairGainPhase_AllSub.fig'];
                savefig(fig1,fname1)
                saveas(fig1,strrep(fname1,'.fig','.png'))
                %% Over Freq Grouped Sub (Change groupings of subjects as needed)
                each_unique_fig_text = {'MVI1-5','MVI6-10','MVI11-13'};
                each_sub_num = [{1:5},{6:10},{11:13}];
                for f = 1:length(each_unique_fig_text)
                    unique_fig_text = each_unique_fig_text{f};
                    sub_num = each_sub_num{f}; %Subjects to graph (array of #)
                    %Make fig
                    spread = 0.05;
                    fig2 = figure;
                    set(fig2,'Units','inches','Position',[1 1 5 8],'Color',[1,1,1]);
                    xmin = 0.09;
                    xmax = 0.98;
                    xspc = 0.05;
                    ymin = 0.05;
                    ymax = 0.97;
                    yspc = 0.01;
                    %Initialize axes
                    N = length(sub_num);
                    ha = gobjects(N,2);
                    h1 = gobjects(6,1);
                    %Set positions
                    xwid = (xmax-xmin-xspc)/2; %set for 2 cols
                    x = xmin:(xwid+xspc):xmax;
                    ywid = (ymax-ymin-(N-1)*yspc)/N;
                    y = ymin:(ywid+yspc):ymax;
                    %Plot
                    for i = 1:N
                        ha(i,1) = subplot(N,2,2*i-1);
                        plot(NaN,NaN)
                        hold on
                        h1(6) = fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.gain-norm_dat.gain_std,fliplr(norm_dat.gain+norm_dat.gain_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
                        h1(5) = plot(norm_dat.freq,norm_dat.gain,'k--','LineWidth',2);
                        h1(1) = errorbar((1-1.5*spread)*freqs,gain.pre(sub_num(i),:),gain_sd.pre(sub_num(i),:),'k-','LineWidth',1.5);
                        h1(2) = errorbar((1-0.5*spread)*freqs,gain.post(sub_num(i),:),gain_sd.post(sub_num(i),:),'-.','Color',0.5*[1,1,1],'LineWidth',1.5);
                        h1(4) = errorbar((1+1.5*spread)*freqs,gain.mmo(sub_num(i),:),gain_sd.mmo(sub_num(i),:),'r-','LineWidth',1.5);
                        h1(3) = errorbar((1+0.5*spread)*freqs,gain.cro(sub_num(i),:),gain_sd.cro(sub_num(i),:),'r:','LineWidth',1);
                        if sub_num(i)==1 %Add MVI001 pre-op done elsewhere as black X
                            freqs2 = [0.01;0.02;0.04;0.08;0.16;0.32;0.64];
                            gains2 = [0.0083;0.0250;0.0292;0.0458;0.0500;0.0333;0.0625];
                            plot(freqs2,gains2,'kx','LineWidth',2)
                        end
                        hold off
                        axis([0.041 1.2 -0.025 0.79])
                        ylabel(all_subs{sub_num(i)}(1:6),'FontSize',12,'FontWeight','bold')
                        set(gca,'XTick',freqs,'xscale','log','box','on','Layer','top')
                        if i<N
                            set(gca,'XTickLabel',[])
                        else
                            xlabel('Frequency (Hz)')
                        end
                        if i == 1
                            title('Horizontal VOR Gain')
                        end    
                        ha(i,2) = subplot(N,2,2*i);
                        plot(NaN,NaN)
                        hold on
                        fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.phase-norm_dat.phase_std,fliplr(norm_dat.phase+norm_dat.phase_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
                        plot(norm_dat.freq,norm_dat.phase,'k--','LineWidth',2);
                        errorbar((1-1.5*spread)*freqs,phase.pre(sub_num(i),:),phase_sd.pre(sub_num(i),:),'k-','LineWidth',1.5)
                        errorbar((1-0.5*spread)*freqs,phase.post(sub_num(i),:),phase_sd.post(sub_num(i),:),'-.','Color',0.5*[1,1,1],'LineWidth',1.5)
                        errorbar((1+1.5*spread)*freqs,phase.mmo(sub_num(i),:),phase_sd.mmo(sub_num(i),:),'r-','LineWidth',1.5)
                        errorbar((1+0.5*spread)*freqs,phase.cro(sub_num(i),:),phase_sd.cro(sub_num(i),:),'r:','LineWidth',1)
                        hold off
                        axis([0.041 1.2 -75 135])
                        set(gca,'XTick',freqs,'xscale','log','box','on','Layer','top')
                        if i<N
                            set(gca,'XTickLabel',[])
                        else
                            xlabel('Frequency (Hz)')
                        end
                        if i == 1
                            title('Phase Lead (deg)')
                        end
                    end
                    %Sizing figure and figure letter labels
                    annot_wid_x = 0.047;
                    annot_wid_y = 0.035;
                    annot_pos_x = x+0.009;
                    annot_pos_y = fliplr(y)+ywid-0.04;                     
                    for i = 1:N
                        ha(i,1).Position = [x(1),y(N+1-i),xwid,ywid];
                        ha(i,2).Position = [x(2),y(N+1-i),xwid,ywid];
                        annotation('textbox',[annot_pos_x(1),annot_pos_y(i),annot_wid_x,annot_wid_y],...
                            'String',char(i+64),'HorizontalAlignment','center',...
                            'VerticalAlignment','middle','FontWeight','bold','FontSize',20,...
                            'Fitboxtotext','off','BackgroundColor',[1,1,1]);
                        annotation('textbox',[annot_pos_x(2),annot_pos_y(i),annot_wid_x,annot_wid_y],...
                            'String',char((N+i)+64),'HorizontalAlignment','center',...
                            'VerticalAlignment','middle','FontWeight','bold','FontSize',20,...
                            'Fitboxtotext','off','BackgroundColor',[1,1,1]);
                    end
                    %Add legend
                    leg = legend(ha(1,2),h1,{'Pre-op','Post-op','MVI OFF','MVI ON','Norm Mean','Norm±SD'},'NumColumns',2,'Location','northeast');
                    leg.ItemTokenSize(1) = 15;
                    title(leg,'Condition')                    
                    fname2 = ['Summary Figures',filesep,datestr(now,'yyyymmdd'),'_SummaryRotaryChairGainPhase_EachSub_',unique_fig_text,'.fig'];
                    savefig(fig2,fname2)
                    saveas(fig2,strrep(fname2,'fig','png'))
                end
                %% Each Sub, All Parameters Over Time
                
                
            case 'Candidate: Rotary Chair/Sine'
                results_mat = dir([candidate,filesep,'Visit 0',filesep,'Rotary Chair',filesep,'*Results.mat']);
                if isempty(results_mat)
                    error('No Results files found in the candidate Rotary Chair directory')
                end
                load([candidate,filesep,'Visit 0',filesep,'Rotary Chair',filesep,results_mat(end).name],'all_results')
                all_results = all_results(contains(all_results.Experiment,'RotaryChair')&contains(all_results.Type,'Sine')&~contains(all_results.Condition,'Light'),:);
                cand_gain_L = NaN(1,length(freqs));
                cand_gain_R = NaN(1,length(freqs));
                cand_phase = NaN(1,length(freqs));
                for j = 1:length(freqs)
                    subtab_L = all_results(find(all_results.Frequency==freqs(j)&...
                        contains(all_results.AxisName,'L'),1,'last'),:);
                    subtab_R = all_results(find(all_results.Frequency==freqs(j)&...
                        contains(all_results.AxisName,'R'),1,'last'),:);
                    if ~isempty(subtab_L)
                        cand_gain_L(j) = subtab_L.Gain;
                        cand_phase(j) = subtab_L.Phase;
                    end
                    if ~isempty(subtab_R)
                        cand_gain_R(j) = subtab_R.Gain;
                        cand_phase(j) = subtab_R.Phase;
                    end
                end
                %% Over Freq
                spread = 0.05;
                mark_size = 7;
                h1 = gobjects(6,1);
                h2 = gobjects(length(all_subs)+2,1);
                fig1 = figure;
                set(fig1,'Units','inches','Position',[1 1 6 4],'Color',[1 1 1])
                ha(1) = subplot(2,1,1);
                plot(NaN,NaN)
                hold on
                h1(6) = fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.gain-norm_dat.gain_std,fliplr(norm_dat.gain+norm_dat.gain_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
                h1(5) = plot(norm_dat.freq,norm_dat.gain,'k--','LineWidth',2);
                h1(1) = errorbar((1-1.5*spread)*freqs,gain_mean.pre,gain_std.pre,'k-','LineWidth',1.5);
                h1(2) = errorbar((1-0.5*spread)*freqs,gain_mean.post,gain_std.post,'-.','Color',0.5*[1,1,1],'LineWidth',1);
                h1(3) = errorbar((1+0.5*spread)*freqs,gain_mean.cro,gain_std.cro,'r:','LineWidth',1);
                h1(4) = errorbar((1+1.5*spread)*freqs,gain_mean.mmo,gain_std.mmo,'r-','LineWidth',1.5);
                for i = 1:length(all_subs)
                    plot((1-0.5*spread)*freqs,gain.post(i,:),'.','Color',0.5*[1,1,1],'Marker',sub_mark(i),'MarkerSize',mark_size)
                    plot((1-1.5*spread)*freqs,gain.pre(i,:),'k.','Marker',sub_mark(i),'MarkerSize',mark_size,'MarkerFaceColor','k')
                    plot((1+1.5*spread)*freqs,gain.mmo(i,:),'r.','Marker',sub_mark(i),'MarkerSize',mark_size,'MarkerFaceColor','r')
                    plot((1+0.5*spread)*freqs,gain.cro(i,:),'r.','Marker',sub_mark(i),'MarkerSize',mark_size)
                end
                plot((1-3*spread)*freqs,cand_gain_L,'b*','MarkerSize',2*mark_size,'MarkerFaceColor','b')
                plot((1-3*spread)*freqs,cand_gain_R,'m*','MarkerSize',2*mark_size,'MarkerFaceColor','m')
                hold off
                set(gca,'xscale','log','xminortick','off','XTick',freqs,'XTickLabel',[],'box','on','Layer','top')
                axis([0.041 1.2 -0.025 0.79])
                ylabel('Horizontal VOR Gain')
                leg1 = legend(h1,{'Pre-Op','Post-Op','MVI OFF','MVI ON','Normal Mean','Normal±SD'},'Location','north','NumColumns',3);
                title(leg1,'Condition')
                leg1.ItemTokenSize(1) = 15;
                ha(2) = subplot(2,1,2);
                plot(NaN,NaN)
                hold on
                fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.phase-norm_dat.phase_std,fliplr(norm_dat.phase+norm_dat.phase_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
                plot(norm_dat.freq,norm_dat.phase,'k--','LineWidth',2);
                errorbar((1-0.5*spread)*freqs,phase_mean.post,phase_std.post,'-.','Color',0.5*[1,1,1],'LineWidth',1);
                errorbar((1-1.5*spread)*freqs,phase_mean.pre,phase_std.pre,'k-','LineWidth',1.5);
                errorbar((1+1.5*spread)*freqs,phase_mean.mmo,phase_std.mmo,'r-','LineWidth',1.5);
                errorbar((1+0.5*spread)*freqs,phase_mean.cro,phase_std.cro,'r:','LineWidth',1);
                for i = 1:length(all_subs)
                    plot((1-0.5*spread)*freqs,phase.post(i,:),'.','Color',0.5*[1,1,1],'Marker',sub_mark(i),'MarkerSize',mark_size)
                    plot((1-1.5*spread)*freqs,phase.pre(i,:),'k.','Marker',sub_mark(i),'MarkerSize',mark_size,'MarkerFaceColor','k')
                    plot((1+1.5*spread)*freqs,phase.mmo(i,:),'r.','Marker',sub_mark(i),'MarkerSize',mark_size,'MarkerFaceColor','r')
                    plot((1+0.5*spread)*freqs,phase.cro(i,:),'r.','Marker',sub_mark(i),'MarkerSize',mark_size)
                    h2(i) = plot(NaN,NaN,'k.','Marker',sub_mark(i),'MarkerSize',mark_size,'LineWidth',1);
                end
                h2(i+1) = plot(NaN,NaN,'b*','MarkerSize',mark_size,'LineWidth',1);
                h2(i+2) = plot(NaN,NaN,'m*','MarkerSize',mark_size,'LineWidth',1);
                plot((1-3*spread)*freqs,cand_phase,'b*','MarkerSize',2*mark_size,'MarkerFaceColor','b')
                hold off
                ha(1).Position = [0.09,0.56,0.98-0.09,0.425];
                ha(2).Position = [0.09,0.12,0.98-0.09,0.425];
                set(gca,'xscale','log','xminortick','off','XTick',freqs,'box','on','Layer','top')
                axis([0.041 1.2 -25 135])
                xlabel('Frequency (Hz)')
                ylabel('Phase Lead (deg)')
                leg2 = legend(h2,[split(cellstr(num2str(1:length(all_subs))));{[candidate,'(L)'];[candidate,'(R)']}],'Location','north','NumColumns',length(all_subs)+2);
                leg2.ItemTokenSize(1) = 8;
                title(leg2,'Subjects')
                %Figure letter labels
                annot_wid_x = 0.04;
                annot_wid_y = 0.06;
                annot_pos_x = (0.09+0.009)*ones(1,2);
                annot_pos_y = [0.985, 0.545]-0.075;
                annot_string = {'A','B'};
                for i = 1:length(annot_string)
                    annotation('textbox',[annot_pos_x(i),annot_pos_y(i),annot_wid_x,annot_wid_y],...
                        'String',annot_string{i},'HorizontalAlignment','center',...
                        'VerticalAlignment','middle','FontWeight','bold','FontSize',20,...
                        'Fitboxtotext','off','BackgroundColor',[1,1,1]);
                end 
                fname1 = [candidate,filesep,'Visit 0',filesep,candidate,'_',datestr(now,'yyyymmdd'),'_SummaryRotaryChairGainPhase_AllMVISub.fig'];
                savefig(fig1,fname1)
                saveas(fig1,strrep(fname1,'.fig','.png'))
            case 'vHIT/GNO'
                
            case 'Candidate: vHIT/GNO'
                
        end
    end
end
end