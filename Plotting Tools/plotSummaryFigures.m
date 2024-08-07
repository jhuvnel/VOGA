function plotSummaryFigures(params)
% Note that the mod/tonic rotary chair data plotted are from the most
% recent visit
% Figure options
opts = {'Rotary Chair Sine','vHIT GNO','Candidate: Rotary Chair Sine','Candidate: vHIT GNO'};
% Rotary Chair figure groupings, change as needed with the addition of more subjects
rotchair_each_unique_fig_text = {'MVI1-5','MVI6-10','MVI11-15'};
rotchair_each_sub_num = [{1:5},{6:10},{11:15}];
% Use MVI001's preop rotary chair testing done elsewhere for his figure
MVI001_preop_freqs = [0.01;0.02;0.04;0.08;0.16;0.32;0.64];
MVI001_preop_gains = [0.0083;0.0250;0.0292;0.0458;0.0500;0.0333;0.0625];
%Subject info
sub_info = params.sub_info;
all_subs = sub_info.Subject;
sub_ears = sub_info.Ear;
% Select figure
[ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
    'SelectionMode','single','ListSize',[200 125],'ListString',opts);
if ~tf
    return;
end
sel = opts{ind};
%Select Candidate if appropriate
if contains(sel,'Candidate')
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
% Make plots
if contains(sel,'Rotary')
    %% Setup
    %Load data
    if ~isfile([params.MVIPath,filesep,'ALLMVI-VOGResults.mat'])
        error(['ALLMVI-VOGResults.mat was not found in ',params.MVIPath])
    end
    load([params.MVIPath,filesep,'ALLMVI-VOGResults.mat'],'all_results')
    all_results = sortrows(sortrows(all_results,'Date','ascend'),'Subject','ascend'); %Sort properly
    all_results = all_results(contains(all_results.Experiment,'RotaryChair')&contains(all_results.Type,'Sine'),:);
    if isempty(all_results)
        error('No Rotary Chair Results files found in this directory')
    end
    % Extract relevant normative values
    load('RotaryChairNormativeData.mat','norm_dat')
    % Limit summary stats to 0.05-1Hz (0.05, 0.1, 0.2, 0.5, 1)
    freqs = [0.05,0.1,0.2,0.5,1];
    remove_phase_with_gain_below = 0.025;
    %Summary matrices
    conds = {'pre','post','cro','mmo'};
    colors = [0,0,0;0.5,0.5,0.5;1,0,1;1,0,0];
    cond_logic = [contains(all_results.Visit,'Visit0')&~contains(all_results.Condition,'Light'),...
        contains(all_results.Visit,'Visit3')&contains(all_results.Condition,'NoStim'),...
        contains(all_results.Condition,'Constant'),contains(all_results.Condition,'Motion')];
    gain = NaN(length(conds),length(freqs),length(all_subs));
    gain_sd= NaN(length(conds),length(freqs),length(all_subs));
    phase = NaN(length(conds),length(freqs),length(all_subs));
    phase_sd = NaN(length(conds),length(freqs),length(all_subs));
    for c = 1:length(conds)
        for i = 1:length(all_subs)
            for j = 1:length(freqs)
                subtab = all_results(find(contains(all_results.Subject,all_subs(i))...
                    &all_results.Frequency==freqs(j)&cond_logic(:,c)&...
                    contains(all_results.AxisName,sub_ears(i)),1,'last'),:);
                if ~isempty(subtab)
                    gain(c,j,i) = subtab.Gain;
                    gain_sd(c,j,i) = subtab.Gain_sd;
                    phase(c,j,i) = subtab.Phase;
                    phase_sd(c,j,i) = subtab.Phase_sd;
                end
            end
        end
    end
    phase(gain<remove_phase_with_gain_below) = NaN;
    gain_mean = mean(gain,3,'omitnan');
    gain_std = std(gain,[],3,'omitnan');
    phase_mean = mean(phase,3,'omitnan');
    phase_std = std(phase,[],3,'omitnan');
    cand_gain_L = NaN(1,length(freqs));
    cand_gain_R = NaN(1,length(freqs));
    cand_phase = NaN(1,length(freqs));
    if contains(sel,'Candidate')
        results_mat = dir([params.MVIPath,filesep,candidate,filesep,'Visit 0',filesep,'Rotary Chair',filesep,'*Results.mat']);
        if isempty(results_mat)
            error('No Results files found in the candidate Rotary Chair directory')
        end
        load([params.MVIPath,filesep,candidate,filesep,'Visit 0',filesep,'Rotary Chair',filesep,results_mat(end).name],'all_results')
        all_results = all_results(contains(all_results.Experiment,'RotaryChair')&contains(all_results.Type,'Sine')&~contains(all_results.Condition,'Light'),:);
        for j = 1:length(freqs)
            subtab_L = all_results(find(all_results.Frequency==freqs(j)&contains(all_results.AxisName,'L'),1,'last'),:);
            subtab_R = all_results(find(all_results.Frequency==freqs(j)&contains(all_results.AxisName,'R'),1,'last'),:);
            if ~isempty(subtab_L)
                cand_gain_L(j) = subtab_L.Gain;
                cand_phase(j) = subtab_L.Phase;
            end
            if ~isempty(subtab_R)
                cand_gain_R(j) = subtab_R.Gain;
                cand_phase(j) = subtab_R.Phase;
            end
        end
    end
    %% Over Frequency with all Subjects and Candidate if Selected
    spread = 0.05;
    h1 = gobjects(6,1);
    fig1 = figure;
    set(fig1,'Units','inches','Position',[1 1 6 4],'Color',[1 1 1])
    ha(1) = subplot(2,1,1);
    plot(NaN,NaN)
    hold on
    h1(6) = fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.gain-norm_dat.gain_std,fliplr(norm_dat.gain+norm_dat.gain_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
    h1(5) = plot(norm_dat.freq,norm_dat.gain,'k--','LineWidth',2);
    for c = 1:length(conds)
        h1(c) = errorbar((1+(c-2.5)*spread)*freqs,gain_mean(c,:),gain_std(c,:),'-','Color',colors(c,:),'LineWidth',1.5);
    end
    for i = 1:length(all_subs)
        rel_mat = gain(:,:,i);
        rel_mat(rel_mat>0.79|rel_mat<-0.025) = NaN; %Out of plot limits
        for c = 1:length(conds)
            text((1+(c-2.5)*spread)*freqs,rel_mat(c,:),char(64+i),'Color',colors(c,:),'FontSize',7)
        end
    end
    h2(1) = plot((1-3*spread)*freqs,cand_gain_L,'b*','MarkerSize',10,'MarkerFaceColor','b');
    h2(2) = plot((1-3*spread)*freqs,cand_gain_R,'g*','MarkerSize',10,'MarkerFaceColor','g');
    hold off
    set(gca,'xscale','log','xminortick','off','XTick',freqs,'XTickLabel',[],'box','on','Layer','top')
    axis([0.041 1.2 -0.025 0.79])
    ylabel('Horizontal VOR Gain')
    title('Rotary Chair Horizontal Sinusoids (Most Recent Visit)')
    leg1 = legend(h1,{'Pre-op','Post-op','Placebo','Treatment','Norm Mean','Norm±SD'},'Location','north','NumColumns',length(h1));
    leg1.ItemTokenSize(1) = 15;
    ha(2) = subplot(2,1,2);
    plot(NaN,NaN)
    hold on
    fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.phase-norm_dat.phase_std,fliplr(norm_dat.phase+norm_dat.phase_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
    plot(norm_dat.freq,norm_dat.phase,'k--','LineWidth',2);
    for c = 1:length(conds)
        errorbar((1+(c-2.5)*spread)*freqs,phase_mean(c,:),phase_std(c,:),'-','Color',colors(c,:),'LineWidth',1.5)
    end
    for i = 1:length(all_subs)
        rel_mat = phase(:,:,i);
        rel_mat(rel_mat>135|rel_mat<-25) = NaN; %Out of plot limits
        for c = 1:length(conds)
            text((1+(c-2.5)*spread)*freqs,rel_mat(c,:),char(64+i),'Color',colors(c,:),'FontSize',7)
        end
    end
    plot((1-3*spread)*freqs,cand_phase,'b*','MarkerSize',10,'MarkerFaceColor','b')
    hold off
    set(gca,'xscale','log','xminortick','off','XTick',freqs,'box','on','Layer','top')
    axis([0.041 1.2 -25 135])
    xlabel('Frequency (Hz)')
    ylabel('Phase Lead (deg)')
    ha(1).Position = [0.09,0.54,0.98-0.09,0.41];
    ha(2).Position = [0.09,0.12,0.98-0.09,0.41];
    if contains(sel,'Candidate')
        leg2 = legend(h2,{[candidate,' Left Ear'],[candidate,' Right Ear']},'Location','north','NumColumns',2);
        leg2.ItemTokenSize(1) = 10;
        fname1 = [candidate,filesep,'Visit 0',filesep,candidate,'_',char(datetime('now','Format','yyyyMMdd')),'_SummaryRotaryChairGainPhase_AllSub.fig'];
        savefig(fig1,fname1)
        saveas(fig1,strrep(fname1,'.fig','.png'))
        return;
    end
    fname1 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryRotaryChairGainPhase_AllSub.fig'];
    savefig(fig1,fname1)
    saveas(fig1,strrep(fname1,'.fig','.png'))
    %% Over Frequency Grouped by Sub
    for f = 1:length(rotchair_each_unique_fig_text)
        unique_fig_text = rotchair_each_unique_fig_text{f};
        sub_num = rotchair_each_sub_num{f}; %Subjects to graph (array of #)
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
            for c = 1:length(conds)
                h1(c) = errorbar((1+(c-2.5)*spread)*freqs,gain(c,:,sub_num(i)),gain_sd(c,:,sub_num(i)),'-','Color',colors(c,:),'LineWidth',1.5);
            end
            if sub_num(i)==1 %Add MVI001 pre-op done elsewhere as black X
                plot(MVI001_preop_freqs,MVI001_preop_gains,'kx','LineWidth',2)
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
            for c = 1:length(conds)
                errorbar((1+(c-2.5)*spread)*freqs,phase(c,:,sub_num(i)),phase_sd(c,:,sub_num(i)),'-','Color',colors(c,:),'LineWidth',1.5)
            end
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
        leg = legend(ha(1,2),h1,{'Pre-op','Post-op','Placebo','Treatment','Norm Mean','Norm±SD'},'NumColumns',2,'Location','northeast');
        leg.ItemTokenSize(1) = 15;
        title(leg,'Condition')
        fname2 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryRotaryChairGainPhase_EachSub_',unique_fig_text,'.fig'];
        savefig(fig2,fname2)
        saveas(fig2,strrep(fname2,'fig','png'))
    end
elseif contains(sel,'vHIT')
    %% Setup
    %Load data
    if ~isfile([params.MVIPath,filesep,'vHIT_GNO.mat'])
        error(['vHIT_GNO.mat was not found in ',params.MVIPath])
    end
    load([params.MVIPath,filesep,'vHIT_GNO.mat'],'val_results')
    val_results = sortrows(sortrows(val_results,'Date','ascend'),'Subject','ascend'); %Sort properly
    val_results = val_results(contains(val_results.Side,'ipsi')&~isnan(val_results.Gain),:);
    cond_logic = [strcmp(val_results.Visit,'0'),...
        strcmp(val_results.Visit,'9x')&contains(val_results.Condition,'Motion'),...
        strcmp(val_results.Visit,'10x')&contains(val_results.Condition,'Motion'),...
        strcmp(val_results.Visit,'11x')&contains(val_results.Condition,'Motion'),...
        contains(val_results.Condition,'Motion'),...
        strcmp(val_results.Visit,'9x')&contains(val_results.Condition,'Constant'),...
        strcmp(val_results.Visit,'10x')&contains(val_results.Condition,'Constant'),...
        strcmp(val_results.Visit,'11x')&contains(val_results.Condition,'Constant'),...
        contains(val_results.Condition,'Constant')];
    NC = size(cond_logic,2);
    canals = {'Horizontal','Anterior','Posterior'};
    gain = NaN(NC,length(canals),length(all_subs));
    for c = 1:NC
        for i = 1:length(all_subs)
            for j = 1:length(canals)
                subtab = sortrows(val_results(contains(val_results.Subject,all_subs(i))...
                    &contains(val_results.Canal,canals{j}(1))&cond_logic(:,c),:),'Type','ascend');
                if ~isempty(subtab)
                    gain(c,j,i) = subtab.Gain(end);
                end
            end
        end
    end
    dgain = gain-gain(1,:,:);
    dgain_mean = mean(dgain,3,'omitnan');
    dgain_std = std(dgain,[],3,'omitnan');
    %% Visit 0 Gain Values and Candidate values
    if contains(sel,'Candidate')
        if ~isfile(strrep([params.MVIPath,filesep,candidate,'/Visit 0/vHIT/GNO_Summary.mat'],'/',filesep))
            error('No GNO_Summary.mat file found in the candidate vHIT directory.')
        end
        load(strrep([params.MVIPath,filesep,candidate,'/Visit 0/vHIT/GNO_Summary.mat'],'/',filesep),'tab')
        cand_gain = NaN(2,length(canals));
        for j = 1:length(canals)
            subtab_L = sortrows(tab(strcmp(tab.Canal,['L',canals{j}(1)]),:),'Type','ascend');
            subtab_R = sortrows(tab(strcmp(tab.Canal,['R',canals{j}(1)]),:),'Type','ascend');
            if ~isempty(subtab_L)
                cand_gain(1,j) = subtab_L.Gain(end);
            end
            if ~isempty(subtab_R)
                cand_gain(2,j) = subtab_R.Gain(end);
            end
        end
        rel_mat = permute(gain(1,:,:),[3,2,1]); %v0 of the MVI subjects
        fig1 = figure;
        set(fig1,'Units','inches','Position',[1 1 4 4],'Color',[1 1 1])
        spread = 0.15;
        for i = 1:length(all_subs)
            for c = 1:length(canals)
                text(c+rand*spread,rel_mat(i,c),char(64+i),'FontSize',7)
            end
        end
        hold on
        h2(1) = plot((1:3)-0.1,cand_gain(1,:),'b*','MarkerSize',10,'MarkerFaceColor','b');
        h2(2) = plot((1:3)-0.10,cand_gain(2,:),'g*','MarkerSize',10,'MarkerFaceColor','g');
        xline(1.5)
        xline(2.5)
        hold off
        set(gca,'xminortick','off','XTick',1:3,'XTickLabels',{'Horizontal','Anterior','Posterior'},'box','on')
        axis([0.5 3.5 -0.1 1.1])
        xlabel('Pre-Op Value for Implanted Canals')
        ylabel('Video Head Impulse Gain (GNO)')
        leg2 = legend(h2,{[candidate,' Left Ear'],[candidate,' Right Ear']},'Location','north','NumColumns',2);
        leg2.ItemTokenSize(1) = 10;
        fname1 = [candidate,filesep,'Visit 0',filesep,candidate,'_',char(datetime('now','Format','yyyyMMdd')),'_GNOvHITGain_v0AllSub.fig'];
        savefig(fig1,fname1)
        saveas(fig1,strrep(fname1,'.fig','.png'))
    end
    %% Implanted Canals Over Time with all Subjects
    titles = {'Treatment Mode','Placebo Mode'};
    fig1 = figure;
    set(fig1,'Units','inches','Position',[1 1 6 4],'Color',[1 1 1])
    %Set figure and axes sizing
    nr = length(canals);
    nc = length(titles);
    x_spac = 0.01;
    x_min = 0.10;
    x_max = 1-x_spac;
    y_min = 0.17;
    y_max = 0.95;
    y_spac = 0.02;
    y_wid = (y_max-y_min-y_spac*(nr-1))/nr;
    x_wid = (x_max-x_min-x_spac*(nc-1))/nc;
    x_pos = x_min:(x_wid+x_spac):x_max;
    y_pos = fliplr(y_min:(y_wid+y_spac):y_max);
    spread = 0.1;
    ha = gobjects(nr,nc);
    h1 = gobjects(3,1);
    for i = 1:nr
        for j = 1:nc
            ha(i,j) = subplot(nr,nc,sub2ind([nc,nr],j,i));
            h1(1) = plot(NaN,NaN,'k-','LineWidth',1.5);
            hold on
            h1(3) = fill([0,6,6,0],0.06*[-1,-1,1,1],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
            h1(2) = plot([0,6],[0,0],'k:','LineWidth',2);
            errorbar(1:5,dgain_mean([1,(-2:1)+4*j],i),dgain_std([1,(-2:1)+4*j],i),'k-','LineWidth',1.5);
            for s = 1:length(all_subs)
                text((1:5)+(rand-1.5)*spread,dgain([1,(-2:1)+4*j],i,s),char(64+s),'FontSize',7)
            end
            hold off
        end
        ylabel(ha(i,1),{'\Delta vHIT Gain',canals{i}})
    end
    for j = 1:nc
        for i = 1:nr
            ha(i,j).Position = [x_pos(j),y_pos(i),x_wid,y_wid];
        end
        title(ha(1,j),titles{j})
    end
    set(ha,'xminortick','off','XTick',1:5,'XTickLabels',{'0','0.5','1','2','most recent'},...
        'box','on','XLim',[0.75 5.25],'YLim',[-0.6 0.7],'Layer','top')
    set(ha(:,2:end),'YTickLabel',[])
    set(ha(1:end-1,:),'XTickLabel',[])
    xlabel(ha(end,:),'Years Since Implantation')
    leg1 = legend(ha(1,end),h1,{'Mean±SD Change','No Change','MCID'},'Location','north','NumColumns',length(h1));
    leg1.ItemTokenSize(1) = 10;   
    fname1 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryGNOvHITGain_AllSub.fig'];
    savefig(fig1,fname1)
    saveas(fig1,strrep(fname1,'.fig','.png'))
end
end