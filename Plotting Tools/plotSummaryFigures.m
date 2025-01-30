function plotSummaryFigures(params)
% Note that the mod/tonic rotary chair data plotted are from the most
% recent visit
% Figure options
opts = {'Rotary Chair Sine','vHIT GNO','Candidate: Rotary Chair Sine','Candidate: vHIT GNO'};
% Rotary Chair figure groupings, change as needed with the addition of more subjects
rotchair_each_unique_fig_text = {'MVI1-5','MVI6-10','MVI11-15', 'MVI16-17'};
rotchair_each_sub_num = [{1:5},{6:10},{11:15},{16:17}];
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
    conds = {'pre','post','cr6mo','mm6mo','cro','mmo'};
    colors = [0,0,0;0.5,0.5,0.5;1,0,1;1,0,0];
    cond_logic = [contains(all_results.Visit,'Visit0')&~contains(all_results.Condition,'Light'),...
        contains(all_results.Visit,'Visit3')&contains(all_results.Condition,'NoStim'),...
        contains(all_results.Visit,'Visit9x')&contains(all_results.Condition,'Constant'),...
        contains(all_results.Visit,'Visit9x')&contains(all_results.Condition,'Motion'),...
        contains(all_results.Condition,'Constant'),...
        contains(all_results.Condition,'Motion')];
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

    % median/CI across subjects
    gain_med = median(gain,3,'omitnan'); % median across subjects
    
    sz = size(gain);
    n_gain = nan([sz(1) sz(2)]);
    gain_lowCI = nan([sz(1) sz(2)]);
    gain_highCI = nan([sz(1) sz(2)]);
    gain_lowCIwid = nan([sz(1) sz(2)]);
    gain_highCIwid = nan([sz(1) sz(2)]);
    for i = 1:sz(1)
        for j = 1:sz(2)
            n_gain(i,j) = median95CI(gain(i,j,:),'n'); % n grouping subject
            gain_lowCI(i,j) = median95CI(gain(i,j,:),'lowCI'); % CI for each canal grouping subjects
            gain_highCI(i,j) = median95CI(gain(i,j,:),'highCI');
            
            gain_lowCIwid(i,j) = gain_med(i,j) - gain_lowCI(i,j); % CI width for each canal grouping subjects
            gain_highCIwid(i,j) = gain_highCI(i,j) - gain_med(i,j);
        end
    end

    phase_med = median(phase,3,'omitnan'); % median across subjects
    
    sz = size(phase);
    n_phase = nan([sz(1) sz(2)]);
    phase_lowCI = nan([sz(1) sz(2)]);
    phase_highCI = nan([sz(1) sz(2)]);
    phase_lowCIwid = nan([sz(1) sz(2)]);
    phase_highCIwid = nan([sz(1) sz(2)]);
    for i = 1:sz(1)
        for j = 1:sz(2)
            n_phase(i,j) = median95CI(phase(i,j,:),'n'); % n grouping subject
            phase_lowCI(i,j) = median95CI(phase(i,j,:),'lowCI'); % CI for each canal grouping subjects
            phase_highCI(i,j) = median95CI(phase(i,j,:),'highCI');
            
            phase_lowCIwid(i,j) = phase_med(i,j) - phase_lowCI(i,j); % CI width for each canal grouping subjects
            phase_highCIwid(i,j) = phase_highCI(i,j) - phase_med(i,j);
        end
    end

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
    %% Over Frequency with all Subjects and Candidate if Selected - Visit 9x (6 month)
    condsSubset = [1 2 3 4];% refers to the conditions in conds - [Preop, Postop, 9x constant rate, 9x motion mod]
    spread = 0.05;
    h1 = gobjects(6,1);
    fig1 = figure;
    set(fig1,'Units','inches','Position',[1 1 6 4],'Color',[1 1 1])
    ha(1) = subplot(2,1,1);
    plot(NaN,NaN)
    hold on
    h1(6) = fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.gain-norm_dat.gain_std,fliplr(norm_dat.gain+norm_dat.gain_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
    h1(5) = plot(norm_dat.freq,norm_dat.gain,'k--','LineWidth',2);
    for c = 1:length(condsSubset)
        % h1(c) = errorbar((1+(c-2.5)*spread)*freqs,gain_mean(condsSubset(c),:),gain_std(condsSubset(c),:),'-','Color',colors(c,:),'LineWidth',1.5);
        h1(c) = errorbar((1+(c-2.5)*spread)*freqs,gain_med(condsSubset(c),:),gain_lowCIwid(condsSubset(c),:),gain_highCIwid(condsSubset(c),:),'-','Color',colors(c,:),'LineWidth',1.5);
    end
    for i = 1:length(all_subs)
        rel_mat = gain(:,:,i);
        rel_mat(rel_mat>0.79|rel_mat<-0.025) = NaN; %Out of plot limits
        for c = 1:length(condsSubset)
            text((1+(c-2.5)*spread)*freqs,rel_mat(condsSubset(c),:),char(64+i),'Color',colors(c,:),'FontSize',7)
        end
    end
    h2(1) = plot((1-3*spread)*freqs,cand_gain_L,'b*','MarkerSize',10,'MarkerFaceColor','b');
    h2(2) = plot((1-3*spread)*freqs,cand_gain_R,'g*','MarkerSize',10,'MarkerFaceColor','g');
    hold off
    set(gca,'xscale','log','xminortick','off','XTick',freqs,'XTickLabel',[],'box','on','Layer','top')
    axis([0.041 1.2 -0.025 0.79])
    ylabel('Horizontal VOR Gain')
    title('Rotary Chair Horizontal Sinusoids (6 Months)')
    % leg1 = legend(h1,{'Pre-op','Post-op','Tonic','Treatment','Norm Mean','Norm±SD'},'Location','north','NumColumns',length(h1));
    leg1 = legend(h1,{'Pre-op','Post-op','Tonic','Treatment','Norm Mean','Norm±SD'},'Location','north','NumColumns',length(h1));
    leg1.ItemTokenSize(1) = 15;
    ha(2) = subplot(2,1,2);
    plot(NaN,NaN)
    hold on
    fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.phase-norm_dat.phase_std,fliplr(norm_dat.phase+norm_dat.phase_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
    plot(norm_dat.freq,norm_dat.phase,'k--','LineWidth',2);
    for c = 1:length(condsSubset)
        % errorbar((1+(c-2.5)*spread)*freqs,phase_mean(condsSubset(c),:),phase_std(condsSubset(c),:),'-','Color',colors(c,:),'LineWidth',1.5)
        errorbar((1+(c-2.5)*spread)*freqs,phase_med(condsSubset(c),:),phase_lowCIwid(condsSubset(c),:),phase_highCIwid(condsSubset(c),:),'-','Color',colors(c,:),'LineWidth',1.5);
    end
    for i = 1:length(all_subs)
        rel_mat = phase(:,:,i);
        rel_mat(rel_mat>135|rel_mat<-25) = NaN; %Out of plot limits
        for c = 1:length(condsSubset)
            text((1+(c-2.5)*spread)*freqs,rel_mat(condsSubset(c),:),char(64+i),'Color',colors(c,:),'FontSize',7)
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
        fname1 = [candidate,filesep,'Visit 0',filesep,candidate,'_',char(datetime('now','Format','yyyyMMdd')),'_SummaryRotaryChairGainPhase_AllSub_9x.fig'];
        savefig(fig1,fname1)
        saveas(fig1,strrep(fname1,'.fig','.png'))
    end
    fname1 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryRotaryChairGainPhase_AllSub_9x.fig'];
    savefig(fig1,fname1)
    saveas(fig1,strrep(fname1,'.fig','.svg'))

    %% Over Frequency with all Subjects and Candidate if Selected - Most Recent Visit
    condsSubset = [1 2 5 6]; % refers to the conditions in conds - [Preop, Postop, most recent constant rate, most recent motion mod]
    spread = 0.05;
    h1 = gobjects(6,1);
    fig1 = figure;
    set(fig1,'Units','inches','Position',[1 1 6 4],'Color',[1 1 1])
    ha(1) = subplot(2,1,1);
    plot(NaN,NaN)
    hold on
    h1(6) = fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.gain-norm_dat.gain_std,fliplr(norm_dat.gain+norm_dat.gain_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
    h1(5) = plot(norm_dat.freq,norm_dat.gain,'k--','LineWidth',2);
    for c = 1:length(condsSubset)
        % h1(c) = errorbar((1+(c-2.5)*spread)*freqs,gain_mean(condsSubset(c),:),gain_std(condsSubset(c),:),'-','Color',colors(c,:),'LineWidth',1.5);
        h1(c) = errorbar((1+(c-2.5)*spread)*freqs,gain_med(condsSubset(c),:),gain_lowCIwid(condsSubset(c),:),gain_highCIwid(condsSubset(c),:),'-','Color',colors(c,:),'LineWidth',1.5);
    end
    for i = 1:length(all_subs)
        rel_mat = gain(:,:,i);
        rel_mat(rel_mat>0.79|rel_mat<-0.025) = NaN; %Out of plot limits
        for c = 1:length(condsSubset)
            text((1+(c-2.5)*spread)*freqs,rel_mat(condsSubset(c),:),char(64+i),'Color',colors(c,:),'FontSize',7)
        end
    end
    h2(1) = plot((1-3*spread)*freqs,cand_gain_L,'b*','MarkerSize',10,'MarkerFaceColor','b');
    h2(2) = plot((1-3*spread)*freqs,cand_gain_R,'g*','MarkerSize',10,'MarkerFaceColor','g');
    hold off
    set(gca,'xscale','log','xminortick','off','XTick',freqs,'XTickLabel',[],'box','on','Layer','top')
    axis([0.041 1.2 -0.025 0.79])
    ylabel('Horizontal VOR Gain')
    title('Rotary Chair Horizontal Sinusoids (Most Recent Visit)')
    leg1 = legend(h1,{'Pre-op','Post-op','Tonic','Treatment','Norm Mean','Norm±SD'},'Location','north','NumColumns',length(h1));
    leg1.ItemTokenSize(1) = 15;
    ha(2) = subplot(2,1,2);
    plot(NaN,NaN)
    hold on
    fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.phase-norm_dat.phase_std,fliplr(norm_dat.phase+norm_dat.phase_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
    plot(norm_dat.freq,norm_dat.phase,'k--','LineWidth',2);
    for c = 1:length(condsSubset)
        % errorbar((1+(c-2.5)*spread)*freqs,phase_mean(condsSubset(c),:),phase_std(condsSubset(c),:),'-','Color',colors(c,:),'LineWidth',1.5)
        errorbar((1+(c-2.5)*spread)*freqs,phase_med(condsSubset(c),:),phase_lowCIwid(condsSubset(c),:),phase_highCIwid(condsSubset(c),:),'-','Color',colors(c,:),'LineWidth',1.5);
    end
    for i = 1:length(all_subs)
        rel_mat = phase(:,:,i);
        rel_mat(rel_mat>135|rel_mat<-25) = NaN; %Out of plot limits
        for c = 1:length(condsSubset)
            text((1+(c-2.5)*spread)*freqs,rel_mat(condsSubset(c),:),char(64+i),'Color',colors(c,:),'FontSize',7)
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
        fname1 = [candidate,filesep,'Visit 0',filesep,candidate,'_',char(datetime('now','Format','yyyyMMdd')),'_SummaryRotaryChairGainPhase_AllSub_MostRecent.fig'];
        savefig(fig1,fname1)
        saveas(fig1,strrep(fname1,'.fig','.png'))
        return;
    end
    fname1 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryRotaryChairGainPhase_AllSub_MostRecent.fig'];
    savefig(fig1,fname1)
    saveas(fig1,strrep(fname1,'.fig','.svg'))
    %% Over Frequency Grouped by Sub
    condsSubset = [1 2 5 6];
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
            for c = 1:length(condsSubset)
                h1(c) = errorbar((1+(c-2.5)*spread)*freqs,gain(condsSubset(c),:,sub_num(i)),gain_sd(condsSubset(c),:,sub_num(i)),'-','Color',colors(c,:),'LineWidth',1.5);
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
            for c = 1:length(condsSubset)
                errorbar((1+(c-2.5)*spread)*freqs,phase(condsSubset(c),:,sub_num(i)),phase_sd(condsSubset(c),:,sub_num(i)),'-','Color',colors(c,:),'LineWidth',1.5)
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
        leg = legend(ha(1,2),h1,{'Pre-op','Post-op','Tonic','Treatment','Norm Mean','Norm±SD'},'NumColumns',2,'Location','northeast');
        leg.ItemTokenSize(1) = 15;
        title(leg,'Condition')
        fname2 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryRotaryChairGainPhase_EachSub_',unique_fig_text,'.fig'];
        savefig(fig2,fname2)
        saveas(fig2,strrep(fname2,'fig','svg'))
    end
elseif contains(sel,'vHIT')
    %% Setup
    %Load data
    if ~isfile([params.MVIPath,filesep,'vHIT_GNO.mat'])
        error(['vHIT_GNO.mat was not found in ',params.MVIPath])
    end
    load([params.MVIPath,filesep,'vHIT_GNO.mat'],'val_results')
    val_results = sortrows(sortrows(val_results,'Date','ascend'),'Subject','ascend'); %Sort properly
    val_results = val_results(~isnan(val_results.Gain),:);
    cond_logic = [strcmp(val_results.Visit,'0'),...
        strcmp(val_results.Visit,'3')&contains(val_results.Condition,'NoStim'),...
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
    sides = {'ipsi','contra'};
    gain = NaN(NC,length(canals),length(all_subs), length(sides));
    lat = NaN(NC,length(canals),length(all_subs), length(sides));
    for c = 1:NC
        for i = 1:length(all_subs)
            for j = 1:length(canals)
                for k = 1:length(sides)
                    subtab = sortrows(val_results(contains(val_results.Subject,all_subs(i))...
                        &contains(val_results.Canal,canals{j}(1))&cond_logic(:,c)&contains(val_results.Side,sides(k)),:),'Type','ascend');
                    if ~isempty(subtab)
                        gain(c,j,i,k) = subtab.Gain(end);
                        lat(c,j,i,k) = subtab.Lat(end);
                    end
                end
            end
        end
    end
    dgain = gain-gain(1,:,:,:); % change in gain from preop
    dlat = lat-lat(1,:,:,:); % change in latency from preop

    dgain_mean = mean(dgain,3,'omitnan');
    dgain_mean2 = mean(dgain,[2 3],'omitnan');
    dgain_std = std(dgain,[],3,'omitnan');
    dgain_std2 = std(dgain,[],[2 3],'omitnan');
    dgain_sem = std(dgain,[],3,'omitnan')./sqrt(sum(~isnan(dgain),3));
    dgain_sem2 = std(dgain,[],[2 3],'omitnan')./sqrt(sum(~isnan(dgain),[2 3]));
    
    dlat_mean = mean(dlat,3,'omitnan');
    dlat_mean2 = mean(dlat,[2 3],'omitnan');
    dlat_std = std(dlat,[],3,'omitnan');
    dlat_std2 = std(dlat,[],[2 3],'omitnan');
    dlat_sem = std(dlat,[],3,'omitnan')./sqrt(sum(~isnan(dlat),3));
    dlat_sem2 = std(dlat,[],[2 3],'omitnan')./sqrt(sum(~isnan(dlat),[2 3]));


    % median/CI across subjects only (each canal will be plotted
    % separately)
    dgain_med1 = median(dgain,3,'omitnan'); % median across subjects
    
    sz = size(dgain);
    n_dgain1 = nan([sz(1) sz(2) 1 sz(4)]);
    dgain_lowCI1 = nan([sz(1) sz(2) 1 sz(4)]);
    dgain_highCI1 = nan([sz(1) sz(2) 1 sz(4)]);
    dgain_lowCIwid1 = nan([sz(1) sz(2) 1 sz(4)]);
    dgain_highCIwid1 = nan([sz(1) sz(2) 1 sz(4)]);
    for i = 1:sz(1)
        for j = 1:sz(2)
            for l = 1:sz(4)
                n_dgain1(i,j,1,l) = median95CI(dgain(i,j,:,l),'n'); % n grouping subject
                dgain_lowCI1(i,j,1,l) = median95CI(dgain(i,j,:,l),'lowCI'); % CI for each canal grouping subjects
                dgain_highCI1(i,j,1,l) = median95CI(dgain(i,j,:,l),'highCI');
                
                dgain_lowCIwid1(i,j,1,l) = dgain_med1(i,j,1,l) - dgain_lowCI1(i,j,1,l); % CI width for each canal grouping subjects
                dgain_highCIwid1(i,j,1,l) = dgain_highCI1(i,j,1,l) - dgain_med1(i,j,1,l);
            end
        end
    end

    dlat_med1 = median(dlat,3,'omitnan'); % median across subjects
    
    sz = size(dlat);
    n_dlat1 = nan([sz(1) sz(2) 1 sz(4)]);
    dlat_lowCI1 = nan([sz(1) sz(2) 1 sz(4)]);
    dlat_highCI1 = nan([sz(1) sz(2) 1 sz(4)]);
    dlat_lowCIwid1 = nan([sz(1) sz(2) 1 sz(4)]);
    dlat_highCIwid1 = nan([sz(1) sz(2) 1 sz(4)]);
    for i = 1:sz(1)
        for j = 1:sz(2)
            for l = 1:sz(4)
                n_dlat1(i,j,1,l) = median95CI(dlat(i,j,:,l),'n'); % n grouping subject
                dlat_lowCI1(i,j,1,l) = median95CI(dlat(i,j,:,l),'lowCI'); % CI for each canal grouping subjects
                dlat_highCI1(i,j,1,l) = median95CI(dlat(i,j,:,l),'highCI');
                
                dlat_lowCIwid1(i,j,1,l) = dlat_med1(i,j,1,l) - dlat_lowCI1(i,j,1,l); % CI width for each canal grouping subjects
                dlat_highCIwid1(i,j,1,l) = dlat_highCI1(i,j,1,l) - dlat_med1(i,j,1,l);
            end
        end
    end


    % mean across canals, then take median/CI across subjects
    dgain_mean1 = mean(dgain,2,'omitnan'); % mean across canals
    dgain_med2 = median(dgain_mean1,3,'omitnan'); % mean across subjects after mean across canals

    n_dgain2 = nan([sz(1) 1 1 sz(4)]);
    dgain_lowCI2 = nan([sz(1) 1 1 sz(4)]);
    dgain_highCI2 = nan([sz(1) 1 1 sz(4)]);
    dgain_lowCIwid2 = nan([sz(1) 1 1 sz(4)]);
    dgain_highCIwid2 = nan([sz(1) 1 1 sz(4)]);
    for i = 1:sz(1)
        for l = 1:sz(4)
            n_dgain2(i,1,1,l) = median95CI(dgain_mean1(i,1,:,l),'n'); % n grouping subject
            dgain_lowCI2(i,1,1,l) = median95CI(dgain_mean1(i,1,:,l),'lowCI'); % CI for each canal grouping subjects
            dgain_highCI2(i,1,1,l) = median95CI(dgain_mean1(i,1,:,l),'highCI');
            
            dgain_lowCIwid2(i,1,1,l) = dgain_med2(i,1,1,l) - dgain_lowCI2(i,1,1,l); % CI width for each canal grouping subjects
            dgain_highCIwid2(i,1,1,l) = dgain_highCI2(i,1,1,l) - dgain_med2(i,1,1,l);
        end
    end

    dlat_mean1 = mean(dlat,2,'omitnan'); % mean across canals
    dlat_med2 = median(dlat_mean1,3,'omitnan'); % mean across subjects after mean across canals

    n_dlat2 = nan([sz(1) 1 1 sz(4)]);
    dlat_lowCI2 = nan([sz(1) 1 1 sz(4)]);
    dlat_highCI2 = nan([sz(1) 1 1 sz(4)]);
    dlat_lowCIwid2 = nan([sz(1) 1 1 sz(4)]);
    dlat_highCIwid2 = nan([sz(1) 1 1 sz(4)]);
    for i = 1:sz(1)
        for l = 1:sz(4)
            n_dlat2(i,1,1,l) = median95CI(dlat_mean1(i,1,:,l),'n'); % n grouping subject
            dlat_lowCI2(i,1,1,l) = median95CI(dlat_mean1(i,1,:,l),'lowCI'); % CI for each canal grouping subjects
            dlat_highCI2(i,1,1,l) = median95CI(dlat_mean1(i,1,:,l),'highCI');
            
            dlat_lowCIwid2(i,1,1,l) = dlat_med2(i,1,1,l) - dlat_lowCI2(i,1,1,l); % CI width for each canal grouping subjects
            dlat_highCIwid2(i,1,1,l) = dlat_highCI2(i,1,1,l) - dlat_med2(i,1,1,l);
        end
    end
    
    % median/CI across subjects and canals at the same time
    dgain_med3 = median(dgain, [2 3], 'omitnan'); % median across subjects and canals

    n_dgain3 = nan([sz(1) 1 1 sz(4)]);
    dgain_lowCI3 = nan([sz(1) 1 1 sz(4)]);
    dgain_highCI3 = nan([sz(1) 1 1 sz(4)]);
    dgain_lowCIwid3 = nan([sz(1) 1 1 sz(4)]);
    dgain_highCIwid3 = nan([sz(1) 1 1 sz(4)]);
    for i = 1:sz(1)
        for l = 1:sz(4)
            n_dgain3(i,1,1,l) = median95CI(dgain(i,:,:,l),'n'); % n grouping subject
            dgain_lowCI3(i,1,1,l) = median95CI(dgain(i,:,:,l),'lowCI'); % CI for each canal grouping subjects
            dgain_highCI3(i,1,1,l) = median95CI(dgain(i,:,:,l),'highCI');
            
            dgain_lowCIwid3(i,1,1,l) = dgain_med3(i,1,1,l) - dgain_lowCI3(i,1,1,l); % CI width for each canal grouping subjects
            dgain_highCIwid3(i,1,1,l) = dgain_highCI3(i,1,1,l) - dgain_med3(i,1,1,l);
        end
    end

    dlat_med3 = median(dlat, [2 3], 'omitnan'); % median across subjects and canals

    n_dlat3 = nan([sz(1) 1 1 sz(4)]);
    dlat_lowCI3 = nan([sz(1) 1 1 sz(4)]);
    dlat_highCI3 = nan([sz(1) 1 1 sz(4)]);
    dlat_lowCIwid3 = nan([sz(1) 1 1 sz(4)]);
    dlat_highCIwid3 = nan([sz(1) 1 1 sz(4)]);
    for i = 1:sz(1)
        for l = 1:sz(4)
            n_dlat3(i,1,1,l) = median95CI(dlat(i,:,:,l),'n'); % n grouping subject
            dlat_lowCI3(i,1,1,l) = median95CI(dlat(i,:,:,l),'lowCI'); % CI for each canal grouping subjects
            dlat_highCI3(i,1,1,l) = median95CI(dlat(i,:,:,l),'highCI');
            
            dlat_lowCIwid3(i,1,1,l) = dlat_med3(i,1,1,l) - dlat_lowCI3(i,1,1,l); % CI width for each canal grouping subjects
            dlat_highCIwid3(i,1,1,l) = dlat_highCI3(i,1,1,l) - dlat_med3(i,1,1,l);
        end
    end
    

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
        rel_mat = permute(gain(1,:,:,1),[3,2,1]); %v0 of the MVI subjects
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
    %% Implanted Canals Change from Preop Gain Over Time with all Subjects
    
    titles = {'Treatment Mode','Tonic Mode'};
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
            errorbar(1:4,dgain_med1([1,(-1:1)+4*j],i,1,1),dgain_lowCIwid1([1,(-1:1)+4*j],i,1,1),dgain_highCIwid1([1,(-1:1)+4*j],i,1,1),'k-','LineWidth',1.5);
            errorbar(5,dgain_med1(2+4*j,i,1,1),dgain_lowCIwid1(2+4*j,i,1,1),dgain_highCIwid1(2+4*j,i,1,1),'k','LineWidth',1.5);
            for s = 1:length(all_subs)
                text((1:5)+(rand-1.5)*spread,dgain([1,(-1:2)+4*j],i,s,1),char(64+s),'FontSize',7)
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
    leg1 = legend(ha(1,end),h1,{'Med±95%CI','No Change','MCID'},'Location','north','NumColumns',length(h1));
    leg1.ItemTokenSize(1) = 10;   
    fname1 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryGNOvHITChangeGain_AllSub_SepAxisCanals.fig'];
    savefig(fig1,fname1)
    saveas(fig1,strrep(fname1,'.fig','.svg'))

    %% Implanted Canals (all canals on same axis) Change from Preop Gain Over Time with all Subjects
    % Take mean across all 3 canals for each subject, then take med/CI of
    % the mean vHIT gain for each subject
    
    colors = [0.75,0,0;0,0.25,0;0,0,0.5];
    titles = {'Mean \Delta vHIT Gain'};
    ylabels = {'Treatment Mode','Tonic Mode'};
    fig2 = figure;
    set(fig2,'Units','inches','Position',[1 1 6 4],'Color',[1 1 1])
    %Set figure and axes sizing
    nr = length(ylabels);
    nc = 1;
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
        ha(i,1) = subplot(nr,nc,sub2ind([nr,nc],i,1));
        h1(1) = plot(NaN,NaN,'k-','LineWidth',1.5); % mean line for legend
        hold on
        h1(3) = fill([0,6,6,0],0.06*[-1,-1,1,1],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]); % fill MCID area
        h1(2) = plot([0,6],[0,0],'k:','LineWidth',2); % no change line

        % h1(4) = plot(NaN,NaN,'-','Color',colors(1,:),'LineWidth',2); % canal colors for legend
        % h1(5) = plot(NaN,NaN,'-','Color',colors(2,:),'LineWidth',2);
        % h1(6) = plot(NaN,NaN,'-','Color',colors(3,:),'LineWidth',2);
    
        errorbar(1:4,dgain_med2([1,(-1:1)+4*i],1,1,1),dgain_lowCIwid2([1,(-1:1)+4*i],1,1,1),dgain_highCIwid2([1,(-1:1)+4*i],1,1,1),'k-','LineWidth',1.5);
        errorbar(5,dgain_med2(2+4*i,1,1,1),dgain_lowCIwid2(2+4*i,1,1,1),dgain_highCIwid2(2+4*i,1,1,1),'k','LineWidth',1.5);
        % for j = 1:length(canals)
        for s = 1:length(all_subs)
            text((1:5)+(rand-1.5)*spread,dgain_mean1([1,(-1:2)+4*i],1,s,1),char(64+s),'FontSize',7,'Color','k')
        end
        % end
        hold off
        ylabel(ha(i,1),ylabels{i})
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
    % leg1 = legend(ha(2,end),h1,{'Med±95%CI Change','No Change','MCID','Horizontal','Anterior','Posterior'},'Location','northwest','NumColumns',3,'Orientation','horizontal');
    leg1 = legend(ha(2,end),h1,{'Med±95%CI Change','No Change','MCID'},'Location','northwest','NumColumns',3,'Orientation','horizontal');
    leg1.ItemTokenSize(1) = 10;   
    fname2 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryGNOvHITChangeGain_AllSub_SameAxisCanals_AvgCanals.fig'];
    savefig(fig2,fname2)
    saveas(fig2,strrep(fname2,'.fig','.svg'))
    
    %% Implanted Canals (all canals on same axis) Change from Preop Gain Over Time with all Subjects
    % Take median and CI of all canals and subjects (no mean across canals
    % first)
    
    colors = [0.75,0,0;0,0.25,0;0,0,0.5];
    titles = {'\Delta vHIT Gain'};
    ylabels = {'Treatment Mode','Tonic Mode'};
    fig2 = figure;
    set(fig2,'Units','inches','Position',[1 1 6 4],'Color',[1 1 1])
    %Set figure and axes sizing
    nr = length(ylabels);
    nc = 1;
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
    h1 = gobjects(6,1);
    for i = 1:nr
        ha(i,1) = subplot(nr,nc,sub2ind([nr,nc],i,1));
        h1(1) = plot(NaN,NaN,'k-','LineWidth',1.5); % mean line for legend
        hold on
        h1(3) = fill([0,6,6,0],0.06*[-1,-1,1,1],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]); % fill MCID area
        h1(2) = plot([0,6],[0,0],'k:','LineWidth',2); % no change line

        h1(4) = plot(NaN,NaN,'-','Color',colors(1,:),'LineWidth',2); % canal colors for legend
        h1(5) = plot(NaN,NaN,'-','Color',colors(2,:),'LineWidth',2);
        h1(6) = plot(NaN,NaN,'-','Color',colors(3,:),'LineWidth',2);
    
        errorbar(1:4,dgain_med3([1,(-1:1)+4*i],1,1,1),dgain_lowCIwid3([1,(-1:1)+4*i],1,1,1),dgain_highCIwid3([1,(-1:1)+4*i],1,1,1),'k-','LineWidth',1.5);
        errorbar(5,dgain_med3(2+4*i,1,1,1),dgain_lowCIwid3(2+4*i,1,1,1),dgain_highCIwid3(2+4*i,1,1,1),'k','LineWidth',1.5);
        for j = 1:length(canals)
            for s = 1:length(all_subs)
                text((1:5)+(rand-1.5)*spread,dgain([1,(-1:2)+4*i],j,s,1),char(64+s),'FontSize',7,'Color',colors(j,:))
            end
        end
        hold off
        ylabel(ha(i,1),ylabels{i})
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
    leg1 = legend(ha(2,end),h1,{'Med±95%CI Change','No Change','MCID','Horizontal','Anterior','Posterior'},'Location','northwest','NumColumns',3,'Orientation','horizontal');
    leg1.ItemTokenSize(1) = 10;   
    fname2 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryGNOvHITChangeGain_AllSub_SameAxisCanals_NoAvgCanals.fig'];
    savefig(fig2,fname2)
    saveas(fig2,strrep(fname2,'.fig','.svg'))
    %% Implanted Canals Change from Preop Latency Over Time with all Subjects
    titles = {'Treatment Mode','Tonic Mode'};
    fig3 = figure;
    set(fig3,'Units','inches','Position',[1 1 6 4],'Color',[1 1 1])
    %Set figure and axes sizing
    nr = length(canals);
    nc = length(titles);
    x_spac = 0.01;
    x_min = 0.10;
    x_max = 1-x_spac;
    y_min = 0.17;
    y_max = 0.95;
    y_spac = 0.02;
    
    % this sizes each panel so that the scale is the same but the y limits
    % are bigger or smaller depending on your specifications
    y_lims = [-150 100; -60 120; -60 60]; % manually put in y limits - from top row to bottom row
    y_weights = (y_lims(:,2) - y_lims(:,1))/sum(y_lims(:,2) - y_lims(:,1));

    y_wid = (y_max-y_min-y_spac*(nr-1)).*y_weights;
    x_wid = (x_max-x_min-x_spac*(nc-1))/nc;
    x_pos = x_min:(x_wid+x_spac):x_max;

    y_pos = y_max - cumsum(y_wid + y_spac*[0; ones(length(y_wid)-1,1)]);
    % y_pos = fliplr(y_min:(y_wid+y_spac):y_max);
    spread = 0.1;
    ha = gobjects(nr,nc);
    h1 = gobjects(2,1);
    for i = 1:nr
        for j = 1:nc
            ha(i,j) = subplot(nr,nc,sub2ind([nc,nr],j,i));
            h1(1) = plot(NaN,NaN,'k-','LineWidth',1.5);
            hold on
            h1(2) = plot([0,6],[0,0],'k:','LineWidth',2);
            % errorbar(1:4,dlat_mean([1,(-1:1)+4*j],i,1,1),dlat_sem([1,(-1:1)+4*j],i,1,1),'k-','LineWidth',1.5);
            % errorbar(5,dlat_mean(2+4*j,i,1,1),dlat_sem(2+4*j,i,1,1),'k','LineWidth',1.5);
            errorbar(1:4,dlat_med1([1,(-1:1)+4*j],i,1,1),dlat_lowCIwid1([1,(-1:1)+4*j],i,1,1),dlat_highCIwid1([1,(-1:1)+4*j],i,1,1),'k-','LineWidth',1.5);
            errorbar(5,dlat_med1(2+4*j,i,1,1),dlat_lowCIwid1(2+4*j,i,1,1),dlat_highCIwid1(2+4*j,i,1,1),'k','LineWidth',1.5);
            for s = 1:length(all_subs)
                text((1:5)+(rand-1.5)*spread,dlat([1,(-1:2)+4*j],i,s,1),char(64+s),'FontSize',7)
            end
            hold off
        end
        ylabel(ha(i,1),{'\Delta vHIT Latency [ms]',canals{i}})
    end
    for j = 1:nc
        for i = 1:nr
            ha(i,j).Position = [x_pos(j),y_pos(i),x_wid,y_wid(i)];
            ylim(ha(i,j),y_lims(i,:))
        end
        title(ha(1,j),titles{j})
    end
    % set axis properties (including inverting y axis so up (lower/negative
    % latency) is better
    set(ha,'xminortick','off','XTick',1:5,'XTickLabels',{'0','0.5','1','2','most recent'},...
        'box','on','XLim',[0.75 5.25],'Layer','top', 'YDir', 'reverse')
    set(ha(:,2:end),'YTickLabel',[])
    set(ha(1:end-1,:),'XTickLabel',[])
    xlabel(ha(end,:),'Years Since Implantation')
    leg1 = legend(ha(1,1),h1,{'Med±95%CI Change','No Change'},'Location','south','NumColumns',length(h1));
    leg1.ItemTokenSize(1) = 10;   
    fname3 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryGNOvHITChangeLatency_AllSub_SepAxisCanals.fig'];
    savefig(fig3,fname3)
    saveas(fig3,strrep(fname3,'.fig','.svg'))

    %% Implanted Canals (all canals on same axis) Change from Preop Latency Over Time with all Subjects
    % Plotting mean across canals and then take med/CI across subjects
    colors = [0.75,0,0;0,0.25,0;0,0,0.5];
    titles = {'Mean \Delta vHIT Latency, ms'};
    ylabels = {'Treatment Mode','Tonic Mode'};
    fig4 = figure;
    set(fig4,'Units','inches','Position',[1 1 6 4],'Color',[1 1 1])
    %Set figure and axes sizing
    nr = length(ylabels);
    nc = 1;
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
        ha(i,1) = subplot(nr,nc,sub2ind([nr,nc],i,1));
        h1(1) = plot(NaN,NaN,'k-','LineWidth',1.5); % mean line for legend
        hold on
        h1(3) = fill([0,6,6,0],0.06*[-1,-1,1,1],[1,1,1],'EdgeColor',[1,1,1]); % invisible fill box for legend layout
        h1(2) = plot([0,6],[0,0],'k:','LineWidth',2); % no change line

        % h1(4) = plot(NaN,NaN,'-','Color',colors(1,:),'LineWidth',2); % canal colors for legend
        % h1(5) = plot(NaN,NaN,'-','Color',colors(2,:),'LineWidth',2);
        % h1(6) = plot(NaN,NaN,'-','Color',colors(3,:),'LineWidth',2);
    
        % errorbar(1:4,dlat_mean2([1,(-1:1)+4*i],1,1,1),dlat_sem2([1,(-1:1)+4*i],1,1,1),'k-','LineWidth',1.5);
        % errorbar(5,dlat_mean2(2+4*i,1,1,1),dlat_sem2(2+4*i,1,1,1),'k','LineWidth',1.5);
        errorbar(1:4,dlat_med2([1,(-1:1)+4*i],1,1,1),dlat_lowCIwid2([1,(-1:1)+4*i],1,1,1),dlat_highCIwid2([1,(-1:1)+4*i],1,1,1),'k-','LineWidth',1.5);
        errorbar(5,dlat_med2(2+4*i,1,1,1),dlat_lowCIwid2(2+4*i,1,1,1),dlat_highCIwid2(2+4*i,1,1,1),'k','LineWidth',1.5);
        % for j = 1:length(canals)
        for s = 1:length(all_subs)
            text((1:5)+(rand-1.5)*spread,dlat_mean1([1,(-1:2)+4*i],1,s,1),char(64+s),'FontSize',7,'Color', 'k')
        end
        % end
        hold off
        ylabel(ha(i,1),ylabels{i})
    end
    for j = 1:nc
        for i = 1:nr
            ha(i,j).Position = [x_pos(j),y_pos(i),x_wid,y_wid];
        end
        title(ha(1,j),titles{j})
    end
    set(ha,'xminortick','off','XTick',1:5,'XTickLabels',{'0','0.5','1','2','most recent'},...
        'box','on','XLim',[0.75 5.25],'YLim',[-100 50],'Layer','top', 'YDir', 'reverse')
    set(ha(:,2:end),'YTickLabel',[])
    set(ha(1:end-1,:),'XTickLabel',[])
    xlabel(ha(end,:),'Years Since Implantation')
    % leg1 = legend(ha(2,end),h1,{'Mean±SEM Change','No Change','','Horizontal','Anterior','Posterior'},'Location','northwest','NumColumns',3,'Orientation','horizontal');
    leg1 = legend(ha(2,end),h1,{'Med±95%CI Change','No Change',''},'Location','northwest','NumColumns',3,'Orientation','horizontal');
    leg1.ItemTokenSize(1) = 10;   
    fname4 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryGNOvHITChangeLatency_AllSub_SameAxisCanals_AvgCanals.fig'];
    savefig(fig4,fname4)
    saveas(fig4,strrep(fname4,'.fig','.svg'))
    
    %% Implanted Canals (all canals on same axis) Change from Preop Latency Over Time with all Subjects
    % Med/CI across all subjects and canals, plot each canal individually
    colors = [0.75,0,0;0,0.25,0;0,0,0.5];
    titles = {'\Delta vHIT Latency, ms'};
    ylabels = {'Treatment Mode','Tonic Mode'};
    fig4 = figure;
    set(fig4,'Units','inches','Position',[1 1 6 4],'Color',[1 1 1])
    %Set figure and axes sizing
    nr = length(ylabels);
    nc = 1;
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
    h1 = gobjects(6,1);
    for i = 1:nr
        ha(i,1) = subplot(nr,nc,sub2ind([nr,nc],i,1));
        h1(1) = plot(NaN,NaN,'k-','LineWidth',1.5); % mean line for legend
        hold on
        h1(3) = fill([0,6,6,0],0.06*[-1,-1,1,1],[1,1,1],'EdgeColor',[1,1,1]); % invisible fill box for legend layout
        h1(2) = plot([0,6],[0,0],'k:','LineWidth',2); % no change line

        h1(4) = plot(NaN,NaN,'-','Color',colors(1,:),'LineWidth',2); % canal colors for legend
        h1(5) = plot(NaN,NaN,'-','Color',colors(2,:),'LineWidth',2);
        h1(6) = plot(NaN,NaN,'-','Color',colors(3,:),'LineWidth',2);
    
        % errorbar(1:4,dlat_mean3([1,(-1:1)+4*i],1,1,1),dlat_sem2([1,(-1:1)+4*i],1,1,1),'k-','LineWidth',1.5);
        % errorbar(5,dlat_mean3(2+4*i,1,1,1),dlat_sem2(2+4*i,1,1,1),'k','LineWidth',1.5);
        errorbar(1:4,dlat_med3([1,(-1:1)+4*i],1,1,1),dlat_lowCIwid3([1,(-1:1)+4*i],1,1,1),dlat_highCIwid3([1,(-1:1)+4*i],1,1,1),'k-','LineWidth',1.5);
        errorbar(5,dlat_med3(2+4*i,1,1,1),dlat_lowCIwid3(2+4*i,1,1,1),dlat_highCIwid3(2+4*i,1,1,1),'k','LineWidth',1.5);
        for j = 1:length(canals)
            for s = 1:length(all_subs)
                text((1:5)+(rand-1.5)*spread,dlat([1,(-1:2)+4*i],j,s,1),char(64+s),'FontSize',7,'Color',colors(j,:))
            end
        end
        hold off
        ylabel(ha(i,1),ylabels{i})
    end
    for j = 1:nc
        for i = 1:nr
            ha(i,j).Position = [x_pos(j),y_pos(i),x_wid,y_wid];
        end
        title(ha(1,j),titles{j})
    end
    set(ha,'xminortick','off','XTick',1:5,'XTickLabels',{'0','0.5','1','2','most recent'},...
        'box','on','XLim',[0.75 5.25],'YLim',[-150 150],'Layer','top', 'YDir', 'reverse')
    set(ha(:,2:end),'YTickLabel',[])
    set(ha(1:end-1,:),'XTickLabel',[])
    xlabel(ha(end,:),'Years Since Implantation')
    leg1 = legend(ha(1,end),h1,{'Med±95%CI Change','No Change','','Horizontal','Anterior','Posterior'},'Location','southwest','NumColumns',3,'Orientation','horizontal');
    leg1.ItemTokenSize(1) = 10;   
    fname4 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryGNOvHITChangeLatency_AllSub_SameAxisCanals_NoAvgCanals.fig'];
    savefig(fig4,fname4)
    saveas(fig4,strrep(fname4,'.fig','.svg'))

    %% Gain in all canals for both sides all conditions
    n_subs = length(all_subs);
    %Plot defaults
    sub_mark = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'; % Subject markers
    line_bold = 2;
    line_norm = 1;
    line_thin = 0.5;
    mark_size = 9;
    lab_font = 10;
    %Plot sizes
    nr = 2;
    nc = 1;
    x_plot = [0.07,0.99,0.01]; %min, max, spac
    y_plot = [0.11,0.98,0.01]; %min, max, spac
    xwid = (x_plot(2)-x_plot(1)-(nc-1)*x_plot(3))/nc;
    xpos = x_plot(1):(xwid+x_plot(3)):x_plot(2);
    ywid = (y_plot(2)-y_plot(1)-(nr-1)*y_plot(3))/nr;
    ypos = fliplr(y_plot(1):(ywid+y_plot(3)):y_plot(2));
    %Plot labels
    YLim = [-0.25 1.1];
    ylabs = {'Implanted Ear Gain';'Contralateral Ear Gain'};
    XLim = [0.5 6.5];
    XTick = 1:6;
    XTickLabel = {'','','','','',''};
    colors = [0.75,0,0;0,0.25,0;0,0,0.5];
    % fig = figure(1);
    fig5 = figure;
    clf;
    set(fig5,'Color','w','Units','inches','Position',[1,1,8,7])
    ha = gobjects(nr,nc);
    h1 = gobjects(3,1);
    hs = gobjects(n_subs,1);
    offs = reshape(linspace(-0.25,0.25,n_subs*5),5,n_subs);
    
    % Set up Boxplot Data
    boxplotData = zeros(3*n_subs,6,2);
    mapData2Plot = [1 2 7 3 8 4];
    for i = 1:6
        for j = 1:nr
            temp = gain(mapData2Plot(i),:,:,j);
            boxplotData(:,i,j) = temp(:);
        end
    end
    
    for j = 1:nc
        for i = 1:nr
            ha(i,j) = subplot(nr,nc,sub2ind([nc,nr],j,i));
            plot(NaN,NaN)        
            hold on
            % Plot individual data points
            for s = 1:n_subs   
                for c = 1:length(canals)
                    plot((1:6)+offs(c,s),reshape(gain(mapData2Plot, c, s, i), [1, 6]),':',...
                        'Color',colors(c,:),'LineWidth',line_thin)
                    text((1:6)+offs(c,s),reshape(gain(mapData2Plot, c, s, i), [1, 6]),sub_mark(s),...
                        'Color',colors(c,:),'FontSize',mark_size,'HorizontalAlignment','center')
                end
            end
            
            % Make Boxplots
            boxplot(boxplotData(:,:,i),'Notch','on','Colors','k','Symbol','')
    
        end
        
    end
    for i = 1:nr
        for j = 1:nc
            set(ha(i,j),'Position',[xpos(j) ypos(i) xwid ywid])
        end
    %     ylabel(ha(i,1),ylabs{i},'FontSize',lab_font)
    end
    hold on
    h1(1) = plot(NaN,NaN,'-','Color',colors(1,:),'LineWidth',line_bold);
    h1(2) = plot(NaN,NaN,'-','Color',colors(2,:),'LineWidth',line_bold);
    h1(3) = plot(NaN,NaN,'-','Color',colors(3,:),'LineWidth',line_bold);
    
    hold off
    leg = legend(ha(end,end),h1,{'Horizontal','Anterior','Posterior'},'NumColumns',length(h1),'box','off','Location','north');
    leg.ItemTokenSize(1) = 7;
    annotation('textbox','Position',[0.23,0.915,0.14,0.047],'String',{'Subject:'},...
        'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8,'EdgeColor','none')
    % leg.Position = [0.99-0.33,0.005,0.33,0.047];
    
    % Create subject marker legend
    leg_cell = {'',''};
    for i = 1:n_subs
        if i < n_subs/2 + 1
            leg_cell{1} = [leg_cell{1},sub_mark(i),': ',num2str(i),'  '];
        else
            leg_cell{2} = [leg_cell{2},sub_mark(i),': ',num2str(i)];
            if i ~= n_subs
                leg_cell{2} = [leg_cell{2},'  '];
            end
        end
    end
    
    leg2 = annotation(fig5,'textbox',[0.288, 0.915, 0.34, 0.05],'String',leg_cell,'FontSize',8,'LineStyle','none');
    
    annotation('textbox','Position',[0.345,0.484,0.11,0.047],'String',{'Canal:'},...
        'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8,'EdgeColor','none')
    set(ha,'XLim',XLim,'ygrid','off','xgrid','off','XMinorGrid','off',...
        'XMinorTick','off','XTick',XTick,'XTickLabel',XTickLabel,...
        'XTickLabelRotation',0,'YLim',YLim)
    set(ha(1:end-1,:),'XTickLabel',[]);
    set(ha(:,2:end),'YTickLabel',[]);
    
    % save as .svg and .fig files
    fname5 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryGNOvHITGain_AllSub.fig'];
    savefig(fig5,fname5)
    saveas(fig5,strrep(fname5,'.fig','.svg'))

    %% Latency plot
    YLim = [0 150];
    fig6 = figure;
    clf;
    set(fig6,'Color','w','Units','inches','Position',[1,1,8,7])
    ha = gobjects(nr,nc);
    h1 = gobjects(3,1);
    hs = gobjects(n_subs,1);
    % Set up Boxplot Data
    boxplotData = zeros(3*n_subs,6,2);
    mapData2Plot = [1 2 7 3 8 4];
    for i = 1:6
        for j = 1:nr
            temp = lat(mapData2Plot(i),:,:,j);
            boxplotData(:,i,j) = temp(:);
        end
    end
    
    for j = 1:nc
        for i = 1:nr
            ha(i,j) = subplot(nr,nc,sub2ind([nc,nr],j,i));
            plot(NaN,NaN)        
            hold on
            % Plot individual data points
            for s = 1:n_subs   
                for c = 1:length(canals)
                    plot((1:6)+offs(c,s),reshape(lat(mapData2Plot, c, s, i), [1, 6]),':',...
                        'Color',colors(c,:),'LineWidth',line_thin)
                    text((1:6)+offs(c,s),reshape(lat(mapData2Plot, c, s, i), [1, 6]),sub_mark(s),...
                        'Color',colors(c,:),'FontSize',mark_size,'HorizontalAlignment','center')
                end
            end
            
            % Make Boxplots
            boxplot(boxplotData(:,:,i),'Notch','on','Colors','k','Symbol','')
    
        end
        
    end
    for i = 1:nr
        for j = 1:nc
            set(ha(i,j),'Position',[xpos(j) ypos(i) xwid ywid])
        end
    %     ylabel(ha(i,1),ylabs{i},'FontSize',lab_font)
    end
    hold on
    h1(1) = plot(NaN,NaN,'-','Color',colors(1,:),'LineWidth',line_bold);
    h1(2) = plot(NaN,NaN,'-','Color',colors(2,:),'LineWidth',line_bold);
    h1(3) = plot(NaN,NaN,'-','Color',colors(3,:),'LineWidth',line_bold);
    
    hold off
    leg = legend(ha(end,end),h1,{'Horizontal','Anterior','Posterior'},'NumColumns',length(h1),'box','off','Location','south');
    leg.ItemTokenSize(1) = 7;
    annotation('textbox','Position',[0.28,0.56,0.14,0.047],'String',{'Subject:'},...
        'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8,'EdgeColor','none')
    % leg.Position = [0.99-0.33,0.005,0.33,0.047];
    
    % Create subject marker legend
    leg_cell = {'',''};
    for i = 1:n_subs
        if i < n_subs/2 + 1
            leg_cell{1} = [leg_cell{1},sub_mark(i),': ',num2str(i),'  '];
        else
            leg_cell{2} = [leg_cell{2},sub_mark(i),': ',num2str(i)];
            if i ~= n_subs
                leg_cell{2} = [leg_cell{2},'  '];
            end
        end
    end
    
    leg2 = annotation(fig6,'textbox',[0.335,0.56,0.34,0.05],'String',leg_cell,'FontSize',8,'LineStyle','none');
    
    annotation('textbox','Position',[0.345,0.484,0.11,0.047],'String',{'Canal:'},...
        'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8,'EdgeColor','none')
    set(ha,'XLim',XLim,'ygrid','off','xgrid','off','XMinorGrid','off',...
        'XMinorTick','off','XTick',XTick,'XTickLabel',XTickLabel,...
        'XTickLabelRotation',0,'YLim',YLim, 'YDir', 'reverse')
    set(ha(1:end-1,:),'XTickLabel',[]);
    set(ha(:,2:end),'YTickLabel',[]);
    
    % save as .svg and .fig files
    fname6 = ['Summary Figures',filesep,char(datetime('now','Format','yyyyMMdd')),'_SummaryGNOvHITLatency_AllSub.fig'];
    savefig(fig6,fname6)
    saveas(fig6,strrep(fname6,'.fig','.svg'))

end
end