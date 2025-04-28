close all; clear; clc
%Load tables from experiments A, B, and C (Celia knows what these mean)
files = dir([cd filesep '*' filesep '*VOGResults.mat']);

tab = [];
for i = 1:length(files)
    load([files(i).folder filesep files(i).name],'all_results');
    for j = 1:height(all_results)
        if all_results{j,'CurrentAmp'} > 9
            all_results{j,'CurrentAmp'} = 5*round(all_results{j,'CurrentAmp'}/5);
        end
    end
    tab = [tab; all_results];
end

%Sort table
tab = sortrows(sortrows(sortrows(tab,'Frequency','ascend'),'Date','ascend'),'CurrentAmp','ascend');
tab_LA = tab(2:2:end,{'File','Date','Frequency','CurrentAmp','Cycles','MaxVel_LL','MaxVel_LL_sd','MaxVel_RL','MaxVel_RL_sd','MaxVel_LZ','MaxVel_LZ_sd','MaxVel_RZ','MaxVel_RZ_sd','MaxVel_LR','MaxVel_LR_sd','MaxVel_RR','MaxVel_RR_sd','MaxVel','MaxVel_sd','Align_L','Align_L_sd','Align_R','Align_R_sd','Align','Align_sd'});
tab_RP = tab(1:2:end,{'File','Date','Frequency','CurrentAmp','Cycles','MaxVel_LL','MaxVel_LL_sd','MaxVel_RL','MaxVel_RL_sd','MaxVel_LZ','MaxVel_LZ_sd','MaxVel_RZ','MaxVel_RZ_sd','MaxVel_LR','MaxVel_LR_sd','MaxVel_RR','MaxVel_RR_sd','MaxVel','MaxVel_sd','Align_L','Align_L_sd','Align_R','Align_R_sd','Align','Align_sd'});

%% Maximum velocity with normalized y lims
figure
curr_amps = unique(tab_LA.CurrentAmp);
exps = unique(tab_LA.Date);

cols = ceil(sqrt(length(curr_amps)));
rows = ceil(length(curr_amps) / cols);

t = tiledlayout(rows, cols);
t.TileSpacing = 'compact';
t.Padding = 'tight';

for f = 1:length(curr_amps)
    nexttile
    hold on
    col = {'#b3cde3';'#8c96c6';'#8856a7';'#810f7c'};
    for j = 1:length(exps)
        ind = tab_LA.CurrentAmp==curr_amps(f)&tab_LA.Date==exps(j);
        errorbar(tab_LA.Frequency(ind),tab_LA.MaxVel(ind),tab_LA.MaxVel_sd(ind),'Color',col{j},'LineWidth',2)
    end
    legend(datestr(exps,'yyyymmdd'),'Location','northwest')
    title(['Left Anterior Stimulation at ',num2str(curr_amps(f)),'\muA'])
    xlim([0 11])
    ylim([0 max(tab_LA.MaxVel)+100])
    xticks([0.5 1 2 3 8 10])
    grid on
end
sgtitle('Maximum velocity with normalized y lims')
xlabel(t,'Frequency (Hz)')
ylabel(t,'Maximum Eye Velocity (dps)')

%% Maximum velocity
figure
curr_amps = unique(tab_LA.CurrentAmp);
exps = unique(tab_LA.Date);

cols = ceil(sqrt(length(curr_amps)));
rows = ceil(length(curr_amps) / cols);

t = tiledlayout(rows, cols);
t.TileSpacing = 'compact';
t.Padding = 'tight';

for f = 1:length(curr_amps)
    nexttile
    hold on
    col = {'#b3cde3';'#8c96c6';'#8856a7';'#810f7c'};
    for j = 1:length(exps)
        ind = tab_LA.CurrentAmp==curr_amps(f)&tab_LA.Date==exps(j);
        errorbar(tab_LA.Frequency(ind),tab_LA.MaxVel(ind),tab_LA.MaxVel_sd(ind),'Color',col{j},'LineWidth',2)
    end
    legend(datestr(exps,'yyyymmdd'),'Location','northwest')
    title(['Left Anterior Stimulation at ',num2str(curr_amps(f)),'\muA'])
    xlim([0 11])
%     ylim([0 max(tab_LA.MaxVel)+100])
    xticks([0.5 1 2 3 8 10])
    grid on
end
sgtitle('Maximum velocity')
xlabel(t,'Frequency (Hz)')
ylabel(t,'Maximum Eye Velocity (dps)')
%% Max response misalignment
figure
curr_amps = unique(tab_LA.CurrentAmp);
exps = unique(tab_LA.Date);

cols = ceil(sqrt(length(curr_amps)));
rows = ceil(length(curr_amps) / cols);

t = tiledlayout(rows, cols);
t.TileSpacing = 'compact';
t.Padding = 'tight';

for f = 1:length(curr_amps)
    nexttile
    hold on
    col = {'#b3cde3';'#8c96c6';'#8856a7';'#810f7c'};
    for j = 1:length(exps)
        ind = tab_LA.CurrentAmp==curr_amps(f)&tab_LA.Date==exps(j);
        errorbar(tab_LA.Frequency(ind),tab_LA.Align(ind),tab_LA.Align_sd(ind),'Color',col{j},'LineWidth',2)
    end
    legend(datestr(exps,'yyyymmdd'),'Location','northwest')
    title(['Left Anterior Stimulation at ',num2str(curr_amps(f)),'\muA'])
    xlim([0 11])
%     ylim([0 max(tab_LA.MaxVel)+100])
    xticks([0.5 1 2 3 8 10])
    grid on
end
sgtitle('Max response misalignment')
xlabel(t,'Frequency (Hz)')
ylabel(t,'Misalignment (deg)')
%% Maximum Left Eye LARP velocity with normalized y lims
figure
curr_amps = unique(tab_LA.CurrentAmp);
exps = unique(tab_LA.Date);

cols = ceil(sqrt(length(curr_amps)));
rows = ceil(length(curr_amps) / cols);

t = tiledlayout(rows, cols);
t.TileSpacing = 'compact';
t.Padding = 'tight';

for f = 1:length(curr_amps)
    nexttile
    hold on
    col = {'#b3cde3';'#8c96c6';'#8856a7';'#810f7c'};
    for j = 1:length(exps)
        ind = tab_LA.CurrentAmp==curr_amps(f)&tab_LA.Date==exps(j);
        errorbar(tab_LA.Frequency(ind),tab_LA.MaxVel_LL(ind),tab_LA.MaxVel_LL_sd(ind),'Color',col{j},'LineWidth',2)
    end
    legend(datestr(exps,'yyyymmdd'),'Location','northwest')
    title(['Left Anterior Stimulation at ',num2str(curr_amps(f)),'\muA'])
    xlim([0 11])
    ylim([0 max(tab_LA.MaxVel)+100])
    xticks([0.5 1 2 3 8 10])
    grid on
end
sgtitle('Maximum Left Eye LARP velocity with normalized y lims')
xlabel(t,'Frequency (Hz)')
ylabel(t,'Maximum Eye Velocity (dps)')
%% Maximum Left Eye LARP velocity
figure
curr_amps = unique(tab_LA.CurrentAmp);
exps = unique(tab_LA.Date);

cols = ceil(sqrt(length(curr_amps)));
rows = ceil(length(curr_amps) / cols);

t = tiledlayout(rows, cols);
t.TileSpacing = 'compact';
t.Padding = 'tight';

for f = 1:length(curr_amps)
    nexttile
    hold on
    col = {'#b3cde3';'#8c96c6';'#8856a7';'#810f7c'};
    for j = 1:length(exps)
        ind = tab_LA.CurrentAmp==curr_amps(f)&tab_LA.Date==exps(j);
        errorbar(tab_LA.Frequency(ind),tab_LA.MaxVel_LL(ind),tab_LA.MaxVel_LL_sd(ind),'Color',col{j},'LineWidth',2)
    end
    legend(datestr(exps,'yyyymmdd'),'Location','northwest')
    title(['Left Anterior Stimulation at ',num2str(curr_amps(f)),'\muA'])
    xlim([0 11])
%     ylim([0 max(tab_LA.MaxVel)+100])
    xticks([0.5 1 2 3 8 10])
    grid on
end
sgtitle('Maximum Left Eye LARP velocity')
xlabel(t,'Frequency (Hz)')
ylabel(t,'Maximum Eye Velocity (dps)')

%% Left Eye misalignment
figure
curr_amps = unique(tab_LA.CurrentAmp);
exps = unique(tab_LA.Date);

cols = ceil(sqrt(length(curr_amps)));
rows = ceil(length(curr_amps) / cols);

t = tiledlayout(rows, cols);
t.TileSpacing = 'compact';
t.Padding = 'tight';

for f = 1:length(curr_amps)
    nexttile
    hold on
    col = {'#b3cde3';'#8c96c6';'#8856a7';'#810f7c'};
    for j = 1:length(exps)
        ind = tab_LA.CurrentAmp==curr_amps(f)&tab_LA.Date==exps(j);
        errorbar(tab_LA.Frequency(ind),tab_LA.Align_L(ind),tab_LA.Align_L_sd(ind),'Color',col{j},'LineWidth',2)
    end
    legend(datestr(exps,'yyyymmdd'),'Location','northwest')
    title(['Left Anterior Stimulation at ',num2str(curr_amps(f)),'\muA'])
    xlim([0 11])
%     ylim([0 max(tab_LA.MaxVel)+100])
    xticks([0.5 1 2 3 8 10])
    grid on
end
sgtitle('Left Eye misalignment')
xlabel(t,'Frequency (Hz)')
ylabel(t,'Misalignment (deg)')
%% Maximum Right Eye LARP velocity with normalized y lims
figure
curr_amps = unique(tab_LA.CurrentAmp);
exps = unique(tab_LA.Date);

cols = ceil(sqrt(length(curr_amps)));
rows = ceil(length(curr_amps) / cols);

t = tiledlayout(rows, cols);
t.TileSpacing = 'compact';
t.Padding = 'tight';

for f = 1:length(curr_amps)
    nexttile
    hold on
    col = {'#b3cde3';'#8c96c6';'#8856a7';'#810f7c'};
    for j = 1:length(exps)
        ind = tab_LA.CurrentAmp==curr_amps(f)&tab_LA.Date==exps(j);
        errorbar(tab_LA.Frequency(ind),tab_LA.MaxVel_RL(ind),tab_LA.MaxVel_RL_sd(ind),'Color',col{j},'LineWidth',2)
    end
    legend(datestr(exps,'yyyymmdd'),'Location','northwest')
    title(['Left Anterior Stimulation at ',num2str(curr_amps(f)),'\muA'])
    xlim([0 11])
    ylim([0 max(tab_LA.MaxVel)+100])
    xticks([0.5 1 2 3 8 10])
    grid on
end
sgtitle('Maximum Right Eye LARP velocity with normalized y lims')
xlabel(t,'Frequency (Hz)')
ylabel(t,'Maximum Eye Velocity (dps)')
%% Maximum Right Eye LARP velocity
figure
curr_amps = unique(tab_LA.CurrentAmp);
exps = unique(tab_LA.Date);

cols = ceil(sqrt(length(curr_amps)));
rows = ceil(length(curr_amps) / cols);

t = tiledlayout(rows, cols);
t.TileSpacing = 'compact';
t.Padding = 'tight';

for f = 1:length(curr_amps)
    nexttile
    hold on
    col = {'#b3cde3';'#8c96c6';'#8856a7';'#810f7c'};
    for j = 1:length(exps)
        ind = tab_LA.CurrentAmp==curr_amps(f)&tab_LA.Date==exps(j);
        errorbar(tab_LA.Frequency(ind),tab_LA.MaxVel_RL(ind),tab_LA.MaxVel_RL_sd(ind),'Color',col{j},'LineWidth',2)
    end
    legend(datestr(exps,'yyyymmdd'),'Location','northwest')
    title(['Left Anterior Stimulation at ',num2str(curr_amps(f)),'\muA'])
    xlim([0 11])
%     ylim([0 max(tab_LA.MaxVel)+100])
    xticks([0.5 1 2 3 8 10])
    grid on
end
sgtitle('Maximum Right Eye LARP velocity')
xlabel(t,'Frequency (Hz)')
ylabel(t,'Maximum Eye Velocity (dps)')

%% Right Eye misalignment
figure
curr_amps = unique(tab_LA.CurrentAmp);
exps = unique(tab_LA.Date);

cols = ceil(sqrt(length(curr_amps)));
rows = ceil(length(curr_amps) / cols);

t = tiledlayout(rows, cols);
t.TileSpacing = 'compact';
t.Padding = 'tight';

for f = 1:length(curr_amps)
    nexttile
    hold on
    col = {'#b3cde3';'#8c96c6';'#8856a7';'#810f7c'};
    for j = 1:length(exps)
        ind = tab_LA.CurrentAmp==curr_amps(f)&tab_LA.Date==exps(j);
        errorbar(tab_LA.Frequency(ind),tab_LA.Align_R(ind),tab_LA.Align_R_sd(ind),'Color',col{j},'LineWidth',2)
    end
    legend(datestr(exps,'yyyymmdd'),'Location','northwest')
    title(['Left Anterior Stimulation at ',num2str(curr_amps(f)),'\muA'])
    xlim([0 11])
%     ylim([0 max(tab_LA.MaxVel)+100])
    xticks([0.5 1 2 3 8 10])
    grid on
end
sgtitle('Right Eye misalignment')
xlabel(t,'Frequency (Hz)')
ylabel(t,'Misalignment (deg)')

