function plotRotaryChairDSMB
%% Make Rotary Chair Results mat
if ispc %AIA lab machine
    Path = 'Z:\Study Subjects';
else %AIA Mac
    Path = '/Volumes/vnelhuman$/MVI/Study Subjects';
end
sub_folds = extractfield(dir(Path),'name',extractfield(dir(Path),'isdir')&contains(extractfield(dir(Path),'name'),'MVI'));
tabs = cell(length(sub_folds),1);
for i = 1:length(sub_folds)
    rel_fold = extractfield(dir([Path,filesep,sub_folds{i},filesep,'*Results.mat']),'name');
    if ~isempty(rel_fold)
        %Load most recent version
        fname = [Path,filesep,sub_folds{i},filesep,rel_fold{end}]; %get the most recent item
        load(fname,'all_results')
        %Run again
        tabs{i} = all_results;
    end
end
if ~isempty(vertcat(tabs{:}))
    all_results_allexp = vertcat(tabs{:});
    delete([Path,filesep,'*Results.mat']) %Remove outdated versions
    exp_types = unique(all_results_allexp.Experiment,'stable');
    for i = 1:length(exp_types)
        all_results = all_results_allexp(contains(all_results_allexp.Experiment,exp_types{i}),:);
        save([Path,filesep,datestr(now,'yyyymmdd_HHMMSS'),'_AllSubjects_',strrep(exp_types{i},' ',''),'Results.mat'],'all_results')
    end
end
all_results_MVI = all_results_allexp(contains(all_results_allexp.Experiment,'Rotary')&contains(all_results_allexp.Type,'Sine'),:);
[filename,pathname] = uigetfile('RotaryChairResults.mat', 'To plot a candidate on the summary figure, select their Results file.');
if all(pathname==0)
    all_results = [];
else
    load([pathname,filesep,filename],'all_results')
    all_results = all_results(contains(all_results.Type,'Sine')&contains(all_results.Condition,'NoStim'),:);
end
load('RotaryChairNormativeData.mat','norm_dat')
%Make sure it's sorted by Subject first and then by date
all_results_MVI = sortrows(sortrows(all_results_MVI,'Date','ascend'),'Subject','ascend');
%% Extract relevant values
% Limit to 0.05-1Hz (0.05, 0.1, 0.2, 0.5, 1)
remove_phase_with_gain_below = 0.025;
freqs = [0.05,0.1,0.2,0.5,1];
all_subs = unique(all_results_MVI.Subject);
%Make structs for each condition
conds = {'pre','mmo'};
cond_logic = [contains(all_results_MVI.Visit,'Visit0')&contains(all_results_MVI.Condition,'NoStim'),...
    contains(all_results_MVI.Condition,'Motion')];
%Initialize arrays
for c = 1:length(conds)
    gain.(conds{c}) = NaN(length(all_subs),length(freqs));
    gain_sd.(conds{c}) = NaN(length(all_subs),length(freqs));
    phase.(conds{c}) = NaN(length(all_subs),length(freqs));
    phase_sd.(conds{c}) = NaN(length(all_subs),length(freqs));
end
%Populate arrays
for i = 1:length(all_subs) 
    if contains(all_subs(i),{'MVI001','MVI002',',MVI003','MVI004','MVI007','MVI009'}) %L ear
        tab_field1 = 'Gain_LZ_HIGH';
        tab_field2 = 'Gain_RZ_HIGH';        
    else %R ear
        tab_field1 = 'Gain_LZ_LOW';
        tab_field2 = 'Gain_RZ_LOW';  
    end
    for j = 1:length(freqs)
        for c = 1:length(conds)
            subtab = all_results_MVI(contains(all_results_MVI.Subject,all_subs(i))&all_results_MVI.('Frequency(Hz)')==freqs(j)&cond_logic(:,c),:);
            if ~isempty(subtab)
                subtab = subtab(end,:);
                if subtab.(tab_field1)>subtab.(tab_field2) %Use L eye
                    gain.(conds{c})(i,j) = subtab.(tab_field1);
                    gain_sd.(conds{c})(i,j) = subtab.([tab_field1,'_sd']);
                    phase.(conds{c})(i,j) = subtab.Phase_L;
                    phase_sd.(conds{c})(i,j) = subtab.Phase_L_sd;
                else
                    gain.(conds{c})(i,j) = subtab.(tab_field2);
                    gain_sd.(conds{c})(i,j) = subtab.([tab_field2,'_sd']);
                    phase.(conds{c})(i,j) = subtab.Phase_R;
                    phase_sd.(conds{c})(i,j) = subtab.Phase_R_sd;
                end
            end
        end
    end
end
% Make summary lines
for c = 1:length(conds)
    phase.(conds{c})(gain.(conds{c})<remove_phase_with_gain_below) = NaN;
    gain_mean.(conds{c}) = mean(gain.(conds{c}),1,'omitnan');
    gain_std.(conds{c}) = std(gain.(conds{c}),[],1,'omitnan');
    phase_mean.(conds{c}) = mean(phase.(conds{c}),1,'omitnan');
    phase_std.(conds{c}) = std(phase.(conds{c}),[],1,'omitnan');
end
% See if there is anything to plot for a prospective subject
if ~isempty(all_results)
    cand_gain_L = NaN(1,length(freqs));
    cand_gain_R = NaN(1,length(freqs));
    cand_phase = NaN(1,length(freqs));
    for j = 1:length(freqs)
        subtab = all_results(all_results.('Frequency(Hz)')==freqs(j),:);
        if ~isempty(subtab)
            subtab = subtab(end,:);
            if subtab.Gain_LZ_HIGH>subtab.Gain_RZ_HIGH %Use L eye
                cand_gain_L(1,j) = subtab.Gain_LZ_HIGH;
            else
                cand_gain_L(1,j) = subtab.Gain_RZ_HIGH;
            end
            if subtab.Gain_LZ_LOW>subtab.Gain_RZ_LOW %Use L eye
                cand_gain_R(1,j) = subtab.Gain_LZ_LOW;
            else
                cand_gain_R(1,j) = subtab.Gain_RZ_LOW;
            end  
            if (subtab.Gain_LZ_HIGH+subtab.Gain_LZ_LOW)>(subtab.Gain_RZ_HIGH+subtab.Gain_RZ_LOW) %Use L eye
                cand_phase(1,j) = subtab.Phase_L;
            else
                cand_phase(1,j) = subtab.Phase_R;
            end
        end
    end
    cand_phase(cand_gain_L<remove_phase_with_gain_below&cand_gain_R<remove_phase_with_gain_below)= NaN; 
end
%% Plot
sub_mark = 'xdo^ps+hv<';
spread = 0.05;
mark_size = 7;
h2 = gobjects(length(all_subs),1);
figure('Units','inches','Position',[-7 1 6 4],'Color',[1 1 1])
ha(1) = subplot(2,1,1);
ha(2) = subplot(2,1,2);
ha(1).Position = [0.09,0.56,0.98-0.09,0.425];
ha(2).Position = [0.09,0.12,0.98-0.09,0.425];
axes(ha(1))
plot(NaN,NaN)
hold on
h1(2) = fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.gain-norm_dat.gain_std,fliplr(norm_dat.gain+norm_dat.gain_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
h1(1) = plot(norm_dat.freq,norm_dat.gain,'k--','LineWidth',2);
h1(3) = errorbar((1-0.5*spread)*freqs,gain_mean.pre,gain_std.pre,'k-','LineWidth',1.5);
h1(4) = errorbar((1+0.5*spread)*freqs,gain_mean.mmo,gain_std.mmo,'r-','LineWidth',1.5);
for i = 1:length(all_subs)
    plot((1-0.5*spread)*freqs,gain.pre(i,:),'k.','Marker',sub_mark(i),'MarkerSize',mark_size,'MarkerFaceColor','k')
    plot((1+0.5*spread)*freqs,gain.mmo(i,:),'r.','Marker',sub_mark(i),'MarkerSize',mark_size,'MarkerFaceColor','r')
end
if ~isempty(all_results)
    h1(5) = plot((1-0.75*spread)*freqs,cand_gain_L,'b*--','MarkerSize',mark_size);
    h1(6) = plot((1-0.75*spread)*freqs,cand_gain_R,'b*:','MarkerSize',mark_size);
    leg1_lab = {'Normal Mean','Normal±SD','Pre-Op','MVI ON','Left Ear','Right Ear'};   
else
    leg1_lab = {'Normal Mean','Normal±SD','Pre-Op','MVI ON'};
end
hold off
set(gca,'xscale','log','xminortick','off','XTick',freqs,'XTickLabel',[],'box','off')
axis([0.041 1.2 -0.025 0.79])
text(0.037,0.75,'A','FontSize',14)
ylabel('Horizontal VOR Gain')
leg1 = legend(h1,leg1_lab,'Location','northwest','NumColumns',3);
title(leg1,'Condition')
axes(ha(2))
plot(NaN,NaN)
hold on
fill([norm_dat.freq,fliplr(norm_dat.freq)],[norm_dat.phase-norm_dat.phase_std,fliplr(norm_dat.phase+norm_dat.phase_std)],0.85*[1,1,1],'EdgeColor',0.85*[1,1,1]);
plot(norm_dat.freq,norm_dat.phase,'k--','LineWidth',2);
errorbar((1-0.5*spread)*freqs,phase_mean.pre,phase_std.pre,'k-','LineWidth',1.5);
errorbar((1+0.5*spread)*freqs,phase_mean.mmo,phase_std.mmo,'r-','LineWidth',1.5);
for i = 1:length(all_subs)
    plot((1-0.5*spread)*freqs,phase.pre(i,:),'k.','Marker',sub_mark(i),'MarkerSize',mark_size,'MarkerFaceColor','k')
    plot((1+0.5*spread)*freqs,phase.mmo(i,:),'r.','Marker',sub_mark(i),'MarkerSize',mark_size,'MarkerFaceColor','r')
    h2(i) = plot(NaN,NaN,'k.','Marker',sub_mark(i),'MarkerSize',mark_size,'LineWidth',1);
end
if ~isempty(all_results)
    h2(i+1) = plot((1-0.75*spread)*freqs,cand_phase,'b*','MarkerSize',mark_size);
    leg2_lab = [split(cellstr(num2str(1:length(all_subs))));all_results.Subject(1)];
else
    leg2_lab = split(cellstr(num2str(1:length(all_subs))));
end
hold off
set(gca,'xscale','log','xminortick','off','XTick',freqs,'box','off')
axis([0.041 1.2 -25 135])
text(0.037,125,'B','FontSize',14)
xlabel('Frequency (Hz)')
ylabel('Phase Lead (deg)')
leg2 = legend(h2,leg2_lab,'Location','northwest','NumColumns',length(leg2_lab));
leg2.ItemTokenSize(1) = 8;
title(leg2,'Subjects')
end