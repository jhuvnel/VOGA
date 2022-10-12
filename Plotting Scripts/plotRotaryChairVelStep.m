%% Plot the Gain Tau Product Over Time for Each Subject
%Find the directory for each MVI Subject
allMVIpath = uigetdir('','Select the path with the study subject folders.');
sub_paths = strcat(allMVIpath,filesep,extractfield(dir([allMVIpath,filesep,'MVI*R*']),'name',extractfield(dir([allMVIpath,filesep,'MVI*R*']),'isdir')));
savefigpath = uigetdir('','Select where to save the figures.');
%Load items
load('VNELcolors.mat','colors')
%Load the subject info file
warning('off')
sub_info = readtable('SubjectInfo.xlsx');
warning('on')
Subs = sub_info{:,1};
Ears = sub_info{:,2};
for s = 1:length(sub_paths)
    %Load, find rotary chair exponentials and error handle
    tab_name = extractfield(dir([sub_paths{s},filesep,'*Results.mat']),'name');
    if isempty(tab_name)
        error(['No subject-wide results table exists yet in ',cd])
    end
    load([sub_paths{s},filesep,tab_name{end}],'all_results')
    load(strrep([sub_paths{s},filesep,tab_name{end}],'Results','CycParam'),'cyc_params')
    %For the purpose of this figure, only look at +/- 240dps
    [~,sub_num] = ismember(all_results.Subject,Subs);
    wanted_amp = 240*ones(length(all_results.Subject),1); %Assume left ear and correct
    wanted_amp(strcmp(Ears(sub_num),'R')) = -240;
    stim_ear_bool = all_results.('Amplitude(dps)')==wanted_amp&contains(all_results.Experiment,'RotaryChair')&contains(all_results.Type,'Exp');
    contra_ear_bool = all_results.('Amplitude(dps)')==-wanted_amp&contains(all_results.Experiment,'RotaryChair')&contains(all_results.Type,'Exp');
    if ~any(stim_ear_bool)
        error(['No rotary chair velocity step experiments detected in the results table ',tab_name{end}])
    end
    stim_ear_params = cyc_params(stim_ear_bool,:);
    contra_ear_params = cyc_params(contra_ear_bool,:);
    contra_ear_tab = all_results(contra_ear_bool,[1:5,8,13]);
    stim_ear_tab = all_results(stim_ear_bool,[1:5,8,13]);
    stim_ear_tab.K1 = NaN(size(stim_ear_tab,1),1);
    stim_ear_tab.Tau1 = NaN(size(stim_ear_tab,1),1);
    stim_ear_tab.K2 = NaN(size(stim_ear_tab,1),1);
    stim_ear_tab.Tau2 = NaN(size(stim_ear_tab,1),1);
    contra_ear_tab.K1 = NaN(size(contra_ear_tab,1),1);
    contra_ear_tab.Tau1 = NaN(size(contra_ear_tab,1),1);
    contra_ear_tab.K2 = NaN(size(contra_ear_tab,1),1);
    contra_ear_tab.Tau2 = NaN(size(contra_ear_tab,1),1);
    for i = 1:size(stim_ear_tab,1)
        lz_vals = stim_ear_params{i,2}.exp_ord1_const_high_params.lz;
        rz_vals = stim_ear_params{i,2}.exp_ord1_const_high_params.rz;
        gains = -[lz_vals(1)+lz_vals(3), rz_vals(1)+rz_vals(3)]./stim_ear_tab.('Amplitude(dps)')(i);
        taus = [lz_vals(2), rz_vals(2)];
        [K,eye] = max(gains); %use the bigger gain and associated time constant
        stim_ear_tab.K1(i) = K;
        stim_ear_tab.Tau1(i) = taus(eye);
        lz_vals = stim_ear_params{i,2}.exp_ord1_high_params.lz;
        rz_vals = stim_ear_params{i,2}.exp_ord1_high_params.rz;
        gains = -[lz_vals(1), rz_vals(1)]./stim_ear_tab.('Amplitude(dps)')(i);
        taus = [lz_vals(2), rz_vals(2)];
        [K,eye] = max(gains); %use the bigger gain and associated time constant
        stim_ear_tab.K2(i) = K;
        stim_ear_tab.Tau2(i) = taus(eye);
    end
    stim_ear_tab.GTP1 = stim_ear_tab.K1.*stim_ear_tab.Tau1;
    stim_ear_tab.GTP2 = stim_ear_tab.K2.*stim_ear_tab.Tau2;
    for i = 1:size(contra_ear_tab,1)
        lz_vals = contra_ear_params{i,2}.exp_ord1_const_high_params.lz;
        rz_vals = contra_ear_params{i,2}.exp_ord1_const_high_params.rz;
        gains = -[lz_vals(1)+lz_vals(3), rz_vals(1)+rz_vals(3)]./contra_ear_tab.('Amplitude(dps)')(i);
        taus = [lz_vals(2), rz_vals(2)];
        [K,eye] = max(gains); %use the bigger gain and associated time constant
        contra_ear_tab.K1(i) = K;
        contra_ear_tab.Tau1(i) = taus(eye);
        lz_vals = contra_ear_params{i,2}.exp_ord1_high_params.lz;
        rz_vals = contra_ear_params{i,2}.exp_ord1_high_params.rz;
        gains = -[lz_vals(1), rz_vals(1)]./contra_ear_tab.('Amplitude(dps)')(i);
        taus = [lz_vals(2), rz_vals(2)];
        [K,eye] = max(gains); %use the bigger gain and associated time constant
        contra_ear_tab.K2(i) = K;
        contra_ear_tab.Tau2(i) = taus(eye);
    end
    contra_ear_tab.GTP1 = contra_ear_tab.K1.*contra_ear_tab.Tau1;
    contra_ear_tab.GTP2 = contra_ear_tab.K2.*contra_ear_tab.Tau2;
    %Partition by condition
    no_stim_bool = contains(stim_ear_tab.Condition,'NoStim');
    mot_mod_bool = contains(stim_ear_tab.Condition,'MotionMod');
    const_rate_bool = contains(stim_ear_tab.Condition,'ConstantRate');
    t_yr = days(stim_ear_tab.Date - sub_info{sub_num(1),3})/365.25; %in years
    %% Plot over time
    fig_name = [stim_ear_tab.Subject{1},' Velocity Step Gain Tau Product'];
    fig1 = figure(1);
    set(fig1,'Color',[1,1,1])
    plot(NaN,NaN)
    hold on
    h1 = plot(t_yr(no_stim_bool),stim_ear_tab.GTP2(no_stim_bool),'ko');
    plot(t_yr(no_stim_bool),stim_ear_tab.GTP2(no_stim_bool),'k--')
    h2 = plot(t_yr(mot_mod_bool),stim_ear_tab.GTP2(mot_mod_bool),'go');
    plot(t_yr(mot_mod_bool),stim_ear_tab.GTP2(mot_mod_bool),'g--')
    h3 = plot(t_yr(const_rate_bool),stim_ear_tab.GTP2(const_rate_bool),'ro');
    plot(t_yr(const_rate_bool),stim_ear_tab.GTP2(const_rate_bool),'r--')
    hold off
    xlabel('Time (yr)')
    ylabel('Gain Tau Product (s)')
    title(fig_name)
    leg1 = legend([h1(1),h2(1),h3(1)],{'Off','On','Tonic'},'location','northwest','NumColumns',3);
    leg1.ItemTokenSize(1) = 6;
    title(leg1,'Condition')
    set(gca,'XLim',[-0.5 6.5],'YLim',[-0.5 12])
    savefig(fig1,[savefigpath,filesep,date,'_',strrep(fig_name,' ',''),'.fig'])
    saveas(fig1,[savefigpath,filesep,date,'_',strrep(fig_name,' ',''),'.png'])
    %% Plot all time points as dots (K vs. tau) for stim ear
    % fig2 = figure(2);
    % set(fig2,'Color',[1,1,1])
    % plot(NaN,NaN)
    % hold on
    % h1 = plot(stim_ear_tab.Tau2(no_stim_bool),stim_ear_tab.K2(no_stim_bool),'ko');
    % h2 = plot(stim_ear_tab.Tau2(mot_mod_bool),stim_ear_tab.K2(mot_mod_bool),'go');
    % h3 = plot(stim_ear_tab.Tau2(const_rate_bool),stim_ear_tab.K2(const_rate_bool),'ro');
    % hold off
    % xlabel('Tau (s)')
    % ylabel('Gain')
    % title([stim_ear_tab.Subject(1),' Velocity Step Results'])
    % leg1 = legend([h1(1),h2(1),h3(1)],{'Off','On','Tonic'},'location','southeast');
    % title(leg1,'Condition')
end