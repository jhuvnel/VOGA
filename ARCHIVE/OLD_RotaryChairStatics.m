%% Rotary Chair Statistics.m 
% Run ART in R first
% Needs to be editted before use
%% Get information of what to plot
[path2,path1] = uigetfile('*.mat','Select Parameter Fit File');
load([path1,filesep,path2],'all_results')
subjects = unique(all_results.Subject);
sub_num = nmlistdlg('PromptString','Subject:','SelectionMode','single','ListString',subjects);
subject = subjects(sub_num);
tab = sortrows(all_results(contains(all_results.Subject,subject),:),'Date','ascend');
visits = unique(tab.Visit,'stable');
freqs = {'0.01Hz','0.02Hz','0.05Hz','0.1Hz','0.2Hz','0.5Hz','1Hz','2Hz'};
%% Ranksum to isolate by freq and condition type
p_gain = NaN(4,5);
z_gain = NaN(4,5);
p_phase = NaN(4,5);
z_phase = NaN(4,5);
for i = 1:5
    MMO_gain = s_gain_cyc{4,i};
    MMO_phase = s_phase_cyc{4,i};    
    CRO_gain = s_gain_cyc{3,i};
    CRO_phase = s_phase_cyc{3,i};   
    PreOp_gain = s_gain_cyc{1,i};
    PreOp_phase = s_phase_cyc{1,i};    
    PostOp_gain = s_gain_cyc{2,i};
    PostOp_phase = s_phase_cyc{2,i};    
    %Motion Mod vs. Preop
    if ~isempty(MMO_gain) && ~isempty(PreOp_gain)
        [p_gain(1,i),~,stat] = ranksum(MMO_gain,PreOp_gain,'method','approximate');
        z_gain(1,i) = stat.zval;
        if ~isempty(MMO_phase) && ~isempty(PreOp_phase)
            [p_phase(1,i),~,stat] = ranksum(MMO_phase,PreOp_phase,'method','approximate');
            z_phase(1,i) = stat.zval;
        end
    end
    %Motion Mod vs. ConstantRate
    if ~isempty(MMO_gain) && ~isempty(CRO_gain)
        [p_gain(2,i),~,stat] = ranksum(MMO_gain,CRO_gain,'method','approximate');
        z_gain(2,i) = stat.zval;
        if ~isempty(MMO_phase) && ~isempty(CRO_phase)
            [p_phase(2,i),~,stat] = ranksum(MMO_phase,CRO_phase,'method','approximate');
            z_phase(2,i) = stat.zval;
        end
    end
    %Post-op vs. Pre-op
    if ~isempty(PreOp_gain) && ~isempty(PostOp_gain)
        [p_gain(3,i),~,stat] = ranksum(PostOp_gain,PreOp_gain,'method','approximate');
        z_gain(3,i) = stat.zval;
        if ~isempty(PreOp_phase) && ~isempty(PostOp_phase)
            [p_phase(3,i),~,stat] = ranksum(PostOp_phase,PreOp_phase,'method','approximate');
            z_phase(3,i) = stat.zval;
        end
    end
    %Post-op vs. Motion Mod On
    if ~isempty(PostOp_gain) && ~isempty(MMO_gain)
        [p_gain(4,i),~,stat] = ranksum(PostOp_gain,MMO_gain,'method','approximate');
        z_gain(4,i) = stat.zval;
        if ~isempty(MMO_phase) && ~isempty(PostOp_phase)
            [p_phase(4,i),~,stat] = ranksum(PostOp_phase,MMO_phase,'method','approximate');
            z_phase(4,i) = stat.zval;
        end
    end
end

% Now test just along condition type
p_gain_cond = NaN(4,1);
ind1 = [4,4,2,2];
ind2 = [1,3,1,4];
for i = 1:4
    if ~all(isnan(gain_mat(ind1(i),:,1))|isnan(gain_mat(ind2(i),:,1)))
        p_gain_cond(i,1) = signrank(gain_mat(ind1(i),:,1),gain_mat(ind2(i),:,1),'method','approximate');
    end
end

stat_tab = table();
stat_tab(1:4,1) = cell2table({'Gain: Motion Mod vs. PreOp';...
                              'Gain: Motion Mod vs. Constant Rate';...
                              'Gain: Post-Op vs. Pre-Op';...
                              'Gain: Post-Op vs. Motion Mod'});
stat_tab(1:4,2:7) = array2table([p_gain_cond,p_gain]);
tab_labs = {'Test','Overall','p_0p1Hz','p_0p2Hz','p_0p5Hz','p_1Hz','p_2Hz'};
stat_tab.Properties.VariableNames = tab_labs;
disp('NO p-value corrections')
disp(stat_tab)

figure;
% Get the table in string form.
TString = evalc('disp(stat_tab)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);