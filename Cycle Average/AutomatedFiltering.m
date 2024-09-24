function [CycAvg,analyzed] = AutomatedFiltering(Data,mode)
%% Initialize
% Plot defaults
% Set colors and line widths
load('VNELcolors.mat','colors')
colors.cyc_keep = [0.85 0.85 0.85]; % Fill colors for cycle selection
colors.cyc_rm = [1 1 1];
plot_info.colors = colors;
plot_info.line_wid.norm = 0.5;
plot_info.line_wid.bold = 2;
% Set Experimentor/version
VOGA_VerInfo = rows2vars(readtable([userpath,filesep,'VOGA_VerInfo.txt'],'ReadVariableNames',false,'ReadRowNames',true));
Data.info.Analyzer = VOGA_VerInfo.Experimenter{:};
Data.info.ver = VOGA_VerInfo.Version{:};
% Get or set the file name
if ~isfield(Data.info,'name')
    Data.info.name = [Data.info.subject,'-',Data.info.visit,'-',Data.info.exp_date,'-',Data.info.goggle_ver,'-',Data.info.dataType,'.mat'];
end
fname = Data.info.name;
% Load default filters
filt_params = VOGA__saveLastUsedParams;
%For speed during autoscan analysis, do not make cycle average figures
buffer_pix = [25,50,450,100]; %left, bottom, right, and top # of pixel buffer
plot_info.screen_size = get(0,"ScreenSize");
plot_info.fig_space = [buffer_pix(1:2),plot_info.screen_size(3:4)-buffer_pix(3:4)-buffer_pix(1:2)];
plot_info.menu_space = [sum(plot_info.fig_space([1,3]))+buffer_pix(1),...
    plot_info.fig_space(2),buffer_pix(3)-2*buffer_pix(1),plot_info.fig_space(4)];
fig_title = strrep(strrep(strrep(fname,'_',' '),'-',' '),'.mat','');
fig_title(strfind(fig_title,'['):strfind(fig_title,']')) = strrep(fig_title(strfind(fig_title,'['):strfind(fig_title,']')),' ','-'); 
close all;
figure('Color','w','Position',plot_info.fig_space)
annotation('textbox',[0 .9 1 .1],'String',fig_title,'FontSize',14,'HorizontalAlignment','center','EdgeColor','none');
all_traces = {'LX','RX','LY','RY','LZ','RZ','LLARP','RLARP','LRALP','RRALP'};
traces_pos1 = all_traces(1:6); %XYZ pos
if contains(fname,{'X','Y'})
    traces_vel1 = all_traces(1:6);
elseif contains(fname,{'GNO','ESC3'})
    traces_vel1 = all_traces(4+2*find(cellfun(@(x) contains(fname,x),{{'LH','RH'};{'LA','RP'};{'RA','LP'}})));
elseif contains(fname,'ESC')
    traces_vel1 = all_traces(3+2*find(cellfun(@(x) contains(fname,x),{{'LH','RH'};{'LA','RP'};{'RA','LP'}})));
else
    traces_vel1 = all_traces(5:end);
end
type_ord = [3,4,2,1]; %Decision tree of which type it is
type_logic = [contains(fname,{'GNO','ESC'})&&~contains(fname,'ESC3'),...
    contains(fname,'Impulse'),contains(fname,{'Activation','Step'}),true];
type = type_ord(find(type_logic,1,'first'));
Data.info.type = type;
plot_info.traces_pos = traces_pos1;
plot_info.traces_vel = traces_vel1;
plot_info.YLim.Pos = [NaN NaN];
plot_info.YLim.Vel = [NaN NaN];
Data.info.TriggerShift2 = 0;
Data = MakeCycAvg__alignCycles(Data); % Cycle Align  
traces_pos = plot_info.traces_pos;
traces_vel = plot_info.traces_vel;
tr2cyc = @(tr) tr(Data.keep_inds);
T = size(Data.keep_inds,1);
te = Data.te;
ts = Data.ts;
t_snip = Data.t_snip;
stim = Data.stim;
stims = Data.stims;
fname = Data.info.name;
keep_inds = Data.keep_inds;
keep_tr = true(1,size(keep_inds,2));
filt = filt_params.filt1;
filt.t_interp = [];
filt.keep_tr = keep_tr;
%% Position
pos_thresh = 50; %nothing should be out of these bounds
dropout_thresh = 0.2;
error_thresh = 0.01;
spline_param = [1,1-5*10.^(-11:-1)]; %spline values to try
Data_pos_filt = Data;
[~,Data_pos] = angpos2angvel(Data);
rel_inds = @(trace) find((diff(trace(1:end-1))>dropout_thresh&diff(trace(2:end))<-dropout_thresh)|...
        (diff(trace(1:end-1))<-dropout_thresh&diff(trace(2:end))>dropout_thresh))+1;
for i = 1:length(traces_pos)
    raw_tr = Data.([traces_pos{i}(1),'E_Position_',traces_pos{i}(2)]);
    raw_tr(abs(raw_tr)>pos_thresh|isnan(raw_tr)) = 0;
    tr = raw_tr;
    % Check for dropouts by looking for indeces with rapid changes in a direction and then back. 
    % Iteratively 1-D median filter using a frame length of 3 over those indeces.
    iter = 1;
    while any(rel_inds(tr))&&iter<10
        iter = iter+1;
        tr1 = medfilt1(tr,3);
        tr(rel_inds(tr)) = tr1(rel_inds(tr));
        filt.pos.median(traces_pos{i}) = 3;
    end
    %Remove point-by-point noise with spline
    spline_tr = NaN(length(tr),length(spline_param));
    for j = 1:length(spline_param)
        spline_tr(:,j) = filterTrace('spline',tr,spline_param(j),Data.te,Data.te);
    end
    % Metrics are bumpiness (zero crossings in the derivative), and root mean square error.
    [~,counts] = zerocrossrate(gradient(spline_tr));
    counts = counts+sum(diff(spline_tr)==0);
    err = [counts/length(te);sqrt(mean((tr-spline_tr).^2))/sqrt(mean((raw_tr-mean(raw_tr)).^2))]';
    [~,temp] = zerocrossrate(gradient(tr));
    raw_bump = sum(diff(tr)==0)+temp;
    if median(abs(diff(tr)))>0.05&&raw_bump/length(te)>0.3 %Prioritize filtering
        [~,ind] = min(sum(err,2));
        ind = ind-1;
    else %Prioritize not overfiltering
        ind = find(err(:,1)<=err(1,1)&err(:,2)<=error_thresh,1,'last');
    end
    filt.pos.spline(traces_pos{i}) = spline_param(ind);
    tr = spline_tr(:,ind);
    Data_pos_filt.([traces_pos{i}(1),'E_Position_',traces_pos{i}(2)]) = tr;
    if mode
        plot(ts,stim,'k',te,raw_tr,'b',te,tr,'r'); pause;
    end
end
eye_max = NaN(length(plot_info.traces_pos),2);
for t = 1:length(plot_info.traces_pos)
    tr = strrep(strrep([plot_info.traces_pos{t}(1),'E_Position_',plot_info.traces_pos{t}(2)],'_L','_LARP'),'_R','_RALP');
    if isfield(Data_pos_filt,tr)
        eye_max(t,:) = [min(Data_pos_filt.(tr)),max(Data_pos_filt.(tr))];
    end
end
eye_max(eye_max<0) = 10*floor(eye_max(eye_max<0)/10);
eye_max(eye_max>=0) = 10*ceil(eye_max(eye_max>=0)/10);
plot_info.YLim.Pos = [max([min(eye_max(:,1)),-40]),min([max(eye_max(:,2)),40])]; %Keep bounded by ±40
%% Velocity
Data_vel = angpos2angvel(Data_pos_filt);
Data_vel.t = ts;
Data_vel.stim = stim;
Data_vel_filt.t = ts;
Data_vel_filt.stim = stim;
Data_cyc.t = t_snip;
Data_cyc.stim = stims;
Data_cyc.keep_inds = keep_inds;
spline_param = [1,1-5*10.^(-11:-1),0]; %spline values to try
ip = round(T*0.16);
error_thresh = 0.95;
for i = 1:length(traces_vel)
    name_tr = strrep(strrep([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2)],'_L','_LARP'),'_R','_RALP');
    full_tr = Data_vel.(name_tr);
    full_tr(isnan(full_tr)) = 0;
    cyc_tr = tr2cyc(full_tr);
    med_tr = median(tr2cyc(full_tr),2,'omitnan');
    is_sac = false(size(cyc_tr));
    for j = 1:size(cyc_tr,2)
        is_sac(:,j) = islocalmax(abs(cyc_tr(:,j)-med_tr))&abs(cyc_tr(:,j)-med_tr)>3*iqr(cyc_tr,2);
    end
    sac_ind = Data.keep_inds(is_sac);
    med_max = max(abs(med_tr));
    raw_bump = sum(islocalmax(full_tr)|islocalmin(full_tr));
    %Desaccade with irlssmooth if possible    
    if length(full_tr)>5000||med_max<10||raw_bump/length(full_tr)>0.5
       %Use this heuristic if there is a small response, high noise, or long trace
        irlsp = ip;    
    else
        err = NaN(ip+1,2);
        irls_tr = NaN(length(full_tr),ip+1);
        for j = 1:ip+1
            irls_tr(:,j) = irlssmooth(full_tr,j-1);
            err(j,1) = max(abs(median(tr2cyc(irls_tr(:,j)),2,'omitnan')))./med_max; %Don't decrease max value
            err(j,2) = 1-median(irls_tr(sac_ind,j)./full_tr(sac_ind)); %Do decrease saccade height
        end
        err(err(:,1)<error_thresh,:) = NaN;
        [~,ind] = max(sum(err,2));
        irlsp = ind-1;
    end
    filt.vel.irlssmooth(traces_vel{i}) = irlsp;
    filt_tr = irlssmooth(full_tr,irlsp);
    %Remove point-by-point noise with spline
    spline_tr = NaN(length(filt_tr),length(spline_param));
    err = NaN(length(spline_param),2);
    filt_bump = sum(islocalmax(full_tr)|islocalmin(full_tr));
    med_max = max(abs(median(tr2cyc(filt_tr),2,'omitnan')));
    for j = 1:length(spline_param)
        spline_tr(:,j) = filterTrace('spline',filt_tr,spline_param(j),Data.te,Data.te);
        err(j,1) = sum(islocalmax(spline_tr(:,j))|islocalmin(spline_tr(:,j)))/filt_bump;
        err(j,2) = 1-max(abs(median(tr2cyc(spline_tr(:,j)),2,'omitnan')))./med_max; %Don't decrease max value
    end
    %err(:,2) = sqrt(mean((filt_tr-spline_tr).^2))/sqrt(mean((filt_tr-mean(filt_tr)).^2));
    %[~,ind] = min(sum(err,2)); %Minimize # of peaks and total RMSE
    err(err(:,2)>(1-error_thresh),:) = NaN;
    [~,ind] = min(sum(err,2)); %Minimize # of peaks and total RMSE
    filt.vel.spline(traces_vel{i}) = spline_param(ind);
    filt_tr = filterTrace('spline',filt_tr,spline_param(ind),Data.te,Data.te);
    Data_vel_filt.(name_tr) = filt_tr;
    Data_cyc.(name_tr) = Data_vel_filt.(name_tr)(keep_inds);
    if mode
        plot(ts,stim,'k',ts,full_tr,'b',ts,filt_tr,'r')
        set(gca,'YLim',300*[-1,1])
        pause;
    end
end
%Set Velocity YLim
eye_max = NaN(length(plot_info.traces_vel),2);
for t = 1:length(plot_info.traces_vel)
    tr = strrep(strrep([plot_info.traces_vel{t}(1),'E_Vel_',plot_info.traces_vel{t}(2)],'_L','_LARP'),'_R','_RALP');
    if isfield(Data_cyc,tr)
        eye_max(t,:) = [min(1.5*median(Data_cyc.(tr),2)),max(1.5*median(Data_cyc.(tr),2))];
    end
end
if contains(Data.info.name,'Impulse')
    eye_max(t+1,:) = [min(Data_vel_filt.stim),max(Data_vel_filt.stim)];
end
eye_max(eye_max<0) = 25*floor(eye_max(eye_max<0)/25);
eye_max(eye_max>=0) = 25*ceil(eye_max(eye_max>=0)/25);
plot_info.YLim.Vel = [min(eye_max(:,1)),max(eye_max(:,2))];
if plot_info.YLim.Vel(1)==plot_info.YLim.Vel(2)
    plot_info.YLim.Vel(2) = plot_info.YLim.Vel(1)+1;
end
%% Make the CycAvg Struct
%Data Traces
CycAvg.t = Data_cyc.t;
if all(size(Data_cyc.stim)>1) %multiple head traces
    CycAvg.stim_cyc = Data_cyc.stim(:,keep_tr);
    CycAvg.stim_cycavg = mean(Data_cyc.stim(:,keep_tr),2,'omitnan')';
    CycAvg.stim_cycstd = std(Data_cyc.stim(:,keep_tr),0,2,'omitnan')';
    CycAvg.stim = CycAvg.stim_cyc;
else
    CycAvg.stim = Data_cyc.stim';
end
traces = traces_vel;
for i = 1:length(traces)
    trac = lower(traces{i}(1:2));
    var_n = [traces{i}(1),'E_Vel_',traces{i}(2:end)];
    if isfield(Data_cyc,var_n)
        if all(size(Data_cyc.(var_n)(:,keep_tr))>1)
            CycAvg.([trac,'_cycavg']) = mean(Data_cyc.(var_n)(:,keep_tr),2,'omitnan')';
            CycAvg.([trac,'_cycstd']) = std(Data_cyc.(var_n)(:,keep_tr),0,2,'omitnan')';
            CycAvg.([trac,'_cyc']) = Data_cyc.(var_n)(:,keep_tr)';
        else
            CycAvg.([trac,'_cyc']) = Data_cyc.(var_n)';
        end
        if contains(fname,{'Activation','Step'})
            CycAvg.([trac,'_cyc_prefilt']) = Data_vel.(var_n);
        elseif contains(fname,'Impulse')&&isfield(Data_vel,var_n)
            CycAvg.([trac,'_cyc_prefilt']) = Data_vel.(var_n)(Data_cyc.keep_inds(:,keep_tr));
            CycAvg.([trac,'_cyc_QP']) = Data_cyc.([var_n,'_smooth'])(:,keep_tr);
            long_t = repmat(Data_cyc.t,1,sum(keep_tr));
            CycAvg.([trac,'_saccade_time']) = long_t(logical(reshape(Data_cyc.([var_n,'_saccade'])(:,keep_tr),[],1)));
        end
    end
end
%File Information and intermediate steps
%Other relevant items
CycAvg.Fs = Data.Fs;
CycAvg.info = Data.info;
CycAvg.name = ['CycAvg_',fname];
CycAvg.filt = filt;
CycAvg.cyclist = find(keep_tr);
CycAvg.keep_tr = keep_tr;
CycAvg.detec_tr = Data.detec_tr;
CycAvg.t_interp = filt.t_interp;
%Steps in Data Analsis
CycAvg.Data = Data; %Direct output of segment file
CycAvg.Data_rawpos = Data_pos;
CycAvg.Data_filtpos = Data_pos_filt;
CycAvg.Data_rawvel = Data_vel;
CycAvg.Data_filtvel = Data_vel_filt;
CycAvg.Data_allcyc = Data_cyc;
%% Cycle Selection
if length(filt.keep_tr)>5
    %Extract eye and stim data
    fields = fieldnames(Data_cyc);
    traces_vel1 = traces_vel;
    for i = 1:length(traces_vel)
        traces_vel1{i} = [traces_vel1{i}(1),'E_Vel_',traces_vel{i}(2:end)];
    end
    eye_fields = fields(contains(fields,traces_vel1)&~contains(fields,{'saccade','smooth'}));
    all_trac = cell(length(eye_fields),1);
    for i = 1:length(eye_fields)
        all_trac{i} = Data_cyc.(eye_fields{i});
    end
    trac_cyc = vertcat(all_trac{:});
    %Use the median 
    [~,ind] = sort(sqrt(mean((trac_cyc-median(trac_cyc,2,'omitnan')).^2)));
    mean_med = NaN(sum(keep_tr),1);
    for i = 3:sum(keep_tr)
        mean_med(i) = sqrt(mean((mean(trac_cyc(:,ind(1:i)),2,'omitnan')-median(trac_cyc(:,ind(1:i)),2,'omitnan')).^2))/sqrt(i);
    end
    [~,ind2] = min(mean_med);
    keep_tr(ind(ind2+1:end)) = 0;
    trac_cyc(:,~keep_tr) = NaN;
    [m_val,m_ind] = max(abs(median(trac_cyc,2,'omitnan')));    
    [~,ind] = sort(abs(trac_cyc(m_ind,:)-m_val));
    max_sem = NaN(sum(keep_tr),1);
    for i = 3:sum(keep_tr)
        max_sem(i) = max(std(trac_cyc(1:m_ind,ind(1:i)),[],2))/i^2;
    end
    [~,ind3] = min(max_sem);
    keep_tr(ind(ind3+1:end)) = 0;
    %Recalculate and plot
    CycAvg.keep_tr = keep_tr;
    CycAvg = MakeCycAvg__filterTraces([],[],[],CycAvg);
end
MakeCycAvg__plotFullCycAvg([],CycAvg,plot_info);
CycAvg = ParameterizeCycAvg(CycAvg);
filt_params.filt.pos = filt.pos;
filt_params.filt.vel = filt.vel;
filt_params.YLim = plot_info.YLim;
VOGA__saveLastUsedParams(filt_params);
analyzed = 1;
end