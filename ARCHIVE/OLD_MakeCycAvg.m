%% MakeCycAvg.m
% This function was created to automate the data analysis pipeline for the
% LDVOG data sets as much as possible.
% Cycle average before filtering
% 
% Version from 2020-10-29 before altering
function done = MakeCycAvg_v1(path,Seg_Path,Cyc_Path)
clc;       
% Standardize Colors
% Normal colors
colors.l_x = [237,150,33]/255;
colors.l_y = [125,46,143]/255;
colors.l_z = [1 0 0];
colors.l_l = [0,128,0]/255;
colors.l_r = [0 0 1];
colors.r_x = [237,204,33]/255;
colors.r_y = [125,46,230]/255;
colors.r_z = [1,0,1];
colors.r_l = [0 1 0];
colors.r_r = [64,224,208]/255;
% Faded colors
colors.l_x_s = colors.l_x + 0.5*(1-colors.l_x);
colors.l_y_s = colors.l_y + 0.5*(1-colors.l_y);
colors.l_z_s = colors.l_z + 0.5*(1-colors.l_z);
colors.l_l_s = colors.l_l + 0.5*(1-colors.l_l);
colors.l_r_s = colors.l_r + 0.5*(1-colors.l_r);
colors.r_x_s = colors.r_x + 0.5*(1-colors.r_x);
colors.r_y_s = colors.r_y + 0.5*(1-colors.r_y);
colors.r_z_s = colors.r_z + 0.5*(1-colors.r_z);
colors.r_l_s = colors.r_l + 0.5*(1-colors.r_l);
colors.r_r_s = colors.r_r + 0.5*(1-colors.r_r);
%Numbers for graphing
x_pos = 0.1;
y_pos1 = 0.08;
y_pos2 = 0.53;
x_wid = 0.85;
y_height = 0.4;
%% Load in data
progress_tab = assessProgress(path);
progress_i = [find(~progress_tab{:,2}&~progress_tab{:,3});find(progress_tab{:,2}|progress_tab{:,3})]; %put unanalyzed files at the top
set(0,'units','pixels')  
Pix_SS = get(0,'MonitorPositions');
pix_wid = Pix_SS(1,3);
pix_height = Pix_SS(1,4);
f = figure;
set(f,'Units','normalized','Position',[1 - 500/pix_wid, 0, 500/pix_wid, 1])
uit = uitable(f,'Data',table2cell(progress_tab(progress_i,:)),'ColumnName',progress_tab.Properties.VariableNames);
uit.Units = 'normalized';
uit.Position = [20/500, 20/pix_height, 1 - 40/500, 1 - 40/pix_height];
uit.ColumnWidth = [{320},{55},{50}];
[indx,tf] = nmlistdlg('PromptString','Select a file to analyze:',...
                       'SelectionMode','single',...
                       'ListSize',[350 300],...
                       'ListString',table2cell(progress_tab(progress_i,1)));  
close(f)
%If the user selects cancel
if ~tf
    disp('Operation Ended.')
    done = true;
    return; 
end
done = false;
a = dir([Seg_Path,filesep,'*.mat']);
a = {a.name}';
a = a(progress_i);
In_FileName = a{indx};
load([Seg_Path,filesep,In_FileName],'Data');
%% Extract and plot raw position data
te = Data.Time_Eye - Data.Time_Eye(1);
ts = Data.Time_Stim - Data.Time_Stim(1);
info = Data.info;
dType = strrep(strrep(info.dataType,'_',' '),'-',' ');
name1 = [info.subject,'-',info.visit,'-',info.exp_date,'-',info.dataType,'.mat'];
Fs = Data.Fs;
%Trigger multiplier
if contains(info.dataType,'RotaryChair')
    if isfield(Data,'HeadVel_Z')
         stim = Data.HeadVel_Z;
    else
        stim = Data.HeadMPUVel_Z; 
    end   
    sm = 1;
else
    stim = Data.Trigger; 
    sm = 20;
end       
%Fix huge number of NaN values in torsion traces of NKI traces
if contains(info.goggle_ver,'NKI')
    if sum(isnan(Data.LE_Position_X)) > 0.9*length(te) %less than 10% data integrity
        Data.LE_Position_X = zeros(length(te),1); %set to 0 so no torsion
    else
        Data.LE_Position_X = spline(te(~isnan(Data.LE_Position_X)),Data.LE_Position_X(~isnan(Data.LE_Position_X)),te);
    end    
    if sum(isnan(Data.RE_Position_X)) > 0.9*length(te)
        Data.RE_Position_X = zeros(length(te),1);
    else
        Data.RE_Position_X = spline(te(~isnan(Data.RE_Position_X)),Data.RE_Position_X(~isnan(Data.RE_Position_X)),te);
    end    
end

fig = figure(1);
delete(findall(gcf,'type','annotation')) %in case there are leftover anotations
fig.Units = 'inches';
fig.Position = [0 0 7 10];
%Title
annotation('textbox',[0 .9 1 .1],'String',dType,'FontSize',14,...
    'HorizontalAlignment','center','EdgeColor','none');
subplot(1,1,1)
plot(ts,sm*stim,'k')
hold on
plot(te,Data.LE_Position_X,'Color',colors.l_x)
plot(te,Data.LE_Position_Y,'Color',colors.l_y)
plot(te,Data.LE_Position_Z,'Color',colors.l_z)
plot(te,Data.RE_Position_X,'Color',colors.r_x)
plot(te,Data.RE_Position_Y,'Color',colors.r_y)
plot(te,Data.RE_Position_Z,'Color',colors.r_z)
hold off
axis([0 te(end) -20 20])
title('Raw Angular Position')
xlabel('Time(s)')
ylabel('Position (deg)')
legend('Stim','L X','L Y','L Z','R X','R Y','R Z')
%% Determine if data are analyzeable
analyze = nmquestdlg('Are these data analyzeable?','Assess data','Yes','No','Auto','Auto');
%Only cycle average if you can analyze the data
if strcmp(analyze,'No') 
    name = ['NotAnalyzeable_',name1];
    save([Cyc_Path,filesep,name],'Data')
    %If the file was previously designated as analyzeable,
    %remove that file
    d1 = dir([Cyc_Path,filesep,'*.mat']);
    if ~isempty(d1)
        d1 = {d1.name};
        if ismember(['CycAvg_',name1],d1)
            delete([Cyc_Path,filesep,'CycAvg_',name1])
        end
    end
    return;
end
%% Cycle Average
if (contains(info.dataType,'RotaryChair')||contains(info.dataType,'aHIT'))&&contains(info.dataType,'Sine')
    fparts = split(info.dataType,'-');
    freqs = fparts(contains(fparts,'Hz'));
    freq = zeros(1,length(freqs));
    for i = 1:length(freqs)
        freq(i) = str2double(strrep(freqs(i),'Hz',''));
    end
    amp = str2double(strrep(fparts{contains(fparts,'dps')},'dps',''));
    snip_len = floor(Fs/min(freq));
    template = zeros(length(freq),snip_len);  
    for i = 1:length(freq)
        template(i,:) = sin(2*pi*freq(i)*ts(1:snip_len));
    end
    template = sum(template,1);
    template = amp/max(template)*template;
    if size(stim,2)==1
        template = template';
    end
    %Find the mismatch between signal and template
    errors = NaN(1,length(stim));
    sub_i = floor(linspace(1,length(stim)-snip_len,1000));
    for i = 1:1000  
        errors(sub_i(i)) = sum((template-stim(sub_i(i):sub_i(i)+snip_len-1)).^2);
    end  
    errors = (errors - min(errors))/(max(errors)-min(errors));
    TF = find(islocalmin(errors));
    starts = TF(errors(TF)<0.1); 
    %Now find zeros crossings on the stimulus trace itself
    neg_stim = find(stim < 0);
    pos_stim = find(stim > 0);
    for i = 1:length(starts)
        if stim(starts(i)) < 0 %look for next positive value
            l_i = pos_stim(find(pos_stim>starts(i),1,'first'));
            if isempty(l_i)
                l_i = length(starts);
            end 
            poss_i = [l_i-1 l_i];
        else %look for previous negative value
            l_i = neg_stim(find(neg_stim<starts(i),1,'last'));
            if isempty(l_i)
                l_i = 1;
            end 
            poss_i = [l_i l_i+1];
        end
        [~,ind] = min(abs(stim(poss_i)));
        starts(i) = poss_i(ind); 
    end
    starts = unique(starts); %remove duplicates
    if length(starts) > 1
        snip_len1 = min(diff(starts));
        snip_len2 = length(stim) - starts(end);
        if abs(snip_len2-snip_len)/snip_len < 0.01 %Less than 1% off of the expected cycle length
            snip_len = min([snip_len1 snip_len2]);
        else
            snip_len = snip_len1; 
        end
    end 
    ends = starts + snip_len - 1;
    %Delete incomplete cycles
    starts(ends>length(stim)) = [];
    ends(ends>length(stim)) = [];
    all_stim = zeros(snip_len,length(starts));
    keep_inds = zeros(snip_len,length(starts));
    for i = 1:length(starts)
        keep_inds(:,i) = starts(i):ends(i);
        all_stim(:,i) = stim(starts(i):ends(i));
    end
    %Manually remove any erroneous motion traces
    tol = 0.2; %Amplitude can be 20% wrong and still be tolerates
    rm_tr = abs(max(all_stim)-amp)/amp > tol | abs(min(all_stim)+amp)/amp > tol;
    plot(1:snip_len,all_stim(:,~rm_tr),'k') %Samples for trace selection
    xlabel('Samples')
    ylabel('Angular Velocity (dps)')
    title(['Aligned Stimulus Traces: ',num2str(sum(~rm_tr)),' cycles'])
    remove = nmquestdlg('Remove an erroneous stimulus trace by clicking on plot?','','Yes','No','No');
    while strcmp(remove,'Yes')
        [x,y] = ginput(1);
        [~,i1] = min(abs(all_stim(round(x),:) - y));
        rm_tr(i1) = true;
        plot(1:snip_len,all_stim(:,~rm_tr),'k') 
        xlabel('Samples')
        ylabel('Angular Velocity (dps)')
        title(['Aligned Stimulus Traces: ',num2str(sum(~rm_tr)),' cycles'])
        remove = nmquestdlg('Remove an erroneous stimulus trace by clicking on plot?','','Yes','No','No');
    end
    starts(rm_tr) = [];
    ends(rm_tr) = [];
    all_stim(:,rm_tr) = [];
    keep_inds(:,rm_tr) = [];
    stims = mean(all_stim,2);         
elseif contains(info.dataType,'RotaryChair') && contains(info.dataType,'Step')
    t_snip = ts - ts(1);
    stims = stim;
    starts = 1;
    ends = ts(end);
    snip_len = length(ts);
    keep_inds = 1:length(ts);
elseif contains(info.dataType,'65Vector')||contains(info.dataType,'MultiVector')
    %The trigger is actually showing when the trapezoids start and end. There
    %are only 20 cycles of the stimulus applied and there are 40 trigger
    %toggles. To make the stimulus trace, I assumed the trigger was high when
    %the virtual velocity was linear and that the max virtual velocity is 50
    %dps (estimated from PJB's 2019 manuscript).
    %The above was confirmed by the the VOG code.
    %I also want both excitation and inhibition phases of the stimulus during
    %alignemnt.
    %Find time window for alignment
    trig = diff(stim);
    starts = find(trig==1)-1;
    starts = starts(1:2:end);
    snip_len = min(diff(starts));
    ends = starts + snip_len;
    keep_inds = zeros(snip_len,length(starts));
    for i = 1:length(starts)
        keep_inds(:,i) = starts(i):ends(i);
    end
    %Create model stimulus trace
    stims = stim(starts(1):ends(1));
    ind = find(stims==1);
    end1 = ind(diff(ind)>1);
    start2 = ind(find(diff(ind)>1)+1);
    stims(ind(1):end1) = linspace(0,50,length(ind(1):end1));
    stims(end1+1:start2-1) = 50*ones(length(end1+1:start2-1),1);
    stims(start2:ind(end)) = linspace(50,0,length(start2:ind(end)));
else %pulse train and frequency sweep
    if contains(info.dataType,'Sine') %sine (toggle = new cycle)
        trig = abs(diff(stim));
        starts = find(trig==1);
        snip_len = min(diff(starts));
        ends = starts + snip_len - 1;
        keep_inds = zeros(snip_len,length(starts));
        for i = 1:length(starts)
            keep_inds(:,i) = starts(i):ends(i);
        end
        fparts = split(info.dataType,'-');
        freq = str2double(strrep(fparts{contains(fparts,'Hz')},'Hz',''));
        amp = str2double(strrep(fparts{contains(fparts,'dps')},'dps',''));
        stims = amp*sin(2*pi*freq*(ts(1:snip_len)));
    else %pulse train or autoscan (high = on, low = off)
        trig = (diff(stim));
        starts = find(trig==1);
        snip_len = min(diff(starts));
        ends = starts + snip_len - 1;
        keep_inds = zeros(snip_len,length(starts));
        for i = 1:length(starts)
            keep_inds(:,i) = starts(i):ends(i);
        end
        stims = 50*stim(starts(1):ends(1));
    end
end 
%Remove any unnecessary trace the start and end
te = te(starts(1):ends(end)) - te(starts(1));
ts = ts(starts(1):ends(end)) - ts(starts(1));
Data.LE_Position_X = Data.LE_Position_X(starts(1):ends(end));
Data.LE_Position_Y = Data.LE_Position_Y(starts(1):ends(end));
Data.LE_Position_Z = Data.LE_Position_Z(starts(1):ends(end));
Data.RE_Position_X = Data.RE_Position_X(starts(1):ends(end));
Data.RE_Position_Y = Data.RE_Position_Y(starts(1):ends(end));
Data.RE_Position_Z = Data.RE_Position_Z(starts(1):ends(end));
keep_inds = keep_inds - starts(1)+1;
%% Filter Traces
%Offer a median and spline filter for each eye's position trace
%For all traces but rotary chair step and activation data, offer 
%irlssmooth and spline for velocity traces

%Initialize parameters
pos_med = [11;11;3;3;3;3];
pos_spline = [0.99995;0.99995;0.999995;0.999995;0.9999995;0.9999995];
vel_smooth = 30;
vel_spline = 1;
vel_acc = 5;   
init_filt_params = [pos_med;vel_smooth;vel_spline;pos_spline;vel_acc];
filt_params = init_filt_params;
while ~isempty(filt_params)
    [filt_params,filt,Data_calc] = filter_traces(filt_params,info,te,ts,stims,Data,snip_len,keep_inds,colors,analyze);
end    
%% Select Cycles      
if contains(info.dataType,'RotaryChair') && contains(info.dataType,'Step')
    LZ = qp_LZ;
    RZ = qp_RZ;
    base_labs = {'GyroZ','LX','RX','LY','RY','LZ','RZ','LLARP','RLARP','LRALP','RRALP'};
    proceed = 'Fit Again';
    while ~strcmp(proceed,'Proceed')
        h(1) = plot(ts,stim,'k');
        hold on
        h(2) = plot(ts,LE_Vel_X,'.','Color',colors.l_x);
        h(4) = plot(ts,LY,'.','Color',colors.l_y);
        h(3) = plot(ts,RX,'.','Color',colors.r_x);
        h(5) = plot(ts,RY,'.','Color',colors.r_y);
        h(8) = plot(ts,LL,'.','Color',colors.l_l);
        h(10) = plot(ts,LR,'.','Color',colors.l_r);
        h(9) = plot(ts,RL,'.','Color',colors.r_l);
        h(11) = plot(ts,RR,'.','Color',colors.r_r);
        h(6) = plot(ts,LZ,'.','Color',colors.l_z);
        h(7) = plot(ts,RZ,'.','Color',colors.r_z);
        hold off
        legend(h,base_labs)
        indx = nmlistdlg('PromptString','Select which traces should be fit to an exponential:',...
                                   'ListSize',[300 300],...
                                   'ListString',{'X','Y','LHRH','LARP','RALP'});
        eyes = questdlg('Which eye?','','Both','Left','Right','Both');
        uiwait(msgbox('Select the start and the end points for the exponential fit'))
        [x,~] = ginput(2);
        [~,i1] = min(abs(ts - x(1)));
        [~,i2] = min(abs(ts - x(2)));
        t_sub = ts(i1:i2);
        nonan_l = find(~isnan(LZ)==1);
        is_l = nonan_l(find(nonan_l>=i1,1,'first'));
        nonan_r = find(~isnan(RZ)==1);
        is_r = nonan_r(find(nonan_r>=i1,1,'first'));
        exp_fit = @(tt,p) p(1)*exp(-(tt-t_sub(1))/p(2)) + p(3);
        all_fit = zeros(5,length(ts));
        all_params = NaN(3,5);
        L = [LE_Vel_X,LY,LZ,LL,LR];
        R = [RX,RY,RZ,RL,RR];
        for i = 1:5
            if ismember(indx,i) 
                switch eyes
                    case 'Both'
                        LSCF = @(p) sum((exp_fit(t_sub,p) - L(i1:i2,i)).^2,'omitnan') + sum((exp_fit(t_sub,p) - R(i1:i2,i)).^2,'omitnan');  
                        all_params(:,i) = fminsearch(LSCF,[mean([L(is_l,i);R(is_r,i)]);10;0]);
                    case 'Left'
                        LSCF = @(p) sum((exp_fit(t_sub,p) - L(i1:i2,i)).^2,'omitnan');  
                        all_params(:,i) = fminsearch(LSCF,[mean(L(is_l,i));10;0]);
                    case 'Right'
                        LSCF = @(p) sum((exp_fit(t_sub,p) - R(i1:i2,i)).^2,'omitnan');  
                        all_params(:,i) = fminsearch(LSCF,[mean(R(is_r,i));10;0]);
                end
                all_fit(i,:) = exp_fit(ts,all_params(:,i)); 
            else
                all_fit(i,:) = mean((L(:,i) + R(:,i))/2,'omitnan')*ones(1,length(ts));
                all_params(:,i) = [0;0;mean((L(:,i) + R(:,i))/2,'omitnan')];
            end
            all_fit(i,1:i1-1) = NaN;
        end
        X_fit = all_fit(1,:);
        Y_fit = all_fit(2,:);
        LHRH_fit = all_fit(3,:);
        LARP_fit = all_fit(4,:);
        RALP_fit = all_fit(5,:);
        %Check fits
        h2(1) = plot(ts,stim,'k');
        hold on
        h2(2) = plot(ts,LE_Vel_X,'.','Color',colors.l_x);
        h2(4) = plot(ts,LY,'.','Color',colors.l_y);
        h2(3) = plot(ts,RX,'.','Color',colors.r_x);
        h2(5) = plot(ts,RY,'.','Color',colors.r_y);
        h2(8) = plot(ts,LL,'.','Color',colors.l_l);
        h2(10) = plot(ts,LR,'.','Color',colors.l_r);
        h2(9) = plot(ts,RL,'.','Color',colors.r_l);
        h2(11) = plot(ts,RR,'.','Color',colors.r_r);
        h2(6) = plot(ts,LZ,'.','Color',colors.l_z);
        h2(7) = plot(ts,RZ,'.','Color',colors.r_z);
        h2(12) = plot(ts,X_fit,'Color',colors.l_x,'LineWidth',3);
        h2(13) = plot(ts,Y_fit,'Color',colors.l_y,'LineWidth',3);
        h2(15) = plot(ts,LARP_fit,'Color',colors.l_l,'LineWidth',3);
        h2(16) = plot(ts,RALP_fit,'Color',colors.l_r,'LineWidth',3);
        h2(14) = plot(ts,LHRH_fit,'Color',colors.l_z,'LineWidth',3);
        hold off
        legend(h2,[base_labs,{'X Fit','Y Fit','Z Fit','LARP Fit','RALP Fit'}])
        %Create the params table
        params = array2table(all_params);
        params.Properties.VariableNames = {'X','Y','LHRH','LARP','RALP'};
        params.Properties.RowNames = {'Magnitude','Time_Constant','Offset'};
        disp(params)
        proceed = nmquestdlg('Proceed or Fit Again?','','Proceed','Fit Again','Proceed');
    end
    %Set mean as exp fit (same for both eyes)
    m_LE_V_LHRH = LHRH_fit';
    m_LE_V_LARP = LARP_fit';
    m_LE_V_RALP = RALP_fit';
    m_LE_V_X = X_fit';
    m_LE_V_Y = Y_fit';
    m_RE_V_LHRH = LHRH_fit';
    m_RE_V_LARP = LARP_fit';
    m_RE_V_RALP = RALP_fit';
    m_RE_V_X = X_fit';
    m_RE_V_Y = Y_fit';
    %Make standard deviations empty
    std_LE_V_LHRH = [];
    std_LE_V_LARP = [];
    std_LE_V_RALP = [];
    std_LE_V_X = [];
    std_LE_V_Y = [];
    std_RE_V_LHRH = [];
    std_RE_V_LARP = [];
    std_RE_V_RALP = [];
    std_RE_V_X = [];
    std_RE_V_Y = [];
    %Set 
    LE_V_LHRH_cyc = Data_calc.LE_Vel_Z;
    RE_V_LHRH_cyc = Data_calc.RE_Vel_Z;
    LE_V_LARP_cyc = Data_calc.LE_Vel_LARP;
    RE_V_LARP_cyc = Data_calc.RE_Vel_LARP;
    LE_V_RALP_cyc = Data_calc.LE_Vel_RALP;
    RE_V_RALP_cyc = Data_calc.RE_Vel_RALP;
    LE_V_X_cyc = Data_calc.LE_Vel_X;
    RE_V_X_cyc = Data_calc.RE_Vel_X;
    LE_V_Y_cyc = Data_calc.LE_Vel_Y;
    RE_V_Y_cyc = Data_calc.RE_Vel_Y;
else   
    if strcmp(analyze,'Auto') %currently removed this option
        %Create new time vec starting at 0 that is the length of the segment
        t_snip = ts(1:length(stims))-ts(1);
        t_snip = reshape(t_snip,1,[]);
        %Initialize vectors with the cycles
        LE_V_LHRH = Data_calc.LE_Vel_Z(keep_inds);
        LE_V_LARP = Data_calc.LE_Vel_LARP(keep_inds);
        LE_V_RALP = Data_calc.LE_Vel_RALP(keep_inds);
        LE_V_X = Data_calc.LE_Vel_X(keep_inds);
        LE_V_Y = Data_calc.LE_Vel_Y(keep_inds);
        RE_V_LHRH = Data_calc.RE_Vel_Z(keep_inds);
        RE_V_LARP = Data_calc.RE_Vel_LARP(keep_inds);
        RE_V_RALP = Data_calc.RE_Vel_RALP(keep_inds);
        RE_V_X = Data_calc.RE_Vel_X(keep_inds);
        RE_V_Y = Data_calc.RE_Vel_Y(keep_inds); 
        keep_tr = true(1,length(starts)); 
        LE_V_LHRH_cyc = LE_V_LHRH(:,keep_tr);
        LE_V_LARP_cyc = LE_V_LARP(:,keep_tr);
        LE_V_RALP_cyc = LE_V_RALP(:,keep_tr);
        LE_V_X_cyc = LE_V_X(:,keep_tr);
        LE_V_Y_cyc = LE_V_Y(:,keep_tr);
        RE_V_LHRH_cyc = RE_V_LHRH(:,keep_tr);
        RE_V_LARP_cyc = RE_V_LARP(:,keep_tr);
        RE_V_RALP_cyc = RE_V_RALP(:,keep_tr);
        RE_V_X_cyc = RE_V_X(:,keep_tr);
        RE_V_Y_cyc = RE_V_Y(:,keep_tr);
    else                 
        %Different ways of cycle analysis
        cyc_rem = nmquestdlg({'How would you like to do cycle selection?';...
            '1) Select On the Graph';'2) Cycle by Cycle with keyboard shortcuts'},'','1','2','2');
        confirm_cyc = 'No';        
        while(strcmp(confirm_cyc,'No'))
            %Create new time vec starting at 0 that is the length of the segment
            t_snip = ts(1:length(stims))-ts(1);
            t_snip = reshape(t_snip,1,[]);
            %Initialize vectors with the cycles
            LE_V_LHRH = Data_calc.LE_Vel_Z(keep_inds);
            LE_V_LARP = Data_calc.LE_Vel_LARP(keep_inds);
            LE_V_RALP = Data_calc.LE_Vel_RALP(keep_inds);
            LE_V_X = Data_calc.LE_Vel_X(keep_inds);
            LE_V_Y = Data_calc.LE_Vel_Y(keep_inds);
            RE_V_LHRH = Data_calc.RE_Vel_Z(keep_inds);
            RE_V_LARP = Data_calc.RE_Vel_LARP(keep_inds);
            RE_V_RALP = Data_calc.RE_Vel_RALP(keep_inds);
            RE_V_X = Data_calc.RE_Vel_X(keep_inds);
            RE_V_Y = Data_calc.RE_Vel_Y(keep_inds);                        
            if(strcmp(cyc_rem,'1'))
                keep_tr = true(1,length(starts)); 
                %Remove traces in each direction
                for i = 1:3
                    rm_tr = 'Yes';
                    while(strcmp(rm_tr,'Yes'))
                        switch i
                            case 1 %LHRH
                                LE = LE_V_LHRH;
                                RE = RE_V_LHRH;
                                col_l = colors.l_z;
                                col_r = colors.r_z;
                                plot_t = 'LHRH';
                            case 2 %LARP
                                LE = LE_V_LARP;
                                RE = RE_V_LARP;
                                col_l = colors.l_l;
                                col_r = colors.r_l;
                                plot_t = 'LARP';
                            case 3 %RALP
                                LE = LE_V_RALP;
                                RE = RE_V_RALP;
                                col_l = colors.l_r;
                                col_r = colors.r_r;
                                plot_t = 'RALP';
                        end
                        %Plot so the user can decide to delete more traces or
                        %not
                        ha1 = subplot(2,1,1);
                        ha1.Position = [x_pos y_pos2 x_wid y_height];
                        plot(t_snip,stims,'k')
                        hold on
                        %Plot deleted cycles lightly
                        if(any(~keep_tr))
                            plot(t_snip,LE_V_LHRH(:,~keep_tr),'Color',colors.l_z_s)
                            plot(t_snip,LE_V_LARP(:,~keep_tr),'Color',colors.l_l_s)
                            plot(t_snip,LE_V_RALP(:,~keep_tr),'Color',colors.l_r_s)
                            plot(t_snip,RE_V_LHRH(:,~keep_tr),'Color',colors.r_z_s)
                            plot(t_snip,RE_V_LARP(:,~keep_tr),'Color',colors.r_l_s)
                            plot(t_snip,RE_V_RALP(:,~keep_tr),'Color',colors.r_r_s)
                        end
                        %Plot kept cycles normally
                        if(any(keep_tr))
                            plot(t_snip,LE_V_LHRH(:,keep_tr),'Color',colors.l_z)
                            plot(t_snip,LE_V_LARP(:,keep_tr),'Color',colors.l_l)
                            plot(t_snip,LE_V_RALP(:,keep_tr),'Color',colors.l_r)
                            plot(t_snip,RE_V_LHRH(:,keep_tr),'Color',colors.r_z)
                            plot(t_snip,RE_V_LARP(:,keep_tr),'Color',colors.r_l)
                            plot(t_snip,RE_V_RALP(:,keep_tr),'Color',colors.r_r)
                        end
                        hold off
                        title('All Cycles')
                        axis([t_snip(1),t_snip(end),-110,110])
                        ha2 = subplot(2,1,2);
                        ha2.Position = [x_pos y_pos1 x_wid y_height];
                        plot(t_snip,stims,'k')
                        hold on
                        plot(t_snip,LE(:,keep_tr),'Color',col_l)
                        plot(t_snip,RE(:,keep_tr),'Color',col_r)
                        title(plot_t)
                        hold off
                        axis([t_snip(1),t_snip(end),-110,110])
                        rm_tr = nmquestdlg(['Remove a trace from the ',plot_t,' traces?'],'','Yes','No','No');
                        if(strcmp(rm_tr,'Yes'))
                            [x,y] = ginput(1);
                            %Find the trace the person was trying to remove
                            %from the valid traces left
                            val_tr = find(keep_tr == 1);
                            [~,x_ind] = min(abs(t_snip - x));
                            dist_tr = ([LE(x_ind,keep_tr);RE(x_ind,keep_tr)] - y).^2;
                            [~,n_ind] = min(min(dist_tr));
                            m_ind = val_tr(n_ind); 
                            %New plot to confirm whether to delete that trace
                            ha1 = subplot(2,1,1);
                            ha1.Position = [x_pos y_pos2 x_wid y_height];
                            plot(t_snip,stims,'k')
                            hold on
                            %Plot deleted cycles lightly
                            if(any(~keep_tr))
                                plot(t_snip,LE_V_LHRH(:,~keep_tr),'Color',colors.l_z_s)
                                plot(t_snip,LE_V_LARP(:,~keep_tr),'Color',colors.l_l_s)
                                plot(t_snip,LE_V_RALP(:,~keep_tr),'Color',colors.l_r_s)
                                plot(t_snip,RE_V_LHRH(:,~keep_tr),'Color',colors.r_z_s)
                                plot(t_snip,RE_V_LARP(:,~keep_tr),'Color',colors.r_l_s)
                                plot(t_snip,RE_V_RALP(:,~keep_tr),'Color',colors.r_r_s)
                            end
                            %Plot kept cycles normally
                            if(any(keep_tr))
                                plot(t_snip,LE_V_LHRH(:,keep_tr),'Color',colors.l_z)
                                plot(t_snip,LE_V_LARP(:,keep_tr),'Color',colors.l_l)
                                plot(t_snip,LE_V_RALP(:,keep_tr),'Color',colors.l_r)
                                plot(t_snip,RE_V_LHRH(:,keep_tr),'Color',colors.r_z)
                                plot(t_snip,RE_V_LARP(:,keep_tr),'Color',colors.r_l)
                                plot(t_snip,RE_V_RALP(:,keep_tr),'Color',colors.r_r)
                            end
                            %Plot trace in question in bold on both plots
                            plot(t_snip,LE_V_LHRH(:,m_ind),'Color',colors.l_z,'LineWidth',2)
                            plot(t_snip,LE_V_LARP(:,m_ind),'Color',colors.l_l,'LineWidth',2)
                            plot(t_snip,LE_V_RALP(:,m_ind),'Color',colors.l_r,'LineWidth',2)
                            plot(t_snip,RE_V_LHRH(:,m_ind),'Color',colors.r_z,'LineWidth',2)
                            plot(t_snip,RE_V_LARP(:,m_ind),'Color',colors.r_l,'LineWidth',2)
                            plot(t_snip,RE_V_RALP(:,m_ind),'Color',colors.r_r,'LineWidth',2)
                            hold off
                            title('All Cycles')
                            axis([t_snip(1),t_snip(end),-110,110])
                            ha2 = subplot(2,1,2);
                            ha2.Position = [x_pos y_pos1 x_wid y_height];
                            plot(t_snip,stims,'k')
                            hold on
                            plot(t_snip,LE(:,keep_tr),'Color',col_l)
                            plot(t_snip,RE(:,keep_tr),'Color',col_r)
                            plot(t_snip,LE(:,m_ind),'Color',col_l,'LineWidth',2)
                            plot(t_snip,RE(:,m_ind),'Color',col_r,'LineWidth',2)
                            hold off
                            title(plot_t)
                            axis([t_snip(1),t_snip(end),-110,110])
                            del_tr = nmquestdlg('Delete this cycle?','','Yes','No','No');
                            if(strcmp(del_tr,'Yes'))
                                keep_tr(m_ind) = false; 
                            end
                        end
                    end
                end
            else   
                %Remove by manual keyboard cycle selection
                keep_tr = true(1,length(starts));
                for i = 1:length(starts) 
                    subplot(1,1,1)
                    plot(t_snip,stims,'k')
                    hold on
                    plot(t_snip,LE_V_LHRH(:,keep_tr),'Color',colors.l_z)
                    plot(t_snip,LE_V_LARP(:,keep_tr),'Color',colors.l_l)
                    plot(t_snip,LE_V_RALP(:,keep_tr),'Color',colors.l_r)
                    plot(t_snip,RE_V_LHRH(:,keep_tr),'Color',colors.r_z)
                    plot(t_snip,RE_V_LARP(:,keep_tr),'Color',colors.r_l)
                    plot(t_snip,RE_V_RALP(:,keep_tr),'Color',colors.r_r)
                    %Cycle of interest plotted normally
                    plot(t_snip,LE_V_LHRH(:,i),'Color',colors.l_z,'LineWidth',2)
                    plot(t_snip,LE_V_LARP(:,i),'Color',colors.l_l,'LineWidth',2)
                    plot(t_snip,LE_V_RALP(:,i),'Color',colors.l_r,'LineWidth',2)
                    plot(t_snip,RE_V_LHRH(:,i),'Color',colors.r_z,'LineWidth',2)
                    plot(t_snip,RE_V_LARP(:,i),'Color',colors.r_l,'LineWidth',2)
                    plot(t_snip,RE_V_RALP(:,i),'Color',colors.r_r,'LineWidth',2)
                    hold off       
                    axis([t_snip(1),t_snip(end),-110,110])
                    keep = input('Keep trace? k = keep rest (y/n/k) ','s');
                    while ~strcmp(keep,'y') && ~strcmp(keep,'n') && ~strcmp(keep,'k')
                        disp('Only the characters "y", "n", and "k" are accepted as a response.')
                        keep = input('Keep trace? (y/n) ','s');
                    end
                    if strcmp(keep,'k')
                        break;
                    elseif strcmp(keep,'n')
                        keep_tr(i) = false; 
                    end
                end
            end
            LE_V_LHRH_cyc = LE_V_LHRH(:,keep_tr);
            LE_V_LARP_cyc = LE_V_LARP(:,keep_tr);
            LE_V_RALP_cyc = LE_V_RALP(:,keep_tr);
            LE_V_X_cyc = LE_V_X(:,keep_tr);
            LE_V_Y_cyc = LE_V_Y(:,keep_tr);
            RE_V_LHRH_cyc = RE_V_LHRH(:,keep_tr);
            RE_V_LARP_cyc = RE_V_LARP(:,keep_tr);
            RE_V_RALP_cyc = RE_V_RALP(:,keep_tr);
            RE_V_X_cyc = RE_V_X(:,keep_tr);
            RE_V_Y_cyc = RE_V_Y(:,keep_tr);
            ha1 = subplot(2,1,1);
            ha1.Position = [x_pos y_pos2 x_wid y_height];
            plot(t_snip,stims,'k')
            hold on
            plot(t_snip,LE_V_LHRH,'Color',colors.l_z)
            plot(t_snip,LE_V_LARP,'Color',colors.l_l)
            plot(t_snip,LE_V_RALP,'Color',colors.l_r)
            plot(t_snip,RE_V_LHRH,'Color',colors.r_z)
            plot(t_snip,RE_V_LARP,'Color',colors.r_l)
            plot(t_snip,RE_V_RALP,'Color',colors.r_r)
            hold off
            title('No Cycles Removed')
            ylabel('Angular Velocity (dps)')
            axis([t_snip(1),t_snip(end),-110,110])
            ha2 = subplot(2,1,2);
            ha2.Position = [x_pos y_pos1 x_wid y_height];
            plot(t_snip,stims,'k')
            hold on
            plot(t_snip,LE_V_LHRH_cyc,'Color',colors.l_z)
            plot(t_snip,LE_V_LARP_cyc,'Color',colors.l_l)
            plot(t_snip,LE_V_RALP_cyc,'Color',colors.l_r)
            plot(t_snip,RE_V_LHRH_cyc,'Color',colors.r_z)
            plot(t_snip,RE_V_LARP_cyc,'Color',colors.r_l)
            plot(t_snip,RE_V_RALP_cyc,'Color',colors.r_r)
            hold off
            title('Only Selected Cycles')
            xlabel('Time (s)')
            ylabel('Angular Velocity (dps)')
            axis([t_snip(1),t_snip(end),-110,110])
            confirm_cyc = nmquestdlg({[num2str(size(LE_V_LHRH_cyc,2)),' of ',num2str(size(LE_V_LHRH,2)),' cycles kept.'];...
                'Confirm cycle selection?'},'','Yes','No','Filter Again','Yes');
            if strcmp(confirm_cyc,'Filter Again')
                filt_params = init_filt_params;
                while ~isempty(filt_params)
                    [filt_params,filt,Data_calc] = filter_traces(filt_params,info,te,ts,stims,Data,snip_len,keep_inds,colors,analyze);
                end    
                confirm_cyc = 'No';
            end
        end
    end
    %% Cycle Average
    %Calculate means
    m_LE_V_LHRH = mean(LE_V_LHRH_cyc,2)';
    m_LE_V_LARP = mean(LE_V_LARP_cyc,2)';
    m_LE_V_RALP = mean(LE_V_RALP_cyc,2)';
    m_LE_V_X = mean(LE_V_X_cyc,2)';
    m_LE_V_Y = mean(LE_V_Y_cyc,2)';
    m_RE_V_LHRH = mean(RE_V_LHRH_cyc,2)';
    m_RE_V_LARP = mean(RE_V_LARP_cyc,2)';
    m_RE_V_RALP = mean(RE_V_RALP_cyc,2)';
    m_RE_V_X = mean(RE_V_X_cyc,2)';
    m_RE_V_Y = mean(RE_V_Y_cyc,2)';
    %Make standard deviations
    std_LE_V_LHRH = std(LE_V_LHRH_cyc,0,2)';
    std_LE_V_LARP = std(LE_V_LARP_cyc,0,2)';
    std_LE_V_RALP = std(LE_V_RALP_cyc,0,2)';
    std_LE_V_X = std(LE_V_X_cyc,0,2)';
    std_LE_V_Y = std(LE_V_Y_cyc,0,2)';
    std_RE_V_LHRH = std(RE_V_LHRH_cyc,0,2)';
    std_RE_V_LARP = std(RE_V_LARP_cyc,0,2)';
    std_RE_V_RALP = std(RE_V_RALP_cyc,0,2)';
    std_RE_V_X = std(RE_V_X_cyc,0,2)';
    std_RE_V_Y = std(RE_V_Y_cyc,0,2)';
end
%% Make CycAvg struct
CycAvg.old_Fs = Data.Fs;
CycAvg.Fs = Fs;
CycAvg.t = t_snip';
CycAvg.stim = stims';
%Averages
CycAvg.lz_cycavg = m_LE_V_LHRH;
CycAvg.rz_cycavg = m_RE_V_LHRH;
CycAvg.ll_cycavg = m_LE_V_LARP;
CycAvg.rl_cycavg = m_RE_V_LARP;
CycAvg.lr_cycavg = m_LE_V_RALP;
CycAvg.rr_cycavg = m_RE_V_RALP;
CycAvg.lx_cycavg = m_LE_V_X;
CycAvg.rx_cycavg = m_RE_V_X;
CycAvg.ly_cycavg = m_LE_V_Y;
CycAvg.ry_cycavg = m_RE_V_Y;
%StandardDeviations
CycAvg.lz_cycstd = std_LE_V_LHRH;
CycAvg.rz_cycstd = std_RE_V_LHRH;
CycAvg.ll_cycstd = std_LE_V_LARP;
CycAvg.rl_cycstd = std_RE_V_LARP;
CycAvg.lr_cycstd = std_LE_V_RALP;
CycAvg.rr_cycstd = std_RE_V_RALP;
CycAvg.lx_cycstd = std_LE_V_X;
CycAvg.rx_cycstd = std_RE_V_X;
CycAvg.ly_cycstd = std_LE_V_Y;
CycAvg.ry_cycstd = std_RE_V_Y;
%All cycles
CycAvg.lz_cyc = LE_V_LHRH_cyc';
CycAvg.rz_cyc = RE_V_LHRH_cyc';
CycAvg.ll_cyc = LE_V_LARP_cyc';
CycAvg.rl_cyc = RE_V_LARP_cyc';
CycAvg.lr_cyc = LE_V_RALP_cyc';
CycAvg.rr_cyc = RE_V_RALP_cyc';
CycAvg.lx_cyc = LE_V_X_cyc';
CycAvg.rx_cyc = RE_V_X_cyc';
CycAvg.ly_cyc = LE_V_Y_cyc';
CycAvg.ry_cyc = RE_V_Y_cyc';
%Raw data too
CycAvg.raw_Data = Data;
%Other relevant items
info.colors = colors;
CycAvg.info = info;
CycAvg.filt = filt;
cyc_numbers = 1:length(starts);
CycAvg.cyclist = cyc_numbers(keep_tr);
if ~(contains(info.dataType,'RotaryChair') && contains(info.dataType,'Step'))
    %% Plot
    subplot(1,1,1)
    if length(CycAvg.t) > 1000
        s = round(linspace(1,length(CycAvg.t),1000));
    else
        s = 1:length(CycAvg.t);
    end
    h = [];
    h(1) = plot(CycAvg.t(s),CycAvg.stim(s),'k');
    hold on
    %Now add the fills and standard deviations and means
    %LE-LHRH
    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s)))],colors.l_z_s)
    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s)))],colors.l_z_s)
    plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
    plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
    h(2) = plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
    %RE-LHRH
    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s)))],colors.r_z_s)
    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s)))],colors.r_z_s)
    plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
    plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
    h(3) = plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
    %LE-LARP
    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ll_cycavg(s),fliplr((CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s)))],colors.l_l_s)
    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ll_cycavg(s),fliplr((CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s)))],colors.l_l_s)
    plot(CycAvg.t(s),CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s),'Color',colors.l_l)
    plot(CycAvg.t(s),CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s),'Color',colors.l_l)
    h(4) = plot(CycAvg.t(s),CycAvg.ll_cycavg(s),'Color',colors.l_l,'LineWidth',2);
    %RE-LARP
    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rl_cycavg(s),fliplr((CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s)))],colors.r_l_s)
    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rl_cycavg(s),fliplr((CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s)))],colors.r_l_s)
    plot(CycAvg.t(s),CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s),'Color',colors.r_l)
    plot(CycAvg.t(s),CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s),'Color',colors.r_l)
    h(5) = plot(CycAvg.t(s),CycAvg.rl_cycavg(s),'Color',colors.r_l,'LineWidth',2);
    %LE_RALP
    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lr_cycavg(s),fliplr((CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s)))],colors.l_r_s)
    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lr_cycavg(s),fliplr((CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s)))],colors.l_r_s)
    plot(CycAvg.t(s),CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s),'Color',colors.l_r)
    plot(CycAvg.t(s),CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s),'Color',colors.l_r)
    h(6) = plot(CycAvg.t(s),CycAvg.lr_cycavg(s),'Color',colors.l_r,'LineWidth',2);
    %RE-RALP
    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rr_cycavg(s),fliplr((CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s)))],colors.r_r_s)
    fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rr_cycavg(s),fliplr((CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s)))],colors.r_r_s)
    plot(CycAvg.t(s),CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s),'Color',colors.r_r)
    plot(CycAvg.t(s),CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s),'Color',colors.r_r)
    h(7) = plot(CycAvg.t(s),CycAvg.rr_cycavg(s),'Color',colors.r_r,'LineWidth',2);
    hold off
    xlabel('Time (s)')
    ylabel('Angular Velocity (dps)')
    title('Cycle Average')
    legend(h,'Stim','Left Z','Right Z','Left LARP','Right LARP','Left RALP','Right RALP')  
else
    CycAvg.expfit = params;
end
%% Save
save_cycavg = nmquestdlg('Save the cycle average file?','','Yes','No','No');
if strcmp(save_cycavg,'Yes')
    name1 = [info.subject,'-',info.visit,'-',info.exp_date,'-',info.dataType,'.mat'];
    name = ['CycAvg_',name1];
    save([Cyc_Path,filesep,name],'CycAvg')
    %If the file was previously designated as un-analyzeable,
    %remove that file
    d1 = dir([Cyc_Path,filesep,'*.mat']);
    if ~isempty(d1)
        d1 = {d1.name};
        if ismember(['NotAnalyzeable_',name1],d1)
            delete([Cyc_Path,filesep,'NotAnalyzeable_',name1])
        end
    end
end
%% Make a Sphere Plot Using Peter's script
% if ~contains(info.dataType,'RotaryChair')
%     MakeSpherePlot(2,CycAvg,1,1,0);
%     pause;
% end
end

%% Functions

%Filter position and velcity
function [filt_params,filt,Data_calc] = filter_traces(filt_params,info,te,ts,stims,Data,snip_len,keep_inds,colors,analyze)
    %Deal with bad parameters
    bad_med = (filt_params<1|mod(filt_params,2)~=1)&[true(6,1);false(9,1)];
    bad_smooth = (filt_params<0|mod(filt_params,1)~=0)&[false(6,1);true;false(8,1)];
    bad_spline = (filt_params>1|filt_params<0)&[false(7,1);true(7,1);false];
    bad_acc = (filt_params<0)&[false(14,1);true];
    filt_params(bad_med) = floor(filt_params(bad_med)/2)+1;
    filt_params(bad_spline) = 1;
    filt_params(bad_smooth) = abs(floor(filt_params(bad_smooth)));
    filt_params(bad_acc) = abs(filt_params(bad_acc));
    %Save to structures
    filt.pos.l_x.medianfilt = filt_params(1);
    filt.pos.r_x.medianfilt = filt_params(2);
    filt.pos.l_y.medianfilt = filt_params(3);
    filt.pos.r_y.medianfilt = filt_params(4);
    filt.pos.l_z.medianfilt = filt_params(5);
    filt.pos.r_z.medianfilt = filt_params(6);
    filt.vel.irlssmooth = filt_params(7);
    filt.vel.spline = filt_params(8);
    filt.pos.l_x.spline = filt_params(9);
    filt.pos.r_x.spline = filt_params(10);
    filt.pos.l_y.spline = filt_params(11);
    filt.pos.r_y.spline = filt_params(12);
    filt.pos.l_z.spline = filt_params(13);
    filt.pos.r_z.spline = filt_params(14);
    filt.vel.accelcutoff = filt_params(15);
    %Filter position
    LE_X = medfilt1(Data.LE_Position_X,filt.pos.l_x.medianfilt,'omitnan');
    RE_X = medfilt1(Data.RE_Position_X,filt.pos.r_x.medianfilt,'omitnan');
    LE_Y = medfilt1(Data.LE_Position_Y,filt.pos.l_y.medianfilt,'omitnan');
    RE_Y = medfilt1(Data.RE_Position_Y,filt.pos.r_y.medianfilt,'omitnan');
    LE_Z = medfilt1(Data.LE_Position_Z,filt.pos.l_z.medianfilt,'omitnan');
    RE_Z = medfilt1(Data.RE_Position_Z,filt.pos.r_z.medianfilt,'omitnan');
    LE_P_X_f = spline_filt(te,LE_X,ts,filt.pos.l_x.spline);
    RE_P_X_f = spline_filt(te,RE_X,ts,filt.pos.r_x.spline);
    LE_P_Y_f = spline_filt(te,LE_Y,ts,filt.pos.l_y.spline);
    RE_P_Y_f = spline_filt(te,RE_Y,ts,filt.pos.r_y.spline);
    LE_P_Z_f = spline_filt(te,LE_Z,ts,filt.pos.l_z.spline);
    RE_P_Z_f = spline_filt(te,RE_Z,ts,filt.pos.r_z.spline);
    % Recalculate Angular Velocity
    Data_In.Data_LE_Pos_X = LE_P_X_f;
    Data_In.Data_LE_Pos_Y = LE_P_Y_f;
    Data_In.Data_LE_Pos_Z = LE_P_Z_f;
    Data_In.Data_RE_Pos_X = RE_P_X_f;
    Data_In.Data_RE_Pos_Y = RE_P_Y_f;
    Data_In.Data_RE_Pos_Z = RE_P_Z_f;
    Data_In.Fs = Data.Fs;
    %First param says no initial coordinate transforms, second is a 
    %struct with the filtered position traces.
    Data_cal = angpos2angvel(1,Data_In);  
    if ~strcmp(analyze,'Auto')
        % Make the plot
        x_pos = 0.1;
        y_pos1 = 0.08;
        y_pos2 = 0.53;
        x_wid = 0.85;
        y_height = 0.4;
        ha1 = subplot(2,1,1);
        ha1.Position = [x_pos y_pos2 x_wid y_height];        
        h(1) = plot(ts(1:snip_len),stims,'k');
        hold on
        %Raw Data
        plot(te(1:snip_len),Data.LE_Position_X(keep_inds),'Color',colors.l_x_s)
        plot(te(1:snip_len),Data.LE_Position_Y(keep_inds),'Color',colors.l_y_s)
        plot(te(1:snip_len),Data.LE_Position_Z(keep_inds),'Color',colors.l_z_s)
        plot(te(1:snip_len),Data.RE_Position_X(keep_inds),'Color',colors.r_x_s)
        plot(te(1:snip_len),Data.RE_Position_Y(keep_inds),'Color',colors.r_y_s)
        plot(te(1:snip_len),Data.RE_Position_Z(keep_inds),'Color',colors.r_z_s)
        %Filtered Data
        plot(ts(1:snip_len),LE_P_X_f(keep_inds),'Color',colors.l_x);
        plot(ts(1:snip_len),LE_P_Y_f(keep_inds),'Color',colors.l_y);
        plot(ts(1:snip_len),LE_P_Z_f(keep_inds),'Color',colors.l_z);
        plot(ts(1:snip_len),RE_P_X_f(keep_inds),'Color',colors.r_x);
        plot(ts(1:snip_len),RE_P_Y_f(keep_inds),'Color',colors.r_y);
        plot(ts(1:snip_len),RE_P_Z_f(keep_inds),'Color',colors.r_z);
        h(2) = plot(NaN,NaN,'Color',colors.l_x);
        h(3) = plot(NaN,NaN,'Color',colors.l_y);
        h(4) = plot(NaN,NaN,'Color',colors.l_z);
        h(5) = plot(NaN,NaN,'Color',colors.l_l);
        h(6) = plot(NaN,NaN,'Color',colors.l_r);
        h(7) = plot(NaN,NaN,'Color',colors.r_x);
        h(8) = plot(NaN,NaN,'Color',colors.r_y);
        h(9) = plot(NaN,NaN,'Color',colors.r_z);  
        h(10) = plot(NaN,NaN,'Color',colors.r_l);
        h(11) = plot(NaN,NaN,'Color',colors.r_r);
        hold off
        axis([0 ts(snip_len) -20 20])
        title('Angular Position')
        ylabel('Position (deg)')
        set(ha1,'XTickLabel',[])
        legend(h,{'Trigger','L X','L Y','L Z','L L','L R','R X','R Y','R Z','R L','R R'})
        %Filter and plot velocity
        if contains(info.dataType,'RotaryChair') && contains(info.dataType,'Step') %Add case for activation too
            %Turn points into dots and desaccade that way
            Data_calc = Data_cal;
            %Left Eye QPR
            d_LZ = [0;diff(Data_calc.LE_Vel_Z)];
            L_up_sacc = find(d_LZ>filt.vel.accelcutoff);
            L_up_sacc_s = L_up_sacc([true;diff(L_up_sacc) > 1]);
            for i = 1:length(L_up_sacc_s)
                Data_calc.LE_Vel_Z(find(d_LZ(1:L_up_sacc_s(i)-1)<0,1,'last'):find(d_LZ(L_up_sacc_s(i)+1:end)<0,1,'first')+L_up_sacc_s(i)) = NaN;
            end
            L_down_sacc = find(d_LZ<-filt.vel.accelcutoff);
            L_down_sacc_s = L_down_sacc([true;diff(L_down_sacc) > 1]);
            for i = 1:length(L_down_sacc_s)
                Data_calc.LE_Vel_Z(find(d_LZ(1:L_down_sacc_s(i)-1)>0,1,'last'):find(d_LZ(L_down_sacc_s(i)+1:end)>0,1,'first')+L_down_sacc_s(i)) = NaN;
            end
            LZ_nan = find(isnan(Data_calc.LE_Vel_Z));
            small_spac = find(diff(LZ_nan) > 1 & diff(LZ_nan) < 3);
            for i = 1:length(small_spac)
                Data_calc.LE_Vel_Z(LZ_nan(small_spac(i)):LZ_nan(small_spac(i)+1)) = NaN;
            end
            Data_calc.LE_Vel_LARP(isnan(Data_calc.LE_Vel_Z)) = NaN;
            Data_calc.LE_Vel_RALP(isnan(Data_calc.LE_Vel_Z)) = NaN;
            Data_calc.LE_Vel_X(isnan(Data_calc.LE_Vel_Z)) = NaN;
            Data_calc.LE_Vel_Y(isnan(Data_calc.LE_Vel_Z)) = NaN;
            %Right eye QPR
            d_RZ = [0;diff(Data_calc.RE_Vel_Z)];
            R_up_sacc = find(d_RZ>filt.vel.accelcutoff);
            R_up_sacc_s = R_up_sacc([true;diff(R_up_sacc) > 1]);
            for i = 1:length(R_up_sacc_s)
                Data_calc.RE_Vel_Z(find(d_RZ(1:R_up_sacc_s(i)-1)<0,1,'last'):find(d_RZ(R_up_sacc_s(i)+1:end)<0,1,'first')+R_up_sacc_s(i)) = NaN;
            end
            R_down_sacc = find(d_RZ<-filt.vel.accelcutoff);
            R_down_sacc_s = R_down_sacc([true;diff(R_down_sacc) > 1]);
            for i = 1:length(R_down_sacc_s)
                Data_calc.RE_Vel_Z(find(d_RZ(1:R_down_sacc_s(i)-1)>0,1,'last'):find(d_RZ(R_down_sacc_s(i)+1:end)>0,1,'first')+R_down_sacc_s(i)) = NaN;
            end
            RZ_nan = find(isnan(Data_calc.RE_Vel_Z));
            small_spac = find(diff(RZ_nan) > 1 & diff(RZ_nan) < 3);
            for i = 1:length(small_spac)
                Data_calc.LE_Vel_Z(RZ_nan(small_spac(i)):RZ_nan(small_spac(i)+1)) = NaN;
            end
            Data_calc.RE_Vel_LARP(isnan(Data_calc.RE_Vel_Z)) = NaN;
            Data_calc.RE_Vel_RALP(isnan(Data_calc.RE_Vel_Z)) = NaN;
            Data_calc.RE_Vel_X(isnan(Data_calc.RE_Vel_Z)) = NaN;
            Data_calc.RE_Vel_Y(isnan(Data_calc.RE_Vel_Z)) = NaN;
            %Plot
            ha2 = subplot(2,1,2);
            ha2.Position = [x_pos y_pos1 x_wid y_height];
            plot(ts(1:snip_len),stims,'k')
            hold on
            plot(ts(1:snip_len),Data_cal.LE_Vel_LARP(keep_inds),'.','Color',colors.l_l_s)
            plot(ts(1:snip_len),Data_cal.LE_Vel_RALP(keep_inds),'.','Color',colors.l_r_s)
            plot(ts(1:snip_len),Data_cal.RE_Vel_LARP(keep_inds),'.','Color',colors.r_l_s)
            plot(ts(1:snip_len),Data_cal.RE_Vel_RALP(keep_inds),'.','Color',colors.r_r_s)
            plot(ts(1:snip_len),Data_cal.LE_Vel_Z(keep_inds),'.','Color',colors.l_z_s)
            plot(ts(1:snip_len),Data_cal.RE_Vel_Z(keep_inds),'.','Color',colors.r_z_s)
            plot(ts(1:snip_len),Data_calc.LE_Vel_LARP(keep_inds),'.','Color',colors.l_l)
            plot(ts(1:snip_len),Data_calc.LE_Vel_RALP(keep_inds),'.','Color',colors.l_r)
            plot(ts(1:snip_len),Data_calc.RE_Vel_LARP(keep_inds),'.','Color',colors.r_l)
            plot(ts(1:snip_len),Data_calc.RE_Vel_RALP(keep_inds),'.','Color',colors.r_r)
            plot(ts(1:snip_len),Data_calc.LE_Vel_Z(keep_inds),'.','Color',colors.l_z)
            plot(ts(1:snip_len),Data_calc.RE_Vel_Z(keep_inds),'.','Color',colors.r_z)
            hold off
            axis([0 ts(snip_len) -110 110])
            xlabel('Time(s)')
            ylabel('Velocity (dps)')
            title('Angular Velocity')
            linkaxes([ha1,ha2],'x')             
        else %All other data types
            Data_calc.LE_Vel_X = spline_filt(ts,irlssmooth(Data_cal.LE_Vel_X,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_X = spline_filt(ts,irlssmooth(Data_cal.RE_Vel_X,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.LE_Vel_Y = spline_filt(ts,irlssmooth(Data_cal.LE_Vel_Y,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_Y = spline_filt(ts,irlssmooth(Data_cal.RE_Vel_Y,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.LE_Vel_Z = spline_filt(ts,irlssmooth(Data_cal.LE_Vel_Z,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_Z = spline_filt(ts,irlssmooth(Data_cal.RE_Vel_Z,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.LE_Vel_LARP = spline_filt(ts,irlssmooth(Data_cal.LE_Vel_LARP,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_LARP = spline_filt(ts,irlssmooth(Data_cal.RE_Vel_LARP,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.LE_Vel_RALP = spline_filt(ts,irlssmooth(Data_cal.LE_Vel_RALP,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_RALP = spline_filt(ts,irlssmooth(Data_cal.RE_Vel_RALP,filt.vel.irlssmooth),ts,filt.vel.spline);
            %Plot
            ha2 = subplot(2,1,2);
            ha2.Position = [x_pos y_pos1 x_wid y_height];
            plot(ts(1:snip_len),stims,'k')
            hold on
            plot(ts(1:snip_len),Data_cal.LE_Vel_LARP(keep_inds),'Color',colors.l_l_s)
            plot(ts(1:snip_len),Data_cal.LE_Vel_RALP(keep_inds),'Color',colors.l_r_s)
            plot(ts(1:snip_len),Data_cal.RE_Vel_LARP(keep_inds),'Color',colors.r_l_s)
            plot(ts(1:snip_len),Data_cal.RE_Vel_RALP(keep_inds),'Color',colors.r_r_s)
            plot(ts(1:snip_len),Data_cal.LE_Vel_Z(keep_inds),'Color',colors.l_z_s)
            plot(ts(1:snip_len),Data_cal.RE_Vel_Z(keep_inds),'Color',colors.r_z_s)
            plot(ts(1:snip_len),Data_calc.LE_Vel_LARP(keep_inds),'Color',colors.l_l)
            plot(ts(1:snip_len),Data_calc.LE_Vel_RALP(keep_inds),'Color',colors.l_r)
            plot(ts(1:snip_len),Data_calc.RE_Vel_LARP(keep_inds),'Color',colors.r_l)
            plot(ts(1:snip_len),Data_calc.RE_Vel_RALP(keep_inds),'Color',colors.r_r)
            plot(ts(1:snip_len),Data_calc.LE_Vel_Z(keep_inds),'Color',colors.l_z)
            plot(ts(1:snip_len),Data_calc.RE_Vel_Z(keep_inds),'Color',colors.r_z)
            hold off
            axis([0 ts(snip_len) -110 110])
            xlabel('Time(s)')
            ylabel('Velocity (dps)')
            title('Angular Velocity')
            linkaxes([ha1,ha2],'x')        
        end    
        %Get new parameter values
        prompt = {['Press OK to refilter.',newline,newline,'Position Filters:',newline,newline,'Median Filter (odd)',newline,'Left Torsion:'],...
            'Right Torsion:','Left Vertical:','Right Vertical:','Left Horizontal:','Right Horizontal:',...
            ['Velocity Filters:',newline,newline,'Irlssmooth (>=0):'],'Spline (0-1):',...
            ['Press CANCEL to confirm.',newline,newline,newline,newline,'Spline (0-1)',newline,'Left Torsion:'],'Right Torsion:',...
            'Left Vertical:','Right Vertical:','Left Horizontal:','Right Horizontal:',...
            'Acceleration Threshold:'};
        dlgtitle = 'Filter parameters';
        dims = [1 25];
        definput = cellfun(@(x) num2str(x,10),num2cell(filt_params),'UniformOutput',false);
        filt_params = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,dims,definput,'on',2));
    else
        filt_params = [];
        %Filter and plot velocity
        if contains(info.dataType,'RotaryChair') && contains(info.dataType,'Step') %Add case for activation too
            %Turn points into dots and desaccade that way
            Data_calc = Data_cal;
            %Left Eye QPR
            d_LZ = [0;diff(Data_calc.LE_Vel_Z)];
            L_up_sacc = find(d_LZ>filt.vel.accelcutoff);
            L_up_sacc_s = L_up_sacc([true;diff(L_up_sacc) > 1]);
            for i = 1:length(L_up_sacc_s)
                Data_calc.LE_Vel_Z(find(d_LZ(1:L_up_sacc_s(i)-1)<0,1,'last'):find(d_LZ(L_up_sacc_s(i)+1:end)<0,1,'first')+L_up_sacc_s(i)) = NaN;
            end
            L_down_sacc = find(d_LZ<-filt.vel.accelcutoff);
            L_down_sacc_s = L_down_sacc([true;diff(L_down_sacc) > 1]);
            for i = 1:length(L_down_sacc_s)
                Data_calc.LE_Vel_Z(find(d_LZ(1:L_down_sacc_s(i)-1)>0,1,'last'):find(d_LZ(L_down_sacc_s(i)+1:end)>0,1,'first')+L_down_sacc_s(i)) = NaN;
            end
            LZ_nan = find(isnan(Data_calc.LE_Vel_Z));
            small_spac = find(diff(LZ_nan) > 1 & diff(LZ_nan) < 3);
            for i = 1:length(small_spac)
                Data_calc.LE_Vel_Z(LZ_nan(small_spac(i)):LZ_nan(small_spac(i)+1)) = NaN;
            end
            Data_calc.LE_Vel_LARP(isnan(Data_calc.LE_Vel_Z)) = NaN;
            Data_calc.LE_Vel_RALP(isnan(Data_calc.LE_Vel_Z)) = NaN;
            Data_calc.LE_Vel_X(isnan(Data_calc.LE_Vel_Z)) = NaN;
            Data_calc.LE_Vel_Y(isnan(Data_calc.LE_Vel_Z)) = NaN;
            %Right eye QPR
            d_RZ = [0;diff(Data_calc.RE_Vel_Z)];
            R_up_sacc = find(d_RZ>filt.vel.accelcutoff);
            R_up_sacc_s = R_up_sacc([true;diff(R_up_sacc) > 1]);
            for i = 1:length(R_up_sacc_s)
                Data_calc.RE_Vel_Z(find(d_RZ(1:R_up_sacc_s(i)-1)<0,1,'last'):find(d_RZ(R_up_sacc_s(i)+1:end)<0,1,'first')+R_up_sacc_s(i)) = NaN;
            end
            R_down_sacc = find(d_RZ<-filt.vel.accelcutoff);
            R_down_sacc_s = R_down_sacc([true;diff(R_down_sacc) > 1]);
            for i = 1:length(R_down_sacc_s)
                Data_calc.RE_Vel_Z(find(d_RZ(1:R_down_sacc_s(i)-1)>0,1,'last'):find(d_RZ(R_down_sacc_s(i)+1:end)>0,1,'first')+R_down_sacc_s(i)) = NaN;
            end
            RZ_nan = find(isnan(Data_calc.RE_Vel_Z));
            small_spac = find(diff(RZ_nan) > 1 & diff(RZ_nan) < 3);
            for i = 1:length(small_spac)
                Data_calc.LE_Vel_Z(RZ_nan(small_spac(i)):RZ_nan(small_spac(i)+1)) = NaN;
            end
            Data_calc.RE_Vel_LARP(isnan(Data_calc.RE_Vel_Z)) = NaN;
            Data_calc.RE_Vel_RALP(isnan(Data_calc.RE_Vel_Z)) = NaN;
            Data_calc.RE_Vel_X(isnan(Data_calc.RE_Vel_Z)) = NaN;
            Data_calc.RE_Vel_Y(isnan(Data_calc.RE_Vel_Z)) = NaN;             
        else %All other data types
            Data_calc.LE_Vel_X = spline_filt(ts,irlssmooth(Data_cal.LE_Vel_X,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_X = spline_filt(ts,irlssmooth(Data_cal.RE_Vel_X,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.LE_Vel_Y = spline_filt(ts,irlssmooth(Data_cal.LE_Vel_Y,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_Y = spline_filt(ts,irlssmooth(Data_cal.RE_Vel_Y,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.LE_Vel_Z = spline_filt(ts,irlssmooth(Data_cal.LE_Vel_Z,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_Z = spline_filt(ts,irlssmooth(Data_cal.RE_Vel_Z,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.LE_Vel_LARP = spline_filt(ts,irlssmooth(Data_cal.LE_Vel_LARP,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_LARP = spline_filt(ts,irlssmooth(Data_cal.RE_Vel_LARP,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.LE_Vel_RALP = spline_filt(ts,irlssmooth(Data_cal.LE_Vel_RALP,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_RALP = spline_filt(ts,irlssmooth(Data_cal.RE_Vel_RALP,filt.vel.irlssmooth),ts,filt.vel.spline);       
        end 
    end
end