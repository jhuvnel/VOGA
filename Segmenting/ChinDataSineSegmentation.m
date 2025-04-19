%% Chin Data Sine Segmentation
%This function takes in the mat files created from EMA, segments them into
%different frequencies and then alligns them for further processing.
%It takes in the arguments of chin ID, a struct called SineData, a boolean
%that allows the user to save each of the segments of the file, and a 
%boolean that decides whether to preview the data during segmentation.
%It outputs a struct SegDat that contains the segmented data.
%The input file must have a struct called SineData with fields Time, Chair,
%Lz, and Rz.

%Written by Celia Fernandez Brillet based on code from Andrianna Ayiotis

function ChinDataSineSegmentation(date, rawfile, chin,in_fname,SineData,save_data)
disp(chin)
disp(in_fname)
%% Get data from file
% First get rid of rows where all traces have nan values (oversampled data
% from EMA)
theta_L_cam = SineData(1,1);
theta_R_cam = SineData(1,2);
% idx_missing=find(isnan(SineData(6,:))&isnan(SineData(7,:))&isnan(SineData(8,:))&isnan(SineData(9,:))&isnan(SineData(10,:))&isnan(SineData(11,:)));
% SineData(:,idx_missing)=[];
% Now save different arrays from table and (linearly) interpolate across missing data
% time_L = SineData(2,:);
% time = time - time(1);
time = SineData(2,:);
Lx = SineData(6,:);
% Lx_interp = interp1(time(~isnan(Lx)),Lx(~isnan(Lx)),time);
Ly = SineData(7,:);
% Ly_interp = interp1(time(~isnan(Ly)),Ly(~isnan(Ly)),time);
Lz = SineData(8,:);
% Lz_interp = interp1(time(~isnan(Lz)),Lz(~isnan(Lz)),time);
Rx = SineData(9,:);
% Rx_interp = interp1(time(~isnan(Rx)),Rx(~isnan(Rx)),time);
Ry = SineData(10,:);
% Ry_interp = interp1(time(~isnan(Ry)),Ry(~isnan(Ry)),time);
Rz = SineData(11,:);
% Rz_interp = interp1(time(~isnan(Rz)),Rz(~isnan(Rz)),time);
if contains(lower(in_fname),'dps')
    original_stim = SineData(4,:);
elseif contains(lower(in_fname),'ua')
    original_stim = 13.5*SineData(5,:);
else
    disp('Check file name or data type. Should either contain mechanical motion or electrical stimulation data.');
end
stim = round(original_stim); % doing this for easier segmenting
maxamp = max(stim);
% Now convert from rotation to Fick -- only necessary when exporting from
% EMA in rotation vector rather than Fick
% fick_L = rot2fick([Lz_interp Ly_interp Lx_interp]);
% Lx = fick_L(:,1);
% Ly = fick_L(:,2);
% Lz = fick_L(:,3);
% fick_R = rot2fick([Rz_interp Ry_interp Rx_interp]);
% Rx = fick_R(:,1);
% Ry = fick_R(:,2);
% Rz = fick_R(:,3);
% Lx = Lx_interp;
% Ly = Ly_interp;
% Lz = Lz_interp;
% Rx = Rx_interp;
% Ry = Ry_interp;
% Rz = Rz_interp;
%% Segment into different frequencies using the chair velocity
l = length(time);
% This accounts for whether there is zero padding at the beginning of the
% data or not
if stim(1) == 0
    starts = 1;
else
   starts = []; 
end
stops = [];
%Find the string of 0s in the array
for i = 2:l
    if (stim(i) == 0 && stim(i-1) ~= 0) %Beginning of a sequence of 0s
        starts = [starts;i];
    elseif (stim(i) ~= 0 && stim(i-1) == 0) %End of a sequence of 0s
        stops = [stops;i-1];
    end
end
%%Accounts for 0 padding at the end of an array if needed
if length(starts) == length(stops)
    res = [starts,stops,stops-starts+1];
elseif length(starts) == length(stops)+1
    stops = [stops;l];
    res = [starts,stops,stops-starts+1];
else
    disp('Discrepancy between number of start and stop points.')
end
%Find only the places that have more than 100 consecutive 0s
breaks = res((res(:,3)>50),:);
blen = size(breaks,1);
%From the breaks, I can find the segments of data
%If no 0 padding at start
if breaks(1,1)-50 > 0
    segs(1,1) = 1;
    segs(1,2) = breaks(1,1)-1;
else
    segs = [];
end
for i = 1:blen-1
    if length([breaks(i,2)+1:breaks(i+1,1)-1])>50
        segs = [segs; breaks(i,2)+1,breaks(i+1,1)-1];
    end
end
%If no 0 padding at end
if l - breaks(blen,2) > 100
    segs = [segs; breaks(blen,2)+1,l];
end

%Make sure the cycle starts on a 0
segs(:,1) = segs(:,1)-1;
if(segs(1,1)==0)
   segs(1,1) = 1; 
end
p2 = plot(time,Lx,'r');
hold on
p3 = plot(time,Rx,'m');
p4 = plot(time,Ly,'Color',"#77AC30");
p5 = plot(time,Ry,'g');
p6 = plot(time,Lz,'Color',"#0072BD");
p7 = plot(time,Rz,'c');
p1 = plot(time,stim,'k','LineWidth',3);
hold off
close
p1=plot(time,stim,'k','LineWidth',3);
for i = 1:size(segs,1)
    h1 = line([time(segs(i,1)) time(segs(i,1))],[-1.1*maxamp 1.1*maxamp]);
    h2 = line([time(segs(i,2)) time(segs(i,2))],[-1.1*maxamp 1.1*maxamp]);
    % Set properties of lines
    set([h1 h2],'Color','k','LineWidth',2)
    % Add a patch
    patch([time(segs(i,1)) time(segs(i,2)) time(segs(i,2)) time(segs(i,1))],1000*[-1.1*maxamp -1.1*maxamp 1.1*maxamp 1.1*maxamp],'g')
end
set(gca,'children',flipud(get(gca,'children')))
title('Check segments of data in the file')
legend(p1,{'Stim'})
xlabel('Time (s)')
ylabel('Deg/s')
axis([0 time(end) -1.1*maxamp 1.1*maxamp])
pause;
% plot(time(segs(2,1):segs(2,2)),15*sin(2*pi*1*(time(segs(2,1):segs(2,2))-time(segs(2,1)))),'r','LineWidth',3)
%% Analyze the data separately and save each aligned segment
slen = size(segs,1);
%Ensure t vector is rounded correctly (artifact of data input)
dt = mean(diff(time));
% time = (0:dt:(length(time)-1)*dt)';
for i = 1:slen
    t = time(segs(i,1):segs(i,2));
    cdat = stim(segs(i,1):segs(i,2));
    lxdat = Lx(segs(i,1):segs(i,2));
    lydat = Ly(segs(i,1):segs(i,2));
    lzdat = Lz(segs(i,1):segs(i,2));
    rxdat = Rx(segs(i,1):segs(i,2));
    rydat = Ry(segs(i,1):segs(i,2));
    rzdat = Rz(segs(i,1):segs(i,2));
    %%%%%%%%%%%%%%%%% Uncomment the following and Fs=1/0.001; if I want to
    %%%%%%%%%%%%%%%%% oversample again
    % Interpolate data points to allign
%     tt_c = t(1):0.001:t(end);
%     cdats_c = spline(t,cdat,tt_c);
%     lxdats_c = spline(t,lxdat,tt_c);
%     lydats_c = spline(t,lydat,tt_c);
%     lzdats_c = spline(t,lzdat,tt_c);
%     rxdats_c = spline(t,rxdat,tt_c);
%     rydats_c = spline(t,rydat,tt_c);
%     rzdats_c = spline(t,rzdat,tt_c);
    Fs = 1/dt;
    tt_c = t;
    cdats_c = cdat;
    lxdats_c = lxdat;
    lydats_c = lydat;
    lzdats_c = lzdat;
    rxdats_c = rxdat;
    rydats_c = rydat;
    rzdats_c = rzdat;
    %% Find the native frequency of the chair velocity
%     Fs=1/0.001;
    N=length(cdats_c);
    fx=fft(cdats_c)/N;
    f=(0:N-1)*Fs/N;
    x = f(1:floor(end/2));
    y = 10*log10(abs(fx(1:floor(end/2))));
    [~,ind] = max(y);
    freq = x(ind);
    if(freq > 0.9) %This should get frequencies 10, 5, 2, and 1 Hz.
        round_freq = round(freq);
    elseif (freq>0.08&&freq<=0.9) %This should get 0.5, 0.2, and 0.1 Hz.
        round_freq = round(freq,1);
    else %This should get 0.05 and 0.02 Hz.
        round_freq = round(freq,2);
    end
    if contains(lower(in_fname),'dps')
        cdat = -maxamp*sin(2*pi*round_freq*(time(segs(i,1):segs(i,2))-time(segs(i,1))));
    elseif contains(lower(in_fname),'ua')
        cdat = maxamp*sin(2*pi*round_freq*(time(segs(i,1):segs(i,2))-time(segs(i,1))));
    end
    cdats_c = spline(t,cdat,tt_c);
    cdats_c = spline(t,original_stim(segs(i,1):segs(i,2)),tt_c);
     %% Catch multiple versions (up to 3)
    if contains(lower(in_fname),'dps')
        expt_unit = 'dps';
    elseif contains(lower(in_fname),'ua')
        expt_unit = 'uA';
    end
    in_fname_elements = split(string(in_fname(1:end-4)),'_');
    if sum(contains(in_fname_elements,'SINE'))
        expt_type = 'SINE';
        if sum(contains(in_fname_elements,'REDO')) % should improve coding this so it's not hardcoded and can accept more arguments
            expt_type = 'SINE-REDO';
        end
        out_file = [date,'-',chin,'-Moogles-',expt_type,'-',num2str(maxamp),expt_unit,'-',strrep(num2str(round_freq),'.','p'),'HzSegmented.mat'];
        if isfile(out_file)
            out_file = [date,'-',chin,'-Moogles-',expt_type,'-',num2str(maxamp),expt_unit,'-',strrep(num2str(round_freq),'.','p'),'HzSegmented_v2.mat'];
            if isfile(out_file)
                out_file = [date,'-',chin,'-Moogles-',expt_type,'-',num2str(maxamp),expt_unit,'-',strrep(num2str(round_freq),'.','p'),'HzSegmented_v3.mat'];
                disp('Now three versions of this file exist. More versions will write over v3.')
            end
        end
        info.dataType = ['Sine-DC-LARP-',num2str(round_freq),'Hz-',num2str(maxamp),expt_unit];
    elseif sum(contains(in_fname_elements,'STEP'))
        out_file = [date,'-',chin,'-Moogles-',char(in_fname_elements(end)),'-',num2str(maxamp),expt_unit,'-','Segmented.mat'];
        if isfile(out_file)
            out_file = [date,'-',chin,'-Moogles-',char(in_fname_elements(end)),'-',num2str(maxamp),expt_unit,'-','Segmented_v2.mat'];
            if isfile(out_file)
                out_file = [date,'-',chin,'-Moogles-',char(in_fname_elements(end)),'-',num2str(maxamp),expt_unit,'-','Segmented_v3.mat'];
                disp('Now three versions of this file exist. More versions will write over v3.')
            end
        end
        info.dataType = ['Step-DC-LARP-',num2str(maxamp),'uA'];
    end
    %Make struct to save
    info.rawfile = rawfile;
    info.rawnotes = rawfile;
    info.subject = chin;
    info.ear = 'L';
    info.visit = 'NA';
    info.exp_date = date;
    info.goggle_ver = 'Moogles';
    info.goggle_reorient_ang = 0;
    info.TriggerShift = 0;
    %info.dataType = ['Sine-DC-LARP-',num2str(round_freq),'Hz-',num2str(maxamp),'dps'];
    info.stim_axis = [-1,1,0];
    info.Name = out_file;
    info.theta_L_cam = theta_L_cam;
    info.theta_R_cam = theta_R_cam;
%     info.round_freq = round_freq;
%     info.true_freq = freq;
%     info.maxvel = maxamp;
%     info.fname = out_file;
    Data.Fs = Fs; 
    Data.Time_Eye = tt_c';
    Data.Time_Stim = tt_c';
    Data.raw_start_t = nan;
    Data.raw_end_t = nan;
    Data.Trigger = cdats_c;
    Data.LE_Position_X = lxdats_c';
    Data.RE_Position_X = rxdats_c';
    Data.LE_Position_Y = lydats_c';
    Data.RE_Position_Y = rydats_c';
    Data.LE_Position_Z = lzdats_c';
    Data.RE_Position_Z = rzdats_c';
    Data.HeadVel_X = nan;
    Data.HeadVel_Y = nan;
    Data.HeadVel_Z = nan;
    Data.HeadAccel_X = nan;
    Data.HeadAccel_Y = nan;
    Data.HeadAccel_Z = nan;
    Data.rawfile = rawfile;
    Data.info = info;   
    
    if(save_data)
       save(out_file,'Data') 
    end
end
cd ../
end
