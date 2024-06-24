%% Segment
%Make a segmenting pipeline
%Parameters assume all LDVOG data was recorded AFTER 2016-09-26, consistent
%with recorded dates for the MVI trial
%Works for LDVOG, NKI, GNO, and ESC. Limited workability for Moog
function Segment(In_Path,Seg_Path)
%% Find the Notes File
slash = find(In_Path == filesep,1,'last');
dot = find(In_Path == '.',1,'last');
rawfile = In_Path(slash+1:dot-1);
notesfile = extractfield(dir([In_Path(1:min([strfind(In_Path,'_Updated')-1,dot-1])),'*-Notes.txt']),'name');
if isempty(notesfile)
    disp(['No notes files were found for this file: ',In_Path])
    return;
elseif length(notesfile)>1
    disp(['Too many notes files were found for this file: ',In_Path])
    return;
end
fileinfo = table2cell(readtable([In_Path(1:slash),notesfile{:}],'ReadVariableNames',false,'Format','%s %s'));
if isempty(fileinfo{7,2})
    fileinfo = table2cell(readtable([In_Path(1:slash),notesfile{:}],'ReadVariableNames',false,'Delimiter',' '));
end
%% Parse Notes File
%Change this to reflect new notes file types when they are created
info.rawfile = In_Path;
info.rawnotes = [In_Path(1:slash),notesfile{:}];
info.subject = fileinfo{1,end};
info.ear = fileinfo{2,end};
info.visit = strrep(fileinfo{3,end},' ','');
info.exp_date = fileinfo{4,end};
info.goggle_ver = fileinfo{5,end}; %This should say NKI or not
info.goggle_reorient_ang = str2double(fileinfo{6,end});
if length(info.exp_date)==14 %Just missing '-' in between date and time
    info.exp_date = [info.exp_date(1:8),'-',info.exp_date(9:end)];
elseif length(info.exp_date)==8 %Just date, add time
    underscore = strfind(rawfile,'_');
    if contains(info.goggle_ver,{'NKI','NL'})
        VOG_time = char(datetime(rawfile(underscore(4)+1:underscore(7)-1),'InputFormat','hh_mm_a'),'HHmmSS');
    elseif contains(info.goggle_ver,'LDVOG')
        if ~isempty(underscore)
            rawfile = rawfile(1:underscore(1));
        end
        rawparts = split(strrep(rawfile,'_','-'),'-');
        VOG_time = rawparts{3};
    elseif contains(info.goggle_ver,'GNO')
        VOG_time = strrep(rawfile(underscore(5)+1:underscore(8)-1),'_','');
    else
        VOG_time = '000000';
    end
    info.exp_date = [info.exp_date,'-',VOG_time];
end
%% Load Data Types In
if contains(info.goggle_ver,{'NKI','NL'})
    info.TriggerShift = 0; %measure and change
    %Suppress the warning that one of the columns is not a proper column name for a MATLAB table so it got renamed
    warning('off')
    data = readtable(In_Path,'ReadVariableNames',true); %.dat file with uncalibrated eye data
    data2 = table2array(readtable(strrep(info.rawnotes,'-Notes.txt','.csv'))); %.csv with calibrated eye data
    warning('on')
    len = min([size(data,1),size(data2,1)]);
    data = data(1:len,:);
    data2 = data2(1:len,:);
    % Data Index
    HLeftIndex = 2;
    VLeftIndex = 4;
    TLeftIndex = 10;
    HRightIndex = 3;
    VRightIndex = 12;
    TRightIndex = 13;
    Time_Eye = data2(:,1);
    Time_Stim = data2(:,1);
    Stim = zeros(length(Time_Eye),1);
    Stim(data.EventCode ~= 0) = 1;
    % Load Gyroscope and accelerometer readings
    XAxisVelHead = data.GyroY - median(data.GyroY);
    YAxisVelHead = -(data.GyroX - median(data.GyroX));
    ZAxisVelHead = -(data.GyroZ - median(data.GyroZ));
    XAxisAccelHead = data.AccelY - median(data.AccelY);
    YAxisAccelHead = -(data.AccelX - median(data.AccelX));
    ZAxisAccelHead = -(data.AccelZ - median(data.AccelZ));
    Fs = 1/median(abs(diff(Time_Eye)));
    % Load raw eye position data in Fick coordinates [degrees] but
    % adjust for the reverse of the X/Y/Z axes direction. This is validated by
    % experiments done on the aHIT in light on normal subjects:
    % https://docs.google.com/document/d/1euk4V0fbnbhg_paUZOaO32s_FEhStmtCbrHJV6WQLu8/edit
    Horizontal_LE_Position = -data2(:,HLeftIndex);
    Vertical_LE_Position = -data2(:,VLeftIndex);
    Torsion_LE_Position = data2(:,TLeftIndex);
    Horizontal_RE_Position = -data2(:,HRightIndex);
    Vertical_RE_Position = -data2(:,VRightIndex);
    Torsion_RE_Position = data2(:,TRightIndex);
    %Remove repeated values in torsion
    [gr,gc] = groupcounts(diff(find(diff([Torsion_LE_Position;Torsion_RE_Position])~=0&...
        ~isnan([Torsion_LE_Position(2:end);Torsion_RE_Position]))));
    reps = gc(gr==max(gr));
    if reps>1
        M = [movmean(Torsion_LE_Position,reps);zeros(floor(reps/2)-1,1)];
        Torsion_LE_Position = [M(round(reps/2):end);zeros(round(reps/2)-1,1)];
        M = [movmean(Torsion_RE_Position,reps);zeros(floor(reps/2)-1,1)];
        Torsion_RE_Position = [M(round(reps/2):end);zeros(round(reps/2)-1,1)];        
    end    
    GyroX = medfilt1(reshape(XAxisVelHead,[],1),3);
    GyroY = medfilt1(reshape(YAxisVelHead,[],1),3);
    GyroZ = medfilt1(reshape(ZAxisVelHead,[],1),3);
    GyroLARP = (GyroX - GyroY)/sqrt(2);
    GyroRALP = (GyroX + GyroY)/sqrt(2);
elseif contains(info.goggle_ver,'LDVOG')
    % Load data, determine VOG type
    data1 = readtable(In_Path);
    data = data1{:,contains(varfun(@class,data1,'OutputFormat','cell'),'double')};
    % Generate Time_Eye vector
    Time = data(:,2);
    % The time vector recorded by the VOG goggles resets after it reaches a value of 128.
    Time_Eye = (0:length(Time)-1)*median(diff(Time));
    % Generate the time vector for the MPU9250 Data
    Head_Sensor_Latency = 0.047; % From Mehdi Rahman bench tests, the data acquisition of the MPU9250 leads the LD VOG Goggles by 47ms
    Time_Stim = Time_Eye - Head_Sensor_Latency;
    % Data Index
    HLeftIndex = 40;
    VLeftIndex = 41;
    TLeftIndex = 42;
    HRightIndex = 43;
    VRightIndex = 44;
    TRightIndex = 45;
    XvelHeadIndex = 30;
    YvelHeadIndex = 29;
    ZvelHeadIndex = 28;
    XaccelHeadIndex = 27;
    YaccelHeadIndex = 26;
    ZaccelHeadIndex = 25;
    StimIndex = 35;
    % Around 2018-04, PJB noticed that the LD VOG system appeared to record the
    % PCU trigger LATE relative to the collected eye movement data. PJB, MR,
    % and NV performed some experiments outlined here:
    % https://docs.google.com/document/d/16EppbHOlsSabjR61xbDi798p7ldrqqTV6tOXjH1KZbI/edit#heading=h.xh2azhw75tvn
    % These experiments resulted in the following trigger recording delay in
    % samples (where positive values indicate the trigger was recorded AFTER
    % the eye movement).
    switch info.goggle_ver
        case {'1','LDVOG1'}
            TriggerDelay = 3;
            XvelHeadOffset = 2.7084;
            YvelHeadOffset = -0.5595;
            ZvelHeadOffset = 0.7228;
        case {'2','LDVOG2'}
            TriggerDelay = 5;
            XvelHeadOffset = 2.3185;
            YvelHeadOffset = -1.5181;
            ZvelHeadOffset = -1.0424;
        case {'3','LDVOG3'}
            TriggerDelay = 2;
            XvelHeadOffset = 1.9796;
            YvelHeadOffset = 0.1524;
            ZvelHeadOffset = -0.5258;
        otherwise
            TriggerDelay = 0;
            XvelHeadOffset = 0;
            YvelHeadOffset = 0;
            ZvelHeadOffset = 0;
    end
    % Load raw eye position data in Fick coordinates [degrees]
    Horizontal_LE_Position = data(:,HLeftIndex);
    Vertical_LE_Position = data(:,VLeftIndex);
    Torsion_LE_Position = data(:,TLeftIndex);
    Horizontal_RE_Position = data(:,HRightIndex);
    Vertical_RE_Position = data(:,VRightIndex);
    Torsion_RE_Position = data(:,TRightIndex);
    % Index for the VOG GPIO line
    Stim = data(1:length(Time_Eye),StimIndex);
    Stim = [Stim(TriggerDelay + 1:end) ; Stim(end)*ones(TriggerDelay,1)];
    info.TriggerShift = TriggerDelay;
    % We need to correct each gyroscope signal by subtracting the
    % correct device-specifc MPU9250 gyroscope offset. Each offset for
    % each VOG goggle ID was measured by Mehdi Rahman and posted on the
    % Google Doc located here: https://docs.google.com/a/labyrinthdevices.com/document/d/1UlZpovNkwer608aswJWdkLhF0gF-frajAdu1qgMJt9Y/edit?usp=sharing
    XvelHeadRaw = data(1:length(Time_Eye),XvelHeadIndex) + XvelHeadOffset;
    YvelHeadRaw = data(1:length(Time_Eye),YvelHeadIndex) + YvelHeadOffset;
    ZvelHeadRaw = data(1:length(Time_Eye),ZvelHeadIndex) + ZvelHeadOffset;
    XaccelHeadRaw = data(1:length(Time_Eye),XaccelHeadIndex);
    YaccelHeadRaw = data(1:length(Time_Eye),YaccelHeadIndex);
    ZaccelHeadRaw = data(1:length(Time_Eye),ZaccelHeadIndex);
    % Based on the orientation of the MPU9250 module relative to the
    % patients head, a passive (coordinate system) -150deg rotation is
    % necesary to align the coordinate system of MPU to the coordinate
    % system of the patient's head with the VOG goggles on.
    % When the patient is tested with the Automatic Head Impulse Test
    % (aHIT) device, the individual's head is pitched down 20deg to
    % (roughly) align the patient's +LHRH SCC axis with the +Z 'world'
    % axis that the aHIT rotates about. This means we need to apply a
    % (-150deg + -20deg) to orient the +Z axis of the MPU9250 seated in
    % the VOG goggles with the +LHRH axis of SCC coordinate system.
    phi = info.goggle_reorient_ang;
    Rotation_Head = [cosd(phi),0,sind(phi);0,1,0;-sind(phi),0,cosd(phi)];
    % NOTE: We are transposing the rotation matrix in order to apply a
    % PASSIVE (i.e., a coordinate system) transformation
    A = Rotation_Head' * [XvelHeadRaw' ; YvelHeadRaw' ; ZvelHeadRaw'];
    XAxisVelHead = A(1,:);
    YAxisVelHead = A(2,:);
    ZAxisVelHead = A(3,:);
    B = Rotation_Head' * [XaccelHeadRaw' ; YaccelHeadRaw' ; ZaccelHeadRaw'];
    XAxisAccelHead = B(1,:);
    YAxisAccelHead = B(2,:);
    ZAxisAccelHead = B(3,:);
    temp_SampPer = median(abs(diff(Time)));
    temp_Fs = 1/temp_SampPer;
    expected_Fs = [100,150,250]; %these are the expected computer frame rates we've seen for the LDVOG system
    comp = abs(expected_Fs - temp_Fs) < expected_Fs*0.1;
    if length(find(comp)) ~= 1
        answer = questdlg(['The computed sample rate produced an unexpected value.' newline 'Expected Sample Rate: ' num2str(expected_Fs) 'Hz' newline ...
            'Mean Computed Sample Rate: ' num2str(1/abs(mean(diff(Time)))) 'Hz' newline 'Median Computed Sample Rate: ' num2str(temp_Fs) 'Hz' newline ...
            'How would you like to proceed?'], ...
            'Sample Rate Error', ...
            'Let me define the sample rate',['Force the expected rate:' num2str(expected_Fs) 'Hz'],['Force the expected rate:' num2str(expected_Fs) 'Hz']);
        % Handle response
        switch answer
            case 'Let me define the sample rate'
                prompt = {'Enter Sample Rate:'};
                title = 'Input';
                dims = [1 35];
                definput = {'150','hsv'};
                answer = str2double(inputdlg(prompt,title,dims,definput));
                Fs = answer;
            case ['Force the expected rate:' num2str(expected_Fs) 'Hz']
                Fs = expected_Fs;
        end
    else
        Fs = expected_Fs(comp==1);
    end
    GyroX = reshape(XAxisVelHead,[],1);
    GyroY = reshape(YAxisVelHead,[],1);
    GyroZ = reshape(ZAxisVelHead,[],1);
    GyroLARP = (GyroX - GyroY)/sqrt(2);
    GyroRALP = (GyroX + GyroY)/sqrt(2);
elseif contains(info.goggle_ver,'GNO') % GNO File
    %Load the values
    data = table2array(readtable(In_Path));
    Time_Eye = (data(:,1) - data(1,1))/10e6;
    Time_Stim = Time_Eye;
    GyroZ = data(:,4);
    GyroLARP = -data(:,3); %Invert to make it aligned with +L canal, wasn't negative here until 01/10/22
    GyroRALP = data(:,2);
    Fs = round(1/median(diff(Time_Eye)));    
    Horizontal_RE_Velocity = data(:,5);
    Vertical_RE_Velocity = data(:,6);  
    info.TriggerShift = 0;
    % Find the accepted traces
    if isfile(strrep(In_Path,'.txt','.csv'))  
        opts = detectImportOptions(strrep(In_Path,'.txt','.csv'),'ReadVariableNames',false);
        GNO_CSV = table2cell(readtable(strrep(In_Path,'.txt','.csv'),setvartype(opts,'char')));
        if size(GNO_CSV,2)>1
            GNO_CSV = join(GNO_CSV,',');
        end
        imp_ind = find(contains(GNO_CSV,{'Direction,L','Direction,R'}));
        i_end = [imp_ind(2:end)-1;length(GNO_CSV)];
        imp_dat = cell(length(imp_ind),2);
        for i = 1:length(imp_ind)
            head_ind = find(contains(GNO_CSV(imp_ind(i):i_end(i)),'Head'));
            if length(head_ind)==1
                imp_dat(i,1) = extract(GNO_CSV{imp_ind(i)},"L"|"R");
                vec = str2double(split(strrep(GNO_CSV(imp_ind(i)-1+head_ind),'Head',''),','));
                imp_dat{i,2} = vec(~isnan(vec));
            end
        end
        if contains(In_Path,'Lateral')
            DetectedTraces_HeadVel = [horzcat(imp_dat{contains(imp_dat(:,1),'L'),2}),-horzcat(imp_dat{contains(imp_dat(:,1),'R'),2})];
        elseif contains(In_Path,{'LARP','RALP'})
            DetectedTraces_HeadVel = [-horzcat(imp_dat{contains(imp_dat(:,1),'L'),2}),horzcat(imp_dat{contains(imp_dat(:,1),'R'),2})];
        end
    else
        GNO_CSV = [];
        DetectedTraces_HeadVel=[];
    end
    %Put xml data in the segment too
    if isfile(strrep(In_Path,'.txt','.xml'))
        GNO_XML = cellstr(readlines(strrep(In_Path,'.txt','.xml')));
    else
        GNO_XML = [];
    end
elseif contains(info.goggle_ver,{'ESC1','ESC2'}) % ESC File with .mat
    data1 = load(In_Path);
    if isfield(data1,'content_RAW')
        phi = -20;
        data = data1.content_RAW.Data;
        if isfield(data1.content_CAL,'R')
            eye_cal_gyro = ([cosd(phi),0,sind(phi);0,1,0;-sind(phi),0,cosd(phi)]'*data1.content_CAL.R);
            eye_cal_accel = ([cosd(phi),0,sind(phi);0,1,0;-sind(phi),0,cosd(phi)]'*data1.content_CAL.R);
        else %Missing calibration 
            eye_cal_gyro = eye(3);
            eye_cal_accel = eye(3);
        end
        labs = cellstr(data1.content_RAW.DataNames);
    elseif isfield(data1,'Data')
        data = data1.Data;
        labs = cellstr(data1.DataNames);
        eye_cal_gyro = eye(3);
        eye_cal_accel = eye(3);   
    end 
    EyeTimeIndex = find(contains(labs,'SystemTime'));
    if ~isempty(EyeTimeIndex) % Most cases
        old_Time_Eye = data(:,EyeTimeIndex) - data(1,EyeTimeIndex);
        keep_eye = [true;diff(old_Time_Eye)>0];
        Fs = 1/median(diff(old_Time_Eye(keep_eye))); 
        Time_Eye = 0:(1/Fs):old_Time_Eye(end);
        Time_Stim = Time_Eye;
        warning('off')
        Horizontal_LE_Velocity = spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'EyeVelZ')),Time_Eye);
        Vertical_LE_Velocity = spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'EyeVelY')),Time_Eye);
        Torsion_LE_Velocity = spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'EyeVelX')),Time_Eye); 
        GyroX = spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'InertialVelX')),Time_Eye); 
        GyroY = spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'InertialVelY')),Time_Eye); 
        GyroZ = spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'InertialVelZ')),Time_Eye); 
        A = eye_cal_gyro*[GyroX;GyroY;GyroZ];
        GyroX = A(1,:);
        GyroY = A(2,:);
        GyroZ = A(3,:);
        GyroLARP = (GyroX - GyroY)/sqrt(2);
        GyroRALP = (GyroX + GyroY)/sqrt(2); 
        AccelX = spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'InertialAccelX')),Time_Eye); 
        AccelY = spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'InertialAccelY')),Time_Eye); 
        AccelZ = spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'InertialAccelZ')),Time_Eye); 
        B = eye_cal_accel*[AccelX;AccelY;AccelZ];
        AccelX = B(1,:);
        AccelY = B(2,:);
        AccelZ = B(3,:);
    elseif any(contains(labs,'RealTime'))
        %Old version missing much data
        %Time vec is in msec
        t1 = data(:,strcmp(labs,'Time')); %in msec
        nan_time = isnan(t1);
        len = size(data,1)-sum(nan_time);
        Time_Eye = (0:(len-1))*10-3; %now in sec
        Time_Stim = Time_Eye;
        Horizontal_LE_Velocity = NaN(len,1);
        Vertical_LE_Velocity = NaN(len,1);
        Torsion_LE_Velocity = NaN(len,1);
        GyroX = NaN(len,1);
        GyroY = NaN(len,1);
        GyroZ = NaN(len,1);
        GyroLARP = NaN(len,1);
        GyroRALP = NaN(len,1);
        AccelX = NaN(len,1);
        AccelY = NaN(len,1);
        AccelZ = NaN(len,1);
        if any(contains(reshape(fileinfo,[],1),'LHRH'))
            Horizontal_LE_Velocity = data(~nan_time,contains(labs,'EyeVel'));
            GyroZ = data(~nan_time,contains(labs,'HeadVel'));
        elseif any(contains(reshape(fileinfo,[],1),'LARP'))
            Vertical_LE_Velocity = data(~nan_time,contains(labs,'EyeVel'));
            GyroLARP = data(~nan_time,contains(labs,'HeadVel'));
        elseif any(contains(reshape(fileinfo,[],1),'RALP'))
            Vertical_LE_Velocity = data(~nan_time,contains(labs,'EyeVel'));
            GyroRALP = data(~nan_time,contains(labs,'HeadVel'));
        end       
    else
        disp(['Not segmented: ',In_Path])
        return;
    end
    Fs = 1/median(diff(Time_Eye)); 
    info.TriggerShift = 0;
    if isfile(strrep(strrep(In_Path,'_export',''),'.mat','.xls'))
        fname = strrep(strrep(In_Path,'_export',''),'.mat','.xls');
        movefile(fname,strrep(fname,'.xls','.csv'))
        ESC_CSV = cellstr(readlines(strrep(fname,'.xls','.csv')));
        ESC_CSV(cellfun(@isempty,ESC_CSV)) = [];
        ESC_CSV = split(ESC_CSV,',');
        movefile(strrep(fname,'.xls','.csv'),fname)
    else
        ESC_CSV = [];
    end    
    %plot(Time_Eye,GyroZ,'k-',Time_Eye,GyroY,'g-',Time_Eye,GyroX,'b-')
elseif contains(info.goggle_ver,'ESC3') % ESC File with .csv data
    ESC_CSV = cellstr(readlines(strrep(In_Path,'ImuData','Numerical')));
    ESC_CSV(cellfun(@isempty,ESC_CSV)) = [];
    max_comma = max(count(ESC_CSV,','));
    comma_cnt = count(ESC_CSV,',');
    while any(comma_cnt~=max_comma)
        ESC_CSV(comma_cnt~=max_comma) = strcat(ESC_CSV(comma_cnt~=max_comma),',');        
        comma_cnt = count(ESC_CSV,',');
    end
    ESC_CSV = split(ESC_CSV,',');
    info.TriggerShift = 0;
    data = readtable(In_Path);
    data1 = data;
    Time_Eye = data.LeftTime - data.LeftTime(1);
    Time_Stim = Time_Eye; %Update later as needed
    Stim = zeros(length(Time_Eye),1); %No external trigger yet
    % Load Gyroscope and accelerometer readings
    % Load raw eye position data in Fick coordinates [degrees] with no torsion.
    % Gyros appear to be properly aligned with +XYZ. This is shown here:
    % https://docs.google.com/document/d/1euk4V0fbnbhg_paUZOaO32s_FEhStmtCbrHJV6WQLu8/edit
    GyroX = data.HeadInertialVelX;
    GyroY = data.HeadInertialVelY;
    GyroZ  = data.HeadInertialVelZ;
    GyroLARP = (GyroX - GyroY)/sqrt(2);
    GyroRALP = (GyroX + GyroY)/sqrt(2); 
    AccelX = data.HeadInertialAccelX;
    AccelY = data.HeadInertialAccelY;
    AccelZ = data.HeadInertialAccelZ;
    Fs = 1/median(abs(diff(Time_Eye)));
    Horizontal_RE_Position = data.LeftPupilCol;
    Vertical_RE_Position = data.LeftPupilRow;
    Torsion_RE_Position = 0*Time_Eye; %NO TORSION TRACKING 
end
%% Figure out how many experiments there are
if size(fileinfo,2)==2 %New way w/ 2 columns
    stim_info = fileinfo(7:end,2);
else %old way
    if length(strsplit(fileinfo{7},' ')) > 1
        stim_info = strrep(fileinfo(8:end),'"','');
    else %Only one experiment type
        stim_info = strrep(strcat(fileinfo{7},'-',fileinfo(8:end)),'"','');
    end
end
%% Segment
if all(contains(stim_info,{'RotaryChair','aHIT','manual','Manual','trash'})) %Fit on motion trace
    stim_info = strrep(stim_info,'manual','Manual'); %in case I forgot to capitalize
    GyroAll = GyroZ+GyroLARP+GyroRALP;
    if any(contains(lower(stim_info),'velstep'))
        thresh = 50; %Adjust as needed
        thresh2 = 5; %Counts as 0.
        approx0 = find(abs(GyroAll)< thresh2);
        start1 = find(abs(GyroAll) > thresh);
        start2 = start1([true;diff(start1)>Fs*10]); %make sure it's at least 45s long
        [~,inds] = min(abs(repmat(approx0',length(start2),1) - start2),[],2);
        seg_start = approx0(inds);
        seg_start = seg_start([true;diff(seg_start)>Fs*10]); %make sure it's at least 45s long
        seg_end = [seg_start(2:end)-1;length(Time_Eye)];
        too_long = (seg_end-seg_start)>round(130*Fs);
        seg_end(too_long) = seg_start(too_long)+130*Fs; %make sure segment is not longer than 2.5 min
    elseif any(contains(stim_info,{'Sine','Impulse'}))&&length(stim_info) > 1 %multiple sine stimuli
        %Segment and allow user to select which segments are "real"
        thresh = 20; %change as needed but works for 0.01Hz - 2Hz
        abov_thresh = find(abs(GyroAll) > thresh);
        zero_cross = [1;find(GyroAll(2:end).*GyroAll(1:end-1)<0);length(GyroAll)];
        starts1 = abov_thresh([true;diff(abov_thresh)>1]);
        ends1 = abov_thresh([diff(abov_thresh)>1;true]);
        for i = 1:length(starts1)
            starts1(i) = max([zero_cross(find(zero_cross<starts1(i),1,'last'))-10,1]);
            ends1(i) = min([zero_cross(find(zero_cross>ends1(i),1,'first'))+10,length(GyroAll)]);
        end
        % Combine segments with only a few samples between them
        for i = 1:length(starts1)-1
            if starts1(i+1)-ends1(i) < 2*Fs
                starts1(i+1) = NaN;
                ends1(i) = NaN;
            end
        end
        starts1(isnan(starts1)) = [];
        ends1(isnan(ends1)) = [];
        seg_start = starts1;
        seg_end = ends1;
        if length(seg_start)~=length(seg_end)
            error('Length of start and stop are unequal. Please manually segment')
        elseif length(seg_start)<length(stim_info)
            error('Less segments found than in the notes file. Please manually segment')
        elseif length(seg_start)>length(stim_info)
            small_seg = seg_end-seg_start < 2*Fs;
            if (length(seg_start)-sum(small_seg))>=length(stim_info)
                seg_start(small_seg) = [];
                seg_end(small_seg) = [];
            end
        end
    else %Only one segment
        seg_start = 1;
        seg_end = length(Time_Eye);
    end   
    seg_start = round(seg_start,0);
    seg_end = round(seg_end,0);
    if length(seg_start) ~= length(stim_info)
        keep = false(1,length(seg_start));
        plot(NaN,NaN)
        hold on
        for j = 1:length(seg_start) %Now plot all fills
            fill(Time_Eye([seg_start(j),seg_end(j),seg_end(j),seg_start(j)])',500*[1,1,-1,-1]',0.85*[1,1,1]);
        end
        plot(Time_Eye,GyroLARP,'k:',Time_Eye,GyroRALP,'k--',Time_Eye,GyroZ,'k-')
        hold off
        uiwait(msgbox('Click on all valid segments on the figure.'))
        for i = 1:length(stim_info)
            [x,~] = ginput(1);
            j = find(Time_Eye(seg_start)-x<0,1,'last');
            keep(j) = true;
            hold on
            fill(Time_Eye([seg_start(j),seg_end(j),seg_end(j),seg_start(j)])',500*[1,1,-1,-1]','g');
            plot(Time_Eye,GyroLARP,'k:',Time_Eye,GyroRALP,'k--',Time_Eye,GyroZ,'k-')
            hold off
        end
        seg_start = seg_start(keep);
        seg_end = seg_end(keep);
    end
    start = round(seg_start,0);
    ends = round(seg_end,0);
    plot(NaN,NaN)
    hold on
    %Now plot all fills
    for j = 1:length(start)
        fill([Time_Eye(start(j)),Time_Eye(ends(j)),Time_Eye(ends(j)),Time_Eye(start(j))]',[500,500,-500,-500]','g');
    end
    plot(Time_Eye,GyroLARP,'k:',Time_Eye,GyroRALP,'k--',Time_Eye,GyroZ,'k-')
    hold off
    pause(1)    
    if any(contains(stim_info,{'Impulse','Gaussian'}))
        %Add a stim_info entry for each canal (LHRH -> LH and RH)
        rep_ind = sort([(1:length(start))';find(contains(stim_info,{'Impulse','Gaussian'}))]);
        start = start(rep_ind);
        ends = ends(rep_ind); 
        stim_info2 = repmat(stim_info,1,2);
        stim_info2(~contains(stim_info2(:,2),{'Impulse','Gaussian'}),2) = {''};
        %Make the canals one-sided here
        stim_info2(contains(stim_info2(:,1),{'Impulse','Gaussian'}),1) = strrep(strrep(strrep(stim_info2(contains(stim_info2(:,1),{'Impulse','Gaussian'}),1),'LHRH','LH'),'LARP','RP'),'RALP','RA');
        stim_info2(contains(stim_info2(:,1),{'Impulse','Gaussian'}),2) = strrep(strrep(strrep(stim_info2(contains(stim_info2(:,1),{'Impulse','Gaussian'}),2),'LHRH','RH'),'LARP','LA'),'RALP','LP');
        stim_info = reshape(stim_info2',[],1);
        stim_info(cellfun(@isempty,stim_info)) = [];                        
    end 
elseif all(contains(stim_info,{'eeVOR','trash'})) %Good for all externally triggered stimuli for now
    if any(contains(stim_info,'Activation'))
        %Use this to make sure the segments don't start or end on a
        %long blink (which will mess with filtering)
        all_traces = medfilt1([Torsion_LE_Position,Torsion_RE_Position,...
            Vertical_LE_Position,Vertical_RE_Position,...
            Horizontal_LE_Position,Horizontal_RE_Position],5);
        if all(diff(Stim)==0) %Issue w/ the software trigger
            %Fit from long blinks and recorded transition times
            nan_trace = find(all(isnan(all_traces),2));
            nan_trace_s = nan_trace([true;diff(nan_trace)>1]);
            nan_trace_e = nan_trace([diff(nan_trace)>1;true]);
            %Make sure there's at least 1 second of eyes closed
            long_e = find((nan_trace_e-nan_trace_s)>1*Fs);
            Time_e = Time_Eye(nan_trace_s(long_e))/60;
            proj_t = [2,6,7,11,12,16,17]+[20,10,34,29,44,37,32]/60;
            inds = NaN*proj_t;
            for i = 1:length(inds)
                [~,inds(i)] = min(abs(Time_e - proj_t(i)));
            end
            inds2 = [1;nan_trace_e(long_e(inds));length(Stim)];
            %Write the new stim trace
            Stim = 0*Stim;
            for i = 2:2:length(inds2)
                Stim(inds2(i):inds2(i+1)) = 1;
            end
            %                    %Check plot
            plot(Time_Eye/60,all_traces,'k',Time_Eye/60,Stim,'b',Time_Eye(inds2)/60,20*ones(1,length(inds2)),'g*')
            xlabel('Time (min)')
            ylabel('Position (dps)')
            set(gca,'XMinorGrid','on','YLim',[-30,30])
        end
        if length(stim_info) > 1 %each light/dark cycle is seperate
            %Every trigger toggle is a cycle
            temp = find(abs(diff(Stim))==1);
            start = [1;temp+1];
            ends = [temp;length(Stim)];
            val_trace = find(all(~isnan(all_traces),2));
            %constrain the segments by non-nan values
            for i = 1:length(ends)
                start(i) = val_trace(find(val_trace>=start(i),1,'first'));
                ends(i) = val_trace(find(val_trace<=ends(i),1,'last'));
            end
        else %one big trace
            start = 1;
            ends = length(Time_Eye);
        end
    elseif any(contains(stim_info,'MultiVector'))
        %These files have trigger toggles for the "ramp" of every cycle
        if Stim(1)==1 %If Stim started "high", the first ramp is missing so add it back in
            temp = diff(find(abs(diff(Stim))==1));
            ramp_dur = temp(2);
            Stim(1:(find(diff(Stim)==-1,1,'first')-ramp_dur)) = 0;
        end
        temp = [find(diff(Stim)==1)-1;length(Stim)];
        temp2 = diff(temp);
        [~,ind] = sort(temp2,'descend');
        ind2 = sort(ind(1:length(stim_info)));
        start = temp([0;ind2(1:end-1)]+1)-1;
        ends = temp(ind2)+sum(temp2(ind2-(1:2)),2);
    elseif any(contains(stim_info,'VelStep'))
        %These files have trigger toggles for the "ramp" of every cycle
        if Stim(1)==1 %If Stim started "high", the first ramp is missing so add it back in
            temp = diff(find(abs(diff(Stim))==1));
            ramp_dur = temp(2);
            Stim(1:(find(diff(Stim)==-1,1,'first')-ramp_dur)) = 0;
        end
        temp = [find(diff(Stim)==1)-1;length(Stim)];
        all_starts = temp(1:2:end);
        all_ends = [temp(3:2:end);length(Stim)];
        temp2 = all_ends-all_starts;
        [~,ind] = sort(temp2,'descend');
        ind2 = sort(ind(1:length(stim_info)));
        start = all_starts([1;ind2(1:end-1)+1]);
        ends = all_ends(ind2);  
    elseif any(contains(stim_info,'Sine'))
        %Every trigger toggle is a cycle
        thresh_tol = Fs*0.4; %tolerance for differences in cyc length (always at least 400ms of break)
        temp = find(abs(diff(Stim))==1);
        temp2 = diff(diff(temp));
        all_ends = unique(temp([find(temp2<-thresh_tol);find(temp2>thresh_tol)+1;end]));
        all_starts = unique(temp([1;find(temp2<-thresh_tol)+1;find(temp2>thresh_tol)+2]));        
        start = all_starts-5;
        ends = all_ends+5;        
    else %Pulse Train and Autoscan--high period is stimulation and low is break.
        temp = find(diff(Stim)==1)-1;
        temp2 = diff(temp);
        all_starts = temp;
        all_ends = all_starts+median(temp2);
        [~,ind] = sort(temp2,'descend');
        ind2 = sort(ind(1:length(stim_info)-1));
        start = all_starts([1;ind2+1])-round(0.5*median(diff(temp)));
        ends = all_ends([ind2;end])+round(0.5*median(diff(temp)));
    end
    stim = Stim;
    start(start<0) = 1;
    ends(ends>length(Time_Stim)) = length(Time_Stim);
    start = round(start);
    ends = round(ends);
    plot(Time_Stim,stim,'k',Time_Stim(start),stim(start),'r*',Time_Stim(ends),stim(ends),'b*')
    xlabel('Time (s)')
    ylabel('Trigger')
    legend('Stimulus','Start','End')
    pause(1)
else
    disp(stim_info)
    error('Unrecognized experiment type during segmenting.')
end
%% Save
%Remove any traces marked as "trash" or "calibrations"
rm = contains(lower(stim_info),{'trash','calib'});
stim_info(rm) = [];
start(rm) = [];
ends(rm) = [];
if ~isempty(stim_info)
    if length(start) ~= length(ends)
        disp('Something went awry when detecting segments')
        disp('There are unequal numbers of starting and ending points. Please adjust or segment the following file manually:')
        disp(In_Path)
        savefig([In_Path(1:end-4),'SegERR.fig'])
    elseif length(stim_info) > length(start)
        disp('There are more experiments in the notes file than there are detected segments. Please adjust the notes for the following file:')
        disp(In_Path)
        savefig([In_Path(1:end-4),'SegERR.fig'])
    elseif length(stim_info) < length(start)
        disp('There are more detected segments than experiments in the notes file. Please adjust the notes for the following file:')
        disp(In_Path)
        savefig([In_Path(1:end-4),'SegERR.fig'])
    else %Can automatically segment
        for i = 1:length(start)
            info.dataType = stim_info{i};
            %Fix this to be more general for pulse trains
            if contains(stim_info{i},'65Vector') %stim vec in the info file
                info.stim_axis = str2double(info.dataType(strfind(info.dataType,'['):strfind(info.dataType,']')));
            elseif contains(stim_info{i},{'Impulse','Gaussian'})
                %Add code here
                if contains(stim_info{i},'LH')
                    info.stim_axis = [0,0,1];
                elseif contains(stim_info{i},'RH')
                    info.stim_axis = [0,0,-1];
                elseif contains(stim_info{i},'RP')
                    info.stim_axis = [1,0,0];
                elseif contains(stim_info{i},'LA')
                    info.stim_axis = [-1,0,0];
                elseif contains(stim_info{i},'RA')
                    info.stim_axis = [0,1,0];
                elseif contains(stim_info{i},'LP')
                    info.stim_axis = [0,-1,0];
                else
                    info.stim_axis = [0,0,0];
                end
            elseif contains(stim_info{i},'RotaryChair') %horizontal only
                info.stim_axis = [0,0,1]; %Moves to the left first
            elseif contains(stim_info{i},{'PulseTrain','Autoscan'})
                if contains(stim_info{i},{'RALP','LP','RA'})
                    if strcmp(info.ear,'L')
                        info.stim_axis = [0,-1,0];
                    elseif strcmp(info.ear,'R')
                        info.stim_axis = [0,1,0];
                    else
                        info.stim_axis = [0,0,0];
                    end
                elseif contains(stim_info{i},{'LARP','LA','RP'})
                    if strcmp(info.ear,'L')
                        info.stim_axis = [-1,0,0];
                    elseif strcmp(info.ear,'R')
                        info.stim_axis = [1,0,0];
                    else
                        info.stim_axis = [0,0,0];
                    end
                elseif contains(stim_info{i},{'LHRH','LH','RH'})
                    if strcmp(info.ear,'L')
                        info.stim_axis = [0,0,1];
                    elseif strcmp(info.ear,'R')
                        info.stim_axis = [0,0,-1];
                    else
                        info.stim_axis = [0,0,0];
                    end
                else
                    info.stim_axis = [0,0,0];
                end
            elseif contains(stim_info{i},'Sine')
                if contains(stim_info{i},'RALP')
                    info.stim_axis = [0,1,0];
                elseif contains(stim_info{i},'LHRH')
                    info.stim_axis = [0,0,1];
                elseif contains(stim_info{i},'LARP')
                    info.stim_axis = [1,0,0];
                elseif contains(stim_info{i},'X')
                    info.stim_axis = [0.707,0.707,0];
                elseif contains(stim_info{i},'Y')
                    info.stim_axis = [-0.707,0.707,0];
                else
                    info.stim_axis = [0,0,0];
                end    
            else
                info.stim_axis = [0,0,0];
            end
            i1 = start(i);
            i2 = ends(i);
            Data.info = info;
            Data.Fs = Fs;
            Data.Time_Eye = Time_Eye(i1:i2);
            Data.Time_Stim = Time_Stim(i1:i2);
            Data.raw_start_t = Time_Stim(i1);
            Data.raw_end_t = Time_Stim(i2);
            if contains(info.goggle_ver,'GNO')
                Data.RE_Vel_Y = Vertical_RE_Velocity(i1:i2);
                Data.RE_Vel_Z = Horizontal_RE_Velocity(i1:i2);
                Data.HeadVel_Z = GyroZ(i1:i2);
                Data.HeadVel_L = GyroLARP(i1:i2);
                Data.HeadVel_R = GyroRALP(i1:i2);
                Data.DetectedTraces_HeadVel = DetectedTraces_HeadVel;
                Data.CSVData = GNO_CSV;
                Data.XMLData = GNO_XML;
            elseif contains(info.goggle_ver,'ESC')
                if contains(info.goggle_ver,{'ESC3'})
                    Data.RE_Position_Y = Vertical_RE_Position(i1:i2);
                    Data.RE_Position_Z = Horizontal_RE_Position(i1:i2);
                    Data.RE_Position_X = Torsion_RE_Position(i1:i2);
                else
                    Data.LE_Vel_Y = Vertical_LE_Velocity(i1:i2);
                    Data.LE_Vel_Z = Horizontal_LE_Velocity(i1:i2);
                    Data.LE_Vel_X = Torsion_LE_Velocity(i1:i2);
                end
                Data.HeadVel_X = GyroX(i1:i2);
                Data.HeadVel_Y = GyroY(i1:i2);
                Data.HeadVel_Z = GyroZ(i1:i2);
                Data.HeadVel_L = GyroLARP(i1:i2);
                Data.HeadVel_R = GyroRALP(i1:i2);
                Data.HeadAccel_X = AccelX(i1:i2);
                Data.HeadAccel_Y = AccelY(i1:i2);
                Data.HeadAccel_Z = AccelZ(i1:i2);
                Data.AllData = data1;
                Data.CSVData = ESC_CSV;
            else %LDVOG and NKI           
                Data.Trigger = Stim(i1:i2); % computer trigger
                Data.LE_Position_X = Torsion_LE_Position(i1:i2);
                Data.LE_Position_Y = Vertical_LE_Position(i1:i2);
                Data.LE_Position_Z = Horizontal_LE_Position(i1:i2);
                Data.RE_Position_X = Torsion_RE_Position(i1:i2);
                Data.RE_Position_Y = Vertical_RE_Position(i1:i2);
                Data.RE_Position_Z = Horizontal_RE_Position(i1:i2);
                Data.HeadVel_X = GyroX(i1:i2);
                Data.HeadVel_Y = GyroY(i1:i2);
                Data.HeadVel_Z = GyroZ(i1:i2);
                Data.HeadAccel_X = XAxisAccelHead(i1:i2);
                Data.HeadAccel_Y = YAxisAccelHead(i1:i2);
                Data.HeadAccel_Z = ZAxisAccelHead(i1:i2);
            end
            Data.rawfile = {info.rawfile};
            save_flag = 1;
            fname = [info.subject,'-',info.visit,'-',info.exp_date,'-',info.goggle_ver,'-',info.dataType];
            %Save but make a new ending if there are multiple segments with
            %the same information
            if exist([Seg_Path,filesep,fname,'.mat'],'file') %Already an instance of this file
                Data2 = Data; %set the segment to Data 2 to check against current file
                load([Seg_Path,filesep,fname,'.mat'],'Data')
                if ~any(ismember(Data.rawfile,Data2.rawfile)&Data.raw_start_t==Data2.raw_start_t)
                    disp([fname,' already exists in this folder and they were combined.'])
                    Data = CombineSegments(Data,Data2);
                    delete([Seg_Path,filesep,fname,'.mat'])
                else
                    disp([fname,' already exists in this folder and was ignored.'])
                    save_flag = 0;
                end
            end
            Data.info.name = [fname,'.mat'];
            if save_flag
                save([Seg_Path,filesep,fname,'.mat'],'Data')
                %Plot and save figure of the segment
                fig = plotSegment(Data);
                savefig(fig,[Seg_Path,filesep,fname,'.fig'])
                close;
            end
        end
    end
end
end