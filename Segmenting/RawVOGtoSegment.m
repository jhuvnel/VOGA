%% Raw VOG to Segment
%
% Works for LDVOG, NKI, GNO, and ESC. Limited workability for Moog. This is
% function that takes the output of the VOG recording system and turns it
% into a Data struct. This has not yet partitioned it by experiment.
%
function [Data,stim_info] = RawVOGtoSegment(In_Path)
%% Load and Parse Notes File
% Find and load in the the Notes File 
% Sometimes the Notes file name won't exactly match the VOG file name
% because of extensions added to the VOG file to signify updating the
% trigger or something similar. All notes are after an _Updated in the file
% name so ignore text after that if necessary.
slash = find(In_Path == filesep,1,'last');
dot = find(In_Path == '.',1,'last');
rawfile = In_Path(slash+1:dot-1);
notesfile = extractfield(dir([In_Path(1:min([strfind(In_Path,'_Updated')-1,dot-1])),'*-Notes.txt']),'name');
if isempty(notesfile)
    error(['No notes files were found for this file: ',In_Path])
elseif length(notesfile)>1
    error(['Too many notes files were found for this file: ',In_Path])
end
fileinfo = table2cell(readtable([In_Path(1:slash),notesfile{:}],'ReadVariableNames',false,'Format','%s %s'));
if isempty(fileinfo{7,2})
    fileinfo = table2cell(readtable([In_Path(1:slash),notesfile{:}],'ReadVariableNames',false,'Delimiter',' '));
end
% Figure out how many experiments there are and save them into a stim_info
% variable.
if size(fileinfo,2)==2 %New way w/2 columns
    stim_info = fileinfo(7:end,2);
elseif length(strsplit(fileinfo{7},' ')) > 1
    stim_info = strrep(fileinfo(8:end),'"','');
else %Only one experiment type
    stim_info = strrep(strcat(fileinfo{7},'-',fileinfo(8:end)),'"','');
end
stim_info = strrep(stim_info,'manual','Manual'); %in case I forgot to capitalize
%Put all of the relevant info into the info struct
info.rawfile = In_Path;
info.rawnotes = [In_Path(1:slash),notesfile{:}];
info.subject = fileinfo{1,end};
info.ear = fileinfo{2,end};
info.visit = strrep(fileinfo{3,end},' ','');
info.exp_date = fileinfo{4,end};
info.goggle_ver = fileinfo{5,end}; %This should say NKI or not
info.goggle_reorient_ang = str2double(fileinfo{6,end});
info.TriggerShift = 0; %this is the sample # offset between goggle and MVI fitting software trigger, should be measured for each goggle set that does eeVOR
% Fix the date and time in YYYYMMDD-hhmmss format if it's not already there
if length(info.exp_date)==14 %Just missing '-' in between date and time
    info.exp_date = [info.exp_date(1:8),'-',info.exp_date(9:end)];
elseif length(info.exp_date)==8 %Just date, add time
    VOG_time = '000000'; %Default
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
    end
    info.exp_date = [info.exp_date,'-',VOG_time];
end
%% Load Data In for Each Type of Goggle
% For each goggle, create these variables in a struct named Data
% if available: Fs, Time_Eye, Time_Stim, Trigger, 
% LE_Position_X, LE_Position_Y, LE_Position_Z,
% RE_Position_X, RE_Position_Y, RE_Position_Z, 
% HeadVel_X, HeadVel_Y, HeadVel_Z, 
% HeadAccel_X, HeadAccel_Y, and HeadAccel_Z
Data.info = info;
if contains(info.goggle_ver,{'NKI','NL'})
    %Load in both the .dat and .csv files for a Neurolign file. The .csv
    %has the properly calibrated eye data but the .dat file has more
    %reliable gyro and trigger line data
    warning('off') %Suppress the warning that one of the columns is not a proper column name for a MATLAB table so it got renamed
    data = readtable(In_Path,'ReadVariableNames',true); %.dat file with uncalibrated eye data
    data2 = table2array(readtable(strrep(info.rawnotes,'-Notes.txt','.csv'))); %.csv with calibrated eye data
    warning('on')
    % Truncate the longer file if they aren't the same length
    len = min([size(data,1),size(data2,1)]);
    data = data(1:len,:); data2 = data2(1:len,:);
    % Set the values of Data
    Data.Fs = 1/median(abs(diff(data2(:,1))));
    Data.Time_Eye = data2(:,1); Data.Time_Stim = data2(:,1);
    Data.Trigger = zeros(length(Data.Time_Eye),1);
    Data.Trigger(data.EventCode ~= 0) = 1; %Any time the EventCode field isnt' 0, it's high
    % Load gyroscope, acceleration, and raw eye position data in Fick
    % coordinates [degrees].
    % Adjust for the camera values being reversed. This is validated by
    % Experiments done on the aHIT in light on normal subjects can be found here:
    % https://docs.google.com/document/d/1euk4V0fbnbhg_paUZOaO32s_FEhStmtCbrHJV6WQLu8/edit
    % Shows that camera Y and Z values are both mirrored.
    % Remove repeated values in torsion by detecting how many samples in a
    % row are part of the first order hold and doing a moving average
    % instead.
    Torsion_LE_Position = data2(:,10); Torsion_RE_Position = data2(:,13);
    [gr,gc] = groupcounts(diff(find(diff([Torsion_LE_Position;Torsion_RE_Position])~=0&...
        ~isnan([Torsion_LE_Position(2:end);Torsion_RE_Position]))));
    reps = gc(gr==max(gr));
    if reps>1
        M = [movmean(Torsion_LE_Position,reps);zeros(floor(reps/2)-1,1)];
        Torsion_LE_Position = [M(round(reps/2):end);zeros(round(reps/2)-1,1)];
        M = [movmean(Torsion_RE_Position,reps);zeros(floor(reps/2)-1,1)];
        Torsion_RE_Position = [M(round(reps/2):end);zeros(round(reps/2)-1,1)];        
    end 
    % The gryo/accelerometers are oriented with the left hand rule
    % and therefore a passive rotation cannot not be used to compensate
    % unless we already switch the gyros/accelerometers labeled "X" and "Y"
    psi = 180; %Based on +y pointing out of the right ear and +z pointing down
    Rotation_Head = [1,0,0;0,cosd(psi),-sind(psi);0,sind(psi),cosd(psi)];
    % NOTE: We are transposing the rotation matrix in order to apply a
    % PASSIVE (i.e., a coordinate system) transformation
    A = Rotation_Head'*[(data.GyroY-median(data.GyroY))';(data.GyroX-median(data.GyroX))';(data.GyroZ-median(data.GyroZ))'];
    B = Rotation_Head'*[(data.AccelY-median(data.AccelY))';(data.AccelX-median(data.AccelX))';(data.AccelZ-median(data.AccelZ))'];
    % We want to rotate our data from an [X,Y,Z] coordinate system,
    % into a [LARP,RALP,LHRH] coordinate system, where Z = LHRH. To
    % accomplish this, we will perform a PASSIVE (i.e., 'alias' or
    % 'coordinate system') rotation of -45 degrees. In order to realize
    % this rotation, I will generate a -45deg rotation matrix and
    % RIGHT multiple our data by the TRANSPOSE of this rotation
    % matrix. Note that rotation matrices are orthonormal, and
    % their inverses are equivalent to their transpose. -PJB
    C = (rotZ3deg(-45)'*A); 
    Data.LE_Position_X = Torsion_LE_Position;
    Data.LE_Position_Y = -data2(:,4);
    Data.LE_Position_Z = -data2(:,2);
    Data.RE_Position_X = Torsion_RE_Position;
    Data.RE_Position_Y = -data2(:,12);
    Data.RE_Position_Z = -data2(:,3);   
    Data.HeadVel_X = medfilt1(reshape(A(1,:),[],1),3); %remove outliers/droupouts
    Data.HeadVel_Y = medfilt1(reshape(A(2,:),[],1),3);
    Data.HeadVel_Z = medfilt1(reshape(A(3,:),[],1),3);
    Data.HeadVel_L = medfilt1(reshape(C(1,:),[],1),3);
    Data.HeadVel_R = medfilt1(reshape(C(2,:),[],1),3);
    Data.HeadAccel_X = medfilt1(reshape(B(1,:),[],1),3); 
    Data.HeadAccel_Y = medfilt1(reshape(B(2,:),[],1),3); 
    Data.HeadAccel_Z = medfilt1(reshape(B(3,:),[],1),3); 
elseif contains(info.goggle_ver,'LDVOG')
    % Load data, remove timestamps if they are being read as doubles
    data1 = readtable(In_Path);
    data = data1{:,contains(varfun(@class,data1,'OutputFormat','cell'),'double')};
    % Around 2018-04, PJB noticed that the LD VOG system appeared to record the
    % PCU trigger LATE relative to the collected eye movement data. PJB, MR,
    % and NV performed some experiments outlined here:
    % https://docs.google.com/document/d/16EppbHOlsSabjR61xbDi798p7ldrqqTV6tOXjH1KZbI/edit#heading=h.xh2azhw75tvn
    % These experiments resulted in the following trigger recording delay in
    % samples (where positive values indicate the trigger was recorded AFTER
    % the eye movement).
    % We also need to correct each gyroscope signal by subtracting the
    % correct device-specifc MPU9250 gyroscope offset. Each offset for
    % each VOG goggle ID was measured by Mehdi Rahman and posted on the
    % Google Doc located here: https://docs.google.com/a/labyrinthdevices.com/document/d/1UlZpovNkwer608aswJWdkLhF0gF-frajAdu1qgMJt9Y/edit?usp=sharing
    TriggerDelay = 0; XvelHeadOffset = 0; YvelHeadOffset = 0; ZvelHeadOffset = 0; % Defaults
    switch info.goggle_ver
        case {'1','LDVOG1'}
            TriggerDelay = 3; XvelHeadOffset = 2.7084; YvelHeadOffset = -0.5595; ZvelHeadOffset = 0.7228;
        case {'2','LDVOG2'}
            TriggerDelay = 5; XvelHeadOffset = 2.3185; YvelHeadOffset = -1.5181; ZvelHeadOffset = -1.0424;
        case {'3','LDVOG3'}
            TriggerDelay = 2; XvelHeadOffset = 1.9796; YvelHeadOffset = 0.1524; ZvelHeadOffset = -0.5258;           
    end
    Head_Sensor_Latency = 0.047; % From Mehdi Rahman bench tests, the data acquisition of the MPU9250 leads the LD VOG Goggles by 47ms
    % Set Data values
    Time = data(:,2);
    Data.Fs = 1/median(abs(diff(Time)));
    % The time vector recorded by the VOG goggles resets after it reaches a value of 128.
    Data.Time_Eye = (0:length(Time)-1)*median(diff(Time));
    Data.Time_Stim = Data.Time_Eye - Head_Sensor_Latency;
    Stim = data(1:length(Time),35);
    Data.Trigger = [Stim(TriggerDelay + 1:end) ; Stim(end)*ones(TriggerDelay,1)];
    Data.info.TriggerShift = TriggerDelay;
    % Load raw eye position data in Fick coordinates [degrees]
    Data.LE_Position_X = data(:,42); Data.LE_Position_Y = data(:,41); Data.LE_Position_Z = data(:,40);
    Data.RE_Position_X = data(:,45); Data.RE_Position_Y = data(:,44); Data.RE_Position_Z = data(:,43);   
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
    A = Rotation_Head'*[(data(1:length(Time),30)+XvelHeadOffset)';(data(1:length(Time),29)+YvelHeadOffset)';(data(1:length(Time),28)+ZvelHeadOffset)'];
    B = Rotation_Head'*[data(1:length(Time),27)';data(1:length(Time),26)';data(1:length(Time),25)'];
    % We want to rotate our data from an [X,Y,Z] coordinate system,
    % into a [LARP,RALP,LHRH] coordinate system, where Z = LHRH. To
    % accomplish this, we will perform a PASSIVE (i.e., 'alias' or
    % 'coordinate system') rotation of -45 degrees. In order to realize
    % this rotation, I will generate a -45deg rotation matrix and
    % RIGHT multiple our data by the TRANSPOSE of this rotation
    % matrix. Note that rotation matrices are orthonormal, and
    % their inverses are equivalent to their transpose. -PJB
    C = (rotZ3deg(-45)'*A); 
    Data.HeadVel_X = reshape(A(1,:),[],1); Data.HeadVel_Y = reshape(A(2,:),[],1); Data.HeadVel_Z = reshape(A(3,:),[],1);
    Data.HeadVel_L = reshape(C(1,:),[],1); Data.HeadVel_R = reshape(C(2,:),[],1);  
    Data.HeadAccel_X = reshape(B(1,:),[],1); Data.HeadAccel_Y = reshape(B(2,:),[],1); Data.HeadAccel_Z = reshape(B(3,:),[],1); 
elseif contains(info.goggle_ver,'GNO') % GNO File
    %Load the data
    data = table2array(readtable(In_Path));
    Data.Fs = round(1/median(diff(data(:,1))))*10e6;
    Data.Time_Eye = (data(:,1) - data(1,1))/10e6; Data.Time_Stim = Data.Time_Eye;
    % The only eye traces available from GNO
    Data.RE_Vel_Y = data(:,6); Data.RE_Vel_Z = data(:,5);
    %Inverted LARP gyro to adjust for +L canal, wasn't negative here until 01/10/22
    Data.HeadVel_Z = data(:,4); Data.HeadVel_L = -data(:,3); Data.HeadVel_R = data(:,2); 
    % Add the information from the other files to segment in case you need
    % to reference it later
    Data.DetectedTraces_HeadVel = []; Data.CSVData = []; Data.XMLData = []; %Defaults
    if isfile(strrep(In_Path,'.txt','.csv'))  % Find the accepted traces
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
            if isscalar(head_ind)
                imp_dat(i,1) = extract(GNO_CSV{imp_ind(i)},"L"|"R");
                vec = str2double(split(strrep(GNO_CSV(imp_ind(i)-1+head_ind),'Head',''),','));
                imp_dat{i,2} = vec(~isnan(vec));
            end
        end
        Data.CSVData = GNO_CSV;
        Data.DetectedTraces_HeadVel = [horzcat(imp_dat{contains(imp_dat(:,1),'L'),2}),-horzcat(imp_dat{contains(imp_dat(:,1),'R'),2})];
        if ~contains(In_Path,'Lateral')
            Data.DetectedTraces_HeadVel = -Data.DetectedTraces_HeadVel; %Flip them to be in the correct +LRZ orientation 
        end 
    end
    if isfile(strrep(In_Path,'.txt','.xml')) %Put xml data in the segment too
        Data.XMLData = cellstr(readlines(strrep(In_Path,'.txt','.xml')));
    end
elseif contains(info.goggle_ver,{'ESC1','ESC2'}) % ESC File with .mat
    data1 = load(In_Path);
    %Get values to calibrate the gyro
    phi = -20; motion_cal = eye(3); %default
    if isfield(data1,'content_RAW')
        data = data1.content_RAW.Data;        
        if isfield(data1.content_CAL,'R')
            motion_cal = ([cosd(phi),0,sind(phi);0,1,0;-sind(phi),0,cosd(phi)]'*data1.content_CAL.R);
        end
        labs = cellstr(data1.content_RAW.DataNames);
    elseif isfield(data1,'Data')
        data = data1.Data;
        labs = cellstr(data1.DataNames);  
    end 
    EyeTimeIndex = find(contains(labs,'SystemTime'));
    if ~isempty(EyeTimeIndex) % Most cases
        old_Time_Eye = data(:,EyeTimeIndex) - data(1,EyeTimeIndex);
        keep_eye = [true;diff(old_Time_Eye)>0];
        Data.Fs = 1/median(diff(old_Time_Eye(keep_eye))); 
        Data.Time_Eye = 0:(1/Data.Fs):old_Time_Eye(end);
        Data.Time_Stim = Data.Time_Eye;
        warning('off')
        Data.LE_Position_X = spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'EyeVelX')),Data.Time_Eye); 
        Data.LE_Position_Y = spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'EyeVelY')),Data.Time_Eye);
        Data.LE_Position_Z = spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'EyeVelZ')),Data.Time_Eye);
        A = motion_cal*[spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'InertialVelX')),Data.Time_Eye);... 
            spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'InertialVelY')),Data.Time_Eye);... 
            spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'InertialVelZ')),Data.Time_Eye)];         
        B = motion_cal*[spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'InertialAccelX')),Data.Time_Eye);... 
            spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'InertialAccelY')),Data.Time_Eye);... 
            spline(old_Time_Eye(keep_eye),data(keep_eye,contains(labs,'InertialAccelZ')),Data.Time_Eye)]; 
        % We want to rotate our data from an [X,Y,Z] coordinate system,
        % into a [LARP,RALP,LHRH] coordinate system, where Z = LHRH. To
        % accomplish this, we will perform a PASSIVE (i.e., 'alias' or
        % 'coordinate system') rotation of -45 degrees. In order to realize
        % this rotation, I will generate a -45deg rotation matrix and
        % RIGHT multiple our data by the TRANSPOSE of this rotation
        % matrix. Note that rotation matrices are orthonormal, and
        % their inverses are equivalent to their transpose. -PJB
        C = (rotZ3deg(-45)'*A); 
        Data.HeadVel_X = reshape(A(1,:),[],1); Data.HeadVel_Y = reshape(A(2,:),[],1); Data.HeadVel_Z = reshape(A(3,:),[],1);
        Data.HeadVel_L = reshape(C(1,:),[],1); Data.HeadVel_R = reshape(C(2,:),[],1);
        Data.HeadAccel_X = reshape(B(1,:),[],1); Data.HeadAccel_Y = reshape(B(2,:),[],1); Data.HeadAccel_Z = reshape(B(3,:),[],1); 
    elseif any(contains(labs,'RealTime'))
        %Old version missing much data
        len = sum(~isnan(data(:,strcmp(labs,'Time'))));
        Data.Fs = 10^3; %Time vec is in msec
        Data.Time_Eye = (0:(len-1))*10-3; %now in sec
        Data.Time_Stim = Data.Time_Eye;
        if any(contains(reshape(fileinfo,[],1),'LHRH'))
            Data.LE_Vel_Z = data(~isnan(data(:,strcmp(labs,'Time'))),contains(labs,'EyeVel'));
            Data.HeadVel_Z = data(~isnan(data(:,strcmp(labs,'Time'))),contains(labs,'HeadVel'));
        elseif any(contains(reshape(fileinfo,[],1),'LARP'))
            Data.LE_Vel_Y = data(~isnan(data(:,strcmp(labs,'Time'))),contains(labs,'EyeVel'));
            Data.HeadVel_L = data(~isnan(data(:,strcmp(labs,'Time'))),contains(labs,'HeadVel'));
        elseif any(contains(reshape(fileinfo,[],1),'RALP'))
            Data.LE_Vel_Y = data(~isnan(data(:,strcmp(labs,'Time'))),contains(labs,'EyeVel'));
            Data.HeadVel_R = data(~isnan(data(:,strcmp(labs,'Time'))),contains(labs,'HeadVel'));
        end       
    else
        disp(['Not segmented: ',In_Path])
        return;
    end
    Data.AllData = data1;
    Data.CSVData = [];
    if isfile(strrep(strrep(In_Path,'_export',''),'.mat','.xls'))
        fname = strrep(strrep(In_Path,'_export',''),'.mat','.xls');
        movefile(fname,strrep(fname,'.xls','.csv'))
        ESC_CSV = cellstr(readlines(strrep(fname,'.xls','.csv')));
        ESC_CSV(cellfun(@isempty,ESC_CSV)) = [];
        Data.CSVData = split(ESC_CSV,',');
        movefile(strrep(fname,'.xls','.csv'),fname)
    end    
elseif contains(info.goggle_ver,'ESC3') % ESC File with .csv data
    data = readtable(In_Path);
    Data.Fs = 1/median(diff(data.LeftTime)); 
    Data.Time_Eye = data.LeftTime - data.LeftTime(1);
    Data.Time_Stim = Data.Time_Eye; %Update later as needed
    % Load Gyroscope and accelerometer readings
    % Load raw eye position data in Fick coordinates [degrees] with no torsion.
    % Gyros appear to be properly aligned with +XYZ. This is shown here:
    % https://docs.google.com/document/d/1euk4V0fbnbhg_paUZOaO32s_FEhStmtCbrHJV6WQLu8/edit
    % Set torsion to 0 so it can do 3D vel calculations
    Data.RE_Position_X = 0*data.LeftPupilCol; Data.RE_Position_Y = data.LeftPupilRow; Data.RE_Position_Z = data.LeftPupilCol;
    Data.HeadVel_X = data.HeadInertialVelX; Data.HeadVel_Y = data.HeadInertialVelY; Data.HeadVel_Z = data.HeadInertialVelZ;
    % We want to rotate our data from an [X,Y,Z] coordinate system,
    % into a [LARP,RALP,LHRH] coordinate system, where Z = LHRH. To
    % accomplish this, we will perform a PASSIVE (i.e., 'alias' or
    % 'coordinate system') rotation of -45 degrees. In order to realize
    % this rotation, I will generate a -45deg rotation matrix and
    % RIGHT multiple our data by the TRANSPOSE of this rotation
    % matrix. Note that rotation matrices are orthonormal, and
    % their inverses are equivalent to their transpose. -PJB
    % % % FIX ME To get LARP and RALP gyro
    %C = (rotZ3deg(-45)'*[data.HeadInertialVelX;data.HeadInertialVelY;data.HeadInertialVelZ]')'; 
    %Data.HeadVel_L = reshape(C(1,:),[],1); Data.HeadVel_R = reshape(C(2,:),[],1);
    Data.HeadAccel_X = data.HeadInertialAccelX; Data.HeadAccel_Y = data.HeadInertialAccelY; Data.HeadAccel_Z = data.HeadInertialAccelZ;
    Data.AllData = data;
    ESC_CSV = cellstr(readlines(strrep(In_Path,'ImuData','Numerical')));
    ESC_CSV(cellfun(@isempty,ESC_CSV)) = [];
    max_comma = max(count(ESC_CSV,','));
    comma_count = count(ESC_CSV,',');
    while any(comma_count~=max_comma)
        ESC_CSV(comma_count~=max_comma) = strcat(ESC_CSV(comma_count~=max_comma),',');        
        comma_count = count(ESC_CSV,',');
    end
    Data.CSVData = split(ESC_CSV,',');
end
end