%% Segment
%Make a segmenting pipeline
%Parameters assume all LDVOG data was recorded AFTER 2016-09-26, consistent
%with recorded dates for the MVI trial
%Works for LDVOG, NKI, and GNO. Limited workability for Moog
function Segment(In_Path,Seg_Path)
%% Find the Notes File
slash = find(In_Path == filesep,1,'last');
rawfile = In_Path(slash+1:end-4);
notesfile = extractfield(dir([In_Path(1:min([strfind(In_Path,'_UpdatedTrigger')-1,length(In_Path)-4])),'*-Notes.txt']),'name');
if isempty(notesfile)
    disp(['No notes files were found for this file: ',In_Path])
    return;
elseif length(notesfile)>1
    disp(['Too many notes files were found for this file: ',In_Path])
    return;
end
%fileinfo = table2cell(readtable([In_Path(1:slash),notesfile{:}],'ReadVariableNames',false,'Format','%s %s'));
fileinfo = table2cell(readtable([In_Path(1:slash),notesfile{:}],'ReadVariableNames',false,'Delimiter',' '));
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
    if contains(info.goggle_ver,'NKI')
        VOG_time = datestr(datetime(rawfile(underscore(4)+1:underscore(7)-1),'InputFormat','hh_mm_a'),'HHMMSS');
    elseif contains(info.goggle_ver,'LDVOG')
        if ~isempty(underscore)
            rawfile = rawfile(1:underscore(1));
        end
        rawparts = split(rawfile,'-');
        VOG_time = rawparts{3};
    elseif contains(info.goggle_ver,'GNO')
        VOG_time = strrep(rawfile(underscore(5)+1:underscore(8)-1),'_','');
    else
        VOG_time = '000000';
    end
    info.exp_date = [info.exp_date,'-',VOG_time];
end
%% Load Data Types In
if contains(info.goggle_ver,'NKI')
    info.TriggerShift = 0; %measure and change
    %Suppress the warning that one of the columns is not a proper
    %column name for a MATLAB table so it got renamed
    warning('off')
    data = readtable(In_Path,'ReadVariableNames',true);
    warning('on')
    data.Properties.VariableNames{1} = 'EyeTime';
    Time_Eye = data.EyeTime;
    Time_Stim = data.EyeTime; %Update later as needed
    Stim = zeros(length(Time_Eye),1);
    Stim(data.EventCode ~= 0) = 1;
    % Load Gyroscope and accelerometer readings
    XAxisVelHead = data.GyroY - median(data.GyroY);
    YAxisVelHead  = -(data.GyroX - median(data.GyroX));
    ZAxisVelHead  = -(data.GyroZ - median(data.GyroZ));
    XAxisAccelHead  = data.AccelY - median(data.AccelY);
    YAxisAccelHead = -(data.AccelX - median(data.AccelX));
    ZAxisAccelHead = -(data.AccelZ - median(data.AccelZ));
    Fs = 1/median(abs(diff(Time_Eye)));
    % Load raw eye position data in Fick coordinates [degrees] but
    % adjust for the reverse of the X/Y/Z axes direction. This is validated by
    % experiments done on the aHIT in light on normal subjects:
    % https://docs.google.com/document/d/1euk4V0fbnbhg_paUZOaO32s_FEhStmtCbrHJV6WQLu8/edit
    Horizontal_LE_Position = -data.LeftHoriz;
    Vertical_LE_Position = -data.LeftVert;
    Torsion_LE_Position = data.LeftTorsion;
    Horizontal_RE_Position = -data.RightHoriz;
    Vertical_RE_Position = -data.RightVert;
    Torsion_RE_Position = data.RightTorsion;
    %Remove repeated values in torsion
    X1 = unique(diff(find(diff(Torsion_LE_Position)~=0)))';
    [~,max_L] = max(hist(diff(find(diff(Torsion_LE_Position)~=0)),X1).*X1);
    X2 = unique(diff(find(diff(Torsion_RE_Position)~=0)))';
    [~,max_R] = max(hist(diff(find(diff(Torsion_RE_Position)~=0)),X2).*X2);
    reps = max([X1(max_L),X2(max_R)]);
    %reps = max([median(diff(find(diff(Torsion_LE_Position)==0))),median(diff(find(diff(Torsion_RE_Position)==0)))]);
    if reps > 1
        %Loop through left array and delete duplicates
        i = 1;
        while(i<=length(Torsion_LE_Position))
            if ~isnan(Torsion_LE_Position(i))
                check_i = i+1:i+reps-1;
                check_i(check_i>length(Torsion_LE_Position)) = []; %don't check index bigger than the array
                comp = Torsion_LE_Position(i)==Torsion_LE_Position(check_i);
                if all(comp)
                    Torsion_LE_Position(check_i) = NaN;
                    i = i+reps;
                elseif any(comp)
                    Torsion_LE_Position(check_i(1:find(~comp,1,'first')-1)) = NaN;
                    i = i+find(~comp,1,'first');
                else
                    i = i+1;
                end
            else
                i = i+find(~isnan(Torsion_LE_Position(i+1:end)),1,'first');
            end
        end
        %Loop through right array and delete duplicates
        i = 1;
        while(i<=length(Torsion_RE_Position))
            if ~isnan(Torsion_RE_Position(i))
                check_i = i+1:i+reps-1;
                check_i(check_i>length(Torsion_RE_Position)) = []; %don't check index bigger than the array
                comp = Torsion_RE_Position(i)==Torsion_RE_Position(check_i);
                if all(comp)
                    Torsion_RE_Position(check_i) = NaN;
                    i = i+reps;
                elseif any(comp)
                    Torsion_RE_Position(check_i(1:find(~comp,1,'first')-1)) = NaN;
                    i = i+find(~comp,1,'first');
                else
                    i = i+1;
                end
            else
                i = i+find(~isnan(Torsion_RE_Position(i+1:end)),1,'first');
            end
        end
    end
    GyroX = reshape(XAxisVelHead,[],1);
    GyroY = reshape(YAxisVelHead,[],1);
    GyroZ = reshape(ZAxisVelHead,[],1);
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
    %Extract Eye position vectors
    HLeftIndex = 40;
    VLeftIndex = 41;
    TLeftIndex = 42;
    HRightIndex = 43;
    VRightIndex = 44;
    TRightIndex = 45;
    % Load raw eye position data in Fick coordinates [degrees]
    Horizontal_LE_Position = data(:,HLeftIndex);
    Vertical_LE_Position = data(:,VLeftIndex);
    Torsion_LE_Position = data(:,TLeftIndex);
    Horizontal_RE_Position = data(:,HRightIndex);
    Vertical_RE_Position = data(:,VRightIndex);
    Torsion_RE_Position = data(:,TRightIndex);
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
    % Index for the VOG GPIO line
    StimIndex = 35;
    Stim = data(1:length(Time_Eye),StimIndex);
    Stim = [Stim(TriggerDelay + 1:end) ; Stim(end)*ones(TriggerDelay,1)];
    %info.TriggerShift = ['UpdatedLDVOGTrigger_Shifted' num2str(TriggerDelay) 'SamplesEarlier'];
    info.TriggerShift = TriggerDelay;
    gyroscale = 1;
    accelscale = 1;
    XvelHeadIndex = 30;
    YvelHeadIndex = 29;
    ZvelHeadIndex = 28;
    XaccelHeadIndex = 27;
    YaccelHeadIndex = 26;
    ZaccelHeadIndex = 25;
    % We need to correct each gyroscope signal by subtracting the
    % correct device-specifc MPU9250 gyroscope offset. Each offset for
    % each VOG goggle ID was measured by Mehdi Rahman and posted on the
    % Google Doc located here: https://docs.google.com/a/labyrinthdevices.com/document/d/1UlZpovNkwer608aswJWdkLhF0gF-frajAdu1qgMJt9Y/edit?usp=sharing
    XvelHeadRaw = data(1:length(Time_Eye),XvelHeadIndex)*gyroscale + XvelHeadOffset;
    YvelHeadRaw = data(1:length(Time_Eye),YvelHeadIndex)*gyroscale + YvelHeadOffset;
    ZvelHeadRaw = data(1:length(Time_Eye),ZvelHeadIndex)*gyroscale + ZvelHeadOffset;
    XaccelHeadRaw = data(1:length(Time_Eye),XaccelHeadIndex)*accelscale;
    YaccelHeadRaw = data(1:length(Time_Eye),YaccelHeadIndex)*accelscale;
    ZaccelHeadRaw = data(1:length(Time_Eye),ZaccelHeadIndex)*accelscale;
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
    Rotation_Head = [
        cosd(phi) 0   sind(phi);
        0   1   0;
        -sind(phi)    0   cosd(phi)
        ];
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
    GyroLARP = data(:,3);
    GyroRALP = data(:,2);
    Fs = round(1/median(diff(Time_Eye)));    
    Horizontal_RE_Velocity = data(:,5);
    Vertical_RE_Velocity = data(:,6);  
    info.TriggerShift = 0;
    % Find the accepted traces
    warning('off')
    try
        GNO_processed = readtable(strrep(In_Path,'.txt','.csv'));
        h_ind = find(contains(GNO_processed{:,1},'Head'));
        left_imp = find(contains(GNO_processed{:,1},'Impulse')&contains(GNO_processed{:,2},'L'));
        left_ind = NaN(length(left_imp),1);
        for i = 1:length(left_imp)
            ind = find((h_ind-left_imp(i))>0,1,'first');
            left_ind(i) = h_ind(ind);
        end
        right_imp = find(contains(GNO_processed{:,1},'Impulse')&contains(GNO_processed{:,2},'R'));
        right_ind = NaN(length(right_imp),1);
        for i = 1:length(right_imp)
            ind = find((h_ind-right_imp(i))>0,1,'first');
            right_ind(i) = h_ind(ind);
        end        
        detected_left = str2double(split(strrep(GNO_processed{left_ind,1},',Head,',''),','))';
        detected_right = str2double(split(strrep(GNO_processed{right_ind,1},',Head,',''),','))';    
        if contains(In_Path,'Lateral')
            DetectedTraces_HeadVel = [detected_left,-detected_right];
        elseif contains(In_Path,'LARP')
            DetectedTraces_HeadVel = [-detected_left,detected_right];
        elseif contains(In_Path,'RALP')
            DetectedTraces_HeadVel = [-detected_left,detected_right];
        end
    catch
        DetectedTraces_HeadVel=[];
    end
    warning('on')   
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
    %Check to make sure the right canal is in the notes
    canals = {'LARP','RALP','LHRH'};
    [~,canal_i] = max(max([GyroLARP,GyroRALP,GyroZ]));
    notes_canal = find([any(contains(stim_info,'LARP')),any(contains(stim_info,'RALP')),any(contains(stim_info,'LHRH'))]);
    if canal_i~=notes_canal %Mismatch
        plot(GyroLARP,'Color','g')
        hold on
        plot(GyroRALP,'Color','b')
        plot(GyroZ,'Color','r')
        hold off
        prompt = [In_Path,newline,'Notes file says ',canals{notes_canal},' but the motion traces appear to be in the ',canals{canal_i},' canal.',newline,'Which do you want to use for segmenting?'];
        canal_choice = questdlg(prompt,'',canals{notes_canal},canals{canal_i},'Stop process','Stop process');
        if strcmp(canal_choice,'Stop process')
            error(['Mismatch in canal for file: ',In_Path])
        else
            [~,canal] = ismember(canal_choice,canals);
        end
    else
        canal = notes_canal;
    end
    switch canal
        case 1
            GyroAll = GyroLARP;
        case 2
            GyroAll = GyroRALP;
        case 3
            GyroAll = GyroZ;
    end
    if all(contains(stim_info,'VelStep')) % Velstep
        thresh = 50; %Adjust as needed
        thresh2 = 1; %Counts as 0.
        approx0 = find(abs(GyroAll)< thresh2);
        start1 = find(abs(GyroAll) > thresh);
        start2 = start1([true;diff(start1)>Fs*45]); %make sure it's at least 45s long
        [~,inds] = min(abs(repmat(approx0',length(start2),1) - start2),[],2);
        start = approx0(inds);
        start = start([true;diff(start)>Fs*45]); %make sure it's at least 45s long
        ends = [start(2:end)-1;length(Time_Eye)];
        plot(Time_Stim,GyroAll,'k',Time_Stim(start),GyroAll(start),'r*',Time_Stim(ends),GyroAll(ends),'b*')
        xlabel('Time (s)')
        ylabel('Head Velocity')
        legend('Stimulus','Start','End')
        pause(1)
    else %Sine or Impulses
        if length(stim_info) > 1 %multiple sine stimuli
            %Segment and allow user to select which segments are "real"
            thresh = 20; %change as needed but works for 0.01Hz - 2Hz
            abov_thresh = find(abs(GyroAll) > thresh);
            starts1 = abov_thresh([true;diff(abov_thresh)>1]);
            ends1 = abov_thresh([diff(abov_thresh)>1;true]);
            pos = find(GyroAll >= 0);
            neg = find(GyroAll <= 0);
            for i = 1:length(starts1)
                if GyroAll(starts1(i)) > 0
                    ind1 = neg(find(neg<starts1(i),1,'last'));
                else
                    ind1 = pos(find(pos<starts1(i),1,'last'));
                end
                if GyroAll(ends1(i)) > 0
                    ind2 = neg(find(neg>ends1(i),1,'first'));
                else
                    ind2 = pos(find(pos>ends1(i),1,'first'));
                end
                if isempty(ind1)
                    ind1 = 1;
                end
                if isempty(ind2)
                    ind2 = length(GyroAll);
                end
                if ismember(ind1,starts1)
                    starts1(i) = NaN;
                    ends1(i) = NaN;
                else
                    starts1(i) = ind1;
                    ends1(i) = ind2;
                end
            end
            starts1(isnan(starts1)) = [];
            ends1(isnan(ends1)) = [];
            % Combine segments with only a few samples between them
            for i = 2:length(starts1)
                if abs(ends1(i-1)-starts1(i)) < 2*Fs %number of samples between segments for them to still be counted together. change as needed.
                    starts1(i) = starts1(i-1);
                    starts1(i-1) = NaN;
                    ends1(i-1) = NaN;
                end
            end
            starts1(isnan(starts1)) = [];
            ends1(isnan(ends1)) = [];
            seg_start = max(round(starts1-20,0),1);
            seg_end = min(round(ends1,0)+20,length(GyroAll));
            if(length(seg_start) ~= length(seg_end))
                error('Length of start and stop are unequal. Please manually segment')
            elseif length(seg_start) < length(stim_info)
                error('Less segments found than in the notes file. Please manually segment')
            end
            if length(seg_start) ~= length(stim_info)
                plot(NaN,NaN)
                hold on
                %Now plot all fills
                for j = 1:length(seg_start)
                    fill([Time_Eye(seg_start(j)),Time_Eye(seg_end(j)),Time_Eye(seg_end(j)),Time_Eye(seg_start(j))]',[500,500,-500,-500]',[.85, .85, .85]);
                end
                plot(Time_Eye,GyroLARP,'k:',Time_Eye,GyroRALP,'k--',Time_Eye,GyroZ,'k-')
                hold off
                uiwait(msgbox('Select all valid segments.'))
                keep = false(1,length(seg_start));
                [x,~] = ginput(length(stim_info));
                for i = 1:length(x)
                    t1 = Time_Eye(seg_start) - x(i);
                    i1 = find(t1 < 0);
                    keep(i1(end)) = true;
                end
                start = seg_start(keep);
                ends = seg_end(keep);
            else
                start = seg_start;
                ends = seg_end;
            end
            plot(NaN,NaN)
            hold on
            %Now plot all fills
            for j = 1:length(start)
                fill([Time_Eye(start(j)),Time_Eye(ends(j)),Time_Eye(ends(j)),Time_Eye(start(j))]',[500,500,-500,-500]','g');
            end
            plot(Time_Eye,GyroLARP,'k:',Time_Eye,GyroRALP,'k--',Time_Eye,GyroZ,'k-')
            hold off
            pause(1)
        else % low freq sine
            %Just take the whole thing as a segment
            start = 1;
            ends = length(Time_Eye);
            plot(NaN,NaN)
            hold on
            fill([Time_Eye(start),Time_Eye(ends),Time_Eye(ends),Time_Eye(start)]',[500,500,-500,-500]','g');
            plot(Time_Eye,GyroLARP,'k:',Time_Eye,GyroRALP,'k--',Time_Eye,GyroZ,'k-')
            hold off
            pause(1)
        end
        if any(contains(stim_info,'Impulse'))
            %Add a stim_info entry for each canal (LHRH -> LH and RH)
            rep_ind = sort([(1:length(start))';find(contains(stim_info,'Impulse'))]);
            start = start(rep_ind);
            ends = ends(rep_ind); 
            stim_info2 = repmat(stim_info,1,2);
            stim_info2(~contains(stim_info2(:,2),'Impulse'),2) = {''};
            %Make the canals one-sided here
            stim_info2(contains(stim_info2(:,1),'Impulse')&contains(stim_info2(:,1),'LHRH'),1) = strrep(stim_info2(contains(stim_info2(:,1),'Impulse')&contains(stim_info2(:,1),'LHRH'),1),'LHRH','LH');
            stim_info2(contains(stim_info2(:,1),'Impulse')&contains(stim_info2(:,2),'LHRH'),2) = strrep(stim_info2(contains(stim_info2(:,1),'Impulse')&contains(stim_info2(:,2),'LHRH'),2),'LHRH','RH');
            stim_info2(contains(stim_info2(:,1),'Impulse')&contains(stim_info2(:,1),'LARP'),1) = strrep(stim_info2(contains(stim_info2(:,1),'Impulse')&contains(stim_info2(:,1),'LARP'),1),'LARP','RP');
            stim_info2(contains(stim_info2(:,1),'Impulse')&contains(stim_info2(:,2),'LARP'),2) = strrep(stim_info2(contains(stim_info2(:,1),'Impulse')&contains(stim_info2(:,2),'LARP'),2),'LARP','LA');
            stim_info2(contains(stim_info2(:,1),'Impulse')&contains(stim_info2(:,1),'RALP'),1) = strrep(stim_info2(contains(stim_info2(:,1),'Impulse')&contains(stim_info2(:,1),'RALP'),1),'RALP','RA');
            stim_info2(contains(stim_info2(:,1),'Impulse')&contains(stim_info2(:,2),'RALP'),2) = strrep(stim_info2(contains(stim_info2(:,1),'Impulse')&contains(stim_info2(:,2),'RALP'),2),'RALP','LP');
            stim_info = reshape(stim_info2',[],1);
            stim_info(cellfun(@isempty,stim_info)) = [];                        
        end 
    end
elseif all(contains(stim_info,{'eeVOR','trash'})) %Good for all externally triggered stimuli for now
    if all(contains(stim_info,'Activation'))
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
    else %Should work for Sine, Pulse Trains and Multi Vector
        %Every trigger toggle is a cycle
        temp = find(abs(diff(Stim))==1);
        all_starts = temp(1:end-1);
        all_ends = temp(2:end);
        temp2 = diff(temp);
        start = all_starts;
        ends = all_ends;
        thresh_tol = Fs*0.4; %tolerance for differences in cyc length (always at least 400ms of break)
        %thresh_tol = 50; %old
        for i = 2:length(all_starts)
            if abs(temp2(i)-temp2(i-1))<thresh_tol
                start(i) = start(i-1);
                start(i-1) = NaN;
                ends(i-1) = NaN;
            elseif start(i-1)==all_starts(i-1)
                start(i-1) = NaN;
                ends(i-1) = NaN;
            end
        end
        temp2(isnan(start)) = [];
        start(isnan(start)) = [];
        ends(isnan(ends)) = [];
        start(~ismember(start,ends)) = start(~ismember(start,ends))-temp2(~ismember(start,ends));
        ends = ends+temp2;
        start(start<0) = 1;
        ends(ends>length(Time_Stim)) = length(Time_Stim);
    end
    stim = Stim;
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
            elseif contains(stim_info{i},'Impulse')
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
            else %LDVOG and NKI           
                Data.Trigger = Stim(i1:i2); % computer trigger
                Data.LE_Position_X = Torsion_LE_Position(i1:i2);
                Data.LE_Position_Y = Vertical_LE_Position(i1:i2);
                Data.LE_Position_Z = Horizontal_LE_Position(i1:i2);
                Data.RE_Position_X = Torsion_RE_Position(i1:i2);
                Data.RE_Position_Y = Vertical_RE_Position(i1:i2);
                Data.RE_Position_Z = Horizontal_RE_Position(i1:i2);
                Data.HeadVel_X = XAxisVelHead(i1:i2);
                Data.HeadVel_Y = YAxisVelHead(i1:i2);
                Data.HeadVel_Z = ZAxisVelHead(i1:i2);
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