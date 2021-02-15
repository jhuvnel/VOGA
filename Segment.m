%% Segment
%Make a segmenting pipeline
%Parameters assume all LDVOG data was recorded AFTER 2016-09-26, consistent
%with recorded dates for the MVI trial
function Segment(In_Path,Seg_Path)
    %Make sure there is a real path for the input/output files
    if isnumeric(In_Path) || isempty(In_Path)
        error('Input path to Segment.m not valid.') 
    elseif isnumeric(Seg_Path) || isempty(Seg_Path)
        error('Input path to Segment.m not valid.') 
    end    
    %See if notes file already exists and create one if it doesn't
    try
        fileinfo = importdata([In_Path(1:end-4),'-Notes.txt']);   
    catch
        slash = find(In_Path == filesep,1,'last');
        MakeNotes(In_Path(1:slash-1),{In_Path(slash+1:end)})
        fileinfo = importdata([In_Path(1:end-4),'-Notes.txt']);   
    end
    %Change this to reflect new notes file types when they are created 
    info.rawfile = In_Path;
    info.rawnotes = [In_Path(1:end-4),'-Notes.txt'];
    info.subject = fileinfo{1}(2:end-1);
    info.ear = fileinfo{2}(2:end-1);
    info.visit = strrep(fileinfo{3}(2:end-1),' ','');
    info.exp_date = fileinfo{4}(2:end-1);
    info.goggle_ver = fileinfo{5}(2:end-1); %This should say NKI or not
    info.goggle_reorient_ang = str2double(fileinfo{6}(2:end-1));
    if contains(info.goggle_ver,'NKI')
        %Suppress the warning that one of the columns is not a proper
        %column name for a MATLAB table so it got renamed
        warning('off')
        data = readtable(In_Path,'ReadVariableNames',true);
        warning('on')
        data.Properties.VariableNames{1} = 'EyeTime';
        Time_Eye = data.EyeTime;  
        Time_Stim = data.EyeTime; %Update later as needed
        StimAll = zeros(length(Time_Eye),1);
        StimAll(data.EventCode ~= 0) = 1;
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
        StimAll = data(1:length(Time_Eye),StimIndex);
        StimAll = [StimAll(TriggerDelay + 1:end) ; StimAll(end)*ones(TriggerDelay,1)];
        info.TriggerShift = ['UpdatedLDVOGTrigger_Shifted' num2str(TriggerDelay) 'SamplesEarlier']; 
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
    end    
    %% Figure out how many experiments there are
    exps = strsplit(fileinfo{7}(2:end-1),' ');
    if length(exps) > 1
        plot(XAxisVelHead,'k:')
        hold on
        plot(YAxisVelHead,'k--')
        plot(ZAxisVelHead,'k-')
        plot(100*StimAll,'b')
        hold off
        uiwait(msgbox('Select the boundaries between experiment types.'))
        [x,~] = ginput(length(exps)-1);
        is = [1;round(x,0)];
        ie = [round(x,0);length(Time_Eye)];
        stim_info_all = strrep(fileinfo(8:end),'"','');
        stim_info_ind = ones(length(stim_info_all),1);
        k = 1;
        for j = 1:length(stim_info_ind) %Helps in the case of repeated experiment types
            if ~contains(stim_info_all(j),exps{k})
                k = k+1;
            end
            stim_info_ind(j) = k;
        end
    else %Only one experiment type
        stim_info_all = strrep(strcat(fileinfo{7}(2:end-1),'-',fileinfo(8:end)),'"','');
        stim_info_ind = ones(1,length(stim_info_all));
        is = 1;
        ie = length(Time_Eye);
    end
    %% Segment
    for q = 1:length(exps)
        exp = exps(q);
        stim_info = stim_info_all(stim_info_ind==q);
        GyroX = zeros(length(Time_Eye),1);
        GyroY = zeros(length(Time_Eye),1);
        GyroZ = zeros(length(Time_Eye),1);
        Stim = zeros(length(Time_Eye),1);
        GyroX(is(q):ie(q)) = reshape(XAxisVelHead(is(q):ie(q)),[],1);
        GyroY(is(q):ie(q)) = reshape(YAxisVelHead(is(q):ie(q)),[],1);
        GyroZ(is(q):ie(q)) = reshape(ZAxisVelHead(is(q):ie(q)),[],1);
        Stim(is(q):ie(q)) = StimAll(is(q):ie(q));       
        if contains(exp,{'RotaryChair','aHIT'})           
            %Check to make sure the right canal is in the notes
            canals = {'LARP','RALP','LHRH'};
            GyroLARP = (GyroX - GyroY)/sqrt(2);
            GyroRALP = (GyroX + GyroY)/sqrt(2);            
            [~,canal_i] = max(max([GyroLARP,GyroRALP,GyroZ]));
            notes_canal = find([any(contains(stim_info,'LARP')),any(contains(stim_info,'RALP')),any(contains(stim_info,'LHRH'))]);
            if canal_i~=notes_canal %Mismatch
                plot(GyroLARP,'k:')
                hold on
                plot(GyroRALP,'k--')
                plot(GyroZ,'k-')
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
            if contains(exp,'Sine')           
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
                        if abs(ends1(i-1)-starts1(i)) < Fs %number of samples between segments for them to still be counted together. change as needed.
                            starts1(i) = starts1(i-1);
                            starts1(i-1) = NaN;
                            ends1(i-1) = NaN;
                        end
                    end
                    starts1(isnan(starts1)) = [];
                    ends1(isnan(ends1)) = [];
                    seg_start = max(round(starts1-20,0),1);
                    seg_end = min(round(ends1,0)+20,length(GyroAll));
                    plot(Time_Eye,GyroX,'k:')
                    hold on
                    plot(Time_Eye,GyroY,'k--')
                    plot(Time_Eye,GyroZ,'k-')
                    hold off
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
                        plot(Time_Eye,GyroX,'k:')
                        plot(Time_Eye,GyroY,'k--')
                        plot(Time_Eye,GyroZ,'k-')
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
                    plot(Time_Eye,GyroX,'k:')
                    plot(Time_Eye,GyroY,'k--')
                    plot(Time_Eye,GyroZ,'k-')
                    hold off
                    pause(1)  
                else % low freq sine
                    %Just take the whole thing as a segment
                    start = is(q);
                    ends = ie(q);
                end
            elseif contains(exp,'VelStep') % Velstep               
                thresh = 50; %Adjust as needed
                thresh2 = 1; %Counts as 0.
                approx0 = find(abs(GyroAll)< thresh2);
                start1 = find(abs(GyroAll) > thresh);
                start2 = start1([true;diff(start1)>Fs*45]); %make sure it's at least 45s long
                [~,inds] = min(abs(repmat(approx0',length(start2),1) - start2),[],2);
                start = approx0(inds);
                start = start([true;diff(start)>Fs*45]); %make sure it's at least 45s long
                ends = [start(2:end)-1;ie(q)];
                plot(Time_Stim,GyroAll,'k',Time_Stim(start),GyroAll(start),'r*',Time_Stim(ends),GyroAll(ends),'b*')
                xlabel('Time (s)')
                ylabel('Head Velocity')
                legend('Stimulus','Start','End')
                pause(1)
            elseif contains(exp,'Impulse') %Impulse
                if length(stim_info) >1
                    plot(Time_Eye,GyroAll,'k-')
                    uiwait(msgbox('Select the boundaries between impulse velocities.'))
                    [x,~] = ginput(length(stim_info)-1);
                    start = [is(q);round(x,0)];
                    ends = [round(x,0);ie(q)];
                else %Take the whole thing
                    start = is(q);
                    ends = ie(q);
                end
            end
        elseif contains(exp,'eeVOR') %Good for all externally triggered stimuli for now        
            if contains(exp,{'Sine','Frequency'})
                %Every trigger toggle is a cycle
                temp = find(abs(diff(Stim))==1);
                all_starts = temp(1:end-1);
                all_ends = temp(2:end);
                temp2 = diff(temp);
                start = all_starts;
                ends = all_ends;
                thresh_tol = 30; %tolerance for differences in cyc length
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
                start(isnan(start)) = [];
                ends(isnan(ends)) = [];
                plot(Time_Stim,Stim,'k',Time_Stim(round(start,0)),Stim(round(start,0)),'r*',Time_Stim(round(ends,0)),Stim(round(ends,0)),'b*')
            elseif contains(exp,'Activation')
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
                    start = is(q);
                    ends = ie(q);
                end
            else %Pulse Trains / Multi Vector
                %Sometimes the trigger starts high and needs to be fixed
                if Stim(1)==1
                    %Set all leading 1's to 0
                    Stim(1:find(Stim==0,1,'first')) = 0;                    
                end
                inds = find(Stim==1);
                inds([false;diff(inds)==1]) = []; %Multiple points from same trigger toggle
                pad = median(diff(inds));
                start = inds([true;diff(inds)>2*pad])-pad;
                ends = inds([diff(inds)>2*pad;true])+pad;
            end
            stim = Stim;
            plot(Time_Stim,stim,'k',Time_Stim(round(start,0)),stim(round(start,0)),'r*',Time_Stim(round(ends,0)),stim(round(ends,0)),'b*')
            xlabel('Time (s)')
            ylabel('Trigger')
            legend('Stimulus','Start','End')
            pause(1)
        else
            error('Unrecognized experiment type during segmenting.')
        end
        %% Save
        %Remove any traces marked as "trash"
        rm = contains(stim_info,'trash');
        if any(rm)
            stim_info(rm) = [];
            start(rm) = [];
            ends(rm) = [];
        end    
        if length(start) ~= length(ends)
            disp('Something went awry when detecting segments')
            disp('There are unequal numbers of starting and ending points. Please segment the following file manually in VOMA:')
            disp(In_Path)
        elseif length(stim_info) > length(start)
            disp('There are more experiments in the notes file than there are detected segments. Please adjust the notes for the following file:')
            disp(In_Path)
        elseif length(stim_info) < length(start)
            disp('There are more detected segments than experiments in the notes file. Please adjust the notes for the following file:')
            disp(In_Path)
        else %Can automatically segment
            for i = 1:length(start)
                info.dataType = stim_info{i};
                %Fix this to be more general for pulse trains
                if contains(fileinfo{7}(2:end-1),'65Vector') %stim vec in the info file
                    info.stim_axis = str2double(info.dataType(strfind(info.dataType,'['):strfind(info.dataType,']')));
                elseif contains(fileinfo{7}(2:end-1),'RotaryChair') %horizontal only
                    info.stim_axis = [0,0,1]; %Moves to the left first
                elseif contains(fileinfo{7}(2:end-1),{'PulseTrain','Autoscan'})
                    %first two char are the canal
                    exp = stim_info{i}(2:end-1);
                    canal = exp(1:2);
                    switch canal
                        case 'LP'
                            info.stim_axis = [0,-1,0];
                        case 'LH'
                            info.stim_axis = [0,0,1];
                        case 'LA'
                            info.stim_axis = [-1,0,0];
                        case 'RP'
                            info.stim_axis = [1,0,0];
                        case 'RH'
                            info.stim_axis = [0,0,-1];
                        case 'RA'
                            info.stim_axis = [0,1,0];
                    end
                elseif contains(fileinfo{7}(2:end-1),{'Sine','aHIT'}) %probably not quite right for aHIT
                    exp = stim_info{i}(2:end-1);
                    sep2 = strfind(exp,'_');
                    axis = exp(1:sep2-1);
                    switch axis
                        case 'RALP'
                            info.stim_axis = [0,1,0];
                        case 'LHRH'
                            info.stim_axis = [0,0,1];
                        case 'LARP'
                            info.stim_axis = [1,0,0];
                        case 'X'
                            info.stim_axis = [0.707,0.707,0];
                        case 'Y' 
                            info.stim_axis = [-0.707,0.707,0];
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
                Data.rawfile = {info.rawfile};
                %Save but make a new ending if there are multiple segments with
                %the same information
                fname = [info.subject,'-',info.visit,'-',info.exp_date,'-',info.dataType];
                save_flag = 1;
                if exist([Seg_Path,filesep,fname,'.mat'],'file') %Already an instance of this file
                    Data2 = Data; %set the segment to Data 2 to check against current file
                    load([Seg_Path,filesep,fname,'.mat'],'Data')
                    if ~any(ismember(Data.rawfile,Data2.rawfile)&Data.raw_start_t==Data2.raw_start_t)
                        disp([fname,' already exists in this folder and they were combined.'])
                        %Make the time vectors continuous with the first segment
                        Time_Eye2 = Data2.Time_Eye - Data2.Time_Eye(1) + Data.Time_Eye(end) + mean(diff(Data2.Time_Eye));
                        Time_Stim2 = Data2.Time_Stim - Data2.Time_Stim(1) + Data.Time_Stim(end) + mean(diff(Data2.Time_Stim));
                        %Figure out how many copies of info there are 
                        copynum = sum(contains(fieldnames(Data),'info'))+1;
                        %Combine all of the items
                        Data.(['info',num2str(copynum)]) = info; %as many info copies as needed
                        Data.(['Data',num2str(copynum)]) = Data2;
                        Data.Time_Eye = [reshape(Data.Time_Eye,[],1);reshape(Time_Eye2,[],1)];
                        Data.Time_Stim = [reshape(Data.Time_Stim,[],1);reshape(Time_Stim2,[],1)];
                        Data.raw_start_t = [Data.raw_start_t;Data2.raw_start_t];
                        Data.raw_end_t = [Data.raw_end_t;Data.raw_end_t];
                        Data.Trigger = [Data.Trigger;Data2.Trigger];
                        Data.LE_Position_X = [Data.LE_Position_X;Data2.LE_Position_X];
                        Data.LE_Position_Y = [Data.LE_Position_Y;Data2.LE_Position_Y];
                        Data.LE_Position_Z = [Data.LE_Position_Z;Data2.LE_Position_Z];
                        Data.RE_Position_X = [Data.RE_Position_X;Data2.RE_Position_X];
                        Data.RE_Position_Y = [Data.RE_Position_Y;Data2.RE_Position_Y];
                        Data.RE_Position_Z = [Data.RE_Position_Z;Data2.RE_Position_Z];
                        Data.HeadVel_X = [reshape(Data.HeadVel_X,[],1);reshape(Data2.HeadVel_X,[],1)];
                        Data.HeadVel_Y = [reshape(Data.HeadVel_Y,[],1);reshape(Data2.HeadVel_Y,[],1)];
                        Data.HeadVel_Z = [reshape(Data.HeadVel_Z,[],1);reshape(Data2.HeadVel_Z,[],1)];
                        Data.HeadAccel_X = [reshape(Data.HeadAccel_X,[],1);reshape(Data2.HeadAccel_X,[],1)];
                        Data.HeadAccel_Y = [reshape(Data.HeadAccel_Y,[],1);reshape(Data2.HeadAccel_Y,[],1)];
                        Data.HeadAccel_Z = [reshape(Data.HeadAccel_Z,[],1);reshape(Data2.HeadAccel_Z,[],1)];
                        Data.rawfile = [Data.rawfile;Data2.rawfile];
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
%                    fig = plotSegment(Data);
%                    savefig(fig,[Seg_Path,filesep,fname,'.fig'])
                end
            end
        end
    end
end