In_Path = 'LongTests_4_28_2021_12_12_PM_10Custom Test.dat';
TEMP_In_Path = 'LongTests_4_14_2021_12_20_PM_9Custom Test.dat';
% Standardize Colors
load('VNELColors','colors')
%Load template from file
%Load in file
if isempty(TEMP_In_Path)
    TEMP_Time_Eye = NaN;
    TEMP_Stim = NaN;
else
    if strcmp(TEMP_In_Path(end-2:end),'dat')%NKI
        warning('off')
        data = readtable(TEMP_In_Path,'ReadVariableNames',true);
        warning('on')
        data.Properties.VariableNames{1} = 'EyeTime';
        TEMP_Time_Eye = data.EyeTime;        
        TEMP_Stim = zeros(length(TEMP_Time_Eye),1); 
        TEMP_Stim(data.EventCode ~= 0) = 1;   
    elseif strcmp(TEMP_In_Path(end-2:end),'txt') %LDVOG
        data = readtable(TEMP_In_Path);
        % Generate Time_Eye vector
        Time = data{:,2};
        TEMP_Time_Eye = (0:length(Time)-1)'*median(diff(Time));
        % Index for the VOG GPIO line
        StimIndex = 35; 
        TEMP_Stim = data{1:length(TEMP_Time_Eye),StimIndex};
    else
        TEMP_Time_Eye = NaN;
        TEMP_Stim = NaN;
    end
end  
% Load in file
if strcmp(In_Path(end-2:end),'dat')%NKI
    warning('off')
    data = readtable(In_Path,'ReadVariableNames',true);
    warning('on')
    data.Properties.VariableNames{1} = 'EyeTime';
    Time_Eye = data.EyeTime;        
    %Transform coordinates to be in standard canal coordinates (X,Y,Z)
    GX = data.GyroX - median(data.GyroX);
    GY = data.GyroY - median(data.GyroY);
    GZ = data.GyroZ - median(data.GyroZ);
    GyroX = GY; 
    GyroY = -GX; 
    GyroZ = -GZ; 
    Stim = zeros(length(Time_Eye),1); 
    Stim(data.EventCode ~= 0) = 1;
    LZ = -data.LeftHoriz;
    LY = -data.LeftVert;
    LX = data.LeftTorsion;
    RZ = -data.RightHoriz;
    RY = -data.RightVert;
    RX = data.RightTorsion;    
elseif strcmp(In_Path(end-2:end),'txt') %LDVOG
    data = readtable(In_Path);
    % Generate Time_Eye vector
    Time = data{:,2};
    Time_Eye = (0:length(Time)-1)'*median(diff(Time));
    % Index for the VOG GPIO line
    StimIndex = 35; 
    XvelHeadIndex = 30;
    YvelHeadIndex = 29;
    ZvelHeadIndex = 28;
    Stim = data{1:length(Time_Eye),StimIndex};
    %Transform coordinates to be in standard canal coordinates (X,Y,Z)
    GyroX = data{1:length(Time_Eye),XvelHeadIndex};
    GyroY = data{1:length(Time_Eye),YvelHeadIndex};
    GyroZ = -data{1:length(Time_Eye),ZvelHeadIndex};  
    HLeftIndex = 40;
    VLeftIndex = 41;
    TLeftIndex = 42;
    HRightIndex = 43;
    VRightIndex = 44;
    TRightIndex = 45;       
    % Load raw eye position data in Fick coordinates [degrees]
    LZ = data{:,HLeftIndex};
    LY = data{:,VLeftIndex};
    LX = data{:,TLeftIndex};
    RZ = data{:,HRightIndex};
    RY = data{:,VRightIndex};
    RX = data{:,TRightIndex}; 
else
    Time_Eye = NaN;
    Stim = NaN;
end
Fs = 1/mean(diff(Time_Eye));
% Plot
plot(NaN,NaN)
hold on
plot(Time_Eye,LX,'Color',colors.l_x);plot(Time_Eye,RX,'Color',colors.r_x);plot(Time_Eye,LY,'Color',colors.l_y);
plot(Time_Eye,RY,'Color',colors.r_y);plot(Time_Eye,LZ,'Color',colors.l_z);plot(Time_Eye,RZ,'Color',colors.r_z);
plot(Time_Eye,GyroX,'k:',Time_Eye,GyroY,'k--',Time_Eye,GyroZ,'k-')
plot(Time_Eye,Stim,'b'); plot(TEMP_Time_Eye,TEMP_Stim,'g');
hold off
xlabel('Time (s)')
ylabel('Velocity (dps)')
title(strrep(In_Path,'_',' '))
axis([Time_Eye(1) Time_Eye(end) -1 10])
%% Template + Stim toggle: Align to stim toggle start or stop
%Must be at least one toggle at the start or end of the sets for alignment
align = 'end'; %can also be start
if length(TEMP_Stim)>length(Stim)
    template = TEMP_Stim(1:length(Stim));
elseif length(TEMP_Stim)<length(Stim)
    template = [TEMP_Stim;TEMP_Stim(end)*ones(length(Stim)-length(TEMP_Stim),1)];
end
if strcmp(align,'end')
    %Flip the template if needed
    if template(end)~=Stim(end)
        template = 1-template;
    end
    k1 = find(Stim==(1-Stim(end)),1,'last');
    k2 = find(template==(1-template(end)),1,'last');
elseif strcmp(align,'start')
    if template(1)~=Stim(1)
        template = 1-template;
    end
    k1 = find(Stim==(1-Stim(1)),1,'first');
    k2 = find(template==(1-template(1)),1,'first');
end
if k1 > k2
    template = [template(1)*ones(k1-k2,1);template];
elseif k1 < k2
    template = [template((k2-k1):end);template(end)*ones(k2-k1,1)];
end
template = template(1:length(Stim));

% Plot
plot(NaN,NaN)
hold on
plot(Time_Eye,LX,'Color',colors.l_x);plot(Time_Eye,RX,'Color',colors.r_x);plot(Time_Eye,LY,'Color',colors.l_y);
plot(Time_Eye,RY,'Color',colors.r_y);plot(Time_Eye,LZ,'Color',colors.l_z);plot(Time_Eye,RZ,'Color',colors.r_z);
plot(Time_Eye,GyroX,'k:',Time_Eye,GyroY,'k--',Time_Eye,GyroZ,'k-')
plot(Time_Eye,Stim,'b',Time_Eye,template,'g');
hold off
xlabel('Time (s)')
ylabel('Velocity (dps)')
title(strrep(In_Path,'_',' '))
axis([Time_Eye(1) Time_Eye(end) -1 10])
%% Template: Adjust Template Zero Padding
TEMP_Stim2 = TEMP_Stim(find(TEMP_Stim==(1-TEMP_Stim(1)),1,'first')-1:end);
l_i = 855;
r_i = 30000;
template = [TEMP_Stim2(1)*ones(l_i,1);TEMP_Stim2;TEMP_Stim2(end)*ones(r_i,1)];
template = template(1:length(Stim));
%template = 1-template;
% Plot
plot(NaN,NaN)
hold on
plot(Time_Eye,LX,'Color',colors.l_x);plot(Time_Eye,RX,'Color',colors.r_x);plot(Time_Eye,LY,'Color',colors.l_y);
plot(Time_Eye,RY,'Color',colors.r_y);plot(Time_Eye,LZ,'Color',colors.l_z);plot(Time_Eye,RZ,'Color',colors.r_z);
plot(Time_Eye,GyroX,'k:',Time_Eye,GyroY,'k--',Time_Eye,GyroZ,'k-')
plot(Time_Eye,Stim,'b',Time_Eye,template,'g');
hold off
xlabel('Time (s)')
ylabel('Velocity (dps)')
title(strrep(In_Path,'_',' '))
%axis([3 18 -1 10])
axis([Time_Eye(1) Time_Eye(end) -1 10])
%% Use Stim Toggle as its own Template
%Frequency Sweep High to Low
%Set time chunk to look at
shift = -7100;
t1 = 120; %s 
t2 = 148; %s 
[~,t1_ind] = min(abs(Time_Eye-t1));
[~,t2_ind] = min(abs(Time_Eye-t2));
template = Stim;
template(t1_ind:t2_ind) = Stim((t1_ind:t2_ind)+shift);
% Plot
plot(NaN,NaN)
hold on
plot(Time_Eye,LX,'Color',colors.l_x);plot(Time_Eye,RX,'Color',colors.r_x);plot(Time_Eye,LY,'Color',colors.l_y);
plot(Time_Eye,RY,'Color',colors.r_y);plot(Time_Eye,LZ,'Color',colors.l_z);plot(Time_Eye,RZ,'Color',colors.r_z);
plot(Time_Eye,GyroX,'k:',Time_Eye,GyroY,'k--',Time_Eye,GyroZ,'k-')
plot(Time_Eye,Stim,'b',Time_Eye,template,'g');
hold off
xlabel('Time (s)')
ylabel('Velocity (dps)')
title(strrep(In_Path,'_',' '))
axis([Time_Eye(1) Time_Eye(end) -1 10])   
%% Keep changes and Save to file
Stim = template;
new_fname = [In_Path(1:end-4),'_UpdatedTrigger_',datestr(now,'yyyymmdd_HHMMSS'),'_FixedShape',In_Path(end-3:end)];
if strcmp(In_Path(end-2:end),'dat')%NKI
    data.EventCode = Stim;
    writetable(data,new_fname,'Delimiter','tab');
elseif strcmp(In_Path(end-2:end),'txt') %LDVOG
    data{:,StimIndex} = reshape(Stim,[],1);
    writetable(data,new_fname);
end