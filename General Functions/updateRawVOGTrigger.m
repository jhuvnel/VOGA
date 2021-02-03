In_Path = 'SESSION-2020Nov18-145320.txt';
% Standardize Colors
    % Normal colors
    colors.l_x = [237,150,33]/255;
    colors.l_y = [125,46,143]/255;
    colors.l_z = [1 0 0];
    %colors.l_l = [0,128,0]/255;
    %colors.l_r = [0 0 1];
    colors.r_x = [237,204,33]/255;
    colors.r_y = [125,46,230]/255;
    colors.r_z = [1,0,1];
    %colors.r_l = [0 1 0];
    %colors.r_r = [64,224,208]/255;
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
end
%% Plot
plot(NaN,NaN)
hold on
plot(Time_Eye,LX,'Color',colors.l_x)
plot(Time_Eye,RX,'Color',colors.r_x)
plot(Time_Eye,LY,'Color',colors.l_y)
plot(Time_Eye,RY,'Color',colors.r_y)
plot(Time_Eye,LZ,'Color',colors.l_z)
plot(Time_Eye,RZ,'Color',colors.r_z)
plot(Time_Eye,GyroX,'k:')
plot(Time_Eye,GyroY,'k--')
plot(Time_Eye,GyroZ,'k-')
plot(Time_Eye,100*Stim,'b')
hold off
xlabel('Time (s)')
ylabel('Velocity (dps)')
legend('GyroX','GyroY','GyroZ','Trigger')
axis([Time_Eye(1) Time_Eye(end) -10 110])
Fs = 1/mean(diff(Time_Eye));
%% OLD Code that worked before
%         % CODED FOR SOME NOISY TRIGGER EXPERIMENTS
%         plot(100*StimAll,'b')
%         noisy = questdlg(['Bad trigger?',newline,In_Path],'','Fine','Noisy','Fine');
%         if strcmp(noisy,'Noisy')
%             if all(StimAll(1:10)) %Just flipped so starts at 0
%                 Stim2 = -StimAll+1; 
%             else
%                 Stim2 = StimAll;
%             end   
%             min_time = 0.1;  %Lets say 100ms is the smallest time frame to expect a change (5Hz)
%             sus_len = floor(min_time*Fs); % minimum number of expected consecutive zeros
%             zero_vals = find(Stim2==0);
%             k = 2;
%             while(k<length(zero_vals))
%                 if zero_vals(k+1)-1 == zero_vals(k)
%                     zero_vals(k) = [];
%                 else
%                     k = k+1;
%                 end
%             end
%             start1 = zero_vals(1:end-1);
%             end1 = zero_vals(2:end);
%             difference = end1 - start1;
%             start1(difference < sus_len) = [];
%             end1(difference < sus_len) = [];
%             Stim3 = ones(length(Stim2),1);
%             for i = 1:length(start1)
%                 Stim3(start1(i):end1(i)) = 0;
%             end
%             StimAll = Stim3;
%         end
%% Templates
%Frequency Sweep High to Low
t1 = 405;
t2 = 440;
cyc_num = 10;
[~,t1_ind] = min(abs(Time_Eye-t1));
[~,t2_ind] = min(abs(Time_Eye-t2));
spike = diff(find(abs(diff(Stim(t1_ind:t2_ind)))>0));
spike(spike==1) = [];
tol = sum(abs((spike/min(spike)-floor(spike/min(spike)))));
if tol < 0.1
    val = min(spike);
else
    val = spike(end);
end
val = floor(Fs/0.1)-62;
changes = [0;val*ones(cyc_num,1)];
k = find(Stim==(1-Stim(t1_ind))&Time_Eye>t1,1,'first');
template = Stim(t1_ind)*ones(1,length(Time));
%changes = [0;...
           %floor(Fs/5)*ones(20,1);floor(Fs*5);...
           %floor(Fs/2)*ones(20,1);floor(Fs*5);...
           %floor(Fs/1)*ones(20,1);floor(Fs*5);...
           %floor(Fs/0.5)*ones(20,1);floor(Fs*5)];...
           %floor(Fs/0.2)*ones(15,1);floor(Fs*10);...
           %floor(Fs/0.1)*ones(10,1)];
for i = 1:length(changes)
    template(k+changes(i):end) = 1-template(k+changes(i):end);
    k = k+changes(i)+1;
end
plot(Time_Eye,GyroX,'k:')
hold on
plot(Time_Eye,LX,'Color',colors.l_x)
plot(Time_Eye,RX,'Color',colors.r_x)
plot(Time_Eye,LY,'Color',colors.l_y)
plot(Time_Eye,RY,'Color',colors.r_y)
plot(Time_Eye,LZ,'Color',colors.l_z)
plot(Time_Eye,RZ,'Color',colors.r_z)
plot(Time_Eye,GyroY,'k--')
plot(Time_Eye,GyroZ,'k-')
plot(Time_Eye,100*Stim,'b')
plot(Time_Eye,100*template,'g')
hold off
axis([t1 t2 -10 110])        
%% Keep changes
Stim(t1_ind:t2_ind) = template(t1_ind:t2_ind);
%Stim(t1_ind:end) = template(t1_ind:end);
Stim(t2_ind:end) = Stim(t2_ind);
%% Save to file
new_fname = [In_Path(1:end-4),'_UpdatedTrigger_',datestr(now,'yyyymmdd_HHMMSS'),'_FixedShape',In_Path(end-3:end)];
if strcmp(In_Path(end-2:end),'dat')%NKI
    data.Stim = Stim;
    %Test the variable names here
elseif strcmp(In_Path(end-2:end),'txt') %LDVOG
    data{:,StimIndex} = Stim;
    writetable(data,new_fname);
end