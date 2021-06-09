function updateRawVOGTrigger(Raw_Path,In_Path,TEMP_In_Path)
if nargin < 1
    Raw_Path = cd;
end
all_files = extractfield(dir(Raw_Path),'name',find(~extractfield(dir(Raw_Path),'isdir')));
LDVOG_files = all_files(contains(all_files,'SESSION')&contains(all_files,'.txt')&~contains(all_files,'Notes.txt'));
NKI_files = all_files(contains(all_files,'.dat'));
VOG_files = [LDVOG_files;NKI_files];
%% Initialize Figure
fig = figure(1);
fig.Units = 'normalized';
fig.Position = [0 0 1 1];
fig.Units = 'inches';
screen_size = fig.Position;
fig.Position = screen_size - [0 0 6 0];
%% Get file names
if nargin < 2
    [ind1,tf1] = nmlistdlg('PromptString','Select an action:',...
                       'SelectionMode','single',...
                       'ListSize',[300 150],...
                       'ListString',VOG_files,...
                       'Position',[screen_size(3)-6,screen_size(4)-3.75,5,3.75]); 
    if ~tf1
        return;
    end
    In_Path = [Raw_Path,filesep,VOG_files{ind1}];
end
if nargin < 3
    [ind2,tf2] = nmlistdlg('PromptString','Select an action:',...
                       'SelectionMode','single',...
                       'ListSize',[300 150],...
                       'ListString',VOG_files,...
                       'Position',[screen_size(3)-6,screen_size(4)-3.75,5,3.75]); 
    if tf2
        TEMP_In_Path = [Raw_Path,filesep,VOG_files{ind2}];
    else
        TEMP_In_Path = [];
    end
end    
%% Load in files
% Standardize Colors
load('VNELColors','colors')  
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
%Load template from file
if isempty(TEMP_In_Path)
    TEMP_Time_Eye = Time_Eye;
    TEMP_Stim = NaN*Stim;
else
    if strcmp(TEMP_In_Path(end-2:end),'dat')%NKI
        warning('off')
        data2 = readtable(TEMP_In_Path,'ReadVariableNames',true);
        warning('on')
        data2.Properties.VariableNames{1} = 'EyeTime';  
        TEMP_Time_Eye = data2.EyeTime;
        TEMP_Stim = zeros(length(TEMP_Time_Eye),1); 
        TEMP_Stim(data2.EventCode ~= 0) = 1; 
    elseif strcmp(TEMP_In_Path(end-2:end),'txt') %LDVOG
        data2 = readtable(TEMP_In_Path);
        % Generate Time_Eye vector
        Time = data2{:,2};
        TEMP_Time_Eye = (0:length(Time)-1)'*median(diff(Time)); 
        % Index for the VOG GPIO line
        StimIndex = 35; 
        TEMP_Stim = data2{1:length(TEMP_Time_Eye),StimIndex};
    else
        TEMP_Time_Eye = Time_Eye;
        TEMP_Stim = NaN*Stim;
    end
end
%% Adjust Trigger
ind = 2; %Start with Initialize
opts = {'Save','Initialize','Set Plot Limits','Set Plot Items','Align From Start','Align From End',' '};
while ~strcmp(opts{ind},'Save')
    if strcmp(opts{ind},'Initialize')
        % Initialize Plot
        template = Stim;
        plot_eyes = 1;
        plot_gyro = 1;
        time1 = Time_Eye(1);
        time2 = TEMP_Time_Eye(1);
        len = min([Time_Eye(end)-Time_Eye(1),TEMP_Time_Eye(end)-TEMP_Time_Eye(1)]);
        [~,t1] = min(abs(Time_Eye-time1));
        [~,t2] = min(abs(Time_Eye-(time1+len)));
        [~,t3] = min(abs(TEMP_Time_Eye-time2));
        [~,t4] = min(abs(TEMP_Time_Eye-(time2+len)));
    elseif strcmp(opts{ind},'Set Plot Limits')
        %Get new parameter values
        prompt = {['Set plot axis limits:',newline,newline,'Sync Time Start:'],...
            'Template Time Start:','Length:'};
        dlgtitle = 'Y-axis Limits';
        definput = cellfun(@(x) num2str(x,10),num2cell([time1,time2,len]),'UniformOutput',false);
        out_nums = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 25],definput,'on',1,[screen_size(3)-6 screen_size(4)-4 3 3]));
        if ~isempty(out_nums)
            time1 = out_nums(1);
            time2 = out_nums(2);
            len = out_nums(3);
        end
        [~,t1] = min(abs(Time_Eye-time1));
        [~,t2] = min(abs(Time_Eye-(time1+len)));
        [~,t3] = min(abs(TEMP_Time_Eye-time2));
        [~,t4] = min(abs(TEMP_Time_Eye-(time2+len)));
        set(gca,'XLim',XLim,'YLim',YLim)
    elseif strcmp(opts{ind},'Set Plot Items')
        plot_opts = {'Plot Eyes','Plot Gyro','Plot Sync'};
        [ind3,tf3] = nmlistdlg('PromptString','Select an action:',...
                           'SelectionMode','multiple',...
                           'ListSize',[150 150],...
                           'ListString',plot_opts,...
                           'Position',[screen_size(3)-6,screen_size(4)-3.75,3,3.75],...
                           'InitialValue',find([plot_eyes,plot_gyro,1])); 
        if tf3
            if ismember(ind3,1)
                plot_eyes = 1;
            else
                plot_eyes = 0;
            end
            if ismember(ind3,2)
                plot_gyro = 1;
            else
                plot_gyro = 0;
            end
        end    
    elseif strcmp(opts{ind},'Align From Start')
        %Flip the template if needed
        template(t1:t2) = TEMP_Stim(t3:t4);
        if template(t1)~=Stim(t1)
            template(t1:t2) = 1-template(t1:t2);
        end
        k1 = find(Stim(t1:t2)==(1-Stim(t1)),1,'first')-1+t1;
        k2 = find(template(t1:t2)==(1-template(t1)),1,'first')-1+t1;
        if k1 > k2
            template1 = [template(t1)*ones(k1-k2,1);template];
        elseif k1 < k2
            template1 = [template((k2-k1):end);template(t2)*ones(k2-k1,1)];
        else
            template1 = template;
        end
        template(t1:t2) = template1(1:length(t1:t2));
    elseif strcmp(opts{ind},'Align From End')
        %Flip the template if needed
        template(t1:t2) = TEMP_Stim(t3:t4);
        if template(t2)~=Stim(t2)
            template(t1:t2) = 1-template(t1:t2);
        end
        k1 = find(Stim(t1:t2)==(1-Stim(t2)),1,'last');
        k2 = find(template(t1:t2)==(1-template(t2)),1,'last');
        if k1 > k2
            template1 = [template(t1)*ones(k1-k2,1);template];
        elseif k1 < k2
            template1 = [template((k2-k1):end);template(t2)*ones(k2-k1,1)];
        else
            template1 = template;
        end
        template(t1:t2) = template1(1:length(t1:t2));
    end
    %Remake Plot
    plot(NaN,NaN)
    hold on
    if plot_eyes
        plot(Time_Eye,LX,'Color',colors.l_x);plot(Time_Eye,RX,'Color',colors.r_x);plot(Time_Eye,LY,'Color',colors.l_y);
        plot(Time_Eye,RY,'Color',colors.r_y);plot(Time_Eye,LZ,'Color',colors.l_z);plot(Time_Eye,RZ,'Color',colors.r_z);
    end
    if plot_gyro
        plot(Time_Eye,GyroX,'k:',Time_Eye,GyroY,'k--',Time_Eye,GyroZ,'k-')
    end
    h(1) = plot(Time_Eye,Stim,'b');
    h(2) = plot(TEMP_Time_Eye,TEMP_Stim-0.25,'g');
    h(3) = plot(Time_Eye,template+0.25,'k');
    plot(Time_Eye(t1),Stim(t1),'b*',Time_Eye(t2),Stim(t2),'b*')
    plot(TEMP_Time_Eye(t3),TEMP_Stim(t3)-0.25,'g*',TEMP_Time_Eye(t4),TEMP_Stim(t4)-0.25,'g*')
    hold off
    xlabel('Time (s)')
    ylabel('Velocity (dps)')
    legend(h,'Original Sync','Template Sync','New Sync')
    %Get next step
    [ind,tf] = nmlistdlg('PromptString','Select an action:',...
                           'SelectionMode','single',...
                           'ListSize',[150 150],...
                           'ListString',opts,...
                           'Position',[screen_size(3)-6,screen_size(4)-3.75,3,3.75]); 
    if ~tf
        return;
    end
end
%% Keep changes and Save to file
new_fname = [In_Path(1:end-4),'_UpdatedTrigger_',datestr(now,'yyyymmdd_HHMMSS'),'_FixedShape',In_Path(end-3:end)];
if strcmp(In_Path(end-2:end),'dat')%NKI
    data.EventCode = template;
    writetable(data,new_fname,'Delimiter','tab');
    if ~any(contains(extractfield(dir(Raw_Path),'name'),'trash')&extractfield(dir(Raw_Path),'isdir'))
        mkdir([Raw_Path,filesep,'trash'])
    end
    movefile(In_Path,[Raw_Path,filesep,'trash'])
elseif strcmp(In_Path(end-2:end),'txt') %LDVOG
    data{:,StimIndex} = reshape(template,[],1);
    writetable(data,new_fname);
    if ~any(contains(extractfield(dir(Raw_Path),'name'),'trash')&extractfield(dir(Raw_Path),'isdir'))
        mkdir([Raw_Path,filesep,'trash'])
    end
    movefile(In_Path,[Raw_Path,filesep,'trash'])
end 
end