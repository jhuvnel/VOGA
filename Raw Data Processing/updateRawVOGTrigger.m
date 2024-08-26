function updateRawVOGTrigger(Raw_Path,In_Path,TEMP_In_Path)
if nargin < 1
    Raw_Path = cd;
    if isfolder([cd,filesep,'Raw Files'])
        Raw_Path = [cd,filesep,'Raw Files'];
    end
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
    [ind1,tf1] = nmlistdlg('PromptString','Select a file to fix:','SelectionMode','single',...
               'ListSize',[300 150],'ListString',VOG_files,'Position',[screen_size(3)-6,screen_size(4)-3.75,5,3.75]); 
    if ~tf1
        return;
    end
    In_Path = [Raw_Path,filesep,VOG_files{ind1}];
elseif ~contains(In_Path,Raw_Path)
    In_Path = [Raw_Path,filesep,In_Path];
end
if nargin < 3
    [ind2,tf2] = nmlistdlg('PromptString','Select a template file:','SelectionMode','single','ListSize',[300 150],...
                       'ListString',VOG_files,'Position',[screen_size(3)-6,screen_size(4)-3.75,5,3.75]); 
    TEMP_In_Path = [];
    if tf2
        TEMP_In_Path = [Raw_Path,filesep,VOG_files{ind2}];
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
   TEMP_In_Path = ''; 
end    
if strcmp(TEMP_In_Path(end-2:end),'dat')%NKI
    warning('off')
    data2 = readtable(TEMP_In_Path,'ReadVariableNames',true);
    warning('on')
    TEMP_Stim = zeros(length(data2.EventCode),1);
    TEMP_Stim(data2.EventCode ~= 0) = 1;    
    small_TEMP = TEMP_Stim(find(TEMP_Stim~=TEMP_Stim(1),1,'first')-1:find(TEMP_Stim~=TEMP_Stim(end),1,'last')+1);
    TEMP_Stim = 0*Stim;
    TEMP_Stim(1:length(small_TEMP)) = small_TEMP;   
elseif strcmp(TEMP_In_Path(end-2:end),'txt') %LDVOG
    data2 = readtable(TEMP_In_Path);
    TEMP_Stim = data2{1:size(data2,1),35};
    small_TEMP = TEMP_Stim(find(TEMP_Stim~=TEMP_Stim(1),1,'first')-1:find(TEMP_Stim~=TEMP_Stim(end),1,'last')+1);
    TEMP_Stim = 0*Stim;
    TEMP_Stim(1:length(small_TEMP)) = small_TEMP;     
else
    TEMP_Stim = NaN*Stim;
end
%% Adjust Trigger
ind = 2; %Start with Initialize
opts = {'Save','Initialize','Set Plot Limits','Set Plot Items','Align From Start','Align From End','Zero Padding'};
while ~strcmp(opts{ind},'Save')
    if strcmp(opts{ind},'Initialize')
        % Initialize Plot        
        plot_eyes = 1;
        plot_gyro = 1;
        rel_shift = 0;
        template = TEMP_Stim;
        plot(NaN,NaN)
        hold on
        plot(Time_Eye',LX','Color',colors.l_x);plot(Time_Eye',RX','Color',colors.r_x);
        plot(Time_Eye,LY,'Color',colors.l_y);plot(Time_Eye,RY,'Color',colors.r_y);
        plot(Time_Eye,LZ,'Color',colors.l_z);plot(Time_Eye,RZ,'Color',colors.r_z);
        plot(Time_Eye,GyroX,'k:',Time_Eye,GyroY,'k--',Time_Eye,GyroZ,'k-')
        h1 = plot(Time_Eye,Stim,'b');
        h2 = plot(Time_Eye,template,'g');
        hold off
        xlabel('Time (s)')
        ylabel('Velocity (dps)')
        legend([h1;h2],'Original Sync','New Sync')        
    elseif strcmp(opts{ind},'Set Plot Items')
        plot_opts = {'Plot Eyes','Plot Gyro'};
        [ind3,tf3] = nmlistdlg('PromptString','Select an action:',...
                           'SelectionMode','multiple',...
                           'ListSize',[150 150],...
                           'ListString',plot_opts,...
                           'Position',[screen_size(3)-6,screen_size(4)-3.75,3,3.75],...
                           'InitialValue',find([plot_eyes,plot_gyro])); 
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
        if template(1)~=Stim(1)
            template = 1-template;
        end
        k1 = find(Stim==(1-Stim(1)),1,'first');
        k2 = find(template==(1-template(1)),1,'first');
        if k1 > k2
            template1 = [template(1)*ones(k1-k2,1);template];
        elseif k1 < k2
            template1 = [template((k2-k1):end);template(1)*ones(k2-k1,1)];
        else
            template1 = template;
        end
        template = template1(1:length(template));
    elseif strcmp(opts{ind},'Align From End')
        %Flip the template if needed
        if template(end)~=Stim(end)
            template = 1-template;
        end
        k1 = find(Stim==(1-Stim(end)),1,'last');
        k2 = find(template==(1-template(end)),1,'last');
        if k1 > k2
            template1 = [template(1)*ones(k1-k2,1);template];
        elseif k1 < k2
            template1 = [template((k2-k1):end);template(1)*ones(k2-k1,1)];
        else
            template1 = template;
        end
        template = template1(1:length(template));
    elseif strcmp(opts{ind},'Zero Padding')
        small_template = TEMP_Stim(find(TEMP_Stim~=TEMP_Stim(1),1,'first')-1:find(TEMP_Stim~=TEMP_Stim(end),1,'last')+1);
        len = length(small_template);
        temp = inputdlg('Zero padding:','Zero padding',1,{num2str(rel_shift)});
        if ~isempty(temp)&&~isnan(str2double(temp{:}))
            rel_shift = str2double(temp{:});
        end
        template = 0*template;
        template((1:len)+rel_shift) = small_template;
    end        
    XLim = get(gca,'XLim');
    YLim = get(gca,'YLim'); 
    %Plot
    plot(NaN,NaN)
    hold on
    if plot_eyes
        plot(Time_Eye',LX','Color',colors.l_x);plot(Time_Eye',RX','Color',colors.r_x);
        plot(Time_Eye,LY,'Color',colors.l_y);plot(Time_Eye,RY,'Color',colors.r_y);
        plot(Time_Eye,LZ,'Color',colors.l_z);plot(Time_Eye,RZ,'Color',colors.r_z);
    end
    if plot_gyro
        plot(Time_Eye,GyroX,'k:',Time_Eye,GyroY,'k--',Time_Eye,GyroZ,'k-')
    end
    h1 = plot(Time_Eye,Stim,'b');
    h2 = plot(Time_Eye,template,'g');
    hold off
    xlabel('Time (s)')
    ylabel('Velocity (dps)')
    legend([h1;h2],'Original Sync','New Sync')
    set(gca,'XLim',XLim','YLim',YLim)
    %Get next step
    [ind,tf] = nmlistdlg('PromptString','Select an action:','SelectionMode','single','ListSize',[150 150],...
                           'ListString',opts,'Position',[screen_size(3)-6,screen_size(4)-3.75,3,3.75]); 
    if ~tf
        return;
    end
end
%% Keep changes and Save to file
new_fname = [In_Path(1:end-4),'_UpdatedTrigger_',char(datetime("now",'Format','yyyyMMdd_HHmmSS')),'_FixedShape',In_Path(end-3:end)];
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