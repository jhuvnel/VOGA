%% Make Notes
% RotaryChair-Sine-Condition-Axis-Freq-Speed
% RotaryChair-VelStep-Condition-Axis-Speed
% RotaryChair-SumSine-Condition-Axis-Freq1-Freq2-Freq3-Speed
% aHIT-Sine-Condition-Axis-Freq-Speed
% aHIT-Impulses-Condition-Axis-Speed
% aHIT-Gaussian-Condition-Axis-Speed
% Manual-Impulses-Condition-Axis-Speed
% eeVOR-Sine-Axis-Freq-Speed
% eeVOR-PulseTrain-PFM/PAM-Axis-pps-uA
% eeVOR-MultiVector-DOMAxis
% eeVOR-Autoscan-Canal/Electrode-pps-us-uA
% eeVOR-Activation-Light/Dark#

function flag = MakeNotes(Raw_Path)
%% Input
if nargin < 1
    if isfolder([cd,filesep,'Raw Files'])
        Raw_Path = [cd,filesep,'Raw Files'];
    else
        Raw_Path = cd;
    end
end
%% Process Log Files first
logtoNotes(Raw_Path)
%% Find Files to Make Notes
VOG_fname_pat = {'SESSION','Lateral.txt','LARP.txt','RALP.txt',...
    '.dat','.mat','ImuData'};
rel_dir = dir(Raw_Path);
rel_dir(extractfield(rel_dir,'isdir')) = [];
file_names = extractfield(rel_dir,'name');
file_date = extractfield(rel_dir,'date');
if isempty(file_names)
    file_names = '';
    file_date = [];
end
Notes_ind = contains(file_names,'-Notes.txt');
VOG_ind = find(contains(file_names,VOG_fname_pat)&~Notes_ind&~contains(file_names,{'Raw','.cal'})); %Raw = LDVOG calibration file
has_notes = contains(file_names(VOG_ind),strrep(file_names(Notes_ind),'-Notes.txt',''));
VOG_files = file_names(VOG_ind(~has_notes));
VOG_files_date = file_date(VOG_ind(~has_notes));
if ~any(VOG_ind)
    flag = ['No VOG files (LDVOG/NL/GNO/ESC) have been detected: ',Raw_Path];
elseif isempty(VOG_files)
    flag = ['All VOG Files have Notes files: ',Raw_Path];
else
    flag = '';
end
if ~isempty(flag)
    disp(flag)
    return;
end
%% Set some parameters that will likely stay the same but can be edited
path_parts = strsplit(strrep(strrep(Raw_Path,'_',''),' ',''),filesep);
if any(contains(path_parts,'MVI')&contains(path_parts,'R')) %subject in expected formatting
    sub = path_parts{contains(path_parts,'MVI')&contains(path_parts,'R')};
    MVI_num = str2double(sub((-3:1:-1)+strfind(sub,'R')));
    if ismember(MVI_num,[1,2,3,4,7,9,11,13]) %EXPAND as needed
        ear = 'L';
    else
        ear = 'R';
    end
else
    sub = '';
    ear = '';
end
if any(contains(path_parts,'Visit'))
    vis = path_parts{contains(path_parts,'Visit')};
else
    vis = '';
end
if isempty(sub)||isempty(ear)||isempty(vis)
    common_notes = inputdlg({'Subject:','Ear:','Visit:'},'Set VOG File Parameters',[1,40],{sub,ear,vis}); 
    if isempty(common_notes)
        return;
    end
    sub = common_notes{1};
    ear = common_notes{2};
    vis = common_notes{3};
end
%% Partition by File Type
if any(contains(VOG_files,{'Lateral.txt','LARP.txt','RALP.txt'})) %GNO    
    ang = '0';
    GNO_files = VOG_files(contains(VOG_files,{'Lateral.txt','LARP.txt','RALP.txt'}));
    GNO_SerialNums = {'129852','5628','989477','3012563'}; %The order for GNO1/2/3/4
    for i = 1:length(GNO_files) 
        fname = GNO_files{i};
        dashes = strfind(fname,'_');
        date = fname;
        date(dashes(5)) = '-';
        date = strrep(date(dashes(2)+1:dashes(8)-1),'_','');
        if contains(fname,'Lateral')
            canal = 'LHRH';
        elseif contains(fname,'LARP')
            canal = 'LARP';
        elseif contains(fname,'RALP')
            canal = 'RALP';
        else
            canal = '';
        end
        %Defaults
        type = 'Manual';
        cond = '';
        %Try to plot the accepted head traces (if none, plot the whole time
        %trace)
        try
            GNO_CSV = readtable([Raw_Path,filesep,fname(1:end-4),'.csv'],'ReadVariableNames',false);
            h_ind = find(contains(GNO_CSV{:,1},'Head'));
            left_imp = find(contains(GNO_CSV{:,1},'Impulse')&contains(GNO_CSV{:,2},'L'));
            left_ind = NaN(length(left_imp),1);
            for ii = 1:length(left_imp)
                ind = find((h_ind-left_imp(ii))>0,1,'first');
                left_ind(ii) = h_ind(ind);
            end
            right_imp = find(contains(GNO_CSV{:,1},'Impulse')&contains(GNO_CSV{:,2},'R'));
            right_ind = NaN(length(right_imp),1);
            for ii = 1:length(right_imp)
                ind = find((h_ind-right_imp(ii))>0,1,'first');
                right_ind(ii) = h_ind(ind);
            end        
            detected_left = str2double(split(strrep(GNO_CSV{left_ind,1},',Head,',''),','))';
            detected_right = str2double(split(strrep(GNO_CSV{right_ind,1},',Head,',''),','))';
            t = (0:size(detected_left)-1)*1000/246;
            figure;
            plot(t,detected_left,'k',t,detected_right,'k')            
        catch
            plotRawVOG([Raw_Path,filesep,fname])
        end
        %Try to load the xml file and automatically make the notes
        try
            fdata = cellstr(readlines([Raw_Path,filesep,fname(1:end-4),'.xml']));
            exp_info = strrep(extractXMLdataline(fdata{contains(fdata,'<Remarks>')}),' ','');
            if contains(lower(exp_info),'ahit')
                type = 'aHIT';
            elseif contains(lower(exp_info),'chair')
                type = 'RotaryChair';
            end
            if contains(lower(exp_info),{'off','nostim','preop','preact','pre-op'})||contains(Raw_Path,'Visit 0')
                cond = 'NoStim';
            elseif contains(lower(exp_info),{'constant','baseline'})
                cond = 'ConstantRate';
            elseif contains(lower(exp_info),{'motionmod','mod','accel'})
                cond = 'MotionMod';
            end
            gog_line = fdata{contains(fdata,'GogglesSN')};
            gog_num = gog_line(ismember(gog_line,'0123456789'));
            gog_detect = find(ismember(GNO_SerialNums,gog_num));
            gog = ['GNO',num2str(gog_detect)];            
        catch
            exp_info = '(No File)';
            gog = 'GNO';
        end
        exp = [type,'-Impulse-',cond,'-',canal,'-150dps'];
        w_notes = {['Subject ',sub];['Ear ',ear];['Visit ',vis];['Date ',date];['Goggle ',gog];['Angle ',ang];['Experiment ',exp]};
        prompt = ['File Name: ',fname,newline,'Experimenter Remarks: ',exp_info,newline,newline,'Check Notes: ',newline,'Ex: aHIT-Impulse-NoStim-LHRH-150dps-pseudorandom'];
        notes_check = inputdlg(prompt,'Set VOG File Parameters',[length(w_notes),80],{strjoin(w_notes,'\n')});
        if ~isempty(notes_check)
            w_notes = cellstr(notes_check{1,1});
            filePh = fopen([Raw_Path,filesep,fname(1:end-4),'-Notes.txt'],'w');
            fprintf(filePh,'%s\n',w_notes{:});
            fclose(filePh);
        end
        close;
    end 
elseif contains(Raw_Path,'ESC') %ESC
    ang = '0';
    ESC_files = find(contains(VOG_files,{'.mat','EyePositionData'}));
    for i = 1:length(ESC_files) 
        fname = VOG_files{ESC_files(i)};
        exp_info = strrep(strrep(strrep(fname,'.',' '),'_',' '),'-',' ');
        date_pat1 = digitsPattern(4)+' '+digitsPattern(2)+' '+digitsPattern(2)+' '+digitsPattern(2)+' '+digitsPattern(2)+' '+digitsPattern(2);
        if any(strfind(exp_info,date_pat1)) %In the expected yyyy-MM-dd_HH.mm.ss file format in file name
            VOG_time = datetime(cell2mat(extract(exp_info,date_pat1)),'InputFormat','yyyy MM dd HH mm ss');
        else
            VOG_time = datetime(VOG_files_date{ESC_files(i)}); %use file creation/saving time since there are variable naming standards
            VOG_time.Format = 'yyyy-MM-dd HH:mm:ss';
        end
        if VOG_time < datetime(2018,1,17,0,0,0) 
            gog = 'ESC1'; %Yuri/Schubert pair? Hard to tell
        elseif VOG_time < datetime(2022,11,2,0,0,0)
            gog = 'ESC2';  %Kheradmand pair 
        else 
            gog = 'ESC3';  %MVI pair
        end        
        date = datestr(VOG_time,'yyyymmdd-HHMMss');
        if contains(exp_info,{'Lateral','LHRH'})
            canal = 'LHRH';
        elseif contains(exp_info,'LARP')
            canal = 'LARP';
        elseif contains(exp_info,'RALP')
            canal = 'RALP';
        else
            canal = '';
        end
        if contains(lower(exp_info),'ahit')
            type = 'aHIT';
        elseif contains(lower(exp_info),'chair')
            type = 'RotaryChair';
        else
            type = 'Manual';
        end
        if contains(lower(exp_info),{'off','nostim','preop','preact','pre-op'})||contains(Raw_Path,'Visit 0')
            cond = 'NoStim';
        elseif contains(lower(exp_info),{'constant','baseline'})
            cond = 'ConstantRate';
        elseif contains(lower(exp_info),{'motionmod','mod','accel'})
            cond = 'MotionMod';
        else
            cond = '';
        end
        exp = [type,'-Impulse-',cond,'-',canal,'-150dps'];
        w_notes = {['Subject ',sub];['Ear ',ear];['Visit ',vis];['Date ',date];['Goggle ',gog];['Angle ',ang];['Experiment ',exp]};
        notes_check = inputdlg([fname,newline,newline,'Check Notes: ',newline,'Ex: aHIT-Impulse-NoStim-LHRH-150dps-pseudorandom'],'Set VOG File Parameters',[length(w_notes),80],{strjoin(w_notes,'\n')});
        if ~isempty(notes_check)
            w_notes = cellstr(notes_check{1,1});
            filePh = fopen([Raw_Path,filesep,fname(1:end-4),'-Notes.txt'],'w');
            fprintf(filePh,'%s\n',w_notes{:});
            fclose(filePh);
        end
    end 
else %LDVOG and NKI
    %% Initialize Figure
    fig = figure(1);
    set(fig,'Color','w','Units','normalized','Position',[0,0,1,1]);
    clf; %in case there are leftover anotations
    fig.Units = 'inches';
    screen_size = fig.Position;
    fig.Position = screen_size - [0 0 6 0];
    ax = subplot(1,1,1);
    xlabel('Time (s)')
    ylabel('Velocity (dps)')
    ax.Position = [0.05 0.1 0.6 0.83];
    notes = annotation('textbox',[0.7 0.1 0.25 0.83],'String','',...
        'FontSize',11,'HorizontalAlignment','left','EdgeColor','none');
    %% Check each VOG file
    for i = 1:length(VOG_files)
        fname = VOG_files{i};
        set(notes,'String','')
        if contains(fname,'SESSION') %LDVOG
            VOG_data = readtable([Raw_Path,filesep,fname]);
            %Make date and other labels
            parts = split(strrep(fname,'_','-'),'-');
            if isduration(VOG_data{end,end})&&~isnan(VOG_data{end,end})
                %Use the timestamps on the file itself
                VOG_times = datetime(strrep(parts{2},'.txt',''),'InputFormat','yyyyMMMdd')+[VOG_data{1,end},VOG_data{end,end}];
            else %Use file creation time
                VOG_times = [datetime(strrep([parts{2},' ',parts{3}],'.txt',''),'InputFormat','yyyyMMMdd HHmmss'),datetime(VOG_files_date{i})];
            end
            VOG_times.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
            date = datestr(VOG_times(1),'yyyymmdd-HHMMss');
            gog = 'LDVOG2';
            ang = '-170';
            %Load items for plotting
            % Generate Time_Eye vector
            Time = VOG_data{:,2};
            Time_Eye = (0:length(Time)-1)'*median(diff(Time));
            % Index for the VOG GPIO line
            StimIndex = 35; 
            XvelHeadIndex = 30;
            YvelHeadIndex = 29;
            ZvelHeadIndex = 28;
            Stim = VOG_data{1:length(Time_Eye),StimIndex};
            %Transform coordinates to be in standard canal coordinates (X,Y,Z)
            GyroX = VOG_data{1:length(Time_Eye),XvelHeadIndex};
            GyroY = VOG_data{1:length(Time_Eye),YvelHeadIndex};
            GyroZ = -VOG_data{1:length(Time_Eye),ZvelHeadIndex};
        elseif contains(fname,'.dat') %NKI
            %Load file
            warning('off')
            VOG_data = readtable([Raw_Path,filesep,fname],'ReadVariableNames',true);
            warning('on')
            VOG_data.Properties.VariableNames{1} = 'EyeTime';
            %Make date and other labels
            VOG_time = datetime(VOG_files_date{i}); %use file creation/saving time since NKI doesn't output time stamps
            VOG_time.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
            VOG_times = [VOG_time-seconds(VOG_data{end,1}) VOG_time];
            date = datestr(VOG_times(1),'yyyymmdd-HHMMss');
            if VOG_times(1) > datetime(2022,07,28)
                gog = 'NL2'; %Changed from NKI1 to NL2 on 2022/07/28
            else
                gog = 'NL1';
            end
            ang = '0';
            %Load items for plotting
            Time_Eye = VOG_data.EyeTime;        
            Stim = zeros(length(Time_Eye),1); 
            Stim(VOG_data.EventCode ~= 0) = 1;
            GyroX = VOG_data.GyroY - median(VOG_data.GyroY); 
            GyroY = -(VOG_data.GyroX - median(VOG_data.GyroX)); 
            GyroZ = -(VOG_data.GyroZ - median(VOG_data.GyroZ)); 
        else%Ignore unknown file type
            %ADD CODE here
        end
        if contains(Raw_Path,'Rotary') %All the normal experiments
            w_notes = {['Subject ',sub];['Ear ',ear];['Visit ',vis];['Date ',date];['Goggle ',gog];['Angle ',ang];...
                'Experiment RotaryChair-VelStep-NoStim-LH-240dps';...
                'Experiment RotaryChair-VelStep-NoStim-RH-240dps';...
                'Experiment RotaryChair-Sine-NoStim-LHRH-0.05Hz-100dps';...
                'Experiment RotaryChair-Sine-NoStim-LHRH-0.1Hz-100dps';...
                'Experiment RotaryChair-Sine-NoStim-LHRH-0.2Hz-100dps';...
                'Experiment RotaryChair-Sine-NoStim-LHRH-0.5Hz-100dps';...
                'Experiment RotaryChair-Sine-NoStim-LHRH-1Hz-100dps';...
                };
        elseif contains(Raw_Path,'aHIT') %All the normal experiments
            w_notes = {['Subject ',sub];['Ear ',ear];['Visit ',vis];['Date ',date];['Goggle ',gog];['Angle ',ang];...
                'Experiment aHIT-Gaussian-LightNoStim-LHRH-150dps';...
                'Experiment aHIT-Sine-LightNoStim-LHRH-0.5Hz-35dps';...
                'Experiment aHIT-Sine-LightNoStim-LHRH-1Hz-70dps';...
                'Experiment aHIT-Sine-LightNoStim-LHRH-2Hz-140dps';...
                'Experiment aHIT-Impulse-LightNoStim-LHRH-150dps';...
                };
        else
            w_notes = {['Subject ',sub];['Ear ',ear];['Visit ',vis];['Date ',date];['Goggle ',gog];['Angle ',ang];'Experiment '};
        end
        sub_i = unique(round(linspace(1,length(Time_Eye),1e6)));
        %Update plots
        plot(Time_Eye(sub_i),100*Stim(sub_i),'Color','b')
        hold on
        plot(Time_Eye(sub_i),GyroX(sub_i),'Color',0.5*[1,1,1])
        plot(Time_Eye(sub_i),GyroY(sub_i),'Color',0.75*[1,1,1])
        plot(Time_Eye(sub_i),GyroZ(sub_i),'Color','k')
        hold off
        title([{Raw_Path};{fname}],'FontSize',12,'interpreter','none')
        legend('GyroX','GyroY','GyroZ','Trigger')
        set(notes,'String',w_notes)
        %Action options
        opts = {'Save','Edit Notes','Fix Trigger','Skip'};
        [ind,tf] = nmlistdlg('PromptString','Select an action:',...
               'SelectionMode','single','ListSize',[100 70],'ListString',opts,...
               'Position',[screen_size(3)-6,screen_size(4)-3,1.5,1.75]); 
        if ~tf
            return;
        end
        while ~strcmp(opts{ind},'Save')&&~strcmp(opts{ind},'Skip')
            if strcmp(opts{ind},'Edit Notes')
                notes_check = inputdlg(['Check Notes: ',newline,'Ex: RotaryChair-Sine-NoStim-LHRH-0.05Hz-100dps'],'Set VOG File Parameters',[length(w_notes),70],{strjoin(w_notes,'\n')}); 
                if ~isempty(notes_check)
                    w_notes = cellstr(notes_check{1,1});
                    set(notes,'String',w_notes)
                else
                    return;
                end
            elseif strcmp(opts{ind},'Fix Trigger')
                updateRawVOGTrigger(Raw_Path,fname);
            end
            [ind,tf] = nmlistdlg('PromptString','Select an action:',...
               'SelectionMode','single','ListSize',[100 70],'ListString',opts,...
               'Position',[screen_size(3)-6,screen_size(4)-3,1.5,1.75]); 
            if ~tf
                return;
            end
        end
        if strcmp(opts{ind},'Save')
            %If you got here it's time to save
            filePh = fopen([Raw_Path,filesep,fname(1:end-4),'-Notes.txt'],'w');
            fprintf(filePh,'%s\n',w_notes{:});
            fclose(filePh);
        end  
    end
end
end