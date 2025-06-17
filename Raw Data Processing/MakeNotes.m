%% Make Notes
%
% This function expects that all of the VOG files are in the Raw_Path
% directory and creates files with -Notes.txt appeneded to the name of
% the VOG file with a description of the experiments containted in the VOG
% file. This script also contains another function LogtoNotes which 
% automatically creates the Notes files based on the log files created by 
% the MVI fitting software. This is especially useful for the eeVOR 
% experiments. The text file generally appears in the following format using 
% spaces and newlines as delimeters. Notes about the items appear in
% brackets. Sometimes it is necessary to add a newline even after the last
% line of the text file.
%
% % Notes File Format:
% 
% Subject [Subject#]
% Ear [Implanted Ear for MVI, U for others]
% Visit [Visit # for MVI, NA for others]
% Date [Date and time in YYYYMMDD-hhmmss format]
% Goggle [VOG system used]
% Angle [Theta (pitch) angle used for gyro reorientation within the goggles]
% Experiment [Experiment type #1]
% Experiment [Experiment type #2]
% ....
% Experiment [Experiment type #N]
%
% % Experiment Types
%
% Experiment types define the type of experiment in a certain order. They 
% can be also added as new experiments are introduced but they generally
% appear in the format of:
%
% Experiment Apparatus (RotaryChair, aHIT, Manual, or eeVOR)
% Type of Waveform (Sine, SumSine, VelStep, Impulses, Gaussian, Autoscan,
%   PulseTrain, MultiVector, or Activation)
% Condition During Stimulation (NoStim, MotionMod, ConstantRate, LightNoStim, Light#, Dark#)
% Axis of Motion or Electrode (LHRH, LARP, RALP, X, Y, LH, RH, LA, RP, RA,
%   LP, LRZ format of [0,0,1] for Multivector, or canal name with E# like
%   LHE6 for Autoscan)
% Frequency (Hz)
% Amplitude (degrees/s)
% Pulse Frequency (pulses/s)
% Phase duration (us/phase)
% Current Amplitude (uA)
% 
% Here are the ones I have encountered:
%
% RotaryChair-Sine-Condition-Axis-Freq-Amp
% RotaryChair-VelStep-Condition-Axis-Amp
% RotaryChair-SumSine-Condition-Axis-Freq1-...-FreqN-Amp
% aHIT-Sine-Condition-Axis-Freq-Amp
% aHIT-Impulses-Condition-Axis-Amp
% aHIT-Gaussian-Condition-Axis-Amp
% Manual-Impulses-Condition-Axis-Speed
% eeVOR-Sine-Axis-Freq-Speed
% eeVOR-VelStep-Axis-Amp
% eeVOR-PulseTrain-Axis-PulseFreq-CurrAmp
% eeVOR-MultiVector-DOMAxis
% eeVOR-Autoscan-Canal/Electrode-PulseFreq-PhaseDur-CurrAmp
% eeVOR-Activation-Light/Dark#

function MakeNotes(Raw_Path)
%% Input
if nargin < 1
    Raw_Path = cd;
end
%Get subject info for which ear to note for an implanted MVI user
VOGA_VerInfo = rows2vars(readtable([userpath,filesep,'VOGA_VerInfo.txt'],...
    'ReadVariableNames',false,'ReadRowNames',true));
ImplantSide = {[],[]};
if isfolder(VOGA_VerInfo.Path{:})
    sub_info = readtable([VOGA_VerInfo.Path{:},filesep,'MVI_Information.xlsx']);
    ImplantSide = {find(strcmp(sub_info.Ear,'L')),find(strcmp(sub_info.Ear,'R'))};
end
% Process log files first
LogtoNotes(Raw_Path,ImplantSide)
%% Find Files to Make Notes
%These are the keywords in the file names that indicate it's a VOG file
VOG_fname_pat = {'SESSION','Lateral.txt','LARP.txt','RALP.txt','.dat','.mat','ImuData','Custom Test'};
rel_dir = dir(Raw_Path);
rel_dir(extractfield(rel_dir,'isdir')) = []; %remove folders
file_names = extractfield(rel_dir,'name');
file_date = extractfield(rel_dir,'date');
if isempty(file_names)
    file_names = ''; file_date = [];
end
Notes_ind = contains(file_names,'-Notes.txt');
%Find the files that fit the VOG file pattern, and aren't notes or
%calibration files
VOG_ind = find(contains(file_names,VOG_fname_pat)&~Notes_ind&~contains(file_names,{'Raw','.cal'})); %Raw = LDVOG calibration file
has_notes = contains(file_names(VOG_ind),strrep(file_names(Notes_ind),'-Notes.txt',''));
VOG_files = file_names(VOG_ind(~has_notes));
VOG_files_date = file_date(VOG_ind(~has_notes));
if ~any(VOG_ind)
    disp(['No VOG files (LDVOG/NL/GNO/ESC) have been detected: ',Raw_Path])
    return;
elseif isempty(VOG_files)
    disp(['All VOG Files have Notes files: ',Raw_Path]);
    return;
end
%% Set some parameters that will likely stay the same but can be edited
path_parts = strsplit(strrep(strrep(Raw_Path,'_',''),' ',''),filesep);
sub = 'Unknown'; ear = 'U'; vis = 'NA'; date = ''; gog = 'NA'; ang = '0';
if any(contains(path_parts,'MVI')&contains(path_parts,'R')) %subject in expected formatting
    sub = path_parts{contains(path_parts,'MVI')&contains(path_parts,'R')};
    MVI_num = str2double(sub((-3:1:-1)+strfind(sub,'R')));
    if ismember(MVI_num,ImplantSide{1})
        ear = 'L';
    elseif ismember(MVI_num,ImplantSide{2})
        ear = 'R';
    end
end
if any(contains(path_parts,'Visit'))
    vis = path_parts{contains(path_parts,'Visit')};
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
        type = 'Manual'; cond = ''; exp_info = '(No File)'; gog = 'GNO'; xml_file = [];
        %Try to plot the accepted head traces (if none, plot the whole time trace)
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
        all_xml = extractfield(dir([Raw_Path,filesep,'*.xml']),'name');
        if isfile([Raw_Path,filesep,fname(1:end-4),'.xml'])
            xml_file = [fname(1:end-4),'.xml'];
        elseif ~isempty(all_xml) %Look for a file with a small time delay
            xml_time = datetime(extract([all_xml;fname],digitsPattern(4)+'_'+digitsPattern(2)+'_'+digitsPattern(2)+'_'+digitsPattern(2)+'_'+digitsPattern(2)+'_'+digitsPattern(2)),'Format','yyyy_MM_dd_HH_mm_ss');
            if any(abs(seconds(xml_time(1:end-1)-xml_time(end)))<120)
                [~,ind] = min(abs(seconds(xml_time(1:end-1)-xml_time(end))));
                xml_file = all_xml{ind};
            end
        end
        if ~isempty(xml_file) %XML file with experiment notes exists!
            fdata = cellstr(readlines([Raw_Path,filesep,xml_file]));
            exp_info = strrep(extractXMLdataline(fdata{contains(fdata,'<Remarks>')}),' ','');
            if contains(lower(exp_info),'ahit') %Experiment
                type = 'aHIT';
            elseif contains(lower(exp_info),'chair')
                type = 'RotaryChair';
            end
            %Set condition
            if contains(lower(exp_info),{'off','nostim','preop','preact','pre-op'})||contains(Raw_Path,'Visit 0')
                cond = 'NoStim';
            elseif contains(lower(exp_info),{'constant','baseline'})
                cond = 'ConstantRate';
            elseif contains(lower(exp_info),{'motionmod','mod','accel'})
                cond = 'MotionMod';
            end
            %Find which goggle set #
            gog_line = fdata{contains(fdata,'GogglesSN')};
            gog = ['GNO',num2str(find(ismember(GNO_SerialNums,gog_line(ismember(gog_line,'0123456789')))))]; 
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
    %Defaults - all the metadata should be in the file name
    ang = '0'; canal = ''; type = 'Manual'; cond = '';
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
        date = char(VOG_time,'yyyyMMdd-HHmmss');
        if contains(exp_info,{'Lateral','LHRH'})
            canal = 'LHRH';
        elseif contains(exp_info,'LARP')
            canal = 'LARP';
        elseif contains(exp_info,'RALP')
            canal = 'RALP';
        end
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
    clf; %in case there are leftover anotations
    %Figure out how big the screen is and then leave 6 inches of space on
    %the right hand side
    set(fig,'Color','w','Units','normalized','Position',[0,0,1,1]);
    fig.Units = 'inches';
    screen_size = fig.Position;
    fig.Position = screen_size - [0 0 6 0];
    ax = subplot(1,1,1);
    xlabel('Time (s)')
    ylabel('Velocity (dps)')
    ax.Position = [0.05 0.1 0.6 0.83]; %Make space for the notes annotations
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
            date = char(VOG_times(1),'yyyyMMdd-HHmmss');
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
        elseif contains(fname,'.dat') || contains(fname,'Custom Test') %NKI
            %Load file
            warning('off')
            % sometimes has an error where it doesn't recognize the 3rd row
            % as column labels
            VOG_data = readtable([Raw_Path,filesep,fname],'ReadVariableNames',true);
            % if contains(fname,'.dat')
            %     VOG_data = readtable([Raw_Path,filesep,fname],'ReadVariableNames',true,'Range',[2 1]);
            % elseif contains(fname,'.csv')
            %     VOG_data = readtable([Raw_Path,filesep,fname],'ReadVariableNames',true,'Range',[3 1]);
            % end
            warning('on')
            VOG_data.Properties.VariableNames{1} = 'EyeTime';
            % VOG_data.Properties.VariableNames{19} = 'EyeTime';
            % VOG_data.Properties.VariableNames{20} = 'EyeTime';
            % VOG_data.Properties.VariableNames{21} = 'EyeTime';
            %Make date and other labels
            VOG_time = datetime(VOG_files_date{i}); %use file creation/saving time since NKI doesn't output time stamps
            VOG_time.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
            VOG_times = [VOG_time-seconds(VOG_data{end,1}) VOG_time];
            date = char(VOG_times(1),'yyyyMMdd-HHmmss');
            if VOG_times(1) > datetime(2024,01,31)
                gog = 'NL3'; %Changed from NL2 to NL3 before testing on 2024/01/31 
            elseif VOG_times(1) > datetime(2022,07,28)
                gog = 'NL2'; %Changed from NKI1 to NL2 on 2022/07/28
            else
                gog = 'NL1';
            end
            ang = '0';
            %Load items for plotting
            Time_Eye = VOG_data.EyeTime;        
            Stim = zeros(length(Time_Eye),1); 
            Stim(VOG_data.EventCode ~= 0) = 1;
            GyroX = medfilt1(VOG_data.GyroY - median(VOG_data.GyroY),3); 
            GyroY = medfilt1(-(VOG_data.GyroX - median(VOG_data.GyroX)),3); 
            GyroZ = medfilt1(-(VOG_data.GyroZ - median(VOG_data.GyroZ)),3); 
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
                'Experiment RotaryChair-Sine-NoStim-LHRH-0.5Hz-35dps';...
                'Experiment RotaryChair-Sine-NoStim-LHRH-1Hz-70dps'};
        elseif contains(Raw_Path,'aHIT') %All the normal experiments
            w_notes = {['Subject ',sub];['Ear ',ear];['Visit ',vis];['Date ',date];['Goggle ',gog];['Angle ',ang];...
                'Experiment aHIT-Sine-NoStim-LHRH-0.5Hz-35dps';...
                'Experiment aHIT-Sine-NoStim-LHRH-1Hz-70dps'};
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
        legend('Trigger','GyroX','GyroY','GyroZ')
        set(notes,'String',w_notes)
        %Action options
        opts = {'Save','Edit Notes','Fix Trigger','Skip'};
        resp = ''; %Inialize for the loop
        while ~contains(resp,{'Save','Skip'}) %Keep running until the user selects "Save" or "Skip"
            if strcmp(resp,'Edit Notes')
                notes_check = inputdlg(['Check Notes: ',newline,...
                    'Ex: RotaryChair-Sine-NoStim-LHRH-0.05Hz-100dps'],...
                    'Set VOG File Parameters',[length(w_notes),70],{strjoin(w_notes,'\n')}); 
                if ~isempty(notes_check)
                    w_notes = cellstr(notes_check{1,1});
                    set(notes,'String',w_notes)
                else
                    return;
                end
            elseif strcmp(resp,'Fix Trigger')
                updateRawVOGTrigger(Raw_Path,fname);
            end
            [ind,tf] = nmlistdlg('PromptString','Select an action:',...
               'SelectionMode','single','ListSize',[100 70],'ListString',opts,...
               'Position',[screen_size(3)-6,screen_size(4)-3,2,2]); 
            if tf
                resp = opts{ind};
            else
                return;
            end
        end
        if strcmp(resp,'Save')
            %If you got here it's time to save
            filePh = fopen([Raw_Path,filesep,fname(1:end-4),'-Notes.txt'],'w');
            fprintf(filePh,'%s\n',w_notes{:});
            fclose(filePh);
        end  
    end
end
end