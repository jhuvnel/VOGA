function MakeNotes(Raw_Path)
%% Notes about the Notes
% fn = {'RotaryChair-Sine','RotaryChair-VelStep','RotaryChair-SumSine',...
%         'aHIT-Sine','aHIT-Impulses',...
%         'eeVOR-Sine','eeVOR-PulseTrain','eeVOR-MultiVector','eeVOR-Autoscan','eeVOR-Activation'};
%     descrip = {'Condition-Axis-Freq-Speed','Condition-Axis-Speed(n=neg velocity)','Condition-Axis-Freq1-Freq2-Freq3-Speed',...
%         'Condition-Axis-Freq-Speed','Condition-Axis-Speed',...
%         'Axis-Freq-Speed','PFM/PAM-Axis-pps-uA','DOMAxis','Canal/Electrode-pps-us-uA','Light/Dark#'};
%% Find Files to Make Notes
rel_dir = [dir([Raw_Path,filesep,'*.txt']);dir([Raw_Path,filesep,'*.dat']);dir([Raw_Path,filesep,'*.mat'])]; %update if more file extensions are added
if isempty(rel_dir)
    disp(['No folders or files have been detected in ',Raw_Path])
    return;
end
file_names = extractfield(rel_dir,'name');
file_date = extractfield(rel_dir,'date');
Notes_ind = contains(file_names,'-Notes.txt');
VOG_ind = contains(file_names,{'SESSION','.dat','Lateral.txt','LARP.txt','RALP.txt','export.mat'})&~Notes_ind&~contains(file_names,'Raw');
VOG_ind_num = find(VOG_ind);
if all(~VOG_ind)
    disp(['No LDVOG, NKI, GNO, or ESC files have been detected in ',Raw_Path])
    return;
end
has_notes = contains(strrep(strrep(strrep(file_names(VOG_ind),'.txt',''),'.dat',''),'.mat',''),strrep(file_names(Notes_ind),'-Notes.txt',''));
VOG_files = file_names(VOG_ind_num(~has_notes));
VOG_files_date = file_date(VOG_ind_num(~has_notes));
%Set some parameters that will likely stay the same but can be edited
path_parts = strsplit(strrep(strrep(Raw_Path,'_',''),' ',''),filesep);
if any(contains(path_parts,'MVI')&contains(path_parts,'R')) %subject in expected formatting
    sub = path_parts{contains(path_parts,'MVI')&contains(path_parts,'R')};
    MVI_num = str2double(sub((-3:1:-1)+strfind(sub,'R')));
    if ismember(MVI_num,[1,2,3,4,7,9]) %EXPAND as needed
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
%% Partition by File Type
if any(contains(VOG_files,{'Lateral.txt','LARP.txt','RALP.txt'})) %GNO
    common_notes = inputdlg({'Subject:','Ear:','Visit:'},'Set VOG File Parameters',[1,40],{sub,ear,vis}); 
    sub = common_notes{1};
    ear = common_notes{2};
    vis = common_notes{3};
    gog = 'GNO';
    ang = '0';
    GNO_files = VOG_files(contains(VOG_files,{'Lateral.txt','LARP.txt','RALP.txt'}));
    for i = 1:length(GNO_files) 
        fname = GNO_files{i};
        dashes = strfind(fname,'_');
        date = fname;
        date(dashes(5)) = '-';
        date = strrep(date(dashes(2)+1:dashes(8)-1),'_','');
        w_notes = {['Subject ',sub];['Ear ',ear];['Visit ',vis];['Date ',date];['Goggle ',gog];['Angle ',ang];'Experiment '};
        notes_check = inputdlg([fname,newline,newline,'Check Notes: ',newline,'Ex: aHIT-Impulse-NoStim-LHRH-150dps-pseudorandom'],'Set VOG File Parameters',[length(w_notes),70],{strjoin(w_notes,'\n')});
        if ~isempty(notes_check)
            w_notes = cellstr(notes_check{1,1});
            filePh = fopen([Raw_Path,filesep,fname(1:end-4),'-Notes.txt'],'w');
            fprintf(filePh,'%s\n',w_notes{:});
            fclose(filePh);
        end
    end 
elseif any(contains(VOG_files,'export.mat')) %ESC
    common_notes = inputdlg({'Subject:','Ear:','Visit:'},'Set VOG File Parameters',[1,40],{sub,ear,vis}); 
    sub = common_notes{1};
    ear = common_notes{2};
    vis = common_notes{3};
    gog = 'ESC';
    ang = '0';
    ESC_files = VOG_files(contains(VOG_files,'export.mat'));
    for i = 1:length(ESC_files) 
        fname = ESC_files{i};
        dashes = strfind(fname,'-');
        date = strrep(strrep(strrep(fname(dashes(1)-4:dashes(1)+14),'-',''),'.',''),'_','-');
        w_notes = {['Subject ',sub];['Ear ',ear];['Visit ',vis];['Date ',date];['Goggle ',gog];['Angle ',ang];'Experiment '};
        notes_check = inputdlg([fname,newline,newline,'Check Notes: ',newline,'Ex: aHIT-Impulse-NoStim-LHRH-150dps-pseudorandom'],'Set VOG File Parameters',[length(w_notes),70],{strjoin(w_notes,'\n')});
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
    fig.Units = 'normalized';
    fig.Position = [0 0 1 1];
    fig.Units = 'inches';
    screen_size = fig.Position;
    fig.Position = screen_size - [0 0 6 0];
    annotation('textbox',[0 .9 1 .1],'String',strrep(Raw_Path,'_',' '),'FontSize',14,...
    'HorizontalAlignment','center','EdgeColor','none');
    %% Check each VOG file
    for i = 1:length(VOG_files)
        fname = VOG_files{i};
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
            gog = 'NKI1';
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
        w_notes = {['Subject ',sub];['Ear ',ear];['Visit ',vis];['Date ',date];['Goggle ',gog];['Angle ',ang]};
        if length(Time_Eye) > 1000000
            sub_i = floor(linspace(1,length(Time_Eye),1000000));
        else
            sub_i = 1:length(Time_Eye);
        end
        plot(Time_Eye(sub_i),GyroX(sub_i),'k:',Time_Eye(sub_i),GyroY(sub_i),'k--',Time_Eye(sub_i),GyroZ(sub_i),'k-',Time_Eye(sub_i),100*Stim(sub_i),'b')
        title(strrep(fname,'_',' '))
        xlabel('Time (s)')
        ylabel('Velocity (dps)')
        legend('GyroX','GyroY','GyroZ','Trigger')
        %Action options
        fig2 = figure(2);
        fig2.Units = 'inches';
        fig2.Position = [screen_size(3)-6,0,5,screen_size(4)-5.50]; 
        annotation('textbox',[0 0 1 1],'String',w_notes,...
            'FontSize',11,'HorizontalAlignment','left','EdgeColor','none');
        opts = {'Save','Edit Notes','Fix Trigger','Skip'};
        [ind,tf] = nmlistdlg('PromptString','Select an action:',...
               'SelectionMode','single',...
               'ListSize',[150 150],...
               'ListString',opts,...
               'Position',[screen_size(3)-6,screen_size(4)-3.75,3,3.75]); 
        if ~tf
            return;
        end
        while ~strcmp(opts{ind},'Save')&&~strcmp(opts{ind},'Skip')
            if strcmp(opts{ind},'Edit Notes')
                notes_check = inputdlg(['Check Notes: ',newline,'Ex: RotaryChair-Sine-NoStim-LHRH-0.05Hz-100dps'],'Set VOG File Parameters',[length(w_notes),70],{strjoin(w_notes,'\n')}); 
                if ~isempty(notes_check)
                    w_notes = cellstr(notes_check{1,1});
                    clf(fig2)
                    annotation('textbox',[0 0 1 1],'String',strrep(w_notes,'_','-'),...
                        'FontSize',11,'HorizontalAlignment','left','EdgeColor','none');
                else
                    return;
                end
            elseif strcmp(opts{ind},'Fix Trigger')
                updateRawVOGTrigger(Raw_Path,fname);
            end
            [ind,tf] = nmlistdlg('PromptString','Select an action:',...
               'SelectionMode','single',...
               'ListSize',[150 150],...
               'ListString',opts,...
               'Position',[screen_size(3)-6,screen_size(4)-3.75,3,3.75]);
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
        close(fig2)    
    end
end
end