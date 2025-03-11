%Takes in the Raw Path to LDHP/LDPC Fitting Software log files and makes the appropriate Notes files
%This ignores VOG files that already have notes files.
function logtoNotes(Raw_Path,ImplantSide)
%Update this line to include more extensions as needed
file_names = extractfield([dir([Raw_Path,filesep,'*.txt']);dir([Raw_Path,filesep,'*.dat']);dir([Raw_Path,filesep,'*.smr'])],'name');
file_dates = extractfield([dir([Raw_Path,filesep,'*.txt']);dir([Raw_Path,filesep,'*.dat']);dir([Raw_Path,filesep,'*.smr'])],'date');
if isempty(file_names)
    file_names = '';
end
Notes_ind = contains(file_names,'-Notes.txt');
Log_ind = contains(file_names,{'LDHP','LDPC','TestingLog','StimLogFile'});
VOG_ind = find(contains(file_names,{'SESSION','.dat','.smr'})&~Notes_ind&~Log_ind);
has_notes = contains(strrep(strrep(strrep(file_names(VOG_ind),'.txt',''),'.dat',''),'.smr',''),strrep(file_names(Notes_ind),'-Notes.txt',''));
VOG_files = file_names(VOG_ind(~has_notes));
possible_log = file_names(Log_ind);
if isempty(VOG_files)||isempty(possible_log) %No files or log files
    return;
end
VOG_files_date = file_dates(VOG_ind(~has_notes));
%% Find all log/autoscan/VOG files
%Start with log files
[indx_l,tf1] = nmlistdlg('PromptString','Select log files:','ListSize',[300 300],'ListString',possible_log);
if ~tf1
    return;
end
logfiles = possible_log(indx_l);
%% Initialize Figure
fig = figure(1);
set(fig,'Color','w','Units','normalized','Position',[0,0,1,1]);
clf; %in case there are leftover anotations
fig.Units = 'inches';
screen_size = fig.Position;
fig.Position = screen_size - [0 0 6 0];
ax = subplot(1,1,1);
xlabel('Time (s)')
ylabel('Stim Signals')
ax.Position = [0.05 0.1 0.6 0.83];
plot_notes = annotation('textbox',[0.7 0.1 0.25 0.83],'String','',...
    'FontSize',11,'HorizontalAlignment','left','EdgeColor','none');
%% Set some parameters that will likely stay the same but can be edited
path_parts = strsplit(strrep(strrep(Raw_Path,'_',''),' ',''),filesep);
sub = '';
ear = '';
vis = '';
if any(contains(path_parts,'MVI')&contains(path_parts,'R')) %subject in expected formatting
    sub = path_parts{contains(path_parts,'MVI')&contains(path_parts,'R')};
    MVI_num = str2double(sub((-3:1:-1)+strfind(sub,'R')));
    if ismember(MVI_num,ImplantSide{1})
        ear = 'L';
    elseif ismember(MVI_num,ImplantSide{2})
        ear = 'R';
    else
        open('MakeNotes.m')
        error('Remember to update the implantation side for this MVI Subject in MakeNotes.m')
    end
elseif any(contains(file_names,'ExperimentInfo'))
    experimentInfo_name = file_names{contains(file_names,'ExperimentInfo')};
    fid = fopen(fullfile(Raw_Path,experimentInfo_name),'r');
    experimentInfoText = textscan(fid,'%s','delimiter','\n');
    experimentInfoText = experimentInfoText{1};
    experimentInfo = split(experimentInfoText,char(9));

    sub = experimentInfo{contains(experimentInfo(:,1),'subject','IgnoreCase',true),2};
    ear = experimentInfo{contains(experimentInfo(:,1),'ear','IgnoreCase',true),2};
    vis = 'NA';

end
if any(contains(path_parts,'Visit'))
    vis = path_parts{contains(path_parts,'Visit')};
end
for f = 1:length(logfiles)
    %% Load and parse log file
    logFile = logfiles{f};
    fullpath = [Raw_Path,filesep,logfiles{f}];

    % Make sure text in log files for Ross/coil data is properly spaced]
    fid = fopen(fullpath);
    logText = textscan(fid,'%s','delimiter','\n');
    logText = logText{1};
    fclose(fid);

    if contains(logFile,'StimLogFile')
        for iLine = 1:length(logText)
            str = extractAfter(logText{iLine},'  ');
            stimFunction = str(isstrprop(str,'alpha'));
            switch stimFunction
                case 'keithleystimdctrapezoid'
                    if isscalar(split(extractAfter(logText{iLine},'  '),char(9)))
                        pat = extractAfter(logText{iLine},'keithley_stim_dc_trapezoid');
                        new_pat = strrep(pat,' ',char(9));
                        logText{iLine} = strrep(logText{iLine},pat,new_pat);
                    end
                case 'keithleystimdcsinusoid'
                    if isscalar(split(extractAfter(logText{iLine},'  '),char(9)))
                        pat = extractAfter(logText{iLine},'keithley_stim_dc_sinusoid');
                        new_pat = strrep(pat,' ',char(9));
                        logText{iLine} = strrep(logText{iLine},pat,new_pat);
                    end
            end
            logText{iLine} = strrep(logText{iLine},':  ',char(9));
        end

        logText = erase(logText,'"');

        fid = fopen(fullpath,'wt');
        fprintf(fid,'%s\n',logText{:});
        fclose(fid);
    end

    %Go through it once just to get dimensions
    fid = fopen(fullpath);
    tline = fgetl(fid);
    k = 0;
    num_cols=1;
    while ischar(tline)
        k = k+1;
        num_cols = max([num_cols,length(split(tline,char(9)))]);
        tline = fgetl(fid);
    end
    fclose(fid);
    %Now get the data
    data = cell(k,num_cols+1); %a little wider than needed for formatting
    fid = fopen(fullpath);
    tline = fgetl(fid);
    k = 0;
    while ischar(tline)
        k = k+1;
        data(k,1:length(split(tline,char(9)))) = split(tline,char(9));
        tline = fgetl(fid);
    end
    fclose(fid);
    %Now make sure it's all in the right format
    %Everything in the first column should be a timestamp (if one is
    %missing, it will be propaged)
    %There should be no empty arrays in the second column
    matchStr = regexp(data(:,1),'\w*:\w*:\w*.\w*','match');
    %Some of these are empty and some don't contain strings--either way,
    %scoot them over to the right one. Assumes first entry is fine
    col1ind = find(cellfun(@isempty,matchStr));
    for i = 1:length(col1ind)
        data(col1ind(i),2:num_cols+1) = data(col1ind(i),1:num_cols);
        data(col1ind(i),1) = data(col1ind(i)-1,1);
    end
    %Now keep scooting rows to the left until the second column has
    %something in it
    col2ind = find(cellfun(@isempty,data(:,2)));
    while ~isempty(col2ind)
        rm_row = false(size(data,1),1);
        for i = 1:length(col2ind)
            if all(cellfun(@isempty,data(col2ind(i),2:end)))
                rm_row(col2ind(i),1) = 1;
            else
                data(col2ind(i),2:num_cols+1) = [data(col2ind(i),3:num_cols+1),{[]}];
            end
        end
        data(rm_row,:) = [];
        col2ind = find(cellfun(@isempty,data(:,2)));
    end
    %Lastly remove any extra columns
    cellempty = all(cellfun(@isempty,data),1);
    data = data(:,~cellempty);
    %% Find relevant text files and make Notes files
    %Now find which text files may correspond to this log File.
    %Txt files must already be in the same folder as the log Files.
    try
        logfile_times = datetime(data(:,1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    catch
        try
            logfile_times = datetime(join([repmat({logFile(1:10)},k,1),data(:,1)]),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
        catch % data collected in lasker system - coils
            data(:,1) = strrep(data(:,1),'/','-');
            logfile_times = datetime(data(:,1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
        end
    end
    logfile_times.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
    %Some log files span several days so find and account for that
    day_inds = find(diff(logfile_times)<0);
    for i = 1:length(day_inds)
        temp = zeros(length(logfile_times),1);
        temp(day_inds(i)+1:end) = 1;
        logfile_times = logfile_times + temp;
    end
    %% Check each VOG file
    for i = 1:length(VOG_files)
        fname = VOG_files{i};
        if contains(fname,'.txt') %LDVOG
            VOG_data = readtable([Raw_Path,filesep,fname]);
            %Make date and other labels
            fname1 = fname;
            if any(strfind(fname,'_'))
                underscore = strfind(fname,'_');
                fname1(underscore(1):end-4) = [];
            end
            parts = split(fname1,'-');
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
        elseif contains(fname,'.dat') %NKI/NL
            %Load file
            warning('off')
            VOG_data = readtable([Raw_Path,filesep,fname],'ReadVariableNames',true);
            warning('on')
            VOG_data.Properties.VariableNames{1} = 'EyeTime';
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
            GyroX = VOG_data.GyroY - median(VOG_data.GyroY); 
            GyroY = -(VOG_data.GyroX - median(VOG_data.GyroX)); 
            GyroZ = -(VOG_data.GyroZ - median(VOG_data.GyroZ)); 
        elseif contains(fname,'.smr') %Spike2 coil data
            % % load eye coil field gains
            % 
            % FieldGainFile = dir(fullfile(Raw_Path,'*Gains*')).name;
            % fieldgainname = fullfile(Raw_Path, FieldGainFile);
            % delimiter = '\t';
            % formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
            % 
            % fileID = fopen(fieldgainname,'r');
            % 
            % dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
            % 
            % fclose(fileID);
            % 
            % FieldGains = [dataArray{1:end-1}];
            % 
            % % have no offsets with a test coil.
            % coilzeros = [0 0 0 0 0 0 0 0 0 0 0 0];
            % % This will use the first data point as a reference position
            % ref = 0;
            % % Apply a -pi/2 YAW reorientation of the raw data
            % data_rot = 2;
            % % Lasker Coil System data acquired using ONLY a CED and Spike 2 software
            % DAQ_code = 3;
            % OutputFormat = 1;
            % Data_In = [];
            % 
            % [d] = voma__processeyemovements(Raw_Path,VOG_files{i},FieldGains,coilzeros,ref,data_rot,DAQ_code,OutputFormat,Data_In);

            date = char(datetime(extract(Raw_Path,digitsPattern(8)),'InputFormat','yyyyMMdd'));
            gog = 'coils';
            ang = '0';

            fid = fopen(fullfile(Raw_Path,VOG_files{i}),'r');
            % ChanList = SONChanList(fid);
            TChannel = SONChannelInfo(fid,1);

            GyroZ = SONGetChannel(fid,2); % lasker system only has one axis of rotation, NOTE: this should be updated if gyro data is available
            GyroX = zeros(length(GyroZ),1);
            GyroY = zeros(length(GyroZ),1);

            elecStimTrig = SONGetChannel(fid,22); % keithley phase marker signal is on channel 22
            elecStimFlag = SONGetChannel(fid,30); % stim text flags are on channel 30
            
            fclose(fid);
          
            %Load items for plotting
            Time_Eye = 1:length(GyroZ);
            Fs = TChannel.idealRate;
            Stim = zeros(length(Time_Eye),1);
            Stim(int32(elecStimTrig*Fs)) = 1;

            flagText = [char(elecStimFlag.text)]';
            pat = digitsPattern(2) + ':' + digitsPattern(2) + ':' + digitsPattern(2);
            flagTimes = string(date) + ' ' + extract(string(flagText),pat);
            VOG_times = datetime(flagTimes,'InputFormat','dd-MMM-yyyy hh:mm:ss');

        else %Ignore unknown file type
            break;
        end
        common_notes = {['Subject ',sub];['Ear ',ear];['Visit ',vis];['Date ',date];['Goggle ',gog];['Angle ',ang]};
        %Find where the time stamps overlap
        rel_inds = logfile_times >= VOG_times(1) & logfile_times <= VOG_times(end);
        rel_col = ~all(cellfun(@isempty,data(rel_inds,:)),1);
        rel_dat = data(rel_inds,rel_col);
        %% Check in case it's in another log file from another computer
        if any(contains(logText,{'teensy','iDC','keithley','ced','labdata'})) && all(~contains(logText,'MVI')) % lasker system experiment- coil data
            % read relevant data for experiment types
            expTypes = unique(rel_dat(:,2));
            expTypes(contains(expTypes,{'pause','delay'})) = [];

            sweeps = rel_dat(matches(rel_dat(:,2),expTypes),:);
            nSweeps = sum(matches(rel_dat(:,2),expTypes));

            experiments = cell(nSweeps,1);

            for iSweep = 1:nSweeps
                
                switch sweeps{iSweep,2}
                    case 'keithley_stim_dc_sinusoid'
                        stimE = sweeps{iSweep,3};
                        refE = sweeps{iSweep,4};
                        stimCan = experimentInfo{contains(experimentInfo(:,1),stimE),2};
                        refCan = experimentInfo{contains(experimentInfo(:,1),refE),2};
                        canal = [stimCan,'-', refCan];

                        nCycles = sweeps{iSweep,5};
                        freq = num2str(round(str2double(sweeps{iSweep,6}),1));
                        amp = num2str(str2double(sweeps{iSweep,7})); 
                        offset = num2str(str2double(sweeps{iSweep,8})); 

                        experiments(iSweep) = {strcat('Experiment iDC-eeVOR-Sine-',canal,'-',freq,'Hz-',amp,'uA','-',offset,'uA-offset','-',nCycles,'cycles')};
                    case 'keithley_stim_dc_trapezoid'
                        keithleyNum = sweeps{iSweep,3};
                        stimE = sweeps{iSweep,4};
                        refE = sweeps{iSweep,5};
                        stimCan = experimentInfo{contains(experimentInfo(:,1),stimE),2};
                        refCan = experimentInfo{contains(experimentInfo(:,1),refE),2};
                        canal = [stimCan,'-', refCan];

                        nCycles = sweeps{iSweep,6};
                        interTrainGap = sweeps{iSweep,7};
                        phase1Dur = sweeps{iSweep,8}; 
                        ipg = sweeps{iSweep,9};
                        phase2Dur = sweeps{iSweep,10}; 
                        if contains(sweeps{iSweep,11},'-')
                            phase1Amp = [num2str(abs(str2double(sweeps{iSweep,11}))),'uA-cathodicFirstPhase']; 
                        else
                            phase1Amp = [num2str(str2double(sweeps{iSweep,11})),'uA-anodicFirstPhase']; 
                        end
                        if contains(sweeps{iSweep,12},'-')
                            phase2Amp = [num2str(abs(str2double(sweeps{iSweep,12}))),'uA-cathodicSecondPhase']; 
                        else
                            phase2Amp = [num2str(str2double(sweeps{iSweep,12})),'uA-anodicSecondPhase']; 
                        end
                        rampTime = sweeps{iSweep,13}; 
                        offset = num2str(str2double(sweeps{iSweep,14})); 

                        experiments(iSweep) = {strcat('Experiment iDC-eeVOR-Trapezoid-','K',keithleyNum,'-',canal,'-', ...
                            phase1Dur,'us-',phase1Amp,'-',phase2Dur,'us-',phase2Amp,'-',interTrainGap,'us-stimGap','-',ipg,'us-ipg','-', ...
                            offset,'uA-offset','-',rampTime,'us-ramp','-',nCycles,'cycles')};
                end
            end
            notes = experiments;
            stimCan = experimentInfo(contains(experimentInfo(:,1),'EyeCoils'),:);
            w_notes = [common_notes;[stimCan{1}, ' ', stimCan{2}];notes];
            %Plot
            sub_i = unique(floor(linspace(1,length(Time_Eye),1000000)));
            %Update plots
            plot(Time_Eye(sub_i),100*Stim(sub_i),'Color','b')
            hold on
            plot(Time_Eye(sub_i),GyroX(sub_i),'Color',0.5*[1,1,1])
            plot(Time_Eye(sub_i),GyroY(sub_i),'Color',0.75*[1,1,1])
            plot(Time_Eye(sub_i),GyroZ(sub_i),'Color','k')
            hold off
            title([{Raw_Path};{fname}],'FontSize',12,'interpreter','none')
            legend('Trigger','GyroX','GyroY','GyroZ')
            set(plot_notes,'String',w_notes)
            filePh = fopen([Raw_Path,filesep,fname(1:end-4),'-Notes.txt'],'w');
            fprintf(filePh,'%s\n',w_notes{:});
            fclose(filePh);
            set(plot_notes,'String','')
        elseif ~isempty(rel_dat)&&(any(contains(rel_dat(:,2),'Electrode characterization.'))||sum(rel_col)>2)
            if any(contains(rel_dat(:,2),'Electrode characterization.')) %Autoscan
                e_i = find(contains(rel_dat(:,2),'Electrode characterization.'));
                if ~any(contains(rel_dat(1:e_i(1),2),'(pps)'))
                    rel_inds = find(contains(data(1:find(rel_inds,1,'first'),2),'(pps)'),1,'last'):find(rel_inds,1,'last');
                    rel_dat = data(rel_inds,rel_col);
                    e_i = find(contains(rel_dat(:,2),'Electrode characterization.'));
                end
                rel_exp = rel_dat(e_i,2);
                %Make cell vector of experiment names
                experiments = cell(length(rel_exp),1);
                for j = 1:length(experiments)
                    line = rel_exp{j};
                    curr = num2str(round(str2double(strrep(line(strfind(line,':')+1:end),' ','')),0));
                    electrode = strrep(line(strfind(line,' E')+1:strfind(line,' E')+3),' ','');
                    rate_line = rel_dat{find(contains(rel_dat(1:e_i(j),2),'(pps)'),1,'last'),2};
                    pps = strrep(rate_line(strfind(rate_line,':')+1:end),' ','');
                    phase_line = rel_dat{find(contains(rel_dat(1:e_i(j),2),[electrode,' Phase Duration']),1,'last'),2};
                    us = strrep(phase_line(strfind(phase_line,':')+1:end),' ','');
                    switch electrode %canal
                        case {'E3','E4','E5'}
                            can = 'P';
                        case {'E6','E7','E8'}
                            can = 'H';
                        case {'E9','E10','E11'}
                            can = 'A';
                    end
                    if contains(common_notes(2),'L')
                        canal = ['L',can,electrode];
                    else
                        canal = ['R',can,electrode];
                    end
                    experiments{j} = [canal,'-',pps,'pps-',us,'us-',curr,'uA']; %rounds to closest uA for file naming
                end
                notes = strcat({'Experiment eeVOR-Autoscan-'},experiments);
            else 
                rel_dat(cellfun(@isempty,rel_dat(:,3)),:) = [];
                %Adjust for misalignment in stim start and file start
                %Scrolls up in the log file until it finds a header
                stim_ind = ~isnan(str2double(rel_dat(:,2)));
                while stim_ind(1)
                    rel_inds(find(rel_inds,1,'first')-1) = 1;
                    rel_dat = data(rel_inds,~all(cellfun(@isempty,data(rel_inds,:))));
                    stim_ind = ~isnan(str2double(rel_dat(:,2)));
                end
                stim_change = diff(stim_ind);
                stim_starts = find(stim_change==1);
                stim_ends = find(stim_change==-1);
                if length(stim_ends) < length(stim_starts)
                    stim_ends = [find(stim_change==-1);find(stim_ind,1,'last')];
                end
                experiments = cell(length(stim_starts),1);
                for j = 1:length(stim_starts)
                    stim_tab = rel_dat(stim_starts(j):stim_ends(j),2:end);
                    col_labs = stim_tab(1,:);
                    col_labs(cellfun(@isempty,col_labs)) = [];
                    %Figure out what type of experiment it is
                    if isempty(col_labs)
                        disp([fname,': Experiment type was not detected.'])
                    elseif any(contains(col_labs,'Depth of Modulation'))&&any(contains(data(rel_inds,2),'VelocityStep-LHRH'))
                        experiments(j) = {[{'Experiment eeVOR-VelStep-LH-240dps'};{'Experiment eeVOR-VelStep-RH-240dps'}]};
                    elseif any(contains(col_labs,'Depth of Modulation'))
                        experiments(j) = {strcat('Experiment eeVOR-MultiVector-[',stim_tab(2:end,2),',',stim_tab(2:end,3),',',stim_tab(2:end,4),']')};
                    elseif any(contains(col_labs,'Frequency'))
                        ax_vec = str2double(stim_tab(2:end,2:4));
                        stim_axis = cell(size(ax_vec,1),1);
                        %Round to the nearest 5 for continuity
                        amp_mat = 5*round(sqrt(ax_vec(:,1).^2+ax_vec(:,2).^2+ax_vec(:,3).^2)/5,0);
                        amp_vec = strrep(cellstr(num2str(amp_mat)),' ','');
                        stim_axis(amp_mat==ax_vec(:,1)) = {'LARP'};
                        stim_axis(amp_mat==ax_vec(:,2)) = {'RALP'};
                        stim_axis(amp_mat==ax_vec(:,3)) = {'LHRH'};
                        stim_axis(ax_vec(:,1)==ax_vec(:,2)&abs(ax_vec(:,1))>0) = {'X'};
                        stim_axis(ax_vec(:,1)==-ax_vec(:,2)&abs(ax_vec(:,1))>0) = {'Y'};
                        combo_vec = cellfun(@isempty,stim_axis);
                        if any(combo_vec)
                            stim_axis(combo_vec) = strcat('[',join(cellfun(@num2str,num2cell(ax_vec(combo_vec,:)./repmat(amp_mat(combo_vec),1,3)),'UniformOutput',false),','),']');
                        end
                        experiments(j) = {strcat('Experiment eeVOR-Sine-',stim_axis,'-',stim_tab(2:end,5),'Hz-',amp_vec,'dps')};
                    elseif any(contains(col_labs,'BSR (pps)'))
                        %Figure out axis and if it's PFM or PAM
                        %Assume that axes are only LHRH, RALP, or LARP and no
                        %combinations
                        p_vec = str2double(stim_tab(2:end,2:13));
                        [~,inds] = max([p_vec(:,1)~=p_vec(:,3),p_vec(:,2)~=p_vec(:,4),...
                            p_vec(:,5)~=p_vec(:,7),p_vec(:,6)~=p_vec(:,8),...
                            p_vec(:,9)~=p_vec(:,11),p_vec(:,10)~=p_vec(:,12)],[],2);
                        possible_exps = {'PFM_LARP','PAM_LARP','PFM_RALP','PAM_RALP','PFM_LHRH','PAM_LHRH'};
                        experiment_type = possible_exps(inds)';
                        pulse_freq = cell(length(experiment_type),1);
                        pulse_freq(inds==1|inds==2) = stim_tab(find(inds==1|inds==2)+1,4);
                        pulse_freq(inds==3|inds==4) = stim_tab(find(inds==3|inds==4)+1,8);
                        pulse_freq(inds==5|inds==6) = stim_tab(find(inds==5|inds==6)+1,12);
                        pulse_amp = cell(length(experiment_type),1);
                        pulse_amp(inds==1|inds==2) = stim_tab(find(inds==1|inds==2)+1,5);
                        pulse_amp(inds==3|inds==4) = stim_tab(find(inds==3|inds==4)+1,9);
                        pulse_amp(inds==5|inds==6) = stim_tab(find(inds==5|inds==6)+1,13);
                        %Round pulse_amp numbers
                        pulse_amp_round = cellfun(@num2str,num2cell(round(str2double(pulse_amp),0)),'UniformOutput',false);
                        experiments(j) = {strcat('Experiment eeVOR-PulseTrain-',experiment_type,'-',pulse_freq,'pps-',pulse_amp_round,'uA')};
                    else
                        disp([fname,': Experiment type was not detected.'])
                    end
                end
                notes = vertcat(experiments{:});
            end
            w_notes = [common_notes;notes];
            %Plot
            sub_i = unique(floor(linspace(1,length(Time_Eye),1000000)));
            %Update plots
            plot(Time_Eye(sub_i),100*Stim(sub_i),'Color','b')
            hold on
            plot(Time_Eye(sub_i),GyroX(sub_i),'Color',0.5*[1,1,1])
            plot(Time_Eye(sub_i),GyroY(sub_i),'Color',0.75*[1,1,1])
            plot(Time_Eye(sub_i),GyroZ(sub_i),'Color','k')
            hold off
            title([{Raw_Path};{fname}],'FontSize',12,'interpreter','none')
            legend('GyroX','GyroY','GyroZ','Trigger')
            set(plot_notes,'String',w_notes)
            filePh = fopen([Raw_Path,filesep,fname(1:end-4),'-Notes.txt'],'w');
            fprintf(filePh,'%s\n',w_notes{:});
            fclose(filePh);
            set(plot_notes,'String','')
        end
    end
end
end