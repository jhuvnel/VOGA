function MakeNotes(Raw_Path,no_notes)
%Set some initial parameters to use in the batch
prompt = {['Set the initial parameters.',newline,newline,'Subject: '];...
    'Implant Ear (L/R): ';'Visit: ';'Date: ';'Goggle Version: ';'Goggle Angle: '};
sub = '';
ear = '';
visit = '';
date = '';
gog = '';
ang = '';
path_parts = strsplit(strrep(strrep(Raw_Path,'_',''),' ',''),filesep);
if any(contains(path_parts,'MVI')&contains(path_parts,'R')) %subject in expected formatting
    sub = path_parts{contains(path_parts,'MVI')&contains(path_parts,'R')};
    MVI_num = str2double(sub((-3:1:-1)+strfind(sub,'R')));
    if ismember(MVI_num,[1,2,3,4,7,9]) %EXPAND as needed
        ear = 'L';
    else
        ear = 'R';
    end
end
if any(contains(path_parts,'Visit'))
   visit = path_parts{contains(path_parts,'Visit')};
end
if strcmp(no_notes{1}(end-2:end),'txt') %LDVOG
    date = datestr(datetime(no_notes{1}(end-19:end-11),'InputFormat','yyyyMMMdd'),'yyyymmdd');
    gog = 'LDVOG2';
    ang = '-170';
elseif strcmp(no_notes{1}(end-2:end),'dat')||strcmp(no_notes{1}(end-2:end),'mat') %NKI
    fparts = strsplit(no_notes{1},'_');
    if length(fparts)>4
        if length(fparts{2}) ==1
            fparts{2} = ['0',fparts{2}];
        end
        if length(fparts{3}) ==1
            fparts{3} = ['0',fparts{3}];
        end
        date = [fparts{4},fparts{2},fparts{3}];
    end
    gog = 'NKI1';
    ang = '0';
end
common_notes = {sub,ear,visit,date,gog,ang};
common_notes = inputdlg(prompt,'Set VOG File Parameters',1,common_notes);        
for i = 1:length(no_notes)
    In_Path = [Raw_Path,filesep,no_notes{i}];
    %Adjust date automatically
    if strcmp(no_notes{i}(end-2:end),'txt') %LDVOG
        common_notes{4} = datestr(datetime(no_notes{i}(end-19:end-11),'InputFormat','yyyyMMMdd'),'yyyymmdd');
    elseif strcmp(no_notes{i}(end-2:end),'dat')||strcmp(no_notes{i}(end-2:end),'mat') %NKI
        fparts = strsplit(no_notes{i},'_');
        if length(fparts)>4
            if length(fparts{2}) ==1
                fparts{2} = ['0',fparts{2}];
            end
            if length(fparts{3}) ==1
                fparts{3} = ['0',fparts{3}];
            end
            common_notes{4} = [fparts{4},fparts{2},fparts{3}];
        end
    end    
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
    elseif strcmp(In_Path(end-2:end),'mat') %NKI, preprocessed
        load(In_Path,'Data')
        data = Data;
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
    end
    if length(Time_Eye) > 1000000
        sub_i = floor(linspace(1,length(Time_Eye),1000000));
    else
        sub_i = 1:length(Time_Eye);
    end
    plot(Time_Eye(sub_i),GyroX(sub_i),'k:',Time_Eye(sub_i),GyroY(sub_i),'k--',Time_Eye(sub_i),GyroZ(sub_i),'k-',Time_Eye(sub_i),100*Stim(sub_i),'b')
    xlabel('Time (s)')
    ylabel('Velocity (dps)')
    legend('GyroX','GyroY','GyroZ','Trigger')
    %Have the user select the experiment type
    fn = {'RotaryChair-Sine','RotaryChair-VelStep','RotaryChair-SumSine',...
        'aHIT-Sine','aHIT-Impulses',...
        'eeVOR-Sine','eeVOR-PulseTrain','eeVOR-MultiVector','eeVOR-Autoscan','eeVOR-Activation'};
    descrip = {'Condition-Axis-Freq-Speed','Condition-Axis-Speed(n=neg velocity)','Condition-Axis-Freq1-Freq2-Freq3-Speed',...
        'Condition-Axis-Freq-Speed','Condition-Axis-Speed',...
        'Axis-Freq-Speed','PFM/PAM-Axis-pps-uA','DOMAxis','Canal/Electrode-pps-us-uA','Light/Dark#'};
    [indx,tf] = nmlistdlg('PromptString','Select the experiment type(s):',...
        'SelectionMode','multiple','ListString',fn);
    if tf
        if contains(fn(indx),'eeVOR-Autoscan')
            descrip1 = 'Canal/Electrode-pps-us';
            descrip2 = 'Current(uA) Level1-Level2-...';
            temp = inputdlg('Number of groups/electrodes tested','',[1 35],{'9'});
            num_el = str2double(temp{:});
            experiments = repmat({''},2*num_el,1);                
            notes1 = [common_notes;fn(indx);experiments];
            prompt = [{[In_Path,newline,newline,'Subject: '];'Implant Ear (L/R): ';...
                'Visit: ';'Date: ';'Goggle Version: ';'Goggle Angle: ';'Experiment Type: '};repmat({['Exp (',descrip1,'): '];descrip2},num_el,1)];
            notes1 = inputdlgcol(prompt,'Set VOG File Parameters',[1 45],notes1,'on',2); 
            subnotes = cell(10,num_el); %initialize to the largest vector needed
            for j = 1:num_el
                str1 = notes1{6+2*j};
                str2 = notes1{7+2*j};
                levels = strsplit(str2,'-')';
                subnotes(1:length(levels),j) = strcat(str1,'-',levels);
            end
            subnotes = reshape(subnotes,[],1);
            subnotes(cellfun(@isempty,subnotes)) = [];
            notes = [notes1(1:7);subnotes];
        else  
            if length(indx) > 1 
                %Select number of segments
                temp = inputdlg(strcat({'Number of segments of '},fn(indx),':'),'',[1 35],repmat({'1'},length(indx),1));
                temp = cellfun(@str2num,temp);
                experiments = cell(sum(temp),1);
                exp_desc = cell(sum(temp),1);
                k = 1;
                for j = 1:length(temp)
                    experiments(k:k+temp(j)-1) = {[fn{indx(j)},'-']};
                    exp_desc(k:k+temp(j)-1) = {['Exp (',descrip{indx(j)},'): ']};
                    k = k+temp(j);
                end             
                notes = [common_notes;strjoin(fn(indx),' ');experiments];
                prompt = [{[In_Path,newline,newline,'Subject: '];'Implant Ear (L/R): ';...
                    'Visit: ';'Date: ';'Goggle Version: ';'Goggle Angle: ';'Experiment Type: '};exp_desc];
                notes = inputdlgcol(prompt,'Set VOG File Parameters',[1 45],notes,'on',2);      
            else
                %Select number of segments
                temp = inputdlg('Number of segments','',[1 35],{'1'});
                experiments = repmat({''},str2double(temp{:}),1);                
                notes = [common_notes;fn(indx);experiments];
                prompt = [{[In_Path,newline,newline,'Subject: '];'Implant Ear (L/R): ';...
                    'Visit: ';'Date: ';'Goggle Version: ';'Goggle Angle: ';'Experiment Type: '};repmat({['Exp (',descrip{indx},'): ']},length(experiments),1)];
                notes = inputdlgcol(prompt,'Set VOG File Parameters',[1 45],notes,'on',2);      
            end
        end
        common_notes = notes(1:6);
        w_notes = strcat('"',notes,'"');
        filePh = fopen([In_Path(1:end-4),'-Notes.txt'],'w');
        fprintf(filePh,'%s\n',w_notes{:});
        fclose(filePh); 
    end
end
end