%Takes in the Raw Path to Fitting Software log files and makes the appropriate Notes files
%This ignores VOG files that already have notes files.
function logtoNotes(Raw_Path)
    rel_dir = [dir([Raw_Path,filesep,'*.txt']);dir([Raw_Path,filesep,'*.dat'])]; %update if more file extensions are added
    file_names = extractfield(rel_dir,'name');
    file_date = extractfield(rel_dir,'date');
    Notes_ind = contains(file_names,'-Notes.txt');
    VOG_ind = contains(file_names,{'SESSION','.dat'})&~Notes_ind;
    VOG_ind_num = find(VOG_ind);
    if all(~VOG_ind)
        error(['No VOG files have been detected in',Raw_Path])
    end
    has_notes = ismember(strrep(strrep(file_names(VOG_ind),'.txt',''),'.dat',''),strrep(file_names(Notes_ind),'-Notes.txt',''));
    VOG_files = file_names(VOG_ind_num(~has_notes));
    VOG_files_date = file_date(VOG_ind_num(~has_notes));    
    if all(VOG_ind|Notes_ind)
        disp('No log files detected:')
        disp(Raw_Path)
        return;
    elseif isempty(VOG_files)
        disp('All VOG files already have Notes files')
        disp(Raw_Path)
        return;
    end    
    possible_log = file_names(~Notes_ind&~VOG_ind);
    %% Find all log/autoscan/VOG files    
    %Start with log files
    [indx_l,tf1] = nmlistdlg('PromptString','Select log files:','ListSize',[300 300],'ListString',possible_log);
    if tf1 == 1
        logfiles = possible_log(indx_l);
        for f = 1:length(logfiles)       
            %% Load and parse log file   
            logFile = logfiles{f};
            fullpath = [Raw_Path,filesep,logfiles{f}];
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
                logfile_times = datetime(join([repmat({logFile(1:10)},k,1),data(:,1)]),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            end
            logfile_times.Format = 'yyyy-MM-dd HH:mm:ss.SSS';    
            %Some log files span several days so find and account for that
            day_inds = find(diff(logfile_times)<0);
            for i = 1:length(day_inds)
                temp = zeros(length(logfile_times),1);
                temp(day_inds(i)+1:end) = 1;
                logfile_times = logfile_times + temp;
            end  
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
            if any(contains(VOG_files,'.txt')) %LDVOG
                common_notes = {sub,ear,vis,datestr(logfile_times(1),'yyyymmdd'),'LDVOG2','-170','200','50'};   
            else
                common_notes = {sub,ear,vis,datestr(logfile_times(1),'yyyymmdd'),'NKI1','0','200','50'};  
            end   
            prompt = {['Set the initial parameters.',newline,newline,'Subject: '];'Implant Ear (L/R): ';'Visit: ';'Date: ';'Goggle Version: ';'Goggle Angle: '};
            common_notes = inputdlg(prompt,'Set VOG File Parameters',1,common_notes);        
            %% For loop
            for i = 1:length(VOG_files)
                fname = VOG_files{i};
                VOG_data = readtable([Raw_Path,filesep,fname]);
                if contains(fname,'.txt') %LDVOG  
                    parts = split(fname,'-');
                    if isduration(VOG_data{2,end})
                        %Use the timestamps on the file itself
                        VOG_times = datetime(strrep(parts{2},'.txt',''),'InputFormat','yyyyMMMMdd')+[VOG_data{1,end},VOG_data{end,end}];
                    else %Use file creation time
                        VOG_times = [datetime(strrep([parts{2},' ',parts{3}],'.txt',''),'InputFormat','yyyyMMMMdd HHmmss'),datetime(VOG_files_date{i})]; 
                    end
                    VOG_times.Format = 'yyyy-MM-dd HH:mm:ss.SSS';  
                else %NKI
                    VOG_time = datetime(VOG_files_date{i}); %use file creation/saving time since NKI doesn't output time stamps
                    VOG_time.Format = 'yyyy-MM-dd HH:mm:ss.SSS'; 
                    VOG_times = [VOG_time-seconds(VOG_data{end,1}) VOG_time];
                end
                if VOG_times(1) >= logfile_times(1) && VOG_times(1) <= logfile_times(end)
                    %Find where the time stamps overlap
                    rel_inds = logfile_times >= VOG_times(1) & logfile_times <= VOG_times(end);
                    rel_dat = data(rel_inds,:);   
                    rel_dat = rel_dat(:,~all(cellfun(@isempty,rel_dat)));
                    if ~isempty(rel_dat) %in case it's in another log file from another computer
                        %Find the number of experiments done and add them
                        %to the experiment cell block by block. Experiment
                        %types can differ.
                        if any(contains(rel_dat(:,2),'Electrode characterization.')) %Autoscan
                            e_i = find(contains(rel_dat(:,2),'Electrode characterization.'));
                            rel_exp = rel_dat(e_i,2);
                            %Make cell vector of experiment names
                            experiments = cell(length(rel_exp),1);
                            for j = 1:length(experiments)
                                line = rel_exp{j};
                                curr = num2str(round(str2num(strrep(line(strfind(line,':')+1:end),' ','')),0));
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
                                if strcmp(common_notes{2},'L')
                                    canal = ['L',can,electrode];
                                else
                                    canal = ['R',can,electrode];
                                end
                                
                                experiments{j} = [canal,'-',pps,'pps-',us,'us-',curr,'uA']; %rounds to closest uA for file naming
                            end
                            notes = [common_notes;{'eeVOR-Autoscan'};experiments];
                        else
                            stim_ind = ~isnan(str2double(rel_dat(:,2)));
                            while stim_ind(1) %Misalignemnt in stim start and file start
                                rel_inds(find(rel_inds,1,'first')-1) = 1;
                                rel_dat = data(rel_inds,:);   
                                rel_dat = rel_dat(:,~all(cellfun(@isempty,rel_dat)));
                                stim_ind = ~isnan(str2double(rel_dat(:,2)));
                            end
                            stim_change = diff(stim_ind);
                            stim_starts = find(stim_change==1);
                            stim_ends = find(stim_change==-1);
                            if length(stim_ends) < length(stim_starts)
                                stim_ends = [find(stim_change==-1);find(stim_ind,1,'last')];
                            end
                            experiments = [];   
                            for j = 1:length(stim_starts)
                                stim_tab = rel_dat(stim_starts(j):stim_ends(j),2:end);
                                %Figure out what type of experiment it is            
                                if any(contains(stim_tab(1,:),'Depth of Modulation'))                
                                    exp_type = 'eeVOR-MultiVector';
                                    experiments = [experiments;strcat('[',stim_tab(2:end,2),',',stim_tab(2:end,3),',',stim_tab(2:end,4),']')];              
                                elseif any(contains(stim_tab(1,:),'Frequency'))
                                    exp_type = 'eeVOR-Sine';
                                    ax_vec = str2double(stim_tab(2:end,2:4));
                                    stim_axis = cell(size(ax_vec,1),1);
                                    %Round to nearest 10 for continuity
                                    amp_mat = round(sqrt(ax_vec(:,1).^2+ax_vec(:,2).^2+ax_vec(:,3).^2),-1);
                                    amp_vec = strrep(cellstr(num2str(amp_mat)),' ','');
                                    stim_axis(amp_mat==ax_vec(:,1)) = {'LARP'};
                                    stim_axis(amp_mat==ax_vec(:,2)) = {'RALP'};
                                    stim_axis(amp_mat==ax_vec(:,3)) = {'LHRH'};
                                    stim_axis(ax_vec(:,1)==ax_vec(:,2)&abs(ax_vec(:,1))>0) = {'X'};
                                    stim_axis(ax_vec(:,1)==-ax_vec(:,2)&abs(ax_vec(:,1))>0) = {'Y'};                            
                                    experiments = [experiments;strcat(stim_axis,'-',stim_tab(2:end,5),'Hz-',amp_vec,'dps')];                                          
                                elseif any(contains(stim_tab(1,:),'BSR (pps)'))
                                    exp_type = 'eeVOR-PulseTrain';
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
                                    experiments = [experiments;strcat(experiment_type,'-',pulse_freq,'pps-',pulse_amp_round,'uA')];    
                                else
                                    disp([fname,': Experiment type was not detected.'])
                                    exp_type = 'Unknown';
                                end
                            end
                            notes = [common_notes;{exp_type};experiments];
                            prompt = [{[fname,newline,newline,'Subject: '];'Implant Ear (L/R): ';...
                                'Visit: ';'Date: ';'Goggle Version: ';'Goggle Angle: ';'Experiment Type: '};repmat({'Exp: '},length(experiments),1)];
                            notes = inputdlgcol(prompt,'Confirm VOG File Parameters',[1 45],notes,'on',2);        
                            common_notes = notes(1:6);
                        end
                        %Now make the file
                        filePh = fopen([Raw_Path,filesep,fname(1:end-4),'-Notes.txt'],'w');
                        w_notes = strcat('"',notes,'"');
                        fprintf(filePh,'%s\n',w_notes{:});
                        fclose(filePh); 
                    end
                end       
            end
        end 
    end
%    %% Autoscan 
%     %Start with log files
%     [indx_a,tf2] = nmlistdlg('PromptString','Select autoscan files:','ListSize',[300 300],'ListString',possible_log);
%     if tf2 == 1
%         logfiles = possible_log(indx_a);
%         % Find the times of the autoscan files (multiple may correspond to
%         % one VOG file)
%         logfile_times = NaT(length(logfiles),1);
%         logfile_times.Format = 'yyyy-MM-dd HH:mm:ss.SSS'; 
%         for f = 1:length(logfiles)
%             logfile_times(f) = datetime(strrep(strrep(logfiles{f},'  DAY @ TIME  ',' '),'.txt',''),'InputFormat','yyyy-MM-dd h-mm-ss-a');    
%         end        
%         %% Find relevant text files and make Notes files
%         %Set some parameters that will likely stay the same but can be edited 
%         path_parts = strsplit(strrep(strrep(Raw_Path,'_',''),' ',''),filesep);
%         if any(contains(path_parts,'MVI')&contains(path_parts,'R')) %subject in expected formatting
%             sub = path_parts{contains(path_parts,'MVI')&contains(path_parts,'R')};
%             MVI_num = str2double(sub((-3:1:-1)+strfind(sub,'R')));
%             if ismember(MVI_num,[1,2,3,4,7,9]) %EXPAND as needed
%                 ear = 'L';
%             else
%                 ear = 'R';
%             end
%         else
%             sub = '';
%             ear = '';
%         end
%         if any(contains(path_parts,'Visit'))
%            vis = path_parts{contains(path_parts,'Visit')};
%         else
%             vis = '';
%         end
%         if any(contains(VOG_files,'.txt')) %LDVOG
%             common_notes = {sub,ear,vis,datestr(logfile_times(1),'yyyymmdd'),'LDVOG2','-170','200','50'};   
%         else
%             common_notes = {sub,ear,vis,datestr(logfile_times(1),'yyyymmdd'),'NKI1','0','200','50'};  
%         end
%         % For loop
%         for i = 1:length(VOG_files)
%             fname = VOG_files{i};
%             if contains(fname,'.txt') %LDVOG  
%                 parts = split(fname,'-');
%                 VOG_times = [datetime(strrep([parts{2},' ',parts{3}],'.txt',''),'InputFormat','yyyyMMMMdd HHmmss'),datetime(VOG_files_date{i})]; 
%                 VOG_times.Format = 'yyyy-MM-dd HH:mm:ss.SSS';  
%             else %NKI
%                 VOG_data = readtable([Raw_Path,filesep,fname]);
%                 VOG_time = datetime(VOG_files_date{i}); %use file creation/saving time since NKI doesn't output time stamps
%                 VOG_time.Format = 'yyyy-MM-dd HH:mm:ss.SSS'; 
%                 VOG_times = [VOG_time-seconds(VOG_data{end,1}) VOG_time];
%             end
%             in_range = VOG_times(1) <= logfile_times & VOG_times(end) >= logfile_times;
%             if any(in_range) %Found corresponding VOG file, this assumes one VOG file to one autoscan file
%                 %% Get notes on this file
%                 prompt = {[fname,newline,newline,'Subject: '];'Implant Ear (L/R): ';'Visit: ';'Date: ';'Goggle Version: ';'Goggle Angle: ';'Burst rate (pps): ';'Phase duration (us): '};
%                 common_notes = inputdlg(prompt,'Set VOG File Parameters',1,common_notes);  
%                 % Load and parse log file   
%                 data = readtable([Raw_Path,filesep,logfiles{in_range}],'PreserveVariableNames',0);
%                 exp_type = 'eeVOR-Autoscan';
%                 experiments = cell(size(data,1),1);
%                 %Make cell vector of experiment names
%                 for j = 1:length(experiments)
%                     switch data{j,1}{:} %canal
%                         case {'E3','E4','E5'}
%                             can = 'P';
%                         case {'E6','E7','E8'}
%                             can = 'H';
%                         case {'E9','E10','E11'}
%                             can = 'A';
%                     end
%                     if strcmp(common_notes{2},'L')
%                         canal = ['L',can,data{j,1}{:}];
%                     else
%                         canal = ['R',can,data{j,1}{:}];
%                     end
%                     experiments{j} = [canal,'-',common_notes{7},'pps-',common_notes{8},'us-',num2str(round(data{j,3},0)),'uA']; %rounds to closest uA for file naming
%                 end               
%                 notes = [common_notes(1:6);{exp_type};experiments];
%                 % Don't show because too many entries
%                 filePh = fopen([Raw_Path,filesep,fname(1:end-4),'-Notes.txt'],'w');
%                 w_notes = strcat('"',notes,'"');
%                 fprintf(filePh,'%s\n',w_notes{:});
%                 fclose(filePh);  
%             end
%        end
%    end
end