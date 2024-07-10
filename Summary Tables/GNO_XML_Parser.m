%% GNO_XML_Parser
% Makes the Data struct and tab table for one folder with xml_files
% Based off of GNO_XML_Creator_v3 in the vHIT analysis folder of the paper
function [tab,Data] = GNO_XML_Parser(Path)
% Most common missing fields
missing_maybe = {'VW_HIImpulse','VW_HIVideo','CalibrationUID','VideoUID','CalibrationGUID','VideoGUID'};
% Values for latency estimation. 
t = (0:174)*1000/245 - 140; %Time is the GNO time scale shown on the PDFs
tt = t(1):1:t(end); %Make it a 1ms step
thresh = 15; %dps threshold for the head/eye start
xml_files = extractfield(dir([Path,filesep,'*.xml']),'name');
if isempty(xml_files)
    disp(['No .xml files were found on the path: ',Path])
    tab = [];
    Data = [];
    return;
end
notes_files = extractfield(dir([Path,filesep,'*-Notes.txt']),'name');
xml_time = datetime(extract(xml_files,digitsPattern(4)+'_'+digitsPattern(2)+'_'+digitsPattern(2)+'_'+digitsPattern(2)+'_'+digitsPattern(2)+'_'+digitsPattern(2)),'Format','yyyy_MM_dd_HH_mm_ss');
notes_time = datetime(extract(notes_files,digitsPattern(4)+'_'+digitsPattern(2)+'_'+digitsPattern(2)+'_'+digitsPattern(2)+'_'+digitsPattern(2)+'_'+digitsPattern(2)),'Format','yyyy_MM_dd_HH_mm_ss');
%Make Data Struct
Data = cell(length(xml_files),1);
for i = 1:length(xml_files)
    notes_file = [];
    if any(seconds(notes_time-xml_time(i))<120)
        [~,ind] = min(abs(seconds(notes_time-xml_time(i))));
        notes_file = notes_files{ind};
    end
    if ~isempty(notes_file)
        n_temp = table2cell(readtable([Path,filesep,notes_file],'ReadVariableNames',false,'Format','%s %s'));
        notes = n_temp{end,end};
        if ~contains(notes,{'trash','RotaryChair'})
            data = readstruct([Path,filesep,xml_files{i}]);
            for j = 1:length(missing_maybe)
                if ~isfield(data,missing_maybe{j})
                    data.(missing_maybe{j}) = [];
                end
            end
            data.FileName = xml_files{i};
            data.Notes = notes;
            data.Subject = n_temp{contains(n_temp(:,1),'Subject'),2};
            data.Visit = strrep(n_temp{contains(n_temp(:,1),'Visit'),2},'Visit','');
            Data{i} = data;
        end
    end
end
Data = vertcat(Data{:});
%Make tab table
n = length(Data);
tab = cell2table(cell(2*n,2),'VariableNames',{'Subject','Visit'});
tab.Subject(:) = {Data.Subject,Data.Subject}';
tab.Visit(:) = {Data.Visit,Data.Visit}';
tab.Date = [Data.StartDateTime,Data.StartDateTime]';
tab.Experiment = strrep({Data.Notes,Data.Notes}','manual','Manual');
tab.Canal = [strrep(extract(cellstr([Data.TestType]'),"Lat"|"LA"|"LP"),'Lat','LH');...
    strrep(extract(cellstr([Data.TestType]'),"Lat"|"RA"|"RP"),'Lat','RH')];
tab.Type = extract(tab.Experiment,"Manual"|"aHIT");
tab.Condition = extract(tab.Experiment,"NoStim"|"ConstantRate"|"MotionMod");
tab.Cycles = [Data.NumAcceptedLeftImpulses,Data.NumAcceptedRightImpulses]';
cyc_var = {'HeadCyc','EyeCyc','Saccade','HeadVelCyc','GainCyc'};
for j = 1:length(cyc_var)
    tab.(cyc_var{j})(:) = {[]}; %Initialize as an empty array
end
for i = 1:n
    if ~isempty(Data(i).VW_HIImpulse)
        subdat = Data(i).VW_HIImpulse;
        %Extract relevant variables and put them into a struct
        HeadVel = cellfun(@str2num,cellstr([subdat.HeadVelocitySamples])','UniformOutput',false);
        HeadVel = vertcat(HeadVel{:});
        EyeVel = cellfun(@str2num,cellstr([subdat.EyeVelocitySamples])','UniformOutput',false);
        EyeVel = vertcat(EyeVel{:});
        temp.HeadCyc = HeadVel;
        temp.EyeCyc = EyeVel;
        temp.Saccade = [subdat.FirstSaccadeLatency;subdat.SecondSaccadeLatency;subdat.ThirdSaccadeLatency]';
        temp.Saccade(temp.Saccade==-999) = NaN;
        temp.HeadVelCyc = [subdat.PeakVelocity]'/10000;
        temp.GainCyc = [subdat.Gain]'/10000;
        for j = 1:length(cyc_var)
            tab.(cyc_var{j}){i} = temp.(cyc_var{j})([subdat.IsDirectionLeft]' == "true",:);
            tab.(cyc_var{j}){i+n} = temp.(cyc_var{j})([subdat.IsDirectionLeft]' ~= "true",:);
        end
    end
end
%Calculate latency for each cycle based on the difference in
%time between the head and eye velocities reaching a threshold
tab.LatCyc(:) = {[]};
has_cyc = find(tab.Cycles>0);
for i = 1:length(has_cyc)
    HeadVel = tab.HeadCyc{has_cyc(i)};
    EyeVel = tab.EyeCyc{has_cyc(i)};
    LatCyc = NaN(size(HeadVel,1),1);
    for j = 1:size(HeadVel,1)
        headvel = spline(t,HeadVel(j,:),tt);
        eyevel = spline(t,EyeVel(j,:),tt);
        %Remove the eye trace before the head starts moving
        eyevel(1:find(headvel(tt<50)<thresh,1,'last')) = []; 
        if any(eyevel>thresh)
            LatCyc(j)= find(eyevel>thresh,1,'first');
        end
    end  
    tab.LatCyc{has_cyc(i)} = LatCyc;
end
% Means
tab.Gain = cellfun(@(x) round(mean(x),2),tab.GainCyc);
tab.Gain_sd = cellfun(@(x) round(std(x),2),tab.GainCyc);
tab.HeadVel = cellfun(@(x) round(mean(x)),tab.HeadVelCyc);
tab.HeadVel_sd = cellfun(@(x) round(std(x)),tab.HeadVelCyc);
tab.Lat = cellfun(@(x) round(mean(x)),tab.LatCyc);
tab.Lat_sd = cellfun(@(x) round(std(x)),tab.LatCyc);
end