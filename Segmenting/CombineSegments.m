function Data = CombineSegments(Data1,Data2)
info_fields = {'info','DetectedTraces_HeadVel','CSVData','XMLData','AllData'};
num_fields = {'Time_Eye','Time_Stim','Trigger','raw_start_t','raw_end_t',...
    'LE_Position_X','LE_Position_Y','LE_Position_Z',...
    'RE_Position_X','RE_Position_Y','RE_Position_Z',...
    'LE_Vel_X','LE_Vel_Y','LE_Vel_Z','RE_Vel_X','RE_Vel_Y','RE_Vel_Z',...
    'HeadVel_Z','HeadVel_L','HeadVel_R','HeadVel_X','HeadVel_Y',...
    'HeadAccel_X','HeadAccel_Y','HeadAccel_Z'};
%Make the time vectors continuous with the first segment
Data2.Time_Eye = Data2.Time_Eye - Data2.Time_Eye(1) + Data1.Time_Eye(end) + mean(diff(Data2.Time_Eye));
Data2.Time_Stim = Data2.Time_Stim - Data2.Time_Stim(1) + Data1.Time_Stim(end) + mean(diff(Data2.Time_Stim));
%Figure out how many copies of info there are 
copynum = num2str(sum(contains(fieldnames(Data1),'info'))+1);
%Combine all of the items
Data = Data1;
for i = 1:length(num_fields)
    if isfield(Data1,num_fields{i})&&isfield(Data2,num_fields{i})
        Data.(num_fields{i}) = [reshape(Data1.(num_fields{i}),[],1);reshape(Data2.(num_fields{i}),[],1)];
    end
end
Data.rawfile = [Data1.rawfile;Data2.rawfile];
for i = 1:length(info_fields)
    if isfield(Data2,info_fields{i})
        Data.([info_fields{i},copynum]) = Data2.(info_fields{i}); %as many info copies as needed
    end
end
Data.(['Data',copynum]) = Data2;
end