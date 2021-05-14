function Data = CombineSegments(Data1,Data2)
    %Make the time vectors continuous with the first segment
    Time_Eye2 = Data2.Time_Eye - Data2.Time_Eye(1) + Data1.Time_Eye(end) + mean(diff(Data2.Time_Eye));
    Time_Stim2 = Data2.Time_Stim - Data2.Time_Stim(1) + Data1.Time_Stim(end) + mean(diff(Data2.Time_Stim));
    %Figure out how many copies of info there are 
    copynum = sum(contains(fieldnames(Data1),'info'))+1;
    %Combine all of the items
    Data = Data1;
    Data.(['info',num2str(copynum)]) = Data2.info; %as many info copies as needed
    Data.(['Data',num2str(copynum)]) = Data2;
    Data.Time_Eye = [reshape(Data1.Time_Eye,[],1);reshape(Time_Eye2,[],1)];
    Data.Time_Stim = [reshape(Data1.Time_Stim,[],1);reshape(Time_Stim2,[],1)];
    Data.raw_start_t = [Data1.raw_start_t;Data2.raw_start_t];
    Data.raw_end_t = [Data1.raw_end_t;Data2.raw_end_t];
    Data.Trigger = [Data1.Trigger;Data2.Trigger];
    Data.LE_Position_X = [Data1.LE_Position_X;Data2.LE_Position_X];
    Data.LE_Position_Y = [Data1.LE_Position_Y;Data2.LE_Position_Y];
    Data.LE_Position_Z = [Data1.LE_Position_Z;Data2.LE_Position_Z];
    Data.RE_Position_X = [Data1.RE_Position_X;Data2.RE_Position_X];
    Data.RE_Position_Y = [Data1.RE_Position_Y;Data2.RE_Position_Y];
    Data.RE_Position_Z = [Data1.RE_Position_Z;Data2.RE_Position_Z];
    Data.HeadVel_X = [reshape(Data1.HeadVel_X,[],1);reshape(Data2.HeadVel_X,[],1)];
    Data.HeadVel_Y = [reshape(Data1.HeadVel_Y,[],1);reshape(Data2.HeadVel_Y,[],1)];
    Data.HeadVel_Z = [reshape(Data1.HeadVel_Z,[],1);reshape(Data2.HeadVel_Z,[],1)];
    Data.HeadAccel_X = [reshape(Data1.HeadAccel_X,[],1);reshape(Data2.HeadAccel_X,[],1)];
    Data.HeadAccel_Y = [reshape(Data1.HeadAccel_Y,[],1);reshape(Data2.HeadAccel_Y,[],1)];
    Data.HeadAccel_Z = [reshape(Data1.HeadAccel_Z,[],1);reshape(Data2.HeadAccel_Z,[],1)];
    Data.rawfile = [Data1.rawfile;Data2.rawfile];
end