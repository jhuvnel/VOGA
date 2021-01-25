function makeVOMAfile(Seg_Path)
    if nargin < 1
        Seg_Path = cd;
    end
    mat_files = dir([Seg_Path,filesep,'*.mat']);
    mat_files = {mat_files.name}';
    [indx,tf] = listdlg('PromptString','Select the .mat files','ListString',mat_files);
    if tf
        files = mat_files(indx);
        for i = 1:length(files)        
            load(files{i},'Data')
            %Build the v struct
            v.Fs = Data.Fs;
            if contains(files{i},'RotaryChair') %use GyroZ
                v.Stim_Trace = Data.HeadVel_Z;
            else
                v.Stim_Trace = Data.Trigger; %use trigger
            end
            v.Stim_t = Data.Time_Stim;
            v.stim_ind = [];
            v.Eye_t = Data.Time_Eye;
            v.Data_LE_Pos_X = Data.LE_Position_X;
            v.Data_LE_Pos_Y = Data.LE_Position_Y;
            v.Data_LE_Pos_Z = Data.LE_Position_Z;
            v.Data_RE_Pos_X = Data.RE_Position_X;
            v.Data_RE_Pos_Y = Data.RE_Position_Y;
            v.Data_RE_Pos_Z = Data.RE_Position_Z;
            %Calculate angular velocity
            %First param says no initial coordinate transforms, second is a 
            %struct with the unfiltered position traces.
            Data_cal = angpos2angvel(1,v); 
            v.Data_LE_Vel_X = Data_cal.LE_Vel_X;
            v.Data_LE_Vel_Y = Data_cal.LE_Vel_Y;
            v.Data_LE_Vel_LARP = Data_cal.LE_Vel_LARP;
            v.Data_LE_Vel_RALP = Data_cal.LE_Vel_RALP;
            v.Data_LE_Vel_Z = Data_cal.LE_Vel_Z;
            v.Data_RE_Vel_X = Data_cal.RE_Vel_X;
            v.Data_RE_Vel_Y = Data_cal.RE_Vel_Y;
            v.Data_RE_Vel_LARP = Data_cal.RE_Vel_LARP;
            v.Data_RE_Vel_RALP = Data_cal.RE_Vel_RALP;
            v.Data_RE_Vel_Z = Data_cal.RE_Vel_Z;
            %build the p struct
            p.Stim_Info.Stim_Type = {Data.info.dataType};
            p.Stim_Info.ModCanal = {''};
            p.Stim_Info.Freq = {''};
            p.Stim_Info.Max_Vel = {''};
            p.Stim_Info.Cycles = {''};
            p.Stim_Info.Notes = {''};
            p.DAQ = Data.info.goggle_ver;
            p.DAQ_code = 5; %technically for LDVOG but needed
            v.Parameters = p;
            %Place in larger struct
            Data_QPR(i).name = files{i};
            Data_QPR(i).VOMA_data = v;
            Data_QPR(i).SoftwareVer = struct([]);
            Data_QPR(i).RawFileName = Data.info.rawfile;
        end
        out_name = inputdlg('Name the .voma file (no suffix)');
        save([Seg_Path,filesep,out_name{:},'.voma'],'Data_QPR')
    end
end