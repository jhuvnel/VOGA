function [Data,info,Fs,te,ts,stim1,type,detec_head,filt1,traces_vel1] = MakeCycAvg_startProcess(Data,info,filt1,all_traces)
info.TriggerShift2 = 0; %Shifting done manually in this file
Fs = Data.Fs;
if contains(info.goggle_ver,'GNO') %No raw position, just velocity
    te = Data.Time_Eye - Data.Time_Eye(1);
    ts = Data.Time_Stim - Data.Time_Stim(1);
    if contains(info.dataType,{'LH','RH'})
        stim1 = Data.HeadVel_Z;
    elseif contains(info.dataType,{'LA','RP'})
        stim1 = Data.HeadVel_L; %was negative here until 01/10/22
    elseif contains(info.dataType,{'RA','LP'})
        stim1 = Data.HeadVel_R;
    else
        stim1 = Data.HeadVel_Z;
    end
elseif contains(info.goggle_ver,'ESC') %No raw position, just velocity
    te = Data.Time_Eye - Data.Time_Eye(1);
    ts = Data.Time_Stim - Data.Time_Stim(1);
    if contains(info.dataType,{'LH','RH'})
        stim1 = Data.HeadVel_Z;
    elseif contains(info.dataType,{'LA','RP'})
        stim1 = Data.HeadVel_L;
    elseif contains(info.dataType,{'RA','LP'})
        stim1 = Data.HeadVel_R;
    else
        stim1 = Data.HeadVel_Z;
    end
else
    if contains(info.dataType,'Activation')
        %preserve the time because this will be used to rejoin them
        te = Data.Time_Eye;
        ts = Data.Time_Stim;
    else
        te = Data.Time_Eye - Data.Time_Eye(1);
        ts = Data.Time_Stim - Data.Time_Stim(1);
    end
    %Assign Trigger
    if contains(info.goggle_ver,'Moogles') %MOOG room coil system
        stim1 = Data.Trigger;
    else %NKI or LDVOG Trigger = internal gyro or comes from the PCU for eeVOR
        if contains(info.dataType,'RotaryChair')
            if isfield(Data,'HeadVel_Z')
                stim1 = Data.HeadVel_Z;
            else
                stim1 = Data.HeadMPUVel_Z;
            end
        elseif contains(info.dataType,'aHIT')
            if contains(info.dataType,'LHRH')
                stim1 = Data.HeadVel_Z;
            elseif contains(info.dataType,'LARP')
                stim1 = (Data.HeadVel_X - Data.HeadVel_Y)/sqrt(2);
            elseif contains(info.dataType,'RALP')
                stim1 = (Data.HeadVel_X + Data.HeadVel_Y)/sqrt(2);
            else
                stim1 = Data.Trigger;
            end
        else
            stim1 = Data.Trigger;
        end
    end
    %Fix huge number of NaN values in torsion traces of NKI traces
    if contains(info.goggle_ver,'NKI')
        if sum(isnan(Data.LE_Position_X)) > 0.9*length(te) %less than 10% data integrity
            Data.LE_Position_X = zeros(length(te),1); %set to 0 so no torsion
        else
            Data.LE_Position_X = spline(te(~isnan(Data.LE_Position_X)),Data.LE_Position_X(~isnan(Data.LE_Position_X)),te);
        end
        if sum(isnan(Data.RE_Position_X)) > 0.9*length(te)
            Data.RE_Position_X = zeros(length(te),1);
        else
            Data.RE_Position_X = spline(te(~isnan(Data.RE_Position_X)),Data.RE_Position_X(~isnan(Data.RE_Position_X)),te);
        end
    end
end
% Assign Type and set detected traces
if contains(info.goggle_ver,{'GNO','ESC'}) %No raw pos traces
    type = 3;
    if isfield(Data,'DetectedTraces_HeadVel')
        detec_head = Data.DetectedTraces_HeadVel;
    else
        detec_head = [];
    end
elseif contains(info.dataType,{'Activation','Step'}) %No cycle averaging
    type = 2;
    detec_head = [];
else
    type = 1;
    detec_head = [];
end
% Set some defaults
% Cycle Align
[~,t_snip] = MakeCycAvg__alignCycles(info,Fs,ts,stim1,[]);
if type == 1
    filt1.vel.irlssmooth(end) = round(length(t_snip)*0.16); %heuristic
end
if contains(info.goggle_ver,'GNO')&&contains(info.dataType,{'LH','RH'})
    traces_vel1 = {'RZ'};
elseif contains(info.goggle_ver,'GNO')&&contains(info.dataType,{'LA','RP'})
    traces_vel1 = {'RLARP'};
elseif contains(info.goggle_ver,'GNO')&&contains(info.dataType,{'RA','LP'})
    traces_vel1 = {'RRALP'};
elseif contains(info.goggle_ver,'ESC')&&contains(info.dataType,{'LH','RH'})
    traces_vel1 = {'LZ'};
elseif contains(info.goggle_ver,'ESC')&&contains(info.dataType,{'LA','RP'})
    traces_vel1 = {'LLARP'};
elseif contains(info.goggle_ver,'ESC')&&contains(info.dataType,{'RA','LP'})
    traces_vel1 = {'LRALP'};
elseif contains(info.dataType,{'X','Y'})
    traces_vel1 = all_traces(1:6); %LRZ vel
else
    traces_vel1 = all_traces(5:10); %LRZ vel
end
end