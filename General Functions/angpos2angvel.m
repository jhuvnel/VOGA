function [Data,Data_In] = angpos2angvel(Data_In)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Adapted from PJB's voma function voma__processeyemovements.m
% Typically, we care more about how the EYES are moving in a HEAD
% coordinate frame (i.e., [Eye_in_Head]. We will need to apply some coordinate system changes
% (also termed 'passive' rotations) to get there.
%% Extract and process the raw data
% LD VOG Goggles - MVI Trial
Fs = Data_In.Fs;
Data_In.LE_Position_X = reshape(Data_In.LE_Position_X,[],1);
Data_In.LE_Position_Y = reshape(Data_In.LE_Position_Y,[],1);
Data_In.LE_Position_Z = reshape(Data_In.LE_Position_Z,[],1);
Data_In.RE_Position_X = reshape(Data_In.RE_Position_X,[],1);
Data_In.RE_Position_Y = reshape(Data_In.RE_Position_Y,[],1);
Data_In.RE_Position_Z = reshape(Data_In.RE_Position_Z,[],1);
rawData_L = [Data_In.LE_Position_X,Data_In.LE_Position_Y,Data_In.LE_Position_Z];
rawData_R = [Data_In.RE_Position_X,Data_In.RE_Position_Y,Data_In.RE_Position_Z];

%Get LARP and RALP in position too
angpos_L = (rotZ3deg(-45)'*rawData_L')'; 
Data_In.LE_Position_LARP = angpos_L(:,1);
Data_In.LE_Position_RALP = angpos_L(:,2);
angpos_R = (rotZ3deg(-45)'*rawData_L')'; 
Data_In.RE_Position_LARP = angpos_R(:,1);
Data_In.RE_Position_RALP = angpos_R(:,2);
%% Compute Angular Velocity
for j=1:2  
    % Lab. Dev. VOG Goggles            
    if j == 1
        rawData = rawData_L;
    elseif j == 2
        rawData = rawData_R;
    end    
    psi = rawData(:,1);
    phi = rawData(:,2);
    theta = rawData(:,3);
    % Computing angular velocity from Fick angular position angles
    angvel_dps_b = [(gradient(psi)*Fs.*cosd(theta).*cosd(phi)) - (gradient(phi)*Fs.*sind(theta)) ...
        (gradient(psi)*Fs.*sind(theta).*cosd(phi)) + (gradient(phi)*Fs.*cosd(theta)) ...
        (gradient(theta)*Fs) - (gradient(psi)*Fs.*sind(phi))];            
    % We want to rotate our data from an [X,Y,Z] coordinate system,
    % into a [LARP,RALP,LHRH] coordinate system, where Z = LHRH. To
    % accomplish this, we will perform a PASSIVE (i.e., 'alias' or
    % 'coordinate system') rotation of -45 degrees. In order to realize
    % this rotation, I will generate a -45deg rotation matrix and
    % RIGHT multiple our data by the TRANSPOSE of this rotation
    % matrix. Note that rotation matrices are orthonormal, and
    % their inverses are equivalent to their transpose. -PJB
    angvel_dps_c = (rotZ3deg(-45)'*angvel_dps_b')';            
    switch j        
        case 1
            % Store Data
            Data.LE_Vel_X = angvel_dps_b(:,1);
            Data.LE_Vel_Y = angvel_dps_b(:,2);
            Data.LE_Vel_Z = angvel_dps_b(:,3);
            Data.LE_Vel_LARP = angvel_dps_c(:,1);
            Data.LE_Vel_RALP = angvel_dps_c(:,2);            
        case 2
            % Store Data
            Data.RE_Vel_X = angvel_dps_b(:,1);
            Data.RE_Vel_Y = angvel_dps_b(:,2);
            Data.RE_Vel_Z = angvel_dps_b(:,3);
            Data.RE_Vel_LARP = angvel_dps_c(:,1);
            Data.RE_Vel_RALP = angvel_dps_c(:,2);
    end
end
if exist('ElecStimTrig','var')
    Data.ElecStimTrig = ElecStimTrig;
end
Data.Fs = Fs;
% Saving both POSITION and RAW data traces
Data.LE_Pos_X = rawData_L(:,1);
Data.LE_Pos_Y = rawData_L(:,2);
Data.LE_Pos_Z = rawData_L(:,3);
Data.RE_Pos_X = rawData_R(:,1);
Data.RE_Pos_Y = rawData_R(:,2);
Data.RE_Pos_Z = rawData_R(:,3);
end