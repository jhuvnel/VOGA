function [Data,Data_In] = angpos2angvel(Data_In)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapted by CFB from EMA by AAM. Important note: this scripts assumes that
% chinch VOG cam phi = 0 and psi = 0 for both eyes, and thus only passes
% theta values to compute rotation matrix to go from camera to head coords.
% If this changes, MAKE SURE TO CHANGE Lcamrotmatrix and Rcamrotmatrix.
%% Extract and process the raw data
Fs = Data_In.Fs;
L_cam_theta = Data_In.info.theta_L_cam;
R_cam_theta = Data_In.info.theta_R_cam;
Lcamrotmatrix = [cosd(L_cam_theta) -sind(L_cam_theta) 0;sind(L_cam_theta) cosd(L_cam_theta) 0;0 0 1];
Rcamrotmatrix = [cosd(R_cam_theta) -sind(R_cam_theta) 0;sind(R_cam_theta) cosd(R_cam_theta) 0;0 0 1];
if isfield(Data_In,'t')
    t = Data_In.t';
elseif isfield(Data_In,'Time_Eye')
    t = Data_In.Time_Eye';
end

% Left eye
LE_rot_head = fick2rot([Data_In.LE_Position_Z Data_In.LE_Position_Y Data_In.LE_Position_X]);
LE_quat_head = rot2quat(LE_rot_head);
LE_quat_head_inv = quatinv(LE_quat_head);
time1 = [t 0 0];
time2 = [0 0 t];
LE_quat_head1 = [LE_quat_head' zeros(4,2)];
LE_quat_head2 = [zeros(4,2) LE_quat_head'];
LE_quat_head_diff = (LE_quat_head1-LE_quat_head2)./repmat(time1-time2,4,1);
LE_quat_head_diff(:,[1 2 end-1 end]) = [];
LE_quat_head_diff = [LE_quat_head_diff(1:4,1) LE_quat_head_diff LE_quat_head_diff(1:4,end)]; % slight diff in time array (suspect it's due to rounding) that carries over
LE_quat_head_mult = (360/pi)*quatmultiply(LE_quat_head_diff',LE_quat_head_inv);
LE_quat_head_mult(:,1)=[];
LE_vel_TVH = nan(length(LE_quat_head_mult),3);
for i = 1:length(LE_quat_head_mult)
    LE_vel_TVH(i,:) = (Lcamrotmatrix*LE_quat_head_mult(i,:)')';
end
LE_vel_TVH = LE_vel_TVH';
LE_vel_LRZ = [cos(-pi/4)*(LE_vel_TVH(1,:) - LE_vel_TVH(2,:));
    cos(-pi/4)*(LE_vel_TVH(1,:) + LE_vel_TVH(2,:));
    LE_vel_TVH(3,:)]; 

% Right eye
RE_rot_head = fick2rot([Data_In.RE_Position_Z Data_In.RE_Position_Y Data_In.RE_Position_X]);
RE_quat_head = rot2quat(RE_rot_head);
RE_quat_head_inv = quatinv(RE_quat_head);
time1 = [t 0 0];
time2 = [0 0 t];
RE_quat_head1 = [RE_quat_head' zeros(4,2)];
RE_quat_head2 = [zeros(4,2) RE_quat_head'];
RE_quat_head_diff = (RE_quat_head1-RE_quat_head2)./repmat(time1-time2,4,1);
RE_quat_head_diff(:,[1 2 end-1 end]) = [];
RE_quat_head_diff = [RE_quat_head_diff(1:4,1) RE_quat_head_diff RE_quat_head_diff(1:4,end)]; % slight diff in time array (suspect it's due to rounding) that carries over
RE_quat_head_mult = (360/pi)*quatmultiply(RE_quat_head_diff',RE_quat_head_inv);
RE_quat_head_mult(:,1)=[];
RE_vel_TVH = nan(length(RE_quat_head_mult),3);
for i = 1:length(RE_quat_head_mult)
    RE_vel_TVH(i,:) = (Rcamrotmatrix*RE_quat_head_mult(i,:)')';
end
RE_vel_TVH = RE_vel_TVH';
RE_vel_LRZ = [cos(-pi/4)*(RE_vel_TVH(1,:) - RE_vel_TVH(2,:));
    cos(-pi/4)*(RE_vel_TVH(1,:) + RE_vel_TVH(2,:));
    RE_vel_TVH(3,:)]; 

% Store Data
Data.LE_Vel_X = LE_vel_TVH(1,:)';
Data.LE_Vel_Y = LE_vel_TVH(2,:)';
Data.LE_Vel_Z = LE_vel_TVH(3,:)';
Data.LE_Vel_LARP = LE_vel_LRZ(1,:)';
Data.LE_Vel_RALP = LE_vel_LRZ(2,:)';
Data.RE_Vel_X = RE_vel_TVH(1,:)';
Data.RE_Vel_Y = RE_vel_TVH(2,:)';
Data.RE_Vel_Z = RE_vel_TVH(3,:)';
Data.RE_Vel_LARP = RE_vel_LRZ(1,:)';
Data.RE_Vel_RALP = RE_vel_LRZ(2,:)';
% Saving both POSITION and RAW data traces
Data.LE_Pos_X = Data_In.LE_Position_X;
Data.LE_Pos_Y = Data_In.LE_Position_Y;
Data.LE_Pos_Z = Data_In.LE_Position_Z;
Data.RE_Pos_X = Data_In.RE_Position_X;
Data.RE_Pos_Y = Data_In.RE_Position_Y;
Data.RE_Pos_Z = Data_In.RE_Position_Z;

if exist('ElecStimTrig','var')
    Data.ElecStimTrig = ElecStimTrig;
end
Data.Fs = Fs;
% -------------------------------------------------------------------------
% Below find AIA's code based on PJB's code, no longer using for chinch VOG
% analysis due to concerns with gimbal lock caused by using Fick coords
% rather than quaternions
% -------------------------------------------------------------------------
% %Has Left Eye Data
% if isfield(Data_In,'LE_Position_Z')&&isfield(Data_In,'LE_Position_Y')&&isfield(Data_In,'LE_Position_X')
%     Data_In.LE_Position_X = reshape(Data_In.LE_Position_X,[],1);
%     Data_In.LE_Position_Y = reshape(Data_In.LE_Position_Y,[],1);
%     Data_In.LE_Position_Z = reshape(Data_In.LE_Position_Z,[],1);
%     rawData = [Data_In.LE_Position_X,Data_In.LE_Position_Y,Data_In.LE_Position_Z];
%     %Get LARP and RALP in position too
%     angpos_L = (rotZ3deg(-45)'*rawData')'; 
%     Data_In.LE_Position_LARP = angpos_L(:,1);
%     Data_In.LE_Position_RALP = angpos_L(:,2);    
%     % Computing angular velocity from Fick angular position angles
%     psi = rawData(:,1);
%     phi = rawData(:,2);
%     theta = rawData(:,3);
%     
%     angvel_dps_a = [(gradient(psi)*Fs.*cosd(theta).*cosd(phi)) - (gradient(phi)*Fs.*sind(theta)) ...
%         (gradient(psi)*Fs.*sind(theta).*cosd(phi)) + (gradient(phi)*Fs.*cosd(theta)) ...
%         (gradient(theta)*Fs) - (gradient(psi)*Fs.*sind(phi))];      
% 
%     Lcamrotmatrix = [cosd(L_cam_theta) -sind(L_cam_theta) 0;sind(L_cam_theta) cosd(L_cam_theta) 0;0 0 1];
%     angvel_dps_b = (Lcamrotmatrix*angvel_dps_a')';
% 
%     % We want to rotate our data from an [X,Y,Z] coordinate system,
%     % into a [LARP,RALP,LHRH] coordinate system, where Z = LHRH. To
%     % accomplish this, we will perform a PASSIVE (i.e., 'alias' or
%     % 'coordinate system') rotation of -45 degrees. In order to realize
%     % this rotation, I will generate a -45deg rotation matrix and
%     % RIGHT multiple our data by the TRANSPOSE of this rotation
%     % matrix. Note that rotation matrices are orthonormal, and
%     % their inverses are equivalent to their transpose. -PJB
%     angvel_dps_c = (rotZ3deg(-45)'*angvel_dps_b')';  
%     % Store Data
%     Data.LE_Vel_X = angvel_dps_b(:,1);
%     Data.LE_Vel_Y = angvel_dps_b(:,2);
%     Data.LE_Vel_Z = angvel_dps_b(:,3);
%     Data.LE_Vel_LARP = angvel_dps_c(:,1);
%     Data.LE_Vel_RALP = angvel_dps_c(:,2); 
%     % Saving both POSITION and RAW data traces
%     Data.LE_Pos_X = rawData(:,1);
%     Data.LE_Pos_Y = rawData(:,2);
%     Data.LE_Pos_Z = rawData(:,3);
% end
% 
% %Has Right Eye Data
% if isfield(Data_In,'RE_Position_Z')&&isfield(Data_In,'RE_Position_Y')&&isfield(Data_In,'RE_Position_X')
%     Data_In.RE_Position_X = reshape(Data_In.RE_Position_X,[],1);
%     Data_In.RE_Position_Y = reshape(Data_In.RE_Position_Y,[],1);
%     Data_In.RE_Position_Z = reshape(Data_In.RE_Position_Z,[],1);
%     rawData = [Data_In.RE_Position_X,Data_In.RE_Position_Y,Data_In.RE_Position_Z];
%     %Get LARP and RALP in position too
%     angpos_R = (rotZ3deg(-45)'*rawData')'; 
%     Data_In.RE_Position_LARP = angpos_R(:,1);
%     Data_In.RE_Position_RALP = angpos_R(:,2);    
%     % Computing angular velocity from Fick angular position angles
%     psi = rawData(:,1);
%     phi = rawData(:,2);
%     theta = rawData(:,3);
%     angvel_dps_a = [(gradient(psi)*Fs.*cosd(theta).*cosd(phi)) - (gradient(phi)*Fs.*sind(theta)) ...
%         (gradient(psi)*Fs.*sind(theta).*cosd(phi)) + (gradient(phi)*Fs.*cosd(theta)) ...
%         (gradient(theta)*Fs) - (gradient(psi)*Fs.*sind(phi))];    
%     Rcamrotmatrix = [cosd(R_cam_theta) -sind(R_cam_theta) 0;sind(R_cam_theta) cosd(R_cam_theta) 0;0 0 1];
%     angvel_dps_b = (Rcamrotmatrix*angvel_dps_a')';
%     % We want to rotate our data from an [X,Y,Z] coordinate system,
%     % into a [LARP,RALP,LHRH] coordinate system, where Z = LHRH. To
%     % accomplish this, we will perform a PASSIVE (i.e., 'alias' or
%     % 'coordinate system') rotation of -45 degrees. In order to realize
%     % this rotation, I will generate a -45deg rotation matrix and
%     % RIGHT multiple our data by the TRANSPOSE of this rotation
%     % matrix. Note that rotation matrices are orthonormal, and
%     % their inverses are equivalent to their transpose. -PJB
%     angvel_dps_c = (rotZ3deg(-45)'*angvel_dps_b')';  
%     % Store Data
%     Data.RE_Vel_X = angvel_dps_b(:,1);
%     Data.RE_Vel_Y = angvel_dps_b(:,2);
%     Data.RE_Vel_Z = angvel_dps_b(:,3);
%     Data.RE_Vel_LARP = angvel_dps_c(:,1);
%     Data.RE_Vel_RALP = angvel_dps_c(:,2); 
%     % Saving both POSITION and RAW data traces
%     Data.RE_Pos_X = rawData(:,1);
%     Data.RE_Pos_Y = rawData(:,2);
%     Data.RE_Pos_Z = rawData(:,3);
% end
end