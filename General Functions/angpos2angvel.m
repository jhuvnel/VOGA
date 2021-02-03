function Data = angpos2angvel(data_rot,Data_In)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Adapted from PJB's voma function voma__processeyemovements.m
%
%   data_rot:
%       1: Apply no coordinate system transformations to raw data
%       2: Apply a -pi/2 YAW reorientation of the raw data
%       3: Apply a -pi/4 YAW reorientation of the raw data
%       4: Apply a -3*pi/4 YAW reorientation of the raw data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here, we are setting up the rotation needed to convert eye movement data
% presented in a COIL FRAME coordinate system, into a HEAD coordinate
% system. In particular for oculomotor data collected using a search coil
% system, in general the data computed from the raw coil signals provides
% data representing the rotations of the COILS on the eye in a COIL FRAME
% coordinate system.
% In square bracket rotation vector notation, we can say the system
% provides:
% [Coil_in_Frame].
% Typically, we care more about how the EYES are moving in a HEAD
% coordinate frame (i.e., [Eye_in_Head]. We will need to apply some coordinate system changes
% (also termed 'passive' rotations) to get there.
switch data_rot
    case 1
        % In cases where the coil system uses a HEAD FIXED coil frame, we often
        % assume that the HEAD coordinate system is identical to the FRAME
        % coordinate system. (Thus, we apply NO additional coordinate system
        % changes to the [Coil_in_Frame] data.
        radH = 0;
        radV = 0;
        radT = 0;
    case 2
        % Some cases in VNEL where this is NOT true:
        
        % In Ross710, the coil frame cardinal axes are defined as +X =
        % Left +Y = Occipital, +Z = Up convention (which we will refer to as
        % 'Frame')
        % We would like to use the standard coordinate system of +X = Nasal,
        % +Y = Left, +Z = Up (which we will refer to as 'Head')
        %
        % Thus, we can thus express Ross710 data in
        % 'Head' (the standard X-Y-Z convention) coordiantes by applying a -pi/2 Z-axis rotation
        %
        % To make this correction, we will convert the Coil-in-Frame rotation
        % vectors we get from the raw coil system voltages into Coil-in-Head
        % rotation vectors by applying a rotation describing
        % Head-in-Frame (described below in Fick angles)
        radH = -pi/2;
        radV = 0;
        radT = 0;
    case 3
        % When providing LARP or RALP mechanical rotations for a monkey in
        % the Ross710 Monkey Chamber, we must re-orient the monkey to be on
        % its SIDE (or its BACK, depending on how you seat the monkey in
        % the chair) and then physically apply a +pi/4 yaw rotation in HEAD
        % coordinates (so rotate the monkey's headcap within the 'ring'
        % headcap holder). For this system, our new [Head_in_Frame]
        % conversion is a -pi/4 yaw re-orientation to align the coordinate
        % systems.
        % NOTE: The correct head reorientation needed CHANGES if the monkey
        % is seated with its LEFT EAR to the chamber door, or if the BACK
        % OF THE HEAD is to the chamber door.
        % This case is used for:
        % Left Ear to Door: RALP rotation
        % Back of Head to Door: LARP rotation
        radH = -pi/4;
        radV = 0;
        radT = 0;
    case 4
        % When providing LARP or RALP mechanical rotations for a monkey in
        % the Ross710 Monkey Chamber, we must re-orient the monkey to be on
        % its SIDE (or its BACK, depending on how you seat the monkey in
        % the chair) and then physically apply a -pi/4 yaw rotation in HEAD
        % coordinates (so rotate the monkey's headcap within the 'ring'
        % headcap holder). For this system, our new [Head_in_Frame]
        % conversion is a -3pi/4 yaw re-orientation to align the coordinate
        % systems.
        % NOTE: The correct head reorientation needed CHANGES if the monkey
        % is seated with its LEFT EAR to the chamber door, or if the BACK
        % OF THE HEAD is to the chamber door.
        % This case is used for:
        % Left Ear to Door: LARP rotation
        % Back of Head to Door: RALP rotation
        radH = -3*pi/4;
        radV = 0;
        radT = 0;
end
HinFfick = [(radH*180/pi) (radV*180/pi) (radT*180/pi)]; % Fick angle in [H V T] convention
%% Extract and process the raw data
% LD VOG Goggles - MVI Trial
Fs = Data_In.Fs;
if (isrow(Data_In.Data_LE_Pos_X))
    Data_In.Data_LE_Pos_X=Data_In.Data_LE_Pos_X';
end
if (isrow(Data_In.Data_LE_Pos_Y))
    Data_In.Data_LE_Pos_Y=Data_In.Data_LE_Pos_Y';
end
if (isrow(Data_In.Data_LE_Pos_Z))
    Data_In.Data_LE_Pos_Z=Data_In.Data_LE_Pos_Z';
end
if (isrow(Data_In.Data_RE_Pos_X))
    Data_In.Data_RE_Pos_X=Data_In.Data_RE_Pos_X';
end
if (isrow(Data_In.Data_RE_Pos_Y))
    Data_In.Data_RE_Pos_Y=Data_In.Data_RE_Pos_Y';
end
if (isrow(Data_In.Data_RE_Pos_Z))
    Data_In.Data_RE_Pos_Z=Data_In.Data_RE_Pos_Z';
end
rawData_L = [ Data_In.Data_LE_Pos_X  Data_In.Data_LE_Pos_Y  Data_In.Data_LE_Pos_Z ];
rawData_R = [ Data_In.Data_RE_Pos_X  Data_In.Data_RE_Pos_Y  Data_In.Data_RE_Pos_Z ];
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

function ang=rot2angvel(rot)
dr = diff(rot,1,1);
r = rot(1:(end-1),:) + dr./2;
denom = (1 + dot(r,r,2) - dot(dr,dr,2));
ang = 2*(dr + cross(r,dr,2)) ./ denom(:,[1 1 1]);

function rot=fick2rot(fick)
% rot = fick2rot(fick)
% Generate nx3 rotation vectors from the given nx3 horizontal, vertical
% and torsional Fick angles in degrees.  +h left, +v down, +t clockwise.
% Get the tangents of half of the angles.
fick = tan(fick.*(pi/180/2));
% We make rotation vectors in the same way we would make rotation matricies,
% by simply applying the three Fick rotations (Euler "gimbal" angles) in the
% proper order.  The "straight ahead" rotation vector is just [0 0 0], and
% for a rotation about a single axis you just put the tan of half the angle
% of rotation in the X, Y, or Z component of the rotation vector.  So we
% create these three simple, single rotations, and combine them.
% We need a vector of zeros with the same number of rows as our input.
zr = zeros(size(fick,1),1);
% The X (torsion) rotation is the first and innermost rotation.
% The Y (vertical) is the middle rotation.
% The Z (horizontal) is the last and outermost rotation.
rot = rot2rot([zr zr fick(:,1)], rot2rot([zr fick(:,2) zr], [fick(:,3) zr zr]));

function EinH = raw2rot(data, gains, ref, option, HinF)
%RAW2ROT.M
%100896 version
%Purpose	Produces rotations vectors from digitized data.
%
%PJB 2017-06-15
%I am adding code to correctly process coil signals acquired where the HEAD
%and FRAME coordinate systems are not equivalent. In VNEL, the Lasker
%system head frame has been setup with a non-traditional coordinate system:
%+X = Coming out of the LEFT ear
%+Y = Coming out the BACK of the subject's head
%+Z = Coming out the TOP of the subject's head
%Thus, we will need to perform the following operations:
%
% Assumptions:
% - From the coil system, we get [Coil_in_Frame]
% - [Head] |= [Frame] (i.e, the head is NOT aligned with the frame
% coordinate system
% - [Head_in_Frame] = [0 0 -1]
% - [Coil_in_Eye] = [Coil_in_Head]_@t=0 (i.e., we are assuming the subject
% is looking forward at t=0)
%
% [Coil_in_Head] = [Head_in_Frame]^T[Coil_in_Frame]
% [Eye_in_Head] = [Coil_in_Head][Coil_in_Eye]^T (i.e., 'taking a reference
% position')
% If the user does not input a 'HinF'
%
%Description	Remark
%		The below description uses the right-hand rule with
%		the following Cartesian coordinate system:
%		x: about nasooccipital pointing forward
%		y: horizontal (interaural) pointing leftward
%		z: perpendicular to x and y pointing upward
%
%		Data format
%		Matrix with 6 cols:
%		(1) torsional field direction coil (x)
%		(2) vertical field direction coil (y)
%		(3) horizontal (interaural) field direction coil (z)
%		(4) torsional field torsion coil (x)
%		(5) vertical field torsion coil (y)
%		(6) horizontal (interaural) field torsion coil (z)
%
%		Data is A/D converted by a 12-bit converter. Numbers range
%		from 0 to 4096. Therefore 2048 have to be subtracted from
%		each data channel.
%
%		Gains format
%		Vector with 6 elements.
%		During in-vitro calibrations the signals are inverted such
%		that positive signals appear leftward, upward, and with
%		positive torsion (extorsion of the right eye, intorsion of
%		the left eye).
%		(1) Maximal signal with direction coil sensitivity vector
%		along x
%		(2) along y
%		(3) along z
%		(4) Maximal signal with torsion coil sensitivity vector
%		along x
%		(5) along y
%		(6) along z
%
%Input		data
%		gains
%		ref: reference position = 6 element vector or 0 (= first sample taken)
%		option (optional): 0 = default (compute, rotate, and plot), 1 = compute only.
%
%Output		rot:	rotation vectors (see Haustein 1989): 3 cols (x, y, z - components).
%		reg:	4 col vector (yslope, zslope, offset, stdev).
%
%Call		r = raw2rot(data, gains, option, ref_option)
%
('Checking arguments...');
if nargin < 2
    gains = [1 1 1 1 1 1];
    %  disp('  Bad argument number.');
    %  return;
end
if nargin < 3
    option = 0;
    ref = 0;
    HinF = 0;
end
if nargin < 4
    option = 0;
    HinF = 0;
end
if nargin < 5
    HinF = 0;
end
%disp('Checking input formats...');
[~,dc] = size(data);
if dc ~= 6
    disp('  Data matrix must have 6 cols.');
    return;
end
gr = size(gains,1);
if (gr~= 1)
    gains = gains';
end
[~,gc] = size(gains);
if (gc ~= 6) && (gr ~= 1)
    disp('  Gains vector must have 6 elements...');
    return;
end
%disp('Constructing reference sample');
if (length(ref) == 6)
    data = [ref; data];
end
%disp('Dividing all data channels with gains (sign correction included)...');
dxv = (data(:,1)) ./ gains(1);
dyv = (data(:,2)) ./ gains(2);
dzv = (data(:,3)) ./ gains(3);
txv = (data(:,4)) ./ gains(4);
tyv = (data(:,5)) ./ gains(5);
tzv = (data(:,6)) ./ gains(6);
clear data;
%Length of coil vectors
dlv = sqrt(dxv.^2 + dyv.^2 + dzv.^2);
%tlv = sqrt(txv.^2 + tyv.^2 + tzv.^2);
%Inproduct of both coil vectors
inpv = dxv.*txv + dyv.*tyv + dzv.*tzv;
%Orthogonalize torsional vector
toxv = txv - ((inpv .* dxv) ./ (dlv.^2));
clear txv;
toyv = tyv - ((inpv .* dyv) ./ (dlv.^2));
clear tyv;
tozv = tzv - ((inpv .* dzv) ./ (dlv.^2));
clear tzv inpv;
%Length of orthogonalized torsional vector
tolv = sqrt(toxv.^2 + toyv.^2 + tozv.^2);
%Construct instantaneous new orthogonal coordinate system of coils
f11 = dxv ./ dlv; clear dxv;
f12 = dyv ./ dlv; clear dyv;
f13 = dzv ./ dlv; clear dzv dlv;
f21 = toxv ./ tolv; clear toxv;
f22 = toyv ./ tolv; clear toyv;
f23 = tozv ./ tolv; clear tozv tolv;
f31 = f12.*f23-f13.*f22;
f32 = f13.*f21-f11.*f23;
f33 = f11.*f22-f12.*f21;
% Create rotation vectors from rotation matrix.
denom = 1+f11+f22+f33;
rot = [(f23-f32) (f31-f13) (f12-f21)] ./ denom(:,[1 1 1]);
% NOTE: At this point, we have constructed our time series of rotation
% vectors for the COIL in the FRAME, or [Coil_in_Frame]
CinF = rot;
if length(HinF) == 3
    % If the user sent a [3x1] rotation vector describing [Head_in_Frame],
    % apply that rotation.
    CinH = rot2rot(-HinF,CinF);
else
    % Used in cases where the user is assuming the [Frame] coordinate
    % system is identical to the [Head] coordinate system.
    CinH = CinF;
end
% Apply reference position (taken from first sample).
% i.e., Compute [Eye_in_Head] from [Coil_in_Head][Coil_in_Eye]^T
%
% Where:
% - [Coil_in_Eye] = [Coil_in_Head]_@t=0 (i.e., we are assuming the subject
% is looking forward at t=0)
CinE = CinH(1,:);
EinH = rot2rot(CinH,-CinE);
if (length(ref) == 6)
    EinH = EinH(2:length(EinH),:);
end

function nrot = rot2rot(rot1, rot2)
%ROT2ROT.M
%Order for creating eye in head = -head in space (tiltrot) o eye in space (rot)
%DSt
%
%Purpose
%Rotates rotation vectors by a constant rotation
%
%Call
%nrot = rot2rot(rot1, rot2);
%
%Arguments
%nrot: rotated rotation vectors
%rot1: input rotation vector or constant (single row) vector
%rot2: input rotation vector or constant (single row) vector
%
%  The rotation may be applied in either order, so one argument
% is the data vector (many rows) and the other is the constant
% rotation to be applied (single row).
%
%  DCR 7/22/97 Both arguments may now be many rows (# rows must
%              be equal in both args).  So it can, for instance,
%              be used to calculate eye in head from head in
%              space and eye in space.

if nargin ~= 2
    disp('Must have two arguments!')
    return
end
[~,cc] = size(rot1);
if cc ~= 3
    disp('First arg must have 3 columns');
    return;
end
[~,cc] = size(rot2);
if cc ~= 3
    disp('Second arg must have 3 columns');
    return
end
% Whichever argument is just a single row, expand it to the
% same number of rows as the other argument.
if size(rot1,1) == 1
    rot1 = ones(size(rot2,1),1) * rot1;
elseif size(rot2,1) == 1
    rot2 = ones(size(rot1,1),1) * rot2;
end
cr = cross(rot1',rot2')';
ir = dot(rot1', rot2')' * [1 1 1];
nrot = (rot1+rot2+cr)./(1-ir);

function fick = rot2fick(rot)
%Purpose Converts rotation vectors into Fick angles
%Algorithm from a C-program written by Th. Haslwanter in 1990.
% Essentially does rot2mat() then mat2fick(), but we only
% compute the 3 elements of the mat that we need for Fick.
x = rot(:,1);
y = rot(:,2);
z = rot(:,3);
clear rot;
scalar = sqrt(x.*x + y.*y + z.*z);
alpha = 2*atan(scalar);
% Avoid divide-by-zero errors.
zr = scalar == 0;
scalar(zr) = 1;
nx = x ./ scalar; %normalize
ny = y ./ scalar;
nz = z ./ scalar;
clear x y z
sn = sin(alpha);
cs = cos(alpha);
clear alpha
f31 = nz .* nx .* (1-cs) - ny .* sn;
f21 = ny .* nx .* (1-cs) + nz .* sn;
f32 = nz .* ny .* (1-cs) + nx .* sn;
clear sn cs
clear nx ny nz
% See Haslwanter et al., 1995 Eq. A4
vert = asin(f31);
horz = asin(f21 ./ cos(vert));
tors = asin(f32 ./ cos(vert));
fick = [horz -vert tors] .* 180 / pi;