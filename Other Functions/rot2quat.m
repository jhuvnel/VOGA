function q=rot2quat(rot)

% quat = rot2quat(rot)
%
%  Convert 3-element rotation vectors to 4-element quaternions.
% rot must be an Nx3 array of rotation vectors, and quat will be an Nx4
% array of quaternions.

% Get the magnitude of each rotation vector.
mag = sqrt(dot(rot,rot,2));

% The magnitude of each rotation vector is just the tangent of half
% of the rotation angle.  Normally we would multiply by 2 to get the
% angle, but quaternions also use the half angle, so we just leave it
% as half the angle here, and take its cosine.
cos_half = cos(atan(mag));

% q0 is just the cosine of the half-angle, and the magnitude of the
% vector part of the quaternion ([q1 q2 q3]) is just the sine of the
% half-angle.  Since the rotation vector's magnitude is the tangent
% of the half-angle, we just multiply it by the cosine to get the
% sine.  i.e., sin(a) = tan(a) * cos(a).
q = [cos_half rot.*cos_half(:,[1 1 1])];
