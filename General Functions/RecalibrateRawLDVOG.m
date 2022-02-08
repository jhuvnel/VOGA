function RecalibrateRawLDVOG
Raw_Path = [cd,filesep,'Raw Files'];
file_names = extractfield(dir(Raw_Path),'name');
file_names(~contains(file_names,'.txt')) = [];
if isempty(file_names)
    file_names = {''};
end
% Load LE and RE calibration files, and distance to the wall
[LEcalibFileName,LEcalibPathName] = uigetfile('.txt','Select LEFT EYE calibration file.');
[REcalibFileName,REcalibPathName] = uigetfile('.txt','Select RIGHT EYE calibration file.');
calib_file_fullpath_left = strcat(LEcalibPathName, LEcalibFileName);
calib_file_fullpath_right = strcat(REcalibPathName, REcalibFileName);
prompt = {'Enter the distance from the subject to the calibration grid [in]'};
title = 'Distance To Wall (in.)';
definput = {'55'};
opts.Interpreter = 'tex';
dist_to_wall = inputdlg(prompt,title,[1 40],definput,opts);
DistanceToWall = str2double(dist_to_wall);
Polynomials_Left = VOG_Calibration_9_Points_AAedit20220208(DistanceToWall,calib_file_fullpath_left);
H_Left_Coeffs = Polynomials_Left(:,1);
V_Left_Coeffs = Polynomials_Left(:,2);
Polynomials_Right = VOG_Calibration_9_Points_AAedit20220208(DistanceToWall,calib_file_fullpath_right);
H_Right_Coeffs = Polynomials_Right(:,1);
V_Right_Coeffs = Polynomials_Right(:,2);
% These are the column indices of the relevant parameters saved to file
TIndex = 2;
HLeftIndex_Pix = 3;
VLeftIndex_Pix = 4;
HRightIndex_Pix = 15;
VRightIndex_Pix = 16;
HLeftIndex_Deg = 40;
VLeftIndex_Deg = 41;
HRightIndex_Deg = 43;
VRightIndex_Deg = 44;
% Select which file(s) to recalibrate
indx = nmlistdlg('PromptString','Select files to segment:','ListSize',[300 300],'ListString',file_names);
sel_files = file_names(indx);
for f = 1:length(sel_files)
    % Load Data from File
    data1 = readtable([Raw_Path,filesep,sel_files{f}]);
    vogdata = data1{:,contains(varfun(@class,data1,'OutputFormat','cell'),'double')};
    % load eye positions in degrees
    H_LE_pix = vogdata(:,HLeftIndex_Pix);
    V_LE_pix = vogdata(:,VLeftIndex_Pix);
    H_RE_pix = vogdata(:,HRightIndex_Pix);
    V_RE_pix = vogdata(:,VRightIndex_Pix);
    H_LE_deg = vogdata(:,HLeftIndex_Deg);
    V_LE_deg = vogdata(:,VLeftIndex_Deg);
    H_RE_deg = vogdata(:,HRightIndex_Deg);
    V_RE_deg = vogdata(:,VRightIndex_Deg);
    N = length(vogdata(:,TIndex));
    H_LE_deg_new = zeros(size(H_LE_deg));
    V_LE_deg_new = zeros(size(H_LE_deg));
    H_RE_deg_new = zeros(size(H_LE_deg));
    V_RE_deg_new = zeros(size(H_LE_deg));
    for i = 1:N
        H_LE_deg_new(i) = Polynomial_Surface_Mult(H_Left_Coeffs, H_LE_pix(i), V_LE_pix(i));
        V_LE_deg_new(i) = Polynomial_Surface_Mult(V_Left_Coeffs, H_LE_pix(i), V_LE_pix(i));
        H_RE_deg_new(i) = Polynomial_Surface_Mult(H_Right_Coeffs, H_RE_pix(i), V_RE_pix(i));
        V_RE_deg_new(i) = Polynomial_Surface_Mult(V_Right_Coeffs, H_RE_pix(i), V_RE_pix(i));
    end
    % The raw data in pixels will saturate when the subject blinks/tracking
    % is lost. The VOG software will replace these data points w/ NaNs
    % since there is no real tracking measurements. We will replace these
    % saturated values in our newly calibrated data with NaNs.
    H_LE_deg_new(isnan(H_LE_deg)) = nan(length(H_LE_deg_new(isnan(H_LE_deg))),1);
    V_LE_deg_new(isnan(V_LE_deg)) = nan(length(V_LE_deg_new(isnan(V_LE_deg))),1);
    H_RE_deg_new(isnan(H_RE_deg)) = nan(length(H_RE_deg_new(isnan(H_RE_deg))),1);
    V_RE_deg_new(isnan(V_RE_deg)) = nan(length(V_RE_deg_new(isnan(V_RE_deg))),1);
    vogdata(:,HLeftIndex_Deg) = H_LE_deg_new;
    vogdata(:,VLeftIndex_Deg) = V_LE_deg_new;
    vogdata(:,HRightIndex_Deg) = H_RE_deg_new;
    vogdata(:,VRightIndex_Deg) = V_RE_deg_new;
    fname_new = [Raw_Path,filesep,sel_files{f}(1:end-4),'_UpdatedVOGCalib_',num2str(DistanceToWall),'in__',datestr(now,'yyyy-mm-dd'),'.txt'];
    dlmwrite(fname_new, vogdata, 'delimiter',' ','precision','%.4f');
end
end


