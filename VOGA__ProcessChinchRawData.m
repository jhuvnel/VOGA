function VOGA__ProcessChinchRawData
savefile=1;
[chairfiles,filepaths] = uigetfile('*CHAIR*','Select Files With Chair Data.','MultiSelect','on');
if ischar(chairfiles)
    chairfiles = {chairfiles};  % Convert to cell array if only one file is selected
end
for numfiles = 1:length(chairfiles)
    chair_data = table2array(readtable([filepaths chairfiles{numfiles}]));
    chair_data(:,1:2) = chair_data(:,1:2)/16; %CDS122917 for EMAv66; this change accounts for the fact that Dale and Peter Boutros changed the line count/rev on May 2015
    
    vog_data = table2array(readtable([filepaths strrep(chairfiles{numfiles},'CHAIR','VOG')]));
    vog_data = [zeros(7,26); vog_data];

    fprintf('\n:::::::::::::::::::::::::::::::::::::::\n%s\n%s',[filepaths chairfiles{numfiles}],[filepaths strrep(chairfiles{numfiles},'CHAIR','VOG')])
    
    prompt = {'Left eye theta:','Left eye y_0:','Left eye z_0:','Right eye theta:','Right eye y_0:','Right eye z_0:'};
    dlgtitle = 'Please input calibration parameters';
    fieldsize = [1 45; 1 45];
    try
        try
            load(dir([strrep(chairfiles{numfiles},'CHAIR','') '_CAL']).name,'cal_data')
        catch
            if exist('cal_data','var')
                definput = {num2str(cal_data(1,3)),num2str(cal_data(1,1)),num2str(cal_data(1,2)),num2str(cal_data(2,3)),num2str(cal_data(2,1)),num2str(cal_data(2,2)),};
            else
                load(dir('*_CAL.mat').name,'cal_data')
            end
        end
    catch
        definput = {'65','0','0','-65','0','0'};
    end
    if exist('cal_data','var')
        definput = {num2str(cal_data(1,3)),num2str(cal_data(1,1)),num2str(cal_data(1,2)),num2str(cal_data(2,3)),num2str(cal_data(2,1)),num2str(cal_data(2,2)),};   
    end
    opts.Interpreter = 'tex';
    dlganswer = inputdlgcol(prompt,dlgtitle,1,definput,'on',2,[],[],opts);
    
    cal_data = [str2num(dlganswer{2}) str2num(dlganswer{3}) str2num(dlganswer{1}); str2num(dlganswer{5}) str2num(dlganswer{6}) str2num(dlganswer{4})];
    if savefile
        save([strrep(chairfiles{numfiles},'CHAIR','') '_CAL'],'cal_data')
    end
    fprintf('\n%s',[num2str(cal_data(1,:)) ' ' num2str(cal_data(2,:))])

    % Subtract y0, z0
    % Left eye
    vog_data(:,[8 10 12]) = vog_data(:,[8 10 12]) - cal_data(1,1);
    vog_data(:,[9 11 13]) = vog_data(:,[9 11 13]) - cal_data(1,2);
    % Right eye
    vog_data(:,[21 23 25]) = vog_data(:,[21 23 25]) - cal_data(2,1);
    vog_data(:,[22 24 26]) = vog_data(:,[22 24 26]) - cal_data(2,2);
    
    % Get real-time A values for each eye
    LE_A = vog_data(:,4);
    RE_A = vog_data(:,17);
    %% Select time reference
    % Build arrays of YZ coordinates for each point
    % Left eye
    LE_Point3 = [vog_data(:,12) vog_data(:,13)]';
    LE_Point2 = [vog_data(:,10) vog_data(:,11)]';
    LE_Point1 = [vog_data(:,8) vog_data(:,9)]';
    LE_time = vog_data(:,1);
    % Right eye
    RE_Point3 = [vog_data(:,25) vog_data(:,26)]';
    RE_Point2 = [vog_data(:,23) vog_data(:,24)]';
    RE_Point1 = [vog_data(:,21) vog_data(:,22)]';
    RE_time = vog_data(:,14);

    % Camera coordinates into head coordinates
    % Left eye
    Lthetad = cal_data(1,3);
    Lphir = 0;
    Lpsir = 0;
    Lthetar = Lthetad*pi/180;
    LCam2Head = [cos(Lthetar)*cos(Lphir), sin(Lthetar)*cos(Lphir), -sin(Lphir)
        cos(Lthetar)*sin(Lphir)*sin(Lpsir) - sin(Lthetar)*cos(Lpsir), sin(Lthetar)*sin(Lphir)*sin(Lpsir) + cos(Lthetar)*cos(Lpsir), cos(Lphir)*sin(Lpsir)
        cos(Lthetar)*sin(Lphir)*cos(Lpsir) + sin(Lthetar)*sin(Lpsir), sin(Lthetar)*sin(Lphir)*cos(Lpsir) - cos(Lthetar)*sin(Lpsir), cos(Lphir)*cos(Lpsir)]';
    % Right eye
    Rthetad = cal_data(2,3);
    Rphir = 0;
    Rpsir = 0;
    Rthetar = Rthetad*pi/180;
    RCam2Head = [cos(Rthetar)*cos(Rphir), sin(Rthetar)*cos(Rphir), -sin(Rphir)
        cos(Rthetar)*sin(Rphir)*sin(Rpsir) - sin(Rthetar)*cos(Rpsir), sin(Rthetar)*sin(Rphir)*sin(Rpsir) + cos(Rthetar)*cos(Rpsir), cos(Rphir)*sin(Rpsir)
        cos(Rthetar)*sin(Rphir)*cos(Rpsir) + sin(Rthetar)*sin(Rpsir), sin(Rthetar)*sin(Rphir)*cos(Rpsir) - cos(Rthetar)*sin(Rpsir), cos(Rphir)*cos(Rpsir)]';


    % Select reference (routine adapted from CDS SADA)
    eyeball_mm = 6.775; % chinchilla eyeball radius in mm
    dots_dist = 1; % minimum distance between dots in mm

    % Left eye
    reffig = figure;
    set(reffig, 'Position', get(0,'Screensize') - [0 0 0 80]);
    plot(LE_time,LE_Point1,'.',LE_time,LE_Point2,'.',LE_time,LE_Point3,'.');
    legend({'Point 1 Y','Point 1 Z','Point 2 Y','Point 2 Z','Point 3 Y','Point 3 Z'})
    hold on;
    xlabel('time (sec)');
    ylabel('realtime YZ coordinates for L eye');
    title(['LEFT EYE: ' strrep(chairfiles{numfiles},'_CHAIR','')], 'Press z to zoom and p to pan. click on plot with R mouse to select ref and middle scroll button to quit', 'Interpreter', 'none');
    enableDefaultInteractivity(gca);
    YL = get(gca, 'YLim');
    button = 1;
    ylim manual
    ylim(YL)
    existingref = 0;
    while button %continue until it is set to zero by case(2), indicating middle button hit
        [xi,~,button] = ginputzp(1);
        switch(button)
            case(3)
                if existingref
                    delete(refline)
                end
                existingref = 1;
                refline = plot([xi xi],[-1000 1000],'b-');
                [~,Lrefindx]=min(abs(LE_time-xi));
            case(2)
                if ~isempty(Lrefindx)
                    button=0;
                else
                    disp('have not yet hit button');
                end
        end
    end
    close(reffig)

    Lrefval = LE_time(Lrefindx);
    LE_Points = [LE_Point1; LE_Point2; LE_Point3];


    % Right eye
    reffig = figure;
    set(reffig, 'Position', get(0,'Screensize') - [0 0 0 80]);
    plot(RE_time,RE_Point1,'.',RE_time,RE_Point2,'.',RE_time,RE_Point3,'.');
    legend({'Point 1 Y','Point 1 Z','Point 2 Y','Point 2 Z','Point 3 Y','Point 3 Z'})
    hold on;
    xlabel('time (sec)');
    ylabel('realtime YZ coordinates for R eye');
    title(['RIGHT EYE: ' strrep(chairfiles{numfiles},'_CHAIR','')], 'Press z to zoom and p to pan. click on plot with R mouse to select ref and middle scroll button to quit', 'Interpreter', 'none');
    enableDefaultInteractivity(gca);
    YL = get(gca, 'YLim');
    button = 1;
    ylim manual
    ylim(YL)
    existingref = 0;
    while button %continue until it is set to zero by case(2), indicating middle button hit
        [xi,~,button] = ginputzp(1);
        switch(button)
            case(3)
                if existingref
                    delete(refline)
                end
                existingref = 1;
                refline = plot([xi xi],[-1000 1000],'b-');
                [~,Rrefindx]=min(abs(RE_time-xi));
            case(2)
                if ~isempty(Rrefindx)
                    button=0;
                else
                    disp('have not yet hit button');
                end
        end
    end
    close(reffig)

    Rrefval = RE_time(Rrefindx);
    RE_Points = [RE_Point1; RE_Point2; RE_Point3];
    %% Now calculate position
    % Left eye
    LE_eyeball_pix = eyeball_mm*min([sqrt(sum((LE_Point1(:,Lrefindx)-LE_Point3(:,Lrefindx)).^2)) sqrt(sum((LE_Point2(:,Lrefindx)-LE_Point3(:,Lrefindx)).^2))])/dots_dist;
    LE_pupil_baseline = LE_Points(:,Lrefindx);

    % Sort dot positions and get 3D rotation
    % LE_previous = LE_Points(:,Lrefindx); % note - this is how it was programmed in EMA v63, but does not seem correct. Instead, I think we should be calculating distance traveled from one time frame to next.
    LE_present_new = nan(size(LE_Points));
    LE_ValidCalc = nan(1,length(LE_Points));
    LE_T = nan(1,length(LE_Points));
    LE_V = nan(1,length(LE_Points));
    LE_H = nan(1,length(LE_Points));
    LE_present_new(:,1) = LE_Points(:,1);
    for i = 2:length(LE_Points)
%         LE_previous = LE_present_new(:,i-1);
        LE_previous = LE_Points(:,Lrefindx);
        LE_present = LE_Points(:,i);
        LE_lse_array = [norm(LE_previous(1:2)-LE_present(1:2))^2 + norm(LE_previous(3:4)-LE_present(3:4))^2 + norm(LE_previous(5:6)-LE_present(5:6))^2
            norm(LE_previous(1:2)-LE_present(1:2))^2 + norm(LE_previous(5:6)-LE_present(3:4))^2 + norm(LE_previous(3:4)-LE_present(5:6))^2
            norm(LE_previous(3:4)-LE_present(1:2))^2 + norm(LE_previous(1:2)-LE_present(3:4))^2 + norm(LE_previous(5:6)-LE_present(5:6))^2
            norm(LE_previous(3:4)-LE_present(1:2))^2 + norm(LE_previous(5:6)-LE_present(3:4))^2 + norm(LE_previous(1:2)-LE_present(5:6))^2
            norm(LE_previous(5:6)-LE_present(1:2))^2 + norm(LE_previous(3:4)-LE_present(3:4))^2 + norm(LE_previous(1:2)-LE_present(5:6))^2
            norm(LE_previous(5:6)-LE_present(1:2))^2 + norm(LE_previous(1:2)-LE_present(3:4))^2 + norm(LE_previous(3:4)-LE_present(5:6))^2];
        [~,LE_min_perm_idx]=min(LE_lse_array);
        switch LE_min_perm_idx
            case 1
                LE_present_new(1:2,i) = LE_present(1:2);
                LE_present_new(3:4,i) = LE_present(3:4);
                LE_present_new(5:6,i) = LE_present(5:6);
            case 2
                LE_present_new(1:2,i) = LE_present(1:2);
                LE_present_new(3:4,i) = LE_present(5:6);
                LE_present_new(5:6,i) = LE_present(3:4);
            case 3
                LE_present_new(1:2,i) = LE_present(3:4);
                LE_present_new(3:4,i) = LE_present(1:2);
                LE_present_new(5:6,i) = LE_present(5:6);
            case 4
                LE_present_new(1:2,i) = LE_present(5:6);
                LE_present_new(3:4,i) = LE_present(1:2);
                LE_present_new(5:6,i) = LE_present(3:4);
            case 5
                LE_present_new(1:2,i) = LE_present(5:6);
                LE_present_new(3:4,i) = LE_present(3:4);
                LE_present_new(5:6,i) = LE_present(1:2);
            case 6
                LE_present_new(1:2,i) = LE_present(3:4);
                LE_present_new(3:4,i) = LE_present(5:6);
                LE_present_new(5:6,i) = LE_present(1:2);
        end
        LE_previous = LE_Points(:,Lrefindx);
        LE_previous_eyeball = [sqrt(LE_eyeball_pix^2 - norm(LE_previous(1:2))^2);
            sqrt(LE_eyeball_pix^2 - norm(LE_previous(3:4))^2);
            sqrt(LE_eyeball_pix^2 - norm(LE_previous(5:6))^2)];
        LE_previous_eyeball_mat = [LE_previous_eyeball reshape(LE_previous,2,3)'];
        LE_cond1 = sum(sum(((LE_previous_eyeball_mat'*LE_previous_eyeball_mat)*pinv(LE_previous_eyeball_mat'*LE_previous_eyeball_mat) - eye(3,3)).^2)) < 0.01;

        LE_present_eyeball = [sqrt(LE_eyeball_pix^2 - norm(LE_present_new(1:2,i))^2);
            sqrt(LE_eyeball_pix^2 - norm(LE_present_new(3:4,i))^2);
            sqrt(LE_eyeball_pix^2 - norm(LE_present_new(5:6,i))^2)];
        LE_present_eyeball_mat = [LE_present_eyeball reshape(LE_present_new(:,i),2,3)'];
        LE_temp_mat = (LE_present_eyeball_mat'*LE_previous_eyeball_mat)*pinv(LE_previous_eyeball_mat'*LE_previous_eyeball_mat);
        LE_T(i) =180*asin(LE_temp_mat(3,2)/cos(asin(-LE_temp_mat(3,1))))/pi;
        LE_V(i) =180*asin(-LE_temp_mat(3,1))/pi;
        LE_H(i) =180*asin(LE_temp_mat(2,1)/cos(asin(-LE_temp_mat(3,1))))/pi;
        
        if imag(LE_T(i)) ~= 0
            LE_T(i) = nan;
        end
        if imag(LE_V(i)) ~= 0
            LE_V(i) = nan;
        end
        if imag(LE_H(i)) ~= 0
            LE_H(i) = nan;
        end

        LE_cond2 = LE_H(i)>-45 & LE_H(i)<45;
        LE_cond3 = LE_V(i)>-45 & LE_V(i)<45;
        LE_cond4 = LE_T(i)>-45 & LE_T(i)<45;
        LE_ValidCalc(i) = double(LE_cond1)*double(LE_cond2)*double(LE_cond3)*double(LE_cond4);
    end

    LE_T_head = nan(1,length(LE_T));
    LE_V_head = nan(1,length(LE_T));
    LE_H_head = nan(1,length(LE_T));
    for i = 1:length(LE_H)
        Lthetar2 = pi*LE_H(i)/180;
        Lphi2 = pi*LE_V(i)/180;
        Lpsi2 = pi*LE_T(i)/180;
        Lrotmat_cam = [cos(Lthetar2)*cos(Lphi2), sin(Lthetar2)*cos(Lphi2), -sin(Lphi2)
            cos(Lthetar2)*sin(Lphi2)*sin(Lpsi2) - sin(Lthetar2)*cos(Lpsi2), sin(Lthetar2)*sin(Lphi2)*sin(Lpsi2) + cos(Lthetar2)*cos(Lpsi2), cos(Lphi2)*sin(Lpsi2)
            cos(Lthetar2)*sin(Lphi2)*cos(Lpsi2) + sin(Lthetar2)*sin(Lpsi2), sin(Lthetar2)*sin(Lphi2)*cos(Lpsi2) - cos(Lthetar2)*sin(Lpsi2), cos(Lphi2)*cos(Lpsi2)]';
        Lrotcam_head = Lrotmat_cam*LCam2Head;
        LE_T_head(i) =180*asin(Lrotcam_head(3,2)/cos(asin(-Lrotcam_head(3,1))))/pi;
        LE_V_head(i) =180*asin(-Lrotcam_head(3,1))/pi;
        LE_H_head(i) =180*asin(Lrotcam_head(2,1)/cos(asin(-Lrotcam_head(3,1))))/pi;
    end

    LE_pos_HVT_head = [LE_H_head; LE_V_head; LE_T_head]';
    LE_pos_TVH_head = [LE_T_head; LE_V_head; LE_H_head]'; % this is the left eye position array

    %% Right eye
    RE_eyeball_pix = eyeball_mm*min([sqrt(sum((RE_Point1(:,Rrefindx)-RE_Point3(:,Rrefindx)).^2)) sqrt(sum((RE_Point2(:,Rrefindx)-RE_Point3(:,Rrefindx)).^2))])/dots_dist;
    RE_pupil_baseline = RE_Points(:,Rrefindx);

    % Sort dot positions and get 3D rotation
    % RE_previous = RE_Points(:,Rrefindx); % note - this is how it was programmed in EMA v63, but does not seem correct. Instead, I think we should be calculating distance traveled from one time frame to next.
    RE_present_new = nan(size(RE_Points));
    RE_ValidCalc = nan(1,length(RE_Points));
    RE_T = nan(1,length(RE_Points));
    RE_V = nan(1,length(RE_Points));
    RE_H = nan(1,length(RE_Points));
    RE_present_new(:,1) = RE_Points(:,1);
    for i = 2:length(RE_Points)
%         RE_previous = RE_present_new(:,i-1);
        RE_previous = RE_Points(:,Rrefindx);
        RE_present = RE_Points(:,i);
        RE_lse_array = [norm(RE_previous(1:2)-RE_present(1:2))^2 + norm(RE_previous(3:4)-RE_present(3:4))^2 + norm(RE_previous(5:6)-RE_present(5:6))^2
            norm(RE_previous(1:2)-RE_present(1:2))^2 + norm(RE_previous(5:6)-RE_present(3:4))^2 + norm(RE_previous(3:4)-RE_present(5:6))^2
            norm(RE_previous(3:4)-RE_present(1:2))^2 + norm(RE_previous(1:2)-RE_present(3:4))^2 + norm(RE_previous(5:6)-RE_present(5:6))^2
            norm(RE_previous(3:4)-RE_present(1:2))^2 + norm(RE_previous(5:6)-RE_present(3:4))^2 + norm(RE_previous(1:2)-RE_present(5:6))^2
            norm(RE_previous(5:6)-RE_present(1:2))^2 + norm(RE_previous(3:4)-RE_present(3:4))^2 + norm(RE_previous(1:2)-RE_present(5:6))^2
            norm(RE_previous(5:6)-RE_present(1:2))^2 + norm(RE_previous(1:2)-RE_present(3:4))^2 + norm(RE_previous(3:4)-RE_present(5:6))^2];
        [~,RE_min_perm_idx]=min(RE_lse_array);
        switch RE_min_perm_idx
            case 1
                RE_present_new(1:2,i) = RE_present(1:2);
                RE_present_new(3:4,i) = RE_present(3:4);
                RE_present_new(5:6,i) = RE_present(5:6);
            case 2
                RE_present_new(1:2,i) = RE_present(1:2);
                RE_present_new(3:4,i) = RE_present(5:6);
                RE_present_new(5:6,i) = RE_present(3:4);
            case 3
                RE_present_new(1:2,i) = RE_present(3:4);
                RE_present_new(3:4,i) = RE_present(1:2);
                RE_present_new(5:6,i) = RE_present(5:6);
            case 4
                RE_present_new(1:2,i) = RE_present(5:6);
                RE_present_new(3:4,i) = RE_present(1:2);
                RE_present_new(5:6,i) = RE_present(3:4);
            case 5
                RE_present_new(1:2,i) = RE_present(5:6);
                RE_present_new(3:4,i) = RE_present(3:4);
                RE_present_new(5:6,i) = RE_present(1:2);
            case 6
                RE_present_new(1:2,i) = RE_present(3:4);
                RE_present_new(3:4,i) = RE_present(5:6);
                RE_present_new(5:6,i) = RE_present(1:2);
        end
        RE_previous = RE_Points(:,Rrefindx);
        RE_previous_eyeball = [sqrt(RE_eyeball_pix^2 - norm(RE_previous(1:2))^2);
            sqrt(RE_eyeball_pix^2 - norm(RE_previous(3:4))^2);
            sqrt(RE_eyeball_pix^2 - norm(RE_previous(5:6))^2)];
        RE_previous_eyeball_mat = [RE_previous_eyeball reshape(RE_previous,2,3)'];
        RE_cond1 = sum(sum(((RE_previous_eyeball_mat'*RE_previous_eyeball_mat)*pinv(RE_previous_eyeball_mat'*RE_previous_eyeball_mat) - eye(3,3)).^2)) < 0.01;

        RE_present_eyeball = [sqrt(RE_eyeball_pix^2 - norm(RE_present_new(1:2,i))^2);
            sqrt(RE_eyeball_pix^2 - norm(RE_present_new(3:4,i))^2);
            sqrt(RE_eyeball_pix^2 - norm(RE_present_new(5:6,i))^2)];
        RE_present_eyeball_mat = [RE_present_eyeball reshape(RE_present_new(:,i),2,3)'];
        RE_temp_mat = (RE_present_eyeball_mat'*RE_previous_eyeball_mat)*pinv(RE_previous_eyeball_mat'*RE_previous_eyeball_mat);
        RE_T(i) =180*asin(RE_temp_mat(3,2)/cos(asin(-RE_temp_mat(3,1))))/pi;
        RE_V(i) =180*asin(-RE_temp_mat(3,1))/pi;
        RE_H(i) =180*asin(RE_temp_mat(2,1)/cos(asin(-RE_temp_mat(3,1))))/pi;

        if imag(RE_T(i)) ~= 0
            RE_T(i) = nan;
        end
        if imag(RE_V(i)) ~= 0
            RE_V(i) = nan;
        end
        if imag(RE_H(i)) ~= 0
            RE_H(i) = nan;
        end

        RE_cond2 = RE_H(i)>-45 & RE_H(i)<45;
        RE_cond3 = RE_V(i)>-45 & RE_V(i)<45;
        RE_cond4 = RE_T(i)>-45 & RE_T(i)<45;
        RE_ValidCalc(i) = double(RE_cond1)*double(RE_cond2)*double(RE_cond3)*double(RE_cond4);
    end

    RE_T_head = nan(1,length(RE_T));
    RE_V_head = nan(1,length(RE_T));
    RE_H_head = nan(1,length(RE_T));
    for i = 1:length(RE_H)
        Rthetar2 = pi*RE_H(i)/180;
        Rphi2 = pi*RE_V(i)/180;
        Rpsi2 = pi*RE_T(i)/180;
        Rrotmat_cam = [cos(Rthetar2)*cos(Rphi2), sin(Rthetar2)*cos(Rphi2), -sin(Rphi2)
            cos(Rthetar2)*sin(Rphi2)*sin(Rpsi2) - sin(Rthetar2)*cos(Rpsi2), sin(Rthetar2)*sin(Rphi2)*sin(Rpsi2) + cos(Rthetar2)*cos(Rpsi2), cos(Rphi2)*sin(Rpsi2)
            cos(Rthetar2)*sin(Rphi2)*cos(Rpsi2) + sin(Rthetar2)*sin(Rpsi2), sin(Rthetar2)*sin(Rphi2)*cos(Rpsi2) - cos(Rthetar2)*sin(Rpsi2), cos(Rphi2)*cos(Rpsi2)]';
        Rrotcam_head = Rrotmat_cam*RCam2Head;
        RE_T_head(i) =180*asin(Rrotcam_head(3,2)/cos(asin(-Rrotcam_head(3,1))))/pi;
        RE_V_head(i) =180*asin(-Rrotcam_head(3,1))/pi;
        RE_H_head(i) =180*asin(Rrotcam_head(2,1)/cos(asin(-Rrotcam_head(3,1))))/pi;
    end

    RE_pos_HVT_head = [RE_H_head; RE_V_head; RE_T_head]';
    RE_pos_TVH_head = [RE_T_head; RE_V_head; RE_H_head]'; % this is the left eye position array


    %% Need to double check interp routine
    neg_LE_time = find(diff(LE_time)<0);
    neg_RE_time = find(diff(RE_time)<0);

    if neg_LE_time > 0
        LE_time(neg_LE_time+1:end) = [];
        LE_pos_HVT_head(neg_LE_time+1:end,:) = [];
        LE_pos_TVH_head(neg_LE_time+1:end,:) = [];
    end

    if neg_RE_time > 0
        RE_time(neg_RE_time+1:end) = [];
        RE_pos_HVT_head(neg_RE_time+1:end,:) = [];
        RE_pos_TVH_head(neg_RE_time+1:end,:) = [];
    end

    zero_LE_time = find(diff(LE_time)==0);
    zero_RE_time = find(diff(RE_time)==0);
    
    if ~isempty(zero_LE_time)
        LE_time(zero_LE_time) = [];
        LE_pos_HVT_head(zero_LE_time,:) = [];
        LE_pos_TVH_head(zero_LE_time,:) = [];
    end

    if ~isempty(zero_RE_time)
        RE_time(zero_RE_time) = [];
        RE_pos_HVT_head(zero_RE_time,:) = [];
        RE_pos_TVH_head(zero_RE_time,:) = [];
    end

    % Now interpolate data
    interp_time = interp1(chair_data(:,6),chair_data(:,6),0:0.0055:chair_data(end,6)); % First interpolate time array from chair to have fixed sampling frequency -- note this is equivalent to just setting interp_time = 0:0.0055:chair_data(end,6)
    interp_position_RE = interp1(RE_time, RE_pos_TVH_head,interp_time,'spline');
    interp_position_LE = interp1(LE_time, LE_pos_TVH_head,interp_time,'spline');
    interp_chairpos = interp1(chair_data(:,6), chair_data(:,1),interp_time,'spline');
    interp_chairvel = interp1(chair_data(:,6), chair_data(:,2),interp_time,'spline');
    interp_Rsense = interp1(chair_data(:,6), chair_data(:,3),interp_time,'spline');
    
    if savefile
        % Now export data
        SineData = nan(11,length(interp_time));
        SineData(1,1:4) = [Lthetad cal_data(1,1) cal_data(1,2) Lrefindx];
        SineData(2,1:4) = [Rthetad cal_data(2,1) cal_data(2,2) Rrefindx];
        SineData(2,:) = interp_time;
        SineData(3,:) = interp_chairpos; % first col of chair data, chair pos
        SineData(4,:) = interp_chairvel; % second col of chair data, chair vel
        SineData(5,:) = interp_Rsense; % third col of chair data, Rsense
        SineData(6:8,:) = interp_position_LE'; % LE XYZ angular position
        SineData(9:11,:) = interp_position_RE'; % RE XYZ angular position

        mkdir('Presegmented Files')
        cd 'Presegmented Files'
        filename = [strrep(chairfiles{numfiles},'_CHAIR',''),'.mat'];
        if ~isempty(dir(filename))
            for napp=2:100
                filename = [strrep(filename,'.mat','')];
                newfilename = sprintf('%s_v%d.mat',filename,napp);
                if isempty(dir(newfilename))
                    save(newfilename,'SineData');
                    break;
                end
            end
        else
            save(filename,'SineData');
        end
        cd ../
    end
end

% %% Now calculate velocity
% 
% % Left eye
% LE_rot_head = fick2rot(interp1(LE_time, LE_pos_HVT_head,interp_time,'spline'));
% LE_quat_head = rot2quat(LE_rot_head);
% LE_quat_head_inv = quatinv(LE_quat_head);
% LE_time1 = [interp_time 0 0];
% LE_time2 = [0 0 interp_time];
% LE_quat_head1 = [LE_quat_head' zeros(4,2)];
% LE_quat_head2 = [zeros(4,2) LE_quat_head'];
% LE_quat_head_diff = (LE_quat_head1-LE_quat_head2)./repmat(LE_time1-LE_time2,4,1);
% LE_quat_head_diff(:,[1 2 end-1 end]) = [];
% LE_quat_head_diff = [LE_quat_head_diff(1:4,1) LE_quat_head_diff LE_quat_head_diff(1:4,end)]; % slight diff in time array (suspect it's due to rounding) that carries over
% LE_quat_head_mult = (360/pi)*quatmultiply(LE_quat_head_diff',LE_quat_head_inv);
% LE_quat_head_mult(:,1)=[];
% LE_vel_TVH = nan(length(LE_quat_head_mult),3);
% for i = 1:length(LE_quat_head_mult)
%     LE_vel_TVH(i,:) = (LCam2Head*LE_quat_head_mult(i,:)')';
% end
% LE_vel_TVH = LE_vel_TVH';
% LE_vel_LRZ = [cos(-pi/4)*(LE_vel_TVH(1,:) - LE_vel_TVH(2,:));
%     cos(-pi/4)*(LE_vel_TVH(1,:) + LE_vel_TVH(2,:));
%     LE_vel_TVH(3,:)]; 
% 
% % Right eye
% RE_rot_head = fick2rot(interp1(RE_time, RE_pos_HVT_head,interp_time,'spline'));
% RE_quat_head = rot2quat(RE_rot_head);
% RE_quat_head_inv = quatinv(RE_quat_head);
% RE_time1 = [interp_time 0 0];
% RE_time2 = [0 0 interp_time];
% RE_quat_head1 = [RE_quat_head' zeros(4,2)];
% RE_quat_head2 = [zeros(4,2) RE_quat_head'];
% RE_quat_head_diff = (RE_quat_head1-RE_quat_head2)./repmat(RE_time1-RE_time2,4,1);
% RE_quat_head_diff(:,[1 2 end-1 end]) = [];
% RE_quat_head_diff = [RE_quat_head_diff(1:4,1) RE_quat_head_diff RE_quat_head_diff(1:4,end)]; % slight diff in time array (suspect it's due to rounding) that carries over
% RE_quat_head_mult = (360/pi)*quatmultiply(RE_quat_head_diff',RE_quat_head_inv);
% RE_quat_head_mult(:,1)=[];
% RE_vel_TVH = nan(length(RE_quat_head_mult),3);
% for i = 1:length(RE_quat_head_mult)
%     RE_vel_TVH(i,:) = (RCam2Head*RE_quat_head_mult(i,:)')';
% end
% RE_vel_TVH = RE_vel_TVH';
% RE_vel_LRZ = [cos(-pi/4)*(RE_vel_TVH(1,:) - RE_vel_TVH(2,:));
%     cos(-pi/4)*(RE_vel_TVH(1,:) + RE_vel_TVH(2,:));
%     RE_vel_TVH(3,:)]; 
end