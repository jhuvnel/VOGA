function plotRawVOG(Raw_Path,plot_eyes,xyz_lrz)
    %% Standardize Colors
    % Normal colors
    colors.l_x = [237,150,33]/255;
    colors.l_y = [125,46,143]/255;
    colors.l_z = [1 0 0];
    %colors.l_l = [0,128,0]/255;
    %colors.l_r = [0 0 1];
    colors.r_x = [237,204,33]/255;
    colors.r_y = [125,46,230]/255;
    colors.r_z = [1,0,1];
    %colors.r_l = [0 1 0];
    %colors.r_r = [64,224,208]/255;
    %% Select File
    if nargin < 3
        xyz_lrz = 'xyz';
    end
    if nargin < 2
        plot_eyes = 0;
    end   
    if nargin < 1
        Raw_Path = cd;
    end
    temp = dir(Raw_Path);    
    all_files = {temp.name}';
    all_files = all_files(~logical([temp.isdir]));
    VOG_files = all_files(((contains(all_files,'SESSION')&contains(all_files,'.txt'))|contains(all_files,'.dat')|contains(all_files,'.mat'))&~contains(all_files,'Notes'));
    [indx,tf] = nmlistdlg('PromptString','Select files to plot:',...
                               'ListSize',[300 300],...
                               'ListString',VOG_files,...
                               'SelectionMode','Single');
    if tf~=1
        return;
    end
    %% Make Plot
    file = [Raw_Path,filesep,VOG_files{indx}];
    if contains(file,'.dat') %NKI
        warning('off')
        data = readtable(file,'ReadVariableNames',true);
        warning('on')
        data.Properties.VariableNames{1} = 'EyeTime';
        Time_Eye = data.EyeTime;    
        Stim = zeros(length(Time_Eye),1);
        Stim(data.EventCode ~= 0) = 1;
        LZ = -data.LeftHoriz;
        LY = -data.LeftVert;
        LX = data.LeftTorsion;
        RZ = -data.RightHoriz;
        RY = -data.RightVert;
        RX = data.RightTorsion;
        GyroX = data.GyroY - median(data.GyroY);
        GyroY = -(data.GyroX- median(data.GyroX));
        GyroZ = -(data.GyroZ- median(data.GyroZ));
    elseif contains(file,'.mat') %NKI
        load(file,'Data')
        data = Data;
        Time_Eye = data.EyeTime;    
        Stim = zeros(length(Time_Eye),1);
        Stim(data.EventCode ~= 0) = 1;
        LZ = -data.LeftHoriz;
        LY = -data.LeftVert;
        LX = data.LeftTorsion;
        RZ = -data.RightHoriz;
        RY = -data.RightVert;
        RX = data.RightTorsion;
        GyroX = data.GyroY - median(data.GyroY);
        GyroY = -(data.GyroX- median(data.GyroX));
        GyroZ = -(data.GyroZ- median(data.GyroZ));
    else %LDVOG
        data = readtable(file);
        % Generate Time_Eye vector
        Time = data{:,2};
        Time_Eye = (0:length(Time)-1)'*median(diff(Time));
        % Index for the LDVOG lines
        StimIndex = 35; %GPIO line
        XvelHeadIndex = 30;
        YvelHeadIndex = 29;
        ZvelHeadIndex = 28;
        Stim = data{1:length(Time_Eye),StimIndex};
        GyroX = data{1:length(Time_Eye),XvelHeadIndex} - median(data{1:length(Time_Eye),XvelHeadIndex});
        GyroY = data{1:length(Time_Eye),YvelHeadIndex} - median(data{1:length(Time_Eye),YvelHeadIndex});
        GyroZ = data{1:length(Time_Eye),ZvelHeadIndex} - median(data{1:length(Time_Eye),ZvelHeadIndex});
        phi = -170;
        Rotation_Head = [
                cosd(phi) 0   sind(phi);
                0   1   0;
                -sind(phi)    0   cosd(phi)
                ];        
        % NOTE: We are transposing the rotation matrix in order to apply a
        % PASSIVE (i.e., a coordinate system) transformation
        A = Rotation_Head' * [GyroX' ; GyroY' ; GyroZ'];
        GyroX = A(1,:);
        GyroY = A(2,:);
        GyroZ = A(3,:);
        HLeftIndex = 40;
        VLeftIndex = 41;
        TLeftIndex = 42;
        HRightIndex = 43;
        VRightIndex = 44;
        TRightIndex = 45;       
        % Load raw eye position data in Fick coordinates [degrees]
        LZ = data{:,HLeftIndex};
        LY = data{:,VLeftIndex};
        LX = data{:,TLeftIndex};
        RZ = data{:,HRightIndex};
        RY = data{:,VRightIndex};
        RX = data{:,TRightIndex}; 
    end
    GyroL = (GyroX - GyroY)/sqrt(2);
    GyroR = (GyroX + GyroY)/sqrt(2);
    figure;
    if any(strcmp(xyz_lrz,{'xyz','XYZ'}))
        plot(Time_Eye,GyroX,'k:',Time_Eye,GyroY,'k--',Time_Eye,GyroZ,'k-',Time_Eye,40*Stim-30,'b')
        leg_text = {'GyroX','GyroY','GyroZ'};
    else
        leg_text = {'GyroL','GyroR','GyroZ'};
        plot(Time_Eye,GyroL,'k:',Time_Eye,GyroR,'k--',Time_Eye,GyroZ,'k-',Time_Eye,40*Stim-30,'b')
    end
    if plot_eyes
        hold on
        plot(Time_Eye,LX,'Color',colors.l_x)
        plot(Time_Eye,RX,'Color',colors.r_x)
        plot(Time_Eye,LY,'Color',colors.l_y)
        plot(Time_Eye,RY,'Color',colors.r_y)
        plot(Time_Eye,LZ,'Color',colors.l_z)
        plot(Time_Eye,RZ,'Color',colors.r_z)
        hold off
        legend([leg_text,{'Trigger','LX','RX','LY','RY','LZ','RZ'}])
    else
        legend([leg_text,{'Trigger'}])
    end        
    xlabel('Time (s)')
    ylabel('Velocity (dps)')
    title(strrep(strrep(VOG_files{indx},'_',' '),'-',' '))        
end