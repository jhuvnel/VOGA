%% Static Tilt Measurements from VOG
% Standardize Colors
load('VNELcolors.mat','colors')
clc;
%% Select File
Raw_Path = [cd,filesep,'Raw Files'];
all_files = extractfield(dir(Raw_Path),'name',find(~extractfield(dir(Raw_Path),'isdir')));
LDVOG_files = all_files(contains(all_files,'SESSION')&contains(all_files,'.txt')&~contains(all_files,'Notes.txt'));
NKI_files = all_files(contains(all_files,'.dat'));
GNO_files = all_files(contains(all_files,{'Lateral.txt','LARP.txt','RALP.txt'}));
VOG_files = [LDVOG_files;NKI_files;GNO_files];
[indx,tf] = nmlistdlg('PromptString','Select files to analyze:',...
    'ListSize',[300 300],...
    'ListString',VOG_files,...
    'SelectionMode','multiple');
if tf~=1
    return;
end
for i = 1:length(indx)
    %% Make Plot
    figure('Color',[1 1 1],'units','inches','Position',[1 1 8 8]);
    file = [Raw_Path,filesep,VOG_files{indx(i)}];
    if ismember(VOG_files{indx(i)},NKI_files)
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
        plot(Time_Eye,50*(Stim-0.5),'b')
        hold on
        plot(Time_Eye,LX,'Color',colors.l_x)
        plot(Time_Eye,RX,'Color',colors.r_x)
        plot(Time_Eye,LY,'Color',colors.l_y)
        plot(Time_Eye,RY,'Color',colors.r_y)
        plot(Time_Eye,LZ,'Color',colors.l_z)
        plot(Time_Eye,RZ,'Color',colors.r_z)
        hold off
        legend({'Trigger','LX','RX','LY','RY','LZ','RZ'})
    elseif ismember(VOG_files{indx(i)},LDVOG_files)
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
        GyroL = (GyroX - GyroY)/sqrt(2);
        GyroR = (GyroX + GyroY)/sqrt(2);
        if any(strcmp(lrz_xyz,{'xyz','XYZ'}))
            plot(Time_Eye,GyroX,'k:',Time_Eye,GyroY,'k--',Time_Eye,GyroZ,'k-',Time_Eye,40*Stim-30,'b')
            leg_text = {'GyroX','GyroY','GyroZ'};
        else
            leg_text = {'GyroL','GyroR','GyroZ'};
            plot(Time_Eye,GyroL,'k:',Time_Eye,GyroR,'k--',Time_Eye,GyroZ,'k-',Time_Eye,40*Stim-30,'b')
        end
            hold on
            plot(Time_Eye,LX,'Color',colors.l_x)
            plot(Time_Eye,RX,'Color',colors.r_x)
            plot(Time_Eye,LY,'Color',colors.l_y)
            plot(Time_Eye,RY,'Color',colors.r_y)
            plot(Time_Eye,LZ,'Color',colors.l_z)
            plot(Time_Eye,RZ,'Color',colors.r_z)
            hold off
            legend([leg_text,{'Trigger','LX','RX','LY','RY','LZ','RZ'}])
    elseif ismember(VOG_files{indx(i)},GNO_files)
        data = table2array(readtable(file));
        Time_Eye = (data(:,1) - data(1,1))/10e6;
        GyroZ = data(:,4);
        GyroL = data(:,3);
        GyroR = data(:,2);
        if contains(VOG_files{indx(i)},'LARP')
            RZ = data(:,5);
            RL = data(:,6);
            RR = NaN*data(:,6);
            RY = NaN*data(:,6);
        elseif contains(VOG_files{indx(i)},'RALP')
            RZ = data(:,5);
            RR = data(:,6);
            RL = NaN*data(:,6);
            RY = NaN*data(:,6);
        elseif contains(VOG_files{indx(i)},'Lateral')
            RZ = data(:,5);
            RY = data(:,6);
            RR = NaN*data(:,6);
            RL = NaN*data(:,6);
        end
        leg_text = {'GyroL','GyroR','GyroZ'};
        plot(Time_Eye,GyroL,'k:',Time_Eye,GyroR,'k--',Time_Eye,GyroZ,'k-')
            hold on
            plot(Time_Eye,RY,'Color',colors.r_y)
            plot(Time_Eye,RL,'Color',colors.r_l)
            plot(Time_Eye,RR,'Color',colors.r_r)
            plot(Time_Eye,RZ,'Color',colors.r_z)
            hold off
            legend([leg_text,{'RY','RL','RR','RZ'}])
        set(gca,'YLim',[-300 300])
    end
    xlabel('Time (s)')
    ylabel('Position (deg)')
    title(strrep(strrep(VOG_files{indx(i)},'_',' '),'-',' '))
    %% Find the median eye position during the time specified
    q_another = questdlg('Select another range?','Select another range?','Yes','No','Yes');
    while(strcmp(q_another,'Yes'))
        [xpos,~] = ginput(2);
        rel_inds = Time_Eye>xpos(1)&Time_Eye<xpos(2);
        disp(['Time(s): ',num2str(xpos(1)),'-',num2str(xpos(2))]) 
        disp(['LZ: ',num2str(median(LZ(rel_inds),'omitnan'))])
        disp(['RZ: ',num2str(median(RZ(rel_inds),'omitnan'))])
        disp(['LY: ',num2str(median(LY(rel_inds),'omitnan'))])
        disp(['RY: ',num2str(median(RY(rel_inds),'omitnan'))])
        disp(['LX: ',num2str(median(LX(rel_inds),'omitnan'))])
        disp(['RX: ',num2str(median(RX(rel_inds),'omitnan'))])
        disp(newline)
        q_another = questdlg('Select another range?','Select another range?','Yes','No','Yes');
    end
end