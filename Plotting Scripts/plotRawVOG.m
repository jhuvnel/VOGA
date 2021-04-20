function plotRawVOG(Raw_Path,plot_eyes,lrz_xyz)
    % Standardize Colors
    load('VNELcolors.mat','colors')
    %% Select File
    if nargin < 3
        lrz_xyz = 'xyz';
    end
    if nargin < 2
        plot_eyes = 0;
    end   
    if nargin < 1
        Raw_Path = cd;
    end
    all_files = extractfield(dir(Raw_Path),'name',find(~extractfield(dir(Raw_Path),'isdir')));
    LDVOG_files = all_files(contains(all_files,'SESSION')&contains(all_files,'.txt')&~contains(all_files,'Notes.txt'));
    NKI_files = all_files(contains(all_files,'.dat'));
    GNO_files = all_files(contains(all_files,{'Lateral.txt','LARP.txt','RALP.txt'}));
    VOG_files = [LDVOG_files;NKI_files;GNO_files];   
    [indx,tf] = nmlistdlg('PromptString','Select files to plot:',...
                               'ListSize',[300 300],...
                               'ListString',VOG_files,...
                               'SelectionMode','multiple');
    if tf~=1
        return;
    end
    for i = 1:length(indx)
        %% Make Plot
        file = [Raw_Path,filesep,VOG_files{indx(i)}];
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
        elseif contains(file,'.mat') %NKI preprocessed/uncommon
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
        if any(strcmp(lrz_xyz,{'xyz','XYZ'}))
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
        title(strrep(strrep(VOG_files{indx(i)},'_',' '),'-',' '))  
    end
end