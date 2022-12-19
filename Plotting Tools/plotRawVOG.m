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
    if isfile(Raw_Path)
        rel_files = {Raw_Path};
        indx = 1;
    else
        VOG_files = extractfield([dir([Raw_Path,filesep,'SESSION*.txt']);dir([Raw_Path,filesep,'*.dat'])...
            dir([Raw_Path,filesep,'*Lateral.txt']);dir([Raw_Path,filesep,'LARP*.txt']);dir([Raw_Path,filesep,'RALP*.txt'])],'name'); 
        VOG_files(contains(VOG_files,'-Notes')) = [];
        [indx,tf] = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',VOG_files,'SelectionMode','multiple');
        if tf~=1
            return;
        end
        rel_files = strcat(Raw_Path,filesep,VOG_files(indx));
    end
    for i = 1:length(rel_files)
        %% Make Plot
        figure;
        file = rel_files{i};
        if contains(file,'.dat')
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
            GyroL = (GyroX - GyroY)/sqrt(2);
            GyroR = (GyroX + GyroY)/sqrt(2);
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
        elseif contains(file,'SESSION')
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
        elseif contains(file,{'Lateral.txt','LARP.txt','RALP.txt'})
            data = table2array(readtable(file));
            Time_Eye = (data(:,1) - data(1,1))/10e6;
            GyroZ = data(:,4);
            GyroL = data(:,3);
            GyroR = data(:,2);
            if contains(rel_files{indx(i)},'LARP')
                RZ = data(:,5);
                RL = data(:,6); 
                RR = NaN*data(:,6);
                RY = NaN*data(:,6);
            elseif contains(rel_files{indx(i)},'RALP')
                RZ = data(:,5);
                RR = data(:,6); 
                RL = NaN*data(:,6);
                RY = NaN*data(:,6);
            elseif contains(rel_files{indx(i)},'Lateral')
                RZ = data(:,5);
                RY = data(:,6);
                RR = NaN*data(:,6);
                RL = NaN*data(:,6);
            end
            leg_text = {'GyroL','GyroR','GyroZ'};
            plot(Time_Eye,GyroL,'k:',Time_Eye,GyroR,'k--',Time_Eye,GyroZ,'k-')
            if plot_eyes
                hold on
                plot(Time_Eye,RY,'Color',colors.r_y)
                plot(Time_Eye,RL,'Color',colors.r_l)
                plot(Time_Eye,RR,'Color',colors.r_r)
                plot(Time_Eye,RZ,'Color',colors.r_z)
                hold off
                legend([leg_text,{'RY','RL','RR','RZ'}])
            else
                legend(leg_text)
            end 
            set(gca,'YLim',[-300 300])
        end      
        xlabel('Time (s)')
        ylabel('Velocity (dps)')
        title(strrep(strrep(rel_files{indx(i)},'_',' '),'-',' '),'interpreter','none') 
    end
end