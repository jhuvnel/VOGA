load('VNELcolors.mat','colors')
lrz_xyz = 'xyz';
plot_eyes = 1;
Raw_Path = cd;
all_files = extractfield(dir(Raw_Path),'name',find(~extractfield(dir(Raw_Path),'isdir')));
LDVOG_files = all_files(contains(all_files,'SESSION')&contains(all_files,'.txt')&~contains(all_files,'Notes.txt'));
NKI_files = all_files(contains(all_files,'.dat'));
GNO_files = all_files(contains(all_files,{'Lateral.txt','LARP.txt','RALP.txt'}));
VOG_files = [LDVOG_files;NKI_files;GNO_files];
[indx,tf] = nmlistdlg('PromptString','Select files to plot:',...
    'ListSize',[300 300],...
    'ListString',VOG_files,...
    'SelectionMode','multiple');
figure;
file = [Raw_Path,filesep,VOG_files{indx(1)}];
%%
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
AccelX = data.AccelX;
AccelY = data.AccelY;
AccelZ = data.AccelZ;
%%
%plot(Time_Eye,GyroX,'k:',Time_Eye,GyroY,'k--',Time_Eye,GyroZ,'k-')
%leg_text = {'GyroX','GyroY','GyroZ'};
plot(Time_Eye,AccelX,'g:',Time_Eye,AccelY,'b--',Time_Eye,AccelZ,'k-')
leg_text = strrep({'GyroX','GyroY','GyroZ'},'Gyro','Accel');
legend(leg_text)
%hold on
% plot(Time_Eye,LX,'Color',colors.l_x)
% plot(Time_Eye,RX,'Color',colors.r_x)
% plot(Time_Eye,LY,'Color',colors.l_y)
% plot(Time_Eye,RY,'Color',colors.r_y)
% plot(Time_Eye,LZ,'Color',colors.l_z)
% plot(Time_Eye,RZ,'Color',colors.r_z)
%hold off
%legend([leg_text,{'LX','RX','LY','RY','LZ','RZ'}])