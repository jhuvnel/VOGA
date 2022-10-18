dth1 = 30; %deg
dth2 = 0; %set to 0 for only one spacing
theta = [0:dth1:360,0:dth2:360];
phi = [0:dth1:180,0:dth2:180];
vects = NaN(length(theta)*length(phi),3);
k = 1;
for i = 1:length(theta)
    for j = 1:length(phi)
        vects(k,:) = round([sind(phi(j))*cosd(theta(i)),sind(phi(j))*sind(theta(i)),cosd(phi(j))],3);
        k = k+1;
    end
end
vects = fliplr(unique(vects,'rows'));
vects = vects((length(vects)/2+1):end,:); %Take the positive half
%Remove +- LRZXY (already tested in eeVOR set)
cardinal_ax = [1,   0,  0;...
               -1,  0,  0;...
               0,   1,  0;...
               0,   -1, 0;...
               0,   0,  1;...
               0,   0,  -1;...
               0.707,0.707,0;...
               -0.707,-0.707,0;...
               0.707,-0.707,0;...
               -0.707,0.707,0];
vects(ismember(vects,cardinal_ax,'rows'),:) = [];                   
%Make plot
[x,y,z]=sphere();
h=surf(0.5*x,0.5*y,0.5*z);
set(h,'FaceColor','white')
axis off
% Align axis for viewing in 3D
axis vis3d
hold on
plot3vect([1,0,0],'.','--','g');
plot3(1,0,0,'.','MarkerSize',18,'LineWidth',2,'Color','g')
plot3vect([0,1,0],'','--','b');
plot3(0,1,0,'.','MarkerSize',18,'LineWidth',2,'Color','b')
plot3vect([0,0,1],'','--','r');
plot3(0,0,1,'.','MarkerSize',18,'LineWidth',2,'Color','r')
for i = 1:length(vects)
    plot3vect(vects(i,:),'','--','k');
    plot3vect(-vects(i,:),'','--',.5*[1,1,1]);
end
hold off
%% Write txt files
filenum = 3; %set manually--it will try to evenly distribute
%Find the number of stimuli per file
num_stim = floor(size(vects,1)/filenum)*ones(1,filenum);
%Fix if not an even distribution
left = zeros(1,filenum);
left(1:(size(vects,1)-sum(num_stim)))=1;
num_stim = num_stim+left;
stim_starts = cumsum(num_stim(1:end-1))+1;
ind_num = (0:size(vects,1)-1)';
for i = 1:length(stim_starts)
   ind_num(stim_starts(i):end) = ind_num(stim_starts(i):end)-ind_num(stim_starts(i));
end
labs = {'Stimulus Number','LARP Velocity (Deg/s)','RALP Velocity (Deg/s)','LHRH Velocity (Deg/s)',...
    'Frequency(Hz)','Freq Index','Number of Cycles','Time to Next Stimulus(ms)'};
full_array = [ind_num,round(500*vects),repmat([2,7,20,5000],size(vects,1),1)];
starts = find(ind_num==0);
ends = [starts(2:end)-1;length(ind_num)];
for f = 1:filenum
    T = cell2table([labs;num2cell(full_array(starts(f):ends(f),:))]);
    writetable(T,['File ',num2str(f),' MultiVector Sinusoids.txt'],'Delimiter','\t','WriteVariableNames',false);
end