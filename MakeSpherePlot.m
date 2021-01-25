function hg = MakeSpherePlot(CycAvg,inputfig,Function,plotstimaxis,plotelecaxis,normlen,plot_colors,stim_ear)
%% Adpated from PJB's MVIPlotVect_CurrFit on vnelhuman/Figure/Scripts
% inputfig
% the handle to an existing figure if one exists--neccessary for plotting
% results from multiple CycAvg files

% Function (default 1)
% 1 | Plot Mean VOR axis
% 2 | Plot Mean VOR axis and a cone representing the covariance matrix of
%     the saved cycles
% 3 | Plot Mean VOR axis, the axis of each cycle, and a cone representing
%     the covariance matrix of the saved cycles

% plotstimaxis (default 1)
% 0 | do not plot the intended axis of eye movement
% 1 | plot the intended axis of eye movement

% plotelecaxis (default (1)
% 0 | do not plot the axes of eye movements that would be created by pure 
%     excitatory stimulation of the implanted canals
% 1 | plot the axes of eye movements that would be created by pure 
%     excitatory stimulation of the implanted canals

% normlen (default 1)
% 0 | plot the eye response vectors with different lengths depending on
%     the magnitude of eye response (100 deg/s = length of axes)
% 1 | plot all the eye response vectors with a normalized length of 1

% plot_colors (default black ([0,0,0]) left eye and gray (0.5*[1,1,1]) right
% eye)
% Can enter colors for the output axis with the left eye in the first row 
% and the right eye in second row. Input is a 2x3 matrix.

% stim_ear (default left)
% Use this parameter to enter the implant ear if there is no
% CycAvg.info.ear variable. It will accept L/R and Left/Right in lowercase
% as well.
%% Set defaults 
%inputfig
if nargin<2 || (nargin >=2 && isempty(inputfig))
    hg = figure;
else
    hg = inputfig;
    figure(hg)
end
%Function
if nargin<3 || (nargin >=3 && isempty(Function))
    Function = 1;
end
%plotstimaxis
if nargin<4 || (nargin >=4 && isempty(plotstimaxis))
    plotstimaxis = 1;
end
%plotelecaxis
if nargin<5 || (nargin >=5 && isempty(plotelecaxis))
    plotelecaxis = 1;
end
%normlen
if nargin<6 || (nargin >=6 && isempty(normlen))
    normlen = 1;
end
%plot_colors
if nargin<7 || (nargin >=7 && isempty(plot_colors))
    plot_colors = [0,0,0;0.5,0.5,0.5];
end
%stim_ear
if nargin<8 || (nargin >=8 && isempty(stim_ear))
    if ismember('info',fieldnames(CycAvg)) && strcmp(CycAvg.info.ear,'R') %Right
        w = [-1,-1,1];
        view_ang = [-45,15];
    else %Left is the default
        w = [1,1,-1];
        view_ang = [135,-15];
    end
else
    if ismember(stim_ear,{'R','r','right','Right'})
        w = [-1,-1,1];
        view_ang = [-45,15];
    else
        w = [1,1,-1];
        view_ang = [135,-15];
    end
end
%% Plot
ondur = 200/1000; %on duration in seconds 
%Standardize colors
colors.l_x = [255 140 0]/255;
colors.l_y = [128 0 128]/255;
colors.l_z = [1 0 0];
colors.l_l = [0,128,0]/255;
colors.l_r = [0 0 1];
colors.r_x = [238 238 0]/255;
colors.r_y = [138 43 226]/255;
colors.r_z = [255,0,255]/255;
colors.r_l = [0 1 0];
colors.r_r = [64,224,208]/255;

[x,y,z]=sphere();
h=surf(0.5*x,0.5*y,0.5*z);
set(h,'FaceColor','white')
axis off
% Align axis for viewing in 3D
axis vis3d
hold on
%Plot electrode stimulation axis
if(plotelecaxis)
    plot3vect([1*w(1);0;0],'','--',colors.l_l);
    plot3(1*w(1),0,0,'','MarkerSize',8,'LineWidth',2,'Color',colors.l_l)
    plot3vect([0;1*w(2);0],'','--',colors.l_r);
    plot3(0,1*w(2),0,'','MarkerSize',8,'LineWidth',2,'Color',colors.l_r)
    plot3vect([0;0;1*w(3)],'','--',colors.l_z);
    plot3(0,0,1*w(3),'','MarkerSize',8,'LineWidth',2,'Color',colors.l_z)
end
%Plot intended axis (opposite of virtual head velocity)
if(plotstimaxis)
    stimaxis = -CycAvg.info.stim_axis;
    plot3vect(stimaxis','','k--');
    plot3(stimaxis(1),stimaxis(2),stimaxis(3),'','MarkerSize',8,'LineWidth',2,'Color','k')
end

Fs = CycAvg.Fs;
ondur_samp = ondur*Fs;
start1 = 1;
end1 = ondur_samp; 

%All Plot Calculations

%Left eye
all_max_l = zeros(1,3);
all_inds_l = zeros(1,3);
[all_max_l(1),all_inds_l(1)] = max(abs(CycAvg.lz_cycavg(start1:end1)));
[all_max_l(2),all_inds_l(2)] = max(abs(CycAvg.lr_cycavg(start1:end1)));
[all_max_l(3),all_inds_l(3)] = max(abs(CycAvg.ll_cycavg(start1:end1)));
[~,i_m] = max(all_max_l);
Imax_l = all_inds_l(i_m) + start1 -1 ;
VORdat_l = [CycAvg.ll_cycavg(Imax_l) CycAvg.lr_cycavg(Imax_l) CycAvg.lz_cycavg(Imax_l)];
VORdatnorm_l = norm(VORdat_l(end,:),2);
if(normlen)
    VORdatvect_l = VORdat_l(end,:)/VORdatnorm_l(end);
else
    VORdatvect_l = VORdat_l(end,:)/VORdatnorm_l(end)*(VORdatnorm_l(end)/100 + 0.5);
end
VORcycvals_l = [CycAvg.ll_cyc(:,Imax_l),CycAvg.lr_cyc(:,Imax_l),CycAvg.lz_cyc(:,Imax_l)];
dataToUse_l = VORcycvals_l./repmat(sqrt((VORcycvals_l(:,1).^2) + (VORcycvals_l(:,2).^2) + (VORcycvals_l(:,3).^2)),1,3);
%Find covariance matrix
C = cov(dataToUse_l);
% Find the eigenvalue decomposition
[W,D]  = eig(C);
N = 1; % choose how many 'eigenvalues' to plot as your vector length
radii = N*sqrt(diag(D));
[xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));
a = kron(W(:,1),xc-xc); b = kron(W(:,2),yc); c = kron(W(:,3),zc);
data = a+b+c; n = size(data,2);
x = data(1:n,:)+(VORdatvect_l(end,1)); y = data(n+1:2*n,:)+(VORdatvect_l(end,2)); z = data(2*n+1:end,:)+(VORdatvect_l(end,3));
N=80;
zerotemp = zeros(size(x,1),size(x,2));
%x
dx = x/(N-1);
X2 = repmat(dx,1,1,N);
X2(:,:,1) = zerotemp;
X2 = cumsum(X2,3);
%y
dy = y/(N-1);
Y2 = repmat(dy,1,1,N);
Y2(:,:,1) = zerotemp;
Y2 = cumsum(Y2,3);
%z
dz = z/(N-1);
Z2 = repmat(dz,1,1,N);
Z2(:,:,1) = zerotemp;
Z2 = cumsum(Z2,3);        
            
%Right eye
all_max_r = zeros(1,3);
all_inds_r = zeros(1,3);
[all_max_r(1),all_inds_r(1)] = max(abs(CycAvg.rl_cycavg(start1:end1)));
[all_max_r(2),all_inds_r(2)] = max(abs(CycAvg.rz_cycavg(start1:end1)));
[all_max_r(3),all_inds_r(3)] = max(abs(CycAvg.rr_cycavg(start1:end1)));
[~,i_m] = max(all_max_r);
Imax_r = all_inds_r(i_m) + start1 - 1;
VORdat_r = [CycAvg.rl_cycavg(Imax_r) CycAvg.rr_cycavg(Imax_r) CycAvg.rz_cycavg(Imax_r)];
VORdatnorm_r = norm(VORdat_r(end,:),2);
if(normlen)
    VORdatvect_r = VORdat_r(end,:)/VORdatnorm_r(end);
else
    VORdatvect_r = VORdat_r(end,:)/VORdatnorm_r(end)*(VORdatnorm_r(end)/100 + 0.5);
end
VORcycvals_r = [CycAvg.rl_cyc(:,Imax_r),CycAvg.rr_cyc(:,Imax_r),CycAvg.rz_cyc(:,Imax_r)]; 
dataToUse_r = VORcycvals_r./repmat(sqrt((VORcycvals_r(:,1).^2) + (VORcycvals_r(:,2).^2) + (VORcycvals_r(:,3).^2)),1,3);
plot3vect(VORdatvect_r(end,:),'','k:');
%Find covariance matrix
C = cov(dataToUse_r);
% Find the eigenvalue decomposition
[W,D]  = eig(C);
N = 1; % choose how many 'eigenvalues' to plot as your vector length
radii = N*sqrt(diag(D));
[xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));
a = kron(W(:,1),xc-xc); b = kron(W(:,2),yc); c = kron(W(:,3),zc);
data = a+b+c; n = size(data,2);
x = data(1:n,:)+(VORdatvect_r(end,1)); y = data(n+1:2*n,:)+(VORdatvect_r(end,2)); z = data(2*n+1:end,:)+(VORdatvect_r(end,3));
N=80;
zerotemp = zeros(size(x,1),size(x,2));
%x
dx = x/(N-1);
X2r = repmat(dx,1,1,N);
X2r(:,:,1) = zerotemp;
X2r = cumsum(X2r,3);
%y
dy = y/(N-1);
Y2r = repmat(dy,1,1,N);
Y2r(:,:,1) = zerotemp;
Y2r = cumsum(Y2r,3);
%z
dz = z/(N-1);
Z2r = repmat(dz,1,1,N);
Z2r(:,:,1) = zerotemp;
Z2r = cumsum(Z2r,3);    
       
%Plot        
plot3vect(VORdatvect_l(end,:),'','-',plot_colors(1,:));
plot3vect(VORdatvect_r(end,:),'','-',plot_colors(2,:));
if Function == 2
    for ll=1:N
        sc = surf(X2(:,:,ll),Y2(:,:,ll),Z2(:,:,ll));
        set(sc,'FaceColor','none','FaceAlpha',0.1,'EdgeColor',0.5*(1+plot_colors(1,:)),'EdgeAlpha',0.1)
        sc = surf(X2r(:,:,ll),Y2r(:,:,ll),Z2r(:,:,ll));
        set(sc,'FaceColor','none','FaceAlpha',0.1,'EdgeColor',0.5*(1+plot_colors(2,:)),'EdgeAlpha',0.1)
    end
end
if Function == 3
    for i=1:size(VORcycvals_l,1)
        %Left
        temp_vects = VORcycvals_l(i,:);
        norm_vects = sqrt(temp_vects(:,1).^2 + temp_vects(:,2).^2 + temp_vects(:,3).^2);
        if(normlen)
            plot3vect((temp_vects/norm_vects),'','-',plot_colors(1,:),0.5);
        else
            plot3vect((temp_vects/norm_vects)*(norm_vects/100+0.5),'','-',plot_colors(1,:),0.5);
        end
        %Right
        temp_vects = VORcycvals_r(i,:);
        norm_vects = sqrt(temp_vects(:,1).^2 + temp_vects(:,2).^2 + temp_vects(:,3).^2);
        if(normlen)
            plot3vect((temp_vects/norm_vects),'',':',plot_colors(2,:),0.5);
        else
            plot3vect((temp_vects/norm_vects)*(norm_vects/100+0.5),'','-',plot_colors(2,:),0.5);
        end
    end  
end   
view(view_ang)
end