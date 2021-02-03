function [params,fits] = ExponentialFit(CycAvg)
%% Plot Data and decide on analysis
t = CycAvg.t;
stim = CycAvg.stim;
colors = CycAvg.info.colors;

LZ = CycAvg.lz_cyc;
RZ = CycAvg.rz_cyc;
LL = CycAvg.ll_cyc;
RL = CycAvg.rl_cyc;
LR = CycAvg.lr_cyc;
RR = CycAvg.rr_cyc;
LX = CycAvg.lx_cyc;
RX = CycAvg.rx_cyc;
LY = CycAvg.ly_cyc;
RY = CycAvg.ry_cyc;

base_labs = {'GyroZ','LX','RX','LY','RY','LZ','RZ','LLARP','RLARP','LRALP','RRALP'};
h(1) = plot(t,stim,'k');
hold on
h(2) = plot(t,LX,'.','Color',colors.l_x);
h(4) = plot(t,LY,'.','Color',colors.l_y);
h(3) = plot(t,RX,'.','Color',colors.r_x);
h(5) = plot(t,RY,'.','Color',colors.r_y);
h(8) = plot(t,LL,'.','Color',colors.l_l);
h(10) = plot(t,LR,'.','Color',colors.l_r);
h(9) = plot(t,RL,'.','Color',colors.r_l);
h(11) = plot(t,RR,'.','Color',colors.r_r);
h(6) = plot(t,LZ,'.','Color',colors.l_z);
h(7) = plot(t,RZ,'.','Color',colors.r_z);
hold off
legend(h,base_labs)

indx = nmlistdlg('PromptString','Select which traces should be fit to an exponential:',...
                           'ListSize',[300 300],...
                           'ListString',{'X','Y','LHRH','LARP','RALP'});
uiwait(msgbox('Select the start and the end points for the exponential fit'))
[x,~] = ginput(2);
[~,i1] = min(abs(t - x(1)));
[~,i2] = min(abs(t - x(2)));

t_sub = t(i1:i2);
nonan_l = find(~isnan(LZ)==1);
is_l = nonan_l(find(nonan_l>=i1,1,'first'));
nonan_r = find(~isnan(RZ)==1);
is_r = nonan_r(find(nonan_r>=i1,1,'first'));
%% Combination of both eyes
exp_fit = @(tt,p) p(1)*exp(-(tt-t_sub(1))/p(2)) + p(3);

if ismember(indx,1) %X
    LSCF = @(p) sum((exp_fit(t_sub,p) - LX(i1:i2)').^2 + (exp_fit(t_sub,p) - RX(i1:i2)').^2,'omitnan');  
    X_params = fminsearch(LSCF,[mean([LX(is_l);RX(is_r)]);10;0]);
    X_fit = exp_fit(t,X_params); 
else
    X_fit = mean((LX + RX)/2,'omitnan')*ones(1,length(t));
    X_params = [0;0;mean((LX + RX)/2,'omitnan')];
end
X_fit(1:i1-1) = NaN;
if ismember(indx,2) %Y
    LSCF = @(p) sum((exp_fit(t_sub,p) - LY(i1:i2)').^2 + (exp_fit(t_sub,p) - RY(i1:i2)').^2,'omitnan');  
    Y_params = fminsearch(LSCF,[mean([LY(is_l);RY(is_r)]);10;0]);
    Y_fit = exp_fit(t,Y_params);
else
    Y_fit = mean((LY + RY)/2,'omitnan')*ones(1,length(t));
    Y_params = [0;0;mean((LY + RY)/2,'omitnan')];
end
Y_fit(1:i1-1) = NaN;
if ismember(indx,3) %LHRH
    LSCF = @(p) sum((exp_fit(t_sub,p) - LZ(i1:i2)').^2 + (exp_fit(t_sub,p) - RZ(i1:i2)').^2,'omitnan');  
    LHRH_params = fminsearch(LSCF,[mean([LZ(is_l);RZ(is_r)]);10;0]);
    LHRH_fit = exp_fit(t,LHRH_params);
else
    LHRH_fit = mean((LZ + RZ)/2,'omitnan')*ones(1,length(t));
    LHRH_params = [0;0;mean((LZ + RZ)/2,'omitnan')];
end
LHRH_fit(1:i1-1) = NaN;
if ismember(indx,4) %LARP
    LSCF = @(p) sum((exp_fit(t_sub,p) - LL(i1:i2)').^2 + (exp_fit(t_sub,p) - RL(i1:i2)').^2,'omitnan');  
    LARP_params = fminsearch(LSCF,[mean([LL(is_l);RL(is_r)]);10;0]);
    LARP_fit = exp_fit(t,LARP_params);
else
    LARP_fit = mean((LL + RL)/2,'omitnan')*ones(1,length(t));
    LARP_params = [0;0;mean((LL + RL)/2,'omitnan')];
end
LARP_fit(1:i1-1) = NaN;
if ismember(indx,5) %RALP
    LSCF = @(p) sum((exp_fit(t_sub,p) - LR(i1:i2)').^2 + (exp_fit(t_sub,p) - RR(i1:i2)').^2,'omitnan');  
    RALP_params = fminsearch(LSCF,[mean([LR(is_l);RR(is_r)]);10;0]);
    RALP_fit = exp_fit(t,RALP_params);
else
    RALP_fit = mean((LR + RR)/2,'omitnan')*ones(1,length(t));
    RALP_params = [0;0;mean((LR + RR)/2,'omitnan')];
end
RALP_fit(1:i1-1) = NaN;

h2(1) = plot(t,stim,'k');
hold on
h2(2) = plot(t,LX,'.','Color',colors.l_x);
h2(4) = plot(t,LY,'.','Color',colors.l_y);
h2(3) = plot(t,RX,'.','Color',colors.r_x);
h2(5) = plot(t,RY,'.','Color',colors.r_y);
h2(8) = plot(t,LL,'.','Color',colors.l_l);
h2(10) = plot(t,LR,'.','Color',colors.l_r);
h2(9) = plot(t,RL,'.','Color',colors.r_l);
h2(11) = plot(t,RR,'.','Color',colors.r_r);
h2(6) = plot(t,LZ,'.','Color',colors.l_z);
h2(7) = plot(t,RZ,'.','Color',colors.r_z);
h2(12) = plot(t,X_fit,'Color',colors.l_x,'LineWidth',3);
h2(13) = plot(t,Y_fit,'Color',colors.l_y,'LineWidth',3);
h2(15) = plot(t,LARP_fit,'Color',colors.l_l,'LineWidth',3);
h2(16) = plot(t,RALP_fit,'Color',colors.l_r,'LineWidth',3);
h2(14) = plot(t,LHRH_fit,'Color',colors.l_z,'LineWidth',3);
hold off
legend(h2,[base_labs,{'X Fit','Y Fit','Z Fit','LARP Fit','RALP Fit'}])

params = array2table([X_params,Y_params,LHRH_params,LARP_params,RALP_params]);
params.Properties.VariableNames = {'X','Y','LHRH','LARP','RALP'};
params.Properties.RowNames = {'Magnitude','Time_Constant','Offset'};

fits = array2table([X_fit',Y_fit',LHRH_fit',LARP_fit',RALP_fit']);
fits.Properties.VariableNames = {'X','Y','LHRH','LARP','RALP'};
end