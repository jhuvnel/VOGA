%Takes in a cycle average file that should have a some burst of stimulation and
%finds the maximum eye velocity and misalignment

% REDUNDANT WITH ParameterizeCycAvg.m--this version not used
function results = ParameterizePulseStim(fname,CycAvg)
%% Extract important information out of the filename
slashes = strfind(fname,filesep);
if isempty(slashes)
    fparts = split(fname,'-');
else
    fparts = split(fname(slashes(end)+1:end),'-');
end
% Subject
if any(contains(fparts,'MVI'))
    subject = strrep(fparts{contains(fparts,'MVI')},'CycAvg_','');
elseif any(contains(fparts,'R'))
    subject = strrep(fparts{contains(fparts,'R')},'CycAvg_','');
else
    subject = '';
end
% Visit
visit = strrep(fparts{contains(fparts,'Visit')},' ','');
% Visit Date
date = datetime(fparts{cellfun(@str2double,fparts)>20160000},'InputFormat','yyyyMMdd');
%Experiment
experiment = strrep(strjoin(fparts(find(contains(fparts,'eeVOR'))+1:end)),'.mat','');
%Ear
if strcmp(CycAvg.info.ear,'L')
    w = 1;
else
    w = -1;
end
%Stim Vector
if contains(experiment,{'LA','RP'})
    stim_vect = w*[-1 0 0];
elseif contains(experiment,{'LP','RA'})
    stim_vect = w*[0 -1 0];
elseif contains(experiment,{'LH','RH'})
    stim_vect = w*[0 0 1];
elseif contains(experiment,'X')
    stim_vect = w*[-0.707 -0.707 0];
elseif contains(experiment,'Y')
    stim_vect = w*[0.707 -0.707 0];
else %MultiVector
    b1 = strfind(experiment,'[');
    b2 = strfind(experiment,']');
    stim_vect = str2num(experiment(b1+1,b2-1));
end
%% Process Inputs
%Make sure stim trace is time points long (1 x nt)
[nc,nt] = size(CycAvg.ll_cyc);
[ab,ac] = size(CycAvg.stim);
if ac == 1 || ac == nc 
    CycAvg.stim = CycAvg.stim';
    [ab,ac] = size(CycAvg.stim);
end
if ab == 1 && ac == nt %already in the right format
    Stim_CycAvg = CycAvg.stim;
elseif ab==1 && ac > nt %One long stim trace
    Stim_CycAvg = mean(reshape(CycAvg.stim',nt,[]),2)';
elseif ab > 1 %One line per stim trace instead of a mean
    Stim_CycAvg = mean(CycAvg.stim,1);
else
    error('Stim trace not in a recognized format.')
end
%% Find maximum velocity (for all directions) and misalignment
%Find the index with the maximum magnitude in the intended direction on the
%average trace
stim_high  = 1:find(diff(Stim_CycAvg)<0);
[max_l,max_l_i] = max(-stim_vect*[CycAvg.ll_cycavg(stim_high);CycAvg.lr_cycavg(stim_high);CycAvg.lz_cycavg(stim_high)]);
[max_r,max_r_i] = max(-stim_vect*[CycAvg.rl_cycavg(stim_high);CycAvg.rr_cycavg(stim_high);CycAvg.rz_cycavg(stim_high)]);
if max_l > max_r
    I = max_l_i;
else
    I = max_r_i;
end
max_vel = [CycAvg.ll_cyc(:,I),CycAvg.lr_cyc(:,I),CycAvg.lz_cyc(:,I),CycAvg.ly_cyc(:,I),CycAvg.lx_cyc(:,I),...
    CycAvg.rl_cyc(:,I),CycAvg.rr_cyc(:,I),CycAvg.rz_cyc(:,I),CycAvg.ry_cyc(:,I),CycAvg.rx_cyc(:,I)];
max_vel_avg = mean(max_vel);
max_vel_std = std(max_vel);
[mis_l_avg,mis_l_std] = calc_misalignment(stim_vect,[CycAvg.ll_cyc(:,I) CycAvg.lr_cyc(:,I) CycAvg.lz_cyc(:,I)]);
[mis_r_avg,mis_r_std] = calc_misalignment(stim_vect,[CycAvg.rl_cyc(:,I) CycAvg.rr_cyc(:,I) CycAvg.rz_cyc(:,I)]);
%% Make Summary Table
dir = {'LARP_','RALP_','LHRH_','Y_','X_'};
labs = [strcat('L_',dir,'MaxVel'),strcat('R_',dir,'MaxVel'),...
    strcat('L_',dir,'MaxVel_sd'),strcat('R_',dir,'MaxVel_sd'),...
    {'L_Align','R_Align','L_Align_sd','R_Align_sd'}]; 
%% Compile into a table
tab1 = cell2table({fname,subject,visit,date,experiment});
tab1.Properties.VariableNames = {'File','Subject','Visit','Date','Experiment'};
tab2 = array2table([nc,max_vel_avg,max_vel_std,mis_l_avg,mis_r_avg,mis_l_std,mis_r_std]);
tab2.Properties.VariableNames = [{'n_cycles'},labs];
results = [tab1,tab2];
end