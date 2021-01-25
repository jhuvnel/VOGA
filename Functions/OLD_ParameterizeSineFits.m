%Takes in a cycle average file that should have a sinusoidal output and
%finds the gain, phase, and misalignment

%COMPLETELY OBSOLETE WITH ParameterizeCycAvg.m
function [results,CycAvg] = ParameterizeSineFits(fname,CycAvg)
%% Extract important information out of the filename
slashes = strfind(fname,filesep);
if isempty(slashes)
    fparts = split(fname,'-');
else
    fparts = split(fname(slashes(end)+1:end),'-');
end
% Subject
if contains(fparts,'MVI')
    subject = strrep(fparts{contains(fparts,'MVI')},'CycAvg_','');
elseif any(~isnan(str2double(strrep(strrep(fparts,'CycAvg_',''),'R','')))&contains(fparts,'R')) %in form R#
    subject = strrep(fparts{~isnan(str2double(strrep(strrep(fparts,'CycAvg_',''),'R','')))&contains(fparts,'R')},'CycAvg_','');
else
    subject = '';
end
% Visit
visit = strrep(fparts{contains(fparts,'Visit')},' ','');
% Visit Date
date = datetime(fparts{cellfun(@str2double,fparts)>20160000},'InputFormat','yyyyMMdd');
% Frequency
freqs = fparts(contains(fparts,'Hz'));
nf = length(freqs);
freq = zeros(1,nf);
for i = 1:nf
    freqs{i} = [strrep(strrep(strrep(freqs{i},'p','.'),'.mat',''),'Hz',''),'Hz'];
    freq(i) = str2double(strrep(freqs{i},'Hz',''));
end
%Condition
if contains(visit,'Visit0')
    condition = 'Pre-Op';
elseif any(contains(fparts,'MotionMod')) || any(contains(fparts,'ModON'))
    condition = 'MotionMod';
elseif any(contains(fparts,'Baseline')) || any(contains(fparts,'Constant'))
    condition = 'ConstantRate';
elseif any(contains(fparts,'Dark')) || any(contains(fparts,'NoStim')) || any(contains(fparts,'NOSTIM'))
    condition = 'NoStim';
else
    condition = '';
end
%Axis
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
ChairAmp = (max(Stim_CycAvg) - min(Stim_CycAvg))/2;
Fs = CycAvg.Fs;
tt = 0:1/Fs:(length(Stim_CycAvg)-1)/Fs;
if length(tt) > 100
    t = tt(round(linspace(1,length(Stim_CycAvg),100)));
else
    t = tt;
end
%% Find gain, phase, and misalignment for each trace
traces = {'lz','rz','ll','rl','lr','rr','lx','rx','ly','ry'};
options = optimset('Display','off');
makefit = @(trace) fminsearchbnd(@(p) sum((sine_fit(t,freq,p)-trace).^2,'omitnan'),repmat([abs(min(trace));abs(max(trace));0],nf,1),repmat([0;0;-180],1,nf),repmat([inf,inf,180],1,nf),options); 
All_Cyc = NaN(nc,3*length(traces),nf);
cyc_fit = NaN(nc,length(tt),length(traces));
for i = 1:nc
    for j = 1:length(traces)
        All_Cyc(i,(3*j-2):3*j,:) = reshape(makefit(spline(tt,CycAvg.([traces{j},'_cyc'])(i,:),t)),1,3,nf);
        cyc_fit(i,:,j) = sine_fit(tt,freq,makefit(spline(tt,CycAvg.([traces{j},'_cyc'])(i,:),t)));
    end
end
%% Make Summary Table
half_cyc = {'L_POS','R_POS','L_NEG','R_NEG'};
labs = [strcat([strcat('LHRH_',half_cyc),strcat('LARP_',half_cyc),strcat('RALP_',half_cyc),strcat('X_',half_cyc),strcat('Y_',half_cyc)],'_Gain'),{'L_Phase','R_Phase'}]; 
Cycles = NaN(nc,22,nf);
val_mean = NaN(nf,22);
val_std = NaN(nf,22);
for i = 1:nf
    %Take the phase from the axis with the most response
    sub_Cyc = All_Cyc(:,:,i);
    [~,ax_ind] = max([mean(mean(sub_Cyc(:,[1,2,4,5]))),mean(mean(sub_Cyc(:,[1,2,4,5]+6))),mean(mean(sub_Cyc(:,[1,2,4,5]+12))),mean(mean(sub_Cyc(:,[1,2,4,5]+18))),mean(mean(sub_Cyc(:,[1,2,4,5]+24)))]);
    Cycles(:,:,i) = [sub_Cyc(:,[1,4,2,5,7,10,8,11,13,16,14,17,19,22,20,23,25,28,26,29])/ChairAmp,sub_Cyc(:,3+(ax_ind-1)*6),sub_Cyc(:,6+(ax_ind-1)*6)];
    val_mean(i,:) = mean(Cycles(:,:,i),1);
    val_std(i,:) = std(Cycles(:,:,i),1);
    if length(val_std)==1
        val_std = zeros(1,size(val_mean,2));
    end
end
%% Misalignment
stim_vect = [0 0 1]; %ideally perfectly horizontal
ang_mean = zeros(1,4);
ang_std = zeros(1,4);
[~,i_lmin] = min(CycAvg.lz_cycavg);
[~,i_lmax] = max(CycAvg.lz_cycavg);
[~,i_rmin] = min(CycAvg.rz_cycavg);
[~,i_rmax] = max(CycAvg.rz_cycavg);
[ang_mean(1), ang_std(1)] = calc_misalignment(stim_vect,[CycAvg.ll_cyc(:,i_lmin),CycAvg.lr_cyc(:,i_lmin),CycAvg.lz_cyc(:,i_lmin)]);
[ang_mean(2), ang_std(2)] = calc_misalignment(stim_vect,[CycAvg.rl_cyc(:,i_rmin),CycAvg.rr_cyc(:,i_rmin),CycAvg.rz_cyc(:,i_rmin)]);
[ang_mean(3), ang_std(3)] = calc_misalignment(-stim_vect,[CycAvg.ll_cyc(:,i_lmax),CycAvg.lr_cyc(:,i_lmax),CycAvg.lz_cyc(:,i_lmax)]);
[ang_mean(4), ang_std(4)] = calc_misalignment(-stim_vect,[CycAvg.rl_cyc(:,i_rmax),CycAvg.rr_cyc(:,i_rmax),CycAvg.rz_cyc(:,i_rmax)]);
%% Compile into a table
tab1 = cell2table([repmat({fname,subject,visit,date,condition},nf,1),freqs]);
tab1.Properties.VariableNames = {'File','Subject','Visit','Date','Condition','Frequency'};
tab2 = array2table([repmat([nc,ChairAmp],nf,1),val_mean,val_std,repmat([ang_mean,ang_std],nf,1)]);
tab2.Properties.VariableNames = [{'n_cycles','ChairAmp'},labs,strcat(labs,'_sd'),strcat(half_cyc,'_Align'),strcat(half_cyc,'_Align_sd')];
results = [tab1,tab2];
CycAvg.parameterized = results;
%% Add fits to the CycAvg Struct
for i = 1:length(traces)
    CycAvg.([traces{i},'_cycavg_fit']) = sine_fit(tt,freq,makefit(spline(tt,CycAvg.([traces{i},'_cycavg']),t)));
    CycAvg.([traces{i},'_cyc_fit']) = reshape(cyc_fit(:,:,i),nc,nt);
end
end