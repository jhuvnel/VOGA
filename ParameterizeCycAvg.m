%% ParameterizeCycAvg.m 
%This function takes in a CycAvg file (either PJB with added info or AA)
%and parameterizes accordingly. It adds the results table and any
%appropriate fits to the CycAvg struct.
%Type tells you what kind of output to expect:
%1. Sine Fit (Half-Cycle Magnitude, Gain, Phase, Misalignment)
%2. Exponential Fit (Tau, Magnititude, Misalignment)
%3. Impulse (Magnitude, Gain (Position in High, Accel in Low), Latency)
%4. Magnitude (Max Magnitude, Misalignment)

function [CycAvg,type] = ParameterizeCycAvg(CycAvg)
%% Extract the relevant descriptive parameters
% Subject, Visit, Date, 
% Experiment (RotaryChair/eeVOR/aHIT/Manual), 
% Type(Sine/Exponential/Impulse/PulseTrain),
% Condition (MotionMod/ConstantRate/NoStim/NoStimLight/eeVOR/Autoscan),
% axis_name (X/Y/LHRH/LARP/RALP/''),
% stim_vect ([L,R,Z] coordinates)
% If applicable:
% Frequency (Hz)
% Amplitude (dps)
% Electrode (e.g. LPE3)
% Pulses (pps)
% PhaseDur (us)
% CurrAmp (uA)
fname = strrep(strrep(strrep(CycAvg.name,'CycAvg_',''),'.mat',''),'_CycAvg','');
if contains(fname,'_Updated')
    und = strfind(fname,'_');
    fname = fname(1:und(1)-1);
end
fparts = split(strrep(fname,'NA',''),'-');
fparts(cellfun(@isempty,fparts)) = [];
% Subject
if contains(fname,'MVI')
    subject = fparts{contains(fparts,'MVI')};
    fparts(contains(fparts,'MVI'))= [];
elseif any(~isnan(str2double(strrep(fparts,'R','')))&contains(fparts,'R')) %in form R#
    subject = fparts{~isnan(str2double(strrep(fparts,'R','')))&contains(fparts,'R')};
    fparts(~isnan(str2double(strrep(fparts,'R','')))&contains(fparts,'R')) = [];
else
    subject = fparts{1};
    fparts(1)= [];
end
% Visit
if any(contains(fparts,'Visit'))
    visit = strrep(fparts{contains(fparts,'Visit')},' ','');
    fparts(contains(fparts,'Visit')) = [];
else
    visit = 'NA';
end
% Date
date_ind = find(cellfun(@str2double,fparts)>20160000);
if ~isnan(str2double(fparts{date_ind+1})) %Time
    date = datetime([fparts{date_ind},'-',fparts{date_ind+1}],'InputFormat','yyyyMMdd-HHmmss');
    fparts(date_ind:date_ind+1) = [];
else %no time
    date = datetime(fparts{date_ind},'InputFormat','yyyyMMdd');
    fparts(date_ind) = [];
end
% Experiment
known_exps = {'RotaryChair','aHIT','eeVOR','Manual'};
fname = strrep(strrep(fname,'CurrentFitting','eeVOR'),'ElectricalOnly','eeVOR'); %For old autoscan files
if contains(fname,known_exps)
    experiment = known_exps(cellfun(@(x) contains(fname,x),known_exps));
    fparts(contains(fparts,experiment)) = [];
else
    disp('Cannot parameterize and assign experiment type to:')
    disp(fname)
    return;
end   
%Stim Vector
if contains(fname,{'LARP','RP'})
    stim_vect = [1 0 0];
    axis_name = 'LARP';
elseif contains(fname,'LA')
    stim_vect = [-1 0 0];  
    axis_name = 'LARP';
elseif contains(fname,{'RALP','RA'})
    stim_vect = [0 1 0];
    axis_name = 'RALP';
elseif contains(fname,'LP')
    stim_vect = [0 -1 0];
    axis_name = 'RALP';
elseif contains(fname,{'LHRH','LH','RotaryChair'})
    stim_vect = [0 0 1];
    axis_name = 'LHRH';
elseif contains(fname,'RH')
    stim_vect = [0 0 -1];
    axis_name = 'LHRH';
elseif contains(fname,'X')
    stim_vect = [0.707 0.707 0];
    axis_name = 'X';
elseif contains(fname,'Y')
    stim_vect = [-0.707 0.707 0];
    axis_name = 'Y';
elseif contains(fname,'[')&&contains(fname,']') %MultiVector, SineMultivector
    b1 = strfind(fname,'[');
    b2 = strfind(fname,']');
    stim_vect = str2num(strrep(fname(b1+1:b2-1),' ','-'));
    axis_name = '';
else %Including Activation
    stim_vect = [0,0,0];
    axis_name = '';
end
% Electrode #
if any(contains(fparts,{'LAE','LHE','LPE','RAE','RHE','RPE'}))
    electrode = fparts{contains(fparts,{'LAE','LHE','LPE','RAE','RHE','RPE'})};
else
    electrode = '';
end
fparts(contains(fparts,{'LA','LH','LP','RA','RH','RP','['})) = [];
%Amplitude
amp_cond = ~isnan(str2double(strrep(strrep(fparts,'dps',''),'n','-')))&contains(fparts,'dps');
if any(amp_cond)
    amp = str2double(strrep(strrep(fparts{amp_cond},'n','-'),'dps',''));
    fparts(amp_cond) = [];
else
    amp = NaN;
end
% Pulses per second
if any(~isnan(str2double(strrep(fparts,'pps','')))&contains(fparts,'pps'))
    pps = str2double(strrep(fparts{~isnan(str2double(strrep(fparts,'pps','')))&contains(fparts,'pps')},'pps',''));
    fparts(~isnan(str2double(strrep(fparts,'pps','')))&contains(fparts,'pps')) = [];
else
    pps = NaN;
end
%PhaseDur
if any(~isnan(str2double(strrep(strrep(fparts,'us',''),'uS','')))&contains(fparts,{'us','uS'}))
    phase_dur = str2double(strrep(strrep(fparts{~isnan(str2double(strrep(strrep(fparts,'us',''),'uS','')))&contains(fparts,{'us','uS'})},'us',''),'uS',''));
    fparts(~isnan(str2double(strrep(strrep(fparts,'us',''),'uS','')))&contains(fparts,{'us','uS'})) = [];
else
    phase_dur = NaN;
end
%Current Amplitude
if any(~isnan(str2double(strrep(strrep(fparts,'ua',''),'uA','')))&contains(fparts,{'ua','uA'}))
    curr_amp = str2double(strrep(strrep(fparts{~isnan(str2double(strrep(strrep(fparts,'ua',''),'uA','')))&contains(fparts,{'ua','uA'})},'ua',''),'uA',''));
    fparts(~isnan(str2double(strrep(strrep(fparts,'ua',''),'uA','')))&contains(fparts,{'ua','uA'})) = [];
else
    curr_amp = NaN;
end
%% Type/Condition
Types = {'Sine','Exponential','Impulse','PulseTrain'};
if contains(fname,{'Sine','Sinusoid','Sin'}) %Sine fit
    type = 1;
    %Frequency
    freqs = strrep(strrep(fparts(contains(fparts,'Hz')),'p','.'),'.mat','');
    if isempty(freqs) %Nothing contains Hz so look for # after sinusoid (don't worry about multiple frequencies)
        freqs = fparts(find(contains(fparts,{'Sin','Sine','Sinusoid'}))+1);
        fparts(find(contains(fparts,{'Sin','Sine','Sinusoid'}))+1) = [];
    end
    freqs = str2double(strrep(freqs,'Hz','')); 
    fparts(contains(fparts,'Hz')) = [];
    fparts(contains(fparts,{'Sine','Sinusoid','Sin'})) = [];
elseif contains(fname,{'VelStep','Activation'}) %Exponential fit
    type = 2;
    fparts(contains(fparts,{'VelStep','Activation'})) = [];
    %Frequency
    freqs = NaN;
elseif contains(fname,{'Impulse'}) %Impulse w/ Head Motion
    type = 3;
    fparts(contains(fparts,{'Impulse'})) = [];
    %Frequency
    freqs = NaN;   
elseif contains(fname,'eeVOR') %Pulse Stim, Autoscan
    type = 4;
    fparts(contains(fparts,{'eeVOR'})) = [];
    %Frequency
    freqs = NaN;
else 
    disp('Cannot parameterize and assign experiment type to:')
    disp(fname)
    return;
end
if ~isfield(CycAvg,'info')
   goggle = 'LDVOG';
elseif ~isfield(CycAvg.info,'goggle_ver')
   goggle = 'LDVOG';
else
    if isnumeric(CycAvg.info.goggle_ver)
        goggle = ['LDVOG',num2str(CycAvg.info.goggle_ver)];
    else
        goggle = CycAvg.info.goggle_ver;
    end
end
fparts(contains(fparts,goggle)) = [];
Type = Types{type};
nf = length(freqs);
%Condition
if contains(fname,'eeVOR')&&type==1
    condition = 'eeVOR';
elseif contains(fname,{'Baseline','Constant'})
    condition = 'ConstantRate';
elseif contains(fname,{'Motion','Mod'}) %may have numbers or words after this
    if contains(fname,{'Mod-,','on-','ON-'})
        condition = 'MotionMod';
    else
        condition = fparts{contains(fparts,'Mod')};        
    end
elseif contains(fname,{'Light','LIGHT'})
    condition = 'Light';
elseif contains(fname,{'NoStim','Dark','DARK'})
    condition = 'NoStim';
else
    condition = strrep(strjoin(fparts,' '),'_','');
    if isempty(condition) %Likely a normal subject
        condition = 'NoStim';
    end
end
if ~isfield(CycAvg,'lx_cyc')&&isfield(CycAvg,'CycAvg.ll_cyc')
    CycAvg.lx_cyc = (CycAvg.ll_cyc+CycAvg.lr_cyc)/sqrt(2);
    CycAvg.lx_cycavg = (CycAvg.ll_cycavg+CycAvg.lr_cycavg)/sqrt(2);
    CycAvg.rx_cyc = (CycAvg.rl_cyc+CycAvg.rr_cyc)/sqrt(2);
    CycAvg.rx_cycavg = (CycAvg.rl_cycavg+CycAvg.rr_cycavg)/sqrt(2);    
    CycAvg.ly_cyc = (-CycAvg.ll_cyc+CycAvg.lr_cyc)/sqrt(2);
    CycAvg.ly_cycavg = (-CycAvg.ll_cycavg+CycAvg.lr_cycavg)/sqrt(2);
    CycAvg.ry_cyc = (-CycAvg.rl_cyc+CycAvg.rr_cyc)/sqrt(2);
    CycAvg.ry_cycavg = (-CycAvg.rl_cycavg+CycAvg.rr_cycavg)/sqrt(2);
end
traces = {'lz','rz','ll','rl','lr','rr','lx','rx','ly','ry'};
options = optimset('Display','off'); %suppress output from fminsearch
%initialize vectors for the table
MaxVel = NaN(nf,length(traces)*4);
Gain = NaN(nf,length(traces)*4);
Tau = NaN(nf,length(traces)*4);
RMSE = NaN(nf,length(traces)*2);
Latency = NaN(nf,length(traces)*2);
Phase = NaN(nf,4);
Align = NaN(nf,8);
Disc = NaN(nf,4);
%% Create a table for each type
switch type 
    case 1 
        %% Sine 
        %Make sure stim trace is time points long (1 x nt)
        [nc1,nt] = size(CycAvg.ll_cyc);
        nc2 = size(CycAvg.rl_cyc,1);
        if nc1 > nc2 %more L than R
            CycAvg.rl_cyc = [CycAvg.rl_cyc;NaN(nc1-nc2,nt)];
            CycAvg.rr_cyc = [CycAvg.rr_cyc;NaN(nc1-nc2,nt)];
            CycAvg.rz_cyc = [CycAvg.rz_cyc;NaN(nc1-nc2,nt)];
            CycAvg.rx_cyc = [CycAvg.rx_cyc;NaN(nc1-nc2,nt)];
            CycAvg.ry_cyc = [CycAvg.ry_cyc;NaN(nc1-nc2,nt)];
        elseif nc1 < nc2 % more R than L
            CycAvg.ll_cyc = [CycAvg.ll_cyc;NaN(nc2-nc1,nt)];
            CycAvg.lr_cyc = [CycAvg.lr_cyc;NaN(nc2-nc1,nt)];
            CycAvg.lz_cyc = [CycAvg.lz_cyc;NaN(nc2-nc1,nt)];
            CycAvg.lx_cyc = [CycAvg.lx_cyc;NaN(nc2-nc1,nt)];
            CycAvg.ly_cyc = [CycAvg.ly_cyc;NaN(nc2-nc1,nt)];
        end
        nc = nc1;    
        [ab,ac] = size(CycAvg.stim);
        if ac == 1 || ac == nc1 || ac==nc2 
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
        % Find gain, phase, and misalignment for each trace
        makefit = @(trace) fminsearchbnd(@(p) sum((sine_fit(t,freqs,p)-trace).^2,'omitnan'),repmat([abs(min(trace));abs(max(trace));0],nf,1),repmat([0;0;-180],1,nf),repmat([inf,inf,180],1,nf),options); 
        All_Cyc = NaN(nc,3*length(traces),nf);
        cyc_fit = NaN(nc,length(tt),length(traces));
        for i = 1:nc
            for j = 1:length(traces)
                if ~all(isnan(CycAvg.([traces{j},'_cyc'])(i,:)))
                    All_Cyc(i,(3*j-2):3*j,:) = reshape(makefit(spline(tt,CycAvg.([traces{j},'_cyc'])(i,:),t)),1,3,nf);
                    cyc_fit(i,:,j) = sine_fit(tt,freqs,makefit(spline(tt,CycAvg.([traces{j},'_cyc'])(i,:),t)));
                end
            end
        end
        for i = 1:nf
            %Take the phase from the axis with the most response
            sub_Cyc = All_Cyc(:,:,i);
            [~,ax_ind] = max([mean(mean(sub_Cyc(:,[1,2,4,5]))),mean(mean(sub_Cyc(:,[1,2,4,5]+6))),mean(mean(sub_Cyc(:,[1,2,4,5]+12))),mean(mean(sub_Cyc(:,[1,2,4,5]+18))),mean(mean(sub_Cyc(:,[1,2,4,5]+24)))]);
            phases = [sub_Cyc(:,3+(ax_ind-1)*6),sub_Cyc(:,6+(ax_ind-1)*6)];
            %MaxVel
            MaxVel(i,1:4:end) = mean(sub_Cyc(:,1:3:end),1);
            MaxVel(i,2:4:end) = std(sub_Cyc(:,1:3:end),1);
            MaxVel(i,3:4:end) = mean(sub_Cyc(:,2:3:end),1);
            MaxVel(i,4:4:end) = std(sub_Cyc(:,2:3:end),1);
            %Phase
            Phase(i,1:2:end) = mean(phases);
            Phase(i,2:2:end) = std(phases);
            if nc==1 %If only 1 cycle, set sd to 0
                MaxVel(i,2:4:end) = 0;
                MaxVel(i,4:4:end) = 0;
                Phase(i,2:2:end) = 0;
            end
        end
        %Gain
        Gain = MaxVel/ChairAmp;
        %RMSE
        for i = 1:length(traces)
            RMSE(:,[i,i+length(traces)]) = sum((CycAvg.([traces{i},'_cycavg'])-sine_fit(tt,freqs,makefit(spline(tt,CycAvg.([traces{i},'_cycavg']),t)))).^2)/length(tt);
        end
        % Misalignment
        A = NaN(1,8);
        [~,i_lmax] = max(stim_vect*[CycAvg.ll_cycavg;CycAvg.lr_cycavg;CycAvg.lz_cycavg]);
        [~,i_rmax] = max(stim_vect*[CycAvg.rl_cycavg;CycAvg.rr_cycavg;CycAvg.rz_cycavg]);
        [~,i_lmin] = min(stim_vect*[CycAvg.ll_cycavg;CycAvg.lr_cycavg;CycAvg.lz_cycavg]);
        [~,i_rmin] = min(stim_vect*[CycAvg.rl_cycavg;CycAvg.rr_cycavg;CycAvg.rz_cycavg]);
        [A(1),A(2)] = calc_misalignment(-stim_vect,[CycAvg.ll_cyc(:,i_lmin),CycAvg.lr_cyc(:,i_lmin),CycAvg.lz_cyc(:,i_lmin)]);
        [A(3),A(4)] = calc_misalignment(stim_vect,[CycAvg.ll_cyc(:,i_lmax),CycAvg.lr_cyc(:,i_lmax),CycAvg.lz_cyc(:,i_lmax)]);
        [A(5),A(6)] = calc_misalignment(-stim_vect,[CycAvg.rl_cyc(:,i_rmin),CycAvg.rr_cyc(:,i_rmin),CycAvg.rz_cyc(:,i_rmin)]);
        [A(7),A(8)] = calc_misalignment(stim_vect,[CycAvg.rl_cyc(:,i_rmax),CycAvg.rr_cyc(:,i_rmax),CycAvg.rz_cyc(:,i_rmax)]);
        Align = repmat(A,nf,1);
        %Disconjugacy
        B = NaN(1,4);        
        [B(1),B(2)] = calc_misalignment([CycAvg.ll_cyc(:,i_lmin),CycAvg.lr_cyc(:,i_lmin),CycAvg.lz_cyc(:,i_lmin)],[CycAvg.rl_cyc(:,i_lmin),CycAvg.rr_cyc(:,i_lmin),CycAvg.rz_cyc(:,i_lmin)]);
        [B(3),B(4)] = calc_misalignment([CycAvg.ll_cyc(:,i_lmax),CycAvg.lr_cyc(:,i_lmax),CycAvg.lz_cyc(:,i_lmax)],[CycAvg.rl_cyc(:,i_lmax),CycAvg.rr_cyc(:,i_lmax),CycAvg.rz_cyc(:,i_lmax)]);
        Disc = repmat(B,nf,1);
        % Add fits to the CycAvg Struct
        for i = 1:length(traces)
            CycAvg.([traces{i},'_cycavg_fit']) = sine_fit(tt,freqs,makefit(spline(tt,CycAvg.([traces{i},'_cycavg']),t)));
            CycAvg.([traces{i},'_cyc_fit']) = reshape(cyc_fit(:,:,i),nc,nt);
        end
        cycle_params = [];
    case 2 
        %% Exponential
        nc = 1;
        t = reshape(CycAvg.t,1,[]);
        %Make sure stim trace is time points long (1 x nt)
        if max(abs(CycAvg.stim)) > 1 %Triggered by motion. Should be one high/low and then a period of 0
            stim = reshape(CycAvg.stim,1,[]);
            max_vel = amp;
            Stim = NaN(1,length(t));
            Stim(abs(stim)>0.9*abs(max_vel)) = 1;
            after_motion = zeros(1,length(Stim));
            after_motion(1,find(Stim==1,1,'last'):end) = 1;
            Stim(find(abs(stim)<1&after_motion,1,'first'):end) = 0; 
        else
            Stim = reshape(abs(CycAvg.stim),1,[]);
        end
        if length(t) > 1000
            sub_i = false(1,length(t));
            sub_i(floor(linspace(1,length(t),1000))) = true;
        else
            sub_i = true(1,length(t));
        end
        tt = t(sub_i);
        high_i = t(find(Stim==1,1,'first'));
        low_i = t(find(Stim==0,1,'first'));
        % Find exponential fits
        %2nd order exponential fit same syntax as matlab's built-in 
        makefit = @(tt,p) p(:,1).*exp(p(:,2).*tt)+p(:,3).*exp(p(:,4).*tt);         
        all_fits = NaN(length(traces),length(t));
        params = NaN(2*length(traces),4);
        dconfint = NaN(2*length(traces),4);
        %1order exponential fit with decay to 0
        exp_fit1 = @(tt,p) p(1).*exp(-tt./p(2));
        makefit2 = @(trace,T) fminsearch(@(p) sum((exp_fit1(T,p)-trace).^2,'omitnan')+sum(p.^2),[median(trace(1));1],options); 
        all_fits2 = NaN(length(traces),length(t));
        params2 = NaN(2*length(traces),2);
        %1order exponential fit with decay not forced to 0. 
        exp_fit2 = @(tt,p) p(1).*exp(-tt./p(2)) + p(3);
        makefit3 = @(trace,T) fminsearch(@(p) sum((exp_fit2(T,p)-trace).^2,'omitnan')+sum(p.^2),[median(trace);1;0],options); 
        all_fits3 = NaN(length(traces),length(t));
        params3 = NaN(2*length(traces),3);
        for i = 1:length(traces)               
            long_dat = CycAvg.([traces{i},'_cyc']);
            dat = reshape(long_dat(sub_i),[],1);
            t_high = tt(Stim(sub_i)'==1&~isnan(dat))'-high_i;
            t_low = tt(Stim(sub_i)'==0&~isnan(dat))'-low_i;
            dat_high = dat(Stim(sub_i)'==1&~isnan(dat));
            dat_low = dat(Stim(sub_i)'==0&~isnan(dat));
            %Stim is high
            [fitobj,gof] = fit(t_high,dat_high,'exp2');
            RMSE(1,i) = gof.rmse;
            params(i,:) = coeffvalues(fitobj);
            dconfint(i,:) = mean(abs(confint(fitobj)-mean(confint(fitobj))));
            params2(i,:) = makefit2(dat_high,t_high);
            params3(i,:) = makefit3(dat_high,t_high);
            all_fits(i,Stim==1) = makefit(t(Stim==1)-high_i,params(i,:));
            all_fits2(i,Stim==1) = exp_fit1(t(Stim==1)-high_i,params2(i,:));
            all_fits3(i,Stim==1) = exp_fit2(t(Stim==1)-high_i,params3(i,:));
            %Stim is low
            if length(t_low)>4
                [fitobj,gof] = fit(t_low,dat_low,'exp2');
                RMSE(1,length(traces)+i) = gof.rmse;
                params(length(traces)+i,:) = coeffvalues(fitobj);
                dconfint(length(traces)+i,:) = mean(abs(confint(fitobj)-mean(confint(fitobj))));
                params2(length(traces)+i,:) = makefit2(dat_low,t_low);
                params3(length(traces)+i,:) = makefit3(dat_low,t_low);                
                all_fits(i,Stim==0) = makefit(t(Stim==0)-low_i,params(length(traces)+i,:));
                all_fits2(i,Stim==0) = exp_fit1(t(Stim==0)-low_i,params2(length(traces)+i,:));
                all_fits3(i,Stim==0) = exp_fit2(t(Stim==0)-low_i,params3(length(traces)+i,:));
            end
        end    
        tau_all = -1./params(:,[2,4]);
        tau_all(tau_all < 0) = NaN;  
        [tau,tau_i] = min(tau_all(:,1:2),[],2);
        tau_conf = dconfint(:,2);
        tau_conf(tau_i==2) = dconfint(tau_i==2,4);
        tau(isnan(tau_conf)) = NaN;
        tau_conf(isnan(tau)) = NaN;
        Tau = reshape([tau';tau_conf'],1,[]);
        max_vel = sum(params(:,[1,3]),2);
        max_vel_conf = dconfint(:,1)+dconfint(:,3);
        max_vel(isnan(max_vel_conf)) = NaN;
        max_vel_conf(isnan(max_vel)) = NaN;
        MaxVel = reshape([max_vel';max_vel_conf'],1,[]);     
        % Misalignment
        A = NaN(1,8);
        all_high_L = all_fits([3,5,1],Stim==1);
        all_low_L = all_fits([3,5,1],Stim==0);
        all_high_R = all_fits([4,6,2],Stim==1);
        all_low_R = all_fits([4,6,2],Stim==0);
        [A(1),A(2)] = calc_misalignment(-stim_vect,all_high_L');
        [A(3),A(4)] = calc_misalignment(stim_vect,all_low_L');
        [A(5),A(6)] = calc_misalignment(-stim_vect,all_high_R');
        [A(7),A(8)] = calc_misalignment(stim_vect,all_low_R');
        Align = repmat(A,nf,1);
        %Disconjugacy
        B = NaN(1,4);        
        [B(1),B(2)] = calc_misalignment(all_high_L',all_high_R');
        [B(3),B(4)] = calc_misalignment(all_low_L',all_low_R');
        Disc = repmat(B,nf,1);
        % Add fits to the CycAvg Struct
        for i=1:length(traces)
            CycAvg.([traces{i},'_cycavg_fit']) = all_fits(i,:);
            CycAvg.([traces{i},'_cyc_fit']) = all_fits(i,:);
        end       
        cycle_params.exp_ord2_high_params = array2table(params(1:length(traces),:)','VariableNames',traces);
        cycle_params.exp_ord2_low_params = array2table(params(length(traces)+1:end,:)','VariableNames',traces);
        cycle_params.exp_ord2_fits = all_fits;
        cycle_params.exp_ord1_high_params = array2table(params2(1:length(traces),:)','VariableNames',traces);
        cycle_params.exp_ord1_low_params = array2table(params2(length(traces)+1:end,:)','VariableNames',traces);
        cycle_params.exp_ord1_fits = all_fits2;
        cycle_params.exp_ord1_const_high_params = array2table(params3(1:length(traces),:)','VariableNames',traces);
        cycle_params.exp_ord1_const_low_params = array2table(params3(length(traces)+1:end,:)','VariableNames',traces);
        cycle_params.exp_ord1_const_fits = all_fits3;        
    case 3 
        %% Impulse
        % Process Inputs
        %Make sure stim trace is time points long (1 x nt)
        if isfield(CycAvg,'rz_cyc')
            [nc,nt] = size(CycAvg.rz_cyc);
        else
            [nc,nt] = size(CycAvg.lz_cyc);
        end
        [ab,ac] = size(CycAvg.stim);
        if ab == nt && ac == nc 
            CycAvg.stim = CycAvg.stim';
            [ab,ac] = size(CycAvg.stim);
        end
        if ab == nc && ac == nt %already in the right format
            Stim_All = CycAvg.stim;
        else
            error('Stim trace not in a recognized format.')
        end
        t = CycAvg.t;
        t_upsamp = t(1):0.0001:t(end);
        max_vel = NaN(nc,length(traces));
        gain_AUC = NaN(nc,length(traces));
        lat = NaN(nc,length(traces)); 
        gain_Ga = NaN(nc,length(traces));
        avg_max_vel = NaN(1,length(traces));
        avg_gain_AUC = NaN(1,length(traces));
        avg_lat = NaN(1,length(traces)); 
        avg_gain_Ga = NaN(1,length(traces));
        head_width = NaN(nc,1);
        head_maxvel = NaN(nc,1);
        head_accel = NaN(nc,1);
        saccade_timing = cell(nc,2);
        time_to_target = NaN(nc,2);
        %Where to look for saccades
        if contains(CycAvg.name,{'LH','RH'})
            rel_tr = {'lz','rz'};
        elseif contains(CycAvg.name,{'LA','RP'})
            rel_tr = {'ll','rl'};
        elseif contains(CycAvg.name,{'LP','RA'})
            rel_tr = {'lr','rr'};
        else
            rel_tr = [];
        end
        if -min(CycAvg.stim_cycavg)>max(CycAvg.stim_cycavg)    
            neg_flag = -1;
        else
            neg_flag = 1;
        end
        for tr = 1:length(traces)
            if isfield(CycAvg,[traces{tr},'_cyc'])&&sum(~isnan(CycAvg.([traces{tr},'_cycavg'])))>1
                %Now find maxvel, gain and latency
                eye_cyc = CycAvg.([traces{tr},'_cyc']);
                %define start and end as 10dps
                warning('off')
                head = neg_flag*spline(t,CycAvg.stim_cycavg,t_upsamp);
                eye = -neg_flag*spline(t,CycAvg.([traces{tr},'_cycavg']),t_upsamp);
                warning('on')
                [~,h_max] = max(head);
                h_start = find(head<10&t_upsamp<t_upsamp(h_max),1,'last');
                %Incase head zero crossing is weird, make the end the same size as start to max
                h_end = min([find(head<10&t_upsamp>t_upsamp(h_max),1,'first'),(h_max-h_start)*2+h_start]);
                head_start = t_upsamp(h_start);
                head_end = t_upsamp(h_end);
                avg_max_vel(tr) = max(eye(h_start:((h_max-h_start)*3+h_start)));
                %Gain by taking the ratio of area under the curve
                h_AUC = trapz(head(h_start:h_end))*median(diff(t_upsamp));
                e_AUC = trapz(eye(h_start:h_end))*median(diff(t_upsamp)); 
                avg_gain_AUC(tr) = e_AUC/h_AUC;
                %Eye Latency
                e_start = find(eye>(10+eye(h_start))&t_upsamp>=t_upsamp(h_start),1,'first');
                if isempty(e_start) || e_start > h_end
                    avg_lat(tr) = 1000*(t_upsamp(h_end)-t_upsamp(h_start)); %The full trace
                    e_start = h_start;
                else
                    avg_lat(tr) = 1000*(t_upsamp(e_start)-t_upsamp(h_start)); %in ms
                end
                %Gain by taking the the ratio of head accel/eye accel
                rel_eye = eye(e_start:h_end);
                [~,max_eye] = max(rel_eye);
                rel_eye = rel_eye(1:max_eye);
                ts_e = t_upsamp(1:length(rel_eye));
                eye_lin_fit = [ones(length(ts_e),1),ts_e']\rel_eye';
                rel_head = head(h_start:h_max);
                ts_h = t_upsamp(1:length(rel_head));
                head_lin_fit = [ones(length(ts_h),1),ts_h']\rel_head';
                avg_gain_Ga(tr) = eye_lin_fit(2)/head_lin_fit(2);
                %Look for saccades
                if any(contains(rel_tr,traces{tr})) %Look for saccades
                    eye_cyc_prefilt = -neg_flag*spline(t,CycAvg.([traces{tr},'_cyc_prefilt'])',t_upsamp);
                end
                for i = 1:nc
                    %define start and end as 10dps
                    warning('off')
                    head = neg_flag*spline(t,Stim_All(i,:),t_upsamp);
                    eye = -neg_flag*spline(t,eye_cyc(i,:),t_upsamp);
                    warning('on')
                    [max_head,h_max] = max(head);
                    h_start = find(head<10&t_upsamp<t_upsamp(h_max),1,'last');
                    %Incase head zero crossing is weird, make the end the same size as start to max
                    h_end = min([find(head<10&t_upsamp>t_upsamp(h_max),1,'first'),(h_max-h_start)*2+h_start]);
                    max_vel(i,tr) = max(eye(h_start:((h_max-h_start)*3+h_start)));
                    head_maxvel(i,1) = max_head;
                    head_width(i,1) = 1000*(t_upsamp(h_end)-t_upsamp(h_start));
                    %Gain by taking the ratio of area under the curve
                    h_AUC = trapz(head(h_start:h_end))*median(diff(t_upsamp));
                    e_AUC = trapz(eye(h_start:h_end))*median(diff(t_upsamp)); 
                    gain_AUC(i,tr) = e_AUC/h_AUC;
                    %Eye Latency
                    e_start = find(eye>(10+eye(h_start))&t_upsamp>=t_upsamp(h_start),1,'first');
                    if isempty(e_start) || e_start > h_end
                        lat(i,tr) = 1000*(t_upsamp(h_end)-t_upsamp(h_start)); %The full trace
                        e_start = h_start;
                    else
                        lat(i,tr) = 1000*(t_upsamp(e_start)-t_upsamp(h_start)); %in ms
                    end
                    %Gain by taking the the ratio of head accel/eye accel
                    rel_eye = eye(e_start:h_end);
                    [~,max_eye] = max(rel_eye);
                    rel_eye = rel_eye(1:max_eye);
                    ts_e = t_upsamp(1:length(rel_eye));
                    eye_lin_fit = [ones(length(ts_e),1),ts_e']\rel_eye';
                    rel_head = head(h_start:h_max);
                    ts_h = t_upsamp(1:length(rel_head));
                    head_lin_fit = [ones(length(ts_h),1),ts_h']\rel_head';
                    gain_Ga(i,tr) = eye_lin_fit(2)/head_lin_fit(2);
                    head_accel(i) = head_lin_fit(2);
                    if any(contains(rel_tr,traces{tr})) %Look for saccades
                        Ts = median(diff(t_upsamp));
                        eye_pf = eye_cyc_prefilt(i,:);
                        pos_eye = cumtrapz(eye_pf)*Ts;
                        pos_head = cumtrapz(head)*median(diff(t_upsamp));
                        % Saccade detection method
                        prom = 0.7;
                        step = 1e-3;
                        wind1 = floor(10e-3/Ts);
                        wind2 = floor(step/Ts);
                        t2 = (wind1+1):wind2:length(pos_eye);
                        dpos_eye = spline_filt(t_upsamp(t2),pos_eye(t2)-pos_eye(t2-wind1),t_upsamp,0.9999995);
                        [pk,loc,wid,prom2] = findpeaks(dpos_eye,t_upsamp,'MinPeakProminence',prom,'WidthReference','halfprom','MinPeakDistance',50e-3); %At least 25ms apart 
                        starts = zeros(1,length(loc));
                        for ii = 1:length(starts)
                            rel_thresh = dpos_eye<(pk(ii)-0.5*prom2(ii));
                            rel_thresh(t_upsamp>loc(ii)) = false;
                            starts(ii) = t_upsamp(find(rel_thresh,1,'last'));
                        end                                            
                        TF = starts(pk>prom&wid>10e-3&wid<150e-3&starts>t_upsamp(h_max)); %saccade between 10-150ms long
                        saccade_t = (TF-t_upsamp(h_start))*1000;                     
                        saccade_timing{i,contains(rel_tr,traces{tr})} = saccade_t';
                        if any(pos_eye(h_end:end)>pos_head(h_end:end))
                            time_to_target(i,contains(rel_tr,traces{tr})) = t_upsamp(find(pos_eye(h_end:end)>pos_head(h_end:end),1,'first'))*1000;             
                        end                          
                        % Plot and check
%                         subplot(2,1,1)
%                         plot(t_upsamp,pos_head,'k')
%                         hold on
%                         plot(t_upsamp,pos_eye,'b')
%                         plot(saccade_t/1000+t_upsamp(h_start),zeros(length(saccade_t),1),'go')
%                         hold off
%                         set(gca,'Xlim',[0 0.7])
%                         subplot(2,1,2)
%                         plot([0 0.7],prom*[1 1],'k:')
%                         hold on
%                         findpeaks(dpos_eye,t_upsamp,'MinPeakProminence',prom,'MinPeakDistance',50e-3,'WidthReference','halfprom','Annotate','extents')  
%                         hold off
%                         set(gca,'Xlim',[0 0.7])
%                         pause;
                    end
                end               
            end    
        end  
        keep_inds = CycAvg.Data_allcyc.keep_inds(:,CycAvg.keep_tr);
        head_3D = [CycAvg.Data.HeadVel_L,CycAvg.Data.HeadVel_R,CycAvg.Data.HeadVel_Z];
        h_i = NaN(nc,1);       
        for i = 1:nc
            [~,head_i] = max(abs(Stim_All(i,:)));
            h_i(i) = keep_inds(head_i,i);
        end     
        [~,~,head_misc] = calc_misalignment(stim_vect,head_3D(h_i,:));
        MaxVel = reshape([avg_max_vel;std(max_vel);avg_max_vel;std(max_vel)],1,[]);
        Gain = reshape([avg_gain_AUC;std(gain_AUC);avg_gain_Ga;std(gain_Ga)],1,[]);
        Latency = reshape([avg_lat;std(lat,'omitnan')],1,[]);
        cycle_params.Gain = array2table(gain_AUC,'VariableNames',traces);
        cycle_params.GainGa = array2table(gain_Ga,'VariableNames',traces);
        cycle_params.Latency = array2table(lat,'VariableNames',traces); 
        cycle_params.StimLen = head_width;
        cycle_params.StimMaxVel = head_maxvel;
        cycle_params.StimLinAccel = head_accel;
        %Find Latency bounds
        CycAvg.stim_start = head_start;       
        CycAvg.stim_end = head_end;
        switch axis_name
            case 'LHRH'
                cols = [1,2];                
            case 'RALP'
                cols = [5,6];
            case 'LARP'
                cols = [3,4];
            otherwise
                cols = 1:10;
        end     
        CycAvg.trace_start = CycAvg.stim_start+mean(reshape(lat(:,cols),[],1),'omitnan')/1000;
        cycle_params.EyeGain = max(gain_AUC(:,cols),[],2);
        cycle_params.EyeLat = mean(lat(:,cols),2,'omitnan');
        if ~isempty(CycAvg.detec_tr)
            cycle_params.GNODetec = length(CycAvg.detec_tr);
            cycle_params.PercDetec = length(CycAvg.detec_tr)/length(CycAvg.keep_tr);
        else
            cycle_params.GNODetec = 0;
            cycle_params.PercDetec = NaN;
        end
        cycle_params.TimeElapsed = CycAvg.Data_rawvel.t(end); 
        cycle_params.Saccades = vertcat(saccade_timing{:});
        cycle_params.TimeToTarget = mean(time_to_target,2,'omitnan');
        cycle_params.HeadMisc = head_misc;
    case 4 
        %% Pulse Train
        % Process Inputs
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
        % Find maximum velocity (for all directions) and misalignment
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
        max_vel = NaN(nc,length(traces));
        for i = 1:length(traces)
            max_vel(:,i) = CycAvg.([traces{i},'_cyc'])(:,I);
        end
        MaxVel = reshape([mean(max_vel);std(max_vel);mean(max_vel);std(max_vel)],1,[]);
        % Misalignment
        A = NaN(1,8);
        [A(1),A(2)] = calc_misalignment(-stim_vect,[CycAvg.ll_cyc(:,I) CycAvg.lr_cyc(:,I) CycAvg.lz_cyc(:,I)]);
        [A(3),A(4)] = calc_misalignment(-stim_vect,[CycAvg.ll_cyc(:,I) CycAvg.lr_cyc(:,I) CycAvg.lz_cyc(:,I)]);
        [A(5),A(6)] = calc_misalignment(-stim_vect,[CycAvg.rl_cyc(:,I) CycAvg.rr_cyc(:,I) CycAvg.rz_cyc(:,I)]);
        [A(7),A(8)] = calc_misalignment(-stim_vect,[CycAvg.rl_cyc(:,I) CycAvg.rr_cyc(:,I) CycAvg.rz_cyc(:,I)]);
        Align = repmat(A,nf,1);
        %Disconjugacy
        B = NaN(1,4);        
        [B(1),B(2)] = calc_misalignment([CycAvg.ll_cyc(:,I) CycAvg.lr_cyc(:,I) CycAvg.lz_cyc(:,I)],[CycAvg.rl_cyc(:,I) CycAvg.rr_cyc(:,I) CycAvg.rz_cyc(:,I)]);
        [B(3),B(4)] = calc_misalignment([CycAvg.ll_cyc(:,I) CycAvg.lr_cyc(:,I) CycAvg.lz_cyc(:,I)],[CycAvg.rl_cyc(:,I) CycAvg.rr_cyc(:,I) CycAvg.rz_cyc(:,I)]);
        Disc = repmat(B,nf,1);
        % Add fits to the CycAvg Struct
        for i = 1:length(traces)
            CycAvg.([traces{i},'_cycavg_fit']) = mean(max_vel(:,i))*ones(1,nt);
            CycAvg.([traces{i},'_cyc_fit']) = max_vel(:,i)*ones(1,nt);
        end
        cycle_params = [];
end
%% Make the table
labs = [{'Frequency(Hz)';'Amplitude(dps)';'PulseFreq(pps)';'PhaseDur(us)';'CurrentAmp(uA)';'Cycles'};...
    reshape(strcat('MaxVel_',repmat(upper(traces),4,1),repmat({'_HIGH';'_HIGH_sd';'_LOW';'_LOW_sd'},1,length(traces))),[],1);...
    reshape(strcat('Gain_',repmat(upper(traces),4,1),repmat({'_HIGH';'_HIGH_sd';'_LOW';'_LOW_sd'},1,length(traces))),[],1);...
    reshape(strcat('Tau_',repmat(upper(traces),4,1),repmat({'_HIGH';'_HIGH_sd';'_LOW';'_LOW_sd'},1,length(traces))),[],1);...
    reshape(strcat('RMSE_',repmat(upper(traces),2,1),repmat({'_HIGH';'_LOW'},1,length(traces))),[],1);...
    reshape(strcat('Latency_',repmat(upper(traces),2,1),repmat({'';'_sd'},1,length(traces))),[],1);...
    {'Phase_L';'Phase_L_sd';'Phase_R';'Phase_R_sd'};...
    {'Align_L_HIGH';'Align_L_HIGH_sd';'Align_L_LOW';'Align_L_LOW_sd';'Align_R_HIGH';'Align_R_HIGH_sd';'Align_R_LOW';'Align_R_LOW_sd'};...
    {'Disc_HIGH';'Disc_HIGH_sd';'Disc_LOW';'Disc_LOW_sd'}];
tab1 = cell2table(repmat({CycAvg.name,subject,visit,date,goggle,experiment,Type,condition,axis_name,{stim_vect},electrode},nf,1));
tab1.Properties.VariableNames = {'File','Subject','Visit','Date','Goggle','Experiment','Type','Condition','AxisName','StimAxis','Electrode'};
tab2 = array2table([freqs,repmat([amp,pps,phase_dur,curr_amp,nc],nf,1),MaxVel,Gain,Tau,RMSE,Latency,Phase,Align,Disc]);
tab2.Properties.VariableNames = labs;
results = [tab1,tab2];
CycAvg.parameterized = results;
CycAvg.cycle_params = cycle_params;
end