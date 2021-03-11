%% ParameterizeCycAvg.m 
%This function takes in a CycAvg file (either PJB with added info or AA)
%and parameterizes accordingly. It adds the results table and any
%appropriate fits to the CycAvg struct.
%Type tells you what kind of output to expect
%1. Sine Fit (Gain, Phase, Misalignment)
%2. Exponential Fit (Tau, Magnititude)
%3. Linear Fit (Gain, Latency, Saccade Analysis)
%4. Magnitude (Max Magnitude, Misalignment)
%These tables can ONLY be concatenated to other tables of the same type
%because they have different labels.

function [CycAvg,type] = ParameterizeCycAvg(CycAvg,fname)
%% Extract the relevant descriptive parameters
% This assumes the file is named properly. Please reference the
% makeNotes.m for proper naming conventions and use fixFileName.m as needed
if nargin < 2
    fname = CycAvg.name;
end
fparts = split(fname,'-');
% Subject
if contains(fname,'MVI')
    subject = strrep(fparts{contains(fparts,'MVI')},'CycAvg_','');
elseif any(~isnan(str2double(strrep(strrep(fparts,'CycAvg_',''),'R','')))&contains(fparts,'R')) %in form R#
    subject = strrep(fparts{~isnan(str2double(strrep(strrep(fparts,'CycAvg_',''),'R','')))&contains(fparts,'R')},'CycAvg_','');
else
    subject = '';
end
% Visit
if any(contains(fparts,'Visit'))
    visit = strrep(fparts{contains(fparts,'Visit')},' ','');
else
    visit = 'NA';
end
% Date
date = datetime(fparts{cellfun(@str2double,fparts)>20160000},'InputFormat','yyyyMMdd');
% Experiment, Condition, and Stim Axes
known_exps = {'RotaryChair','aHIT','eeVOR'};
if contains(fname,known_exps)
    experiment = known_exps(cellfun(@(x) contains(fname,x),known_exps));
    condition = strrep(strjoin(fparts(find(contains(fparts,experiment))+1:end)),'.mat','');        
else
    experiment = '';
    condition = '';
end   
%Stim Vector
if contains(condition,{'LARP','RP'})
    stim_vect = [1 0 0];
elseif contains(condition,'LA')
    stim_vect = [-1 0 0];  
elseif contains(condition,{'RALP','RA'})
    stim_vect = [0 1 0];
elseif contains(condition,'LP')
    stim_vect = [0 -1 0];
elseif contains(condition,{'LHRH','LH'})
    stim_vect = [0 0 1];
elseif contains(condition,'RH')
    stim_vect = [0 0 -1];
elseif contains(condition,'X')
    stim_vect = [0.707 0.707 0];
elseif contains(condition,'Y')
    stim_vect = [-0.707 0.707 0];
elseif contains(condition,'Activation')
    stim_vect = [0,0,0];
elseif contains(condition,'Vect') %MultiVector
    b1 = strfind(condition,'[');
    b2 = strfind(condition,']');
    stim_vect = str2num(strrep(condition(b1+1:b2-1),' ','-'));
else
    stim_vect = [0,0,0];
end    
%Make the type switch
if contains(condition,{'Sine','Sinusoid'}) %Sine fit
    type = 1;
    % Frequency
    freqs = fparts(contains(fparts,'Hz'));
    if isempty(freqs) %Nothing contains Hz so look for # after sinusoid (don't worry about multiple frequencies)
        freqs = fparts(find(contains(fparts,{'Sin','Sine','Sinusoid'}))+1);
    end
    nf = length(freqs);
    freq = zeros(1,nf);
    for i = 1:nf
        freqs{i} = [strrep(strrep(strrep(freqs{i},'p','.'),'.mat',''),'Hz',''),'Hz'];
        freq(i) = str2double(strrep(freqs{i},'Hz',''));
    end  
elseif contains(condition,{'VelStep','Activation'}) %Exponential fit
    type = 2;
    nf = 1;
    freqs = {''};
elseif contains(condition,{'Impulse'}) %Linear Fit
    type = 3;
    nf = 1;
    freqs = {''};
elseif contains(experiment,'eeVOR') %Pulse Stim
    type = 4;
    nf = 1;
    freqs = {''};
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
if ~isfield(CycAvg,'lx_cyc')
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
        makefit = @(trace) fminsearchbnd(@(p) sum((sine_fit(t,freq,p)-trace).^2,'omitnan'),repmat([abs(min(trace));abs(max(trace));0],nf,1),repmat([0;0;-180],1,nf),repmat([inf,inf,180],1,nf),options); 
        All_Cyc = NaN(nc,3*length(traces),nf);
        cyc_fit = NaN(nc,length(tt),length(traces));
        for i = 1:nc
            for j = 1:length(traces)
                if ~all(isnan(CycAvg.([traces{j},'_cyc'])(i,:)))
                    All_Cyc(i,(3*j-2):3*j,:) = reshape(makefit(spline(tt,CycAvg.([traces{j},'_cyc'])(i,:),t)),1,3,nf);
                    cyc_fit(i,:,j) = sine_fit(tt,freq,makefit(spline(tt,CycAvg.([traces{j},'_cyc'])(i,:),t)));
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
            RMSE(:,[i,i+length(traces)]) = sum((CycAvg.([traces{i},'_cycavg'])-sine_fit(tt,freq,makefit(spline(tt,CycAvg.([traces{i},'_cycavg']),t)))).^2)/length(tt);
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
            CycAvg.([traces{i},'_cycavg_fit']) = sine_fit(tt,freq,makefit(spline(tt,CycAvg.([traces{i},'_cycavg']),t)));
            CycAvg.([traces{i},'_cyc_fit']) = reshape(cyc_fit(:,:,i),nc,nt);
        end
    case 2 
        %% Exponential
        nc = 1;
        t = reshape(CycAvg.t,1,[]);
        %Make sure stim trace is time points long (1 x nt)
        if max(abs(CycAvg.stim)) > 1 %Triggered by motion. Should be one high/low and then a period of 0
            stim = reshape(CycAvg.stim,1,[]);
            max_vel = str2double(strrep(strrep(fparts{contains(fparts,'dps')},'dps.mat',''),'n','-'));
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
        % Find exponential fits
        makefit = @(tt,p) p(:,1).*exp(p(:,2).*tt)+p(:,3).*exp(p(:,4).*tt); %2nd order exponential fit same syntax as matlab's built-in         
        all_fits = NaN(length(traces),length(t));
        params = NaN(2*length(traces),4);
        dconfint = NaN(2*length(traces),4);
        high_i = t(find(Stim==1,1,'first'));
        low_i = t(find(Stim==0,1,'first'));
        for i = 1:length(traces)
            tt = t(sub_i);
            long_dat = CycAvg.([traces{i},'_cyc']);
            dat = long_dat(sub_i);
            %Stim is high
            t_high = tt(Stim(sub_i)==1&~isnan(dat));
            [fitobj,gof] = fit(t_high'-high_i,dat(Stim(sub_i)==1&~isnan(dat))','exp2');
            RMSE(1,i) = gof.rmse;
            params(i,:) = coeffvalues(fitobj);
            dconfint(i,:) = mean(abs(confint(fitobj)-mean(confint(fitobj))));
            %Stim is low
            t_low = tt(Stim(sub_i)==0&~isnan(dat));
            [fitobj,gof] = fit(t_low'-low_i,dat(Stim(sub_i)==0&~isnan(dat))','exp2');
            RMSE(1,length(traces)+i) = gof.rmse;
            params(length(traces)+i,:) = coeffvalues(fitobj);
            dconfint(length(traces)+i,:) = mean(abs(confint(fitobj)-mean(confint(fitobj))));
            %Trace fit
            high_fit = makefit(t-high_i,params(i,1:4));
            low_fit = makefit(t-low_i,params(length(traces)+i,1:4));
            all_fits(i,Stim==1) = high_fit(Stim==1);
            all_fits(i,Stim==0) = low_fit(Stim==0);
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
    case 3 
        %% Linear -- IN PROGRESS

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
end
%% Make the table
labs = [{'Cycles'};...
    reshape(strcat('MaxVel_',repmat(upper(traces),4,1),repmat({'_HIGH';'_HIGH_sd';'_LOW';'_LOW_sd'},1,length(traces))),[],1);...
    reshape(strcat('Gain_',repmat(upper(traces),4,1),repmat({'_HIGH';'_HIGH_sd';'_LOW';'_LOW_sd'},1,length(traces))),[],1);...
    reshape(strcat('Tau_',repmat(upper(traces),4,1),repmat({'_HIGH';'_HIGH_sd';'_LOW';'_LOW_sd'},1,length(traces))),[],1);...
    reshape(strcat('RMSE_',repmat(upper(traces),2,1),repmat({'_HIGH';'_LOW'},1,length(traces))),[],1);...
    reshape(strcat('Latency_',repmat(upper(traces),2,1),repmat({'';'_sd'},1,length(traces))),[],1);...
    {'Phase_L';'Phase_L_sd';'Phase_R';'Phase_R_sd'};...
    {'Align_L_HIGH';'Align_L_HIGH_sd';'Align_L_LOW';'Align_L_LOW_sd';'Align_R_HIGH';'Align_R_HIGH_sd';'Align_R_LOW';'Align_R_LOW_sd'};...
    {'Disc_HIGH';'Disc_HIGH_sd';'Disc_LOW';'Disc_LOW_sd'}];
tab1 = cell2table([repmat({fname,subject,visit,date,goggle,experiment,condition},nf,1),freqs]);
tab1.Properties.VariableNames = {'File','Subject','Visit','Date','Goggle','Experiment','Condition','Frequency'};
tab2 = array2table([nc*ones(nf,1),MaxVel,Gain,Tau,RMSE,Latency,Phase,Align,Disc]);
tab2.Properties.VariableNames = labs;
results = [tab1,tab2];
CycAvg.parameterized = results;  
end