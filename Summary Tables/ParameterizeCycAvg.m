%% ParameterizeCycAvg.m
%This function takes in a CycAvg file (either PJB with added info or AA)
%and parameterizes accordingly. It adds the results table and any
%appropriate fits to the CycAvg struct.

%Type tells you what kind of output to expect:
%1. Sine Fit (Magnitude, Half-Cycle Fit Gain, Half-Cycle Fit Phase, Misalignment)
%2. Exponential Fit (Magnititude, Tau, Misalignment)
%3. Impulse (Magnitude, Gain, Latency)
%4. Pulse Train (Magnitude, Misalignment)

%Impulses and Pulse Trains already have a directionality to them. N = 1.
%Sine and exponential fit have half-cycles or high/low periods,
%so each of those files will generate a table where N=2.
%Sum of sine files are no longer supported.

function [CycAvg,type] = ParameterizeCycAvg(CycAvg)
%% Initialize table
options = optimset('Display','off'); %suppress output from fminsearch
traces = {'lz','rz','ll','rl','lr','rr','lx','rx','ly','ry'};
cell_labs = {'File','Subject','Visit','Date','Goggle','Experiment',...
    'Type','Condition','AxisName','AxisLetter','StimAxis','Electrode'};
param_tr = {'MaxVel','Gain','Tau','Latency','Phase','RMSE'};
sub_num_labs = [reshape(strcat(repmat(param_tr,length(traces),1),'_',...
    repmat(upper(traces)',1,length(param_tr))),[],1)',...
    {'Align_L','Align_R','Disc'},param_tr,{'Align'}];
num_labs = [{'Frequency';'Amplitude';'PulseFreq';'PhaseDur';'CurrentAmp';'Cycles'};...
    reshape([sub_num_labs;strcat(sub_num_labs,'_sd')],[],1)]';
%Load in file name and use to figure out how many rows are needed in the table
old_fname = CycAvg.name;
fname = strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(old_fname,'CurrentFitting','eeVOR'),...
    'PreAct_ElectrodeAutoscan','eeVOR'),'ElectricalOnly','eeVOR'),'ElectricOnly','eeVOR'),'CycAvg_',''),'.mat',''),'_CycleAvg',''),'NA','');
fname = fname(1:min([strfind(fname,'_'),length(fname)+1])-1);
fparts = split(strrep(fname,'.mat',''),'-');
fparts(cellfun(@isempty,fparts)) = [];
N = double(contains(lower(fname),{'sin','velstep','activation'}))+1; %Pos/Neg Half-cycles;
% Initialize table and structure
cycle_params = struct;
results = [cell2table(repmat({''},N,length(cell_labs)),'VariableNames',cell_labs),...
    array2table(NaN(N,length(num_labs)),'VariableNames',num_labs)];
results.Date = NaT(N,1); %Fix date into a datetime
results.Properties.VariableUnits = repmat({''},1,length(results.Properties.VariableNames));
results.Properties.VariableUnits(contains(results.Properties.VariableNames,'StimAxis')) = {'[L,R,Z] Normalized rotation vector'};
results.Properties.VariableUnits(contains(results.Properties.VariableNames,'Frequency')) = {'Hz'};
results.Properties.VariableUnits(contains(results.Properties.VariableNames,'PulseFreq')) = {'pulse/s'};
results.Properties.VariableUnits(contains(results.Properties.VariableNames,'PhaseDur')) = {'us/phase'};
results.Properties.VariableUnits(contains(results.Properties.VariableNames,'CurrentAmp')) = {'uA'};
results.Properties.VariableUnits(contains(results.Properties.VariableNames,'Tau')) = {'s'};
results.Properties.VariableUnits(contains(results.Properties.VariableNames,'Latency')) = {'ms'};
results.Properties.VariableUnits(contains(results.Properties.VariableNames,{'Amplitude','MaxVel','RMSE'})) = {'deg/s'};
results.Properties.VariableUnits(contains(results.Properties.VariableNames,{'Phase','Align','Disc'})) = {'deg'};
%% Extract the relevant descriptive parameters
results.File(:) = {old_fname};
% Subject
sub_pat = "R"+digitsPattern(3);
ind = 1;
if contains(fname,sub_pat)
    ind = contains(fparts,extract(fname,sub_pat));
end
results.Subject(:) = fparts(ind);
fparts(ind) = [];
% Visit
results.Visit(:) = {'NA'};
if any(contains(fparts,'Visit'))
    results.Visit(:) = strrep(fparts(contains(fparts,'Visit')),' ','');
    fparts(contains(fparts,'Visit')) = [];
end
% Date
date_ind = find(cellfun(@str2double,fparts)>20160000);
if ~isnan(str2double(fparts{date_ind+1})) %Time
    results.Date(:) = datetime([fparts{date_ind},'-',fparts{date_ind+1}],'InputFormat','yyyyMMdd-HHmmss');
    fparts(date_ind:date_ind+1) = [];
elseif ~isempty(date_ind)
    results.Date(:) = datetime(fparts{date_ind},'InputFormat','yyyyMMdd');
    fparts(date_ind) = [];
else %Wrong/missing date format
    error(['Improper/missing date format. Expecting yyyyMMdd. : ',old_fname])
end
%Goggle
if ~isfield(CycAvg,'info')||~isfield(CycAvg.info,'goggle_ver')
    goggle = {'LDVOG'};
elseif isnumeric(CycAvg.info.goggle_ver)
    goggle = {['LDVOG',num2str(CycAvg.info.goggle_ver)]};
else
    goggle = {CycAvg.info.goggle_ver};
end
results.Goggle(:) = goggle;
fparts(contains(fparts,goggle)) = [];
% Experiment
known_exps = {'RotaryChair','aHIT','eeVOR','Manual'};
if contains(fname,known_exps)
    experiment = known_exps(cellfun(@(x) contains(fname,x),known_exps));
    results.Experiment(:) = experiment;
    fparts(contains(fparts,experiment)) = [];
else
    error(['Cannot parameterize and assign experiment type to: ',old_fname])
end
%Amplitude
if any(contains(fparts,'dps'))
    results.Amplitude(:) = str2double(strrep(strrep(fparts{contains(fparts,'dps')},'n','-'),'dps',''));
    fparts(contains(fparts,'dps')) = [];
end
%Pulses Frequency
if any(contains(fparts,'pps'))
    results.PulseFreq(:) = str2double(strrep(fparts{contains(fparts,'pps')},'pps',''));
    fparts(contains(fparts,'pps')) = [];
end
%PhaseDur
if any(contains(fparts,{'us','uS'}))
    results.PhaseDur(:) = str2double(strrep(lower(fparts{contains(fparts,{'us','uS'})}),'us',''));
    fparts(contains(fparts,{'us','uS'})) = [];
end
%Current Amplitude
if any(contains(fparts,{'ua','uA'}))
    results.CurrentAmp(:) = str2double(strrep(lower(fparts{contains(fparts,{'ua','uA'})}),'ua',''));
    fparts(contains(fparts,{'ua','uA'})) = [];
end
%Stim Vector and Axis Name
n = [1,2];
if contains(fname,{'LA','LP','RH'})&&~contains(fname,{'LARP','LHRH','RALP','X','Y'}) %Negative
    n = [2,1];
end
Rotations = table();
Rotations.Name1 = {{'RP','LA'},{'RA','LP'},{'LH','Rotary','RH'},'X','Y'}';
Rotations.Name2 = {{'RP','LA'},{'RA','LP'},{'LH','RH'},{'+X','-X'},{'+Y','-Y'}}';
Rotations.AxisLetter = {'L','R','Z','X','Y'}';
Rotations.Val = [1 0 0;0 1 0;0 0 1;0.707 0.707 0;-0.707 0.707 0];
is_common_vect = find(cellfun(@(x) contains(fname,x),Rotations.Name1),1,'first');
results.StimAxis(:) = {[0,0,0]};
if ~isempty(is_common_vect)    
    rel_names = Rotations.Name2{is_common_vect}(n);
    rel_vals = [Rotations.Val(is_common_vect,:);-1*Rotations.Val(is_common_vect,:)];
    rel_vals = rel_vals(n,:);
    for i = 1:size(results,1)
       results.AxisName{i} = rel_names{i};
        results.StimAxis{i} = rel_vals(i,:);
    end
    results.AxisLetter(:) = Rotations.AxisLetter(is_common_vect);
elseif contains(fname,'[')&&contains(fname,']') %MultiVector, SineMultivector
    results.StimAxis(:) = {str2double(split(fname(strfind(fname,'[')+1:strfind(fname,']')-1),','))'};
elseif contains(fname,'Activation')&&contains(fname,{'LA','LP','LH'})
    results.StimAxis(:) = {[-1,-1,1]/sqrt(3)};
elseif contains(fname,'Activation')&&contains(fname,{'RA','RP','RH'})
    results.StimAxis(:) = {[1,1,-1]/sqrt(3)};
end
% Electrode #
if contains(fname,{'LAE','LHE','LPE','RAE','RHE','RPE'})
    results.Electrode(:) = fparts(contains(fparts,{'LAE','LHE','LPE','RAE','RHE','RPE'}));
end
fparts(contains(fparts,{'LA','LH','LP','RA','RH','RP','['})) = [];
% Type/Condition
Types = {'Sine','Exponential','Impulse','PulseTrain','Activation'};
is_known_type = [contains(fname,'Sin'),contains(fname,'VelStep'),...
    contains(fname,{'Impulse','Gaussian'}),contains(fname,'eeVOR')];
if ~any(is_known_type)&&~contains(fname,'Activation')
    error(['Unknown experiemnt type. Cannot parameterize and assign experiment type to: ',old_fname])
end
type = 5; %Activation
if ~contains(fname,'Activation')
    type = find(is_known_type,1,'first'); %type 4 should be eeVOR is isn't of another type already, type 5 is activation
end
if type == 1
    %Frequency
    freqs = strrep(strrep(fparts(contains(fparts,'Hz')),'p','.'),'.mat','');
    if isempty(freqs) %Nothing contains Hz so look for # after sinusoid (don't worry about multiple frequencies)
        freqs = fparts(find(contains(fparts,'Sin'))+1);
        fparts(find(contains(fparts,{'Sin','Sine','Sinusoid'}))+1) = [];
    end
    freqs = str2double(strrep(freqs,'Hz',''));
    results.Frequency(:) = freqs;
elseif type == 2 && contains(fname,'VelStep')
    %On and "off" period for exponential
    results.StimAxis{2} = 0*results.StimAxis{2};
    results.AxisName{2} = '';
end
fparts(contains(fparts,{'Hz','Sin','VelStep','Activation','Impulse','Gaussian','eeVOR'})) = [];
results.Type(:) = Types(type);
%Condition
if contains(lower(fname),'autoscan')
    condition = 'Autoscan';
elseif contains(lower(fname),'eevor')
    condition = 'eeVOR';
elseif contains(lower(fname),{'baseline','constant'})
    condition = 'ConstantRate';
elseif contains(lower(fname),{'motion','mod'})
    condition = 'MotionMod';
elseif contains(lower(fname),{'nostim','dark','off'})||isempty(fparts) %if isempty, probably a normal subject
    condition = 'NoStim';
else
    condition = strrep(strjoin(fparts,' '),'_','');
end
if contains(lower(fname),'light')
    results.Condition(:) = {[condition,'Light']};
else
    results.Condition(:) = {condition};
end
%% Create a table for each type
switch type
    case 1
        %% Sine
        % Ensure equal cycle size by padding with NaNs if needed
        [nc1,nt] = size(CycAvg.lz_cyc);
        nc2 = size(CycAvg.rz_cyc,1);
        nc = max([nc1,nc2]);
        results.Cycles(:) = nc;
        %Figure out which traces are used
        rm_tr = true(1,length(traces));
        for tr = 1:length(traces)
            if isfield(CycAvg,[traces{tr},'_cyc'])&&all(any(~isnan(CycAvg.([traces{tr},'_cyc']))))
                CycAvg.([traces{tr},'_cyc']) = [CycAvg.([traces{tr},'_cyc']);NaN(nc-size(CycAvg.([traces{tr},'_cyc']),1),nt)];
                rm_tr(tr) = 0;
            end
        end
        traces = traces(~rm_tr);
        %Make sure stim trace is time points long (1 x nt)
        if ismember(size(CycAvg.stim,2),[1,nc1,nc2])
            CycAvg.stim = CycAvg.stim';
        end
        [ab,ac] = size(CycAvg.stim);
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
        if ~isfield(CycAvg,'t')
            CycAvg.t = tt;
        end
        t = tt(unique(round(linspace(1,length(Stim_CycAvg),100))));
        % Find gain, phase, and misalignment for each trace
        makefit = @(trace) fminsearchbnd(@(p) sum((sine_fit(t,freqs,p)-trace).^2,'omitnan'),...
            [min(trace)-median(trace);max(trace)-median(trace);0;median(trace)],...
            [-inf;-inf;-180;-inf],[inf;inf;180;inf],options);        
        % Find index of interest for maximum and minimum
        L_Eye = permute(reshape([CycAvg.ll_cyc,CycAvg.lr_cyc,CycAvg.lz_cyc],nc,[],3),[3,2,1]);
        R_Eye = permute(reshape([CycAvg.rl_cyc,CycAvg.rr_cyc,CycAvg.rz_cyc],nc,[],3),[3,2,1]);
        [~,I1] = max(max([-results.StimAxis{1}*mean(L_Eye,3,'omitnan');-results.StimAxis{1}*mean(R_Eye,3,'omitnan')]));
        [~,I2] = max(max([-results.StimAxis{2}*mean(L_Eye,3,'omitnan');-results.StimAxis{2}*mean(R_Eye,3,'omitnan')]));
        for tr = 1:length(traces)
            cyc_param = NaN(nc,4);
            cyc_fit = NaN(nc,length(tt));
            cyc_tr = CycAvg.([traces{tr},'_cyc']); 
            cyc_avg = mean(cyc_tr,1,'omitnan');
            for i = 1:nc
                if sum(~isnan(cyc_tr(i,:))) > 2
                    cyc_param(i,:) = makefit(spline(tt,cyc_tr(i,:),t));
                    cyc_fit(i,:) = sine_fit(tt,freqs,cyc_param(i,:));
                end
            end            
            cycavg_param = makefit(spline(tt,cyc_avg,t));
            tr_fit = sine_fit(tt,freqs,cycavg_param);
            err_tr_fit = (tr_fit-cyc_avg).^2;
            err_cyc_fit = (cyc_fit-cyc_tr).^2;
            t0 = sin(2*pi*freqs*t+cycavg_param(3)*pi/180)>=0; %pos half-cycle
            results.(['MaxVel_',upper(traces{tr})]) = [max(abs(cyc_avg(t0)));max(abs(cyc_avg(~t0)))];
            results.(['MaxVel_',upper(traces{tr}),'_sd']) = [std(max(abs(cyc_tr(:,t0)),[],2),'omitnan');std(max(abs(cyc_tr(:,~t0)),[],2),'omitnan')];
            results.(['Gain_',upper(traces{tr})]) = [-1;1].*cycavg_param(1:2)/ChairAmp;
            results.(['Gain_',upper(traces{tr}),'_sd']) = [std(cyc_param(:,1),'omitnan');std(cyc_param(:,2),'omitnan')]/ChairAmp;
            results.(['Phase_',upper(traces{tr})])(:) = cycavg_param(3);
            results.(['Phase_',upper(traces{tr}),'_sd'])(:) = std(cyc_param(:,3),'omitnan');
            results.(['RMSE_',upper(traces{tr})])(:) = [sqrt(mean(err_tr_fit(:,t0),'omitnan'));sqrt(mean(err_tr_fit(:,~t0),'omitnan'))];
            results.(['RMSE_',upper(traces{tr}),'_sd'])(:) = [std(sqrt(mean(err_cyc_fit(:,t0),2,'omitnan')),'omitnan');std(sqrt(mean(err_cyc_fit(:,~t0),2,'omitnan')),'omitnan')];
            CycAvg.([traces{tr},'_cyc_fit']) = cyc_fit;
            CycAvg.([traces{tr},'_cycavg_fit']) = tr_fit;
        end
        % Misalignment
        [results.Align_L(1),results.Align_L_sd(1)] = calc_misalignment(-results.StimAxis{1},permute(L_Eye(:,I1,:),[3,1,2]));
        [results.Align_R(1),results.Align_R_sd(1)] = calc_misalignment(-results.StimAxis{1},permute(R_Eye(:,I1,:),[3,1,2]));
        [results.Align_L(2),results.Align_L_sd(2)] = calc_misalignment(-results.StimAxis{2},permute(L_Eye(:,I2,:),[3,1,2]));
        [results.Align_R(2),results.Align_R_sd(2)] = calc_misalignment(-results.StimAxis{2},permute(R_Eye(:,I2,:),[3,1,2]));
        %Disconjugacy
        [results.Disc(1),results.Disc_sd(1)] = calc_misalignment(permute(L_Eye(:,I1,:),[3,1,2]),permute(R_Eye(:,I1,:),[3,1,2]));
        [results.Disc(2),results.Disc_sd(2)] = calc_misalignment(permute(L_Eye(:,I2,:),[3,1,2]),permute(R_Eye(:,I2,:),[3,1,2]));        
        %Put the median analysis in cycle params file
        cycle_params.t = CycAvg.t;
        cycle_params.stim = CycAvg.stim;
        if isfield(CycAvg,'Data')
            Data_med = angpos2angvel(CycAvg.Data);
            cyc_inds = CycAvg.Data_allcyc.keep_inds;
            for tr = 1:length(traces)
                cycle_params.([traces{tr},'_cycavg']) = CycAvg.([traces{tr},'_cycavg']);
                cycle_params.([traces{tr},'_cycstd']) = CycAvg.([traces{tr},'_cycstd']);
                cycle_params.([traces{tr},'_cyc']) = CycAvg.([traces{tr},'_cyc']);
                long_name = [upper(traces{tr}(1)),'E_Vel_',upper(traces{tr}(2))];   
                long_name = strrep(strrep(long_name,'_L','_LARP'),'_R','_RALP');
                cycle_params.([traces{tr},'_cycmed']) = median(Data_med.(long_name)(cyc_inds),2,'omitnan');
            end
        end
    case 2
        %% Exponential
        results.Cycles(:) = 1;
        t = reshape(CycAvg.t,1,[]);
        %Make sure stim trace is time points long (1 x nt)
        if max(abs(CycAvg.stim)) > 1 %Triggered by motion. Should be one high/low and then a period of 0
            stim = reshape(CycAvg.stim,1,[]);
            max_vel = results.Amplitude(1);
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
        % Find exponential fits: Different Model Types
        %2nd order exponential fit same syntax as matlab's built-in
        makefit = @(tt,p) p(:,1).*exp(p(:,2).*tt)+p(:,3).*exp(p(:,4).*tt);
        %1order exponential fit with decay to 0
        exp_fit1 = @(tt,p) p(1).*exp(-tt./p(2));
        makefit2 = @(trace,T) fminsearch(@(p) sum((exp_fit1(T,p)-trace).^2,'omitnan')+sum(p.^2),[median(trace(1));1],options);
        %1order exponential fit with decay not forced to 0.
        exp_fit2 = @(tt,p) p(1).*exp(-tt./p(2)) + p(3);
        makefit3 = @(trace,T) fminsearch(@(p) sum((exp_fit2(T,p)-trace).^2,'omitnan')+sum(p.^2),[median(trace);1;0],options);
        %Intialize arrays for the different fits
        all_fits = NaN(length(traces),length(t));
        all_fits2 = NaN(length(traces),length(t));
        all_fits3 = NaN(length(traces),length(t));
        %Intialize tables for the different model types
        params = array2table(NaN(8,length(traces)),'VariableNames',traces,'RowNames',{'K1_High','Tau1_High','K2_High','Tau2_High','K1_Low','Tau1_Low','K2_Low','Tau2_Low'});
        params2 = array2table(NaN(4,length(traces)),'VariableNames',traces,'RowNames',{'K_High','Tau_High','K_Low','Tau_Low'});
        params3 = array2table(NaN(6,length(traces)),'VariableNames',traces,'RowNames',{'K_High','Tau_High','C_High','K_Low','Tau_Low','C_Low'});
        mult = [-1;1];
        if contains(results.AxisName(1),{'RH','LA','LP'}) %Negative directions
            mult = [1;-1];
        end        
        for i = 1:length(traces)
            if isfield(CycAvg,[traces{i},'_cyc'])
                results.(['RMSE_',upper(traces{i}),'_sd'])(:) = 0; %Only 1 cyc
                results.(['Tau_',upper(traces{i}),'_sd'])(:) = 0; %Only 1 cyc
                results.(['MaxVel_',upper(traces{i}),'_sd'])(:) = 0; %Only 1 cyc
                results.(['Gain_',upper(traces{i}),'_sd'])(:) = 0; %Only 1 cyc
                long_dat = reshape(CycAvg.([traces{i},'_cyc']),[],1);
                dat = long_dat(sub_i);
                t_high = tt(Stim(sub_i)'==1&~isnan(dat))'-high_i;
                t_low = tt(Stim(sub_i)'==0&~isnan(dat))'-low_i;
                dat_high = dat(Stim(sub_i)'==1&~isnan(dat));
                dat_low = dat(Stim(sub_i)'==0&~isnan(dat));
                %Stim is high
                [fitobj_h,gof_h] = fit(t_high,dat_high,'exp2', 'Lower', [0 0 0 0]);
                params{1:4,i} = coeffvalues(fitobj_h)';
                %dconfint_h = mean(abs(confint(fitobj_h)-mean(confint(fitobj_h))));
                params2{1:2,i} = makefit2(dat_high,t_high);
                params3{1:3,i} = makefit3(dat_high,t_high);
                all_fits(i,Stim==1) = makefit(t(Stim==1)-high_i,params{1:4,i}');
                all_fits2(i,Stim==1) = exp_fit1(t(Stim==1)-high_i,params2{1:2,i}');
                all_fits3(i,Stim==1) = exp_fit2(t(Stim==1)-high_i,params3{1:3,i}');
                gof_h.rmse = sqrt(mean((all_fits2(i,Stim==1)-long_dat(Stim==1)').^2,'omitnan'));
                %Stim is low
                %dconfint_l = NaN*dconfint_h;
                %gof_l = gof_h;
                gof_l.rmse = NaN;
                if length(t_low)>4
                    [fitobj_l,gof_l] = fit(t_low,dat_low,'exp2', 'Lower', [0 0 0 0]);
                    params{5:8,i} = coeffvalues(fitobj_l)';
                    %dconfint_l = mean(abs(confint(fitobj_l)-mean(confint(fitobj_l))));
                    params2{3:4,i} = makefit2(dat_low,t_low);
                    params3{4:6,i} = makefit3(dat_low,t_low);
                    all_fits(i,Stim==0) = makefit(t(Stim==0)-low_i,params{5:8,i}');
                    all_fits2(i,Stim==0) = exp_fit1(t(Stim==0)-low_i,params2{3:4,i}');
                    all_fits3(i,Stim==0) = exp_fit2(t(Stim==0)-low_i,params3{4:6,i}');
                    gof_l.rmse = sqrt(mean((all_fits2(i,Stim==0)-long_dat(Stim==0)').^2,'omitnan'));
                end                
                results.(['RMSE_',upper(traces{i})]) = [gof_h.rmse;gof_l.rmse];
                results.(['MaxVel_',upper(traces{i})]) = params2{[1,3],i};
                results.(['Tau_',upper(traces{i})]) = params2{[2,4],i};
                results.(['Gain_',upper(traces{i})]) = results.(['MaxVel_',upper(traces{i})])./(results.Amplitude(1)*mult);                                 
                CycAvg.([traces{i},'_cycavg_fit']) = all_fits2(i,:);
                CycAvg.([traces{i},'_cyc_fit']) = all_fits2(i,:);
            end
        end
        % Misalignment
        all_high_L = all_fits([3,5,1],Stim==1);
        all_low_L = all_fits([3,5,1],Stim==0);
        all_high_R = all_fits([4,6,2],Stim==1);
        all_low_R = all_fits([4,6,2],Stim==0);
        [results.Align_L(1),results.Align_L_sd(1)] = calc_misalignment(-results.StimAxis{1},all_high_L');
        [results.Align_L(2),results.Align_L_sd(2)] = calc_misalignment(-results.StimAxis{2},all_low_L');
        [results.Align_R(1),results.Align_R_sd(1)] = calc_misalignment(-results.StimAxis{1},all_high_R');
        [results.Align_R(2),results.Align_R_sd(2)] = calc_misalignment(-results.StimAxis{2},all_low_R');
        %Disconjugacy
        [results.Disc(1),results.Disc_sd(1)] = calc_misalignment(all_high_L',all_high_R');
        [results.Disc(2),results.Disc_sd(2)] = calc_misalignment(all_low_L',all_low_R');
        cycle_params.exp_ord2_params = params;
        cycle_params.exp_ord2_fits = all_fits;
        cycle_params.exp_ord1_params = params2;
        cycle_params.exp_ord1_fits = all_fits2;
        cycle_params.exp_ord1_const_params = params3;
        cycle_params.exp_ord1_const_fits = all_fits3;
    case 3
        %% Impulse
        % Process Inputs
        stim_vect = results.StimAxis{:};
        %Make sure stim trace is time points long (1 x nt)
        if isfield(CycAvg,'rz_cyc')
            [nc,nt] = size(CycAvg.rz_cyc);
        else
            [nc,nt] = size(CycAvg.lz_cyc);
        end
        results.Cycles(:) = nc;
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
                avg_max_vel(tr) = max(abs(eye(h_start:((h_max-h_start)*3+h_start))));
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
               end
                %Update table for trace-dependent things
                results.(['MaxVel_',upper(traces{tr})]) = avg_max_vel(tr);
                results.(['MaxVel_',upper(traces{tr}),'_sd']) = std(max_vel(:,tr),'omitnan');
                results.(['Gain_',upper(traces{tr})]) = avg_gain_AUC(tr);
                results.(['Gain_',upper(traces{tr}),'_sd']) = std(gain_AUC(:,tr),'omitnan');
                results.(['Latency_',upper(traces{tr})]) = avg_lat(tr);
                results.(['Latency_',upper(traces{tr}),'_sd']) = std(lat(:,tr),'omitnan');
            end
        end
        fields = fieldnames(CycAvg);
        % Misalignment
        % Only calculate if there are 3D traces for either eye
        if all(ismember({'ll_cyc','lr_cyc','lz_cyc'},fields))
            LEye = -neg_flag*[CycAvg.ll_cycavg;CycAvg.lr_cycavg;CycAvg.lz_cycavg];
            [~,i_lmax] = max(stim_vect*LEye);
            [results.Align_L(1),results.Align_L_sd(1)] = calc_misalignment(stim_vect,LEye(:,i_lmax)');
        end
        if all(ismember({'rl_cyc','rr_cyc','rz_cyc'},fields))
            REye = -neg_flag*[CycAvg.rl_cycavg;CycAvg.rr_cycavg;CycAvg.rz_cycavg];
            [~,i_rmax] = max(stim_vect*REye);
            [results.Align_R(1),results.Align_R_sd(1)] = calc_misalignment(stim_vect,REye(:,i_rmax)');
        end
        % Disconjugacy
        % Only calculate if there are 3D traces for BOTH eyes
        if all(ismember({'ll_cyc','lr_cyc','lz_cyc','rl_cyc','rr_cyc','rz_cyc'},fields))
            nc = min([size(LEye,1),size(REye,1)]); %in case of unequal cycles for left and right eye
            [results.Disc(1),results.Disc_sd(1)] = calc_misalignment(LEye(1:nc,i_lmax),REye(1:nc,i_rmax));
        end
        keep_inds = CycAvg.Data_allcyc.keep_inds(:,CycAvg.keep_tr);
        head_3D = [reshape(CycAvg.Data.HeadVel_L,[],1),reshape(CycAvg.Data.HeadVel_R,[],1),reshape(CycAvg.Data.HeadVel_Z,[],1)];
        h_i = NaN(nc,1);
        for i = 1:nc
            [~,head_i] = max(abs(Stim_All(i,:)));
            h_i(i) = keep_inds(head_i,i);
        end
        [~,~,head_misc] = calc_misalignment(stim_vect,head_3D(h_i,:));
        %Cycle Params Cell
        cycle_params.Gain = array2table(gain_AUC,'VariableNames',traces);
        cycle_params.GainGa = array2table(gain_Ga,'VariableNames',traces);
        cycle_params.Latency = array2table(lat,'VariableNames',traces);
        cycle_params.StimLen = head_width;
        cycle_params.StimMaxVel = head_maxvel;
        cycle_params.StimLinAccel = head_accel;
        %Find Latency bounds
        CycAvg.stim_start = head_start;
        CycAvg.stim_end = head_end;
        switch results.AxisName{1}
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
        %See if there was a CSV file with the gain values computed by the
        %system (GNO)
        if contains(results.Goggle(1),'GNO')&&~isempty(CycAvg.Data.CSVData)
            CSVData = CycAvg.Data.CSVData;
            if istable(CSVData)
                CSVData = table2cell(CSVData);
            end
            if any(any(cellfun(@isnumeric,CSVData)))
                rep_gains = CSVData(contains(CSVData(:,2),'Gain'),3);
            else
                if size(CSVData,2)>1                
                    CSVData = join(CSVData,',');
                end
                rep_gains = str2double(strrep(strrep(CSVData(contains(CSVData,'Gain')),'Gain',''),',',''));
                rep_dir = CSVData(find(contains(CSVData,'Gain'))-1);
            end
            if contains(fname,{'LA','LH','LP'})
                rel_canals = {'LA','LH','LP','Left'};
            else
                rel_canals = {'RA','RH','RP','Right'};
            end
            cycle_params.GoggleGain_Avg = round(mean(rep_gains(contains(rep_dir,rel_canals)),'omitnan'),2);
            cycle_params.GoggleGain_SD = round(std(rep_gains(contains(rep_dir,rel_canals)),'omitnan'),2);
            cycle_params.GoggleGain_All = rep_gains(contains(rep_dir,rel_canals));
        elseif contains(results.Goggle(1),'ESC') %UPDATE FOR THE NEW ESC VERSION
            if contains(fname,{'LA','LH','LP'})
                rel_side = 'Left';
            else
                rel_side = 'Right';
            end
            if isstruct(CycAvg.Data.AllData)
                cycle_params.GoggleGain_Avg = round(CycAvg.Data.AllData.content_HIT.Result(contains(cellstr(CycAvg.Data.AllData.content_HIT.ResultNames),['Gain',rel_side,'Mean4'])),2); %4 is the average over 0-100ms of the HIT
                cycle_params.GoggleGain_SD = round(CycAvg.Data.AllData.content_HIT.Result(contains(cellstr(CycAvg.Data.AllData.content_HIT.ResultNames),['Gain',rel_side,'Std4'])),2);
            elseif contains(results.Goggle(1),'ESC3')
                cycle_params.GoggleGain_Avg = round(str2double(CycAvg.Data.CSVData{contains(CycAvg.Data.CSVData(:,1),rel_side)&contains(CycAvg.Data.CSVData(:,2),'100 ms'),3}),2); %The average over 0-100ms of the HIT
                cycle_params.GoggleGain_SD = round(str2double(CycAvg.Data.CSVData{contains(CycAvg.Data.CSVData(:,1),rel_side)&contains(CycAvg.Data.CSVData(:,2),'100 ms'),4}),2); 
            end
        end
    case 4
        %% Pulse Train
        % Process Inputs
        stim_vect = results.StimAxis{:};
        %Make sure stim trace is time points long (1 x nt)
        [nc,nt] = size(CycAvg.ll_cyc);
        results.Cycles(:) = nc;
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
        stim_high  = 1:find(diff(Stim_CycAvg)<0, 1, 'first');
        [~,I] = max(max([-stim_vect*[CycAvg.ll_cycavg(stim_high);CycAvg.lr_cycavg(stim_high);CycAvg.lz_cycavg(stim_high)];...
            -stim_vect*[CycAvg.rl_cycavg(stim_high);CycAvg.rr_cycavg(stim_high);CycAvg.rz_cycavg(stim_high)]]));
        for tr = 1:length(traces)
            if isfield(CycAvg,[traces{tr},'_cyc'])&&sum(~isnan(CycAvg.([traces{tr},'_cycavg'])))>1
                max_vel = CycAvg.([traces{tr},'_cyc'])(:,I);
                results.(['MaxVel_',upper(traces{tr})]) = abs(mean(max_vel,'omitnan'));
                results.(['MaxVel_',upper(traces{tr}),'_sd']) = std(max_vel,'omitnan');
                CycAvg.([traces{tr},'_cyc_fit']) = max_vel*ones(1,nt);
                CycAvg.([traces{tr},'_cycavg_fit']) = mean(max_vel,'omitnan')*ones(1,nt);
            end
        end
        %Misalignment
        n_lr = [0,0];
        if isfield(CycAvg,'ll_cyc')
            n_lr(1) = size(CycAvg.ll_cyc,1);
            [results.Align_L(1),results.Align_L_sd(1)] = calc_misalignment(-stim_vect,[CycAvg.ll_cyc(:,I),CycAvg.lr_cyc(:,I),CycAvg.lz_cyc(:,I)]);
        else
            results{1,{'Align_L','Align_L_sd'}} = 0;
        end
        if isfield(CycAvg,'rl_cyc')
            n_lr(2) = size(CycAvg.rl_cyc,1);
            [results.Align_R(1),results.Align_R_sd(1)] = calc_misalignment(-stim_vect,[CycAvg.rl_cyc(:,I),CycAvg.rr_cyc(:,I),CycAvg.rz_cyc(:,I)]);
        else
            results{1,{'Align_R','Align_R_sd'}} = 0;
        end
        %Disconjugacy
        if all(n_lr>0)
            nlr = min(n_lr);
            [results.Disc(1),results.Disc_sd(1)] = calc_misalignment([CycAvg.ll_cyc(1:nlr,I),CycAvg.lr_cyc(1:nlr,I),CycAvg.lz_cyc(1:nlr,I)],[CycAvg.rl_cyc(1:nlr,I),CycAvg.rr_cyc(1:nlr,I),CycAvg.rz_cyc(1:nlr,I)]);
        else
            results{1,{'Disc','Disc_sd'}} = 0;
        end
        cycle_params = [];
end
%Chose the main value for each parameter as the one with the larger eye
%response in the intended canal if applicable.
eyes = {'L','R'};
all_canals = {'LA','LP','LH','RP','RA','RH','+X','-X','+Y','-Y'}; %Canals
var = 'MaxVel';
for i = 1:size(results,1)
    if ismember(results.AxisName(i),all_canals)
        [~,eye] = max(abs([results.([var,'_L',results.AxisLetter{i}])(i),results.([var,'_R',results.AxisLetter{i}])(i)]));
        results.MaxVel(i) = results.(['MaxVel_',eyes{eye},results.AxisLetter{i}])(i);
        results.MaxVel_sd(i) = results.(['MaxVel_',eyes{eye},results.AxisLetter{i},'_sd'])(i);
        results.Phase(i) = results.(['Phase_',eyes{eye},results.AxisLetter{i}])(i);
        results.Phase_sd(i) = results.(['Phase_',eyes{eye},results.AxisLetter{i},'_sd'])(i);
        results.Gain(i) = results.(['Gain_',eyes{eye},results.AxisLetter{i}])(i);
        results.Gain_sd(i) = results.(['Gain_',eyes{eye},results.AxisLetter{i},'_sd'])(i);
        results.Latency(i) = results.(['Latency_',eyes{eye},results.AxisLetter{i}])(i);
        results.Latency_sd(i) = results.(['Latency_',eyes{eye},results.AxisLetter{i},'_sd'])(i);
        results.Tau(i) = results.(['Tau_',eyes{eye},results.AxisLetter{i}])(i);
        results.Tau_sd(i) = results.(['Tau_',eyes{eye},results.AxisLetter{i},'_sd'])(i);
        results.RMSE(i) = results.(['RMSE_',eyes{eye},results.AxisLetter{i}])(i);
        results.RMSE_sd(i) = results.(['RMSE_',eyes{eye},results.AxisLetter{i},'_sd'])(i);
        results.Align(i) = results.(['Align_',eyes{eye}])(i);
        results.Align_sd(i) = results.(['Align_',eyes{eye},'_sd'])(i);
    end
end  
% Make the table and cell
CycAvg.parameterized = results;
CycAvg.cycle_params = cycle_params;
end