%Takes in a cycle average file that should have an exponential output and
%finds the exponential fit for all the relevant sections
function [results,CycAvg] = ParameterizeExponentialFits(fname,CycAvg)
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
%Condition
if contains(fname,{'Light','LIGHT'}) %Always mark the ones in light
    condition = 'Light';
elseif contains(visit,'Visit0') 
    condition = 'Pre-Op';
elseif contains(fname,{'MotionMod','ModOn'})
    condition = 'MotionMod';
elseif contains(fname,{'Baseline','Constant'})
    condition = 'ConstantRate';
elseif contains(fname,'Activation')
    condition = 'Activation';
elseif contains(eeVOR,'eeVOR')
    
elseif contains(fname,{'Dark','DARK','NoStim','NOSTIM'})
    condition = 'NoStim';
else
    condition = '';
end
%% Process Inputs
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
%% Find exponential fits
exp_fit2 = @(tt,p) p(1)*exp(p(2)*tt)+p(3)*exp(p(4)*tt); %2nd order exponential fit same syntax as matlab's built-in
if any(contains(fparts,'RotaryChair')) %Each change is a segment 
    all_traces = [CycAvg.lz_cyc;CycAvg.rz_cyc;...
        CycAvg.ll_cyc;CycAvg.rl_cyc;CycAvg.lr_cyc;CycAvg.rr_cyc;...
        CycAvg.lx_cyc;CycAvg.rx_cyc;CycAvg.ly_cyc;CycAvg.ry_cyc];      
    all_fits = NaN(size(all_traces));
    params = NaN(10,8);
    RMSE = NaN(1,20);
    for i = 1:10
        %Stim is high
        tt = t(Stim==1);
        traces = all_traces(:,Stim==1);
        [fitobj,gof] = fit(tt(~isnan(traces(i,:)))'-tt(1),traces(i,~isnan(traces(i,:)))','exp2');
        RMSE(1,i) = gof.rmse;
        params(i,1:4) = coeffvalues(fitobj);
        all_fits(i,Stim==1) = exp_fit2(tt-tt(1),params(i,1:4));
        %Stim is low
        tt = t(Stim==0);
        traces = all_traces(:,Stim==0);
        [fitobj,gof] = fit(tt(~isnan(traces(i,:)))'-tt(1),traces(i,~isnan(traces(i,:)))','exp2');
        RMSE(1,i+10) = gof.rmse;
        params(i,5:8) = coeffvalues(fitobj);
        all_fits(i,Stim==0) = exp_fit2(tt-tt(1),params(i,5:8));
    end    
    tau_on = -1./params(:,[2,4]);
    tau_on(tau_on < 0) = NaN;
    tau_off = -1./params(:,[6,8]);
    tau_off(tau_off < 0) = NaN;   
    tau = [min(tau_on,[],2)',min(tau_off,[],2)'];
    starts = [all_fits(:,find(Stim==1,1,'first'));all_fits(:,find(Stim==0,1,'first'))]';
    ends = [all_fits(:,find(Stim==1,1,'last'));all_fits(:,find(Stim==0,1,'last'))]';
    all_params = reshape([tau;starts;ends;RMSE],[],1)';
    labs = strcat(repmat(reshape(repmat({'L_','R_'},4,1),[],1)',1,10),repmat(reshape(repmat({'LHRH_','LARP_','RALP_','X_','Y_'},8,1),[],1)',1,2),reshape(repmat({'Motion_','NoMotion_'},40,1),[],1)',repmat({'Tau','Start','End','RMSE'},1,20));
    % Compile into a table
    tab1 = cell2table({fname,subject,visit,date,condition});
    tab1.Properties.VariableNames = {'File','Subject','Visit','Date','Condition'};
    tab2 = array2table([max_vel,all_params]);
    tab2.Properties.VariableNames = [{'Amplitude'},labs];
    results = [tab1,tab2];    
    %Save in the struct
    CycAvg.parameterized = results;
    CycAvg.lz_cycavg_fit = all_fits(1,:);
    CycAvg.rz_cycavg_fit = all_fits(2,:);
    CycAvg.ll_cycavg_fit = all_fits(3,:);
    CycAvg.rl_cycavg_fit = all_fits(4,:);
    CycAvg.lr_cycavg_fit = all_fits(5,:);
    CycAvg.rr_cycavg_fit = all_fits(6,:);
    CycAvg.lx_cycavg_fit = all_fits(7,:);
    CycAvg.rx_cycavg_fit = all_fits(8,:);
    CycAvg.ly_cycavg_fit = all_fits(9,:);
    CycAvg.ry_cycavg_fit = all_fits(10,:);   
elseif contains(fparts('Activation')) %Each trigger is light vs. dark
    %Add code here
    
    
else
    error('Unrecognized experiment type based on file name')
end
end