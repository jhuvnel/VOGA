% manualExponentialFit.m
% Do exponential fit by manually selecting points on velstep cycle average
% plot.
% This will change the Gain and Tau fields in the cycle average parameters
% (summary table) for that cycle average file. You must "Parametrize and
% Plot" in VOGA afterwards to update the summary table. Only updates those
% two fields, no others (such as the max vel and tau for each trace), so
% the generated plot will not have the new fit that you generated.
% close all
files = dir('Cycle Averages\*.mat');
filenames = {files(contains({files.name}, 'VelStep')).name};

[indx_exp, tf] = listdlg('ListString', filenames, 'PromptString', 'Choose velstep file to do manual exponential fit for.', 'SelectionMode', 'single', 'ListSize', [600 300]);

if tf
    figName = filenames{indx_exp};
    figName = figName(1:end-4);
    load(['Cycle Averages\',figName,'.mat']) % open cycle average .mat file
    
    fig = openfig(['Cycle Averages\',figName,'.fig']);
    ax = findall(fig, 'Type', 'axes');
    vel_axis = ax(4);
    axes(vel_axis)
    
    
    amp = regexp(figName,'\d+dps','match');
    max_stim = str2double(amp{1}(1:end-3));
    Stim_start_idx = find(abs(CycAvg.stim)>0.9*abs(max_stim),1,'first');
    Stim_start = CycAvg.t(Stim_start_idx); % time the "stim" started (as in velocity reached >90% of programmed max)
    Stim_end_idx = find(abs(CycAvg.stim(CycAvg.t>Stim_start))<0.9*abs(max_stim),1,'first');
    Stim_end = CycAvg.t(Stim_end_idx); % time the "stim" end (as in velocity reached <90% of programmed max)
    
    mult = [-1;1];
    if contains(CycAvg.parameterized.AxisName(1),{'RH','LA','LP'}) %Negative directions
        mult = [1;-1];
    end   
    
    disp('Select any number of points (from left to right) for the true slow phase velocity response in order to make an exponential fit.')
    [x, y] = ginput;
    y(x<=Stim_start) = [];
    x(x<=Stim_start) = [];
    y(x>=Stim_end) = [];
    x(x>=Stim_end) = [];
    x = x - Stim_start; % set t = 0 to be start of stim
    
    options = optimset('Display','off'); %suppress output from fminsearch
    %1order exponential fit with decay to 0
    exp_fit1 = @(tt,p) p(1).*exp(-tt./p(2));
    makefit2 = @(trace,T) fminsearch(@(p) sum((exp_fit1(T,p)-trace).^2,'omitnan')+sum(p.^2),[median(trace(1));1],options);

    %1order exponential fit with decay not forced to 0.
    exp_fit2 = @(tt,p) p(1).*exp(-tt./p(2)) + p(3);
    % makefit3 = @(trace,T) fminsearch(@(p) sum((exp_fit2(T,p)-trace).^2,'omitnan')+sum(p.^2),[median(trace);1;0],options);
    % makefit3 = @(trace,T) fminsearch(@(p) sum((exp_fit2(T,p)-trace).^2,'omitnan'),[max(abs(trace))*mult(1);10;median(trace)],options);
    % makefit3 = @(trace,T) fminsearch(@(p) sum((exp_fit2(T,p)-trace).^2 + sum(p),'omitnan'),[median(trace(1:round(length(trace)/10)));10;median(trace(end-round(length(trace)/4):end))],options);
    makefit3 = @(trace,T) fminsearch(@(p) sum((exp_fit2(T,p)-trace).^2 + sum(abs(p)),'omitnan'),...
        [median(trace(1:round(length(trace)/10))) - median(trace(end-round(length(trace)/4):end));...
        10;median(trace(end-round(length(trace)/4):end))],options);

    params2 = makefit2(y, x);

    params3 = makefit3(y, x);
    
    % tau = params3(2);
    % max_vel = params3(1) + params3(3);

    tau = params2(2);
    max_vel = params2(1);
    
    % save tau and gain in CycAvg struct
    CycAvg.parameterized.Tau(1) = tau;
    CycAvg.parameterized.Gain(1) = max_vel./max_stim*mult(1);

    hold on
    plot(x+Stim_start,y,'k.','MarkerSize',10)
    t_plot = CycAvg.t(Stim_start_idx:Stim_end_idx) - Stim_start;
    yfit_plot2 = exp_fit1(t_plot,params2);
    plot(Stim_start + t_plot,yfit_plot2,'k','LineWidth',1.5)
    % yfit_plot3 = exp_fit2(t_plot,params3);
    % plot(Stim_start + t_plot,yfit_plot3,'k--','LineWidth',1.5)
    
    archiveFolder = 'Cycle Averages\Archive';
    if ~isfolder(archiveFolder)
        mkdir(archiveFolder)
    end
    movefile(['Cycle Averages\',filenames{indx_exp}],archiveFolder)
    movefile(['Cycle Averages\',figName,'.fig'],archiveFolder)
    % re-save CycAvg struct
    save(['Cycle Averages\',figName,'.mat'], 'CycAvg')
    saveas(fig, ['Cycle Averages\',figName],'fig')
end