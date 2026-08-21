% exponentialFitTesting.m
% Tests an alternative exponential fit function on data loaded from a
% velstep CycAvg .mat file.
% Currently this script does not save any data.

% close all
load('VNELcolors.mat','colors')
files = dir('Cycle Averages\*.mat');
filenames = {files(contains({files.name}, 'VelStep')).name};

[indx_exp, tf] = listdlg('ListString', filenames, 'PromptString', 'Choose velstep file to test exponential fits for.', 'SelectionMode', 'single', 'ListSize', [600 300]);

if tf
    figName = filenames{indx_exp};
    figName = figName(1:end-4);
    load(['Cycle Averages\',figName,'.mat']) % open cycle average .mat file
    amp = regexp(figName,'\d+dps','match');
    max_stim = str2double(amp{1}(1:end-3));

    [CycAvgOut, type] = ParameterizeCycAvg(CycAvg);
    t = CycAvgOut.t;
    sz = size(t);
    if sz(2) > 10
        t = t';
    end


    Stim_start_idx = find(abs(CycAvg.stim)>0.9*abs(max_stim),1,'first');
    Stim_start = CycAvg.t(Stim_start_idx); % time the "stim" started (as in velocity reached >90% of programmed max)

    % Stim_start_idx = find(CycAvg.t>=13,1,'first'); % cut off first 5 seconds
    % Stim_start = CycAvg.t(Stim_start_idx);
    Stim_end_idx = find(abs(CycAvg.stim(CycAvg.t>Stim_start))<0.9*abs(max_stim),1,'first');
    Stim_end = CycAvg.t(Stim_end_idx); % time the "stim" end (as in velocity reached <90% of programmed max)

    t_start = t(Stim_start_idx);
    ind_high = Stim_start_idx:Stim_end_idx;
    t_high = t(ind_high);
    t_high = t_high - t_start;



    % plotSegment(CycAvgOut.Data)
    % hold on
    % cycle_params.exp_ord1_const_params
    % plotCycAvg(CycAvgOut, 1)
    
    fig = figure('Units','inches','Position',[0 0 7 10]);
    sz = size(CycAvgOut.Data_filtvel.LE_Vel_LARP);
    if sz(2) > 10
        LE_Vel_LARP = CycAvgOut.Data_filtvel.LE_Vel_LARP(1,:)';
        LE_Vel_RALP = CycAvgOut.Data_filtvel.LE_Vel_RALP(1,:)';
        LE_Vel_Z = CycAvgOut.Data_filtvel.LE_Vel_Z(1,:)';
        RE_Vel_LARP = CycAvgOut.Data_filtvel.RE_Vel_LARP(1,:)';
        RE_Vel_RALP = CycAvgOut.Data_filtvel.RE_Vel_RALP(1,:)';
        RE_Vel_Z = CycAvgOut.Data_filtvel.RE_Vel_Z(1,:)';
    else
        LE_Vel_LARP = CycAvgOut.Data_filtvel.LE_Vel_LARP(:,1);
        LE_Vel_RALP = CycAvgOut.Data_filtvel.LE_Vel_RALP(:,1);
        LE_Vel_Z = CycAvgOut.Data_filtvel.LE_Vel_Z(:,1);
        RE_Vel_LARP = CycAvgOut.Data_filtvel.RE_Vel_LARP(:,1);
        RE_Vel_RALP = CycAvgOut.Data_filtvel.RE_Vel_RALP(:,1);
        RE_Vel_Z = CycAvgOut.Data_filtvel.RE_Vel_Z(:,1);
    end
    plot(t,LE_Vel_LARP,'Color',colors.l_l)
    hold on
    plot(t,LE_Vel_RALP,'Color',colors.l_r)
    plot(t,LE_Vel_Z,'Color',colors.l_z)
    plot(t,RE_Vel_LARP,'Color',colors.r_l)
    plot(t,RE_Vel_RALP,'Color',colors.r_r)
    plot(t,RE_Vel_Z,'Color',colors.r_z)
    plot(t,CycAvgOut.stim(:,1),'-k')


    % plotInds = ~isnan(CycAvgOut.ll_cyc_fit(:,1));
    plot(t,CycAvgOut.ll_cyc_fit(1,:)','--','Color',colors.l_l,'LineWidth',1)
    plot(t,CycAvgOut.lr_cyc_fit(1,:)','--','Color',colors.l_r,'LineWidth',1)
    plot(t,CycAvgOut.lz_cyc_fit(1,:)','--','Color',colors.l_z,'LineWidth',1)
    plot(t,CycAvgOut.rl_cyc_fit(1,:)','--','Color',colors.r_l,'LineWidth',1)
    plot(t,CycAvgOut.rr_cyc_fit(1,:)','--','Color',colors.r_r,'LineWidth',1)
    plot(t,CycAvgOut.rz_cyc_fit(1,:)','--','Color',colors.r_z,'LineWidth',1)


    dat_high = RE_Vel_Z(ind_high);
    sz = size(dat_high);

    if  sz(2) > 10
        dat_high = dat_high';
    end
    % find expected eye vel direction
    if CycAvgOut.info.stim_axis(3) == 1 % positive LHRH motion
        direction = -1;
    elseif CycAvgOut.info.stim_axis(3) == -1
        direction = 1;
    else
        direction = 0;
    end
    % [expc_params, yfit, gof, exp1_c_fitobj] = expcfit(t_high, dat_high, direction);
    % exp1_c_fitobj
    % plot(t_high + t_start, yfit, '-.','Color','k','LineWidth',2)

    % axis([t(1) t(end) min() 20])
    title('Filtered Angular Velocity')
    xlabel('Time (s)')
    ylabel('Velocity (deg/s)')
    annotation('textbox', [0.1 0.9 0.5 0.1], 'String', CycAvg.info.name, 'FitBoxToText', 'on', 'LineStyle', 'none')
    % legend('LE LARP','LE RALP','LE Z','RE LARP','RE RALP','RE Z','Stim', 'LE LARP Fit', 'LE RALP Fit', 'LE Z Fit', 'RE LARP Fit', 'RE RALP Fit', 'RE Z Fit')
    

end

function [params, yfit, gof, exp1_c_fitobj] = expcfit(x, y, direction)
% Fits to the function y = k.*exp(-x./tau) + c.
% Parameter bounds and initial conditions are based on expected velocity
% step experiment results. direction argument specifies which direction you
% expect the eye velocity response to be in, i.e. if you have an LH motion,
% which is a positive LHRH head motion, the expected LHRH eye movement
% would be negative, and the direction argument would be -1. If you expect
% a positive slow-phase eye velocity, input 1. If 0 or no direction is
% input for the direction, then k and c are not constrained to be positive
% or negative.
% Returns the parameters ([k, tau, c]), fitted curve over inputs x,
% goodness of fit metrics, and the fit object.
    
    if nargin < 3
        direction = 0;
    end
    exp1_c_fit = @(k,tau,c,x) k.*exp(-x./tau) + c;
    
    % lower and upper bounds for parameters
    if direction == 1
        pl = [0 0.0001 0];
        pu = [300 1000 300];
    elseif direction == -1
        pl = [-300 0.0001 -300];
        pu = [0 1000 0];
    else % if no direction provided, do not limit bounds for parameters to be positive or negative
        pl = [-300 0.0001 -300];
        pu = [300 1000 300];
    end
    
    % initial conditions set to:
        % [k0 (exponential constant) = median of first 10th of trace - c0, 
        % tau0 (time constant) = 10, 
        % c0 (plateau) = median of last 4th of trace]
     p0 = [median(y(1:round(length(y)/10))) - median(y(end-round(length(y)/4):end)),...
        10,median(y(end-round(length(y)/4):end))];
     if sign(direction) ~= sign(p0(1))
        p0(1) = -1*p0(1);
     end

     [exp1_c_fitobj, gof] = fit(x,y, exp1_c_fit, 'Lower', pl, 'Upper', pu, 'StartPoint', p0);
    
    yfit = exp1_c_fitobj(x);
    params = nan(3,1);
    params(1) = exp1_c_fitobj.k;
    params(2) = exp1_c_fitobj.tau;
    params(3) = exp1_c_fitobj.c;

end