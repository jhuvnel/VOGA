function filt = MakeCycAvg__autoFilter(Data,filt,plot_info)
type = Data.info.type;
%More filtering
traces_pos = plot_info.traces_pos;     
plot_info.traces_vel = plot_info.traces_pos;
%Reset filters
filt.pos{:,:} = NaN; 
filt.vel{:,:} = NaN;

% CycAvg = MakeCycAvg__filterTraces(Data,filt);
% ha = MakeCycAvg__plotFullCycAvg([],CycAvg,plot_info); 
% 
% t = Data.Time_Eye;
% LZ = Data.LE_Position_Z;
% RZ = Data.RE_Position_Z;
% LX = Data.LE_Position_X;
% RX = Data.RE_Position_X;
% LY = Data.LE_Position_Y;
% RY = Data.RE_Position_Y;
% 
% 
% noise = @(trace) median(abs(gradient(gradient(trace,median(diff(t))),median(diff(t)))),2);
% rmse = @(fit_trace,og_trace) sqrt(mean((fit_trace-og_trace).^2));
% disp(noise([LX,LY,LZ,RX,RY,RZ]'))
% %Median filter is always good for dropouts
% filt.pos.median(traces_pos) = 3; 
% plot_info.traces_pos = {'LZ'};
% plot_info.traces_vel = {'LZ'};

% median_cyc_vel = NaN(length(traces_vel),length(CycAvg.t));
% for i = 1:length(traces_vel)
%    median_cyc_vel(i,:) = median(CycAvg.([lower(traces_vel{i}(1:2)),'_cyc']),'omitnan'); 
% end
% pred_max = prctile(max(abs(median_cyc_vel),[],'omitnan'),75)/.75;
% int_vel = cumsum(median_cyc_vel);
filt.pos.median(end) = 3;
filt.pos.spline({'LX','RX'}) = 1-5e-5;
if type == 1
    filt.vel.irlssmooth(end) = round(length(Data.t_snip)*0.16); %heuristic
elseif type == 2    
    filt.vel.irlssmooth(end) = 200; %heuristic
end
end