function filt_params = VOGA__saveLastUsedParams(filt_params)
if nargin > 0 && ~isempty(filt_params) %Valid Filter Parameters passed through
    save([userpath,filesep,'VOGA_LastUsedFilterParams.mat'],'filt_params')
elseif ~isfile([userpath,filesep,'VOGA_LastUsedFilterParams.mat']) %Does not exist yet
    %Default table, no presets
    trace_names = {'LX','RX','LY','RY','LZ','RZ','LLARP','RLARP','LRALP','RRALP','ALL'};
    pos_filts = {'median','mean','sgolay','spline'};
    vel_filts = {'median','sgolay','accel','irlssmooth','spline'};   
    filt_params.filt1.pos = array2table(NaN(length(trace_names),length(pos_filts)),'VariableNames',pos_filts,'RowNames',trace_names);
    filt_params.filt1.vel = array2table(NaN(length(trace_names),length(vel_filts)),'VariableNames',vel_filts,'RowNames',trace_names);
    %Make filt_params struct
    filt_params.YLim.Pos = [-30,30];
    filt_params.YLim.Vel = [-100,100];  
    save([userpath,filesep,'VOGA_LastUsedFilterParams.mat'],'filt_params')
else
    load([userpath,filesep,'VOGA_LastUsedFilterParams.mat'],'filt_params')
end
end