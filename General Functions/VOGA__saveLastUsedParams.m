function VOGA__saveLastUsedParams(filt_params)
if ~any(contains(extractfield(dir(userpath),'name'),'VOGA_LastUsedFilterParams.mat')) %first time running
    %Default table, no presets
    trace_names = {'LX','RX','LY','RY','LZ','RZ','LLARP','RLARP','LRALP','RRALP','ALL'};
    pos_filts = {'median','spline','sgolay1','sgolay2'}; %add more as needed
    vel_filts = {'accel','median','spline','sgolay1','sgolay2','irlssmooth','slip'}; %add more as needed
    filt1.pos = array2table(NaN(length(trace_names),length(pos_filts)));
    filt1.pos.Properties.VariableNames = pos_filts;
    filt1.pos.Properties.RowNames = trace_names;
    filt1.vel = array2table(NaN(length(trace_names),length(vel_filts)));
    filt1.vel.Properties.VariableNames = vel_filts;
    filt1.vel.Properties.RowNames = trace_names;
    %Make filt_params struct
    filt_params.filt1 = filt1;
    filt_params.YLim.Pos = [-100,100];
    filt_params.YLim.Vel = [-100,100];  
end
save([userpath,filesep,'VOGA_LastUsedFilterParams.mat'],'filt_params') 
end