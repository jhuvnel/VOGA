%% Filter position and velocity
function [filt,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc] = MakeCycAvg__filterTraces(filt,keep_inds,te,ts,t_snip,stim,stims,Data,t_interp)
    %% Validate Parameters
    %Set the all values
    filt.pos(:,~isnan(filt.pos{'ALL',:})) = repmat(filt.pos(end,~isnan(filt.pos{'ALL',:})),11,1);
    filt.vel(:,~isnan(filt.vel{'ALL',:})) = repmat(filt.vel(end,~isnan(filt.vel{'ALL',:})),11,1);       
    %Position
    filt.pos.median = MakeCycAvg__filterParameterCheck(filt.pos.median,'med');
    filt.pos.spline = MakeCycAvg__filterParameterCheck(filt.pos.spline,'spline');
    sgolay_f = MakeCycAvg__filterParameterCheck([filt.pos.sgolay1,filt.pos.sgolay2],'sgolay');
    filt.pos.sgolay1 = sgolay_f(:,1);
    filt.pos.sgolay2 = sgolay_f(:,2);
    %Velocity
    filt.vel.median = MakeCycAvg__filterParameterCheck(filt.vel.median,'med');
    filt.vel.spline = MakeCycAvg__filterParameterCheck(filt.vel.spline,'spline');
    filt.vel.irlssmooth = MakeCycAvg__filterParameterCheck(filt.vel.irlssmooth,'smooth');
    sgolay_f = MakeCycAvg__filterParameterCheck([filt.vel.sgolay1,filt.vel.sgolay2],'sgolay');
    filt.vel.sgolay1 = sgolay_f(:,1);
    filt.vel.sgolay2 = sgolay_f(:,2);
    filt.vel.accel = MakeCycAvg__filterParameterCheck(filt.vel.accel,'accel');
    %% Filter Position
    if contains(Data.info.goggle_ver,'GNO') %No position traces
        Data_pos = [];
        Data_pos_filt = [];
        %Make Data_vel struct        
        Data_vel.t = ts;
        Data_vel.stim = stim;
        Data_vel.RE_Vel_Z = -Data.RE_Vel_Z;  
        Data_vel.stim_cond = stim(reshape(keep_inds,[],1));
        t_cond = median(diff(ts))*(0:length(Data_vel.stim_cond)-1)';
        Data_vel.t_cond = t_cond;
        Data_vel.RE_Vel_Z_cond = -Data.RE_Vel_Z(reshape(keep_inds,[],1)); 
        if contains(Data.info.dataType,{'RA','LP'})
            Data_vel.RE_Vel_RALP = -Data.RE_Vel_Y; 
            Data_vel.RE_Vel_RALP_cond = -Data.RE_Vel_Y(reshape(keep_inds,[],1)); 
        elseif contains(Data.info.dataType,{'LA','RP'})
            Data_vel.RE_Vel_LARP = Data.RE_Vel_Y;
            Data_vel.RE_Vel_LARP_cond = Data.RE_Vel_Y(reshape(keep_inds,[],1)); 
        else
            Data_vel.RE_Vel_Y = -Data.RE_Vel_Y;
            Data_vel.RE_Vel_Y_cond = -Data.RE_Vel_Y(reshape(keep_inds,[],1)); 
        end
    else
        Data_pos = Data;
        Data_pos.te = te;
        Data_pos.ts = ts;
        Data_pos.stim = stim;
        [~,Data_pos] = angpos2angvel(Data_pos); %Just to get LARP/RALP transformation
        traces = filt.pos.Properties.RowNames;    
        for i = 1:length(traces)
            var_n = [traces{i}(1),'E_Position_',traces{i}(2:end)];
            if isfield(Data,var_n)
                Data_pos_filt.(var_n) = spline_filt(te,sgolay_filt(med_filt(Data.(var_n),filt.pos.median(traces{i})),[filt.pos.sgolay1(traces{i}),filt.pos.sgolay2(traces{i})]),ts,filt.pos.spline(traces{i}));
            end
        end
        Data_pos_filt.t = ts;
        Data_pos_filt.Fs = Data.Fs;
        Data_pos_filt.stim = stim;
        [Data_vel,Data_pos_filt] = angpos2angvel(Data_pos_filt); %Only uses X, Y, Z data right now
        Data_vel.t = ts;
        Data_vel.stim = stim;
    end
    %% Filter Velocity
    traces = filt.vel.Properties.RowNames;   
    Data_vel_filt.t = ts;
    Data_vel_filt.stim = stim;
    if contains(Data.info.goggle_ver,'GNO')
        Data_vel_filt.t_cond = Data_vel.t_cond;
        Data_vel_filt.stim_cond = stim(reshape(keep_inds,[],1));
    end
    Data_cyc.t = t_snip;
    Data_cyc.stim = stims;      
    for i = 1:length(traces)            
        canal = traces{i}(2:end); 
        var_n = [traces{i}(1),'E_Vel_',canal];
        if isfield(Data_vel,var_n)
            if contains(Data.info.goggle_ver,'GNO')
                rel_t = t_cond;                
                long_i = reshape(1:(size(keep_inds,1)*size(keep_inds,2)),size(keep_inds,1),[]);
                rel_trace = Data_vel.([var_n,'_cond']);
                temp1 = med_filt(rel_trace,filt.vel.median(traces{i}));
                temp2 = sgolay_filt(temp1,[filt.vel.sgolay1(traces{i}),filt.vel.sgolay2(traces{i})]);
                temp3 = accel_QPR(rel_t,temp2,filt.vel.accel(traces{i}));
                temp4 = irls_filt(temp3,filt.vel.irlssmooth(traces{i}));   
                temp5 = temp4;
                temp5(long_i(t_interp,:)) = NaN; %Now remove trace if needed
                temp5 = interp1(rel_t(~isnan(temp5)),temp5(~isnan(temp5)),rel_t);
                temp6 = spline_filt(rel_t,temp5,rel_t,filt.vel.spline(traces{i}));
                Data_vel_filt.(var_n) = temp6;
                Data_cyc.(var_n) = reshape(Data_vel_filt.(var_n),size(keep_inds,1),[]);
            else
                rel_t = ts;
                rel_trace = Data_vel.(var_n);
                temp1 = med_filt(rel_trace,filt.vel.median(traces{i}));
                temp2 = sgolay_filt(temp1,[filt.vel.sgolay1(traces{i}),filt.vel.sgolay2(traces{i})]);
                temp3 = accel_QPR(rel_t,temp2,filt.vel.accel(traces{i}));
                temp4 = irls_filt(temp3,filt.vel.irlssmooth(traces{i}));
                temp5 = temp4;
                temp5(keep_inds(t_interp,:)) = NaN; %Now remove trace if needed
                temp5 = interp1(rel_t(~isnan(temp5)),temp5(~isnan(temp5)),rel_t);
                temp6 = spline_filt(rel_t,temp5,rel_t,filt.vel.spline(traces{i}));
                Data_vel_filt.(var_n) = temp6;
                Data_cyc.(var_n) = Data_vel_filt.(var_n)(keep_inds);
            end
        end
    end  
    Data_cyc.keep_inds = keep_inds;        
    %Clear the 'ALL' value on the filter
    filt.pos{end,:} = NaN; 
    filt.vel{end,:} = NaN;    
end