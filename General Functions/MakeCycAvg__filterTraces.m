%% Filter position and velocity
function [filt,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc] = MakeCycAvg__filterTraces(filt,keep_inds,te,ts,t_snip,stim,stims,Data)
    %% Validate Parameters
    %Position
    filt.pos.median = MakeCycAvg__filterParameterCheck(filt.pos.median,'med');
    filt.pos.spline = MakeCycAvg__filterParameterCheck(filt.pos.spline,'spline');
    sgolay_f = MakeCycAvg__filterParameterCheck([filt.pos.sgolay1,filt.pos.sgolay2],'sgolay');
    filt.pos.sgolay1 = sgolay_f(:,1);
    filt.pos.sgolay2 = sgolay_f(:,2);
    %Set the all value
    filt.pos(:,~isnan(filt.pos{'ALL',:})) = repmat(filt.pos(end,~isnan(filt.pos{'ALL',:})),11,1);
    %Velocity
    filt.vel.median = MakeCycAvg__filterParameterCheck(filt.vel.median,'med');
    filt.vel.spline = MakeCycAvg__filterParameterCheck(filt.vel.spline,'spline');
    filt.vel.irlssmooth = MakeCycAvg__filterParameterCheck(filt.vel.irlssmooth,'smooth');
    sgolay_f = MakeCycAvg__filterParameterCheck([filt.vel.sgolay1,filt.vel.sgolay2],'sgolay');
    filt.vel.sgolay1 = sgolay_f(:,1);
    filt.vel.sgolay2 = sgolay_f(:,2);
    %Set the all value
    filt.vel(:,~isnan(filt.vel{'ALL',:})) = repmat(filt.vel(end,~isnan(filt.vel{'ALL',:})),11,1);       
    %% Filter Position
    if contains(Data.info.goggle_ver,'GNO') %No position traces
        Data_pos = [];
        Data_pos_filt = [];
        %Make Data_vel struct
        Data_vel.t = ts;
        Data_vel.stim = stim;
        Data_vel.RE_Vel_Z = -Data.RE_Vel_Z;  
        if contains(Data.info.dataType,{'RA','LP'})
            Data_vel.RE_Vel_RALP = -Data.RE_Vel_Y;            
        elseif contains(Data.info.dataType,{'LA','RP'})
            Data_vel.RE_Vel_LARP = Data.RE_Vel_Y;
        else
            Data_vel.RE_Vel_Y = -Data.RE_Vel_Y;
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
    Data_cyc.t = t_snip;
    Data_cyc.stim = stims;   
    for i = 1:length(traces)            
        canal = traces{i}(2:end); 
        var_n = [traces{i}(1),'E_Vel_',canal];
        if isfield(Data_vel,var_n)
            if contains(Data.info.goggle_ver,'GNO') %Save a desaccaded version
                Data_vel_filt.(var_n) = spline_filt(ts,sgolay_filt(med_filt(Data_vel.(var_n),filt.vel.median(traces{i})),[filt.vel.sgolay1(traces{i}),filt.vel.sgolay2(traces{i})]),ts,filt.vel.spline(traces{i}));
                if filt.vel.irlssmooth(traces{i})==1
                    Data_vel_filt.([var_n,'_QPR']) = Data_vel_filt.(var_n);
                else    
                    Data_vel_filt.([var_n,'_QPR']) = irls_filt(Data_vel.(var_n),filt.vel.irlssmooth(traces{i}));  
                end
                Data_cyc.(var_n) = Data_vel_filt.(var_n)(keep_inds);
                Data_cyc.([var_n,'_QPR']) = Data_vel_filt.([var_n,'_QPR'])(keep_inds);
            else
                Data_vel_filt.(var_n) = spline_filt(ts,sgolay_filt(med_filt(Data_vel.(var_n),filt.vel.median(traces{i})),[filt.vel.sgolay1(traces{i}),filt.vel.sgolay2(traces{i})]),ts,filt.vel.spline(traces{i}));
                Data_cyc.(var_n) = Data_vel_filt.(var_n)(keep_inds);
            end
        end
    end  
    Data_cyc.keep_inds = keep_inds;        
    %Clear the 'ALL' value on the filter
    filt.pos{end,:} = NaN; 
    filt.vel{end,:} = NaN;    
end