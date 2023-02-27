%% Filter position and velocity
function [filt,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc] = MakeCycAvg__filterTraces(filt,keep_inds,te,ts,t_snip,stim,stims,Data,t_interp)
    %% Validate Parameters
    %Make sure all expected table columns are there
    pos_labs = {'median','spline','sgolay1','sgolay2','lowpass'};
    vel_labs = {'median','spline','irlssmooth','sgolay1','sgolay2','accel'};    
    missing_pos_labs = pos_labs(~ismember(pos_labs,filt.pos.Properties.VariableNames));
    missing_vel_labs = vel_labs(~ismember(vel_labs,filt.vel.Properties.VariableNames));
    for i = 1:length(missing_pos_labs)
        filt.pos.(missing_pos_labs{i}) = NaN(size(filt.pos,1),1);
    end
    for i = 1:length(missing_vel_labs)
        filt.vel.(missing_vel_labs{i}) = NaN(size(filt.vel,1),1);
    end
    %Set the "ALL" values
    filt.pos(:,~isnan(filt.pos{'ALL',:})) = repmat(filt.pos(end,~isnan(filt.pos{'ALL',:})),11,1);
    filt.vel(:,~isnan(filt.vel{'ALL',:})) = repmat(filt.vel(end,~isnan(filt.vel{'ALL',:})),11,1);       
    %Position
    filt.pos.median = MakeCycAvg__filterParameterCheck(filt.pos.median,'med');
    filt.pos.spline = MakeCycAvg__filterParameterCheck(filt.pos.spline,'spline');
    sgolay_f = MakeCycAvg__filterParameterCheck([filt.pos.sgolay1,filt.pos.sgolay2],'sgolay');
    filt.pos.sgolay1 = sgolay_f(:,1);
    filt.pos.sgolay2 = sgolay_f(:,2);
    filt.pos.lowpass = MakeCycAvg__filterParameterCheck(filt.pos.lowpass,'lpass');
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
        if contains(Data.info.dataType,{'RA','LP'})
            Data_vel.RE_Vel_RALP = -Data.RE_Vel_Y; 
        elseif contains(Data.info.dataType,{'LA','RP'})
            Data_vel.RE_Vel_LARP = Data.RE_Vel_Y;
        else
            Data_vel.RE_Vel_Y = -Data.RE_Vel_Y;
        end
    elseif contains(Data.info.goggle_ver,'ESC')&&~contains(Data.info.goggle_ver,'ESC3') %No position traces
        Data_pos = [];
        Data_pos_filt = [];
        %Make Data_vel struct        
        Data_vel.t = ts;
        Data_vel.stim = stim;
        Data_vel.LE_Vel_Z = Data.LE_Vel_Z;  
        if contains(Data.info.dataType,{'RA','LP'})
            Data_vel.LE_Vel_RALP = Data.LE_Vel_Y; 
        elseif contains(Data.info.dataType,{'LA','RP'})
            Data_vel.LE_Vel_LARP = -Data.LE_Vel_Y;
        else
            Data_vel.LE_Vel_Y = Data.LE_Vel_Y;
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
                rel_trace = Data_pos.(var_n);
                temp0 = interp1(te(~isnan(rel_trace)),rel_trace(~isnan(rel_trace)),te);
                temp1 = filterTrace('median',temp0,filt.pos.median(i));
                temp2 = filterTrace('lowpass',temp1,filt.pos.lowpass(i),te);
                temp3 = filterTrace('sgolay',temp2,[filt.pos.sgolay1(i),filt.pos.sgolay2(i)]);
                temp4 = filterTrace('spline',temp3,filt.pos.spline(i),te,ts);
                Data_pos_filt.(var_n) = temp4;
                if contains(Data.info.goggle_ver,'ESC3')
                    Data_pos.([var_n,'_cond']) = Data_pos.(var_n)(reshape(keep_inds,[],1)); 
                    Data_pos_filt.([var_n,'_cond']) = Data_pos_filt.(var_n)(reshape(keep_inds,[],1)); 
                end
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
    if contains(Data.info.name,'Impulse')
        Data_vel.stim_cond = stim(reshape(keep_inds,[],1));
        t_cond = median(diff(ts))*(0:length(Data_vel.stim_cond)-1)';
        Data_vel.t_cond = t_cond;
        Data_vel_filt.t_cond = t_cond;
        Data_vel_filt.stim_cond = stim(reshape(keep_inds,[],1));
    end
    Data_cyc.t = t_snip;
    Data_cyc.stim = stims;  
    for i = 1:length(traces)            
        var_n = [traces{i}(1),'E_Vel_',traces{i}(2:end)];
        if isfield(Data_vel,var_n)
            if contains(Data.info.name,'Impulse')
                rel_t = t_cond;                
                long_i = reshape(1:(size(keep_inds,1)*size(keep_inds,2)),size(keep_inds,1),[]);
                rm_ind = long_i(t_interp,:);
                Data_vel.([var_n,'_cond']) = Data_vel.(var_n)(reshape(keep_inds,[],1)); 
                rel_trace = Data_vel.([var_n,'_cond']);
            else
                rel_t = ts;
                rm_ind = keep_inds(t_interp,:);
                rel_trace = Data_vel.(var_n);
            end
            temp1 = filterTrace('median',rel_trace,filt.vel.median(i));
            temp2 = filterTrace('sgolay',temp1,[filt.vel.sgolay1(i),filt.vel.sgolay2(i)]);
            temp3 = filterTrace('accel',temp2,filt.vel.accel(i),rel_t);
            temp4 = filterTrace('manual_interp',temp3,[],rm_ind,rel_t);
            temp5 = filterTrace('irls',temp4,filt.vel.irlssmooth(i));
            temp6 = filterTrace('spline',temp5,filt.vel.spline(i),rel_t,rel_t);
            Data_vel_filt.(var_n) = temp6;
            if contains(Data.info.name,'Impulse')
                Data_cyc.(var_n) = reshape(Data_vel_filt.(var_n),size(keep_inds,1),[]);
                Data_cyc.([var_n,'_smooth']) = reshape(temp2,size(keep_inds,1),[]);
                warning('off')
                [~,peak_ind] = findpeaks(abs(temp2),'MinPeakDistance',0.05/median(diff(rel_t)),...
                    'MinPeakWidth',0.02/median(diff(rel_t)),...
                    'MinPeakProminence',50,'MinPeakHeight',75);
                warning('on')
                is_sacc = 0*temp2;
                is_sacc(peak_ind) = 1;
                Data_cyc.([var_n,'_saccade']) = reshape(is_sacc,size(keep_inds,1),[]);
            else
                Data_cyc.(var_n) = Data_vel_filt.(var_n)(keep_inds);
            end
        end
    end  
    Data_cyc.keep_inds = keep_inds;
    %Clear the 'ALL' value on the filter
    filt.pos{end,:} = NaN; 
    filt.vel{end,:} = NaN;    
end