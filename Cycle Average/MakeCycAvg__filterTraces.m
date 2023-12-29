%% Filter position and velocity
function [CycAvg,filt] = MakeCycAvg__filterTraces(Data,filt,CycAvg)
if nargin < 3
    %% Feed in variables from "Data" struct
    te = Data.te;
    ts = Data.ts;
    t_snip = Data.t_snip;
    stim = Data.stim;
    stims = Data.stims;
    fname = Data.info.name;
    keep_inds = Data.keep_inds;
    keep_tr = filt.keep_tr;
    %% Validate Parameters
    %Make sure all expected table columns are there
    pos_labs = {'tvd','median','sgolay','spline'};
    vel_labs = {'tvd','sgolay','irlssmooth','spline'};    
    missing_pos_labs = pos_labs(~ismember(pos_labs,filt.pos.Properties.VariableNames));
    missing_vel_labs = vel_labs(~ismember(vel_labs,filt.vel.Properties.VariableNames));
    for i = 1:length(missing_pos_labs)
        filt.pos.(missing_pos_labs{i})(:) = NaN;
    end
    for i = 1:length(missing_vel_labs)
        filt.vel.(missing_vel_labs{i})(:) = NaN;
    end
    %Set the "ALL" values
    filt.pos(:,~isnan(filt.pos{'ALL',:})) = repmat(filt.pos(end,~isnan(filt.pos{'ALL',:})),11,1);
    filt.vel(:,~isnan(filt.vel{'ALL',:})) = repmat(filt.vel(end,~isnan(filt.vel{'ALL',:})),11,1);  
    %Check the parameters for valid values
    filt.pos = MakeCycAvg__filterParameterCheck(filt.pos);
    filt.vel = MakeCycAvg__filterParameterCheck(filt.vel);
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
                temp = interp1(te(~isnan(rel_trace)),rel_trace(~isnan(rel_trace)),te);
                temp(isnan(temp)) = 0;    
                for f = 1:length(pos_labs)
                    temp = filterTrace(pos_labs{f},temp,filt.pos.(pos_labs{f})(i),te,ts);
                end
                Data_pos_filt.(var_n) = temp;
                if contains(Data.info.name,'Impulse')
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
        if isfield(Data_vel,var_n)&&any(~isnan(Data_vel.(var_n)))
            if contains(Data.info.name,'Impulse')
                rel_t = t_cond;                
                Data_vel.([var_n,'_cond']) = Data_vel.(var_n)(reshape(keep_inds,[],1)); 
                rel_trace = Data_vel.([var_n,'_cond']);
            else
                rel_t = ts;
                rel_trace = Data_vel.(var_n);
            end
            temp = interp1(rel_t(~isnan(rel_trace)),rel_trace(~isnan(rel_trace)),rel_t);
            temp(isnan(temp)) = 0;
            for f = 1:length(vel_labs)
                temp = filterTrace(vel_labs{f},temp,filt.vel.(vel_labs{f})(i),rel_t,rel_t);
            end
            Data_vel_filt.(var_n) = temp;
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
    %% Make the CycAvg Struct
    %Data Traces
    CycAvg.t = Data_cyc.t;
    if all(size(Data_cyc.stim)>1) %multiple head traces
        CycAvg.stim_cyc = Data_cyc.stim(:,keep_tr);
        CycAvg.stim_cycavg = mean(Data_cyc.stim(:,keep_tr),2,'omitnan')';
        CycAvg.stim_cycstd = std(Data_cyc.stim(:,keep_tr),0,2,'omitnan')';
        CycAvg.stim = CycAvg.stim_cyc;
    else
        CycAvg.stim = Data_cyc.stim';
    end   
    traces = filt.vel.Properties.RowNames(1:end-1);
    for i = 1:length(traces)
        trac = lower(traces{i}(1:2));
        var_n = [traces{i}(1),'E_Vel_',traces{i}(2:end)];
        if isfield(Data_cyc,var_n)
            if all(size(Data_cyc.(var_n)(:,keep_tr))>1)
                CycAvg.([trac,'_cycavg']) = mean(Data_cyc.(var_n)(:,keep_tr),2,'omitnan')';
                CycAvg.([trac,'_cycstd']) = std(Data_cyc.(var_n)(:,keep_tr),0,2,'omitnan')';
                CycAvg.([trac,'_cyc']) = Data_cyc.(var_n)(:,keep_tr)';
            else
                CycAvg.([trac,'_cyc']) = Data_cyc.(var_n)';
            end
            if contains(fname,{'Activation','Step'})
                CycAvg.([trac,'_cyc_prefilt']) = Data_vel.(var_n); 
            elseif contains(fname,'Impulse')&&isfield(Data_vel,var_n)
                CycAvg.([trac,'_cyc_prefilt']) = Data_vel.(var_n)(Data_cyc.keep_inds(:,keep_tr));      
                CycAvg.([trac,'_cyc_QP']) = Data_cyc.([var_n,'_smooth'])(:,keep_tr); 
                long_t = repmat(Data_cyc.t,1,sum(keep_tr));
                CycAvg.([trac,'_saccade_time']) = long_t(logical(reshape(Data_cyc.([var_n,'_saccade'])(:,keep_tr),[],1))); 
            end
        end
    end         
    %File Information and intermediate steps
    %Other relevant items
    CycAvg.Fs = Data.Fs; 
    CycAvg.info = Data.info;
    CycAvg.name = ['CycAvg_',fname];
    CycAvg.filt = filt;
    CycAvg.cyclist = find(keep_tr);
    CycAvg.keep_tr = keep_tr;
    CycAvg.detec_tr = Data.detec_tr;
    CycAvg.t_interp = filt.t_interp;
    %Steps in Data Analsis
    CycAvg.Data = Data; %Direct output of segment file
    CycAvg.Data_rawpos = Data_pos;
    CycAvg.Data_filtpos = Data_pos_filt;
    CycAvg.Data_rawvel = Data_vel;
    CycAvg.Data_filtvel = Data_vel_filt;
    CycAvg.Data_allcyc = Data_cyc;
else %Just update cycle average values from new CycAvg.keep_tr value
    Data_cyc = CycAvg.Data_allcyc;
    traces = CycAvg.filt.vel.Properties.RowNames(1:end-1);
    fname = CycAvg.info.name;
    Data_vel = CycAvg.Data_rawvel;  
    keep_tr = CycAvg.keep_tr;
    CycAvg.filt.keep_tr = keep_tr;
    filt = CycAvg.filt;
    if all(size(Data_cyc.stim)>1) %multiple head traces
        CycAvg.stim_cyc = Data_cyc.stim(:,keep_tr)';
        CycAvg.stim_cycavg = mean(Data_cyc.stim(:,keep_tr),2,'omitnan')';
        CycAvg.stim_cycstd = std(Data_cyc.stim(:,keep_tr),0,2,'omitnan')';
        CycAvg.stim = CycAvg.stim_cyc;
    else
        CycAvg.stim = Data_cyc.stim';
    end       
    for i = 1:length(traces)
        trac = lower(traces{i}(1:2));
        var_n = [traces{i}(1),'E_Vel_',traces{i}(2:end)];
        if isfield(Data_cyc,var_n)
            if all(size(Data_cyc.(var_n)(:,keep_tr))>1)
                CycAvg.([trac,'_cycavg']) = mean(Data_cyc.(var_n)(:,keep_tr),2,'omitnan')';
                CycAvg.([trac,'_cycstd']) = std(Data_cyc.(var_n)(:,keep_tr),0,2,'omitnan')';
                CycAvg.([trac,'_cyc']) = Data_cyc.(var_n)(:,keep_tr)';
            else
                CycAvg.([trac,'_cyc']) = Data_cyc.(var_n)';
            end
            if contains(fname,{'Activation','Step'})
                CycAvg.([trac,'_cyc_prefilt']) = Data_vel.(var_n); 
            elseif contains(fname,'Impulse')&&isfield(Data_vel,var_n)
                CycAvg.([trac,'_cyc_prefilt']) = Data_vel.(var_n)(Data_cyc.keep_inds(:,keep_tr));      
                CycAvg.([trac,'_cyc_QP']) = Data_cyc.([var_n,'_smooth'])(:,keep_tr); 
                long_t = repmat(Data_cyc.t,1,sum(keep_tr));
                CycAvg.([trac,'_saccade_time']) = long_t(logical(reshape(Data_cyc.([var_n,'_saccade'])(:,keep_tr),[],1))); 
            end
        end
    end         
    CycAvg.cyclist = find(keep_tr);
end
end