%% Filter position and velocity
function [CycAvg,filt,plot_info] = MakeCycAvg__filterTraces(Data,filt,plot_info,CycAvg)
if nargin < 4
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
    pos_labs = filt.pos.Properties.VariableNames;
    vel_labs = filt.vel.Properties.VariableNames;
    %Set the "ALL" values
    filt.pos(plot_info.traces_pos,~isnan(filt.pos{'ALL',:})) = repmat(filt.pos(end,~isnan(filt.pos{'ALL',:})),length(plot_info.traces_pos),1); 
    filt.vel(plot_info.traces_vel,~isnan(filt.vel{'ALL',:})) = repmat(filt.vel(end,~isnan(filt.vel{'ALL',:})),length(plot_info.traces_vel),1);  
    %Clear the 'ALL' value on the filter
    filt.pos{end,:} = NaN; 
    filt.vel{end,:} = NaN; 
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
                if sum(~isnan(Data.(var_n)))>2               
                    temp = interp1(te(~isnan(rel_trace)),rel_trace(~isnan(rel_trace)),te);
                else
                    temp = rel_trace;
                end
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
        Data_pos_filt.info = Data.info;
        [Data_vel,Data_pos_filt] = angpos2angvel(Data_pos_filt); %Only uses X, Y, Z data right now
        Data_vel.t = ts;
        Data_vel.stim = stim;
    end
    if ~isempty(Data_pos_filt)
        eye_max = NaN(length(plot_info.traces_pos),2);
        for t = 1:length(plot_info.traces_pos)
            tr = strrep(strrep([plot_info.traces_pos{t}(1),'E_Position_',plot_info.traces_pos{t}(2)],'_L','_LARP'),'_R','_RALP');
            if isfield(Data_pos_filt,tr)
                eye_max(t,:) = [min(Data_pos_filt.(tr)),max(Data_pos_filt.(tr))];
            end
        end
        eye_max(eye_max<0) = 10*floor(eye_max(eye_max<0)/10);
        eye_max(eye_max>=0) = 10*ceil(eye_max(eye_max>=0)/10);
        plot_info.YLim.Pos = [max([min(eye_max(:,1)),-40]),min([max(eye_max(:,2)),40])]; %Keep bounded by ±40
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
    Data_cyc.keep_inds = keep_inds;
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
    %Set Velocity YLim
    eye_max = NaN(length(plot_info.traces_vel),2);
    for t = 1:length(plot_info.traces_vel)
        tr = strrep(strrep([plot_info.traces_vel{t}(1),'E_Vel_',plot_info.traces_vel{t}(2)],'_L','_LARP'),'_R','_RALP');
        if isfield(Data_cyc,tr)
            eye_max(t,:) = [min(1.2*median(Data_cyc.(tr),2)),max(1.2*median(Data_cyc.(tr),2))];
        end
    end
    if contains(Data.info.name,'Impulse')
        eye_max(t+1,:) = [min(Data_vel_filt.stim),max(Data_vel_filt.stim)];
    end
    eye_max(eye_max<0) = 25*floor(eye_max(eye_max<0)/25);
    eye_max(eye_max>=0) = 25*ceil(eye_max(eye_max>=0)/25);
    plot_info.YLim.Vel = [min(eye_max(:,1)),max(eye_max(:,2))];
    if plot_info.YLim.Vel(1)==plot_info.YLim.Vel(2)
        plot_info.YLim.Vel(2) = plot_info.YLim.Vel(1)+1;
    end
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