%% Make Cycle Average Struct
function CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc)
    %Data Traces
    CycAvg.t = Data_cyc.t;
    if all(size(Data_cyc.stim)>1) %multiple head traces
        CycAvg.stim_cyc = Data_cyc.stim(:,keep_tr);
        CycAvg.stim_cycavg = mean(Data_cyc.stim(:,keep_tr),2)';
        CycAvg.stim_cycstd = std(Data_cyc.stim(:,keep_tr),0,2)';
        CycAvg.stim = CycAvg.stim_cyc;
    else
        CycAvg.stim = Data_cyc.stim';
    end   
    traces = filt.vel.Properties.RowNames(1:end-1);
    for i = 1:length(traces)
        trac = lower(traces{i}(1:2));
        var_n = [traces{i}(1),'E_Vel_',traces{i}(2:end)];
        if isfield(Data_cyc,var_n)
            CycAvg.([trac,'_cycavg']) = mean(Data_cyc.(var_n)(:,keep_tr),2)';
            CycAvg.([trac,'_cycstd']) = std(Data_cyc.(var_n)(:,keep_tr),0,2)';
            CycAvg.([trac,'_cyc']) = Data_cyc.(var_n)(:,keep_tr)';
            if contains(Data.info.goggle_ver,'GNO')&&isfield(Data_cyc,[var_n,'_QPR'])
                CycAvg.([trac,'_cyc_fit']) = Data_cyc.([var_n,'_QPR'])(:,keep_tr)';
                CycAvg.([trac,'_cycavg_fit']) = mean(Data_cyc.([var_n,'_QPR'])(:,keep_tr),2)';
                CycAvg.([trac,'_cycstd_fit']) = std(Data_cyc.([var_n,'_QPR'])(:,keep_tr),0,2)';                
            end
        end
    end         
    %File Information and intermediate steps
    %Other relevant items
    CycAvg.old_Fs = Data.Fs;
    CycAvg.Fs = Fs; 
    CycAvg.info = info;
    CycAvg.name = ['CycAvg_',fname];
    CycAvg.filt = filt;
    CycAvg.cyclist = find(keep_tr);
    CycAvg.keep_tr = keep_tr;
    CycAvg.detec_tr = detec_tr;
    CycAvg.t_interp = t_interp;
    %Steps in Data Analsis
    CycAvg.Data = Data; %Direct output of segment file
    CycAvg.Data_rawpos = Data_pos;
    CycAvg.Data_filtpos = Data_pos_filt;
    CycAvg.Data_rawvel = Data_vel;
    CycAvg.Data_filtvel = Data_vel_filt;
    CycAvg.Data_allcyc = Data_cyc;    
end