% This function can be used to try using different parameters to go from 
% segement to CycAvg 
function CycAvg = plotSeg2CycAvg(CycAvg,Data,params)
    load('VNELcolors.mat','colors')
    % Fill colors for cycle selection
    %Color of cycle on long graph
    colors.cyc_keep = [0.75 1 0.75]; 
    colors.cyc_rm = [1 1 1];
    if ~isempty(CycAvg)
        Data = CycAvg.raw_Data;
        if isempty(params)
            params.filt = CycAvg.filt;
            params.cyclist = CycAvg.cyclist;  
        end
    elseif isempty(Data)||isempty(params)
        return;
    end
    %% Make figure
    fig = figure;
    delete(findall(gcf,'type','annotation')) %in case there are leftover anotations
    fig.Units = 'inches';
    fig.Position = [0 0 11 10];
    %Title
    annotation('textbox',[0 .9 1 .1],'String',strrep(strrep(Data.info.dataType,'_',' '),'-',' '),'FontSize',14,...
        'HorizontalAlignment','center','EdgeColor','none');
    %% Extract raw position data
    info = Data.info;
    info.colors = colors;
    Fs = Data.Fs;
    if contains(info.dataType,'Activation')
        %preserve the time because this will be used to rejoin them
        te = Data.Time_Eye;
        ts = Data.Time_Stim; 
    else
        te = Data.Time_Eye - Data.Time_Eye(1);
        ts = Data.Time_Stim - Data.Time_Stim(1);
    end
    %Trigger multiplier
    if contains(info.dataType,'RotaryChair')
        if isfield(Data,'HeadVel_Z')
            stim = Data.HeadVel_Z;
        else
            stim = Data.HeadMPUVel_Z; 
        end   
    elseif contains(info.dataType,'aHIT')
        if contains(info.dataType,'LHRH')
            stim = Data.HeadVel_Z;
        elseif contains(info.dataType,'LARP')
            stim = (Data.HeadVel_X - Data.HeadVel_Y)/sqrt(2);
        elseif contains(info.dataType,'RALP')
            stim = (Data.HeadVel_X + Data.HeadVel_Y)/sqrt(2);
        else
            stim = Data.Trigger;
        end    
    else
        stim = Data.Trigger; 
    end       
    %Fix huge number of NaN values in torsion traces of NKI traces
    if contains(info.goggle_ver,'NKI')
        if sum(isnan(Data.LE_Position_X)) > 0.9*length(te) %less than 10% data integrity
            Data.LE_Position_X = zeros(length(te),1); %set to 0 so no torsion
        else
            Data.LE_Position_X = spline(te(~isnan(Data.LE_Position_X)),Data.LE_Position_X(~isnan(Data.LE_Position_X)),te);
        end    
        if sum(isnan(Data.RE_Position_X)) > 0.9*length(te)
            Data.RE_Position_X = zeros(length(te),1);
        else
            Data.RE_Position_X = spline(te(~isnan(Data.RE_Position_X)),Data.RE_Position_X(~isnan(Data.RE_Position_X)),te);
        end    
    end
    %% Cycle Align
    %Defines a variable type that determines what type of filtering and
    %analysis needs to be done. type = 1: cycle average, 2: full trace (vel
    %step/activation)
    [type,starts,ends,stims] = MakeCycAvg__alignCycles(info,Fs,ts,stim);
    %Remove any unnecessary trace the start and end
    keep_inds = zeros(ends(1)-starts(1)+1,length(starts));
    for i = 1:length(starts)
        keep_inds(:,i) = starts(i):ends(i);
    end
    keep_inds = keep_inds - starts(1)+1;
    if contains(info.dataType,'Activation')
        te = te(starts(1):ends(end));
        ts = ts(starts(1):ends(end));
        t_snip = reshape(ts(1:length(stims)),1,[]);
    else
        te = te(starts(1):ends(end)) - te(starts(1));
        ts = ts(starts(1):ends(end)) - ts(starts(1));
        t_snip = reshape(ts(1:length(stims))-ts(1),1,[]);
    end
    stim = stim(starts(1):ends(end));
    Data.LE_Position_X = Data.LE_Position_X(starts(1):ends(end));
    Data.LE_Position_Y = Data.LE_Position_Y(starts(1):ends(end));
    Data.LE_Position_Z = Data.LE_Position_Z(starts(1):ends(end));
    Data.RE_Position_X = Data.RE_Position_X(starts(1):ends(end));
    Data.RE_Position_Y = Data.RE_Position_Y(starts(1):ends(end));
    Data.RE_Position_Z = Data.RE_Position_Z(starts(1):ends(end));
    %% Set some Defaults
    keep_tr = ismember(1:length(starts),params.cyclist);
    filt_params_p = NaN(4*6,1);
    traces = {'l_x','r_x','l_y','r_y','l_z','r_z'};
    for i = 1:6
        if isfield(params.filt.pos.(traces{i}),'medianfilt')&&all(~isempty(params.filt.pos.(traces{i}).medianfilt))
            filt_params_p(i) = params.filt.pos.(traces{i}).medianfilt;
        end
        if isfield(params.filt.pos.(traces{i}),'spline')&&all(~isempty(params.filt.pos.(traces{i}).spline))
            filt_params_p(i+6) = params.filt.pos.(traces{i}).spline;
        end
        if isfield(params.filt.pos.(traces{i}),'sgolay')&&all(~isempty(params.filt.pos.(traces{i}).sgolay))
            filt_params_p(i+12) = params.filt.pos.(traces{i}).sgolay(1);
            filt_params_p(i+18) = params.filt.pos.(traces{i}).sgolay(2);
        end
    end
    filt_params_v = NaN(4,1);
    if isfield(params.filt.vel,'irlssmooth')&&all(~isempty(params.filt.vel.irlssmooth))
        filt_params_v(1) = params.filt.vel.irlssmooth;
    end
    if isfield(params.filt.vel,'spline')&&all(~isempty(params.filt.vel.spline))
        filt_params_v(2) = params.filt.vel.spline;
    end
    if isfield(params.filt.vel,'accelcutoff')&&all(~isempty(params.filt.vel.accelcutoff))
        filt_params_v(3) = params.filt.vel.accelcutoff;
    end
    if isfield(params.filt.vel,'medianfilt')&&all(~isempty(params.filt.vel.medianfilt))
        filt_params_v(4) = params.filt.vel.medianfilt;
    end
    %% Make plots      
    YLim_Pos = [-30 30];
    YLim_Vel = [-100 100];
    line_wid.norm = 0.5;
    line_wid.bold = 2;
    [filt,Data_calc,LE_V,RE_V,Data_cal,Data_In] = MakeCycAvg__filterTraces(type,filt_params_p,filt_params_v,te,ts,Data,keep_inds); 
    CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,[info.subject,'-',info.visit,'-',info.exp_date,'-',info.dataType,'.mat']);
    MakeCycAvg__plotFullCycAvg([],type,colors,line_wid,YLim_Pos,YLim_Vel,te,ts,t_snip,stim,stims,Data,Data_In,Data_cal,Data_calc,LE_V,RE_V,CycAvg,keep_inds,keep_tr);   
end