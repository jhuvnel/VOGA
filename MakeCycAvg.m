%% MakeCycAvg.m
% This function was created to automate the data analysis pipeline for the
% LDVOG data sets as much as possible.
% Cycle average before filtering

function [CycAvg,analyzed] = MakeCycAvg(Data,Cyc_Path,from_file)
clc; 
% Plot defaults
all_traces = {'LX','RX','LY','RY','LZ','RZ','LLARP','RLARP','LRALP','RRALP'};
traces_pos1 = all_traces(1:6); %XYZ pos
% Set colors
load('VNELcolors.mat','colors')
% Fill colors for cycle selection
colors.cyc_keep = [0.85 0.85 0.85];
colors.cyc_rm = [1 1 1];
% Set Experimentor/version
if ~any(contains(extractfield(dir(userpath),'name'),'VOGA_VerInfo.txt'))
    VOGA__setVersion;
end
data = readtable('VOGA_VerInfo.txt','ReadVariableNames',false);
version = data{1,2}{:};
Experimenter = data{2,2}{:};
info = Data.info;
info.Analyzer = Experimenter;
info.ver = version;
info.colors = colors;
% Parse inputs and load defaults
try
    fname = Data.info.name;
catch
    fname = [Data.info.subject,'-',Data.info.visit,'-',Data.info.exp_date,'-',Data.info.goggle_ver,'-',Data.info.dataType,'.mat'];
end
%Load default filters for the goggle type
load('VOGA_LastUsedFilterParams.mat','filt_params')
filt1 = filt_params.filt1;
YLim = filt_params.YLim;       
% Initialize Figure
fig = figure(1);
clf; %in case there are leftover anotations
fig.Units = 'normalized';
fig.Position = [0 0 1 1];
fig.Units = 'inches';
screen_size = fig.Position;
fig.Position = screen_size - [0 0 4 0];
%Title
if contains(fname,'[')&&contains(fname,']')
    fig_title = strrep(strrep(strrep(fname,'_',' '),'-',' '),'.mat','');
    b1 = strfind(fig_title,'[');
    b2 = strfind(fig_title,']');
    fig_title(b1:b2) = strrep(fig_title(b1:b2),' ','-');
else
    fig_title = strrep(strrep(strrep(fname,'_',' '),'-',' '),'.mat','');
end
annotation('textbox',[0 .9 1 .1],'String',fig_title,'FontSize',14,...
    'HorizontalAlignment','center','EdgeColor','none');
% Extract raw position data
info.TriggerShift2 = 0; %Shifting done manually in this file
Fs = Data.Fs;
if contains(info.goggle_ver,'GNO') %No raw position, just velocity
    te = Data.Time_Eye - Data.Time_Eye(1);
    ts = Data.Time_Stim - Data.Time_Stim(1);
    if contains(info.dataType,{'LH','RH'})
        stim1 = Data.HeadVel_Z;
    elseif contains(info.dataType,{'LA','RP'})
        stim1 = -Data.HeadVel_L;
    elseif contains(info.dataType,{'RA','LP'})
        stim1 = Data.HeadVel_R;
    else
        stim1 = Data.HeadVel_Z;
    end 
elseif contains(info.goggle_ver,'ESC') %No raw position, just velocity
    te = Data.Time_Eye - Data.Time_Eye(1);
    ts = Data.Time_Stim - Data.Time_Stim(1);
    if contains(info.dataType,{'LH','RH'})
        stim1 = Data.HeadVel_Z;
    elseif contains(info.dataType,{'LA','RP'})
        stim1 = Data.HeadVel_L;
    elseif contains(info.dataType,{'RA','LP'})
        stim1 = Data.HeadVel_R;
    else
        stim1 = Data.HeadVel_Z;
    end    
else
    if contains(info.dataType,'Activation')
        %preserve the time because this will be used to rejoin them
        te = Data.Time_Eye;
        ts = Data.Time_Stim; 
    else
        te = Data.Time_Eye - Data.Time_Eye(1);
        ts = Data.Time_Stim - Data.Time_Stim(1);
    end
    %Assign Trigger
    if contains(info.goggle_ver,'Moogles') %MOOG room coil system
        stim1 = Data.Trigger;
    else %NKI or LDVOG Trigger = internal gyro or comes from the PCU for eeVOR
        if contains(info.dataType,'RotaryChair')
            if isfield(Data,'HeadVel_Z')
                stim1 = Data.HeadVel_Z;
            else
                stim1 = Data.HeadMPUVel_Z; 
            end   
        elseif contains(info.dataType,'aHIT')
            if contains(info.dataType,'LHRH')
                stim1 = Data.HeadVel_Z;
            elseif contains(info.dataType,'LARP')
                stim1 = (Data.HeadVel_X - Data.HeadVel_Y)/sqrt(2);
            elseif contains(info.dataType,'RALP')
                stim1 = (Data.HeadVel_X + Data.HeadVel_Y)/sqrt(2);
            else
                stim1 = Data.Trigger;
            end    
        else
            stim1 = Data.Trigger; 
        end       
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
end
% Assign Type and set detected traces
if contains(info.goggle_ver,{'GNO','ESC'}) %No raw pos traces 
    type = 3;
    if isfield(Data,'DetectedTraces_HeadVel')
        detec_head = Data.DetectedTraces_HeadVel;
    else
        detec_head = [];
    end
elseif contains(info.dataType,{'Activation','Step'}) %No cycle averaging
    type = 2;
    detec_head = [];
else
    type = 1;
    detec_head = [];
end
% Set some defaults
% Cycle Align
[~,t_snip] = MakeCycAvg__alignCycles(info,Fs,ts,stim1,[]);
if type == 1
    filt1.vel.irlssmooth(end) = round(length(t_snip)*0.16); %heuristic
end
line_wid.norm = 0.5;
line_wid.bold = 2;
if contains(info.goggle_ver,'GNO')&&contains(info.dataType,{'LH','RH'})
    traces_vel1 = {'RZ'};
elseif contains(info.goggle_ver,'GNO')&&contains(info.dataType,{'LA','RP'})
    traces_vel1 = {'RLARP'};
elseif contains(info.goggle_ver,'GNO')&&contains(info.dataType,{'RA','LP'})
    traces_vel1 = {'RRALP'};   
elseif contains(info.goggle_ver,'ESC')&&contains(info.dataType,{'LH','RH'})
    traces_vel1 = {'LZ'};
elseif contains(info.goggle_ver,'ESC')&&contains(info.dataType,{'LA','RP'})
    traces_vel1 = {'LLARP'};
elseif contains(info.goggle_ver,'ESC')&&contains(info.dataType,{'RA','LP'})
    traces_vel1 = {'LRALP'};     
elseif contains(info.dataType,{'X','Y'})
    traces_vel1 = all_traces(1:6); %LRZ vel 
else
    traces_vel1 = all_traces(5:10); %LRZ vel 
end
%% No options, just reanalyze file with existing params
if from_file 
    traces_pos = traces_pos1;
    traces_vel = traces_vel1; 
    if isempty(YLim)
        if type == 3 %Impulses
            YLim.Pos = [-30 30];
            YLim.Vel = [-250 250];
        else 
            YLim.Pos = [-30 30];
            YLim.Vel = [-50 50];
        end
    end
    info.TriggerShift2 = 0;
    [stim,t_snip,stims,keep_inds,detec_tr] = MakeCycAvg__alignCycles(info,Fs,ts,stim1,detec_head);
    keep_tr = true(1,size(keep_inds,2));
    file = [Cyc_Path,filesep,'CycAvg_',fname];
    a = load(file);
    CycAvg2 = a.CycAvg;
    if length(CycAvg2.keep_tr) == length(keep_tr) %Only if they are the same size
        filt = CycAvg2.filt;
        keep_tr = CycAvg2.keep_tr;
        if isfield(CycAvg2,'t_interp')
            t_interp = CycAvg2.t_interp;
        else
            t_interp = [];
        end
        [filt,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc] = MakeCycAvg__filterTraces(filt,keep_inds,te,ts,t_snip,stim,stims,Data,t_interp);
        CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);
        ha = MakeCycAvg__plotFullCycAvg([],type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);
    else
        disp('Not a compatible CycAvg file.')
    end 
    % Create to Save
    CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);
    CycAvg = ParameterizeCycAvg(CycAvg);
    filt_params.filt1 = filt;
    filt_params.YLim = YLim;
    VOGA__saveLastUsedParams(filt_params)
    analyzed = 1;
    return;
end
%% Once analyzeable, here is the while loop they stay in until saving or exiting
%You can change the order/existance of these options without ruining
%anything because the comparrisons are all string based
opts = {'Set Y-axis Lim','Choose Coordinates','Filter Position',...
    'Filter Velocity','Manual QPR','Select Cycles','Shift Trigger','Load from File',...
    'Not Analyzeable','Start Over','Save'};
ind = find(contains(opts,'Start Over')); %Run the start procedure first
while ~strcmp(opts{ind},'Save') %Run until it's ready to save or just hopeless
    if strcmp(opts{ind},'Start Over')
        t_interp = [];
        traces_pos = traces_pos1;
        traces_vel = traces_vel1; 
        if isempty(YLim)
            if type == 3 %Impulses
                YLim.Pos = [-30 30];
                YLim.Vel = [-250 250];
            else 
                YLim.Pos = [-30 30];
                YLim.Vel = [-50 50];
            end
        end
        info.TriggerShift2 = 0;
        [stim,t_snip,stims,keep_inds,detec_tr] = MakeCycAvg__alignCycles(info,Fs,ts,stim1,detec_head);
        if contains(info.dataType,'Impulse')&&~isempty(detec_tr) %remove erroneous head traces using auto-detected traces as a template
            head_templ = mean(stims(:,detec_tr),2);
            head_templ_sd = std(stims(:,detec_tr),[],2);            
            out_of_bounds = stims > head_templ+3*head_templ_sd | stims < head_templ-3*head_templ_sd;
            keep_tr = ~any(out_of_bounds(t_snip<0.23&t_snip>0.1,:))&any(stims(t_snip<0.2,:)<10)&any(stims(t_snip<0.4&t_snip>0.1,:)<0);
        else
            keep_tr = true(1,size(keep_inds,2));
        end        
        [filt,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc] = MakeCycAvg__filterTraces(filt1,keep_inds,te,ts,t_snip,stim,stims,Data,t_interp);
        CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);
        ha = MakeCycAvg__plotFullCycAvg([],type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);
    elseif strcmp(opts{ind},'Shift Trigger')
        new_TrigShift = cellfun(@str2double,inputdlgcol('Trigger Shift (samples): ','Shift',[1 15],{num2str(info.TriggerShift2)},'on',1,[screen_size(3)-4 screen_size(4)-1.25 1.75 1.25]));
        if ~isempty(new_TrigShift)
            info.TriggerShift2 = round(new_TrigShift);
            [stim,t_snip,stims,keep_inds,detec_tr] = MakeCycAvg__alignCycles(info,Fs,ts,stim1,detec_head);
            if size(keep_inds,2) > length(keep_tr)
                old_keep_tr = keep_tr;
                keep_tr = [old_keep_tr;true(1,size(keep_inds,2)-length(keep_tr))];
            elseif size(keep_inds,2) < length(keep_tr)
                keep_tr = keep_tr(1:size(keep_inds,2));
            end
            [filt,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc] = MakeCycAvg__filterTraces(filt,keep_inds,te,ts,t_snip,stim,stims,Data,t_interp);
            CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);            
            ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);
        end
    elseif strcmp(opts{ind},'Filter Position')
        %Get new parameter values
        traces = strrep(strrep(filt.pos.Properties.RowNames,'ALP',''),'ARP','');  
        filters = filt.pos.Properties.VariableNames;
        prompt = repmat(traces,1,length(filters));
        prompt(1,:) = strcat(filters,[newline,traces{1}]);
        prompt = reshape(prompt,[],1)';        
        dlgtitle = 'Filter position';
        definput = strrep(cellfun(@(x) num2str(x,10),table2cell(filt.pos),'UniformOutput',false),'NaN','');
        temp_filt_params_p = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 10],definput,'on',length(prompt)/length(traces),[screen_size(3)-3.5 screen_size(4)-7 3.5 7]));
        if ~isempty(temp_filt_params_p) %Didn't hit cancel
            filt.pos{:,:} = reshape(temp_filt_params_p,11,[]);
            [filt,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc] = MakeCycAvg__filterTraces(filt,keep_inds,te,ts,t_snip,stim,stims,Data,t_interp);
            CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);            
            ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);
        end
    elseif strcmp(opts{ind},'Filter Velocity')
        temp_filt_params_v = 'placeholder';
        while ~isempty(temp_filt_params_v) %Didn't hit cancel
        %Get new parameter values
            traces = strrep(strrep(filt.vel.Properties.RowNames,'ALP',''),'ARP','');  
            filters = filt.vel.Properties.VariableNames;
            prompt = repmat(traces,1,length(filters));
            prompt(1,:) = strcat(filters,[newline,traces{1}]);
            prompt = reshape(prompt,[],1)';          
            dlgtitle = 'Filter velocity';
            definput = strrep(cellfun(@(x) num2str(x,10),table2cell(filt.vel),'UniformOutput',false),'NaN','');
            temp_filt_params_v = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 10],definput,'on',length(prompt)/length(traces),[screen_size(3)-5 screen_size(4)-7 5 7]));        
            if ~isempty(temp_filt_params_v)
                filt.vel{:,:} = reshape(temp_filt_params_v,11,[]);
                [filt,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc] = MakeCycAvg__filterTraces(filt,keep_inds,te,ts,t_snip,stim,stims,Data,t_interp);
                CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);                
                ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);        
            end
        end
    elseif strcmp(opts{ind},'Manual QPR') 
        clc;
        disp('Select the start and end time point to linearly interpolate over on the cycle average graph.')
        good_rng = 'Redo';
        if type == 1
            cyc_ax = ha(4);
        elseif type == 3
            cyc_ax = ha(3);
        end
        while strcmp(good_rng,'Redo')       
            [t_spline,~] = ginput(2);
            axes(cyc_ax)
            hold on
            h1 = xline(t_spline(1),'LineWidth',10);
            h2 = xline(t_spline(2),'LineWidth',10);
            good_rng = questdlg('Keep or redo the range?','','Keep','Redo','Redo');  
            delete(h1)
            delete(h2)
            hold off
        end
        if strcmp(good_rng,'Keep')
            t_interp = find(t_snip<t_spline(1),1,'last'):find(t_snip>t_spline(2),1,'first');
            [filt,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc] = MakeCycAvg__filterTraces(filt,keep_inds,te,ts,t_snip,stim,stims,Data,t_interp);
            CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);
            ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);        
        end
    elseif strcmp(opts{ind},'Select Cycles')
        [keep_tr,ha,tf] = MakeCycAvg__selectCycles(ha,type,keep_tr,Data_cyc,screen_size,traces_vel);   
        while tf
            CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);
            ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);
            [keep_tr,ha,tf] = MakeCycAvg__selectCycles(ha,type,keep_tr,Data_cyc,screen_size,traces_vel);   
        end
    elseif strcmp(opts{ind},'Load from File')
        auto = questdlg('Automatically detect file?','','Yes','No','No');
        if ~isempty(auto)
            cyc_files = extractfield(dir([Cyc_Path,filesep,'*.mat']),'name');
            if strcmp(auto,'Yes')&&any(ismember(cyc_files,['CycAvg_',fname]))
                file = [Cyc_Path,filesep,'CycAvg_',fname];
                tf = 1;
            else
                %See what files exist in the Cyc_Path
                cyc_files = extractfield(dir([Cyc_Path,filesep,'*.mat']),'name');
                [indx,tf] = nmlistdlg('PromptString','Select an analyzed file to us:',...
                               'SelectionMode','single',...
                               'ListSize',[500 600],...
                               'ListString',cyc_files); 
                file = [Cyc_Path,filesep,cyc_files{indx}];
            end
            if tf
                a = load(file);
                CycAvg2 = a.CycAvg;
                if length(CycAvg2.keep_tr) == length(keep_tr) %Only if they are the same size
                    filt = CycAvg2.filt;
                    keep_tr = CycAvg2.keep_tr;
                    if isfield(CycAvg2,'t_interp')
                        t_interp = CycAvg2.t_interp;
                    else
                        t_interp = [];
                    end
                    [filt,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc] = MakeCycAvg__filterTraces(filt,keep_inds,te,ts,t_snip,stim,stims,Data,t_interp);
                    CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);
                    ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);
                else
                    disp('Not a compatible CycAvg file.')
                end             
            end
        end
    elseif strcmp(opts{ind},'Set Y-axis Lim')
        %Get new parameter values
        prompt = {['Set Y-axis limits',newline,newline,'Position:',newline,newline,'Lower Limit:'],...
            'Upper Limit:',['Velocity:',newline,newline,'Lower Limit:'],'Upper Limit:'};
        dlgtitle = 'Y-axis Limits';
        definput = cellfun(@(x) num2str(x,10),num2cell([YLim.Pos,YLim.Vel]),'UniformOutput',false);
        out_nums = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 18],definput,'on',2,[screen_size(3)-4 screen_size(4)-2.25 3 2.25]));
        if ~isempty(out_nums)
            %Check to make sure they aren't reversed
            YLim.Pos = sort([out_nums(1),out_nums(2)]);
            YLim.Vel = sort([out_nums(3),out_nums(4)]);
        end
        ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);
    elseif strcmp(opts{ind},'Choose Coordinates') 
        if type ~= 3
            [ind3,tf] = nmlistdlg('PromptString','Select position traces:',...
                               'SelectionMode','multiple',...
                               'InitialValue',find(ismember(all_traces,traces_pos)),...
                               'ListSize',[100 150],...
                               'ListString',all_traces,...
                               'Position',[screen_size(3)-4,screen_size(4)-4,2,4]);           
            if tf
                traces_pos = all_traces(ind3);
            end
        end
        [ind2,tf] = nmlistdlg('PromptString','Select velocity traces:',...
                           'SelectionMode','multiple',...
                           'InitialValue',find(ismember(all_traces,traces_vel)),...
                           'ListSize',[100 150],...
                           'ListString',all_traces,...
                           'Position',[screen_size(3)-4,screen_size(4)-4,2,4]);           
        if tf
            traces_vel = all_traces(ind2);
        end 
        ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);
    elseif strcmp(opts{ind},'Not Analyzeable')
        CycAvg = Data;
        analyzed = 0;
        return;
    else
        disp('This case not yet been defined.')
    end
    [ind,tf2] = nmlistdlg('PromptString','Select an action:',...
                       'SelectionMode','single',...
                       'ListSize',[150 150],...
                       'ListString',opts,...
                       'Position',[screen_size(3)-4,screen_size(4)-3.75,3,3.75]);  
    if tf2 == 0 %Treat this like an exit
        CycAvg = [];
        analyzed = 0;
        return;
    end   
end
%% Create to Save
CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);
CycAvg = ParameterizeCycAvg(CycAvg);
filt_params.filt1 = filt;
filt_params.YLim = YLim;
VOGA__saveLastUsedParams(filt_params)
analyzed = 1;
end