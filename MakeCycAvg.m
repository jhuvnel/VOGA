%% MakeCycAvg.m
% This function was created to automate the data analysis pipeline for the
% LDVOG data sets as much as possible.
% Cycle average before filtering

function done = MakeCycAvg(Path,code_Path,exp_types)
clc;       
% Set colors
load('VNELcolors.mat','colors')
% Fill colors for cycle selection
colors.cyc_keep = [0.85 0.85 0.85];
colors.cyc_bold_k = [0 0.75 0];
colors.cyc_rm = [1 1 1];
colors.cyc_bold_r = [0.85 0.85 0.85];
% Set paths
Seg_Path = [Path,filesep,'Segmented Files'];
Cyc_Path = [Path,filesep,'Cycle Averages'];
% Set Experimentor/version
if ~any(contains(extractfield(dir(code_Path),'name'),'VerInfo.txt'))
    writeInfoFile(code_Path);
end
data = readtable([code_Path,filesep,'VerInfo.txt'],'ReadVariableNames',false);
version = data{1,2}{:};
Experimenter = data{2,2}{:};
%% Load in data
progress_tab = assessProgress(Path);
progress_i = [find(~progress_tab{:,2}&~progress_tab{:,3});find(progress_tab{:,2}|progress_tab{:,3})]; %put unanalyzed files at the top
if ~isempty(exp_types)
    rel_file = false(length(progress_tab{:,1}),length(exp_types));
    for i = 1:length(exp_types)
        components = split(exp_types{i},'-');
        cont = true(length(progress_tab{:,1}),1);
        for j = 1:length(components)
            cont = cont&contains(progress_tab{:,1},components{j});
        end
        rel_file(:,i) = cont;
    end
    progress_i(~any(rel_file,2)) = [];
end
%Change color of option depending on whether it's been analyzed
font_col = repmat({'black'},length(progress_i),1);
font_col(progress_tab{progress_i,2}) = {'green'};
font_col(progress_tab{progress_i,3}) = {'red'};
list = strcat('<HTML><FONT color="',font_col,'">',table2cell(progress_tab(progress_i,1)),'</FONT></HTML>');
[indx,tf] = nmlistdlg('PromptString','Unattempted = Black, Analyzed = Green, Not Analyzeable = Red. Select a file to analyze:',...
                       'SelectionMode','single',...
                       'ListSize',[500 600],...
                       'ListString',list);  
%If the user selects cancel
if ~tf
    disp('Operation Ended.')
    done = true;
    return; 
else
    done = false;
end
In_FileName = progress_tab{progress_i(indx),1}{:};
load([Seg_Path,filesep,In_FileName],'Data');
%% Initialize Figure
fig = figure(1);
delete(findall(gcf,'type','annotation')) %in case there are leftover anotations
fig.Units = 'inches';
fig.Position = [0 0 11 10];
%Title
if contains(In_FileName,'[')&&contains(In_FileName,']')
    fig_title = strrep(strrep(strrep(In_FileName,'_',' '),'-',' '),'.mat','');
    b1 = strfind(fig_title,'[');
    b2 = strfind(fig_title,']');
    fig_title(b1:b2) = strrep(fig_title(b1:b2),' ','-');
else
    fig_title = strrep(strrep(strrep(In_FileName,'_',' '),'-',' '),'.mat','');
end
annotation('textbox',[0 .9 1 .1],'String',fig_title,'FontSize',14,...
    'HorizontalAlignment','center','EdgeColor','none');
%% Extract raw position data
info = Data.info;
info.Analyzer = Experimenter;
info.ver = version;
info.colors = colors;
info.TriggerShift2 = 0; %Shifting done manually in this file
Fs = Data.Fs;
if contains(info.goggle_ver,'GNO') %No raw position, just velocity
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
% Assign Type
if contains(info.goggle_ver,'GNO') %No raw pos traces 
    type = 3;
elseif contains(info.dataType,{'Activation','Step'}) %No cycle averaging
    type = 2;
else
    type = 1;
end
%% Set some Defaults
% Cycle Align
[~,t_snip] = MakeCycAvg__alignCycles(info,Fs,ts,stim1);
%Filter Traces with intial guesses
if type == 1
    if contains(info.goggle_ver,'NKI')
        pos_med = [21;21;3;3;3;3];
        pos_spline = [1;1;0.999995;0.999995;0.9999995;0.9999995];
        pos_sgolay = [2*ones(6,1);5*ones(6,1)];
        if contains(info.dataType,'Sine')
            vel_smooth = round(t_snip(end)*40);
        else
            vel_smooth = 30;
        end
        vel_spline = NaN;
        vel_acc = NaN; 
        vel_med = NaN;
    else %LDVOG, Moogles
        pos_med = [11;11;3;3;3;3];
        pos_spline = [0.99995;0.99995;0.999995;0.999995;0.9999995;0.9999995];
        pos_sgolay = [2*ones(6,1);5*ones(6,1)];
        if contains(info.dataType,'Sine')
            vel_smooth = round(t_snip(end)*40);
        else
            vel_smooth = 30;
        end
        vel_spline = NaN;
        vel_acc = NaN; 
        vel_med = NaN;
    end
elseif type == 2 %velstep and activation
    pos_med = [11;11;5;5;3;3];
    pos_spline = NaN(6,1);
    pos_sgolay = [2*ones(6,1);5*ones(6,1)];
    vel_smooth = NaN;
    vel_spline = NaN;
    vel_acc = 10;
    vel_med = Fs+1-mod(Fs,2); 
elseif type == 3 
    %Coded for GNO right now
    pos_med = NaN(6,1);
    pos_spline = NaN(6,1);
    pos_sgolay = NaN(12,1);
    vel_smooth = 2;
    vel_spline = NaN;
    vel_acc = NaN;
    vel_med = NaN;
end
line_wid.norm = 0.5;
line_wid.bold = 2;
%% Once analyzeable, here is the while loop they stay in until saving or exiting
%You can change the order/existance of these options without ruining
%anything because the comparrisons are all string based
opts = {'Set Y-axis Lim','Filter Position','Filter Velocity','Select Cycles','Shift Trigger','Not Analyzeable','Start Over','Save'};
ind = 7; %Run the start procedure first
while ~strcmp(opts{ind},'Save') %Run until it's ready to save or just hopeless
    if strcmp(opts{ind},'Start Over')
        filt_params_p = [pos_med;pos_spline;pos_sgolay];
        filt_params_v = [vel_smooth;vel_spline;vel_acc;vel_med];        
        YLim_Pos = [-30 30];
        if type == 3
            YLim_Vel = [-50 250];
        else
            YLim_Vel = [-100 100];
        end
        info.TriggerShift2 = 0;
        [stim,t_snip,stims,keep_inds] = MakeCycAvg__alignCycles(info,Fs,ts,stim1);
        keep_tr = true(1,size(keep_inds,2)); 
        [filt,Data_calc,LE_V,RE_V,Data_cal,Data_In] = MakeCycAvg__filterTraces(type,filt_params_p,filt_params_v,te,ts,Data,keep_inds); 
        CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
        ha = MakeCycAvg__plotFullCycAvg([],type,colors,line_wid,YLim_Pos,YLim_Vel,te,ts,t_snip,stim,stims,Data,Data_In,Data_cal,Data_calc,LE_V,RE_V,CycAvg,keep_inds,keep_tr);
    elseif strcmp(opts{ind},'Not Analyzeable')
        save([Cyc_Path,filesep,'NotAnalyzeable_',In_FileName],'Data')
        savefig([Cyc_Path,filesep,'NotAnalyzeable_',In_FileName(1:end-4),'.fig'])
        close;
        %If the file was previously designated as analyzeable, remove that file
        d1 = dir([Cyc_Path,filesep,'*.mat']);
        if ~isempty(d1)
            d1 = {d1.name};
            if ismember(['CycAvg_',In_FileName],d1)
                delete([Cyc_Path,filesep,'CycAvg_',In_FileName])
            end
        end
        return;
    elseif strcmp(opts{ind},'Shift Trigger')
        new_TrigShift = cellfun(@str2double,inputdlgcol('Trigger Shift (samples): ','Shift',[1 15],{num2str(info.TriggerShift2)},'on',1,[11 6.5 1.75 1.25]));
        if ~isempty(new_TrigShift)
            info.TriggerShift2 = round(new_TrigShift);
            [type,stims,t_snip,keep_inds,stim] = MakeCycAvg__alignCycles(info,Fs,ts,stim1);
            if size(keep_inds,2) > length(keep_tr)
                old_keep_tr = keep_tr;
                keep_tr = [old_keep_tr;true(1,size(keep_inds,2)-length(keep_tr))];
            elseif size(keep_inds,2) < length(keep_tr)
                keep_tr = keep_tr(1:size(keep_inds,2));
            end
            [filt,Data_calc,LE_V,RE_V,Data_cal,Data_In] = MakeCycAvg__filterTraces(type,filt_params_p,filt_params_v,te,ts,Data,keep_inds); 
            CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
            ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim_Pos,YLim_Vel,te,ts,t_snip,stim,stims,Data,Data_In,Data_cal,Data_calc,LE_V,RE_V,CycAvg,keep_inds,keep_tr);
        end
    elseif strcmp(opts{ind},'Filter Position')
        %Get new parameter values
        prompt = {['Median',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:',...
            ['Spline',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:',...
            ['Sav-Gol 1',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:',...
             ['Sav-Gol 2',newline,'L X:'],'R X:','L Y:','R Y:','L Z:','R Z:'};
        dlgtitle = 'Filter position';
        definput = strrep(cellfun(@(x) num2str(x,10),num2cell(filt_params_p),'UniformOutput',false),'NaN','');
        temp_filt_params_p = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 10],definput,'on',4,[11 6.5 3.5 4]));
        if ~isempty(temp_filt_params_p) %Didn't hit cancel
            filt_params_p = temp_filt_params_p;
            [filt,Data_calc,LE_V,RE_V,Data_cal,Data_In] = MakeCycAvg__filterTraces(type,filt_params_p,filt_params_v,te,ts,Data,keep_inds); 
            CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
            ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim_Pos,YLim_Vel,te,ts,t_snip,stim,stims,Data,Data_In,Data_cal,Data_calc,LE_V,RE_V,CycAvg,keep_inds,keep_tr);
        end
    elseif strcmp(opts{ind},'Filter Velocity')
        %Get new parameter values
        prompt = {'Irrlsmooth','Spline','Accel Thresh:','Median Filter:'};
        dlgtitle = 'Filter velocity';
        definput = strrep(cellfun(@(x) num2str(x,10),num2cell(filt_params_v),'UniformOutput',false),'NaN','');
        temp_filt_params_v = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 20],definput,'on',1,[11 7.25 2.25 3.25]));
        if ~isempty(temp_filt_params_v) %Didn't hit cancel
            filt_params_v = temp_filt_params_v;
            [filt,Data_calc,LE_V,RE_V,Data_cal,Data_In] = MakeCycAvg__filterTraces(type,filt_params_p,filt_params_v,te,ts,Data,keep_inds); 
            CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
            ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim_Pos,YLim_Vel,te,ts,t_snip,stim,stims,Data,Data_In,Data_cal,Data_calc,LE_V,RE_V,CycAvg,keep_inds,keep_tr);
        end
    elseif strcmp(opts{ind},'Select Cycles')
        [keep_tr,tf] = MakeCycAvg__selectCycles(type,keep_tr,t_snip,stims,LE_V,RE_V);
        while tf
            CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
            ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim_Pos,YLim_Vel,te,ts,t_snip,stim,stims,Data,Data_In,Data_cal,Data_calc,LE_V,RE_V,CycAvg,keep_inds,keep_tr);
            [keep_tr,tf] = MakeCycAvg__selectCycles(type,keep_tr,t_snip,stims,LE_V,RE_V);
        end
    elseif strcmp(opts{ind},'Set Y-axis Lim')
        %Get new parameter values
        prompt = {['Set Y-axis limits',newline,newline,'Position:',newline,newline,'Lower Limit:'],...
            'Upper Limit:',['Velocity:',newline,newline,'Lower Limit:'],'Upper Limit:'};
        dlgtitle = 'Y-axis Limits';
        definput = cellfun(@(x) num2str(x,10),num2cell([YLim_Pos,YLim_Vel]),'UniformOutput',false);
        out_nums = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 18],definput,'on',2,[11 8.25 3 2.25]));
        if ~isempty(out_nums)
            %Check to make sure they aren't reversed
            if out_nums(1)>out_nums(2) 
                YLim_Pos = [out_nums(2),out_nums(1)];
            else
                YLim_Pos = [out_nums(1),out_nums(2)];
            end
            %Check to make sure they aren't reversed
            if out_nums(3)>out_nums(4) 
                YLim_Vel = [out_nums(4),out_nums(3)];
            else
                YLim_Vel = [out_nums(3),out_nums(4)];
            end
        end
        if type==1
            set(ha(1),'YLim',YLim_Pos)
            set(ha(2),'YLim',YLim_Vel)
            set(ha(3),'YLim',YLim_Vel)
            set(ha(4),'YLim',YLim_Vel)
            set(ha(5),'YLim',YLim_Vel)           
        elseif type==2
            set(ha(1),'YLim',YLim_Pos)
            set(ha(2),'YLim',YLim_Vel)
        elseif type==3
            set(ha(1),'YLim',YLim_Vel)
            set(ha(2),'YLim',YLim_Vel)
            set(ha(3),'YLim',YLim_Vel)
        end
    else
        disp('This case not yet been defined.')
    end
    [ind,tf2] = nmlistdlg('PromptString','Select an action:',...
                       'SelectionMode','single',...
                       'ListSize',[100 100],...
                       'ListString',opts,...
                       'Position',[11,7.75,2,2.75]);  
    if tf2 == 0 %Treat this like an exit
        return;
    end
end
%% Save
CycAvg = MakeCycAvg__makeStruct(LE_V,RE_V,keep_tr,Data,Fs,t_snip,stims,info,filt,In_FileName);
savefig([Cyc_Path,filesep,'CycAvg_',In_FileName(1:end-4),'.fig'])
close;
%Save if you reached this point
MakeCycAvg__saveCycAvg(Cyc_Path,In_FileName,CycAvg);
end