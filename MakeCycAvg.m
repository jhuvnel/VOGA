%% MakeCycAvg.m
% This function was created to automate the data analysis pipeline for the
% LDVOG data sets as much as possible.
% Cycle average before filtering

function [CycAvg,analyzed] = MakeCycAvg(Data,Cyc_Path,from_file)
%% General Setup
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
line_wid.norm = 0.5;
line_wid.bold = 2;
%try
    %% Identify relevant traces
    [Data,info,Fs,te,ts,stim1,type,detec_head,filt1,traces_vel1] = MakeCycAvg_startProcess(Data,info,filt1,all_traces);
    %% Reanalyze file with existing params
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
            info.TriggerShift2 = CycAvg2.info.TriggerShift2;
            [filt,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc] = MakeCycAvg__filterTraces(filt,keep_inds,te,ts,t_snip,stim,stims,Data,t_interp);
            CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);
            ha = MakeCycAvg__plotFullCycAvg([],type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);
            % Create to Save
            CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);
            CycAvg = ParameterizeCycAvg(CycAvg);
            filt_params.filt1 = filt;
            filt_params.YLim = YLim;
            VOGA__saveLastUsedParams(filt_params)
            analyzed = 1;
        else
            disp('Not a compatible CycAvg file.')
            CycAvg = Data;
            analyzed = 0;
        end
        return;
    end
    %% Once analyzeable, here is the while loop they stay in until saving or exiting
    %You can change the order/existance of these options without ruining
    %anything because the comparrisons are all string based
    opts = {'Set Y-axis Lim','Choose Coordinates',...
        'Filter Position','Filter Velocity','Manual QPR',...
        'Select Cycles','Shift Trigger','Shorten Segment',...
        'Load from File','Not Analyzeable','Start Over','Save'};
    sel = 'Start Over'; %Run the start procedure first
    while ~strcmp(sel,'Save') %Run until it's ready to save or just hopeless
        switch sel
            case 'Start Over'
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
            case 'Shift Trigger'
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
            case 'Filter Position'
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
            case 'Filter Velocity'
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
            case 'Manual QPR'
                clc;
                disp('Select the start and end time point to linearly interpolate over on the cycle average graph.')
                good_rng = 'Redo';
                if type == 1
                    cyc_ax = ha(4);
                elseif type == 2
                    cyc_ax = ha(2);
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
            case 'Select Cycles'
                [keep_tr,ha,tf] = MakeCycAvg__selectCycles(ha,type,keep_tr,Data_cyc,screen_size,traces_vel);
                while tf
                    CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);
                    ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);
                    [keep_tr,ha,tf] = MakeCycAvg__selectCycles(ha,type,keep_tr,Data_cyc,screen_size,traces_vel);
                end
            case 'Load from File'
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
            case 'Shorten Segment'
                [Data,good_rng] = MakeCycAvg__shortenSegment(ha,te,Data,screen_size);
                if strcmp(good_rng,'Keep') %Did shorten segment
                    save([strrep(Cyc_Path,'Cycle Averages','Segmented Files'),filesep,fname],'Data') 
                    [Data,info,Fs,te,ts,stim1,type,detec_head,filt1,traces_vel1] = MakeCycAvg_startProcess(Data,info,filt1,all_traces);
                    [stim,t_snip,stims,keep_inds,detec_tr] = MakeCycAvg__alignCycles(info,Fs,ts,stim1,detec_head);
                    [filt,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc] = MakeCycAvg__filterTraces(filt1,keep_inds,te,ts,t_snip,stim,stims,Data,t_interp);
                    CycAvg = MakeCycAvg__makeStruct(fname,info,Fs,filt,keep_tr,detec_tr,t_interp,Data,Data_pos,Data_pos_filt,Data_vel,Data_vel_filt,Data_cyc);
                    ha = MakeCycAvg__plotFullCycAvg([],type,colors,line_wid,YLim.Pos,YLim.Vel,traces_pos,traces_vel,CycAvg);
                end
            case 'Set Y-axis Lim'
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
            case 'Choose Coordinates'
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
            case 'Not Analyzeable'
                CycAvg = Data;
                analyzed = 0;
                return;
        end
        [ind,tf2] = nmlistdlg('PromptString','Select an action:',...
            'SelectionMode','single',...
            'ListSize',[150 180],...
            'ListString',opts,...
            'Position',[screen_size(3)-4,screen_size(4)-3.75,3,3.75]);
        sel = opts{ind};
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
% catch e %If you hit an error, mark it as unanalyzeable
%     disp(['Error found in file ',fname])
%     fprintf(1,'The identifier was:\n%s',e.identifier);
%     fprintf(1,'There was an error! The message was:\n%s',e.message);
%     CycAvg = Data;
%     analyzed = 0;
% end
end