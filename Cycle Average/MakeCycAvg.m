%% MakeCycAvg.m
% This function was created to automate the data analysis pipeline for the
% LDVOG data sets as much as possible.
% Cycle average before filtering

function [CycAvg,analyzed] = MakeCycAvg(Data,Cyc_Path,in_opts,has_fig)
%% General Setup
%Loop parameters
opts = {'Set Y-Limit','Set Velocity View','Filter Position',...
    'Filter Velocity','Select Cycles','Advanced','Start Over','Save'};
advanced_opts = {'Set Position View','Shift Trigger','Shorten Segment',...
    'Manual QPR','Load from Selected File','Not Analyzeable','Auto Rerun'};
sel = 'Start Over'; %Run the start procedure first
% Input handling
if nargin < 4
    has_fig = 1;
end
if nargin < 3
    in_opts = {};
end
if ~isempty(in_opts)&&any(~ismember(in_opts,[opts,advanced_opts]))
    disp(in_opts(~ismember(in_opts,[opts,advanced_opts])))
    error('Invalid CycAvg commands.')
end
% Set colors and line widths
load('VNELcolors.mat','colors')
colors.cyc_keep = [0.85 0.85 0.85]; % Fill colors for cycle selection
colors.cyc_rm = [1 1 1];
plot_info.colors = colors;
plot_info.line_wid.norm = 0.5;
plot_info.line_wid.bold = 2;
% Set Experimentor/version
VOGA_VerInfo = rows2vars(readtable([userpath,filesep,'VOGA_VerInfo.txt'],'ReadVariableNames',false,'ReadRowNames',true));
Data.info.Analyzer = VOGA_VerInfo.Experimenter{:};
Data.info.ver = VOGA_VerInfo.Version{:};
% Get or set the file name
if ~isfield(Data.info,'name')
    Data.info.name = [Data.info.subject,'-',Data.info.visit,'-',Data.info.exp_date,'-',Data.info.goggle_ver,'-',Data.info.dataType,'.mat'];
end
fname = Data.info.name;
% Load default filters
filt_params = SaveLastUsedParams;
%For speed during autoscan analysis, do not make cycle average figures
buffer_pix = [25,50,450,100]; %left, bottom, right, and top # of pixel buffer
plot_info.screen_size = get(0,"ScreenSize");
plot_info.fig_space = [buffer_pix(1:2),plot_info.screen_size(3:4)-buffer_pix(3:4)-buffer_pix(1:2)];
plot_info.menu_space = [sum(plot_info.fig_space([1,3]))+buffer_pix(1),...
    plot_info.fig_space(2),buffer_pix(3)-2*buffer_pix(1),plot_info.fig_space(4)];
if has_fig
    % Initialize Figure
    close all;
    figure('Color','w','Position',plot_info.fig_space)
    % Title
    fig_title = strrep(strrep(strrep(fname,'_',' '),'-',' '),'.mat','');
    fig_title(strfind(fig_title,'['):strfind(fig_title,']')) = strrep(fig_title(strfind(fig_title,'['):strfind(fig_title,']')),' ','-');
    annotation('textbox',[0 .9 1 .1],'String',fig_title,'FontSize',14,'HorizontalAlignment','center','EdgeColor','none');
end
% Plot defaults
all_traces = {'LX','RX','LY','RY','LZ','RZ','LLARP','RLARP','LRALP','RRALP'};
traces_pos1 = all_traces(1:6); %XYZ pos
if contains(fname,{'X','Y','Activation'})
    traces_vel1 = all_traces(1:6);
elseif contains(fname,{'GNO','ESC3'})
    traces_vel1 = all_traces(4+2*find(cellfun(@(x) contains(fname,x),{{'LH','RH'};{'LA','RP'};{'RA','LP'}})));
elseif contains(fname,'ESC')
    traces_vel1 = all_traces(3+2*find(cellfun(@(x) contains(fname,x),{{'LH','RH'};{'LA','RP'};{'RA','LP'}})));
else
    traces_vel1 = all_traces(5:end);
end
% Assign Type (Position, Cycle Averaging, Impulse)
% Type 1: Yes, Yes, No (Pulse Train, Sinusoids)
% Type 2: Yes, No, No (Velocity Steps, Activation)
% Type 3: No, Yes, Yes (Impulse from GNO and old ESC)
% Type 4: Yes, Yes, Yes (Impulse from LDVOG and new ESC)
type_ord = [3,4,2,1]; %Decision tree of which type it is
type_logic = [contains(fname,{'GNO','ESC'})&&~contains(fname,'ESC3'),...
    contains(fname,'Impulse'),contains(fname,{'Activation','Step'}),true];
type = type_ord(find(type_logic,1,'first'));
Data.info.type = type;
TriggerShift = 0;
if isfield(Data.info,'TriggerShift2')
    TriggerShift = Data.info.TriggerShift2;
end
%% The Main Loop
while ~strcmp(sel,'Save') %Run until it's ready to save or just hopeless
    if strcmp(sel,'Advanced') %Give the advanced menu and run that selection
        [ind2,tf2] = nmlistdlg('PromptString','Select an action:','ListString',advanced_opts,...
            'SelectionMode','single','ListSize',round(0.9*plot_info.menu_space(3:4)),...
            'Units','pixels','Position',plot_info.menu_space);
        if tf2
            sel = advanced_opts{ind2};
        end
    end
    switch sel
        case 'Start Over'            
            plot_info.traces_pos = traces_pos1;
            plot_info.traces_vel = traces_vel1;
            plot_info.YLim.Pos = [NaN NaN];
            plot_info.YLim.Vel = [NaN NaN];
            Data.info.TriggerShift2 = TriggerShift;
            Data = MakeCycAvg__alignCycles(Data); % Cycle Align            
            filt = MakeCycAvg__autoFilter(Data,filt_params,plot_info); %First pass at filter params
            filt.t_interp = [];
            filt.keep_tr = true(1,size(Data.keep_inds,2));
            [CycAvg,filt,plot_info] = MakeCycAvg__filterTraces(Data,filt,plot_info); 
            if length(filt.keep_tr)>5
                [CycAvg,filt] = MakeCycAvg__selectCycles([],CycAvg,plot_info,'Automatic');
            end
            if has_fig 
                ha = MakeCycAvg__plotFullCycAvg([],CycAvg,plot_info);
            end           
        case {'Filter Position','Filter Velocity'}
            %Load either filt.pos or filt.vel and make a prompt with each
            %of the filters and traces in the table.
            pos_vel = char(lower(extract(sel,"Pos"|"Vel")));
            tab = filt.(pos_vel); 
            traces = strrep(strrep(tab.Properties.RowNames,'ALP',''),'ARP','');
            filters = tab.Properties.VariableNames;
            prompt = reshape([strcat(filters,[newline,traces{1}]);...
                repmat(traces(2:end),1,length(filters))],[],1)';
            position = plot_info.menu_space;
            definput = strrep(cellfun(@(x) num2str(x,10),...
                table2cell(filt.(pos_vel)),'UniformOutput',false),'NaN','');
            temp_filt = cellfun(@str2double,inputdlgcol(prompt,sel,...
                [1 10],definput,'on',length(prompt)/length(traces),...
                position,{'Done','Refilter'},'pixels'));
            %Keep running until the user selects "Done" (the Cancel Button renamed)
            while ~isempty(temp_filt)
                filt.(pos_vel){:,:} = reshape(temp_filt,length(traces),[]);
                CycAvg = MakeCycAvg__filterTraces(Data,filt,plot_info);
                CycAvg = MakeCycAvg__selectCycles(ha,CycAvg,plot_info,'Restart');
                [CycAvg,filt] = MakeCycAvg__selectCycles(ha,CycAvg,plot_info,'Automatic');
                ha = MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
                definput = strrep(cellfun(@(x) num2str(x,10),...
                    table2cell(filt.(pos_vel)),'UniformOutput',false),'NaN','');
                temp_filt = cellfun(@str2double,inputdlgcol(prompt,sel,...
                    [1 10],definput,'on',length(prompt)/length(traces),...
                    position,{'Done','Refilter'},'pixels'));
            end
        case 'Select Cycles'
            [CycAvg,filt] = MakeCycAvg__selectCycles(ha,CycAvg,plot_info);
        case 'Set Y-Limit'
            %Get new parameter values
            prompt = {['Set Y-axis limits',newline,newline,'Position:',newline,newline,'Lower Limit:'],...
                'Upper Limit:',['Velocity:',newline,newline,'Lower Limit:'],'Upper Limit:'};
            dlgtitle = 'Y-axis Limits';
            definput = cellfun(@(x) num2str(x,10),num2cell([plot_info.YLim.Pos,plot_info.YLim.Vel]),'UniformOutput',false);
            position = [plot_info.menu_space(1),sum(plot_info.menu_space([2,4]))-300,plot_info.menu_space(3),300];
            out_nums = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 18],definput,'on',2,position,[],'pixels'));
            if ~isempty(out_nums)
                %Check to make sure they aren't reversed
                plot_info.YLim.Pos = sort([out_nums(1),out_nums(2)]);
                plot_info.YLim.Vel = sort([out_nums(3),out_nums(4)]);
            end
            ha = MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
        case 'Set Position View'
            if type ~= 3
                [ind3,tf] = nmlistdlg('PromptString','Select position traces:',...
                    'InitialValue',find(ismember(all_traces,plot_info.traces_pos)),...
                    'ListString',all_traces,'ListSize',round(0.8*plot_info.menu_space(3:4)),...
                    'Units','pixels','Position',plot_info.menu_space);
                if tf
                    plot_info.traces_pos = all_traces(ind3);
                end
            end
            ha = MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
        case 'Set Velocity View'
            [ind2,tf] = nmlistdlg('PromptString','Select velocity traces:',...
                'InitialValue',find(ismember(all_traces,plot_info.traces_vel)),...
                'ListString',all_traces,'ListSize',round(0.8*plot_info.menu_space(3:4)),...
                'Units','pixels','Position',plot_info.menu_space);
            if tf
                plot_info.traces_vel = all_traces(ind2);
            end
            ha = MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
        case 'Shift Trigger'
            position = [plot_info.menu_space(1),sum(plot_info.menu_space([2,4]))-300,plot_info.menu_space(3),300];
            new_TrigShift = cellfun(@str2double,inputdlgcol('Trigger Shift (samples): ',...
                'Shift',[1 15],{num2str(Data.info.TriggerShift2)},'on',1,position,[],'pixels'));
            while ~isempty(new_TrigShift)
                Data.info.TriggerShift2 = round(new_TrigShift);
                Data = MakeCycAvg__alignCycles(Data);
                if size(Data.keep_inds,2) > length(CycAvg.keep_tr)
                    filt.keep_tr = [CycAvg.keep_tr,true(1,size(Data.keep_inds,2)-length(CycAvg.keep_tr))];
                elseif size(Data.keep_inds,2) < length(CycAvg.keep_tr)
                    filt.keep_tr = CycAvg.keep_tr(1:size(Data.keep_inds,2));
                end
                [CycAvg,filt] = MakeCycAvg__filterTraces(Data,filt,plot_info);
                ha = MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
                position = [plot_info.menu_space(1),sum(plot_info.menu_space([2,4]))-300,plot_info.menu_space(3),300];
                new_TrigShift = cellfun(@str2double,inputdlgcol('Trigger Shift (samples): ',...
                    'Shift',[1 15],{num2str(Data.info.TriggerShift2)},'on',1,position,[],'pixels'));
            end                
        case 'Shorten Segment'
            [Data,good_rng] = MakeCycAvg__shortenSegment(ha,Data,plot_info);
            if strcmp(good_rng,'Keep') %Did shorten segment
                save([strrep(Cyc_Path,'Cycle Averages','Segmented Files'),filesep,fname],'Data')
                Data = MakeCycAvg__alignCycles(Data); % Cycle Align
                filt.keep_tr = true(1,size(Data.keep_inds,2));
                CycAvg = MakeCycAvg__filterTraces(Data,filt,plot_info);
                ha = MakeCycAvg__plotFullCycAvg([],CycAvg,plot_info);
                [CycAvg,filt] = MakeCycAvg__selectCycles(ha,CycAvg,plot_info,'Automatic');
            end
        case 'Load from Selected File'
            cyc_files = extractfield(dir([Cyc_Path,filesep,'*.mat']),'name');
            [indx,tf] = nmlistdlg('PromptString','Select an analyzed file to use:',...
                    'SelectionMode','single','ListSize',round(0.9*plot_info.menu_space(3:4)),...
                    'Units','pixels','Position',plot_info.menu_space,...
                    'InitialValue',find(ismember(cyc_files,['CycAvg_',fname])),...
                    'ListString',cyc_files);
            if tf
                a = load([Cyc_Path,filesep,cyc_files{indx}]);
                CycAvg2 = a.CycAvg;
                if isfield(CycAvg2.info,'TriggerShift2')
                     Data.info.TriggerShift2 = Data.info.TriggerShift2;
                     Data = MakeCycAvg__alignCycles(Data);
                     filt.keep_tr = true(1,size(Data.keep_inds,2));
                end
                if length(CycAvg2.keep_tr) == length(filt.keep_tr) %Only if they are the same size
                    filt = CycAvg2.filt;
                    filt.keep_tr = CycAvg2.keep_tr;
                    if isfield(CycAvg2,'t_interp')
                        filt.t_interp = CycAvg2.t_interp;
                    end
                    [CycAvg,filt] = MakeCycAvg__filterTraces(Data,filt,plot_info);
                    ha = MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
                else
                    disp('Not a compatible CycAvg file: unequal number of cycles detected.')
                end
            end   
        case 'Manual QPR'
            clc;
            uiwait(msgbox('Select the start and end time point to linearly interpolate over on the cycle average graph.'))
            good_rng = 'Redo';
            axes(ha(4))
            hold on
            while strcmp(good_rng,'Redo')                
                [t_spline,~] = ginput(2);
                h1 = xline(t_spline,'LineWidth',10);
                good_rng = questdlg('Keep or redo the range?','','Keep','Redo','Redo');
                delete(h1)
            end
            if strcmp(good_rng,'Keep')
                filt.t_interp = find(Data.t_snip>=t_spline(1)&Data.t_snip<=t_spline(2));
                [CycAvg,filt] = MakeCycAvg__filterTraces(Data,filt,plot_info);
                ha = MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
            end
            hold off    
        case 'Not Analyzeable'
            CycAvg = Data;
            analyzed = 0;
            return;
        case 'Auto Rerun'
            if isfile([Cyc_Path,filesep,'CycAvg_',fname])
                a = load([Cyc_Path,filesep,'CycAvg_',fname]);
                CycAvg2 = a.CycAvg;
                if isfield(CycAvg2.info,'TriggerShift2')
                     Data.info.TriggerShift2 = Data.info.TriggerShift2;     
                end
                Data = MakeCycAvg__alignCycles(Data);
                pos_filt = filt.pos.Properties.VariableNames;
                vel_filt = filt.vel.Properties.VariableNames;
                %Remove unused filters
                CycAvg2.filt.pos = CycAvg2.filt.pos(:,any(~isnan(CycAvg2.filt.pos{:,:})));
                CycAvg2.filt.vel = CycAvg2.filt.vel(:,any(~isnan(CycAvg2.filt.vel{:,:})));
                %Adjust for new Savitsky Golay nomenclature
                if contains(CycAvg2.filt.pos.Properties.VariableNames,'sgolay2')
                    CycAvg2.filt.pos.sgolay = CycAvg2.filt.pos.sgolay2;
                    CycAvg2.filt.pos = removevars(CycAvg2.filt.pos,["sgolay1","sgolay2"]);
                end
                if contains(CycAvg2.filt.vel.Properties.VariableNames,'sgolay2')
                    CycAvg2.filt.vel.sgolay = CycAvg2.filt.vel.sgolay2;
                    CycAvg2.filt.vel = removevars(CycAvg2.filt.vel,["sgolay1","sgolay2"]);
                end
                pos_filt2 = CycAvg2.filt.pos.Properties.VariableNames;
                vel_filt2 = CycAvg2.filt.vel.Properties.VariableNames;
                %Only run if right number of cycles detected and filter types match
                if size(Data.keep_inds,2)==length(CycAvg2.keep_tr)&&...
                        all(ismember(pos_filt2,pos_filt))&&all(ismember(vel_filt2,vel_filt)) 
                    filt.pos(:,pos_filt2) = CycAvg2.filt.pos(:,pos_filt2);
                    filt.vel(:,vel_filt2) = CycAvg2.filt.vel(:,vel_filt2);
                    filt.keep_tr = CycAvg2.keep_tr;
                    if isfield(CycAvg2,'t_interp')
                        filt.t_interp = CycAvg2.t_interp;
                    end
                    [CycAvg,filt] = MakeCycAvg__filterTraces(Data,filt,plot_info);
                    MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info); 
                    break;
                end
            end
    end    
    if ~isempty(in_opts) %Run the next input command
        sel = in_opts{1};
        in_opts(1) = [];
    else
        [ind,tf2] = nmlistdlg('PromptString','Select an action:','ListString',opts,...
            'SelectionMode','single','ListSize',round(0.9*plot_info.menu_space(3:4)),...
            'Units','pixels','Position',plot_info.menu_space);
        if tf2 == 0 %Treat this like an exit
            CycAvg = [];
            analyzed = 0;
            return;
        end
        sel = opts{ind};
    end    
end
%% Create to Save
CycAvg = ParameterizeCycAvg(CycAvg);
filt_params.filt.pos = filt.pos;
filt_params.filt.vel = filt.vel;
filt_params.YLim = plot_info.YLim;
SaveLastUsedParams(filt_params);
analyzed = 1;
end