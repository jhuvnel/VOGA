%% MakeCycAvg.m
% This function was created to automate the data analysis pipeline for the
% LDVOG data sets as much as possible.
% Cycle average before filtering

function [CycAvg,analyzed] = MakeCycAvg(Data,Cyc_Path,in_opts,has_fig)
%% General Setup
%Loop parameters
opts = {'Set Plot View','Filter Position','Filter Velocity',...
    'Select Cycles','Advanced','Start Over','Save'};
advanced_opts = {'Shift Trigger','Shorten Segment',...
    'Load from Selected File','Manual QPR','Not Analyzeable',...
    'Auto Rerun'};
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
filt_params = VOGA__saveLastUsedParams;
filt1 = filt_params.filt1;
YLim = filt_params.YLim; 
%For speed during autoscan analysis, do not make cycle average figures
plot_info.screen_size = [];
if has_fig
    % Initialize Figure
    close all;
    fig = figure('Units','normalized','Position',[0 0 1 1],'Units','inches');
    screen_size = fig.Position;
    fig.Position = screen_size - [-0.5 -0.5 4.5 1.5];
    plot_info.screen_size = screen_size;
    % Title
    fig_title = strrep(strrep(strrep(fname,'_',' '),'-',' '),'.mat','');
    fig_title(strfind(fig_title,'['):strfind(fig_title,']')) = strrep(fig_title(strfind(fig_title,'['):strfind(fig_title,']')),' ','-');
    annotation('textbox',[0 .9 1 .1],'String',fig_title,'FontSize',14,'HorizontalAlignment','center','EdgeColor','none');
end
% Plot defaults
all_traces = {'LX','RX','LY','RY','LZ','RZ','LLARP','RLARP','LRALP','RRALP'};
traces_pos1 = all_traces(1:6); %XYZ pos
if contains(fname,{'X','Y'})
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
if contains(fname,'ESC3')||(contains(fname,'Impulse')&&~contains(fname,{'GNO','ESC'}))
    type = 4;
elseif contains(fname,{'GNO','ESC'})
    type = 3;    
elseif contains(fname,{'Activation','Step'})
    type = 2;
else
    type = 1;
end
Data.info.type = type;
%% The Main Loop
while ~strcmp(sel,'Save') %Run until it's ready to save or just hopeless
    if strcmp(sel,'Advanced') %Give the advanced menu and run that selection
        [ind2,tf2] = nmlistdlg('PromptString','Select an action:',...
            'SelectionMode','single','ListSize',[150 125],'ListString',advanced_opts,...
            'Position',[screen_size(3)-4,screen_size(4)-3.75,3,3.75]);
        if tf2
            sel = advanced_opts{ind2};
        end
    end
    switch sel
        case 'Start Over'
            plot_info.traces_pos = traces_pos1;
            plot_info.traces_vel = traces_vel1;
            plot_info.YLim = YLim;
            Data.info.TriggerShift2 = 0;
            Data = MakeCycAvg__alignCycles(Data); % Cycle Align
            filt = filt1;
            filt.t_interp = [];
            filt.keep_tr = true(1,size(Data.keep_inds,2));
            filt = MakeCycAvg__autoFilter(Data,filt,plot_info);
            CycAvg = MakeCycAvg__filterTraces(Data,filt);
            ha = [];
            if has_fig %Only should not happen during autoscan analysis
                ha = MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
            end
            [CycAvg,filt] = MakeCycAvg__selectCycles(ha,CycAvg,plot_info,'Automatic'); 
        case {'Filter Position','Filter Velocity'}
            %Load either filt.pos or filt.vel and make a prompt with each
            %of the filters and traces in the table.
            pos_vel = char(lower(extract(sel,"Pos"|"Vel")));
            tab = filt.(pos_vel); 
            traces = strrep(strrep(tab.Properties.RowNames,'ALP',''),'ARP','');
            filters = tab.Properties.VariableNames;
            prompt = reshape([strcat(filters,[newline,traces{1}]);...
                repmat(traces(2:end),1,length(filters))],[],1)';
            position = [screen_size(3)-(0.5+0.75*length(filters)),...
                screen_size(4)-(1.0+0.5*length(traces)),...
                0.75*length(filters),0.5+0.5*length(traces)];
            definput = strrep(cellfun(@(x) num2str(x,10),...
                table2cell(filt.(pos_vel)),'UniformOutput',false),'NaN','');
            temp_filt = cellfun(@str2double,inputdlgcol(prompt,sel,...
                [1 10],definput,'on',length(prompt)/length(traces),...
                position,{'Done','Refilter'}));
            %Keep running until the user selects "Done" (the Cancel Button renamed)
            while ~isempty(temp_filt)
                filt.(pos_vel){:,:} = reshape(temp_filt,length(traces),[]);
                CycAvg = MakeCycAvg__filterTraces(Data,filt);
                CycAvg = MakeCycAvg__selectCycles(ha,CycAvg,plot_info,'Restart');
                [CycAvg,filt] = MakeCycAvg__selectCycles(ha,CycAvg,plot_info,'Automatic');
                ha = MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
                definput = strrep(cellfun(@(x) num2str(x,10),...
                    table2cell(filt.(pos_vel)),'UniformOutput',false),'NaN','');
                temp_filt = cellfun(@str2double,inputdlgcol(prompt,sel,...
                    [1 10],definput,'on',length(prompt)/length(traces),...
                    position,{'Done','Refilter'}));
            end
        case 'Select Cycles'
            [CycAvg,filt] = MakeCycAvg__selectCycles(ha,CycAvg,plot_info);
        case 'Set Plot View'
            %Get new parameter values
            prompt = {['Set Y-axis limits',newline,newline,'Position:',newline,newline,'Lower Limit:'],...
                'Upper Limit:',['Velocity:',newline,newline,'Lower Limit:'],'Upper Limit:'};
            dlgtitle = 'Y-axis Limits';
            definput = cellfun(@(x) num2str(x,10),num2cell([plot_info.YLim.Pos,plot_info.YLim.Vel]),'UniformOutput',false);
            out_nums = cellfun(@str2double,inputdlgcol(prompt,dlgtitle,[1 18],definput,'on',2,[screen_size(3)-4 screen_size(4)-4.25 3 2.25]));
            if ~isempty(out_nums)
                %Check to make sure they aren't reversed
                plot_info.YLim.Pos = sort([out_nums(1),out_nums(2)]);
                plot_info.YLim.Vel = sort([out_nums(3),out_nums(4)]);
            end
            if type ~= 3
                [ind3,tf] = nmlistdlg('PromptString','Select position traces:',...
                    'InitialValue',find(ismember(all_traces,plot_info.traces_pos)),...
                    'ListSize',[100 150],'ListString',all_traces,...
                    'Position',[screen_size(3)-4,screen_size(4)-5,2,4]);
                if tf
                    plot_info.traces_pos = all_traces(ind3);
                end
            end
            [ind2,tf] = nmlistdlg('PromptString','Select velocity traces:',...
                'InitialValue',find(ismember(all_traces,plot_info.traces_vel)),...
                'ListSize',[100 150],'ListString',all_traces,...
                'Position',[screen_size(3)-4,screen_size(4)-5,2,4]);
            if tf
                plot_info.traces_vel = all_traces(ind2);
            end
            ha = MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
        case 'Shift Trigger'
            new_TrigShift = cellfun(@str2double,inputdlgcol('Trigger Shift (samples): ','Shift',[1 15],{num2str(Data.info.TriggerShift2)},'on',1,[screen_size(3)-4 screen_size(4)-1.25 1.75 1.25]));
            if ~isempty(new_TrigShift)
                Data.info.TriggerShift2 = round(new_TrigShift);
                Data = MakeCycAvg__alignCycles(Data);
                if size(Data.keep_inds,2) > length(CycAvg.keep_tr)
                    filt.keep_tr = [CycAvg.keep_tr;true(1,size(Data.keep_inds,2)-length(CycAvg.keep_tr))];
                elseif size(Data.keep_inds,2) < length(CycAvg.keep_tr)
                    filt.keep_tr = CycAvg.keep_tr(1:size(Data.keep_inds,2));
                end
                [CycAvg,filt] = MakeCycAvg__filterTraces(Data,filt);
                ha = MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
            end                
        case 'Shorten Segment'
            [Data,good_rng] = MakeCycAvg__shortenSegment(ha,Data,plot_info.screen_size);
            if strcmp(good_rng,'Keep') %Did shorten segment
                save([strrep(Cyc_Path,'Cycle Averages','Segmented Files'),filesep,fname],'Data')
                Data = MakeCycAvg__alignCycles(Data); % Cycle Align
                [CycAvg,filt] = MakeCycAvg__filterTraces(Data,filt);
                ha = MakeCycAvg__plotFullCycAvg([],CycAvg,plot_info);
            end
        case 'Load from Selected File'
            cyc_files = extractfield(dir([Cyc_Path,filesep,'*.mat']),'name');
            [indx,tf] = nmlistdlg('PromptString','Select an analyzed file to use:',...
                    'SelectionMode','single','ListSize',[500 600],...
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
                    [CycAvg,filt] = MakeCycAvg__filterTraces(Data,filt);
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
                [CycAvg,filt] = MakeCycAvg__filterTraces(Data,filt);
                ha = MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
            end
            hold off    
        case 'Not Analyzeable'
            CycAvg = Data;
            analyzed = 0;
            return;
        case 'Auto Rerun'
            a = load([Cyc_Path,filesep,'CycAvg_',fname]);
            CycAvg2 = a.CycAvg;
            if isfield(CycAvg2.info,'TriggerShift2')
                 Data.info.TriggerShift2 = Data.info.TriggerShift2;     
            end
            Data = MakeCycAvg__alignCycles(Data);
            filt = CycAvg2.filt;
            filt.keep_tr = CycAvg2.keep_tr;
            if isfield(CycAvg2,'t_interp')
                filt.t_interp = CycAvg2.t_interp;
            end
            [CycAvg,filt] = MakeCycAvg__filterTraces(Data,filt);
            MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info); 
            break;
    end    
    if ~isempty(in_opts) %Run the next input command
        sel = in_opts{1};
        in_opts{1} = [];
    else
        [ind,tf2] = nmlistdlg('PromptString','Select an action:',...
            'SelectionMode','single','ListSize',[150 180],'ListString',opts,...
            'Position',[screen_size(3)-4,screen_size(4)-3.75,3,3.75]);
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
filt_params.filt1.pos = filt.pos;
filt_params.filt1.vel = filt.vel;
filt_params.YLim = plot_info.YLim;
VOGA__saveLastUsedParams(filt_params);
analyzed = 1;
end