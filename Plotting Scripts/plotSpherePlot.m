function plotSpherePlot(params)
Path = params.Path;
Cyc_Path = params.Cyc_Path;
code_Path = params.code_Path;
version = params.version;
Experimenter = params.Experimenter;
if isfield(params,'annot')
    annot = params.annot;
else
    annot = 1;
end
% Initialize
close all;
load('VNELcolors.mat','colors')
code_name = ['Plotting Scripts',filesep,'plotParamResults.m'];
warning('off')
sub_info = readtable('SubjectInfo.xlsx');
warning('on')
Subs = sub_info{:,1};
Ears = sub_info{:,2};
% Load table in question
res_file = extractfield(dir([Path,filesep,'*Results.mat']),'name')';
if isempty(res_file)
    disp('No table with cycle parameters found on this path. Creating one now.')
    rerun = ~strcmp(questdlg('If a parameter table already exists, use that one or rerun?','','Use existing table','Rerun','Rerun'),'Use existing table');
    MakeCycleSummaryTable(Path,Cyc_Path,rerun);
    res_file = extractfield(dir([Path,filesep,'*Results.mat']),'name')';
end
load(res_file{end},'all_results')
%Pick files to graph, full range to chose any
[indx,tf] = listdlg('ListString',all_results.File,...
    'PromptString','Pick the files to plot','ListSize',[400 300],...
    'SelectionMode','multiple');
if tf == 0
    return;
end
all_results2 = all_results(indx,:);
subject = all_results2.Subject{1};
hg = figure;
if annot
    annotation('textbox',[0 0 1 1],'String',[Path,newline,code_Path,filesep,...
        code_name,newline,...
        'VOGA',version,newline,Experimenter],'FontSize',5,...
        'EdgeColor','none','interpreter','none');
end
Function = 2;
plotstimaxis = 0;
plotelecaxis = 1;
normlen = 1;
stim_ear = Ears{ismember(Subs,subject)};
for i = 1:size(all_results2,1)
    load([Cyc_Path,filesep,all_results2.File{i}],'CycAvg')
    %Determine canal
    if any(contains(all_results2.AxisName{i},{'LP','RA','RALP'})) %RALP
        plot_colors = [colors.l_r;colors.r_r];
    elseif any(contains(all_results2.AxisName{i},{'LH','RH','LHRH'})) %LHRH
        plot_colors = [colors.l_z;colors.r_z];
    elseif any(contains(all_results2.AxisName{i},{'LA','RP','LARP'})) %LARP
        plot_colors = [colors.l_l;colors.r_l];
    else
        plot_colors = [0,0,0;0.5,0.5,0.5]; %black and gray
    end
    hg = MakeSpherePlot(CycAvg,hg,Function,plotstimaxis,plotelecaxis,normlen,plot_colors,stim_ear);
end
end