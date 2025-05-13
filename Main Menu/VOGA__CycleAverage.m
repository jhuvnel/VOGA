%% VOGA__CycleAverage
%
% This function expects that all relevant VOG files are already Segmented.
% From there, this program runs the MakeCycAvg script for filtering,
% cycle selection and parameterization. It requires MATLAB to be run in the
% directory of interest.
%
function VOGA__CycleAverage
if ~MakeFolders
    return;
end
%Expected folder names
Path = cd; Seg_Path = [Path,filesep,'Segmented Files']; Cyc_Path = [Path,filesep,'Cycle Averages'];
% Defaults for the experiment type and to start the loop
tf = 1; indx = 1; done = 0;
% Find unique combinations of experiment type, day, and goggle to analyze
ptab = assessProgress(Path);
all_exp_names = ptab.Segment; %Find all segments
for i = 1:length(all_exp_names)
    n_spl = strsplit(all_exp_names{i},'-');
    %Remove time
    if length(n_spl{3})>8 %If time is included with the date, remove the time
        n_spl{3} = n_spl{3}(1:8);
    elseif ~isnan(str2double(n_spl{4})) %If time is seperate, remove it
        n_spl(4) = [];
    end
    %Take from date to the next thing after experiment type
    all_exp_names(i) = join(n_spl(3:find(contains(n_spl,{'RotaryChair','eeVOR','aHIT','Manual'}))+1),'-');
end
exp_names = unique(all_exp_names);
if length(exp_names)>1 %Choose if there's more than one option
    [indx,tf] = nmlistdlg('PromptString','Select experiment types to analyze:','SelectionMode','multiple','ListSize',[350 300],'ListString',exp_names);
end
if tf == 0 %Selected cancel
    return;
end
exp_types = exp_names(indx);
%% run until the user hits cancel on analyzing a file
while(~done)
    % Select files that have the experiment type indicated above
    ptab = assessProgress(Path);
    rel_file = false(size(ptab,1),1);
    for i = 1:length(exp_types)
        rel_file(contains(ptab.Segment,exp_types{i}(1:8))&contains(ptab.Segment,exp_types{i}(9:end))) = 1;
    end
    ind = [find(~ptab{:,2}&~ptab{:,3});find(ptab{:,2}|ptab{:,3})]; %put unanalyzed files at the top
    % Change color of option depending on whether it's been analyzed
    % Unattempted = Black, Analyzed = Green, Not Analyzeable = Red
    font_color = repmat({'black'},length(ind),1);font_color(ptab{ind,2}) = {'green'};font_color(ptab{ind,3}) = {'red'};
    [indx,tf] = nmlistdlg('PromptString','Unattempted = Black, Analyzed = Green, Not Analyzeable = Red. Select a file to analyze:',...
        'SelectionMode','multiple','ListSize',[500 600],'ListString',strcat('<HTML><FONT color="',font_color,'">',table2cell(ptab(ind,1)),'</FONT></HTML>'));
    if tf % Run unless the user selects cancel
        in_opts = {}; %Default option to run MakeCycAvg
        if length(indx)>1 
            in_opts = {'Auto Rerun'}; %reanalyze files that already exist or analyze multiple files in a row          
        end
        for j = 1:length(indx)
            fname = ptab{ind(indx(j)),1}{:};
            load([Seg_Path,filesep,fname],'Data');
            [CycAvg,analyzed] = MakeCycAvg(Data,Cyc_Path,in_opts,1);
            if analyzed                    
                MakeCycAvg__saveCycAvg(Cyc_Path,strrep(CycAvg.name,'CycAvg_',''),CycAvg,analyzed,1);
                if ~strcmp(CycAvg.name,['CycAvg_',fname]) %Had a different name, delete file
                    delete([Cyc_Path,filesep,'CycAvg_',fname])
                end
            end
        end
    else
        disp('Operation Ended.')
        done = true;
    end
end
end