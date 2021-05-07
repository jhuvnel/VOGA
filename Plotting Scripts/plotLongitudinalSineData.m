%% Load in files
load('VNELColors.mat','colors')
if ispc
    drive_name = 'Z:';
else
    drive_name = '/Volumes/vnelhuman$/MVI';
end
out_path = [drive_name,filesep,'DATA SUMMARY',filesep,'IN PROGRESS',filesep,'VOR - Rotary Chair'];
res_files = extractfield(dir([out_path,filesep,'*RotaryChairResults.mat']),'name');
load([out_path,filesep,res_files{end}],'all_results')
disp(['File: ',out_path,filesep,res_files{end}])
%% See what items are in the file
tab_labs = all_results.Properties.VariableNames;
desc_labs = tab_labs(find(contains(tab_labs,'File'))+1:find(contains(tab_labs,'Cycles'))-1);
desc_labs(contains(desc_labs,'Date')) = [];
for i = 1:length(desc_labs)
    column = all_results.(desc_labs{i});
    out_lab = strrep(strrep(desc_labs{i},'(','_'),')','');
    if isnumeric(column)
        column(isnan(column)) = [];
        contents.(out_lab) = unique(column);
    elseif isnumeric(column{1})
        column = cell2mat(column); 
        contents.(out_lab) = unique(column,'rows');
    else
        contents.(out_lab) = unique(column);
    end    
    disp([desc_labs{i},': '])
    disp(contents.(out_lab))
end
% Select outputs to plot: Magnitude, Gain, Phase, Misalignment, Disc
%all_results = all_results(contains(all_results.Type,'Sine'),:); %Constrain to sinusoids