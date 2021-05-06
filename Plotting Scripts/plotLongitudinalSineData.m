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
    try
        contents.(desc_labs{i}) = unique(all_results.(desc_labs{i}));
    catch %array of numbers
        contents.(desc_labs{i}) = unique(cell2mat(all_results.(desc_labs{i})),'rows');
    end
    disp([desc_labs{i},': '])
    disp(contents.(desc_labs{i}))
end
% Select outputs to plot: Magnitude, Gain, Phase, Misalignment, Disc
%all_results = all_results(contains(all_results.Type,'Sine'),:); %Constrain to sinusoids