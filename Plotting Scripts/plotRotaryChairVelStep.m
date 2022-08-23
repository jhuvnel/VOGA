%% Plot the Magnitude and Time Constant Over Time for One Subject

%Load, find rotary chair exponentials and error handle
tab_name = extractfield(dir([cd,filesep,'*Results.mat']),'name');
if isempty(tab_name)
    error(['No subject-wide results table exists yet in ',cd])
end
load(tab_name{end},'all_results')
all_results = all_results(contains(all_results.Experiment,'RotaryChair')&contains(all_results.Type,'Exp'),:);
if isempty(all_results)
    error(['No rotary chair velocity step experiments detected in the results table ',tab_name{end}])
end
%For the purpose of this figure, only look at +/- 240dps
all_results(~ismember(all_results.('Amplitude(dps)'),[240,-240]),:) = [];
