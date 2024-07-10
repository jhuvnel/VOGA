%% Combine GNO Tables
% Takes the outputs from MakeGNOSummaryTable
function all_results = CombineGNOTables(VOGA_VerInfo,rerun)
% Handle inputs
if nargin < 2
    rerun = 1;
else
    rerun = ~strcmpi(rerun,'load');
end
MVI_path = VOGA_VerInfo.Path{:};
sub_info = readtable([MVI_path,filesep,'MVI_Information.xlsx']);
gno_dir = unique(extractfield(dir(strrep([MVI_path,filesep,'MVI*\Visit*\vHIT'],'\',filesep)),'folder'));
all_results = cell(length(gno_dir),1);
% Make a Summary Table 
for f = 1:length(gno_dir)
    if rerun
        disp([num2str(f),'/',num2str(length(gno_dir))])
    end
    tab = MakeGNOSummaryTable(gno_dir{f},[],rerun);
    all_results{f} = tab;
end
all_results = vertcat(all_results{:});
all_results = sortrows(sortrows(sortrows(sortrows(all_results,'Type','ascend'),...
    'Condition','ascend'),'Date','ascend'),'Subject','ascend');
[~,ind] = ismember(all_results.Subject,sub_info.Subject);
all_results.ImplantTime = years(all_results.Date-sub_info.Surgery(ind));
all_results.ActTime = years(all_results.Date-sub_info.Activation(ind));
all_results.Side(:) = {'contra'};
all_results.Side(contains(strrep(all_results.Canal,sub_info.Ear(ind),'*'),'*')) = {'ipsi'};
all_results = movevars(all_results,{'ImplantTime','ActTime','Side'}, "Before", "Experiment");
all_results = movevars(all_results,{'Gain','Gain_sd','Lat','Lat_sd','HeadVel','HeadVel_sd'}, "After", "Cycles");
save([MVI_path,filesep,'ALL_vHIT_GNO.mat'],'all_results')
val_results = all_results(:,1:find(contains(varfun(@class,all_results,'OutputFormat','cell'),'double'),1,'last'));
save([MVI_path,filesep,'vHIT_GNO.mat'],'val_results')
end