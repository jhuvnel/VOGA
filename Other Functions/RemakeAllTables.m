VOGA_VerInfo = rows2vars(readtable([userpath,filesep,'VOGA_VerInfo.txt'],'ReadVariableNames',false,'ReadRowNames',true));
MVI_path = VOGA_VerInfo.Path{:};
params.version = VOGA_VerInfo.Version{:};
params.Experimenter = VOGA_VerInfo.Experimenter{:};
sub_info = readtable([MVI_path,filesep,'MVI_Information.xlsx']);
params.sub_info = sub_info;
VOG_fnames = dir([MVI_path,filesep,'MVI*',filesep,'Visit*',filesep,'vHIT',filesep,'GNO',filesep,'*Results.mat']);
VOG_fnames(contains(extractfield(VOG_fnames,'name'),'20230224'),:) = [];
%%
for j = 1:length(VOG_fnames)
    disp([num2str(j),'/',num2str(length(VOG_fnames)),': ',VOG_fnames(j).folder])
    VOGA__makeFolders(VOG_fnames(j).folder,1,0);
    try
        MakeCycleSummaryTable(VOG_fnames(j).folder,[VOG_fnames(j).folder,filesep,'Cycle Averages'],1);
        params.Path = VOG_fnames(j).folder;
        params.Cyc_Path = [VOG_fnames(j).folder,filesep,'Cycle Averages'];
        plotParamResults(params);
    catch
        disp(['FIX THIS FOLDER: ',VOG_fnames(j).folder])
    end
end