VOGA_VerInfo = rows2vars(readtable([userpath,filesep,'VOGA_VerInfo.txt'],...
    'ReadVariableNames',false,'ReadRowNames',true));
MVI_path = VOGA_VerInfo.Path{:};
params.version = VOGA_VerInfo.Version{:};
params.Experimenter = VOGA_VerInfo.Experimenter{:};
sub_info = readtable([MVI_path,filesep,'MVI_Information.xlsx']);
params.sub_info = sub_info;
%VOG_fnames = dir([MVI_path,filesep,'MVI*',filesep,'Visit*',filesep,'eeVOR',filesep,'*Results.mat']);
VOG_fnames = unique(strrep(strrep(extractfield(dir([MVI_path,filesep,...
    'MVI*',filesep,'Visit*',filesep,'eeVOR',filesep,'*',filesep,'*Velstep*']),'folder'),...
    [filesep,'Cycle Averages'],''),[filesep,'Segmented Files'],''));
%%
for j = 1:length(VOG_fnames)
    %Path = VOG_fnames(j).folder;
    Path = VOG_fnames{j};
    disp([num2str(j),'/',num2str(length(VOG_fnames)),': ',Path])
    VOGA__makeFolders(Path,1,0);
    try
        MakeCycleSummaryTable(Path,[Path,filesep,'Cycle Averages'],1);
        params.Path = Path;
        params.Cyc_Path = [Path,filesep,'Cycle Averages'];
        plotParamResults(params);
    catch
        disp(['FIX THIS FOLDER: ',Path])
    end
end