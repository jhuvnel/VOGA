function VOGA__SummaryTable
%% Run Once
opts = {'Folder (Load)','Folder (Rerun)',...
    'All Subjects (Load)','Rerun Subject/Experiment',...
    'GNO Folder (Rerun)','GNO All (Load)','GNO All (Rerun)'};
[ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
    'SelectionMode','single',...
    'ListSize',[150 125],...
    'ListString',opts);
VOGA_VerInfo = rows2vars(readtable([userpath,filesep,'VOGA_VerInfo.txt'],'ReadVariableNames',false,'ReadRowNames',true));
if tf
    %Expecting to be in a Visit folder with the right structure
    if contains(opts{ind},'Folder')&&~contains(opts{ind},'GNO')&&~VOGA__makeFolders(cd,0) 
        disp('Expected folder structure not present. Navigate to appropriate directory with "Cycle Average" folder before trying again.')
        return;
    end
    if strcmp(opts{ind},'Folder (Rerun)')
        MakeCycleSummaryTable(cd,[cd,filesep,'Cycle Averages'],'Rerun');
    elseif strcmp(opts{ind},'Folder (Load)')
        MakeCycleSummaryTable(cd,[cd,filesep,'Cycle Averages'],'Load');
    elseif strcmp(opts{ind},'Rerun Subject/Experiment')
        CombineSummaryTables(VOGA_VerInfo,'Rerun');
    elseif strcmp(opts{ind},'All Subjects (Load)')
        CombineSummaryTables(VOGA_VerInfo,'Load');
    elseif strcmp(opts{ind},'GNO Folder (Rerun)')
        MakeGNOSummaryTable(cd,[cd,filesep,'GNO',filesep,'Raw Files'],'Rerun');
    elseif strcmp(opts{ind},'GNO All (Load)')
        CombineGNOTables(VOGA_VerInfo,'Load');
    elseif strcmp(opts{ind},'GNO All (Rerun)')
        CombineGNOTables(VOGA_VerInfo,'Rerun');
    end
end
end