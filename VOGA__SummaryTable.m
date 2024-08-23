function VOGA__SummaryTable(tab_type)
opts = {'Folder (Load)','Folder (Rerun)',...
    'All Subjects (Load)','Rerun Subject/Experiment',...
    'GNO Folder (Rerun)','GNO All (Load)','GNO All (Rerun)'};
if nargin < 1 || ~ismember(tab_type,opts)
    [ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
        'SelectionMode','single','ListSize',[150 125],'ListString',opts);
    tab_type = opts{ind};
else
    tf = 1;
end
VOGA_VerInfo = rows2vars(readtable([userpath,filesep,'VOGA_VerInfo.txt'],'ReadVariableNames',false,'ReadRowNames',true));
if tf
    %Expecting to be in a Visit folder with the right structure
    if contains(tab_type,'Folder')&&~contains(tab_type,'GNO')&&~VOGA__makeFolders(cd,0) 
        disp('Expected folder structure not present. Navigate to appropriate directory with "Cycle Average" folder before trying again.')
        return;
    end
    if strcmp(tab_type,'Folder (Rerun)')
        MakeCycleSummaryTable(cd,[cd,filesep,'Cycle Averages'],'Rerun');
    elseif strcmp(tab_type,'Folder (Load)')
        MakeCycleSummaryTable(cd,[cd,filesep,'Cycle Averages'],'Load');
    elseif strcmp(tab_type,'Rerun Subject/Experiment')
        CombineSummaryTables(VOGA_VerInfo,'Rerun');
    elseif strcmp(tab_type,'All Subjects (Load)')
        CombineSummaryTables(VOGA_VerInfo,'Load');
    elseif strcmp(tab_type,'GNO Folder (Rerun)')
        MakeGNOSummaryTable(cd,[cd,filesep,'GNO',filesep,'Raw Files'],'Rerun');
    elseif strcmp(tab_type,'GNO All (Load)')
        CombineGNOTables(VOGA_VerInfo,'Load');
    elseif strcmp(tab_type,'GNO All (Rerun)')
        CombineGNOTables(VOGA_VerInfo,'Rerun');
    end
end
end