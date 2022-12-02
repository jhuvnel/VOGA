function VOGA__SummaryTable
%% Run Once
opts = {'Folder (Load)','Folder (Rerun)',...
    'All Subjects (Load)','Rerun Subject/Experiment'};
[ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
    'SelectionMode','single',...
    'ListSize',[150 125],...
    'ListString',opts);
if tf
    if contains(opts{ind},'Folder')&&~VOGA__makeFolders(cd,0) %Expecting to be in a Visit folder with the right structure
        disp('Expected folder structure not present. Navigate to appropriate directory with "Cycle Average" folder before trying again.')
        return;
    end
    if strcmp(opts{ind},'Folder (Rerun)')
        MakeCycleSummaryTable(cd,[cd,filesep,'Cycle Averages'],1);
    elseif strcmp(opts{ind},'Folder (Load)')
        MakeCycleSummaryTable(cd,[cd,filesep,'Cycle Averages'],0);
    elseif strcmp(opts{ind},'Rerun Subject/Experiment')
        CombineSummaryTables('Rerun',cd);
    elseif strcmp(opts{ind},'All Subjects (Load)')
        CombineSummaryTables('Load',cd);
    end
end
end