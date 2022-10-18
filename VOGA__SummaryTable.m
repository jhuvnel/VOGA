function VOGA__SummaryTable
%% Run Once
opts = {'Folder (Existing)','Folder (New)',...
    'All Subjects (Existing)','All Subjects (New)'};
[ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
    'SelectionMode','single',...
    'ListSize',[150 125],...
    'ListString',opts);
if tf
    if contains(opts{ind},'Folder')&&~VOGA__checkFolders(0) %Expecting to be in a Visit folder with the right structure
        error('Expected folder structure not present. Navigate to appropriate directory with "Cycle Average" folder before trying again.')
    end
    if strcmp(opts{ind},'Folder (New)')
        MakeCycleSummaryTable(cd,[cd,filesep,'Cycle Averages'],1);
    elseif strcmp(opts{ind},'Folder (Existing)')
        MakeCycleSummaryTable(cd,[cd,filesep,'Cycle Averages'],0);
    elseif strcmp(opts{ind},'All Subjects (New)')
        CombineSummaryTables(cd,'All Subjects (New)');
    elseif strcmp(opts{ind},'All Subjects (Existing)')
        CombineSummaryTables(cd,'All Subjects (Existing)');
    end
end
end