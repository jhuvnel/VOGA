function VOGA__SummaryTable
    %% Run Once
    opts = {'Folder (Existing)','Folder (New)',...
        'All Visits (Existing)','All Visits (New)',...
        'All Subjects (Existing)','All Subjects (New)'};    
    [ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
                       'SelectionMode','single',...
                       'ListSize',[150 125],...
                       'ListString',opts); 
    if tf
        if strcmp(opts{ind},'Folder (New)')
            MakeCycleSummaryTable(cd,[cd,filesep,'Cycle Averages'],1);      
        elseif strcmp(opts{ind},'Folder (Existing)')
            MakeCycleSummaryTable(cd,[cd,filesep,'Cycle Averages'],0);
        elseif strcmp(opts{ind},'All Visits (New)')
            CombineSummaryTables(cd,'All Visits (New)');
        elseif strcmp(opts{ind},'All Visits (Existing)')
            CombineSummaryTables(cd,'All Visits (Existing)');   
        elseif strcmp(opts{ind},'All Subjects (New)')
            CombineSummaryTables(cd,'All Subjects (New)'); 
        elseif strcmp(opts{ind},'All Subjects (Existing)')
            CombineSummaryTables(cd,'All Subjects (Existing)'); 
        end
    end
end