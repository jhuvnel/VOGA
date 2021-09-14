%% VOGA.m
% This script should be used to segment, filter data, select cycles, and
% parameterize batches of files. 
% To initialize after MATLAB has restarted, make sure all folders and
% subfolders within the "VOGA" folder are on the path.
% To analyze data, navivate to the path with the Raw File, Segmented Files,
% and Cycle Average folders and then run "VOGA." Click cancel to end the
% loop.
%
% Written by Andrianna Ayiotis
% Updated on 2021-05-12
opts = {'Initialize','Process Raw Data','Segment','Combine Segments','Cycle Average',...
    'Make New Summary Table','Add to Summary Table','Generate Figures','Set Version'};
ind = 1; %Run the start procedure first
tf1 = 1;
while tf1
    if strcmp(opts{ind},'Initialize')        
        code_Path = [userpath,filesep,'VOGA'];
        addpath(genpath(code_Path));
        if VOGA__checkFolders
            disp('VOGA instance ended.')
            return;
        end
    elseif strcmp(opts{ind},'Process Raw Data') 
        VOGA__ProcessRawData;    
    elseif strcmp(opts{ind},'Segment')
        VOGA__Segment;
    elseif strcmp(opts{ind},'Combine Segments')
        VOGA__combineSegments;
    elseif strcmp(opts{ind},'Cycle Average')
        VOGA__CycAvg;
    elseif strcmp(opts{ind},'Make New Summary Table')
        Path = cd;
        Cyc_Path = [Path,filesep,'Cycle Averages'];
        MakeCycleSummaryTable(Path,Cyc_Path,1);    
    elseif strcmp(opts{ind},'Add to Summary Table')
        Path = cd;
        Cyc_Path = [Path,filesep,'Cycle Averages'];
        MakeCycleSummaryTable(Path,Cyc_Path,0);
    elseif strcmp(opts{ind},'Generate Figures')  
        VOGA__makePlots;
    elseif strcmp(opts{ind},'Set Version')
        VOGA__setVersion;
    end
    %% Poll for new reponse
    [ind,tf1] = nmlistdlg('PromptString','Select an action:',...
                       'SelectionMode','single',...
                       'ListSize',[150 125],...
                       'ListString',opts);    
end
disp('VOGA instance ended.')