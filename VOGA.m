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
% Updated on 2021-03-27
opts = {'Initialize','Make Folders','Segment','Cycle Average',...
    'Summary Table','Generate Figures','Set Version'};
ind = 1; %Run the start procedure first
tf1 = 1;
while tf1
    if strcmp(opts{ind},'Initialize')        
        code_Path = [userpath,filesep,'VOGA'];
        addpath(genpath(code_Path));
    elseif strcmp(opts{ind},'Make Folders') 
        VOGA__makeFolders;
    elseif strcmp(opts{ind},'Segment') 
        flag = VOGA__checkFolders;
        if ~flag
            VOGA__Segment;
        end
    elseif strcmp(opts{ind},'Cycle Average')
        flag = VOGA__checkFolders;
        if ~flag
            VOGA__CycAvg(code_Path);
        end
    elseif strcmp(opts{ind},'Summary Table')
        flag = VOGA__checkFolders;
        if ~flag
            Path = cd;
            Cyc_Path = [Path,filesep,'Cycle Averages'];
            MakeCycleSummaryTable(Path,Cyc_Path);
        end
    elseif strcmp(opts{ind},'Generate Figures')
        flag = VOGA__checkFolders;
        if ~flag  
            Path = cd;
            Seg_Path = [Path,filesep,'Segmented Files'];
            Cyc_Path = [Path,filesep,'Cycle Averages'];
            % Get version and experimenter info from the file
            if ~any(contains(extractfield(dir(code_Path),'name'),'VerInfo.txt'))
                writeInfoFile(code_Path);
            end
            data = readtable([code_Path,filesep,'VerInfo.txt'],'ReadVariableNames',false);
            version = data{1,2}{:};
            Experimenter = data{2,2}{:};
            annot_q = questdlg('Descriptive annotation on figures?','','Yes','No','Yes');
            annot = strcmp(annot_q,'Yes');
        end
        return;
    elseif strcmp(opts{ind},'Set Version')
        VOGA__setVersion(code_Path);
    end
    %% Poll for new reponse
    [ind,tf1] = nmlistdlg('PromptString','Select an action:',...
                       'SelectionMode','single',...
                       'ListSize',[150 125],...
                       'ListString',opts,...
                       'Position',[0,7.75,3,3]);    
end
disp('VOGA instance ended.')