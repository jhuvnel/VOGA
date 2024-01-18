%% VOGA.m
% Video Oculography Analyzer (VOGA)
%
% This script should be used to connect raw files with metadata, segment,
% filter, select cycles, parameterize, and plot batches of files
% containing VOG data. 
% To analyze data, navivate to the path with the Raw File, Segmented Files,
% and Cycle Average folders and then run "VOGA." Click cancel to end the
% loop.
%Current version of VOGA - Changed with each GitHub commit
current_ver = 'v5.6.3';
%VOGA Menu Options
opts = {'Generate Folders','Process Raw Data','Segment','Cycle Average',...
    'CRF','Generate Figures','Advanced'};
advanced_opts = {'Automatic VOG Analysis','Recalibrate LDVOG',...
    'Fix Raw VOG Trigger','Combine Segments','Trim Segment','Rename Files',...
    'Summary Table','Set Version'};
%Run the start procedure first--will make the files needed if they don't
%exist yet
VOGA__setVersion(current_ver,0);
VOGA__saveLastUsedParams;
resp = '';
tf = 1;
while tf
    if strcmp(resp,'Advanced') %Give the advanced menu and run that selection
        [ind2,tf2] = nmlistdlg('SelectionMode','single','ListString',advanced_opts); 
        if tf2
            resp = advanced_opts{ind2};
        end
    end    
    switch resp
        case 'Generate Folders'
            VOGA__makeFolders(cd,1,0);
        case 'Process Raw Data'
            VOGA__ProcessRawData;
        case 'Segment'
            VOGA__Segment;
        case 'Cycle Average'
            VOGA__CycAvg;
        case 'CRF'
            VOGA__CRF;
        case 'Generate Figures'
            VOGA__makePlots;        
        case 'Automatic VOG Analysis'
            VOGA__automaticAnalysis;
        case 'Fix Raw VOG Trigger'
            updateRawVOGTrigger;
        case 'Recalibrate LDVOG'
            RecalibrateRawLDVOG; 
        case 'Combine Segments'
            VOGA__combineSegments;
        case 'Trim Segment'
            VOGA__trimSegment;
        case 'Rename Files'
            VOGA__RenameFiles;
        case 'Summary Table'
            VOGA__SummaryTable;
        case 'Set Version'
            VOGA__setVersion(current_ver,1);                         
    end
    if ismember(resp,advanced_opts) %Stay in advanced menu if started in advanced menu
        resp = 'Advanced';
    else % Poll for new reponse from main menu
        [ind1,tf] = nmlistdlg('SelectionMode','single','ListString',opts);
        if tf
            resp = opts{ind1};
        end
    end    
end
disp('VOGA instance ended.')