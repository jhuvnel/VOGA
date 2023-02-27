%% VOGA.m
% Video Oculography Analyzer (VOGA)
%
% This script should be used to connect raw files with metadata, segment,
% filter, select cycles, parameterize, and plot batches of files
% containing VOG data. 
% To analyze data, navivate to the path with the Raw File, Segmented Files,
% and Cycle Average folders and then run "VOGA." Click cancel to end the
% loop.
%Current version of VOGA - Changed with each GitHub committ
current_ver = 'v5.0.0';
%VOGA Menu Options
opts = {'Generate Folders','Process Raw Data','Segment','Cycle Average',...
    'Summary Table','Generate Figures','CRF','Advanced'};
advanced_opts = {'Automatic VOG Analysis','Recalibrate LDVOG','Combine Segments',...
    'Trim Segment','Rename Files','Set Version'};
resp1 = '';
tf1 = 1;
%Run the start procedure first--will make the files needed if they don't
%exist yet
VOGA__setVersion(current_ver,0);
VOGA__saveLastUsedParams;
while tf1
    switch resp1
        case 'Generate Folders'
            VOGA__makeFolders(cd,1,0);
        case 'Process Raw Data'
            VOGA__ProcessRawData;
        case 'Segment'
            VOGA__Segment;
        case 'Cycle Average'
            VOGA__CycAvg;
        case 'Summary Table'
            VOGA__SummaryTable;
        case 'Generate Figures'
            VOGA__makePlots;
        case 'CRF'
            VOGA__CRF;
        case 'Advanced'
            [ind2,tf2] = nmlistdlg('PromptString','Select an action:','SelectionMode','single',...
                       'ListSize',[150 125],'ListString',advanced_opts); 
            if tf2
                resp2 = advanced_opts{ind2};
                switch resp2
                    case 'Automatic VOG Analysis'
                        VOGA__automaticAnalysis;
                    case 'Recalibrate LDVOG'
                        RecalibrateRawLDVOG; 
                    case 'Combine Segments'
                        VOGA__combineSegments;
                    case 'Trim Segment'
                        VOGA__trimSegment;
                    case 'Rename Files'
                        VOGA__RenameFiles;
                    case 'Set Version'
                        VOGA__setVersion(current_ver,1);
                end
            end           
    end
    % Poll for new reponse
    [ind1,tf1] = nmlistdlg('PromptString','Select an action:','SelectionMode','single',...
                       'ListSize',[150 125],'ListString',opts); 
    if tf1
        resp1 = opts{ind1}; 
    end
end
disp('VOGA instance ended.')