%% VOGA.m
% This script should be used to segment, filter data, select cycles, and
% parameterize batches of files. 
% To initialize after MATLAB has restarted, make sure all folders and
% subfolders within the "VOGA" folder are on the path by running "addVOGA."
% To analyze data, navivate to the path with the Raw File, Segmented Files,
% and Cycle Average folders and then run "VOGA." Click cancel to end the
% loop.
%
% Written by Andrianna Ayiotis
% Updated on 2022-02-14
opts = {'Process Raw Data','Segment','Cycle Average','Summary Table',...
    'Generate Figures','Advanced'};
advanced_opts = {'Recalibrate LDVOG','Combine Segments','Trim Segment','Set Version'};
%Run the start procedure first
code_Path = addVOGA;
if VOGA__checkFolders
    disp('VOGA instance ended.')
    return;
end
resp1 = '';
tf1 = 1;
while tf1
    switch resp1
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
        case 'Advanced'
            [ind2,tf2] = nmlistdlg('PromptString','Select an action:','SelectionMode','single',...
                       'ListSize',[150 125],'ListString',advanced_opts); 
            if tf2
                resp2 = advanced_opts{ind2};
                switch resp2
                    case 'Recalibrate LDVOG'
                        RecalibrateRawLDVOG; 
                    case 'Combine Segments'
                        VOGA__combineSegments;
                    case 'Trim Segment'
                        VOGA__trimSegment;
                    case 'Set Version'
                        VOGA__setVersion;
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