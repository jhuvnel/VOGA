%% VOGA (Video Oculography Analyzer)
%
% The purpose of VOGA is to allow a user to be able to analyze
% video-oculography (VOG) data. This includes associating files with their
% metadata, splitting different experiments into seperate files
% (segmenting), filtering the VOG in position and/or velocity to remove
% high frequency noise, blinks, and quick phases, amalgamating cycles of
% the same experiment together, extracting relevant parameters, making
% summary plots, and generating Case Report Forms (CRFs).
%
% This is the main menu script that allows the user to navigate to the
% desired function. To run VOG, navivate in MATLAB to the path with the 
% VOG data and type "VOGA" into the Command Window. Select the desired 
% functionality and click cancel to end the program. Start with "Generate
% Folders" if the folders of Raw Files, Segmented Files, Cycle
% Averages, and/or Figures don't yet exist. 
%
% This code is available publically on Github at
% https://github.com/jhuvnel/VOGA. It will save a VOGA_VerInfo.txt file and
% a VOGA_LastUsedFilterParams.mat file in MATLAB's default userpath. For
% optimal use, add the folder this code is housed in (and all subfolders)
% to the default path on the computer it is being run on.
%
% Originally created by Andrianna Ayiotis in January of 2021.
%
% % Version Control 
% The standard method for version control numbering is as follows:
% 1. The first number is incremented if the changes are not back-compatible
% with each other. 
% 2. The second number is incremented if a feature is added or removed.
% 3. The third number is incremented for a bug fix.
% 4. When the first or second numbers increment, the numbers after them go
% back to 0.
% Make sure to add the current version into the commit for ease of tracking in the
% future.
current_ver = 'v5.20.0'; %CHANGE ME before each new commit
% % VOGA Menu Options
% It unwieldy to have all of the possible functionalities in the main menu
% so they are split by the options needed in typical folder analysis and 
% the "Advanced" options. More items can be added as needed.
opts = {'Generate Folders','Process Raw Data','Segment','Cycle Average',...
    'Parameterize and Plot','Advanced'};
advanced_opts = {'Automatic Analysis','Set Version','Rename Files',...
    'Recalibrate LDVOG','Update VOG Trigger','Manually Segment',...
    'Combine Segments','Trim Segment','CRF','Summary Table','Make Plots'};
resp = ''; tf = 1; %Starting values so the loop runs
while tf
    if strcmp(resp,'Advanced') %Give the advanced menu and run that selection
        [ind2,tf2] = nmlistdlg('SelectionMode','single','ListString',advanced_opts); 
        resp = '';
        if tf2
            resp = advanced_opts{ind2};
        end
    end    
    switch resp
        case '' %Run the start procedure first--will make the files needed if they don't exist yet
            SetVersion(current_ver,0); %Updates the version # in the VOGA_VerInfo.txt file if different
            SaveLastUsedParams; %Makes sure there is a VOGA_LastUsedFilterParams.mat in the userpath directory
        case 'Set Version' % In case the user wants to make edits to the VOGA_VerInfo.txt file
            SetVersion(current_ver,1); 
        case 'Generate Folders'
            MakeFolders(cd,1,0);
        case 'Manually Segment'
            VOGA__Segment(cd,1);
        case 'Parameterize and Plot'
            VOGA__SummaryTable('Folder (Load)') %Parameterize
            VOGA__MakePlots('Parameterized') %Plot        
        case {'Update VOG Trigger','Recalibrate LDVOG','Rename Files','Trim Segment'}
            run(strrep(resp,' ','')) %Scripts where the menu option is the name of the script
        otherwise
            run(['VOGA__',strrep(resp,' ','')]) %Scripts where VOGA__ and then the menu option is the name of the script
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