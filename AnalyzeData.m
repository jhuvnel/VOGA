%% AnalyzeData.m
% This script should be used to segment, filter data, select cycles, and
% parameterize batches of files. 
% Updated on 2020-11-12
%% Initialize File Tree and Code Path
code_Path = '/Volumes/MVI/DATA SUMMARY/IN PROGRESS/VOG Analysis Scripts/LDVOG_Neurolign';
addpath(genpath(code_Path))
% Assumes you are in the right directory already
path = cd;
%path = uigetdir([],'Select the visit folder where the folder tree should be.');
%cd(path)
%Set paths
Raw_Path = [path,filesep,'Raw Files'];
Seg_Path = [path,filesep,'Segmented Files'];
Cyc_Path = [path,filesep,'Cycle Averages'];
%See if the folders already exist or need to be renamed/created
path_folders = extractfield(dir,'name',find(extractfield(dir,'isdir')));
if any(contains(path_folders,'Raw LD VOG Files'))
    movefile('Raw LD VOG Files','Raw Files')
elseif ~any(contains(path_folders,'Raw Files'))
    mkdir('Raw Files')
end
if any(contains(path_folders,'CycAvg')) 
    movefile('CycAvg','Cycle Averages')
elseif any(contains(path_folders,'Cyc_Avg')) 
    movefile('Cyc_Avg','Cycle Averages')
elseif ~any(contains(path_folders,'Cycle Averages'))
    mkdir('Cycle Averages')
end
if ~any(contains(path_folders,'Segmented Files'))
    mkdir('Segmented Files')
end
%% Prepare to Segment
%Transfer NKI Raw Files from their subfolders if they exist
moveNKIfiles(Raw_Path)
%Detect and process log files and austoscan files
logtoNotes(Raw_Path)
%% Segment
% Select files to segment
all_files = extractfield(dir(Raw_Path),'name',find(~extractfield(dir(Raw_Path),'isdir')));
VOG_files = all_files(((contains(all_files,'SESSION')&contains(all_files,'.txt'))|contains(all_files,'.dat'))...
    &~contains(all_files,'-Notes.txt'));
if isempty(VOG_files) 
    uiwait(msgbox('No VOG files found in the Raw Files folder.'))
else
    [indx,tf] = nmlistdlg('PromptString','Select files to segment:','ListSize',[300 300],'ListString',VOG_files);
    if tf == 1
        sel_files = VOG_files(indx);
        for i = 1:length(sel_files)
            In_Path = [Raw_Path,filesep,sel_files{i}];
            Segment(In_Path,Seg_Path)
        end
    end
end
%% Filter and select cycles
done = false;
Experimenter = 'AIAyiotis';
version = 'VOGA_v1';
while(~done) %run until the user hits cancel on analyzing a file
    done = MakeCycAvg(path,Seg_Path,Cyc_Path,Experimenter,version);
end