%% Segment Data
%This script calls the ChinDataSineSegmentation file to batch segment all
%of the data. 

%Written by Celia Fernandez Brillet based on code from Andrianna Ayiotis

function SegmentData
%% Find names of all files that need to be segmented
%Assuming you're in the folder with the script (Sinusoids)
cd EMAOutput
files=dir(fullfile(cd,'*.csv'));
cd ../
%If a folder called Segmented does not already exist
if ~isfolder('Segmented Files')
    mkdir 'Segmented Files'
    MakePath;
end
cd 'Segmented Files'
fname = {files(:).name}';
fname(contains(fname,'STEP'))=''; % This is so that we only process SINE data, no STEP yet
%% Segment all the files at once
save_data = true;
for i = 1:length(fname)
    SineData = table2array(readtable(fname{i}));
    chin = files(1).name(10:14);
    date = files(1).name(1:8);
    rawfile = [files.folder filesep files.name];
    ChinDataSineSegmentation(date, rawfile, chin,fname{i},SineData,save_data)
end
% %% Deal with multiple versions 
% %Only for case of v2 because no v3 files were created
% %Combine into one file
% dup = dir('*v2.mat');
% dup_name = {dup(:).name}';
% for i = 1:length(dup)
%     f1 = dup_name{i};
%     v2 = load(f1);
%     f2 = [f1(1:end-7),'.mat'];
%     v1 = load(f2);
%     SegDat.info = v1.SegDat.info;
%     SegDat.info2 = v2.SegDat.info;
%     SegDat.t = v2.SegDat.t;
%     %Align the chair traces
%     Chair1 = v1.SegDat.Chair;
%     Chair2 = v2.SegDat.Chair;
%     L1 = v1.SegDat.LEye;
%     L2 = v2.SegDat.LEye;
%     R1 = v1.SegDat.REye;
%     R2 = v2.SegDat.REye;
%     C1 = mean(Chair1);
%     C2 = mean(Chair2);
%     alig = zeros(1,length(C1));
%     for j = 1:length(C1)
%         C = [C1(j:end),C1(1:j-1)];
%         alig(1,j) = sum(abs(C - C2));
%     end
%     [~,i1] = min(alig);
%     [r,c] = size(Chair1);
%     Chair1 = [[Chair1(r,i1:c);Chair1(1:r-1,i1:c)],Chair1(:,1:i1-1)];
%     L1 = [[L1(r,i1:c);L1(1:r-1,i1:c)],L1(:,1:i1-1)];
%     R1 = [[R1(r,i1:c);R1(1:r-1,i1:c)],R1(:,1:i1-1)];
%     
%     SegDat.Chair = [Chair1;Chair2];
%     SegDat.LEye = [L1;L2];
%     SegDat.REye = [R1;R2];
%     save(f2,'SegDat')
% end
% delete *v2.mat
% cd ../
end