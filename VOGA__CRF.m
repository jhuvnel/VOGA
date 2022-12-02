%% CRF Component Creator
% CRFs for eeVOR, Rotary Chair and vHIT need Cycle Average Figures and
% tables with the components of each test

%This script will make:
%(1) .txt file with the information needed underneath the header
%(2) .png files with the relevant figures (need to be copied into the word
%doc CRF)
%(3) .xlsx file with the table of results

%This script requres a Results.mat table to have already been made and
%Cycle Average Figures to already have been generated to fully function. It
%will create some of the files depending on availability of data.

function VOGA__CRF(Path)
%% Initialize
if nargin < 1
    Path = cd;
end
if ~VOGA__makeFolders(Path)
    return;
end
Raw_Path = [Path,filesep,'Raw Files'];
Cyc_Path = [Path,filesep,'Cycle Averages'];
Fig_Path = [Path,filesep,'Figures'];
CRF_Path = [Path,filesep,'CRFs'];
rel_tab = extractfield(dir([Path,filesep,'*Results.mat']),'name');
if isempty(rel_tab)
    error('No Results.mat table has been generated yet. Complete analysis first or create the CRF manually.')
end
load([Path,filesep,rel_tab{end}],'all_results')
all_rawfiles = cell(size(all_results,1),1);
for i = 1:length(all_rawfiles)
    load([Cyc_Path,filesep,all_results.File{i}],'CycAvg')
    all_rawfiles{i} = strrep(CycAvg.info.rawfile,[Raw_Path,filesep],'');
end
rel_figs = extractfield(dir([Fig_Path,filesep,'CycleAverage*.fig']),'name');
%% Text file
sub = strjoin(unique(all_results.Subject),', ');
date = strjoin(unique(cellstr(datestr(all_results.Date,'yyyy-mm-dd'))),', ');
vis = strjoin(unique(all_results.Visit),', ');
test_name = strjoin(unique(strcat(all_results.Experiment,'-',all_results.Type,'-',all_results.Condition,'-',all_results.AxisName)),newline);
gog = strjoin(unique(all_results.Goggle),', ');
raw_files = strjoin(unique(all_rawfiles),newline);
examiner = inputdlg('Name of Experimenter(s): ','Name of Experimenter(s): ',1,{'AIA'});
examiner = examiner{:};
txt = ['Subject ID: ',sub,newline,...
    'Date: ',date,newline,...
    'Visit ID: ',vis,newline,...
    'Test Name: ',test_name,newline,...
    'Goggle ID: ',gog,newline,...
    'Raw Filename(s): ',raw_files,newline,...
    'Examiner(s): ',examiner];
fid = fopen([CRF_Path,filesep,'SubjectInfo.txt'],'w');
fprintf(fid,'%s',txt);
fclose(fid);
%% Png Files
if ~isempty(rel_figs)
    close all;
    for i = 1:length(rel_figs)
        fig1 = open([Fig_Path,filesep,rel_figs{i}]);
        saveas(fig1,[CRF_Path,filesep,strrep(rel_figs{i},'.fig','.png')])
        close;
    end 
else
    disp('No Cycle Average figures found.')
end
%% Xlsx Tables
test_type = cell2table(strrep(strrep(all_results.File,'CycAvg_',''),'.mat',''),'VariableNames',{'TestConditon'});
%Use the indeces for the maximum eye velocity, may need to change #s for
%future versions
LHRH_inds = [18,19,17,22,23,17,20,21,17,24,25,17];
LARP_inds = [26,27,17,30,31,17,28,29,17,32,33,17];
RALP_inds = [34,35,17,38,39,17,36,37,17,40,41,17];
CRF_tabs.LHRH = [test_type,all_results(:,LHRH_inds)];
CRF_tabs.LHRH{:,2:end} = round(CRF_tabs.LHRH{:,2:end},1); %round to 0.1 decimal place
CRF_tabs.LARP = [test_type,all_results(:,LARP_inds)];
CRF_tabs.LARP{:,2:end} = round(CRF_tabs.LARP{:,2:end},1); %round to 0.1 decimal place
CRF_tabs.RALP = [test_type,all_results(:,RALP_inds)];
CRF_tabs.RALP{:,2:end} = round(CRF_tabs.RALP{:,2:end},1); %round to 0.1 decimal place
CRF_tabs.info = [all_results(:,[2,4,3,6,7]),cell2table(all_rawfiles,'VariableNames',{'RawFileName'}),all_results(:,1)];
CRF_tabs.info.Examiner = repmat({examiner},size(all_results,1),1);
save([CRF_Path,filesep,'CRF_Tables.mat'],'CRF_tabs') %save as .mat for easy future loading
writetable(CRF_tabs.LHRH,[CRF_Path,filesep,'CRF_Tables.xlsx'],'Sheet','LHRH');
writetable(CRF_tabs.LARP,[CRF_Path,filesep,'CRF_Tables.xlsx'],'Sheet','LARP');
writetable(CRF_tabs.RALP,[CRF_Path,filesep,'CRF_Tables.xlsx'],'Sheet','RALP');
writetable(CRF_tabs.info,[CRF_Path,filesep,'CRF_Tables.xlsx'],'Sheet','Info');
%% End
disp('Items now in the CRFs folder')
end