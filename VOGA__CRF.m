%% CRF Component Creator
% CRFs for eeVOR, Rotary Chair and vHIT need Cycle Average Figures and
% tables with the components of each test

%This script will make:
%(1) .txt file with the information needed underneath the header
%(2) .xlsx file with the table of results

%This script requres a Results.mat table to have already been made and
%Cycle Average Figures to already have been generated to fully function. It
%will create some of the files depending on availability of data.

function VOGA__CRF(Path)
%% Initialize
if nargin < 1
    Path = cd;
end
rel_tab = extractfield(dir([Path,filesep,'*Results.mat']),'name');
if ~VOGA__makeFolders(Path,0,0)
    disp('Expected file structure missing.')
    return;
elseif isempty(rel_tab)
    error('No Results.mat table has been generated yet. Complete analysis first or create the CRF manually.')
end
CRF_Path = [Path,filesep,'CRFs'];
load([Path,filesep,rel_tab{end}],'all_results')
all_results(cellfun(@(x) all(x==0),all_results.StimAxis),:) = []; %Remove any "rest" experiments
all_results.RawFile(:) = {''};
for i = 1:size(all_results,1)
    load([Path,filesep,'Cycle Averages',filesep,all_results.File{i}],'CycAvg')
    rawfile = split(strrep(strrep(CycAvg.info.rawfile,'\',filesep),'/',filesep),filesep);
    all_results.RawFile(i) = rawfile(end);
end
%% Text file
StimAxStr = strcat('[',strrep(strrep(cellfun(@(x) num2str(x),...
    all_results.StimAxis,'UniformOutput',false),'  ',' '),' ',','),']');
StimAxStr(~cellfun(@isempty,all_results.AxisName)) = all_results.AxisName(~cellfun(@isempty,all_results.AxisName));
date_str = cellstr(datestr(all_results.Date,'yyyy-mm-dd'));
sub = strjoin(unique(all_results.Subject),', ');
date = strjoin(unique(date_str),', ');
vis = strjoin(unique(all_results.Visit),', ');
test_name = strjoin(unique(strcat(all_results.Experiment,'-',...
    all_results.Type,'-',all_results.Condition,'-',StimAxStr),'stable'),newline);
gog = strjoin(unique(all_results.Goggle),', ');
raw_files = strjoin(strrep(extractfield(dir([Path,filesep,'Raw Files',filesep,'*-Notes.txt']),'name'),'-Notes',''),newline);
%examiner = inputdlg('Name of Experimenter(s): ','Name of Experimenter(s): ',1,{'AIA,EOV,RS'});
examiner = 'AIA';
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
%% Xlsx Tables
test_cond = strcat(all_results.Subject,'-',all_results.Visit,'-',date_str,'-',...
    all_results.Experiment,'-',all_results.Type,'-',all_results.Condition,'-',StimAxStr);
vel_labs = {'Test Condition','LE: Vel Max [deg/s]','LE: Vel STD [deg/s]',...
    'LE: # Cycles','RE: Vel Max [deg/s]','RE: Vel STD [deg/s]','RE: # Cycles'};
nc = length(vel_labs);
info_labs = {'Test Condition','Raw File Name','Processed File Name','Examiner(s)'};
info_cell = cell(length(test_cond),nc);
info_cell(:,1:4) = [test_cond,all_results.RawFile,all_results.File,...
    repmat({examiner},length(test_cond),1)];
%Round these all the the nearest 0.1 dps
tab_labs = {'MaxVel_L&','MaxVel_L&_sd','Cycles','MaxVel_R&','MaxVel_R&_sd','Cycles'};
CRF_cell = [{'LHRH Component'},cell(1,nc-1);vel_labs;...
    test_cond,num2cell(round(all_results{:,strrep(tab_labs,'&','Z')},1));...
    {'LARP Component'},cell(1,nc-1);vel_labs;...
    test_cond,num2cell(round(all_results{:,strrep(tab_labs,'&','L')},1));...
    {'RALP Component'},cell(1,nc-1);vel_labs;...
	test_cond,num2cell(round(all_results{:,strrep(tab_labs,'&','R')},1));...
    {'Other Info'},cell(1,nc-1);info_labs,cell(1,nc-length(info_labs));...
    info_cell];
%Make the output cell
writecell(CRF_cell,[CRF_Path,filesep,'CRF_Table.xlsx']);
%% End
disp('Items now in the CRFs folder')
end