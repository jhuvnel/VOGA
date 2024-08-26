%% CRF Component Creator
% Edited on 2024-04-10 to make the VOR and vHIT CRFs (v2024-03-29).
%
% Expects you to be in a Visit path that may have "Rotary Chair", "aHIT",
% "eeVOR", "Autoscan", and "vHIT" folders.
% This script requires the "Raw Files" folder to have been created for each
% of these directories already.

function VOGA__CRF(Path)
%% Defaults
default_vHIT_experimenter = 'CFB';
default_VOG_experimenter = 'EOV';
%% Initialize
if nargin < 1
    Path = cd;
end
if isempty(dir(strrep([Path,'\*\Raw Files'],'\',filesep)))
    disp('No "Raw Files" folders found. Check that you are on a path that includes VOG and vHIT files.')
    return;
end
VOGA_VerInfo = rows2vars(readtable([userpath,filesep,'VOGA_VerInfo.txt'],...
    'ReadVariableNames',false,'ReadRowNames',true));
if ~isfile([char(VOGA_VerInfo.Path),filesep,'MVI_Information.xlsx'])
    error('No MVI_Information.xlsx file found on the MVI Path in VOGA_VerInfo.txt.')
end
sub_info = readtable([char(VOGA_VerInfo.Path),filesep,'MVI_Information.xlsx']);
% Make the figure that will be used for the PDFs (single page)
fig = figure(1);
set(fig,'Units','inches','Position',[1 1 8.5 11],'Color','w');
% Figure out subject and visit
subject = strrep(char(extract(Path,'MVI'+digitsPattern(3)+'_R'+digitsPattern)),'_','');
visit = char(extract(Path,'Visit '+alphanumericsPattern(1,3)));
ind = find(ismember(sub_info.Subject,subject),1);
if ~isempty(ind)
    protocol = char(sub_info.Protocol(ind));
    if strcmp(protocol,'IRB00335294')
        fold = 'IRB00335294 NIDCD';
    elseif strcmp(protocol,'IRB00346924')
        fold = 'IRB00346924 NIA';    
    else
        fold = protocol;
    end
    visit_fold = extractfield(dir([Path(1:(strfind(Path,'Study')-1)),'CRFs',filesep,fold,filesep,subject,filesep,visit,' -*']),'name');
    if isempty(visit_fold) %Non typical visit name found
        visit_fold = {'Visit Nx - (Day XXX) Monitor - X yrs Post-Act - visit applicable only if device still act'};
    end
    CRF_path = [Path(1:(strfind(Path,'Study')-1)),'CRFs',filesep,fold,filesep,subject,filesep,visit_fold{:},filesep];
else
    CRF_path = [Path,filesep,'CRFs',filesep];
end
%% vHIT
out_path = [CRF_path,'14_05 vHIT',filesep,'14_05_CRF_vHIT_',subject,'_',visit,'.pdf'];
examiner = char(inputdlg('Name of vHIT Experimenter(s): ','Name of vHIT Experimenter(s): ',1,{default_vHIT_experimenter}));
source_path = [Path,filesep,'vHIT'];
%Find the first date
GNO_files = extractfield(dir(strrep([Path,filesep,'vHIT\GNO\Raw Files\*.txt'],'\',filesep)),'name');
if ~isempty(GNO_files)
    GNO_times = datetime(extract(GNO_files,digitsPattern(4)+'_'+digitsPattern(2)+'_'+digitsPattern+'_'+digitsPattern+'_'+digitsPattern),'Format','yyyy_MM_dd_HH_mm');
    file_times = sort(datetime(GNO_times,'Format','yyyy-MM-dd HH:mm'));
    % Make text
    CRF_txt = ['Case Report Form Protocol: ',protocol,newline,...
        'Case Report Form Version: 2024-03-29',newline,...
        'Case Report Form Test: Video Head Impulse Testing',newline,...
        'Subject ID: ',subject,newline,'Visit: ',visit,newline,'vHIT: Completed',newline,...
        'Examiners: ',examiner,newline,'Times: ',char(file_times(1)),newline,'Source Data: ',source_path];
    %Save as pdf
    clf;
    annotation('textbox','Position',[0 0 1 1],'EdgeColor','none','String',CRF_txt,'FitBoxToText','on','Interpreter', 'none')
    saveas(fig,out_path)
else
    disp('NO GNO files found. Add and rerun or manually make CRF with deviation.')
end
%% VOR
out_path = [CRF_path,'14_04 VOR',filesep,'14_04_CRF_VOR_',subject,'_',visit,'.pdf'];
examiner = char(inputdlg('Name of VOG Experimenter(s): ','Name of VOG Experimenter(s): ',1,{default_VOG_experimenter}));
exps = {'Rotary Chair','Autoscan','eeVOR','aHIT'};
exp_times = cell(1,length(exps));
for j = 1:length(exps)
    LDVOG = extractfield(dir(strrep([Path,filesep,exps{j},'\Raw Files\SESSION*.txt'],'\',filesep)),'name');
    NL = extractfield(dir(strrep([Path,filesep,exps{j},'\Raw Files\*.dat'],'\',filesep)),'name');
    if ~isempty(LDVOG)
        LDVOG = datetime(extract(LDVOG,digitsPattern(4)+lettersPattern(3)+digitsPattern(2)+'-'+digitsPattern(6)),'Format','yyyyMMMdd-HHmmSS');
    end
    if ~isempty(NL)
        NL = datetime(extract(NL,digitsPattern+'_'+digitsPattern+'_'+digitsPattern+'_'+digitsPattern+'_'+digitsPattern+'_'+lettersPattern(2)),'Format','M_d_yyyy_h_mm_a');
    end
    if ~isempty([LDVOG;NL])
        file_times = sort(datetime([LDVOG;NL],'Format','yyyy-MM-dd HH:mm'));
        exp_times{j} = char(file_times(1));
    end
end
missing_exp = cellfun(@isempty,exp_times);
if ~all(missing_exp)
    source_path = char(join(strcat([Path,filesep],exps(~missing_exp)),[',',newline]));
    datestr = char(join(exp_times(~missing_exp),', '));
    messages = {'Collected','Collected','Collected','Collected'};
    if missing_exp(1)
        messages{1} = char(inputdlg('Explain missing rotary chair data: ','Explain missing rotary chair data: ',1,{'Not collected; protocol deviation'}));
    end
    if strcmp(visit,'Visit 0')
        messages{2} = 'Not collected; test not in protocol for this visit';
        messages{3} = 'Not collected; test not in protocol for this visit';
    else
        if missing_exp(2)
            messages{2} = 'Not collected; test optional--not a protocol deviation';
        end
        if missing_exp(3)
            messages{3} = char(inputdlg('Explain missing eeVOR data: ','Explain missing eeVOR data: ',1,{'Not collected; protocol deviation'}));
        end
    end
    if missing_exp(4)
        messages{4} = 'Not collected; test optional--not a protocol deviation';
    end
    % Make text
    CRF_txt = ['Case Report Form Protocol: ',protocol,newline,...
        'Case Report Form Version: 2024-03-29',newline,...
        'Case Report Form Test: Vestibulo-ocular Reflex (VOR) Testing',newline,...
        'Subject ID: ',subject,newline,'Visit: ',visit,newline,...
        'Rotary Chair: ',messages{1},newline,'Autoscan: ',messages{2},newline,...
        'eeVOR: ',messages{3},newline,'aHIT Sinusoids: ',messages{4},newline,...
        'Examiners: ',examiner,newline,'Times: ',datestr,newline,'Source Data: ',source_path];
    %Save as pdf
    clf;
    annotation('textbox','Position',[0 0 1 1],'EdgeColor','none','String',CRF_txt,'FitBoxToText','on','Interpreter', 'none')
    saveas(fig,out_path)
else
    disp('No VOG files found. Add and rerun or manually make CRF with deviation.')
end
end