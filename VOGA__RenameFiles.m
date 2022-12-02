%% Rename Files
%Adds goggle type to the file name, necessary when goggle is the only
%difference between two segments (same condition/day otherwise) and they
%would therefore have the same name. This script was created to combine
%subfolders under one general 'eeVOR' or 'Rotary Chair' folder and likely
%will not be needed in the future. It looks for the string (second arg) and
%adds the goggle (first arg) in front of that with a - in between so eeVOR
%-> LDVOG2-eeVOR
%
% Some Possible Cases
%RenameFiles('eeVOR','LDVOG2-eeVOR')
%RenameFiles('eeVOR','NKI1-eeVOR')
function VOGA__RenameFiles(str1,str2,Path)
%% Handle Inputs
if nargin < 2
    items = inputdlg({'Item to Change: ','Item to Replace With: '},'');
    if isempty(items)||isempty(items{1})
        error('Process Cancelled.')
    elseif isempty(items{2})
        disp('Second field must contain at least one character.')
        return;
    end
    str1 = items{1};
    str2 = items{2};
end
if nargin < 3
    Path = cd;
end
%% Load .mat files to change the fname/name parameter, otherwise just move the file
all_files = [dir([Path,filesep,'*',filesep,'*']);dir(Path)];
all_files(extractfield(all_files,'isdir')) = [];
if any(contains({all_files.folder},'CRFs'))
    disp('CRF folder detected. Note that the CRF word document will need to updated with new .png files.')
end
for i = 1:length(all_files)  
    old_fname = [all_files(i).folder,filesep,all_files(i).name];
    fname = [all_files(i).folder,filesep,strrep(all_files(i).name,str1,str2)];
    ext = fname(find(fname=='.',1,'last'):end);
    %Move it to its new name (may be the same name if doesn't have the
    %string in it.
    if ~strcmp(old_fname,fname)
        movefile(old_fname,fname)
    end
    %Deal with changing the strings within known .mat files
    if strcmp(ext,'.mat') 
        a = load(fname);
        if isfield(a,'Data') %Segment Struct
            Data = replaceStringStruct(str1,str2,a.Data);
            save(fname,'Data')            
        elseif isfield(a,'CycAvg') %Cycle Average Struct          
            CycAvg = replaceStringStruct(str1,str2,a.CycAvg);
            save(fname,'CycAvg')            
        elseif isfield(a,'cyc_params') %Cyc Param Cell
            cyc_params = a.cyc_params;
            cyc_params(:,1) = strrep(cyc_params(:,1),str1,str2);
            save(fname,'cyc_params')  
        elseif isfield(a,'all_results') %Results Table
            tab = a.all_results;
            for j = 1:size(tab,2)
                if ischar(tab{1,j})||(iscell(tab{1,j})&&ischar(tab{1,j}{:}))
                    tab{:,j} = strrep(tab{:,j},str1,str2);
                end
            end
            all_results = tab;
            save(fname,'all_results')            
        elseif isfield(a,'tab') %Progress Table
            tab = a.tab;
            for j = 1:size(tab,2)
                if ischar(tab{1,j})||(iscell(tab{1,j})&&ischar(tab{1,j}{:}))
                    tab{:,j} = strrep(tab{:,j},str1,str2);
                end
            end
            save(fname,'tab')
        end
    elseif strcmp(ext,'.fig') 
        fig = openfig(fname);
        title = findall(fig,'Type','textbox');
        save_flag = false;
        for j = 1:length(title)
            rel_obj = title(j);
            if any(contains(rel_obj.String,str1)) %Only go through the save process if the title needs to be changed
                rel_obj.String = strrep(rel_obj.String,str1,str2);
                save_flag = true;
            end
        end
        if save_flag
            savefig(fig,fname);
        end
        close all;
    elseif strcmp(ext,'.png')||strcmp(ext,'.xlsx')
        disp('All .png and .xlsx files are being deleted and will need to be regenerated using VOGA.')
        disp(old_fname)
        delete(fname)
    elseif strcmp(ext,'.txt')&&~contains(fname,'LDHP') %Not a log file
        data = table2cell(readtable(fname,'ReadVariableNames',false,'Delimiter','*')); %The * is a delimieter we never use so that each line goes into one cell.
        if any(contains(reshape(data,[],1),str1))
            writecell(strrep(data,str1,str2),fname,'QuoteStrings',0)
        end
    end
end
end