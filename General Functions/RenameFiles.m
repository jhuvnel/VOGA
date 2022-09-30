%% Rename Files
%Adds goggle type to the file name, necessary when goggle is the only
%difference between two segments (same condition/day otherwise) and they
%would therefore have the same name. This script was created to combine
%subfolders under one general 'eeVOR' or 'Rotary Chair' folder and likely
%will not be needed in the future. It looks for the string (second arg) and
%adds the goggle (first arg) in front of that with a - in between so eeVOR
%-> LDVOG2-eeVOR
% Some Possible Cases 
%RenameFiles('eeVOR','LDVOG2-eeVOR','vis')
%RenameFiles('eeVOR','NKI1-eeVOR')
%RenameFiles('Rotary','LDVOG2-Rotary')
%RenameFiles('Rotary','NKI1-Rotary')
function RenameFiles(str1,str2,is_vis)
if strcmp(is_vis,'mat')
    fnames = extractfield([dir(['*',str1,'*.mat']);dir(['*',str1,'*.fig'])],'name');
    for i = 1:length(fnames)
        movefile(fnames{i},strrep(fnames{i},str1,str2))   
    end
elseif strcmp(is_vis,'vis')
    fnames = extractfield(dir(['*',str1,'*']),'name');
    for i = 1:length(fnames)
        if contains(fnames{i},'mat')
            a = load(fnames{i});
            b = fields(a);
            if contains(b,'Data') %Segment
                Data = a.(b{1});
                Data.info.visit = strrep(Data.info.visit,str1,str2);
                Data.info.rawfile = strrep(Data.info.rawfile,str1,str2);
                Data.info.rawnotes = strrep(Data.info.rawnotes,str1,str2);
                Data.info.name = strrep(Data.info.name,str1,str2);
                Data.rawfile = strrep(Data.rawfile,str1,str2);                
                if isfield(Data,'info2')
                    Data.info2.visit = strrep(Data.info2.visit,str1,str2);
                    Data.info2.rawfile = strrep(Data.info2.rawfile,str1,str2);
                    Data.info2.rawnotes = strrep(Data.info2.rawnotes,str1,str2);
                    Data.Data2.info.visit = strrep(Data.Data2.info.visit,str1,str2);
                    Data.Data2.info.rawfile = strrep(Data.Data2.info.rawfile,str1,str2);
                    Data.Data2.info.rawnotes = strrep(Data.Data2.info.rawnotes,str1,str2);
                    Data.Data2.rawfile = strrep(Data.Data2.rawfile,str1,str2);
                end
                save(strrep(fnames{i},str1,str2),'Data')
                delete(fnames{i})
            elseif contains(b,'CycAvg') %Cycle Average
                CycAvg = a.(b{1}); 
                CycAvg.name = strrep(CycAvg.name,str1,str2);
                if isfield(CycAvg,'parameterized')
                    CycAvg.parameterized.Visit = strrep(CycAvg.parameterized.Visit,str1,str2);
                end                
                CycAvg.info.visit = strrep(CycAvg.info.visit,str1,str2);
                CycAvg.info.rawfile = strrep(CycAvg.info.rawfile,str1,str2);
                CycAvg.info.rawnotes = strrep(CycAvg.info.rawnotes,str1,str2);
                CycAvg.info.name = strrep(CycAvg.info.name,str1,str2);
                delete(fnames{i})
                save(strrep(fnames{i},str1,str2),'CycAvg')
            else
                movefile(fnames{i},strrep(fnames{i},str1,str2))  
            end
        else
            movefile(fnames{i},strrep(fnames{i},str1,str2))   
        end
    end
elseif strcmp(is_vis,'all')
    fnames = extractfield(dir(['*',str1,'*']),'name');
    for i = 1:length(fnames)
        movefile(fnames{i},strrep(fnames{i},str1,str2))   
    end
end
end