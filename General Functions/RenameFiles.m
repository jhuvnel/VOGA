%% Rename Files
%Adds goggle type to the file name, necessary when goggle is the only
%difference between two segments (same condition/day otherwise) and they
%would therefore have the same name. This script was created to combine
%subfolders under one general 'eeVOR' or 'Rotary Chair' folder and likely
%will not be needed in the future. It looks for the string (second arg) and
%adds the goggle (first arg) in front of that with a - in between so eeVOR
%-> LDVOG2-eeVOR
% Some Possible Cases 
%RenameFiles('LDVOG2','eeVOR')
%RenameFiles('NKI1','eeVOR')
%RenameFiles('LDVOG2','Rotary')
%RenameFiles('NKI1','Rotary')
function RenameFiles(goggle,str)
fnames = extractfield([dir('*.mat');dir('*.fig')],'name');
for i = 1:length(fnames)
    movefile(fnames{i},strrep(fnames{i},str,[goggle,'-',str]))   
end
end