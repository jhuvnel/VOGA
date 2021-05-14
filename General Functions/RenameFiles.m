%% Rename Segments
goggle = 'LDVOG2';
%goggle = 'NKI1';
str = 'eeVOR';
%str = 'Rotary';
fnames = extractfield([dir('*.mat');dir('*.fig')],'name');
for i = 1:length(fnames)
    movefile(fnames{i},strrep(fnames{i},str,[goggle,'-',str]))   
end