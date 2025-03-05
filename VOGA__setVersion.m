function VOGA__setVersion(current_ver,user_prompt)
%%20230305 CFB test
if nargin < 2
    user_prompt = 1;
end
if nargin < 1
    current_ver = '';
end
prompt = {['Set VOGA parameters.',newline,newline,'Version: '];...
'Experimenter: ';'Path: '};
if isfile([userpath,filesep,'VOGA_VerInfo.txt'])
    data = table2cell(readtable([userpath,filesep,'VOGA_VerInfo.txt'],'ReadVariableNames',false));
    items = data(:,2);
    %Set to prompt user to redo if not the right size, missing data, or outdated version
    if length(items)<3||~strcmp(current_ver,items{1})||isempty(items{3}) 
        user_prompt = 1;
    end
    def_ans = [{current_ver};items(2:3)];
else
    user_prompt = 1;
    def_ans = {current_ver;'';''};
end
if user_prompt
    items = inputdlg(prompt,'Set VOGA Parameters',1,def_ans);
end
if any(cellfun(@isempty,items))
    return;
end
info = strcat({'Version,';'Experimenter,';'Path,'},items);    
filePh = fopen([userpath,filesep,'VOGA_VerInfo.txt'],'w');
fprintf(filePh,'%s\n',info{:});
fclose(filePh);
end