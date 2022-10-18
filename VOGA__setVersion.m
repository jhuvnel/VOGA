function VOGA__setVersion(current_ver,user_prompt)
if nargin < 2
    user_prompt = 1;
end
if nargin < 1
    current_ver = '';
end
prompt = {['Set VOGA parameters.',newline,newline,'Version: '];...
'Experimenter: '};
if isfile([userpath,filesep,'VOGA_VerInfo.txt'])
    data = readtable([userpath,filesep,'VOGA_VerInfo.txt'],'ReadVariableNames',false);
    if ~strcmp(current_ver,data{1,2})||user_prompt 
        %Prompt the user if the version has changed it or the user has
        %flagged that they want to make a change throught the "user_prompt"
        %argument
        items = inputdlg(prompt,'Set VOGA Parameters',1,[{current_ver},data{2,2}]);
    else %User not notified if already up to date
        items = data{:,2};
    end
else %First time running
    %Use the inputed curr_ver if available 
    items = inputdlg(prompt,'Set VOGA Parameters',1,{current_ver;''}); 
end
info = strcat({'Version';'Experimenter'},repmat({' '},2,1),items);    
filePh = fopen([userpath,filesep,'VOGA_VerInfo.txt'],'w');
fprintf(filePh,'%s\n',info{:});
fclose(filePh);
end