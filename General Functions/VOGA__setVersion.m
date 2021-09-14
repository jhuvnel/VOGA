function VOGA__setVersion(new_ver)
    prompt = {['Set VOGA parameters.',newline,newline,'Version: '];...
    'Experimenter: '};
    try 
        data = readtable('VOGA_VerInfo.txt','ReadVariableNames',false);
        if nargin == 1
            items = inputdlg(prompt,'Set VOGA Parameters',1,[{new_ver},data{2,2}]);
        else
            items = inputdlg(prompt,'Set VOGA Parameters',1,data{:,2});
        end
    catch
        if nargin == 1
            items = inputdlg(prompt,'Set VOGA Parameters',1,{new_ver;''});
        else
            items = inputdlg(prompt,'Set VOGA Parameters',1);
        end
    end
    if ~isempty(items)
        info = strcat({'Version';'Experimenter'},repmat({' '},2,1),items);    
        filePh = fopen([userpath,filesep,'VOGA_VerInfo.txt'],'w');
        fprintf(filePh,'%s\n',info{:});
        fclose(filePh); 
    end
end