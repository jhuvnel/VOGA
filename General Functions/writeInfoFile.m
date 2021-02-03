function writeInfoFile(code_Path)
    prompt = {['Set VOGA parameters.',newline,newline,'Version: '];...
    'Experimenter: '};
    try 
        data = readtable([code_Path,filesep,'VerInfo.txt'],'ReadVariableNames',false);
        items = inputdlg(prompt,'Set VOGA Parameters',1,data{:,2});
    catch
        items = inputdlg(prompt,'Set VOGA Parameters',1);
    end
    if ~isempty(items)
        info = strcat({'Version';'Experimenter'},repmat({' '},2,1),items);    
        filePh = fopen([code_Path,filesep,'VerInfo.txt'],'w');
        fprintf(filePh,'%s\n',info{:});
        fclose(filePh); 
    end
end