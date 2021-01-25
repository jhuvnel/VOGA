function writeInfoFile
    prompt = {['Set VOGA parameters.',newline,newline,'Version: '];...
    'Experimenter: ';'Code Path: '};
    try fid = fopen('VerInfo.txt');
        %Now get the data
        data = cell(3,2); %a little wider than needed for formatting
        tline = fgetl(fid);
        k = 0;
        while ischar(tline)
            k = k+1;
            data(k,1:length(split(tline,', '))) = split(tline,', ');
            tline = fgetl(fid);
        end
        fclose(fid);
        items = inputdlg(prompt,'Set VOGA Parameters',1,data(:,2));
    catch
        items = inputdlg(prompt,'Set VOGA Parameters',1);
    end
    info = strcat({'Version';'Experimenter';'Path'},repmat({', '},3,1),items);    
    filePh = fopen('VerInfo.txt','w');
    fprintf(filePh,'%s\n',info{:});
    fclose(filePh); 
end