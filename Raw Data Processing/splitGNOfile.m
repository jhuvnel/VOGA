function splitGNOfile(Path,fname,deidentify)
    %Parse relevant file name parts
    fparts = strsplit(fname,'_');
    if length(fparts)~=8
        return; %Not the right kind of file
    end
    if deidentify
        sub_name = 'DEIDENTIFIED_SUBJECT';
    else
        sub_name = [fparts{1},'_',fparts{2}];
    end
    %Open and read file
    if contains(fname,'.txt')
        fid = fopen([Path,filesep,fname],'r');
        S = textscan(fid,'%s','delimiter',newline);
        fclose(fid); 
        S = S{1};
        %Split into smaller files, renaming as needed with date, test number,
        %and test type
        header = find(contains(S,'Test'));
        ends = [header(2:end)-1;length(S)];
        for i = 1:length(header)
            info = S{header(i)};
            test_info = reshape(strsplit(info(2:end-1),{'>','<'}),3,[]);
            date_time = char(datetime(test_info{2,4},'InputFormat','MM/dd/uuuu hh:mm:ss aa'),'yyyy_MM_dd_HH_mm_SS');
            new_fname = [Path,filesep,sub_name,'_',date_time,'_',strrep(test_info{2,3},'HI - ',''),'.txt'];
            filePh = fopen(new_fname,'w');
            fprintf(filePh,'%s\n',S{header(i):ends(i),:});
            fclose(filePh);
        end
    elseif contains(fname,'.xml')    
        fid = fopen([Path,filesep,fname],'r');
        S = textscan(fid,'%s','delimiter',newline);
        fclose(fid); 
        S = S{1};
        header = find(contains(S,'<VW_HITest'));
        ends = find(contains(S,'</VW_HITest'));
        for i = 1:length(header)
            info = S(header(i):ends(i));
            date_time_str = strsplit(info{find(contains(info,'<StartDateTime'),1,'first')},{'<','>'});
            date_time = strrep(strrep(strrep(date_time_str{1,3}(1:19),'T','_'),':','_'),'-','_');
            test_type = strsplit(info{find(contains(info,'<TestType'),1,'first')},{'<','>'});
            new_fname = [Path,filesep,sub_name,'_',date_time,'_',test_type{1,3},'.xml'];
            filePh = fopen(new_fname,'w');
            fprintf(filePh,'%s\n',S{header(i):ends(i),:});
            fclose(filePh);
        end
    elseif contains(fname,'.csv')
        fid = fopen([Path,filesep,fname],'r');
        S = textscan(fid,'%s','delimiter',newline);
        fclose(fid); 
        S = S{1};
        header = find(contains(S,'Test Date'));
        ends = [header(2:end)-1;length(S)];
        for i = 1:length(header)
            info = S(header(i):ends(i));
            date_time = char(datetime(strrep(info{1},'Test Date,',''),'InputFormat','MM/dd/uuuu hh:mm:ss aa'),'yyyy_MM_dd_HH_mm_SS');
            test_type = strsplit(info{2},' ');
            new_fname = [Path,filesep,sub_name,'_',date_time,'_',test_type{1,4},'.csv'];
            filePh = fopen(new_fname,'w');
            fprintf(filePh,'%s\n',S{header(i):ends(i),:});
            fclose(filePh);
        end
    else
        disp('File type not .xml, .txt, or .csv')
        return;
    end
end