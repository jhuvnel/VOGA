%Run from MVI/Study Subjects
GNO_files = dir([cd,'\MVI*R*\Visit*\vHIT\GNO\Raw Files\*.xml']);
%GNO_files = dir([cd,'\Visit*\vHIT\GNO\Raw Files\*.xml']);
relevant_num = cell(length(GNO_files),1);
file_date = NaT(length(GNO_files),1);
rel_str = 'GogglesSN';
for i = 1:length(GNO_files)
    fname = GNO_files(i).name;
    fparts = split(fname,'_');
    file_date(i) = datetime(str2double(fparts(3:8)'));
    fdata = cellstr(readlines([GNO_files(i).folder,filesep,GNO_files(i).name]));
    if any(contains(fdata,rel_str))
        gog_line = fdata(contains(fdata,rel_str));
        r_num = cell(length(gog_line),1);
        for j = 1:length(gog_line)
            rel_line = gog_line{j};
            try 
                r_num{j} = extractXMLdataline(rel_line);
            catch
                r_num{j} = '';
            end
        end    
        if length(r_num)>1 %Combine same gog #
            r_num = unique(r_num);
        end
        if length(r_num)>1 %If more than one google type, join with a slash
            r_num = {strjoin(r_num,'/')};
        end
        if isempty(r_num{:})
            relevant_num(i) = {'Undetected'};
        else
            relevant_num(i) = r_num;
        end        
    else
        relevant_num(i) = {'Undetected'};
    end
end
full_info = table();
full_info.Name = extractfield(GNO_files,'name');
full_info.Folder = extractfield(GNO_files,'folder');
full_info.Date = file_date;
full_info.Data = relevant_num;
full_info = sortrows(full_info,'Date','ascend');
disp(unique(full_info.Data,'stable'))
%% Folder info
all_vis = unique(full_info.Folder,'stable');
rel_dat = cell(length(all_vis),1);
for i = 1:length(all_vis)
    temp = unique(full_info.Data(strcmp(full_info.Folder,all_vis{i})));
    if length(temp)>1
        temp = {strjoin(temp,'/')};
    end
    rel_dat(i) = temp;
end
vis_info = table();
vis_info.Folder = all_vis;
vis_info.Data = rel_dat;