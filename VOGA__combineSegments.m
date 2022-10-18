function VOGA__combineSegments
Path = cd;
Seg_Path = [cd,filesep,'Segmented Files'];
% File Names
rm_header = @(str) str([strfind(str,'eeVOR'),strfind(str,'aHIT'),strfind(str,'Rotary'),strfind(str,'Manual')]:end);
progress_tab = assessProgress(Path);
all_exp_names = progress_tab.Segment;
[~,sorti] = sort(cellfun(rm_header,all_exp_names,'UniformOutput',false));
all_exp_names = all_exp_names(sorti);
[indx1,tf1] = nmlistdlg('PromptString','Select segments to combine. Segment 1 (combined segment name):',...
    'SelectionMode','single','ListSize',[500 500],'ListString',all_exp_names);
if tf1
    [indx2,tf2] = nmlistdlg('PromptString','Select segments to combine. Segment 2:',...
        'SelectionMode','single','ListSize',[500 500],...
        'ListString',all_exp_names);
end
while(tf1&&tf2)
    if indx1==indx2
        uiwait(msgbox('Cannot combine a segment with itself.'))
    else
        load([Seg_Path,filesep,all_exp_names{indx1}],'Data');
        Data1 = Data;
        load([Seg_Path,filesep,all_exp_names{indx2}],'Data');
        Data2 = Data;
        Data = CombineSegments(Data1,Data2);
        delete([Seg_Path,filesep,all_exp_names{indx1}])
        delete([Seg_Path,filesep,all_exp_names{indx2}])
        save([Seg_Path,filesep,all_exp_names{indx1}],'Data')
    end
    % File Names
    progress_tab = assessProgress(Path);
    all_exp_names = progress_tab.Segment;
    [~,sorti] = sort(cellfun(rm_header,all_exp_names,'UniformOutput',false));
    all_exp_names = all_exp_names(sorti);
    [indx1,tf1] = nmlistdlg('PromptString','Select segments to combine. Segment 1 (combined segment name):',...
        'SelectionMode','single','ListSize',[500 500],'ListString',all_exp_names);
    if tf1
        [indx2,tf2] = nmlistdlg('PromptString','Select segments to combine. Segment 2:',...
            'SelectionMode','single','ListSize',[500 500],...
            'ListString',all_exp_names);
    end
end
end