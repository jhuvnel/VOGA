function OLD_CombineResultsTables(subs)
    %subs is an optional array to combine only a subset of the subjects
    %Update as needed
    all_sub_folders = {'MVI001_R019',...
                       'MVI002_R004',...
                       'MVI003_R140',...
                       'MVI004_R201',...
                       'MVI005_R107',...
                       'MVI006_R296',...
                       'MVI007_R765',...
                       'MVI008_R021'};
    if nargin == 1
        paths = all_sub_folders(subs);
    else
        paths = all_sub_folders;
    end
    in_path = uigetdir('Select the Study Subjects Folder.');
    all_results = table();
    for i = 1:length(paths)
        fstruct = dir([in_path,filesep,paths{i},filesep,'*RotaryChairResults.mat']);
        if ~isempty(fstruct)
            fnames = {fstruct.name}';
            if length(fnames)>1
                indx = listdlg('PromptString',['Select which file to use for ',paths{i}],...
                    'ListString',fnames,'ListSize',[300 300]);
                fnames = fnames(indx);
            end
            a = load([in_path,filesep,paths{i},filesep,fnames{:}],'all_results');
            tab = a.all_results;
            all_results = [all_results;tab];
        else
            disp(['No RotaryChairResults.mat file found in ',in_path,filesep,paths{i}])
        end
    end
    %Save
    out_path = uigetdir('Select where to save this file.');
    fname = inputdlg('File Name (no extension):');
    if ~isempty(fname)
        save([out_path,filesep,fname{:},'.mat'],'all_results')
    end 
end