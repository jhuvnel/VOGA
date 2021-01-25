function cyc_num = howmanycycles
    path = uigetdir('','Select the folder with the Cycle Average files of interest.');
    if ispc
        sep = '\';
    else
        sep = '/';
    end
    %% Load in data
    a = extractfield(dir([path,sep,'*.mat']),'name')';
    if isempty(a)
        error('No segments present.')
    end
    disp(a)
    rm_update = questdlg('Remove files that were not updated?','','Yes','No','No');
    if strcmp(rm_update,'Yes')
        a(~contains(a,'Updated')) = [];
    end
    [indx,tf] = listdlg('PromptString','Select which visits to plot:',...
                               'SelectionMode','multiple',...
                               'ListString',a,...
                               'ListSize',[600 300]);
    if tf == 0
        error('No cycles selected')
    end
    fnames = a(indx);
    cycs = NaN(length(fnames),1);
    for i = 1:length(fnames)
        b = load([path,sep,fnames{i}]);
        c = fieldnames(b);
        d = b.(c{1});
        [n1,n2] = size(d.stim);
        [n3,n4] = size(d.ll_cycavg);
        if n3 ~= 1 && n4 ~= 1
            disp('LL CycAvg was unexpectedly a 2D array with dimensions:')
            disp([n3,n4])
        elseif n3 == 1
            if n1 ~= n4 && n2 ~= n4
                disp(['Expect time array to be ',num2str(n4),' samples long.'])
                disp('Stim array was below dimensions instead.')
                disp([n1,n2])
            elseif n1 == n4
                cycs(i) = n2;
            elseif n2 == n4
                cycs(i) = n1;
            end
        elseif n4 == 1
            if n1 ~= n3 && n2 ~= n3
                disp(['Expect time array to be ',num2str(n3),' samples long.'])
                disp('Stim array was below dimensions instead.')
                disp([n1,n2])
            elseif n1 == n3
                cycs(i) = n2;
            elseif n2 == n3
                cycs(i) = n1;
            end
        end
    end
    cyc_num = [fnames,num2cell(cycs)];
end