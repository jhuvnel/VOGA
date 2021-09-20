function VOGA__CycAvg
    Path = cd;
    Seg_Path = [Path,filesep,'Segmented Files'];
    Cyc_Path = [Path,filesep,'Cycle Averages'];
    %code_Path = [userpath,filesep,'VOGA'];
    done = false;
    % Experiment types
    progress_tab = assessProgress(Path);
    all_exp_names = progress_tab.Segment;
    for i = 1:length(all_exp_names)
        n_spl = strsplit(all_exp_names{i},'-');
        if length(n_spl{3})>8 %If time is included with the date
            n_spl{3} = n_spl{3}(1:8);
        end
        if ~isnan(str2double(n_spl{4})) %If time is seperate
            n_spl(4) = [];
        end
        all_exp_names(i) = join(n_spl(3:find(contains(n_spl,{'RotaryChair','eeVOR','aHIT','Manual'}))+1),'-');
    end
    exp_names = unique(all_exp_names);
    if length(exp_names)>1
        [indx,tf] = nmlistdlg('PromptString','Select experiment types to analyze:',...
                       'SelectionMode','multiple',...
                       'ListSize',[350 300],...
                       'ListString',exp_names); 
    else
        tf = 0;
    end
    if tf == 0
        exp_types = {};
    else
        exp_types = exp_names(indx);
    end
    while(~done) %run until the user hits cancel on analyzing a file
        %% Select file
        progress_tab = assessProgress(Path);
        if ~isempty(exp_types)
            rel_file = false(length(progress_tab{:,1}),length(exp_types));
            for i = 1:length(exp_types)
                components = split(exp_types{i},'-');
                cont = true(length(progress_tab{:,1}),1);
                for j = 1:length(components)
                    cont = cont&contains(progress_tab{:,1},components{j});
                end
                rel_file(:,i) = cont;
            end
            progress_tab(~any(rel_file,2),:) = [];
        end
        progress_i = [find(~progress_tab{:,2}&~progress_tab{:,3});find(progress_tab{:,2}|progress_tab{:,3})]; %put unanalyzed files at the top
        %Change color of option depending on whether it's been analyzed
        font_col = repmat({'black'},length(progress_i),1);
        font_col(progress_tab{progress_i,2}) = {'green'};
        font_col(progress_tab{progress_i,3}) = {'red'};
        list = strcat('<HTML><FONT color="',font_col,'">',table2cell(progress_tab(progress_i,1)),'</FONT></HTML>');
        [indx,tf] = nmlistdlg('PromptString','Unattempted = Black, Analyzed = Green, Not Analyzeable = Red. Select a file to analyze:',...
                               'SelectionMode','single',...
                               'ListSize',[500 600],...
                               'ListString',list);  
        % Run unless the user selects cancel
        if tf
            fname = progress_tab{progress_i(indx),1}{:};
            load([Seg_Path,filesep,fname],'Data');
            [CycAvg,analyzed] = MakeCycAvg(Data);
            if ~isempty(CycAvg)
                MakeCycAvg__saveCycAvg(Cyc_Path,fname,CycAvg,analyzed);
            end
        else 
            disp('Operation Ended.')
            done = true;
        end
    end
end