function VOGA__CycAvg
    Path = cd;
    code_Path = [userpath,filesep,'VOGA'];
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
        done = MakeCycAvg(Path,exp_types);
    end
end