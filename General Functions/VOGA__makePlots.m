function VOGA__makePlots
    code_Path = [userpath,filesep,'VOGA'];
    params.Path = cd;
    params.Raw_Path = [cd,filesep,'Raw Files'];
    params.Seg_Path = [cd,filesep,'Segmented Files'];
    params.Cyc_Path = [cd,filesep,'Cycle Averages'];
    params.code_Path = code_Path;
    % Get version and experimenter info from the file
    if ~any(contains(extractfield(dir(code_Path),'name'),'VerInfo.txt'))
        writeInfoFile(code_Path);
    end
    data = readtable([code_Path,filesep,'VerInfo.txt'],'ReadVariableNames',false);
    params.version = data{1,2}{:};
    params.Experimenter = data{2,2}{:};    
    %% Run until cancel 
    opts = {'Raw VOG','Segment','Cycle Average','Group Cycle Avg','Parameterized'};    
    [ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
                       'SelectionMode','single',...
                       'ListSize',[150 125],...
                       'ListString',opts);  
    while tf
        if strcmp(opts{ind},'Raw VOG')
            in_args = inputdlg({'Plot eye movements (0/1): ','Plot head movement in XYZ or LRZ (xyz/lrz):'},'Plot settings',[1 50],{'0','lrz'});
            if ~isempty(in_args)
                plot_eyes = str2double(in_args{1});
                lrz_xyz = in_args{2};
                plotRawVOG(params.Raw_Path,plot_eyes,lrz_xyz) %Select file inside this function
            end
        elseif strcmp(opts{ind},'Segment')
            %Select files first
            all_files = extractfield(dir([params.Seg_Path,filesep,'*.mat']),'name');
            [indx,tf2] = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',all_files,'SelectionMode','multiple');
            if tf2
               for i = 1:length(indx)
                    load([params.Seg_Path,filesep,all_files{indx(i)}],'Data')
                    plotSegment(Data); 
               end
            end
        elseif strcmp(opts{ind},'Cycle Average')
            in_args = inputdlg({'Plot fits (0/1): ','Plot head movement in XYZ or LRZ (xyz/lrz):'},'Plot settings',[1 50],{'0','lrz'});
            if ~isempty(in_args)
                plot_fits = str2double(in_args{1});
                lrz_xyz = in_args{2};
                %Select files first
                all_files = extractfield(dir([params.Cyc_Path,filesep,'*Cyc*.mat']),'name');
                [indx,tf2] = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',all_files,'SelectionMode','multiple');
                if tf2
                   for i = 1:length(indx)
                        load([params.Cyc_Path,filesep,all_files{indx(i)}],'CycAvg')
                        plotCycAvg(CycAvg,plot_fits,lrz_xyz) ; 
                   end
                end
            end
        elseif strcmp(opts{ind},'Group Cycle Avg')  
            in_args = inputdlg({'Descriptive annotation on figure (0/1): ','Y-axis maximum (dps). Leave empty for default:'},'Plot settings',[1 50],{'1',''});
            if ~isempty(in_args)
                params.annot = str2double(in_args{1});
                params.YMax = str2double(in_args{2});
                plotGroupCycAvg(params);
            end          
        elseif strcmp(opts{ind},'Parameterized')  
            types = {'SineAmpVelLRZ','SineAmpVelXYZ','Autoscan','SpherePlot','SineFreqGainPhase'}; %Update as more things get added
            [indx,tf2] = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',types,'SelectionMode','multiple');
            if tf2
                in_args = inputdlg({'Descriptive annotation on figure (0/1): ','Y-axis maximum (dps). Leave empty for default:'},'Plot settings',[1 50],{'1',''});
                if ~isempty(in_args)
                    params.annot = str2double(in_args{1});
                    params.YMax = str2double(in_args{2});
                    for i = 1:length(indx)
                        plotParamResults(types{indx(i)},params);
                    end
                end
            end 
        end
        [ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
                   'SelectionMode','single',...
                   'ListSize',[150 125],...
                   'ListString',opts);  
    end    
end