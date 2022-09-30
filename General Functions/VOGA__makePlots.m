function VOGA__makePlots
    code_Path = [userpath,filesep,'VOGA'];
    params.Path = cd;
    params.code_Path = code_Path;
    % Get version and experimenter info from the file
    if ~any(contains(extractfield(dir(userpath),'name'),'VOGA_VerInfo.txt'))
        VOGA__setVersion;
    end
    data = readtable('VOGA_VerInfo.txt','ReadVariableNames',false);
    params.version = data{1,2}{:};
    params.Experimenter = data{2,2}{:}; 
    params.annot = 0;
    params.YMax = NaN;
    %% Run until cancel 
    opts = {'Raw VOG','Segment','Cycle Average','Group Cycle Avg',...
        'Parameterized','Across Subjects','Sphere Plot'};    
    [ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
                       'SelectionMode','single',...
                       'ListSize',[150 125],...
                       'ListString',opts);  
    while tf
        %Expect to be in a visit folder with the right structure
        if ~VOGA__checkFolders(0)&&contains(opts{ind},{'Raw VOG','Segment','Cycle Average','Group Cycle Avg','Parameterized','Sphere Plot'})
            error('Expected folder structure not present. Navigate to appropriate folder before trying again.')
        elseif contains(opts{ind},{'Raw VOG','Segment','Cycle Average','Group Cycle Avg','Parameterized','Sphere Plot'})
            params.Raw_Path = [cd,filesep,'Raw Files'];
            params.Seg_Path = [cd,filesep,'Segmented Files'];
            params.Cyc_Path = [cd,filesep,'Cycle Averages'];
        end
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
            in_args = inputdlg({'Plot fits (0/1): ','Plot head movement in XYZ or LRZ (xyz/lrz):'},'Plot settings',[1 50],{'0',''});
            if ~isempty(in_args)
                plot_fits = str2double(in_args{1});
                lrz_xyz = in_args{2};
                %Select files first
                all_files = extractfield(dir([params.Cyc_Path,filesep,'*Cyc*.mat']),'name');
                [indx,tf2] = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',all_files,'SelectionMode','multiple');
                if tf2
                   for i = 1:length(indx)
                        load([params.Cyc_Path,filesep,all_files{indx(i)}],'CycAvg')
                        plotCycAvg(CycAvg,plot_fits,lrz_xyz); 
                   end
                end
            end
        elseif strcmp(opts{ind},'Group Cycle Avg')              
            plotGroupCycAvg(params);            
        elseif strcmp(opts{ind},'Parameterized')  
            plotParamResults(params);
        elseif strcmp(opts{ind},'Across Subjects') 
            plotSummaryFigures(params);
        elseif strcmp(opts{ind},'Sphere Plot')  
            in_args = inputdlg({'Descriptive annotation on figure (0/1): '},'Plot settings',[1 50],{'1'});
            if ~isempty(in_args)
                params.annot = str2double(in_args{1});
                plotSpherePlot(params);
            end    
        end
        [ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
                   'SelectionMode','single',...
                   'ListSize',[150 125],...
                   'ListString',opts);  
    end    
end