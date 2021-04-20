function VOGA__makePlots(code_Path)
    Path = cd;
    Raw_Path = [Path,filesep,'Raw Files'];
    Seg_Path = [Path,filesep,'Segmented Files'];
    Cyc_Path = [Path,filesep,'Cycle Averages'];
    % Get version and experimenter info from the file
    if ~any(contains(extractfield(dir(code_Path),'name'),'VerInfo.txt'))
        writeInfoFile(code_Path);
    end
    data = readtable([code_Path,filesep,'VerInfo.txt'],'ReadVariableNames',false);
    version = data{1,2}{:};
    Experimenter = data{2,2}{:};
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
                plot_eyes = str2num(in_args{1});
                lrz_xyz = in_args{2};
                plotRawVOG(Raw_Path,plot_eyes,lrz_xyz) %Select file inside this function
            end
        elseif strcmp(opts{ind},'Segment')
            %Select files first
            all_files = extractfield(dir([Seg_Path,filesep,'*.mat']),'name');
            [indx,tf2] = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',all_files,'SelectionMode','multiple');
            if tf2
               for i = 1:length(indx)
                    load([Seg_Path,filesep,all_files{indx(i)}],'Data')
                    plotSegment(Data); 
               end
            end
        elseif strcmp(opts{ind},'Cycle Average')
            in_args = inputdlg({'Plot fits (0/1): ','Plot head movement in XYZ or LRZ (xyz/lrz):'},'Plot settings',[1 50],{'0','lrz'});
            if ~isempty(in_args)
                plot_fits = str2num(in_args{1});
                lrz_xyz = in_args{2};
                %Select files first
                all_files = extractfield(dir([Cyc_Path,filesep,'*Cyc*.mat']),'name');
                [indx,tf2] = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',all_files,'SelectionMode','multiple');
                if tf2
                   for i = 1:length(indx)
                        load([Cyc_Path,filesep,all_files{indx(i)}],'CycAvg')
                        plotCycAvg(CycAvg,plot_fits,lrz_xyz) ; 
                   end
                end
            end
        elseif strcmp(opts{ind},'Group Cycle Avg')  
            types = {'SineFreq','SineAmp','SineManual','Autoscan'}; %Update as more things get added
            [indx,tf2] = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',types,'SelectionMode','single');
            if tf2
                in_args = inputdlg({'Descriptive annotation on figure (0/1): ','Y-axis maximum (dps). Leave empty for default:'},'Plot settings',[1 50],{'1',''});
                if ~isempty(in_args)
                    annot = str2num(in_args{1});
                    YMax = str2num(in_args{2});
                    plotGroupCycAvg(types{indx},Path,Cyc_Path,code_Path,version,Experimenter,annot,YMax);
                end
            end            
        elseif strcmp(opts{ind},'Parameterized')  
            types = {'SineAmpVelLRZ','SineAmpVelXYZ','Autoscan','SpherePlot'}; %Update as more things get added
            [indx,tf2] = nmlistdlg('PromptString','Select files to plot:','ListSize',[300 300],'ListString',types,'SelectionMode','single');
            if tf2
                in_args = inputdlg({'Descriptive annotation on figure (0/1): ','Y-axis maximum (dps). Leave empty for default:'},'Plot settings',[1 50],{'1',''});
                if ~isempty(in_args)
                    annot = str2num(in_args{1});
                    YMax = str2num(in_args{2});
                    plotParamResults(types{indx},Path,code_Path,version,Experimenter,annot,YMax);
                end
            end 
        end
        [ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
                   'SelectionMode','single',...
                   'ListSize',[150 125],...
                   'ListString',opts);  
    end    
end