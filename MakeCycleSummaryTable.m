%% MakeCycleSummaryTable
%It can take in arguments of the path to save the table and the path where
%all the cycle averages are located. If not provided, it will allow the
%user to select cycle average directories manually or provide a text file
%with the directories and assume the table should be saved in the current
%directory.
%This version takes in all CycAvg types that have been defined and saves a
%table for each

function all_results = MakeCycleSummaryTable(out_path,cyc_path,rerun)
    if nargin < 1
        out_path = cd;
    end
    files = {};
    if nargin < 2
        % Get all of the paths first
        p_meth = nmquestdlg('Manually select cycle average folders or load text file with paths?','','Manual','File','Manual');
        if strcmp(p_meth,'Manual')
            proceed = 'Yes';
            while strcmp(proceed,'Yes')
                path = uigetdir('','Select the folder with the saved Cycle Averages');
                if isnumeric(path)
                    break; %Assume this means they want to exit the loop
                else
                    fstruct = dir([path,filesep,'*.mat']);
                    if ~isempty(fstruct)
                        fnames = {fstruct.name}';
                        files = [files;strcat(path,filesep,fnames)];
                        disp([num2str(length(fnames)),' files found.'])
                    else
                        disp(['No Cycle Average Files found in this directory: ',path])
                    end
                end
                proceed = nmquestdlg('Add more Cycle Average Files?','','Yes','No','No');
            end
        elseif strcmp(p_meth,'File')
            [path2,path1] = uigetfile('*.txt','Select the text file with all the Cycle Average paths.');
            if filesep=='\'
                paths = strcat(path1,strrep(importdata([path1,filesep,path2]),'/','\'));
            else
                paths = strcat(path1,strrep(importdata([path1,filesep,path2]),'\','/'));
            end
            for i = 1:length(paths)
                fstruct = dir([paths{i},filesep,'*.mat']);
                if ~isempty(fstruct)
                    fnames = {fstruct.name}';
                    files = [files;strcat(paths{i},filesep,fnames)];
                    files(contains(files,'NotAnalyzeable')) = []; %remove un-analyzeable files
                    disp([num2str(length(fnames)),' files found.'])
                else
                    disp(['No Cycle Average Files found in this directory: ',paths{i}])
                end
            end    
        end
    else
        fstruct = dir([cyc_path,filesep,'*.mat']);
        if ~isempty(fstruct)
            fnames = {fstruct.name}';
            files = [files;strcat(cyc_path,filesep,fnames)];
            files(contains(files,'NotAnalyzeable')) = []; %remove un-analyzeable files
            disp([num2str(length(fnames)),' files found.'])
        else
            disp(['No Cycle Average Files found in this directory: ',cyc_path])
            all_results = [];
            return;
        end
    end       
   if nargin < 3
      rerun = 1; 
   end
    %% Now analyze each file and get the table
    tabs = cell(length(files),1);
    for i = 1:length(files)
       disp([cyc_path,': ',num2str(i),'/',num2str(length(files))])
       a = load(files{i});
       b = fieldnames(a);
       CycAvg = a.(b{1});
       slash = strfind(files{i},filesep); 
       fname = files{i}(slash(end)+1:end);
       if ~isfield(CycAvg,'name')||~strcmp(CycAvg.name,fname)||~isfield(CycAvg,'lx_cyc_fit')||rerun
           CycAvg.name = fname;
           CycAvg = ParameterizeCycAvg(CycAvg);
           save(files{i},'CycAvg')
       end
       tabs{i} = CycAvg.parameterized;       
    end  
    all_results = vertcat(tabs{:});
    %% Save
    exp_type = all_results.Experiment{1};
    save([out_path,filesep,datestr(now,'yyyymmdd_HHMMSS'),'_',exp_type,'Results.mat'],'all_results')
end