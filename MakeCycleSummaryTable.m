%% MakeCycleSummaryTable
%It can take in arguments of the path to save the table and the path where
%all the cycle averages are located. If not provided, it will allow the
%user to select cycle average directories manually or provide a text file
%with the directories and assume the table should be saved in the current
%directory.
%This version takes in all CycAvg types that have been defined and saves a
%table for each

function MakeCycleSummaryTable(out_path,cyc_path,rerun)
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
        end
    end       
   if nargin < 3
      rerun = 1; 
   end
    %% Now analyze each file and get the table
    %Initialize table, assumes 1 CycAvg = 1 row so no sumsine
%     all_results = table('Size',[length(files),185],...
%         'VariableTypes',[repmat("cell",3,1);"datetime";repmat("cell",4,1);repmat("double",185-8,1)]);
%     traces = {'lz','rz','ll','rl','lr','rr','lx','rx','ly','ry'};
%     labs = [{'File';'Subject';'Visit';'Date';'Goggle';'Experiment';'Condition';'Frequency';'Cycles'};...
%     reshape(strcat('MaxVel_',repmat(upper(traces),4,1),repmat({'_HIGH';'_HIGH_sd';'_LOW';'_LOW_sd'},1,length(traces))),[],1);...
%     reshape(strcat('Gain_',repmat(upper(traces),4,1),repmat({'_HIGH';'_HIGH_sd';'_LOW';'_LOW_sd'},1,length(traces))),[],1);...
%     reshape(strcat('Tau_',repmat(upper(traces),4,1),repmat({'_HIGH';'_HIGH_sd';'_LOW';'_LOW_sd'},1,length(traces))),[],1);...
%     reshape(strcat('RMSE_',repmat(upper(traces),2,1),repmat({'_HIGH';'_LOW'},1,length(traces))),[],1);...
%     reshape(strcat('Latency_',repmat(upper(traces),2,1),repmat({'';'_sd'},1,length(traces))),[],1);...
%     {'Phase_L';'Phase_L_sd';'Phase_R';'Phase_R_sd'};...
%     {'Align_L_HIGH';'Align_L_HIGH_sd';'Align_L_LOW';'Align_L_LOW_sd';'Align_R_HIGH';'Align_R_HIGH_sd';'Align_R_LOW';'Align_R_LOW_sd'};...
%     {'Disc_HIGH';'Disc_HIGH_sd';'Disc_LOW';'Disc_LOW_sd'}];
%     all_results.Properties.VariableNames = labs;
    % Fill table
    tabs = cell(length(files),1);
    for i = 1:length(files)
       disp(i)
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