%% MakeCycleSineFitTable
%This function makes binocular,
% (1)half-cycle gain fits in all three canal directions,
% (2)phase estimates for the dominant canal direcion,
% (3)half-cycle misalignments, and
% (4)standard deviations for all the above measurements based on individual
% cycle values.
% It can work on single frequency files or sum of sines.
%It can take in arguments of the path to save the table and the path where
%all the cycle averages are located. If not provided, it will allow the
%user to select cycle average directories manually or provide a text file
%with the directories and assume the table should be saved in the current
%directory.

function MakeCycleSineFitTable(out_path,cyc_path)
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
                    fstruct = dir([path,filesep,'*Sin*.mat']);
                    if ~isempty(fstruct)
                        fnames = {fstruct.name}';
                        files = [files;strcat(path,filesep,fnames)];
                        disp([num2str(length(fnames)),' files found.'])
                    else
                        disp('No Cycle Average Files found in this directory')
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
                fstruct = dir([paths{i},filesep,'*Sin*.mat']);
                if ~isempty(fstruct)
                    fnames = {fstruct.name}';
                    files = [files;strcat(paths{i},filesep,fnames)];
                    files(contains(files,'NotAnalyzeable')) = []; %remove un-analyzeable files
                    disp([num2str(length(fnames)),' files found.'])
                else
                    disp('No Cycle Average Files found in this directory')
                end
            end    
        end
    else
        fstruct = dir([cyc_path,filesep,'*Sin*.mat']);
        if ~isempty(fstruct)
            fnames = {fstruct.name}';
            files = [files;strcat(cyc_path,filesep,fnames)];
            files(contains(files,'NotAnalyzeable')) = []; %remove un-analyzeable files
            disp([num2str(length(fnames)),' files found.'])
        else
            disp('No Cycle Average Files found in this directory')
        end
    end   
    %% Now analyze each file and get the table
    all_results = table();
    disp(length(files))
    k = 1;
    for i = 1:length(files)
       disp(i)
       a = load(files{i});
       b = fieldnames(a);
       CycAvg = a.(b{1});
       if ~isfield(CycAvg,'parameterized')||~isfield(CycAvg,'lz_cycavg_fit')
           [results,CycAvg] = ParameterizeSineFits(files{i},CycAvg);
           save(files{i},'CycAvg')
       else
           results = CycAvg.parameterized;
       end
       all_results(k:k-1+size(results,1),:) = results;
       k = k+size(results,1);
    end   
    %% Save
    save([out_path,filesep,datestr(now,'yyyymmdd_HHMMSS'),'_RotaryChairResults.mat'],'all_results')
end