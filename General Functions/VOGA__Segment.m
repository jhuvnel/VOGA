function VOGA__Segment
    Raw_Path = [cd,filesep,'Raw Files'];
    Seg_Path = [cd,filesep,'Segmented Files'];
    %Transfer NKI Raw Files from their subfolders if they exist
    moveNKIfiles(Raw_Path)
    %Split and move GNO files into Raw Files
    moveGNOfiles(cd)
    %Detect and process log files and austoscan files
    logtoNotes(Raw_Path)
    % Select files to segment
    all_files = extractfield(dir(Raw_Path),'name',find(~extractfield(dir(Raw_Path),'isdir')));
    LDVOG_files = all_files(contains(all_files,'SESSION')&contains(all_files,'.txt')&~contains(all_files,'Notes.txt'));
    NKI_files = all_files(contains(all_files,'.dat'));
    GNO_files = all_files(contains(all_files,{'Lateral.txt','LARP.txt','RALP.txt'}));
    VOG_files = [LDVOG_files;NKI_files;GNO_files];
    if isempty(VOG_files) 
        uiwait(msgbox('No VOG files found in the Raw Files folder.'))
    elseif isempty([LDVOG_files;NKI_files]) %only GNO files
        if any(contains(all_files,'FileNotes.txt'))
            data = readtable([Raw_Path,filesep,'FileNotes.txt'],'ReadVariableNames',false);
            GNO_files = data{:,1};
            dataTypes = data{:,2};
        else            
            canals = cell(size(GNO_files));
            canals(contains(GNO_files,'Lateral')) = {'LHRH'};
            canals(contains(GNO_files,'LARP')) = {'LARP'};
            canals(contains(GNO_files,'RALP')) = {'RALP'};
            %open PDF for notes, close manually later
            pdf_files = extractfield(dir([cd,filesep,'*.pdf']),'name');
            for i = 1:length(pdf_files)
                open(pdf_files{i})
            end    
            if(isempty(pdf_files))
                disp('No PDF Reports with Experiment Notes found')
                default_note = strcat({'Manual-Impulse-'},canals,{'-NoNotes'});
            else
                default_note = strcat({'aHIT/Manual-Impulse-'},canals,{'-dpsVelandNotes'});
            end
            dataTypes = inputdlg(GNO_files,'',[1 80],default_note); 
            writetable(cell2table([GNO_files,dataTypes]),[Raw_Path,filesep,'FileNotes.txt'],'WriteVariableNames',false);
        end
        for i = 1:length(GNO_files)
            In_Path = [Raw_Path,filesep,GNO_files{i}];
            Segment(In_Path,Seg_Path,dataTypes{i})
        end
    else
        [indx,tf] = nmlistdlg('PromptString','Select files to segment:','ListSize',[300 300],'ListString',VOG_files);
        if tf == 1
            sel_files = VOG_files(indx);
            for i = 1:length(sel_files)
                In_Path = [Raw_Path,filesep,sel_files{i}];
                Segment(In_Path,Seg_Path)
            end
        end
    end
end