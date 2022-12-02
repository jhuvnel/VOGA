%Add GNO goggle number into old files and resegment all files
GNO_dir = unique(extractfield(dir('\\win.ad.jhu.edu\cloud\vnelhuman$\MVI\Study Subjects\MVI*R*\Visit*\vHIT\GNO'),'folder'));
for i = 1:length(GNO_dir)
    %% Initialize Folder Structure
    Path = GNO_dir{i};
    disp([num2str(i),'/',num2str(length(GNO_dir)),': ',Path])
    VOGA__makeFolders(Path,1,0);
    Raw_Path = [Path,filesep,'Raw Files'];
    Seg_Path = [Path,filesep,'Segmented Files'];
    Cyc_Path = [Path,filesep,'Cycle Averages'];
    CRF_Path = [Path,filesep,'CRFs'];
    Fig_Path = [Path,filesep,'Figures'];
    %% Remove items in folders that will be remade
    delete([Seg_Path,filesep,'*'])
    delete([CRF_Path,filesep,'*'])
    delete([Fig_Path,filesep,'*']) 
    %% Update or create notes files for all VOG files
%     txt_files = extractfield(dir([Raw_Path,filesep,'*.txt']),'name');
%     notes_files = txt_files(contains(txt_files,'-Notes.txt'));
%     VOG_files = txt_files(~contains(txt_files,'-Notes.txt'));
%     VOG_no_notes = VOG_files(~ismember(VOG_files,strrep(notes_files,'-Notes','')));
%     for j = 1:length(notes_files) %Edit existing Notes files
%         fname_notes = [Raw_Path,filesep,notes_files{j}];
%         fname_xml = strrep(fname_notes,'-Notes.txt','.xml');
%         txt_data = table2cell(readtable(fname_notes,'ReadVariableNames',false,'Delimiter','*')); %The * is a delimieter we never use so that each line goes into one cell.
%         if isfile(fname_xml)&&~contains(txt_data{contains(txt_data,'Goggle')},{'1','2','3','4'}) %Already has goggle number            
%             fdata = cellstr(readlines(fname_xml));
%             gog_line = fdata{contains(fdata,'GogglesSN')};
%             gog_num = gog_line(ismember(gog_line,'0123456789'));
%             gog_detect = find(ismember({'129852','5628','989477','3012563'},gog_num));                
%         elseif ~contains(txt_data{contains(txt_data,'Goggle')},{'1','2','3','4'}) %Already has goggle number   
%             disp(['Missing xml file: ',strrep([Raw_Path,filesep,notes_files{j}],'-Notes.txt','.xml')])
%             gog_num = inputdlg('Input Goggle # :','');
%             if isempty(gog_num)
%                 return;
%             end
%             gog_detect = gog_num{:};            
%         end
%         new_gog_line = ['Goggle GNO',num2str(gog_detect)];
%         if ~strcmp(txt_data{contains(txt_data,'Goggle')},new_gog_line)
%             txt_data{contains(txt_data,'Goggle')} = new_gog_line;
%             writecell(txt_data,fname_notes,'QuoteStrings',0)
%         end        
%     end
%     if ~isempty(VOG_no_notes) %Make missing Notes files
%         MakeNotes(Raw_Path)
%         pause; 
%         close all;
%     end    
   %% Segment
   txt_files = extractfield(dir([Raw_Path,filesep,'*.txt']),'name');
   notes_files = txt_files(contains(txt_files,'-Notes.txt'));
   VOG_files = strrep(notes_files,'-Notes','');
    for j = 1:length(VOG_files)
        Segment([Raw_Path,filesep,VOG_files{j}],Seg_Path)
    end
     %% Remake Cycle Averages
    old_cycavg = extractfield(dir([Cyc_Path,filesep,'*.mat']),'name');
    if ~isempty(old_cycavg) %Only updated CycAvg Files at all
        old_cycavg(contains(old_cycavg,{'GNO1','GNO2','GNO3','GNO4'})) = []; 
    end
    if ~isempty(old_cycavg) %Any CycAvg Files at all
        old_cycavgfig = dir([Cyc_Path,filesep,'*.fig']);
        if ~isempty(old_cycavgfig)
            delete([Cyc_Path,filesep,'*.fig']) 
        end
        if any(contains(old_cycavg,'NotAnalyzeable_'))
            delete([Cyc_Path,filesep,'NotAnalyzeable_*.mat'])
            old_cycavg = extractfield(dir([Cyc_Path,filesep,'*.mat']),'name');
        end
        new_seg = extractfield(dir([Seg_Path,filesep,'*.mat']),'name');
        no_num_seg = strrep(strrep(strrep(strrep(new_seg,'GNO4','GNO'),'GNO3','GNO'),'GNO2','GNO'),'GNO1','GNO');
        for j = 1:length(old_cycavg)
            fname = strrep(old_cycavg{j},'CycAvg_','');
            ind = contains(no_num_seg,fname);
            if sum(ind)==1
                load([Seg_Path,filesep,new_seg{ind}],'Data');
                [CycAvg,analyzed] = MakeCycAvg(Data,Cyc_Path,fname); 
                movefile([Cyc_Path,filesep,old_cycavg{j}],[Cyc_Path,filesep,'DELETEME',old_cycavg{j}])
                MakeCycAvg__saveCycAvg(Cyc_Path,strrep(CycAvg.name,'CycAvg_',''),CycAvg,analyzed);
            else
                movefile([Cyc_Path,filesep,old_cycavg{j}],[Cyc_Path,filesep,'DELETEME',old_cycavg{j}])
            end           
        end
        if ~isempty(old_cycavg)
            delete([Cyc_Path,filesep,'DELETEME*']) 
            MakeCycleSummaryTable(Path,Cyc_Path,0);
            VOGA__makePlots('Group Cycle Avg',Path)
            VOGA__makePlots('Parameterized',Path)
            close all;
        end        
    end
end