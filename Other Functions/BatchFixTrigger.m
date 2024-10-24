% Batch trigger fixing
trigger_shift = 15;
Path = cd;
Seg_Path = [Path,filesep,'Segmented Files'];
Cyc_Path = [Path,filesep,'Cycle Averages'];
progress_tab = assessProgress(Path);
%progress_tab(contains(progress_tab.Segment,{'NKI','NL'}),:) = []; %Don't updated the Neurolign files
all_seg = progress_tab.Segment;
for j = 1:length(all_seg)
    fname = all_seg{j};
    load([Seg_Path,filesep,fname],'Data');
    Data.info.TriggerShift2 = trigger_shift;
    save([Seg_Path,filesep,fname],'Data')
    if progress_tab.CycAvgFile(j) %Reanalyze
        [CycAvg,analyzed] = MakeCycAvg(Data,Cyc_Path,{'Auto Rerun'},1);
        if analyzed                    
            MakeCycAvg__saveCycAvg(Cyc_Path,strrep(CycAvg.name,'CycAvg_',''),CycAvg,analyzed,1);
            if ~strcmp(CycAvg.name,['CycAvg_',fname]) %Had a different name, delete file
                delete([Cyc_Path,filesep,'CycAvg_',fname])
            end
        end
    end
end
VOGA__SummaryTable('Folder (Load)')
VOGA__makePlots('Parameterized')