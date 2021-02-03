%% Save
function MakeCycAvg__saveCycAvg(Cyc_Path,In_FileName,CycAvg)
    save([Cyc_Path,filesep,'CycAvg_',In_FileName],'CycAvg')
    %If the file was previously designated as un-analyzeable, remove that file
    d1 = extractfield(dir([Cyc_Path,filesep,'*.mat']),'name');
    if ismember(['NotAnalyzeable_',In_FileName],d1)
        delete([Cyc_Path,filesep,'NotAnalyzeable_',In_FileName])
    end
end