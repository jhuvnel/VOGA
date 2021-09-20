%% Save a CycAvg or NonAnalzeable File
function MakeCycAvg__saveCycAvg(Cyc_Path,fname,CycAvg,analyzed)
if analyzed
    savefig([Cyc_Path,filesep,'CycAvg_',fname(1:end-4),'.fig'])
    close;
    save([Cyc_Path,filesep,'CycAvg_',fname],'CycAvg')
    %If the file was previously designated as un-analyzeable, remove that file
    d1 = extractfield(dir([Cyc_Path,filesep,'*.mat']),'name');
    if ismember(['NotAnalyzeable_',fname],d1)
        delete([Cyc_Path,filesep,'NotAnalyzeable_',fname])
    end
else   
    savefig([Cyc_Path,filesep,'NotAnalyzeable_',fname(1:end-4),'.fig'])
    close;
    Data = CycAvg;
    save([Cyc_Path,filesep,'NotAnalyzeable_',fname],'Data')
    %If the file was previously designated as analyzeable, remove that file
    d1 = extractfield(dir([Cyc_Path,filesep,'*.mat']),'name');
    if ismember(['CycAvg_',fname],d1)
        delete([Cyc_Path,filesep,'CycAvg_',fname])
    end
end
end