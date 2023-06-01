%% Save a CycAvg or NonAnalzeable File
function MakeCycAvg__saveCycAvg(Cyc_Path,fname,CycAvg,analyzed,save_fig)
if analyzed
    filename = [Cyc_Path,filesep,'CycAvg_',fname];
    save([Cyc_Path,filesep,'CycAvg_',fname],'CycAvg')    
    rm_file = {[Cyc_Path,filesep,'NotAnalyzeable_',fname];...
        strrep([Cyc_Path,filesep,'NotAnalyzeable_',fname],'.mat','.fig')};      
else   
    filename = [Cyc_Path,filesep,'NotAnalyzeable_',fname];
    Data = CycAvg;
    save([Cyc_Path,filesep,'NotAnalyzeable_',fname],'Data')
    rm_file = {[Cyc_Path,filesep,'CycAvg_',fname];...
        strrep([Cyc_Path,filesep,'CycAvg_',fname],'.mat','.fig')};  
end
%This will save on top of an existing NonAnalyzable/CycAvg file
%If the file was previously designated as a different type, remove that file
for i = 1:length(rm_file)
    if isfile(rm_file{i})
        delete(rm_file{i})
    end
end
if save_fig
    savefig(strrep(filename,'.mat','.fig'))
    close;
end
end