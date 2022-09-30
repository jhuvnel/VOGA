function plotSummaryFigures(params)
opts = {'Rotary Chair (Sine)','vHIT (GNO)',...
    'Rotary Chair (Sine) + Candidate','vHIT (GNO) + Candidate'};    
[ind,tf] = nmlistdlg('PromptString','Select an plot to make:',...
                   'SelectionMode','multiple',...
                   'ListSize',[200 125],...
                   'ListString',opts);
if tf
    
    for i = 1:length(ind)
        switch opts{ind(i)}
            
        end
    end
end                   
end