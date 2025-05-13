function TrimSegment
Path = cd;
Seg_Path = [Path,filesep,'Segmented Files'];
progress_tab = assessProgress(Path);
exp_names = progress_tab.Segment;
[indx,tf] = nmlistdlg('PromptString','Select experiment types to analyze:','SelectionMode','single','ListSize',[350 300],'ListString',exp_names);  
if ~tf
    return;
end
fname = [Seg_Path,filesep,exp_names{indx}];
load(fname,'Data');
fig = plotSegment(Data);
t = Data.Time_Eye - Data.Time_Eye(1);
tot = t(end)-t(1);
set(fig,'units','inches','position',[0 0 13 9])
set(gca,'XLim',[t(1)-0.1*tot, t(end)+0.1*tot])
good_rng = 'Redo';
def = {num2str(t(1)),num2str(t(end))};
while strcmp(good_rng,'Redo')
    t_bound = inputdlg({'Start Time: ','End Time:'},'Plot settings',[1 50],def); 
    if isempty(t_bound)
       good_rng = 'Done';
    else
        hold on
        h1 = xline(str2num(t_bound{1}),'LineWidth',10);
        h2 = xline(str2num(t_bound{2}),'LineWidth',10);
        good_rng = questdlg('Keep or redo the range?','','Keep','Redo','Redo');  
        delete(h1)
        delete(h2)
        def = t_bound;
        hold off
    end 
end
if strcmp(good_rng,'Keep')
    [~,t1] = min(abs(t-str2num(t_bound{1})));
    [~,t2] = min(abs(t-str2num(t_bound{2})));
    fields = fieldnames(Data);
    old_size = length(t);
    for i = 1:length(fields)
        if any(size(Data.(fields{i}))==old_size)
            Data.(fields{i}) = Data.(fields{i})(t1:t2);
        end
    end
    save(fname,'Data')
end
close(fig);
end