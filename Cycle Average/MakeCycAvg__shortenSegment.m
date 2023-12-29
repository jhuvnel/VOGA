function [Data,good_rng] = MakeCycAvg__shortenSegment(ha,Data,plot_info)
te = Data.te;
def = cellfun(@(x) num2str(x,10),num2cell([te(1),te(end)]),'UniformOutput',false);
good_rng = 'Redo';
axes(ha(1))
while strcmp(good_rng,'Redo')
    position = [plot_info.menu_space(1),sum(plot_info.menu_space([2,4]))-300,plot_info.menu_space(3),300];
    out_nums = cellfun(@str2double,inputdlgcol({'Start Time: ','End Time:'},...
        'X-axis Limits',[1 18],def,'on',1,position,[],'pixels'));
    if isempty(out_nums)
        good_rng = 'Done';
    else
        hold on
        h1 = xline(out_nums(1),'LineWidth',5);
        h2 = xline(out_nums(2),'LineWidth',5);
        good_rng = questdlg('Keep or redo the range?','','Keep','Redo','Redo');
        delete(h1)
        delete(h2)
        hold off
    end
end
if strcmp(good_rng,'Keep')
    [~,t1] = min(abs(te-out_nums(1)));
    [~,t2] = min(abs(te-out_nums(2)));
    fields = fieldnames(Data);
    old_size = length(te);
    for i = 1:length(fields)
        if any(size(Data.(fields{i}))==old_size)
            Data.(fields{i}) = Data.(fields{i})(t1:t2);
        end
    end   
end
end