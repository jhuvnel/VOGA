% Ultimately can be used for a variety of things but now focused on
% adjusting the trigger
close all;
LRZ_XYZ = 1;
fname = 'CycAvg_MVI001R019-Visit13x-20210210-eeVOR-Sine-Y-2Hz-500dps.mat';
load(fname,'CycAvg')
plotCycAvg(CycAvg,1,0,LRZ_XYZ) 
set(gca,'YLim',[-100 100])
[x,~] = ginput(1);
[~,x_i] = min(abs(CycAvg.t - x));

traces = {'lz','ll','lr','lx','ly','rz','rl','rr','rx','ry'};
for i = 1:length(traces)
    CycAvg.([traces{i},'_cycavg']) = CycAvg.([traces{i},'_cycavg'])([x_i:length(CycAvg.t),1:x_i-1]);
    CycAvg.([traces{i},'_cycstd']) = CycAvg.([traces{i},'_cycstd'])([x_i:length(CycAvg.t),1:x_i-1]);
end

plotCycAvg(CycAvg,1,0,LRZ_XYZ) 
set(gca,'YLim',[-100 100])
save(fname,'CycAvg')