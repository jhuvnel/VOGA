%% Autoscan Param Suggestion
%Use the MVI Stim Currents Available as of 2023-05
MVI_curr = [0:1.18:125.08,127.44:2.36:302.8,306.8:4.72:604.16,613.6:9.44:699];
load('AutoscanParameters.mat')
dat = sortrows(sortrows(dat,'MaxScore','descend'),'Canal','ascend');
out_vars = {'BaselineCurr','MaxCurr','StimSteps'};
tab = [dat(:,{'Canal','Cond','MaxScore'}),...
    array2table(NaN(size(dat,1),length(out_vars)),'VariableNames',out_vars)];
tab.MaxScore = round(cell2mat(tab.MaxScore),0);
tab.Cond = strrep(tab.Cond,'\mu','u');
for i = 1:size(tab,1)
    curr = dat.CurrentAmp{i};
    vel = dat.MaxVel{i};
    align = dat.Align{i};
    score = dat.Score{i};
    i_s = find([0;vel]<5,1,'last');
    %Only make predictions if at least two points are above the noise floor
    if sum(vel>5)>1
        %Set the maximum current as the last point within 90% of the
        %velocity and misalignemnt of the point with the highest score
        %Allows for some noise but removes current values with increased
        %curernt spread leading to misalignemnt
        [~,ind] = max(score);
        i_e = find(vel>0.9*vel(ind)&0.9*align<align(ind),1,'last');
        tab.MaxCurr(i) = curr(i_e);
        %Set the baseline current (1/3 of max vel) by linear smoothing
        %first
        coeff = [ones(length(i_s:i_e),1),vel(i_s:i_e)]\curr(i_s:i_e);
        curr_val = coeff(1)+coeff(2)*max(vel(i_s:i_e))/3;
        [~,ind_c] = min(abs(MVI_curr-curr_val));
        [~,ind_m] = min(abs(MVI_curr-curr(i_e)));
        %Find the closes value delieverable by the MVI system
        tab.BaselineCurr(i) = round(MVI_curr(ind_c),0);
        tab.StimSteps(i) = ind_m-ind_c+1;
    end
end
disp(tab)