function keep_tr = MakeCycAvg__automatedCycleRemoval(CycAvg,pass_num)
if nargin < 2
    pass_num = 1;
end
%Load intial values in
keep_tr = CycAvg.keep_tr;
cyclist = CycAvg.cyclist;
%Also isolate the time that the trigger is high for the traces of most
%interest.
rel_t = (CycAvg.stim > 0);
%Turn the 3 traces from each eye into a 3D matrix that is cycles x time x
%trace for easier vector math
all_tr(:,:,1) = CycAvg.lz_cyc;
all_tr(:,:,2) = CycAvg.rz_cyc;
all_tr(:,:,3) = CycAvg.ll_cyc;
all_tr(:,:,4) = CycAvg.rl_cyc;
all_tr(:,:,5) = CycAvg.lr_cyc;
all_tr(:,:,6) = CycAvg.rr_cyc;
NC = size(all_tr,1);
if contains(CycAvg.name,{'LA','RP'})
    rel_tr = all_tr(:,rel_t,3:4);
elseif contains(CycAvg.name,{'RA','LP'})
    rel_tr = all_tr(:,rel_t,5:6);
elseif contains(CycAvg.name,{'LH','RH'})
    rel_tr = all_tr(:,rel_t,1:2);
end
% This block removes blinks and bad tracking across multiple traces.
%Iteratively delete the cycle that is most effecting the group standard
%deviation across traces and determine cycle order to do that
sub_tr = all_tr; %This matrix will have cycles iteratively deleted
sub_cyc = cyclist; 
rm_cyc_ord = NaN(NC,1);
curr_sd = NaN(NC,1); %Save the group sd values at each step
for i = 1:NC-1
    curr_sd(i) = sum(reshape(abs(sub_tr - mean(sub_tr,1,'omitnan')),[],1));
    nc = size(sub_tr,1);
    rm_sd = NaN(nc,1);
    for j = 1:nc
        rm_sd(j) = sum(reshape(abs(sub_tr(~ismember(1:nc,j),:,:) - ...
            mean(sub_tr(~ismember(1:nc,j),:,:),1,'omitnan')),[],1))/curr_sd(i);
    end
    [~,ind] = min(rm_sd);
    rm_cyc_ord(i+1) = sub_cyc(ind);
    sub_cyc(ind) = [];
    sub_tr(ind,:,:) = [];
end
%Minimize the Standard Error of the Mean (balances SD vs cycles)
curr_sem = curr_sd./((NC:-1:1)'.^2);
[~,ind] = min(curr_sem);
if ind ~= 1  %If ind is 1, it's better not to remove any cycles
    keep_tr(rm_cyc_ord(2:ind)) = false;
    rel_tr(~keep_tr,:,:) = NaN;
end
%This block removes clear outliers in the most relevant traces
%Calculate z-scores
if pass_num == 1
    score_thresh = 0.25;
else
    score_thresh = 0.3;
end    
rel_tr = [rel_tr(:,:,1),rel_tr(:,:,2)];
rel_zscore = abs(rel_tr-mean(rel_tr,1,'omitnan'))./std(rel_tr,1,'omitnan');
rel_zscore(rel_zscore <= 1) = 0;
sum_zscore = trapz(mean(diff(CycAvg.t))*(1:2*sum(rel_t)),rel_zscore');
out_rng = sum_zscore>score_thresh;
%Update keep_tr
keep_tr(cyclist(out_rng)) = false;
end