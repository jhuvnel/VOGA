function filt = MakeCycAvg__autoFiltCycSelect(Data,filt,plot_info)
%First time running this
if ~isempty(Data.detec_tr) %remove erroneous head traces using auto-detected traces as a template
    out_of_bounds = abs(Data.stims-mean(Data.stims(:,Data.detec_tr),2))/std(Data.stims(:,Data.detec_tr),[],2)>3;
    filt.keep_tr = ~any(out_of_bounds(Data.t_snip<0.23&Data.t_snip>0.1,:))&any(Data.stims(Data.t_snip<0.2,:)<10)&any(Data.stims(Data.t_snip<0.4&Data.t_snip>0.1,:)<0);
end
%First, filter position
filt.pos.lowpass(end) = 5;
%filt.pos.irlssmooth(end) = 1; %Remove one-point outliers
CycAvg = MakeCycAvg__filterTraces(Data,filt);



ha = MakeCycAvg__plotFullCycAvg([],CycAvg,plot_info);







%%

if type == 1
    filt.vel.irlssmooth(end) = round(length(Data.t_snip)*0.16); %heuristic
end
    trac = {'LZ','RZ','LLARP','RLARP','LRALP','RRALP'};    
    for i = 1:length(trac)
        prc_num = prc_fun(abs(reshape([CycAvg.([lower(trac{i}(1:2)),'_cyc'])],[],1)));
        if prc_num < 20
            filt.vel.irlssmooth(trac{i}) = 25-floor(prc_num);
        else
            filt.vel.spline(trac{i}) = min([1-5*10^(-floor(prc_num)-2),1]);
        end
    end
end