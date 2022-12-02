function filt = MakeCycAvg__automatedFiltering(CycAvg,pos_vel)
filt = CycAvg.filt;
prc_fun = @(trac) prctile(trac-mean(trac,'omitnan'),85);
if strcmp(pos_vel,'pos')
    LX = CycAvg.Data_filtpos.LE_Position_X;
    RX = CycAvg.Data_filtpos.RE_Position_X;        
    spline_lx = min([1-5*10^(-floor(prc_fun(LX))-3),1]);
    spline_rx = min([1-5*10^(-floor(prc_fun(RX))-3),1]);
    filt.pos.spline('LX') = spline_lx;
    filt.pos.spline('RX') = spline_rx;
elseif strcmp(pos_vel,'vel')
    trac = {'LZ','RZ','LLARP','RLARP','LRALP','RRALP'};    
    for i = 1:length(trac)
        prc_num = prc_fun(abs(reshape([CycAvg.([lower(trac{i}(1:2)),'_cyc'])],[],1)));
        %disp([trac{i},': ',num2str(prc_num)])
        if prc_num < 20
            filt.vel.irlssmooth(trac{i}) = 25-floor(prc_num);
        else
            filt.vel.spline(trac{i}) = min([1-5*10^(-floor(prc_num)-2),1]);
        end
    end
end
end