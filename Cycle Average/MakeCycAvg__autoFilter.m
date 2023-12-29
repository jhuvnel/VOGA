function filt = MakeCycAvg__autoFilter(Data,filt,plot_info)
filt.pos{:,:} = NaN;
filt.vel{:,:} = NaN;
type = Data.info.type;
gog = Data.info.goggle_ver;
traces_pos = plot_info.traces_pos;
traces_vel = plot_info.traces_vel;
switch type
    case 1
        filt.vel.irlssmooth(traces_vel) = round(length(Data.t_snip)*0.16);
    case 2
        filt.vel.irlssmooth(traces_vel) = 200;
    case 3
        filt.vel.sgolay(traces_vel) = 5;
        filt.vel.irlssmooth(traces_vel) = 3;
end
end