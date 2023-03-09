function filt = MakeCycAvg__autoFilter(Data,filt)
type = Data.info.type;
if type == 1
    filt.vel.irlssmooth(end) = round(length(Data.t_snip)*0.16); %heuristic
elseif type == 2
    filt.vel.irlssmooth(end) = 400; %heuristic
end
end