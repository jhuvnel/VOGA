function filt = MakeCycAvg__autoFilter(Data,filt)
type = Data.info.type;
if contains(Data.info.name,'Autoscan')
    %More filtering
    eye_trac = {'LX','RX','LY','RY','LZ','RZ'};
    all_trac = [Data.LE_Position_X,Data.RE_Position_X,Data.LE_Position_Y,...
        Data.RE_Position_Y,Data.LE_Position_Z,Data.RE_Position_Z];
    all_trac(isnan(all_trac)) = 0;       
    %Reset filters
    filt.pos{:,:} = NaN; 
    %Median filter is always good for dropouts--higher filter needed for
    %noisier data
    filt.pos.median(eye_trac) = 3;
    all_trac = medfilt1(all_trac,3);    
    %Spline items that are aready smooth
    filt.pos.spline(eye_trac(iqr(abs(diff(all_trac)))<0.5)) = 1-5e-5;
    filt.pos.spline(eye_trac(iqr(abs(diff(all_trac)))<0.2)) = 1-5e-6;
    filt.pos.spline(eye_trac(iqr(abs(diff(all_trac)))<0.1)) = 1-5e-7;
    filt.pos.spline(eye_trac(iqr(abs(diff(all_trac)))<0.05)) = 1-5e-8;
end
if type == 1
    filt.vel.irlssmooth(end) = round(length(Data.t_snip)*0.16); %heuristic
elseif type == 2
    filt.vel{:,:} = NaN;
    filt.vel.irlssmooth(end) = 200; %heuristic
end
end