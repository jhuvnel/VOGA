function filt = MakeCycAvg__autoFilter(Data,filt_params,plot_info)
type = Data.info.type;
gog = Data.info.goggle_ver;
traces_pos = plot_info.traces_pos;
traces_vel = plot_info.traces_vel;
filt = filt_params.filt; %last params for position
filt.vel{:,:} = NaN*filt.vel{:,:}; %Clear the velocity parameters
tr2cyc = @(tr) tr(Data.keep_inds);
T = size(Data.keep_inds,1);
switch type %Update velocity filtering
    case 1        
        if contains(gog,'LDVOG')
            filt.vel.irlssmooth(traces_vel) = round(T*0.16*1.5);
        else
            Data_vel = angpos2angvel(Data);            
            for t = 1:length(traces_vel)
                tr = strrep(strrep([traces_vel{t}(1),'E_Vel_',traces_vel{t}(2)],'_L','_LARP'),'_R','_RALP');
                full_tr = Data_vel.(tr);
                full_tr(isnan(full_tr)) = 0;
                med_tr = median(tr2cyc(full_tr),2);
                med_tr = mean([[NaN;med_tr(1:end-1)],med_tr,[med_tr(2:end);NaN]],2,'omitnan'); %apply a moving average filter with framelength 3
%                 pks= [findpeaks(med_tr,'MinPeakDistance',round(T/2),'NPeaks',1,'SortStr','descend');...
%                     findpeaks(-med_tr,'MinPeakDistance',round(T/2),'NPeaks',1,'SortStr','descend')];            
%                 if max(pks)>100
%                     filt.vel.irlssmooth(traces_vel{t}) = 1;
%                     filt.vel.spline(traces_vel{t}) = 0.9999995;
%                 elseif max(pks)>75
%                     filt.vel.irlssmooth(traces_vel{t}) = 1;
%                     filt.vel.spline(traces_vel{t}) = 0.999995;
%                 elseif max(pks)>10
%                     filt.vel.irlssmooth(traces_vel{t}) = round(T*0.16);
%                 else
%                     filt.vel.irlssmooth(traces_vel{t}) = round(T*0.24);
%                 end
                if contains(lower(Data.info.dataType),'10hz') || contains(lower(Data.info.dataType),'8hz')
                    filt.vel.irlssmooth(traces_vel{t}) = 0;
                    filt.vel.spline(traces_vel{t}) = 1;
                else
                    filt.vel.irlssmooth(traces_vel{t}) = 1;
                    filt.vel.spline(traces_vel{t}) = 1;
                end
                
            end
        end
    case 2
        filt.vel.irlssmooth(traces_vel) = 200;
    case 3
        filt.vel.sgolay(traces_vel) = 5;
        filt.vel.irlssmooth(traces_vel) = 3;
end
end