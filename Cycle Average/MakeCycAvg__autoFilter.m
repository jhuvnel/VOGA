function filt = MakeCycAvg__autoFilter(Data,filt,plot_info)
% Just for debugging
%plot_info.YLim.Pos = 30*[-1,1];
%plot_info.YLim.Vel = 250*[-1,1];
% CycAvg = MakeCycAvg__filterTraces(Data,filt);
% ha = MakeCycAvg__plotFullCycAvg([],CycAvg,plot_info);
filt.pos{:,:} = NaN;
filt.vel{:,:} = NaN;
type = Data.info.type;
%% Filter Velocity
traces_vel = plot_info.traces_vel;
if type==1
%     CycAvg = MakeCycAvg__filterTraces(Data,filt);    
%     med_tr = NaN(length(traces_vel),length(CycAvg.t));
%     med_vel_tr = NaN(length(traces_vel),length(Data.ts));
%     keep_inds = CycAvg.Data_allcyc.keep_inds';
%     irls_p = [0,1,[1,2,5]*round(length(Data.t_snip)*0.16)]; 
%     err_val_irls = NaN(length(traces_vel),length(irls_p));
%     spline_p = [1-5*10.^(-5:-1:-10),1];
%     err_val_spline = NaN(length(traces_vel),length(spline_p));
%     for i = 1:length(traces_vel)
%         trac = traces_vel{i};
%         vel_tr = CycAvg.Data_filtvel.([trac(1),'E_Vel_',trac(2:end)]);
%         med_tr(i,:) = median(CycAvg.([lower(trac(1:2)),'_cyc']));   
%         for j = 1:size(keep_inds,1)
%             med_vel_tr(i,keep_inds(j,:)) = med_tr(i,:);
%         end
%         base_err = mean((vel_tr-med_vel_tr(i,:)).^2,'omitnan');        
%         for j = 1:length(irls_p)
%             filt_tr = reshape(filterTrace('irlssmooth',vel_tr,irls_p(j)),[],1);
%             med_cyc_tr = median(filt_tr(keep_inds));
%             err_val_irls(i,j) = mean(abs(filt_tr-med_vel_tr(i,:)),'omitnan')/base_err-...
%                 max(abs(med_cyc_tr))/max(abs(med_tr(i,:)));
%         end            
%         [~,ind] = min(err_val_irls(i,:));
%         filt.vel.irlssmooth(trac) = irls_p(ind);
%         vel_tr = filterTrace('irlssmooth',vel_tr,irls_p(ind));
%         for j = 1:length(spline_p)
%             filt_tr2 = filterTrace('spline',vel_tr,spline_p(j),Data.te,Data.te);
%             err_val_spline(i,j) = max(filt_tr2)*min(filt_tr2)/max(med_tr(i,:))/min(med_tr(i,:));
%         end   
%         %Spline
%             filt.vel.spline(trac) = spline_p(find(err_val_spline(i,:)>0.9,1,'first'));                
%     end 
    filt.vel.irlssmooth(traces_vel) = round(length(Data.t_snip)*0.16);
elseif type==2
    filt.vel.irlssmooth(traces_vel) = 200;
elseif type==3
    filt.vel.sgolay(traces_vel) = 5;
    filt.vel.irlssmooth(traces_vel) = 3;
end
end