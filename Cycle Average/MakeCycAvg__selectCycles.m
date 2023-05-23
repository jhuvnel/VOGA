function [CycAvg,filt] = MakeCycAvg__selectCycles(ha,CycAvg,plot_info,in_opt)
screen_size = plot_info.screen_size;
traces_vel = plot_info.traces_vel;
Data_cyc = CycAvg.Data_allcyc;
if CycAvg.Data.info.type == 2
    filt = CycAvg.filt;
    return;
end
list = {'Automatic','Click','List','Restart'};
tf = 1;
if nargin <4 ||~ismember(in_opt,list)
    in_opt = '';    
    sel = 'none';
else
    sel = in_opt;
end
while tf
    keep_tr = CycAvg.keep_tr;
    switch sel
        case 'Automatic'
            nc = length(keep_tr);
            %Extract eye and stim data
            fields = fieldnames(Data_cyc);
            traces_vel1 = traces_vel;
            for i = 1:length(traces_vel)
                traces_vel1{i} = [traces_vel1{i}(1),'E_Vel_',traces_vel{i}(2:end)];
            end
            eye_fields = fields(contains(fields,traces_vel1)&~contains(fields,{'saccade','smooth'}));
            all_trac = cell(length(eye_fields),1);
            for i = 1:length(eye_fields)
                all_trac{i} = Data_cyc.(eye_fields{i});
            end
            trac_cyc = vertcat(all_trac{:});
            ind = [];
            val = Inf;
            while val > 1
                keep_tr(ind) = false;
                trac_cyc(:,~keep_tr) = NaN;
                sd_comp = std(trac_cyc,[],2,'omitnan');
                near_dist = NaN*trac_cyc;
                for i = 1:nc
                    near_dist(:,i) = sqrt(sum((trac_cyc-trac_cyc(:,i)).^2,2,'omitnan')/size(trac_cyc,1));
                end
                z_dist = (near_dist-mean(near_dist,2,'omitnan'))./std(near_dist,[],2,'omitnan');
                [~,ind] = max(max(z_dist,[],'omitnan'),[],'omitnan');
                sd_new = std(trac_cyc(:,~ismember(1:nc,ind)),[],2,'omitnan');
                perc_sd = mean((sd_comp-sd_new)./sd_new);
                val = 2*sum(keep_tr)*perc_sd;
            end
            %Recalculate and plot
            CycAvg.keep_tr = keep_tr;
            CycAvg = MakeCycAvg__filterTraces([],[],CycAvg);
            MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
        case 'Click'
            keep_tr1 = keep_tr;
            %Extract eye and stim data
            fields = fieldnames(Data_cyc);
            traces_vel1 = traces_vel;
            for i = 1:length(traces_vel)
                traces_vel1{i} = [traces_vel1{i}(1),'E_Vel_',traces_vel{i}(2:end)];
            end
            eye_fields = fields(contains(fields,traces_vel1)&~contains(fields,{'saccade','smooth'}));
            all_trac = NaN(length(Data_cyc.t),length(keep_tr),length(eye_fields)+1);
            for i = 1:length(eye_fields)
                all_trac(:,:,i) = Data_cyc.(eye_fields{i});
            end
            if all(size(Data_cyc.stim) > 1)
                all_trac(:,:,end) = Data_cyc.stim;
            end
            uiwait(msgbox(['Click near an erroneous trace on the aligned cycle graph.',newline,...
                'Use the up and down arrows to pan through traces.',newline,...
                'Press Enter to delete a trace and ESC to finish the process.']))
            [x,y] = ginput(1); %Assume this is on a cycle graph
            t_ind = find(Data_cyc.t>x,1,'first');
            traces = reshape(all_trac(t_ind,:,:),length(keep_tr),[])';
            [~,trace_i] = sort(min(abs(traces-y)));
            trace_i(ismember(trace_i,find(~keep_tr))) = [];
            ind = 1;
            axes(ha(4))
            hold on
            h = plot(Data_cyc.t,reshape(all_trac(:,trace_i(ind),:),length(Data_cyc.t),[]),'c','LineWidth',2);
            clc;
            k = waitforbuttonpress;
            value = double(get(gcf,'CurrentCharacter'));
            while value ~= 27  %Enter key=13 or ESC = 27, run until time to stop
                if value == 13 %Remove trace
                    hold off
                    keep_tr(trace_i(ind)) = false;
                    CycAvg.keep_tr = keep_tr;
                    CycAvg = MakeCycAvg__filterTraces([],[],CycAvg);
                    MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
                    axes(ha(4))
                    hold on
                    trace_i(ind) = [];
                    if ind>length(trace_i)
                        ind = ind-1;
                    end
                    h = plot(Data_cyc.t,reshape(all_trac(:,trace_i(ind),:),length(Data_cyc.t),[]),'c','LineWidth',2);
                elseif value == 30 %Up key: Trace closer to button press
                    if ind>1
                        delete(h)
                        ind = ind-1;
                        h = plot(Data_cyc.t,reshape(all_trac(:,trace_i(ind),:),length(Data_cyc.t),[]),'c','LineWidth',2);
                    end
                elseif value == 31 %Down key: Trace farther away from button press
                    if ind<length(trace_i)
                        delete(h)
                        ind = ind+1;
                        h = plot(Data_cyc.t,reshape(all_trac(:,trace_i(ind),:),length(Data_cyc.t),[]),'c','LineWidth',2);
                    end
                end
                k = waitforbuttonpress;
                value = double(get(gcf,'CurrentCharacter'));
            end
            delete(h)
            hold off
            disp('Cycles Manually Removed:')
            disp(find(keep_tr1&~keep_tr))
        case 'List'
            cyc_num = 1:length(keep_tr);
            [ind2,tf] = nmlistdlg('PromptString','Select cycles:',...
                'SelectionMode','multiple','InitialValue',cyc_num(keep_tr),...
                'ListSize',[100 250],'ListString',cellstr(num2str(cyc_num')),...
                'Position',[screen_size(3)-4,screen_size(4)-5.25,2,5.25]);
            if ~tf
                filt = CycAvg.filt;
                return;
            end
            keep_tr = false(1,length(keep_tr));
            keep_tr(ind2) = true;
            %Recalculate and plot
            CycAvg.keep_tr = keep_tr;
            CycAvg = MakeCycAvg__filterTraces([],[],CycAvg);
            MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
        case 'Restart'
            keep_tr = true(1,length(keep_tr));
            %Recalculate and plot
            CycAvg.keep_tr = keep_tr;
            CycAvg = MakeCycAvg__filterTraces([],[],CycAvg);
            MakeCycAvg__plotFullCycAvg(ha,CycAvg,plot_info);
    end
    if ~isempty(in_opt)
        break; %Just run once
    end
    %Poll user again
    [ind_l,tf] = nmlistdlg('PromptString','Type of cycles:',...
        'SelectionMode','single','ListSize',[100 80],'ListString',list,...
        'Position',[screen_size(3)-4,screen_size(4)-3,2,3],...
        'CancelString','Done','OKString','Select');
    if tf
        sel = list{ind_l};
    end
end
filt = CycAvg.filt;
end