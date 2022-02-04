function [keep_tr,ha,tf] = MakeCycAvg__selectCycles(ha,type,keep_tr,Data_cyc,screen_size,traces_vel)    
    if type == 2
        disp('No cycle selection for this data type')
        return;
    end
    list = {'Automatic','Click','List'};
    [ind_l,tf] = nmlistdlg('PromptString','Type of cycles:',...
                       'SelectionMode','single',...
                       'ListSize',[100 80],...
                       'ListString',list,...
                       'Position',[screen_size(3)-4,screen_size(4)-3,2,3]); 
    if ~tf
        return;
    end
    if strcmp(list(ind_l),'Automatic')
        keep_tr1 = keep_tr;
        %Extract eye and stim data
        fields = fieldnames(Data_cyc);
        traces_vel1 = traces_vel;
        for i = 1:length(traces_vel)
            if traces_vel1{i}(1) == 'L'
                traces_vel1{i} = ['LE_Vel_',traces_vel{i}(2:end)];
            else
                traces_vel1{i} = ['RE_Vel_',traces_vel{i}(2:end)];
            end
        end
        eye_fields = fields(contains(fields,traces_vel1));
        for i = 1:length(eye_fields)
            templ = mean(Data_cyc.(eye_fields{i})(:,keep_tr),2);
            templ_sd = std(Data_cyc.(eye_fields{i})(:,keep_tr),[],2);            
            out_of_bounds = abs(Data_cyc.(eye_fields{i})-templ) > 2*templ_sd;
            rel_out = ~any(out_of_bounds(Data_cyc.t<0.23&Data_cyc.t>0.1,:))&any(Data_cyc.(eye_fields{i})(Data_cyc.t<0.2,:)<10)&any(Data_cyc.(eye_fields{i})(Data_cyc.t<0.4&Data_cyc.t>0.1,:)<0);
            keep_tr(~rel_out) = false;
        end  
        disp('Cycles Automatically Removed:')
        disp(find(keep_tr1&~keep_tr))
    elseif strcmp(list(ind_l),'Click')
        keep_tr1 = keep_tr;
        %Extract eye and stim data
        fields = fieldnames(Data_cyc);
        traces_vel1 = traces_vel;
        for i = 1:length(traces_vel)
            if traces_vel1{i}(1) == 'L'
                traces_vel1{i} = ['LE_Vel_',traces_vel{i}(2:end)];
            else
                traces_vel1{i} = ['RE_Vel_',traces_vel{i}(2:end)];
            end
        end
        eye_fields = fields(contains(fields,traces_vel1));
        all_trac = NaN(length(Data_cyc.t),length(keep_tr),length(eye_fields)+1);
        for i = 1:length(eye_fields)
            all_trac(:,:,i) = Data_cyc.(eye_fields{i}); 
        end
        if all(size(Data_cyc.stim) > 1)
            all_trac(:,:,end) = Data_cyc.stim;
        end
        [x,y] = ginput(1); %Assume this is on a cycle graph
        t_ind = find(Data_cyc.t>x,1,'first');
        traces = reshape(all_trac(t_ind,:,:),length(keep_tr),[])';
        traces(:,~keep_tr) = NaN;
        [~,trace_i] = sort(min(abs(traces-y)));
        ind = 1;
        if type == 1
            cyc_ax = ha(4);
        elseif type == 3
            cyc_ax = ha(3);
        end
        axes(cyc_ax)
        hold on
        h = plot(Data_cyc.t,reshape(all_trac(:,trace_i(ind),:),length(Data_cyc.t),[]),'c','LineWidth',2);
        clc;
        disp('Hit ESC to stop this process, Enter to delete a trace and use up/down arrows to navigate')
        k = waitforbuttonpress;
        value = double(get(gcf,'CurrentCharacter'));
        while value ~= 27  %Enter key=13 or ESC = 27, run until time to stop
            if value == 13    
                keep_tr(trace_i(ind)) = false;
            elseif value == 30 %Up key
                if ind>1
                    delete(h) 
                    ind = ind -1;
                    h = plot(Data_cyc.t,reshape(all_trac(:,trace_i(ind),:),length(Data_cyc.t),[]),'c','LineWidth',2);
                end
            elseif value == 31 %Down key
                if ind<length(keep_tr)
                    delete(h) 
                    ind = ind+1;
                    h = plot(Data_cyc.t,reshape(all_trac(:,trace_i(ind),:),length(Data_cyc.t),[]),'c','LineWidth',2);
                end
            else
                disp('Unrecognized key press')
            end
            k = waitforbuttonpress;
            value = double(get(gcf,'CurrentCharacter'));
        end
        delete(h)
        hold off
        disp('Cycles Manually Removed:')
        disp(find(keep_tr1&~keep_tr))
    elseif strcmp(list(ind_l),'List')
        cyc_num = 1:length(keep_tr);
        [ind2,tf] = nmlistdlg('PromptString','Select cycles:',...
                           'SelectionMode','multiple',...
                           'InitialValue',cyc_num(keep_tr),...
                           'ListSize',[100 250],...
                           'ListString',cellstr(num2str(cyc_num')),...
                           'Position',[screen_size(3)-4,screen_size(4)-5.25,2,5.25]);           
        if ~tf
            return;
        end
        keep_tr = false(1,length(keep_tr));
        keep_tr(ind2) = true; 
    end
end