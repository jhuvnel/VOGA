function [keep_tr,ha,tf] = MakeCycAvg__selectCycles(ha,type,keep_tr,Data_cyc,screen_size)    
    if type == 2
        disp('No cycle selection for this data type')
        return;
    end
    [ind,tf] = nmlistdlg('PromptString','Type of cycles:',...
                       'SelectionMode','single',...
                       'ListSize',[100 80],...
                       'ListString',{'Click','List'},...
                       'Position',[screen_size(3)-4,screen_size(4)-3,2,3]); 
    if ~tf
        return;
    end
    if ind == 1        
        %Extract eye and stim data
        fields = fieldnames(Data_cyc);
        eye_fields = fields(contains(fields,'Vel'));
        all_trac = NaN(length(Data_cyc.t),length(keep_tr),length(eye_fields)+1);
        for i = 1:length(eye_fields)
            all_trac(:,:,i) = Data_cyc.(eye_fields{i}); 
        end
        all_trac(:,:,end) = Data_cyc.stim;
        
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
    elseif ind ==2
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