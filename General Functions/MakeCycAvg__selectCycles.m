function [keep_tr,tf] = MakeCycAvg__selectCycles(type,keep_tr,Data_cyc,screen_size)    
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
        [x,y] = ginput(1); %Assume this is on a cycle graph
        t_ind = find(Data_cyc.t>x,1,'first');
        %Extract eye data
        fields = fieldnames(Data_cyc);
        eye_fields = fields(contains(fields,'Vel'));
        eyes = NaN(length(eye_fields),length(keep_tr));
        for i = 1:length(eye_fields)
            eyes(i,:) = Data_cyc.(eye_fields{i})(t_ind,:); 
        end
        if all(size(Data_cyc.stim)>1) %multiple head traces too
            traces = [Data_cyc.stim(t_ind,:);eyes];
        else
            traces = eyes;
        end
        traces(:,~keep_tr) = NaN;
        [~,trace_i] = min(min(abs(traces-y)));
        keep_tr(trace_i) = false;
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