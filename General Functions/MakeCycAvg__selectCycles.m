function [keep_tr,tf] = MakeCycAvg__selectCycles(type,keep_tr,t_snip,stims,LE_V,RE_V)    
    if type == 2
        disp('No cycle selection for this data type')
        return;
    end
    [ind,tf] = nmlistdlg('PromptString','Type of cycles:',...
                       'SelectionMode','multiple',...
                       'ListSize',[100 80],...
                       'ListString',{'Click','List'},...
                       'Position',[11,8,2,3]); 
    if ~tf
        return;
    end
    if ind == 1
        [x,y] = ginput(1); %Assume this is on a cycle graph
        t_ind = find(t_snip>x,1,'first');
        if type == 1
            eyes = [LE_V.LHRH(t_ind,:);LE_V.LARP(t_ind,:);LE_V.RALP(t_ind,:);...
                RE_V.LHRH(t_ind,:);RE_V.LARP(t_ind,:);RE_V.RALP(t_ind,:)];
            eyes(:,~keep_tr) = NaN;
            [~,eye_i] = min(min(abs(eyes-y)));
            keep_tr(eye_i) = false;
        elseif type == 3
            stims(t_ind,~keep_tr) = NaN;
            RE_V(t_ind,~keep_tr) = NaN;
            [head_v,head_i] = sort(abs(stims(t_ind,:)-y));
            [eye_v,eye_i] = sort(abs(RE_V(t_ind,:)-y));
            if head_v(1) < eye_v(1)
                keep_tr(head_i(1)) = false;
            else
                keep_tr(eye_i(1)) = false;
            end
        end
    elseif ind ==2
        cyc_num = 1:length(keep_tr);
        [ind2,tf] = nmlistdlg('PromptString','Select cycles:',...
                           'SelectionMode','multiple',...
                           'InitialValue',cyc_num(keep_tr),...
                           'ListSize',[100 250],...
                           'ListString',cellstr(num2str(cyc_num')),...
                           'Position',[11,5.25,2,5.25]);           
        if ~tf
            return;
        end
        keep_tr = false(1,length(keep_tr));
        keep_tr(ind2) = true; 
    end
end