function Data = MakeCycAvg__alignCycles(Data)
fname = Data.info.name;
info = Data.info;
Fs = Data.Fs;
if contains(fname,'Activation') %Preserve time to rejoin them later
    Data.te = Data.Time_Eye;
    Data.ts = Data.Time_Stim;
else %Make 1st point 0
    Data.te = Data.Time_Eye - Data.Time_Eye(1);
    Data.ts = Data.Time_Stim - Data.Time_Stim(1);
end
if isfield(Data,'HeadMPUVel_Z') %True for some rotary chair files
    Data.HeadVel_Z = Data.HeadMPUVel_Z;
end
if ~isfield(Data,'HeadVel_L')&&isfield(Data,'HeadVel_X')&&isfield(Data,'HeadVel_Y')
    Data.HeadVel_L = (Data.HeadVel_X - Data.HeadVel_Y)/sqrt(2);
end
if ~isfield(Data,'HeadVel_R')&&isfield(Data,'HeadVel_X')&&isfield(Data,'HeadVel_Y')
    Data.HeadVel_R = (Data.HeadVel_X + Data.HeadVel_Y)/sqrt(2);
end
if contains(fname,'eeVOR')&&contains(fname,'Step') %Fix the trigger
    %Fix the stimulus trace
    temp_ends = find(diff(Data.Trigger)==-1)-1; %The last index of this is the change to the break.
    Data.stim1 = 0*Data.Trigger;
    Data.stim1(1:(temp_ends(end)-1)) = 1;
elseif contains(fname,{'eeVOR','Moogles'})
    Data.stim1 = Data.Trigger;
else
    dir = ['L','R','Z'];
    Data.stim1 = Data.(['HeadVel_',dir(cellfun(@(x) contains(fname,x),{{'LA','RP'};{'RA','LP'};{'LH','RH','Rotary'}}))]);    
end
%Shift Trigger if needed
ts = Data.ts;
stim1 = Data.stim1;
stim = reshape(stim1,[],1);
len = length(stim);
TrigShift = round(info.TriggerShift2);
if TrigShift > 0
    stim = [repmat(stim(1),TrigShift,1);stim];
    stim = stim(1:len);
elseif TrigShift < 0
    stim = [stim;repmat(stim(end),-TrigShift,1)];
    stim = stim((-TrigShift+1):end);
end
if contains(info.dataType,'Impulse')
    thresh = 50;
    is_pos = (contains(fname,'GNO')&&contains(fname,{'LH','RP','RA'}))||...
        (~contains(fname,'GNO')&&contains(fname,{'LH','RP','LP'}));
    if is_pos
        abov_i = find(stim>thresh);
    else
        abov_i = find(stim<-thresh);
    end
    p_len = 0;
    spike_i = abov_i;
    while p_len~=length(spike_i) %Just to make sure it's really done running
        p_len = length(spike_i);
        for i = 1:p_len
            snip = spike_i(i) + (-floor(0.5*Fs):floor(0.5*Fs)); %1 second apart at least
            snip(snip < 1|snip > length(stim)) = [];
            [~,max_i] = max(abs(stim(snip)));
            spike_i(i) = snip(1)-1+max_i;
        end
        spike_i = unique(spike_i);
    end
    if is_pos
        spike_i(stim(spike_i)<0)=[];
    else
        spike_i(stim(spike_i)>0)=[];
    end
    %Consisitent with GNO's csv, take a 175 sample trace with max at
    %sample 48
    starts = spike_i-floor(0.195*Fs);
    ends = starts+floor(0.71*Fs);
    inv_i = starts<1|ends>length(stim)|(stim(spike_i)-stim(spike_i-floor(0.1*Fs))<50);
    starts(inv_i)= [];
    ends(inv_i) = [];
    t_snip = reshape(median(diff(ts))*(0:ends(1)-starts(1)),1,[]);
    keep_inds = zeros(ends(1)-starts(1)+1,length(starts));
    for i = 1:length(starts)
        keep_inds(:,i) = starts(i):ends(i);
    end
    stims = stim(keep_inds);
elseif contains(info.dataType,{'RotaryChair','aHIT'})||contains(info.goggle_ver,'Moogles') %Align based on real/virtual motion traces
    if contains(info.dataType,'Sine')
        fparts = split(info.dataType,'-');
        freqs = fparts(contains(fparts,'Hz'));
        freq = zeros(1,length(freqs));
        for i = 1:length(freqs)
            freq(i) = str2double(strrep(freqs(i),'Hz',''));
        end
        amp = str2double(strrep(fparts{contains(fparts,'dps')},'dps',''));
        snip_len = floor(Fs/min(freq));
        template = zeros(length(freq),snip_len);
        for i = 1:length(freq)
            template(i,:) = sin(2*pi*freq(i)*ts(1:snip_len));
        end
        template = sum(template,1);
        template = amp/max(template)*template;
        if size(stim,2)==1
            template = template';
        end
        %Find the mismatch between signal and template
        errors = NaN(1,1000);
        sub_i = floor(linspace(1,length(stim)-snip_len,1000));
        for i = 1:1000
            errors(i) = sum((template-stim(sub_i(i):sub_i(i)+snip_len-1)).^2);
        end
        errors = (errors - min(errors))/(max(errors)-min(errors));
        errors2 = errors;
        errors2(errors>0.1)=NaN;
        poss_val = find(~isnan(errors2));
        region_s = [poss_val(1),poss_val([false,diff(poss_val)>1])];
        region_e = [poss_val([diff(poss_val)>1,false]),poss_val(end)];
        starts = region_s;
        for i = 1:length(starts)
            temp_err = NaN*errors;
            temp_err(region_s(i):region_e(i)) = errors(region_s(i):region_e(i));
            [~,m_ind] = min(temp_err);
            starts(i) = sub_i(m_ind);
        end
        %Now find zeros crossings on the stimulus trace itself
        neg_stim = find(stim < 0);
        pos_stim = find(stim > 0);
        for i = 1:length(starts)
            if stim(starts(i)) < 0 %look for next positive value
                l_i = pos_stim(find(pos_stim>starts(i),1,'first'));
                if isempty(l_i)
                    l_i = length(starts);
                end
                poss_i = [l_i-1 l_i];
            else %look for previous negative value
                l_i = neg_stim(find(neg_stim<starts(i),1,'last'));
                if isempty(l_i)
                    l_i = 1;
                end
                poss_i = [l_i l_i+1];
            end
            [~,ind] = min(abs(stim(poss_i)));
            starts(i) = poss_i(ind);
        end
        starts = unique(starts); %remove duplicates
        if length(starts) > 1
            snip_len1 = round(median(diff(starts)));
            snip_len2 = length(stim) - starts(end);
            if abs(snip_len2-snip_len1)/snip_len1 < 0.01 %Less than 1% off of the expected cycle length
                snip_len = min([snip_len1 snip_len2]);
            else
                snip_len = snip_len1;
            end
        end
        ends = starts + snip_len - 1;
        %Delete incomplete cycles
        starts(ends>length(stim)) = [];
        ends(ends>length(stim)) = [];
        all_stim = zeros(snip_len,length(starts));
        for i = 1:length(starts)
            all_stim(:,i) = stim(starts(i):ends(i));
        end
        %Remove any obvious erroneous motion traces
        if contains(info.dataType,{'RotaryChair'})
            tol = 0.2; %Amplitude can be 20% wrong and still be tolerated
            rm_tr = abs(max(all_stim)-amp)/amp > tol | abs(min(all_stim)+amp)/amp > tol;
            starts(rm_tr) = [];
            ends(rm_tr) = [];
            all_stim(:,rm_tr) = [];
        end
        stims = mean(all_stim,2);
    elseif contains(info.dataType,'Step')
        stims = stim;
        starts = 1;
        ends = length(ts);
    else
        error('Unknown Data Type (RotaryChair/aHIT)')
    end
    t_snip = reshape(ts(1:size(stims,1))-ts(1),1,[]);
    keep_inds = zeros(ends(1)-starts(1)+1,length(starts));
    for i = 1:length(starts)
        keep_inds(:,i) = starts(i):ends(i);
    end
elseif contains(info.dataType,'eeVOR') %align using the trigger signal
    if contains(info.dataType,{'65Vector','MultiVector'})
        %The trigger is actually showing when the trapezoids start and end. There
        %are only 20 cycles of the stimulus applied and there are 40 trigger
        %toggles. To make the stimulus trace, I assumed the trigger was high when
        %the virtual velocity was linear and that the max virtual velocity is 50
        %dps (estimated from PJB's 2019 manuscript).
        %The above was confirmed by the the VOG code.
        %I also want both excitation and inhibition phases of the stimulus during
        %alignemnt.
        %Find time window for alignment
        trig = diff(stim);
        starts = find(trig==1)-1;
        starts = starts(1:2:end);
        snip_len = floor(median(diff(starts)));
        ends = starts + snip_len;
        %Create model stimulus trace
        stims = stim(starts(1):ends(1));
        ind1 = find(diff(stims)==1);
        ind2 = find(diff(stims)==-1);
        stims(ind1(1):ind2(1)) = linspace(0,50,length(ind1(1):ind2(1)));
        stims(ind2(1)+1:ind1(2)-1) = 50;
        stims(ind1(2):ind2(2)) = linspace(50,0,length(ind1(2):ind2(2)));
        stims(ind2(2)+1:end) = 0;
    elseif contains(info.dataType,'Sine') %sine (toggle = new cycle), remove last cycle
        trig = abs(diff(stim));
        starts = find(trig==1);
        snip_len = round(median(diff(starts)),0);
        starts(end) = []; %last cycle is the change to rest
        ends = starts + snip_len - 1;
        fparts = split(info.dataType,'-');
        freq = str2double(strrep(fparts{contains(fparts,'Hz')},'Hz',''));
        amp = str2double(strrep(fparts{contains(fparts,'dps')},'dps',''));
        stims = amp*sin(2*pi*freq*(ts(1:snip_len)));
    elseif contains(info.dataType,{'PulseTrain','Autoscan'}) % (high = on, low = off)
        trig = (diff(stim));
        starts = find(trig==1);
        snip_len = round(median(diff(starts)));
        ends = starts + snip_len - 1;
        stims = 0*stim(starts(1):ends(1));
        stims(1+(find(trig(starts(1):ends(1))==1,1,'first'):find(trig(starts(1):ends(1))==-1,1,'first'))) = 50;
    elseif contains(info.dataType,'Activation')
        % for activation data low = dark, high = light but don't cycle average here
        stims = 50*stim;
        starts = 1;
        ends = length(ts);
    elseif contains(info.dataType,'Step')
        fparts = split(info.dataType,'-');
        amp = str2double(strrep(strrep(fparts{contains(fparts,'dps')},'dps',''),'n','-'));
        stims = amp*stim;
        stim = amp*stim;
        starts = 1;
        ends = length(ts);
    else
        error('Unknown Data Type (eeVOR)')
    end
    if ends(end) > length(ts)
        starts(end) = [];
        ends(end) = [];
    end
    t_snip = reshape(ts(1:length(stims))-ts(1),1,[]);
    keep_inds = zeros(ends(1)-starts(1)+1,length(starts));
    for i = 1:length(starts)
        keep_inds(:,i) = starts(i):ends(i);
    end
else
    error('Unknown Data Type')
end
Data.stim = stim;
Data.t_snip = t_snip;
Data.stims = stims;
Data.keep_inds = keep_inds;
%Find the detected impulses in the GNO system
if isfield(Data,'DetectedTraces_HeadVel')&&~isempty(Data.DetectedTraces_HeadVel)
    detec_head = Data.DetectedTraces_HeadVel;
    if max(mean(stims,2))>abs(min(mean(stims,2)))
        detec_head(:,max(detec_head)<abs(min(detec_head))) = [];
    else
        detec_head(:,max(detec_head)>abs(min(detec_head))) = [];
    end
    detec_tr = NaN(1,size(detec_head,2));
    for j = 1:size(detec_head,2)
        [~,detec_tr(:,j)] = min(sum(abs(stims-detec_head(:,j))));
    end
else
    detec_tr = []; %initialize, should only be for GNO traces
end
Data.detec_tr = detec_tr;
end