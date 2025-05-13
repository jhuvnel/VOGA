function Data = MakeCycAvg__alignCycles(Data)
fname = Data.info.name;
info = Data.info;
Fs = Data.Fs;
%Set time values
Data.te = Data.Time_Eye - Data.Time_Eye(1);
Data.ts = Data.Time_Stim - Data.Time_Stim(1);
if contains(fname,'Activation') %Preserve time to rejoin them later
    Data.te = Data.Time_Eye;
    Data.ts = Data.Time_Stim;    
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
if contains(fname,{'eeVOR','Moogles'})
    Data.stim1 = Data.Trigger;
else
    dir = ['L','R','Z'];
    Data.stim1 = Data.(['HeadVel_',dir(cellfun(@(x) contains(fname,x),...
        {{'LA','RP'};{'RA','LP'};{'LH','RH','Rotary'}}))]);    
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
% Set stim and find cycle start/end time indeces
fparts = split(info.dataType,'-');
if contains(info.dataType,'Impulse')
    amp = 50;
    if contains(info.dataType,'dps')
        amp = str2double(strrep(fparts{contains(fparts,'dps')},'dps',''));
    end
    is_pos = (2*double(contains(fname,'GNO')&&contains(fname,{'LH','RP','RA'})||...
        ~contains(fname,'GNO')&&contains(fname,{'LH','RP','LP'}))-1);
    [vals,spike_i,widths] = findpeaks(abs(stim),'MinPeakProminence',amp/2,...
        'MinPeakDistance',round(Fs)); %Impulses will be at least 1 second apart
    spike_i(isoutlier(vals)|isoutlier(widths)) = [];
    spike_i((2*is_pos-1)*stim(spike_i)<0)=[];
    %Consisitent with GNO's csv, take a 175 sample trace with max at sample 48
    starts = spike_i-floor(0.195*Fs);
    ends = starts+floor(0.71*Fs);       
elseif contains(info.dataType,{'RotaryChair','aHIT'})||contains(info.goggle_ver,'Moogles') %Align based on real/virtual motion traces
    if contains(info.dataType,'Sine')
        stim = medfilt1(stim,3);
        pos_stim = (abs(stim)+stim)/2; %Pos half-cycle only
        freq = str2double(strrep(fparts{contains(fparts,'Hz')},'Hz',''));
        amp = str2double(strrep(fparts{contains(fparts,'dps')},'dps',''));
        [vals,spike_i,widths] = findpeaks(pos_stim,'MinPeakProminence',amp/2,...
        'MinPeakDistance',round(0.3*Fs/freq),'Annotate','extents'); %Impulses will be at least 1 second apart
        spike_i(isoutlier(vals,"mean")|isoutlier(widths,"mean")) = [];
        starts = spike_i;
        for i = 1:length(spike_i)
            [min_val,ind] = sort(pos_stim(1:starts(i)));
            starts(i) = ind(find(min_val==min_val(1),1,'last'));            
        end
        snip_len = round(median(diff(starts)));
        ends = starts + snip_len - 1;
        if ends(end)>length(stim)&&(length(stim)-starts(end))/snip_len>0.99
            snip_len = length(stim)-starts(end);
            ends = starts + snip_len - 1;
        end
    elseif contains(info.dataType,'Step')
        stims = stim;
        starts = 1;
        ends = length(ts);
    else
        error('Unknown Data Type (RotaryChair/aHIT)')
    end
elseif contains(info.dataType,'eeVOR') %align using the trigger signal
    if contains(info.dataType,{'65Vector','MultiVector','Step'})
        %The trigger here shows when the stimulus ramps up and down.
        %50dps was chosen as a trigger value based on the Figures in the
        %Boutros 2019 JCI paper.
        %Find time window for alignment
        trig = diff(stim);
        starts = find(trig==1)-1;
        starts = starts(1:2:end);
        if isscalar(starts) %VelStep
            ends = length(ts);
            fparts = split(info.dataType,'-');
            amp = str2double(strrep(strrep(fparts{contains(fparts,'dps')},'dps',''),'n','-'));
            if contains(info.dataType,{'RH','LA','LP'})
                amp = -amp;
            end
        else %Multivector
            snip_len = floor(median(diff(starts)));
            ends = starts + snip_len;
            amp = 50;
        end
        %Create model stimulus trace
        stims = stim(starts(1):ends(1));
        ind1 = find(diff(stims)==1);
        ind2 = find(diff(stims)==-1);
        stims(ind1(1):ind2(1)) = linspace(0,amp,length(ind1(1):ind2(1)));
        stims(ind2(1)+1:ind1(2)-1) = amp;
        stims(ind1(2):ind2(2)) = linspace(amp,0,length(ind1(2):ind2(2)));
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
        temp = find(trig==-1,length(starts),'first');
        l2 = min([length(temp),length(starts)]);
        len2 = round(median(temp(1:l2)-starts(1:l2)));
        if len2 < 0 %cycle not counted
            len2 = round(median(temp(2:end)-starts(1:end-1)));
        end
        stims(2:len2) = 50;
    elseif contains(info.dataType,'Activation')
        % for activation data low = dark, high = light but don't cycle average here
        stims = 50*stim;
        starts = 1;
        ends = length(ts);
    else
        error('Unknown Data Type (eeVOR)')
    end
else
    error('Unknown Data Type')
end
%Remove starts/ends out of bounds
inv_i = starts<1|ends>length(stim);
starts(inv_i)= [];
ends(inv_i) = [];
%Remove repeated starts/ends
starts = unique(starts);
ends = unique(ends);
%Make the keep_inds matrix
keep_inds = zeros(ends(1)-starts(1)+1,length(starts));
for i = 1:length(starts)
    keep_inds(:,i) = starts(i):ends(i);
end
if contains(info.dataType,'Impulse')
    stims = stim(keep_inds);
elseif contains(info.dataType,'Sine')&&~contains(info.dataType,'eeVOR')
    stims = median(stim(keep_inds),2);
elseif contains(info.dataType,'eeVOR')&&contains(info.dataType,'Step')
    stim = stims;
end
Data.stim = stim;
Data.t_snip = reshape(median(diff(ts))*(0:ends(1)-starts(1)),1,[]);
if contains(fname,'Activation') %Preserve time to rejoin them later
    Data.t_snip = Data.te;
end
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