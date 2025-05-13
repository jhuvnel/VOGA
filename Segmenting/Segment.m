%% Segment
%
% Make a segmenting pipeline that can take any Raw VOG file and turn it into
% a Data struct with the eye movement, head movement, and external trigger
% information for each experiment.
% Works for LDVOG, NKI, GNO, and ESC. Other goggle setups should be added 
% to RawVOGtoSegment.
% Works for Rotary Chair, aHIT, Manual (impulses) and eeVOR. Experiment
% types are Sine, VelStep, Impulse, Gaussian, PulseTrain, Activation,
% Autoscan, and MultiVector. Add new experiment types to be segmented into
% this script.
% Any file that fails automatic segmenting gets sent to the ManuallySegment
% script instead.
%
function Segment(In_Path,Seg_Path)
load('VNELcolors.mat','colors')
% Take the Raw VOG files and convert them to a Data struct
[DATA,stim_info] = RawVOGtoSegment(In_Path);
%Extract signals needed for segmenting
Fs = DATA.Fs;
t = DATA.Time_Stim;
Stim = NaN*t;
if isfield(DATA,'Trigger')
    Stim = DATA.Trigger;
end
gyro_fields = {'HeadVel_L','HeadVel_R','HeadVel_Z'};
Gyro = NaN(length(t),length(gyro_fields));
for i = 1:length(gyro_fields)
    if isfield(DATA,gyro_fields{i})
        Gyro(:,i) = DATA.(gyro_fields{i});
    end
end
GyroAll = (2*double((sum(Gyro,2,'omitnan'))>0)-1).*sqrt(sum(Gyro.^2,2,'omitnan'));
% Initialize Figure
fig = figure(1); 
clf; %in case there are leftover anotations
%Figure out how big the screen is and then leave 6 inches of space on
%the right hand side
set(fig,'Color','w','Units','normalized','Position',[0,0,1,1]);
ax = subplot(1,1,1);
plot(t,Gyro(:,1),'k:',t,Gyro(:,2),'k--',t,Gyro(:,3),'k',t,10*Stim,'b')
xlabel(ax,'Time (s)'); ylabel(ax,'Velocity (dps)');
ax.Position = [0.05 0.1 0.7 0.83]; %Make space for the notes annotations
annotation('textbox',[0.8 0.1 0.25 0.83],'String',stim_info,...
    'FontSize',11,'HorizontalAlignment','left','EdgeColor','none');
XLim = get(gca,'XLim'); YLim = get(gca,'YLim');
% Find which experiments are present
exp_types = [extract(stim_info,"RotaryChair"|"aHIT"|"Manual"|"eeVOR"|"trash"),...
    extract(stim_info,"VelStep"|"Sine"|"Impulse"|"Gaussian"|"PulseTrain"|"MultiVector"|"Activation"|"Autoscan"|"trash")];
if any(cellfun(@isempty,exp_types(:,2)))
    error('Unexpected experiment type.')
end
exp_types2 = join(exp_types);
exp_types2(contains(exp_types2,'trash')) = [];
exp_uniq = unique(exp_types2);
starts = []; ends = []; %Initialize segment starts and end indices
%% Find the start and end of each segment based on experiment type
if contains(exp_uniq,'RotaryChair VelStep')
    %Expecting a constant velocity for ~60 seconds followed by no motion for ~60
    %seconds. Often LH followed by RH in the same file.
    head_still = find(abs(GyroAll)<5); %Threshold of no motion.
    head_moving = find(abs(GyroAll)>50); %Threshold of chair moving, change as needed
    head_starts = head_moving([true;diff(head_moving)>Fs*10]); %make sure the break between velstep segments is at least 10s long
    [~,inds] = min(abs(repmat(head_still',length(head_starts),1)-head_starts),[],2); %Find the first no motion point before the motions
    seg_start = head_still(inds);
    %Make sure that detected segments within 10 seconds of each other are
    %combined, which protects against ending a segment early because of a 
    %tracking drop or cough
    seg_start = seg_start([true;diff(seg_start)>Fs*10]); 
    seg_end = [seg_start(2:end)-1;length(t)]; 
    %Make sure any given segment is not longer than 2.25 min by cutting off
    %the segment if it is.
    seg_end((seg_end-seg_start)>round(135*Fs)) = seg_start((seg_end-seg_start)>round(135*Fs))+135*Fs; 
    starts = [starts;seg_start]; ends = [ends;seg_end];
end
if any(contains(exp_uniq,{'aHIT','Manual','RotaryChair Sine',...
        'RotaryChair Impulse','RotaryChair Gaussian'}))
    %Expecting a collection of head motions that cross 0 following by
    %periods of no head motion.
    zero_cross = [1;find(GyroAll(2:end).*GyroAll(1:end-1)<0);length(GyroAll)];
    head_moving = find(abs(GyroAll)>20); %change as needed but works for 0.01Hz - 2Hz
    starts1 = head_moving([true;diff(head_moving)>1]);
    ends1 = head_moving([diff(head_moving)>1;true]);
    for i = 1:length(starts1)
        starts1(i) = max([zero_cross(find(zero_cross<starts1(i),1,'last'))-10,1]);
        ends1(i) = min([zero_cross(find(zero_cross>ends1(i),1,'first'))+10,length(GyroAll)]);
    end
    %Combine segments with only 2 seconds between them
    for i = 1:length(starts1)-1
        if starts1(i+1)-ends1(i)<2*Fs
            starts1(i+1) = NaN; ends1(i) = NaN;
        end
    end
    seg_start = starts1(~isnan(starts1));
    seg_end = ends1(~isnan(ends1));
    %Remove segments that are less than 2 seconds long
    small_seg = seg_end-seg_start<2*Fs;
    seg_start(small_seg) = [];seg_end(small_seg) = [];
    starts = [starts;seg_start]; ends = [ends;seg_end];
end
if contains(exp_uniq,'Activation') % IN PROGRESS %
    %Need to split into light/dark chunks so it's possible to filter but
    %need the chunks to be contiguous so they can be stiched back together
    %In previous versions, there was a software trigger that could be used
    %to indicate Light/Dark cycles but this was not found for activations
    %after August 2018. 
    %This prompts the user to select the start and stop of the cycles by
    %hand
    t2 = t/60;
    cla;
    plot(t2,Gyro(:,1),'k:',t2,Gyro(:,2),'k--',t2,Gyro(:,3),'k',t2,10*Stim,'b')
    hold on
    plot(t2,DATA.LE_Position_X,'Color',colors.l_x)
    plot(t2,DATA.LE_Position_Y,'Color',colors.l_y)
    plot(t2,DATA.LE_Position_Z,'Color',colors.l_z)
    plot(t2,DATA.RE_Position_X,'Color',colors.r_x)
    plot(t2,DATA.RE_Position_Y,'Color',colors.r_y)
    plot(t2,DATA.RE_Position_Z,'Color',colors.r_z)
    hold off
    xlabel(ax,'Time (min)');
    set(gca,'YLim',[-30 30])
    starts = NaN(length(stim_info),1);
    ends = NaN(length(stim_info),1);
    disp('Make sure segments do not start or end on blinks.')
    for i = 1:length(stim_info)
        input(['Hit enter to select the starting and ending point of segment: ',stim_info{i}])
        vals = ginput(2);
        [~,temp] = min(abs(t2-vals(1,1)));
        [~,temp2] = min(abs(t2-vals(2,1)));
        starts(i) = round(temp,0); ends(i) = round(temp2,0);
    end  
    xlabel(ax,'Time (s)'); %Turn it back for the summary graph
end
if any(contains(exp_uniq,{'eeVOR PulseTrain','eeVOR Autoscan'}))
    %Pulse Train and Autoscan--high period is stimulation and low is break.
    temp = find(diff(Stim)==1)-1;
    temp2 = diff(temp);
    all_starts = temp;
    all_ends = all_starts+median(temp2);
    [~,ind] = sort(temp2,'descend');
    ind2 = sort(ind(1:sum(contains(stim_info,{'eeVOR PulseTrain','eeVOR Autoscan'}))-1));
    starts = [starts;all_starts([1;ind2+1])-round(0.5*median(diff(temp)))];
    ends = [ends;all_ends([ind2;end])+round(0.5*median(diff(temp)))];
end
if contains(exp_uniq,'eeVOR Sine')
    %Every trigger toggle is a cycle
    thresh_tol = Fs*0.4; %tolerance for differences in cyc length (always at least 400ms of break)
    temp = find(abs(diff(Stim))==1);
    temp2 = diff(diff(temp));
    all_ends = unique(temp([find(temp2<-thresh_tol);find(temp2>thresh_tol)+1;end]));
    all_starts = unique(temp([1;find(temp2<-thresh_tol)+1;find(temp2>thresh_tol)+2]));        
    starts = [starts;all_starts-5]; ends = [ends;all_ends+5];  
end
if contains(exp_uniq,{'eeVOR VelStep','eeVOR MultiVector'})
    %These files have trigger toggles for the "ramp" of every cycle
    if Stim(1)==1 %If Stim started "high", the first ramp is missing so add it back in
        temp = diff(find(abs(diff(Stim))==1));
        Stim(1:(find(diff(Stim)==-1,1,'first')-temp(2))) = 0;
    end
    temp = find(diff(Stim)==1)-1; %Start of the ramp
    temp2 = find(diff(Stim)==-1)+1; %End of the ramp
    dtemp = diff(temp); %Inter-ramp distance
    %Look at every other ramp because there is a blip for cycle starting and ending
    all_starts = temp(1:2:end); all_ends = temp2(2:2:end); 
    temp2 = [temp(3:2:end);length(Stim)]-all_starts;
    [~,ind] = sort(temp2,'descend');
    ind2 = sort(ind(1:sum(contains(stim_info,{'eeVOR VelStep','eeVOR MultiVector'})))); %Find the breaks that correspond to the number of notes
    cyc_per_seg = diff([0;ind2]);
    seg_starts = [all_starts([1;ind2(1:end-1)+1])-1;length(t)];
    seg_ends = all_ends(ind2);
    for i = 1:length(seg_ends) %Add padding to the last cycle to not cut off those data
        if cyc_per_seg(i) == 1 %Vel Step
            seg_ends(i) = min([seg_ends(i)+dtemp(2*ind2(i)-1)+1,seg_starts(i+1)]); %Make the cycle twice as long as the first part but not so long it clashes with the next one
        else
            seg_ends(i) = min([seg_ends(i)+dtemp(2*ind(i)-2)+1,seg_starts(i+1)]); 
        end
    end
    starts = [starts;seg_starts(1:end-1)]; ends = [ends;seg_ends]; 
end
%% Check to see if it can be saved
%Make sure none of the indices are out of bounds
starts(starts<0) = 1; ends(ends>length(t)) = length(t);
%Remove any duplicates at this stage
starts = unique(starts); ends = unique(ends);
%Sort the indices in case they are out of order;
[starts,sort_ind] = sort(starts); ends = ends(sort_ind);
%Make sure large index numbers are not stored in scientific notation
starts = round(starts,0); ends = round(ends,0);
%Remove any traces marked as "trash" or "calibrations"
rm = contains(lower(stim_info),{'trash','calib'});
starts(rm) = []; ends(rm) = []; stim_info(rm) = [];
if isempty(stim_info)
    return;
end
%If you can't automatically segment, do it manually
if length(starts)~=length(ends)||length(stim_info)>length(starts) 
    disp('There was an issue automatically segmenting.')
    ManuallySegment(In_Path,Seg_Path); 
    return;
elseif length(stim_info)<length(starts) %Extra segments detected
    %Let the user select which ones to use
    keep = false(1,length(starts));
    cla;
    hold on
    fills = gobjects(length(starts),1);
    for j = 1:length(starts) %Now plot all fills
        fills(j) = fill(t([starts(j);ends(j);ends(j);starts(j)]),YLim([2;2;1;1]),0.85*[1,1,1]);
    end
    plot(t,Gyro(:,1),'k:',t,Gyro(:,2),'k--',t,Gyro(:,3),'k',t,10*Stim,'b')
    hold off
    set(gca,'YLim',YLim)
    uiwait(msgbox('Click on all valid segments on the figure.'))
    for i = 1:length(stim_info)
        [x,~] = ginput(1);
        j = find(t(starts)-x<0,1,'last'); %Find closest segment to the click
        keep(j) = true;
        fills(j).FaceColor = 'g'; %Turn that segment green
    end
    starts = starts(keep); ends = ends(keep);
end
%Plot and show the segments
cla;
hold on
for j = 1:length(starts) %Now plot all fills
    fill(t([starts(j);ends(j);ends(j);starts(j)]),YLim([2;2;1;1]),'g');
end
plot(t,Gyro(:,1),'k:',t,Gyro(:,2),'k--',t,Gyro(:,3),'k',t,10*Stim,'b')
hold off
axis([XLim,YLim])
pause(1)
% Make changes to Impulse/Gaussian files
if any(contains(stim_info,{'Impulse','Gaussian'}))
    %Change the stim_info and detected segments to have an entry for each
    %canal. For these experiments, the notes are listed as LHRH but this
    %split that into two entries--one with LH and one with RH
    impgaus_inds = contains(stim_info,{'Impulse','Gaussian'})&contains(stim_info,{'LHRH','LARP','RALP'});
    rep_ind = sort([(1:length(starts))';find(impgaus_inds)]);
    starts = starts(rep_ind); ends = ends(rep_ind); 
    stim_info2 = repmat(stim_info,1,2);
    stim_info2(~impgaus_inds,2) = {''};
    stim_info2(impgaus_ind,:) = [strrep(strrep(strrep(stim_info2(impgaus_inds,1),'LHRH','LH'),'LARP','RP'),'RALP','RA'),...
                                 strrep(strrep(strrep(stim_info2(impgaus_inds,2),'LHRH','RH'),'LARP','LA'),'RALP','LP')];
    stim_info = reshape(stim_info2',[],1); stim_info(cellfun(@isempty,stim_info)) = [];                        
end 
%% Automatic segmenting process
%Define the fields in order that should only be extracted for the relevant
%time of start:end
seg_fields = {'Time_Eye','Time_Stim','Trigger',...
    'LE_Position_X','LE_Position_Y','LE_Position_Z',...
    'RE_Position_X','RE_Position_Y','RE_Position_Z',...
    'HeadVel_X','HeadVel_Y','HeadVel_Z','HeadVel_L','HeadVel_R',...
    'HeadAccel_X','HeadAccel_Y','HeadAccel_Z'};
%These fields should be included in every segment made from the file
seg_fields2 = {'DetectedTraces_HeadVel','AllData','CSVData','XMLData'};
info = DATA.info;
for i = 1:length(starts)
    i1 = starts(i); i2 = ends(i); %Indices of interest
    info.dataType = stim_info{i};
    %Set the stim_axis by the file name or experiment type
    info.stim_axis = double([contains(stim_info{i},{'LA','RP'}),contains(stim_info{i},{'RA','LP'}),contains(stim_info{i},{'LH','RH'})]);            
    if contains(stim_info{i},{'65Vector','MultiVector'}) %Stim vec in the file name
        info.stim_axis = str2double(info.dataType(strfind(info.dataType,'['):strfind(info.dataType,']')));
    elseif contains(stim_info{i},'X')
        info.stim_axis = [0.707,0.707,0];
    elseif contains(stim_info{i},'Y')
        info.stim_axis = [-0.707,0.707,0];
    elseif (contains(stim_info{i},{'Impulse','Gaussian'})&&contains(stim_info{i},'L'))||...
            contains(stim_info{i},{'PulseTrain','Autoscan'})&&strcmp(info.ear,'L')
        info.stim_axis = info.stim_axis.*[-1,-1,1]; %Left ear
    elseif contains(stim_info{i},{'Impulse','Gaussian','PulseTrain','Autoscan'})
         info.stim_axis = info.stim_axis.*[1,1,-1]; %Right ear        
    end
    fname = [info.subject,'-',info.visit,'-',info.exp_date,'-',info.goggle_ver,'-',info.dataType];
    info.name = [fname,'.mat'];
    Data.info = info;
    Data.Fs = DATA.Fs;
    Data.raw_start_t = t(i1);
    Data.raw_end_t = t(i2);
    Data.rawfile = {info.rawfile};
    for ii = 1:length(seg_fields)
        if isfield(DATA,seg_fields{ii})
            Data.(seg_fields{ii}) = DATA.(seg_fields{ii})(i1:i2);
        end
    end
    for ii = 1:length(seg_fields2)
        if isfield(DATA,seg_fields2{ii})
            Data.(seg_fields2{ii}) = DATA.(seg_fields{ii});
        end
    end
    save_flag = 1;
    %Save but make a new ending if there are multiple segments with
    %the same information
    if exist([Seg_Path,filesep,fname,'.mat'],'file') %Already an instance of this file
        Data2 = Data; %set the segment to Data 2 to check against current file
        load([Seg_Path,filesep,fname,'.mat'],'Data')
        if any(ismember(Data.rawfile,Data2.rawfile)&Data.raw_start_t==Data2.raw_start_t)
            % Same rawfile name and segment start time so ignore.
            disp([fname,' already exists in this folder and was ignored.'])
            save_flag = 0;
        else %Combine segments into one
            disp([fname,' already exists in this folder and they were combined.'])
            Data = CombineSegments(Data,Data2);
            delete([Seg_Path,filesep,fname,'.mat']) 
        end
    end
    if save_flag %Plot and save .mat and .fig of the segment
        save([Seg_Path,filesep,fname,'.mat'],'Data')
        fig2 = plotSegment(Data);
        savefig(fig2,[Seg_Path,filesep,fname,'.fig'])
        close(fig2);
    end
end
end