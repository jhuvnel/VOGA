%% Manually Segment
%Supplement to the segmenting pipeline
function ManuallySegment(In_Path,Seg_Path)
% Take the Raw VOG files and convert them to a Data struct
[DATA,stim_info] = RawVOGtoSegment(In_Path);
%Extract signals needed for segmenting
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
YLim = get(gca,'YLim');
%% Manually Select Segments
starts = NaN(1,length(stim_info));
ends = NaN(1,length(stim_info));
for i = 1:length(starts)
    vals = ginput(2);
    [~,temp] = min(abs(Time_Eye-vals(1,1)));
    starts(i) = round(temp,0);
    [~,temp] = min(abs(Time_Eye-vals(2,1)));
    ends(i) = round(temp,0);        
    hold on
    fill(t([starts(i);ends(i);ends(i);starts(i)]),YLim([2;2;1;1]),'g');
    plot(t,Gyro(:,1),'k:',t,Gyro(:,2),'k--',t,Gyro(:,3),'k',t,10*Stim,'b')
    hold off        
end
pause(1)
%% Save
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
%Define the fields in order that should only be extracted for the relevant time of start:end
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
    %Save but make a new ending if there are multiple segments with the same information
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