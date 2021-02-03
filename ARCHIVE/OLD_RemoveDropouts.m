%% RemoveDropouts

%Not currently used but may be used again
function tot_bad_ind = RemoveDropouts(Seg_Path,plot_on)
sep = filesep;
%% Standardize Colors
% Normal colors
colors.l_x = [237,150,33]/255;
colors.l_y = [125,46,143]/255;
colors.l_z = [1 0 0];
colors.l_l = [0,128,0]/255;
colors.l_r = [0 0 1];
colors.r_x = [237,204,33]/255;
colors.r_y = [125,46,230]/255;
colors.r_z = [1,0,1];
colors.r_l = [0 1 0];
colors.r_r = [64,224,208]/255;
% Faded colors
colors.l_x_s = colors.l_x + 0.5*(1-colors.l_x);
colors.l_y_s = colors.l_y + 0.5*(1-colors.l_y);
colors.l_z_s = colors.l_z + 0.5*(1-colors.l_z);
colors.l_l_s = colors.l_l + 0.5*(1-colors.l_l);
colors.l_r_s = colors.l_r + 0.5*(1-colors.l_r);
colors.r_x_s = colors.r_x + 0.5*(1-colors.r_x);
colors.r_y_s = colors.r_y + 0.5*(1-colors.r_y);
colors.r_z_s = colors.r_z + 0.5*(1-colors.r_z);
colors.r_l_s = colors.r_l + 0.5*(1-colors.r_l);
colors.r_r_s = colors.r_r + 0.5*(1-colors.r_r);
%% Find the scores for every file
files = dir([Seg_Path,sep,'*.mat']);
fnames = {files.name}';
perc_drop = zeros(length(fnames),6);
tot_bad_ind = zeros(length(fnames),1);
diff_thresh = 5;
for j = 1:length(fnames)
    load([Seg_Path,sep,fnames{j}],'Data')
    te = Data.Time_Eye;
    ts = Data.Time_Stim;
    stim = Data.Stim_Trig;
    if contains(Data.info.dataType,'RotaryChair')
        sm = 1;
    else
        sm = 20;
    end
    bad_ind = false(length(te),1);
    for i = 1:6
        %Load the raw trace to assess
        switch i
            case 1
                in_tr = Data.LE_Position_X;
            case 2
                in_tr = Data.LE_Position_Y;
            case 3
                in_tr = Data.LE_Position_Z;
            case 4
                in_tr = Data.RE_Position_X;
            case 5
                in_tr = Data.RE_Position_Y;
            case 6
                in_tr = Data.RE_Position_Z;
        end
        raw_tr = in_tr;
        %Set starting values to be overwritten
        prev_rm = 2;
        rm = 1;
        while ~(sum(rm)==0 || prev_rm==sum(rm))
            prev_rm = sum(rm);
            rm1 = isnan(raw_tr);
            stdev = std(raw_tr(~rm1));
            m = mean(raw_tr(~rm1));
            %Interpolate NaN values now    
            raw_tr = reshape(interp1(te(~rm1),raw_tr(~rm1),te),[],1);
            diff_raw_tr = reshape(abs(diff(raw_tr)),[],1);
            %Remove a trace if it's above a velocity threshold and it is 2 stdevs 
            %away from the mean in position or it is 4 stdevs away or is a NaN
            %value
            rm = rm1 | (abs(raw_tr-m)>2*stdev & [true;diff_raw_tr>diff_thresh]) | abs(raw_tr-m)>4*stdev; % 
            filt_tr = interp1(te(~rm),raw_tr(~rm),te);
            diff_filt_tr = abs(diff(filt_tr));
            perc_drop(j,i) = sum(diff_filt_tr>diff_thresh)/length(diff_filt_tr)*100;
            raw_tr = filt_tr;
        end
        %Update the indecies that still have problems
        bad_ind = bad_ind | [false;diff_filt_tr>diff_thresh];
        switch i
            case 1
                Data.LE_Position_X = filt_tr;
                Data.LE_Position_X_raw = in_tr;
            case 2
                Data.LE_Position_Y = filt_tr;
                Data.LE_Position_Y_raw = in_tr;
            case 3
                Data.LE_Position_Z = filt_tr;
                Data.LE_Position_Z_raw = in_tr;
            case 4
                Data.RE_Position_X = filt_tr;
                Data.RE_Position_X_raw = in_tr;
            case 5
                Data.RE_Position_Y = filt_tr;
                Data.RE_Position_Y_raw = in_tr;
            case 6
                Data.RE_Position_Z = filt_tr;
                Data.RE_Position_Z_raw = in_tr;
        end
    end
    tot_bad_ind(j,1) = 100*sum(bad_ind)/length(bad_ind);
    if plot_on
        h7 = plot(ts,sm*stim,'k');
        hold on
        plot(te,Data.LE_Position_X_raw,'Color',colors.l_x_s)
        plot(te,Data.LE_Position_Y_raw,'Color',colors.l_y_s)
        plot(te,Data.LE_Position_Z_raw,'Color',colors.l_z_s)
        plot(te,Data.RE_Position_X_raw,'Color',colors.r_x_s)
        plot(te,Data.RE_Position_Y_raw,'Color',colors.r_y_s)
        plot(te,Data.RE_Position_Z_raw,'Color',colors.r_z_s)
        h1 = plot(te,Data.LE_Position_X,'Color',colors.l_x);
        h2 = plot(te,Data.LE_Position_Y,'Color',colors.l_y);
        h3 = plot(te,Data.LE_Position_Z,'Color',colors.l_z);
        h4 = plot(te,Data.RE_Position_X,'Color',colors.r_x);
        h5 = plot(te,Data.RE_Position_Y,'Color',colors.r_y);
        h6 = plot(te,Data.RE_Position_Z,'Color',colors.r_z);
        hold off
        title(['Filtered Angular Position for ',strrep(fnames{j},'_','-')])
        xlabel('Time(s)')
        ylabel('Position (deg)')
        legend([h7,h1,h2,h3,h4,h5,h6],{'Trigger','L X','L Y','L Z','R X','R Y','R Z'})
        pause;
    end
    save([Seg_Path,sep,fnames{j}],'Data')
end
end