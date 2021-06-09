function fig = plotSegment(Data) 
    % Standardize Colors
    load('VNELcolors.mat','colors')
    %% Extract and plot raw position data
    te = Data.Time_Eye - Data.Time_Eye(1);
    ts = Data.Time_Stim - Data.Time_Stim(1);
    info = Data.info;
    dType = strrep(strrep(info.dataType,'_',' '),'-',' ');
    if contains(dType,'[')&&contains(dType,']')
        dType(strfind(dType,'['):strfind(dType,']')) = strrep(dType(strfind(dType,'['):strfind(dType,']')),' ','-'); %Put negative signs back in vectors
    end
    %Trigger multiplier
    if ~contains(info.goggle_ver,'GNO')
        if contains(info.dataType,'RotaryChair')
            if isfield(Data,'HeadVel_Z')
                 stim = Data.HeadVel_Z;
            else
                stim = Data.HeadMPUVel_Z; 
            end   
            sm = 1;
        else
            stim = Data.Trigger; 
            sm = 20;
        end  
    end
    %Fix huge number of NaN values in torsion traces of NKI traces
    if contains(info.goggle_ver,'NKI')
        if sum(isnan(Data.LE_Position_X)) > 0.9*length(te) %less than 10% data integrity
            Data.LE_Position_X = zeros(length(te),1); %set to 0 so no torsion
        else
            Data.LE_Position_X = spline(te(~isnan(Data.LE_Position_X)),Data.LE_Position_X(~isnan(Data.LE_Position_X)),te);
        end    
        if sum(isnan(Data.RE_Position_X)) > 0.9*length(te)
            Data.RE_Position_X = zeros(length(te),1);
        else
            Data.RE_Position_X = spline(te(~isnan(Data.RE_Position_X)),Data.RE_Position_X(~isnan(Data.RE_Position_X)),te);
        end    
    end    
    fig = figure('Units','inches','Position',[0 0 7 10]);
    %Title
    annotation('textbox',[0 .9 1 .1],'String',dType,'FontSize',14,...
        'HorizontalAlignment','center','EdgeColor','none');
    subplot(1,1,1)
    if ~contains(info.goggle_ver,'GNO')
        plot(ts,sm*stim,'k')
        hold on
        plot(te,Data.LE_Position_X,'Color',colors.l_x)
        plot(te,Data.LE_Position_Y,'Color',colors.l_y)
        plot(te,Data.LE_Position_Z,'Color',colors.l_z)
        plot(te,Data.RE_Position_X,'Color',colors.r_x)
        plot(te,Data.RE_Position_Y,'Color',colors.r_y)
        plot(te,Data.RE_Position_Z,'Color',colors.r_z)
        hold off
        axis([0 te(end) -20 20])
        title('Raw Angular Position')
        xlabel('Time(s)')
        ylabel('Position (deg)')
        legend('Stim','L X','L Y','L Z','R X','R Y','R Z')
    else
        plot(ts,Data.HeadVel_L,'k:')
        hold on
        plot(ts,Data.HeadVel_R,'k--')
        plot(ts,Data.HeadVel_Z,'k-')
        plot(te,Data.RE_Vel_Y,'Color',colors.l_y)
        plot(te,Data.RE_Vel_Z,'Color',colors.l_z)
        hold off
        axis([0 te(end) -300 300])
        title('Angular Velocity')
        xlabel('Time(s)')
        ylabel('Velocity (dps)')
        legend('Head -LARP','Head RALP','Head LHRH','Inv Eye Y','Inv Eye Z')
    end
end