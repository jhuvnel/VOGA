function fig = plotSegment(Data) 
    % Standardize Colors
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
    %% Extract and plot raw position data
    te = Data.Time_Eye - Data.Time_Eye(1);
    ts = Data.Time_Stim - Data.Time_Stim(1);
    info = Data.info;
    dType = strrep(strrep(info.dataType,'_',' '),'-',' ');
    if contains(dType,'[')&&contains(dType,']')
        dType(strfind(dType,'['):strfind(dType,']')) = strrep(dType(strfind(dType,'['):strfind(dType,']')),' ','-'); %Put negative signs back in vectors
    end
    %Trigger multiplier
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

    fig = figure(1);
    delete(findall(gcf,'type','annotation')) %in case there are leftover anotations
    fig.Units = 'inches';
    fig.Position = [0 0 7 10];
    %Title
    annotation('textbox',[0 .9 1 .1],'String',dType,'FontSize',14,...
        'HorizontalAlignment','center','EdgeColor','none');
    subplot(1,1,1)
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
end