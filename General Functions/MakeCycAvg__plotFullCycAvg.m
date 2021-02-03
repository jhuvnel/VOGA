%% Make Full Plot Cyc Avg
function ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim_Pos,YLim_Vel,te,ts,t_snip,stim,stims,Data,Data_In,Data_cal,Data_calc,LE_V,RE_V,CycAvg,keep_inds,keep_tr)
    switch type
        case 1 %Cycle Averaging
            if isempty(ha) %first time running
                ha = gobjects(7,1);
                XLim_Long = [te(1) te(end)];
                XLim_Short = [t_snip(1) t_snip(end)];
                %For plots with cycles to select
                x1 = 0.06;
                x2 = 0.37;
                x3 = 0.68;
                y1 = 0.045;
                y2 = 0.38;
                y3 = 0.67; %0.68;
                wid_x_s = 0.30;
                wid_x_b = 0.92;
                height_y = 0.27;
                ha(1) = subplot(3,3,[1 2 3]);
                ha(2) = subplot(3,3,[4 5 6]);
                ha(3) = subplot(3,3,7);
                ha(4) = subplot(3,3,8);
                ha(5) = subplot(3,3,9);
                ha(1).Position = [x1 y3 wid_x_b height_y];
                ha(2).Position = [x1 y2 wid_x_b height_y];
                ha(3).Position = [x1 y1 wid_x_s height_y];
                ha(4).Position = [x2 y1 wid_x_s height_y];
                ha(5).Position = [x3 y1 wid_x_s height_y];
                %Link Axes
                linkaxes(ha([1,2]),'x')
                linkaxes(ha([2,3]),'y')
                linkaxes(ha([3,4,5]),'xy')
                %Make cycle labels
                cyc_num_labs = 1:length(keep_tr);
                x_tick = te(round(floor(mean(keep_inds)),0));
                %Set Labels
                %ha1
                set(ha(1),'XLim',XLim_Long,'YLim',YLim_Pos,'XTick',x_tick,'XTickLabel',cyc_num_labs,'Xaxislocation','top')
                ylabel(ha(1),'Angular Position (deg)','FontWeight','bold')
                %xlabel(ha(1),'Cycle Number')
                %title(ha(1),'Angular Position')
                %ha2
                set(ha(2),'XLim',XLim_Long,'YLim',YLim_Vel)
                ylabel(ha(2),'Angular Velocity (dps)','FontWeight','bold')
                xlabel(ha(2),'Time (s)')
                %title(ha(2),'Cycle Number')
                %ha3
                set(ha(3),'XLim',XLim_Short,'YLim',YLim_Vel)
                ylabel(ha(3),'Velocity (dps)')
                title(ha(3),'All Filtered Cycles')
                %ha4
                set(ha(4),'XLim',XLim_Short,'YLim',YLim_Vel,'YTickLabel',[])
                xlabel(ha(4),'Time (s)')
                %ha5
                set(ha(5),'XLim',XLim_Short,'YLim',YLim_Vel,'YTickLabel',[])
                title(ha(5),'Cycle Averages')
            end
            %Cycle Aligned Raw and Filtered Position
            axes(ha(1))
            hold on
            cla;
            for j = 1:size(keep_inds,2)
                if keep_tr(j)
                    fill([te(keep_inds(1,j)),te(keep_inds(end,j)),te(keep_inds(end,j)),te(keep_inds(1,j))]',[500,500,-500,-500]',colors.cyc_keep,'Tag',['Cycle_',num2str(j)]);
                else
                    fill([te(keep_inds(1,j)),te(keep_inds(end,j)),te(keep_inds(end,j)),te(keep_inds(1,j))]',[500,500,-500,-500]',colors.cyc_rm,'Tag',['Cycle_',num2str(j)]);
                end
            end
            %Plot Cycles
            h(1) = plot(ts,stim,'k','LineWidth',line_wid.norm);
            %Raw Data
            plot(ts,Data.LE_Position_X,'Color',colors.l_x_s,'LineWidth',line_wid.norm)
            plot(ts,Data.LE_Position_Y,'Color',colors.l_y_s,'LineWidth',line_wid.norm)
            plot(ts,Data.LE_Position_Z,'Color',colors.l_z_s,'LineWidth',line_wid.norm)
            plot(ts,Data.RE_Position_X,'Color',colors.r_x_s,'LineWidth',line_wid.norm)
            plot(ts,Data.RE_Position_Y,'Color',colors.r_y_s,'LineWidth',line_wid.norm)
            plot(ts,Data.RE_Position_Z,'Color',colors.r_z_s,'LineWidth',line_wid.norm)
            %Filtered Data
            h(2) = plot(ts,Data_In.Data_LE_Pos_X,'Color',colors.l_x,'LineWidth',line_wid.norm);
            h(3) = plot(ts,Data_In.Data_LE_Pos_Y,'Color',colors.l_y,'LineWidth',line_wid.norm);
            h(4) = plot(ts,Data_In.Data_LE_Pos_Z,'Color',colors.l_z,'LineWidth',line_wid.norm);
            h(5) = plot(ts,Data_In.Data_RE_Pos_X,'Color',colors.r_x,'LineWidth',line_wid.norm);
            h(6) = plot(ts,Data_In.Data_RE_Pos_Y,'Color',colors.r_y,'LineWidth',line_wid.norm);
            h(7) = plot(ts,Data_In.Data_RE_Pos_Z,'Color',colors.r_z,'LineWidth',line_wid.norm);
            hold off
            leg1 = legend(h,{'Stim','L X','L Y','L Z','R X','R Y','R Z'},'NumColumns',7);
            leg1.ItemTokenSize(1) = 7;
            %All Raw and Filtered Position
            axes(ha(2))
            hold on
            cla; 
            for j = 1:size(keep_inds,2)
                if keep_tr(j)
                    fill([te(keep_inds(1,j)),te(keep_inds(end,j)),te(keep_inds(end,j)),te(keep_inds(1,j))]',[500,500,-500,-500]',colors.cyc_keep,'Tag',['Cycle_',num2str(j)]);
                else
                    fill([te(keep_inds(1,j)),te(keep_inds(end,j)),te(keep_inds(end,j)),te(keep_inds(1,j))]',[500,500,-500,-500]',colors.cyc_rm,'Tag',['Cycle_',num2str(j)]);
                end
            end
            h1(1) = plot(ts,stim,'k','LineWidth',line_wid.norm);
            %Raw Data
            plot(ts,Data_cal.LE_Vel_LARP,'Color',colors.l_l_s,'LineWidth',line_wid.norm)
            plot(ts,Data_cal.LE_Vel_RALP,'Color',colors.l_r_s,'LineWidth',line_wid.norm)
            plot(ts,Data_cal.RE_Vel_LARP,'Color',colors.r_l_s,'LineWidth',line_wid.norm)
            plot(ts,Data_cal.RE_Vel_RALP,'Color',colors.r_r_s,'LineWidth',line_wid.norm)
            plot(ts,Data_cal.LE_Vel_Z,'Color',colors.l_z_s,'LineWidth',line_wid.norm)
            plot(ts,Data_cal.RE_Vel_Z,'Color',colors.r_z_s,'LineWidth',line_wid.norm)
            %Filtered Data
            h1(2) = plot(ts,Data_calc.LE_Vel_LARP,'Color',colors.l_l,'LineWidth',line_wid.norm);
            h1(3) = plot(ts,Data_calc.LE_Vel_RALP,'Color',colors.l_r,'LineWidth',line_wid.norm);
            h1(4) = plot(ts,Data_calc.LE_Vel_Z,'Color',colors.l_z,'LineWidth',line_wid.norm);
            h1(5) = plot(ts,Data_calc.RE_Vel_LARP,'Color',colors.r_l,'LineWidth',line_wid.norm);
            h1(6) = plot(ts,Data_calc.RE_Vel_RALP,'Color',colors.r_r,'LineWidth',line_wid.norm);
            h1(7) = plot(ts,Data_calc.RE_Vel_Z,'Color',colors.r_z,'LineWidth',line_wid.norm);
            hold off
            leg2 = legend(h1,{'Stim','L L','L R','L Z','R L','R R','R Z'},'NumColumns',7);
            leg2.ItemTokenSize(1) = 7;
            axes(ha(3))
            hold on
            cla;
            plot(t_snip,stims,'k','LineWidth',line_wid.norm)
            %Filtered Data Only
            for j = 1:size(keep_inds,2)
                plot(t_snip,LE_V.LARP(:,j),'Color',colors.l_l,'LineWidth',line_wid.norm,'Tag',['Cycle_',num2str(j)]);
                plot(t_snip,LE_V.RALP(:,j),'Color',colors.l_r,'LineWidth',line_wid.norm,'Tag',['Cycle_',num2str(j)]);
                plot(t_snip,RE_V.LARP(:,j),'Color',colors.r_l,'LineWidth',line_wid.norm,'Tag',['Cycle_',num2str(j)]);
                plot(t_snip,RE_V.RALP(:,j),'Color',colors.r_r,'LineWidth',line_wid.norm,'Tag',['Cycle_',num2str(j)]);
                plot(t_snip,LE_V.LHRH(:,j),'Color',colors.l_z,'LineWidth',line_wid.norm,'Tag',['Cycle_',num2str(j)]);
                plot(t_snip,RE_V.LHRH(:,j),'Color',colors.r_z,'LineWidth',line_wid.norm,'Tag',['Cycle_',num2str(j)]);
            end
            hold off
            axes(ha(4))
            hold on
            cla;
            plot(t_snip,stims,'k','LineWidth',line_wid.norm)
            for j = 1:size(keep_inds,2)
                if keep_tr(j)
                    visible = 'on';
                else
                    visible = 'off';
                end
                plot(t_snip,LE_V.LARP(:,j),'Color',colors.l_l,'LineWidth',line_wid.norm,'Visible',visible','Tag',['Cycle_',num2str(j)]);
                plot(t_snip,LE_V.RALP(:,j),'Color',colors.l_r,'LineWidth',line_wid.norm,'Visible',visible','Tag',['Cycle_',num2str(j)]);
                plot(t_snip,RE_V.LARP(:,j),'Color',colors.r_l,'LineWidth',line_wid.norm,'Visible',visible','Tag',['Cycle_',num2str(j)]);
                plot(t_snip,RE_V.RALP(:,j),'Color',colors.r_r,'LineWidth',line_wid.norm,'Visible',visible','Tag',['Cycle_',num2str(j)]);
                plot(t_snip,LE_V.LHRH(:,j),'Color',colors.l_z,'LineWidth',line_wid.norm,'Visible',visible','Tag',['Cycle_',num2str(j)]);
                plot(t_snip,RE_V.LHRH(:,j),'Color',colors.r_z,'LineWidth',line_wid.norm,'Visible',visible','Tag',['Cycle_',num2str(j)]);
            end
            hold off
            title(ha(4),['Accepted Cycles: ',num2str(sum(keep_tr)),' of ',num2str(length(keep_tr))])
            MakeCycAvg__plotCycAvg(ha(5),type,colors,CycAvg);
        case 2 %Just traces
            if isempty(ha) %first time running
                ha = gobjects(2,1);
                XLim_Long = [te(1) te(end)];
                %For plots with cycles to select
                x1 = 0.06;
                y1 = 0.05;
                y2 = 0.52;
                wid_x = 0.90;
                height_y = 0.43;
                ha(1) = subplot(2,1,1);
                ha(2) = subplot(2,1,2);
                ha(1).Position = [x1 y2 wid_x height_y];
                ha(2).Position = [x1 y1 wid_x height_y];
                %Link Axes
                linkaxes(ha,'x')
                %ha1
                set(ha(1),'XLim',XLim_Long,'YLim',YLim_Pos,'XTickLabel',[])
                title(ha(1),'Angular Position')
                ylabel(ha(1),'Position (deg)')
                %ha2
                set(ha(2),'XLim',XLim_Long,'YLim',YLim_Vel)
                title(ha(2),'Angular Velocity')
                xlabel(ha(2),'Time (s)')
                ylabel(ha(2),'Velocity (dps)')
            end
            %Raw and Filtered Position
            axes(ha(1))
            hold on 
            cla;
            h(1) = plot(ts,stim,'k');
            %Raw Data
            plot(ts,Data.LE_Position_X,'Color',colors.l_x_s,'LineWidth',line_wid.norm)
            plot(ts,Data.LE_Position_Y,'Color',colors.l_y_s,'LineWidth',line_wid.norm)
            plot(ts,Data.LE_Position_Z,'Color',colors.l_z_s,'LineWidth',line_wid.norm)
            plot(ts,Data.RE_Position_X,'Color',colors.r_x_s,'LineWidth',line_wid.norm)
            plot(ts,Data.RE_Position_Y,'Color',colors.r_y_s,'LineWidth',line_wid.norm)
            plot(ts,Data.RE_Position_Z,'Color',colors.r_z_s,'LineWidth',line_wid.norm)
            %Filtered Data
            h(2) = plot(ts,Data_In.Data_LE_Pos_X,'Color',colors.l_x,'LineWidth',line_wid.norm);
            h(3) = plot(ts,Data_In.Data_LE_Pos_Y,'Color',colors.l_y,'LineWidth',line_wid.norm);
            h(4) = plot(ts,Data_In.Data_LE_Pos_Z,'Color',colors.l_z,'LineWidth',line_wid.norm);
            h(5) = plot(ts,Data_In.Data_RE_Pos_X,'Color',colors.r_x,'LineWidth',line_wid.norm);
            h(6) = plot(ts,Data_In.Data_RE_Pos_Y,'Color',colors.r_y,'LineWidth',line_wid.norm);
            h(7) = plot(ts,Data_In.Data_RE_Pos_Z,'Color',colors.r_z,'LineWidth',line_wid.norm);
            hold off
            leg1 = legend(h,{'Stim','L X','L Y','L Z','R X','R Y','R Z'},'NumColumns',7);
            leg1.ItemTokenSize(1) = 7;
            axes(ha(2))
            hold on
            cla;
            plot(ts,stim,'k');
            %Raw Data
            plot(ts,Data_cal.LE_Vel_X,'.','Color',colors.l_x_s,'LineWidth',line_wid.norm)
            plot(ts,Data_cal.LE_Vel_Y,'.','Color',colors.l_y_s,'LineWidth',line_wid.norm)
            plot(ts,Data_cal.RE_Vel_X,'.','Color',colors.r_x_s,'LineWidth',line_wid.norm)
            plot(ts,Data_cal.RE_Vel_Y,'.','Color',colors.r_y_s,'LineWidth',line_wid.norm)
            plot(ts,Data_cal.LE_Vel_Z,'.','Color',colors.l_z_s,'LineWidth',line_wid.norm)
            plot(ts,Data_cal.RE_Vel_Z,'.','Color',colors.r_z_s,'LineWidth',line_wid.norm)
            %Filtered Data
            plot(ts,Data_calc.LE_Vel_X,'.','Color',colors.l_x,'LineWidth',line_wid.norm);
            plot(ts,Data_calc.LE_Vel_Y,'.','Color',colors.l_y,'LineWidth',line_wid.norm);
            plot(ts,Data_calc.LE_Vel_Z,'.','Color',colors.l_z,'LineWidth',line_wid.norm);
            plot(ts,Data_calc.RE_Vel_X,'.','Color',colors.r_x,'LineWidth',line_wid.norm);
            plot(ts,Data_calc.RE_Vel_Y,'.','Color',colors.r_y,'LineWidth',line_wid.norm);
            plot(ts,Data_calc.RE_Vel_Z,'.','Color',colors.r_z,'LineWidth',line_wid.norm);
            hold off
    end         
end