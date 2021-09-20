%% Make Full Plot Cyc Avg
function ha = MakeCycAvg__plotFullCycAvg(ha,type,colors,line_wid,YLim_Pos,YLim_Vel,traces_pos,traces_vel,CycAvg)
    switch type
        case 1 %Position, Velocity, and Cycle Averaging
            te = CycAvg.Data_rawpos.te;
            ts = CycAvg.Data_rawpos.ts;
            t_s = CycAvg.t;
            stim = CycAvg.Data_rawpos.stim;
            stims = CycAvg.stim;
            keep_inds = CycAvg.Data_allcyc.keep_inds;
            keep_tr = CycAvg.keep_tr;                       
            if isempty(ha) %first time running
                ha = gobjects(5,1);
                XLim_Long = [te(1) te(end)];
                XLim_Short = [t_s(1) t_s(end)];
                %For plots with cycles to select
                x = [0.06;0.37;0.68];
                y = [0.045;0.38;0.67]; 
                wid_x_s = 0.30;
                wid_x_b = 0.92;
                height_y = 0.27;
                ha(1) = subplot(3,3,[1 2 3]);
                ha(2) = subplot(3,3,[4 5 6]);
                ha(3) = subplot(3,3,7);
                ha(4) = subplot(3,3,8);
                ha(5) = subplot(3,3,9);
                ha(1).Position = [x(1) y(3) wid_x_b height_y];
                ha(2).Position = [x(1) y(2) wid_x_b height_y];
                ha(3).Position = [x(1) y(1) wid_x_s height_y];
                ha(4).Position = [x(2) y(1) wid_x_s height_y];
                ha(5).Position = [x(3) y(1) wid_x_s height_y];
                %Link Axes
                linkaxes(ha([1,2]),'x')
                linkaxes(ha([2,3]),'y')
                linkaxes(ha([3,4,5]),'xy')
                %Set Labels
                set(ha(1),'XLim',XLim_Long,'YLim',YLim_Pos)
                ylabel(ha(1),'Angular Position (deg)','FontWeight','bold')
                set(ha(2),'XLim',XLim_Long,'YLim',YLim_Vel)
                ylabel(ha(2),'Angular Velocity (dps)','FontWeight','bold')
                xlabel(ha(2),'Time (s)')
                set(ha(3),'XLim',XLim_Short,'YLim',YLim_Vel)
                ylabel(ha(3),'Velocity (dps)')
                title(ha(3),'All Filtered Cycles')
                set(ha(4),'XLim',XLim_Short,'YLim',YLim_Vel,'YTickLabel',[])
                xlabel(ha(4),'Time (s)')
                set(ha(5),'XLim',XLim_Short,'YLim',YLim_Vel,'YTickLabel',[])
                title(ha(5),'Cycle Averages')
            end
            %Cycle Aligned Raw and Filtered Position
            axes(ha(1))
            hold on
            cla;
            for j = 1:size(keep_inds,2)
                if keep_tr(j)
                    fill([te(keep_inds(1,j)),te(keep_inds(end,j)),te(keep_inds(end,j)),te(keep_inds(1,j))]',[500,500,-500,-500]',colors.cyc_keep);
                else
                    fill([te(keep_inds(1,j)),te(keep_inds(end,j)),te(keep_inds(end,j)),te(keep_inds(1,j))]',[500,500,-500,-500]',colors.cyc_rm);
                end
            end
            %Make cycle labels
            cyc_num_labs = 1:length(keep_tr);
            x_tick = te(round(floor(mean(keep_inds)),0));
            set(ha(1),'XTick',x_tick,'XTickLabel',cyc_num_labs,'Xaxislocation','top')
            %Plot Cycles
            % Position Data
            h = gobjects(1,length(traces_pos)+1);
            h(1) = plot(ts,stim,'k','LineWidth',line_wid.norm);
            % Raw Traces
            for i = 1:length(traces_pos)
                if isfield(CycAvg.Data_rawpos,[traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end)])
                    plot(te,CycAvg.Data_rawpos.([traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end)]),'Color',colors.([lower(traces_pos{i}(1)),'_',lower(traces_pos{i}(2)),'_s']),'LineWidth',line_wid.norm)
                end
            end
            % Filtered Traces
            for i = 1:length(traces_pos)
                if isfield(CycAvg.Data_filtpos,[traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end)])
                    h(i+1) = plot(ts,CycAvg.Data_filtpos.([traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end)]),'Color',colors.([lower(traces_pos{i}(1)),'_',lower(traces_pos{i}(2))]),'LineWidth',line_wid.norm);
                else
                    h(i+1) = plot(NaN,NaN,'Color',colors.([lower(traces_pos{i}(1)),'_',lower(traces_pos{i}(2))]));
                end
            end
            hold off
            leg1 = legend(h,[{'Stim'};reshape(traces_pos,[],1)],'NumColumns',length(traces_pos)+1);
            leg1.ItemTokenSize(1) = 7;
            % Velocity Data
            axes(ha(2))
            hold on
            cla; 
            for j = 1:size(keep_inds,2)
                if keep_tr(j)
                    fill([ts(keep_inds(1,j)),ts(keep_inds(end,j)),ts(keep_inds(end,j)),ts(keep_inds(1,j))]',[500,500,-500,-500]',colors.cyc_keep,'Tag',['Cycle_',num2str(j)]);
                else
                    fill([ts(keep_inds(1,j)),ts(keep_inds(end,j)),ts(keep_inds(end,j)),ts(keep_inds(1,j))]',[500,500,-500,-500]',colors.cyc_rm,'Tag',['Cycle_',num2str(j)]);
                end
            end
            h1 = gobjects(1,length(traces_vel)+1);
            h1(1) = plot(ts,stim,'k','LineWidth',line_wid.norm);
            % Raw Traces
            for i = 1:length(traces_vel)
                if isfield(CycAvg.Data_rawvel,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
                    plot(ts,CycAvg.Data_rawvel.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)]),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2)),'_s']),'LineWidth',line_wid.norm)
                end
            end
            % Filtered Traces
            for i = 1:length(traces_vel)
                if isfield(CycAvg.Data_filtvel,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
                    h1(i+1) = plot(ts,CycAvg.Data_filtvel.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)]),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm);
                else
                    h1(i+1) = plot(NaN,NaN,'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]));
                end
            end
            hold off
            leg2 = legend(h1,[{'Stim'};reshape(traces_vel,[],1)],'NumColumns',length(traces_vel)+1);
            leg2.ItemTokenSize(1) = 7;
            axes(ha(3))
            hold on
            cla;
            plot(t_s,stims,'k','LineWidth',line_wid.norm)
            %All Filtered Velocity Data
            for i = 1:length(traces_vel)
                if isfield(CycAvg.Data_allcyc,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
                    plot(t_s,CycAvg.Data_allcyc.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)]),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm)
                end
            end
            hold off
            axes(ha(4))
            hold on
            cla;
            plot(t_s,stims,'k','LineWidth',line_wid.norm)
            %Only Selected Filtered Velocity Data
            for i = 1:length(traces_vel)
                if isfield(CycAvg.Data_allcyc,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
                    plot(t_s,CycAvg.Data_allcyc.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])(:,keep_tr),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm)
                end
            end
            hold off
            title(ha(4),['Accepted Cycles: ',num2str(sum(keep_tr)),' of ',num2str(length(keep_tr))])
            axes(ha(5))
            if length(CycAvg.t) > 1000
                s = round(linspace(1,length(CycAvg.t),1000));
            else
                s = 1:length(CycAvg.t);
            end
            hold on
            cla;
            if isfield(CycAvg,'stim_cyc') %Mutiple head traces to show
                trace = 'stim';
                fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[(CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s)),fliplr((CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s)))],0.5*[1,1,1])
                plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s),'Color','k')
                plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s),'Color','k')
                plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s),'Color','k','LineWidth',line_wid.bold);
            else %Only 1
                plot(CycAvg.t(s),CycAvg.stim(s),'k','LineWidth',line_wid.bold);
            end
            %Only Selected Filtered Velocity Data
            for i = 1:length(traces_vel)
                trace = lower(traces_vel{i}(1:2));
                if isfield(CycAvg,[trace,'_cycavg'])&&isfield(CycAvg,[trace,'_cycstd'])
                    fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[(CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s)),fliplr((CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s)))],colors.([trace(1),'_',trace(2),'_s']))
                    plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s),'Color',colors.([trace(1),'_',trace(2)]))
                    plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s),'Color',colors.([trace(1),'_',trace(2)]))
                    plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s),'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',line_wid.bold);
                end
            end    
            hold off
        case 2 %Position and Velocity Traces but NO Cycle Averages
            te = CycAvg.Data_rawpos.te;
            ts = CycAvg.Data_rawpos.ts;
            stim = CycAvg.Data_rawpos.stim;
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
            %Plot Cycles
            % Position Data
            h = gobjects(1,length(traces_pos)+1);
            h(1) = plot(ts,stim,'k','LineWidth',line_wid.norm);
            % Raw Traces
            for i = 1:length(traces_pos)
                if isfield(CycAvg.Data_rawpos,[traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end)])
                    plot(te,CycAvg.Data_rawpos.([traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end)]),'Color',colors.([lower(traces_pos{i}(1)),'_',lower(traces_pos{i}(2)),'_s']),'LineWidth',line_wid.norm)
                end
            end
            % Filtered Traces
            for i = 1:length(traces_pos)
                if isfield(CycAvg.Data_filtpos,[traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end)])
                    h(i+1) = plot(ts,CycAvg.Data_filtpos.([traces_pos{i}(1),'E_Position_',traces_pos{i}(2:end)]),'Color',colors.([lower(traces_pos{i}(1)),'_',lower(traces_pos{i}(2))]),'LineWidth',line_wid.norm);
                else
                    h(i+1) = plot(NaN,NaN,'Color',colors.([lower(traces_pos{i}(1)),'_',lower(traces_pos{i}(2))]));
                end
            end
            hold off
            leg1 = legend(h,[{'Stim'};reshape(traces_pos,[],1)],'NumColumns',length(traces_pos)+1);
            leg1.ItemTokenSize(1) = 7;
            % Velocity Data
            axes(ha(2))
            hold on
            cla; 
            h1 = gobjects(1,length(traces_vel)+1);
            h1(1) = plot(ts,stim,'k','LineWidth',line_wid.norm);
            % Raw Traces
            for i = 1:length(traces_vel)
                if isfield(CycAvg.Data_rawvel,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
                    plot(ts,CycAvg.Data_rawvel.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)]),'.','Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2)),'_s']),'LineWidth',line_wid.norm)
                end
            end
            % Filtered Traces
            for i = 1:length(traces_vel)
                if isfield(CycAvg.Data_filtvel,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
                    h1(i+1) = plot(ts,CycAvg.Data_filtvel.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)]),'.','Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm);
                else
                    h1(i+1) = plot(NaN,NaN,'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'.');
                end
            end
            hold off
            leg2 = legend(h1,[{'Stim'};reshape(traces_vel,[],1)],'NumColumns',length(traces_vel)+1);
            leg2.ItemTokenSize(1) = 7;
        case 3 %No Position Traces, Just Velocity and Cycle Averaging
            te = CycAvg.Data_rawvel.t;
            ts = CycAvg.Data_rawvel.t;
            t_s = CycAvg.t;
            stim = CycAvg.Data_rawvel.stim;
            stims = CycAvg.stim;
            keep_inds = CycAvg.Data_allcyc.keep_inds;
            keep_tr = CycAvg.keep_tr;
            detec_tr = CycAvg.detec_tr;
            if isempty(ha) %first time running
                ha = gobjects(3,1);
                XLim_Long = [te(1) te(end)];
                XLim_Short = [t_s(1) t_s(end)];
                %For plots with cycles to select
                x = [0.06;0.37;0.68];
                y = [0.045;0.51];
                wid_x_s = 0.30;
                wid_x_b = 0.92;
                height_y = 0.43;
                ha(1) = subplot(2,3,[1 2 3]);
                ha(2) = subplot(2,3,4);
                ha(3) = subplot(2,3,5);
                ha(4) = subplot(2,3,6);
                ha(1).Position = [x(1) y(2) wid_x_b height_y];
                ha(2).Position = [x(1) y(1) wid_x_s height_y];
                ha(3).Position = [x(2) y(1) wid_x_s height_y];
                ha(4).Position = [x(3) y(1) wid_x_s height_y];
                %Link Axes
                linkaxes(ha([1,2]),'y')
                linkaxes(ha([2,3,4]),'xy')
                %Make cycle labels
                cyc_num_labs = 1:length(keep_tr);
                x_tick = te(round(floor(mean(keep_inds)),0));
                %Set Labels
                %ha1
                set(ha(1),'XLim',XLim_Long,'YLim',YLim_Vel,'XTick',x_tick,'XTickLabel',cyc_num_labs,'Xaxislocation','top')
                ylabel(ha(1),'Angular Velocity (dps)')
                %ha2
                set(ha(2),'XLim',XLim_Short,'YLim',YLim_Vel)
                ylabel(ha(2),'Velocity (dps)')
                title(ha(2),'All Cycles')
                xlabel(ha(2),'Time (s)')
                %ha3
                set(ha(3),'XLim',XLim_Short,'YLim',YLim_Vel,'YTickLabel',[])
                xlabel(ha(3),'Time (s)')
                %ha4
                set(ha(4),'XLim',XLim_Short,'YLim',YLim_Vel,'YTickLabel',[])
                xlabel(ha(4),'Time (s)')
                title(ha(4),'Cycle Average')
            end
            %Cycle Aligned Raw and Filtered Velocity
            axes(ha(1))
            hold on
            cla;
            for j = 1:size(keep_inds,2)
                if keep_tr(j)
                    fill([te(keep_inds(1,j)),te(keep_inds(end,j)),te(keep_inds(end,j)),te(keep_inds(1,j))]',[500,500,-500,-500]',colors.cyc_keep);
                else
                    fill([te(keep_inds(1,j)),te(keep_inds(end,j)),te(keep_inds(end,j)),te(keep_inds(1,j))]',[500,500,-500,-500]',colors.cyc_rm);
                end 
            end     
            if ~isempty(detec_tr)
                xline(te([keep_inds(1,detec_tr),keep_inds(end,detec_tr)]),'-b','LineWidth',2)     
            end
            h = gobjects(1,length(traces_vel)+1);
            h(1) = plot(ts,stim,'k','LineWidth',line_wid.norm);
            % Raw Traces
            for i = 1:length(traces_vel)
                if isfield(CycAvg.Data_rawvel,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
                    plot(ts,CycAvg.Data_rawvel.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)]),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2)),'_s']),'LineWidth',line_wid.norm)
                end
            end
            % Filtered Traces
            for i = 1:length(traces_vel)
                if isfield(CycAvg.Data_filtvel,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
                    h(i+1) = plot(ts,CycAvg.Data_filtvel.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)]),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm);
                else
                    h(i+1) = plot(NaN,NaN,'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]));
                end
            end
            hold off
            leg2 = legend(h,[{'Stim'};reshape(traces_vel,[],1)],'NumColumns',length(traces_vel)+1);
            leg2.ItemTokenSize(1) = 7;
            axes(ha(2))
            hold on
            cla;
            plot(t_s,stims,'k','LineWidth',line_wid.norm)
            %Cycle Averaged Velocity Data
            for i = 1:length(traces_vel)
                if isfield(CycAvg.Data_allcyc,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
                    plot(t_s,CycAvg.Data_rawvel.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])(CycAvg.Data_allcyc.keep_inds),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2)),'_s']),'LineWidth',line_wid.norm)   
                    plot(t_s,CycAvg.Data_allcyc.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)]),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm)   
                end
            end
            hold off
            axes(ha(3))
            hold on
            cla;
            plot(t_s,stims,'k','LineWidth',line_wid.norm)
            %Only Selected Filtered Velocity Data
            for i = 1:length(traces_vel)
                if isfield(CycAvg.Data_allcyc,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end),'_QPR'])
                    plot(t_s,CycAvg.Data_allcyc.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])(:,keep_tr),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2)),'_s']),'LineWidth',line_wid.norm)
                    plot(t_s,CycAvg.Data_allcyc.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end),'_QPR'])(:,keep_tr),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm)
                elseif isfield(CycAvg.Data_allcyc,[traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])
                    plot(t_s,CycAvg.Data_allcyc.([traces_vel{i}(1),'E_Vel_',traces_vel{i}(2:end)])(:,keep_tr),'Color',colors.([lower(traces_vel{i}(1)),'_',lower(traces_vel{i}(2))]),'LineWidth',line_wid.norm)   
                end
            end
            hold off
            title(ha(3),['Accepted Cycles: ',num2str(sum(keep_tr)),' of ',num2str(length(keep_tr))])
            axes(ha(4))
            if length(CycAvg.t) > 1000
                s = round(linspace(1,length(CycAvg.t),1000));
            else
                s = 1:length(CycAvg.t);
            end
            hold on
            cla;
            if isfield(CycAvg,'stim_cyc') %Mutiple head traces to show
                trace = 'stim';
                fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[(CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s)),fliplr((CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s)))],0.5*[1,1,1])
                plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s),'Color','k')
                plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s),'Color','k')
                plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s),'Color','k','LineWidth',line_wid.bold);
            else %Only 1
                plot(CycAvg.t(s),CycAvg.stim(s),'k','LineWidth',line_wid.bold);
            end
            %Only Selected Filtered Velocity Data
            for i = 1:length(traces_vel)
                trace = lower(traces_vel{i}(1:2));
                if isfield(CycAvg,[trace,'_cyc_fit'])
                    plot(CycAvg.t(s),CycAvg.([trace,'_cyc'])(:,s),'Color',colors.([trace(1),'_',trace(2),'_s']),'LineWidth',0.5);
                    fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[(CycAvg.([trace,'_cycavg_fit'])(s) - CycAvg.([trace,'_cycstd_fit'])(s)),fliplr((CycAvg.([trace,'_cycavg_fit'])(s) + CycAvg.([trace,'_cycstd_fit'])(s)))],colors.([trace(1),'_',trace(2),'_s']))
                    plot(CycAvg.t(s),CycAvg.([trace,'_cycavg_fit'])(s) + CycAvg.([trace,'_cycstd_fit'])(s),'Color',colors.([trace(1),'_',trace(2)]))
                    plot(CycAvg.t(s),CycAvg.([trace,'_cycavg_fit'])(s) - CycAvg.([trace,'_cycstd_fit'])(s),'Color',colors.([trace(1),'_',trace(2)]))
                    plot(CycAvg.t(s),CycAvg.([trace,'_cycavg_fit'])(s),'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',line_wid.bold);
                elseif isfield(CycAvg,[trace,'_cyc'])
                    fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[(CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s)),fliplr((CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s)))],colors.([trace(1),'_',trace(2),'_s']))
                    plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s),'Color',colors.([trace(1),'_',trace(2)]))
                    plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s),'Color',colors.([trace(1),'_',trace(2)]))
                    plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s),'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',line_wid.bold);
                end
            end    
            hold off
    end         
end