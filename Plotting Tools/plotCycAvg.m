function plotCycAvg(CycAvg,plot_fits,lrz_xyz)  
    %First, assign type
    if contains(CycAvg.info.dataType,{'Impulse'}) %All head and cycles
        type = 3;
    elseif contains(CycAvg.info.dataType,{'Activation','Step'}) %No cycle averaging
        type = 2;
    else
        type = 1;
    end
    %plot_fits 0 or no input does not try to plot fits (if they are
    %available)
    %zlr_xyz 0 or no input = Z, LARP, RALP on graph and 1 is X Y Z
    if nargin < 2
        plot_fits = 0;
    end
    if nargin < 3 || isempty(lrz_xyz)
        traces = {'ll','rl','lr','rr','lx','rx','ly','ry','lz','rz'};
    elseif contains(lrz_xyz,{'lrz','LRZ'}) 
        traces = {'ll','rl','lr','rr','lz','rz'};
    else
        traces = {'lx','rx','ly','ry','lz','rz'};
    end
    %Colors
    % Normal colors
    load('VNELcolors.mat','colors')
    figure('Color',[1,1,1]);
    switch type
        case 1 %Sine or other cycle average
            fields = fieldnames(CycAvg);
            if ~ismember('t',fields)
                CycAvg.t = reshape(0:1/CycAvg.Fs:(length(CycAvg.ll_cycavg)-1)/CycAvg.Fs,[],1);
            else
                CycAvg.t = reshape(CycAvg.t,1,[]);
            end
            if length(CycAvg.t) > 1000
                s = round(linspace(1,length(CycAvg.t),1000));
            else
                s = 1:length(CycAvg.t);
            end
            [aa,ab] = size(CycAvg.stim);
            if aa == length(CycAvg.t) && ab ~=1
                CycAvg.stim = mean(CycAvg.stim,2)';
            elseif ab == length(CycAvg.t) && aa ~=1
                CycAvg.stim = mean(CycAvg.stim,1);
            end
            h = gobjects(length(traces)+1,1);
            %h(1) = plot(CycAvg.t(s),CycAvg.stim(s),'k','LineWidth',2);
            h(1) = plot(CycAvg.t(s),CycAvg.stim(s),'--','Color',colors.l_y,'LineWidth',2);
            hold on
            for i = 1:length(traces)
                trace = traces{i};
                h(i+1) = plot(NaN,NaN,'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',2);
                if isfield(CycAvg,[trace,'_cycavg'])&&isfield(CycAvg,[trace,'_cycstd'])
                    fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[(CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s)),fliplr((CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s)))],colors.([trace(1),'_',trace(2),'_s']))
                    plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) + CycAvg.([trace,'_cycstd'])(s),'Color',colors.([trace(1),'_',trace(2)]))
                    plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s) - CycAvg.([trace,'_cycstd'])(s),'Color',colors.([trace(1),'_',trace(2)]))
                    plot(CycAvg.t(s),CycAvg.([trace,'_cycavg'])(s),'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',2);
                end
            end
            if plot_fits
                for i = 1:length(traces)
                    trace = traces{i};
                    if isfield(CycAvg,[trace,'_cycavg_fit'])
                        plot(CycAvg.t(s),CycAvg.([trace,'_cycavg_fit'])(s),'--','Color',colors.([trace(1),'_',trace(2)]),'LineWidth',2)
                    end
                end
            end   
            hold off
            if ismember('name',fields)
                fig_title = strrep(strrep(CycAvg.name,'-',' '),'_',' ');
            else
                fig_title = strrep(CycAvg.info.dataType,'-',' ');
            end
            if contains(fig_title,'[')&&contains(fig_title,']')
                fig_title(strfind(fig_title,'['):strfind(fig_title,']')) = strrep(fig_title(strfind(fig_title,'['):strfind(fig_title,']')),' ','-');
            end
            title(fig_title)
            xlabel('Time (s)')
            ylabel('Angular Velocity (dps)')
            legend(h,[{'Stim'},upper(traces)])
        case 2 %Velocity step or activation data
            ts = CycAvg.t;
            h2 = gobjects(length(traces)+1,1);
            h2(1) = plot(ts,CycAvg.stim,'k');
            hold on
            for i = 1:length(traces)
                trace = traces{i};
                h2(i+1) = plot(NaN,NaN,'.','Color',colors.([trace(1),'_',trace(2)]),'LineWidth',2);
                if isfield(CycAvg,[trace,'_cycavg'])&&isfield(CycAvg,[trace,'_cycstd'])
                    plot(ts,CycAvg.([trace,'_cycavg']),'.','Color',colors.([trace(1),'_',trace(2)]),'LineWidth',2);
                end
            end
            if plot_fits
                for i = 1:length(traces)
                    trace = traces{i};
                    if isfield(CycAvg,[trace,'_cycavg_fit'])
                        plot(ts,CycAvg.([trace,'_cycavg_fit']),'--','Color',colors.([trace(1),'_',trace(2)]),'LineWidth',2)
                    end
                end
            end  
            hold off
            legend(h2,[{'Stim'},upper(traces)])
            fields = fieldnames(CycAvg);
            if ismember('name',fields)
                fig_title = strrep(strrep(CycAvg.name,'-',' '),'_',' ');
            else
                fig_title = strrep(CycAvg.info.dataType,'-',' ');
            end
            if contains(fig_title,'[')&&contains(fig_title,']')
                fig_title(strfind(fig_title,'['):strfind(fig_title,']')) = strrep(fig_title(strfind(fig_title,'['):strfind(fig_title,']')),' ','-');
            end
            title(fig_title)
        case 3 %Impulse data, make the head trace always positive and invert the eye as needed, show all cycles
            fields = fieldnames(CycAvg);
            if ~ismember('t',fields)
                CycAvg.t = reshape(0:1/CycAvg.Fs:(length(CycAvg.ll_cycavg)-1)/CycAvg.Fs,[],1);
            else
                CycAvg.t = reshape(CycAvg.t,1,[]);
            end
            if -min(mean(CycAvg.stim))>max(mean(CycAvg.stim))
                invh = -1;
                inve = 1;
            else
                invh = 1;
                inve = -1;
            end
            h = gobjects(length(traces)+1,1);
            h(1) = plot(NaN,NaN,'k','LineWidth',1);
            hold on
            plot(CycAvg.t,invh*CycAvg.stim,'k','LineWidth',0.5);
            for i = 1:length(traces)                
                trace = traces{i};
                h(i+1) = plot(NaN,NaN,'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',1);
                if plot_fits
                    if isfield(CycAvg,[trace,'_cyc_fit'])
                        plot(CycAvg.t,inve*CycAvg.([trace,'_cyc']),'Color',colors.([trace(1),'_',trace(2),'_s']),'LineWidth',0.5);
                        plot(CycAvg.t,inve*CycAvg.([trace,'_cyc_fit']),'-','Color',colors.([trace(1),'_',trace(2)]),'LineWidth',1)
                    end
                else
                    if isfield(CycAvg,[trace,'_cyc'])
                        plot(CycAvg.t,inve*CycAvg.([trace,'_cyc']),'Color',colors.([trace(1),'_',trace(2)]),'LineWidth',0.5);
                    end
                end
            end   
            hold off
            if ismember('name',fields)
                fig_title = strrep(strrep(CycAvg.name,'-',' '),'_',' ');
            else
                fig_title = strrep(CycAvg.info.dataType,'-',' ');
            end
            if contains(fig_title,'[')&&contains(fig_title,']')
                fig_title(strfind(fig_title,'['):strfind(fig_title,']')) = strrep(fig_title(strfind(fig_title,'['):strfind(fig_title,']')),' ','-');
            end
            title(fig_title)
            xlabel('Time (s)')
            ylabel('Angular Velocity (dps)')
            legend(h,[{'Stim'},upper(traces)])
    end
end