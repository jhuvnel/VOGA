function plotCycAvg(CycAvg,type,plot_fits,lrz_xyz)  
    %type should be 1 for a cycle average and 2 for activation data/vel step
    %plot_fits 0 or no input does not try to plot fits (if they are
    %available)
    %zlr_xyz 0 or no input = Z, LARP, RALP on graph and 1 is X Y Z
    figure;
    if nargin < 3
        plot_fits = 0;
    elseif plot_fits~=0 %Check to see if fits exist 
        fields = fieldnames(CycAvg);
        if ~ismember('lr_cycavg_fit',fields)
            disp('No plot fits found. Plotted without plot fits')
            plot_fits = 0;
        end
    end
    if nargin < 4
        lrz_xyz = 0;
    end
    %figure;
    %Colors
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
    switch type
        case 1 %Sine or other cycle average
            fields = fieldnames(CycAvg);
            if ~ismember('t',fields)
                CycAvg.t = reshape(0:1/CycAvg.Fs:(length(CycAvg.ll_cycavg)-1)/CycAvg.Fs,[],1);
            else
                CycAvg.t = reshape(CycAvg.t,[],1);
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
            h(1) = plot(CycAvg.t(s),CycAvg.stim(s),'k');
            hold on
            %Now add the fills and standard deviations and means
            if ~lrz_xyz  
                %LE-LARP
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ll_cycavg(s),fliplr((CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s)))],colors.l_l_s)
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ll_cycavg(s),fliplr((CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s)))],colors.l_l_s)
                plot(CycAvg.t(s),CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s),'Color',colors.l_l)
                plot(CycAvg.t(s),CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s),'Color',colors.l_l)
                h(2) = plot(CycAvg.t(s),CycAvg.ll_cycavg(s),'Color',colors.l_l,'LineWidth',2);
                %RE-LARP
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rl_cycavg(s),fliplr((CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s)))],colors.r_l_s)
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rl_cycavg(s),fliplr((CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s)))],colors.r_l_s)
                plot(CycAvg.t(s),CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s),'Color',colors.r_l)
                plot(CycAvg.t(s),CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s),'Color',colors.r_l)
                h(3) = plot(CycAvg.t(s),CycAvg.rl_cycavg(s),'Color',colors.r_l,'LineWidth',2);
                %LE_RALP
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lr_cycavg(s),fliplr((CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s)))],colors.l_r_s)
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lr_cycavg(s),fliplr((CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s)))],colors.l_r_s)
                plot(CycAvg.t(s),CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s),'Color',colors.l_r)
                plot(CycAvg.t(s),CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s),'Color',colors.l_r)
                h(4) = plot(CycAvg.t(s),CycAvg.lr_cycavg(s),'Color',colors.l_r,'LineWidth',2);
                %RE-RALP
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rr_cycavg(s),fliplr((CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s)))],colors.r_r_s)
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rr_cycavg(s),fliplr((CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s)))],colors.r_r_s)
                plot(CycAvg.t(s),CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s),'Color',colors.r_r)
                plot(CycAvg.t(s),CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s),'Color',colors.r_r)
                h(5) = plot(CycAvg.t(s),CycAvg.rr_cycavg(s),'Color',colors.r_r,'LineWidth',2);
                leg_lab = {'Stim','Left LARP','Right LARP','Left RALP','Right RALP','Left LHRH','Right LHRH'};
            else
                %LE-X
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lx_cycavg(s),fliplr((CycAvg.lx_cycavg(s) + CycAvg.lx_cycstd(s)))],colors.l_x_s)
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lx_cycavg(s),fliplr((CycAvg.lx_cycavg(s) - CycAvg.lx_cycstd(s)))],colors.l_x_s)
                plot(CycAvg.t(s),CycAvg.lx_cycavg(s) + CycAvg.lx_cycstd(s),'Color',colors.l_x)
                plot(CycAvg.t(s),CycAvg.lx_cycavg(s) - CycAvg.lx_cycstd(s),'Color',colors.l_x)
                h(2) = plot(CycAvg.t(s),CycAvg.lx_cycavg(s),'Color',colors.l_x,'LineWidth',2);
                %RE-X
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rx_cycavg(s),fliplr((CycAvg.rx_cycavg(s) + CycAvg.rx_cycstd(s)))],colors.r_x_s)
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rx_cycavg(s),fliplr((CycAvg.rx_cycavg(s) - CycAvg.rx_cycstd(s)))],colors.r_x_s)
                plot(CycAvg.t(s),CycAvg.rx_cycavg(s) + CycAvg.rx_cycstd(s),'Color',colors.r_x)
                plot(CycAvg.t(s),CycAvg.rx_cycavg(s) - CycAvg.rx_cycstd(s),'Color',colors.r_x)
                h(3) = plot(CycAvg.t(s),CycAvg.rx_cycavg(s),'Color',colors.r_x,'LineWidth',2);
                %LE_RALP
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ly_cycavg(s),fliplr((CycAvg.ly_cycavg(s) + CycAvg.ly_cycstd(s)))],colors.l_y_s)
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ly_cycavg(s),fliplr((CycAvg.ly_cycavg(s) - CycAvg.ly_cycstd(s)))],colors.l_y_s)
                plot(CycAvg.t(s),CycAvg.ly_cycavg(s) + CycAvg.ly_cycstd(s),'Color',colors.l_y)
                plot(CycAvg.t(s),CycAvg.ly_cycavg(s) - CycAvg.ly_cycstd(s),'Color',colors.l_y)
                h(4) = plot(CycAvg.t(s),CycAvg.ly_cycavg(s),'Color',colors.l_y,'LineWidth',2);
                %RE-RALP
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ry_cycavg(s),fliplr((CycAvg.ry_cycavg(s) + CycAvg.ry_cycstd(s)))],colors.r_y_s)
                fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.ry_cycavg(s),fliplr((CycAvg.ry_cycavg(s) - CycAvg.ry_cycstd(s)))],colors.r_y_s)
                plot(CycAvg.t(s),CycAvg.ry_cycavg(s) + CycAvg.ry_cycstd(s),'Color',colors.r_y)
                plot(CycAvg.t(s),CycAvg.ry_cycavg(s) - CycAvg.ry_cycstd(s),'Color',colors.r_y)
                h(5) = plot(CycAvg.t(s),CycAvg.ry_cycavg(s),'Color',colors.r_y,'LineWidth',2);
                leg_lab = {'Stim','Left X','Right X','Left Y','Right Y','Left Z','Right Z'};
            end
            %LE-LHRH
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s)))],colors.l_z_s)
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s)))],colors.l_z_s)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
            h(6) = plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
            %RE-LHRH
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s)))],colors.r_z_s)
            fill([CycAvg.t(s)',fliplr(CycAvg.t(s)')],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s)))],colors.r_z_s)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
            h(7) = plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
            if plot_fits
                if ~lrz_xyz
                    plot(CycAvg.t(s),CycAvg.ll_cycavg_fit(s),'--','Color',colors.l_l,'LineWidth',2)
                    plot(CycAvg.t(s),CycAvg.rl_cycavg_fit(s),'--','Color',colors.r_l,'LineWidth',2)
                    plot(CycAvg.t(s),CycAvg.lr_cycavg_fit(s),'--','Color',colors.l_r,'LineWidth',2)
                    plot(CycAvg.t(s),CycAvg.rr_cycavg_fit(s),'--','Color',colors.r_r,'LineWidth',2)
                else
                    plot(CycAvg.t(s),CycAvg.lx_cycavg_fit(s),'--','Color',colors.l_x,'LineWidth',2)
                    plot(CycAvg.t(s),CycAvg.rx_cycavg_fit(s),'--','Color',colors.r_x,'LineWidth',2)
                    plot(CycAvg.t(s),CycAvg.ly_cycavg_fit(s),'--','Color',colors.l_y,'LineWidth',2)
                    plot(CycAvg.t(s),CycAvg.ry_cycavg_fit(s),'--','Color',colors.r_y,'LineWidth',2)
                end
                plot(CycAvg.t(s),CycAvg.lz_cycavg_fit(s),'--','Color',colors.l_z,'LineWidth',2)
                plot(CycAvg.t(s),CycAvg.rz_cycavg_fit(s),'--','Color',colors.r_z,'LineWidth',2)
            end         
            hold off
            if ismember('name',fields)
                title(strrep(strrep(CycAvg.name,'-',' '),'_',' '))
            else
                title(strrep(CycAvg.info.dataType,'-',' '))
            end
            xlabel('Time (s)')
            ylabel('Angular Velocity (dps)')
            legend(h,leg_lab)
        case 2 %Velocity step or activation data
            ts = CycAvg.t;
            h2(1) = plot(ts,CycAvg.stim,'k');
            hold on
            if ~lrz_xyz
                h2(2) = plot(ts,CycAvg.ll_cyc,'.','Color',colors.l_l);
                h2(3) = plot(ts,CycAvg.rl_cyc,'.','Color',colors.r_l);
                h2(4) = plot(ts,CycAvg.lr_cyc,'.','Color',colors.l_r);
                h2(5) = plot(ts,CycAvg.rr_cyc,'.','Color',colors.r_r);
                leg_lab = {'Stim','Left LARP','Right LARP','Left RALP','Right RALP','Left LHRH','Right LHRH'};
            else
                h2(2) = plot(ts,CycAvg.lx_cyc,'.','Color',colors.l_x);
                h2(3) = plot(ts,CycAvg.rx_cyc,'.','Color',colors.r_x);
                h2(4) = plot(ts,CycAvg.ly_cyc,'.','Color',colors.l_y);
                h2(5) = plot(ts,CycAvg.ry_cyc,'.','Color',colors.r_y);
                leg_lab = {'Stim','Left X','Right X','Left Y','Right Y','Left Z','Right Z'};
            end
            h2(6) = plot(ts,CycAvg.lz_cyc,'.','Color',colors.l_z);
            h2(7) = plot(ts,CycAvg.lz_cyc,'.','Color',colors.r_z);            
            if plot_fits
                if ~lrz_xyz
                    plot(ts,CycAvg.ll_cycavg_fit,'Color',colors.l_l,'LineWidth',2)
                    plot(ts,CycAvg.rl_cycavg_fit,'Color',colors.r_l,'LineWidth',2)
                    plot(ts,CycAvg.lr_cycavg_fit,'Color',colors.l_r,'LineWidth',2)
                    plot(ts,CycAvg.rr_cycavg_fit,'Color',colors.r_r,'LineWidth',2)
                else
                    plot(ts,CycAvg.lx_cycavg_fit,'Color',colors.l_x,'LineWidth',2)
                    plot(ts,CycAvg.rx_cycavg_fit,'Color',colors.r_x,'LineWidth',2)
                    plot(ts,CycAvg.ly_cycavg_fit,'Color',colors.l_y,'LineWidth',2)
                    plot(ts,CycAvg.ry_cycavg_fit,'Color',colors.r_y,'LineWidth',2)
                end
                plot(ts,CycAvg.lz_cycavg_fit,'Color',colors.l_z,'LineWidth',2)
                plot(ts,CycAvg.rz_cycavg_fit,'Color',colors.r_z,'LineWidth',2)
            end             
            hold off
            legend(h2,leg_lab)
            fields = fieldnames(CycAvg);
            if ismember('name',fields)
                title(strrep(strrep(CycAvg.name,'-',' '),'_',' '))
            else
                title(strrep(CycAvg.info.dataType,'-',' '))
            end
    end
end