function MakeCycAvg__plotCycAvg(ax,type,colors,CycAvg)
    switch type
        case 1
            axes(ax)
            if length(CycAvg.t) > 1000
                s = round(linspace(1,length(CycAvg.t),1000));
            else
                s = 1:length(CycAvg.t);
            end
            hold on
            cla;
            plot(CycAvg.t(s),CycAvg.stim(s),'k');
            %Now add the fills and standard deviations and means
            %LE-LHRH
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s)))],colors.l_z_s)
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[CycAvg.lz_cycavg(s),fliplr((CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s)))],colors.l_z_s)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) + CycAvg.lz_cycstd(s),'Color',colors.l_z)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s) - CycAvg.lz_cycstd(s),'Color',colors.l_z)
            plot(CycAvg.t(s),CycAvg.lz_cycavg(s),'Color',colors.l_z,'LineWidth',2);
            %RE-LHRH
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s)))],colors.r_z_s)
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[CycAvg.rz_cycavg(s),fliplr((CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s)))],colors.r_z_s)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) + CycAvg.rz_cycstd(s),'Color',colors.r_z)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s) - CycAvg.rz_cycstd(s),'Color',colors.r_z)
            plot(CycAvg.t(s),CycAvg.rz_cycavg(s),'Color',colors.r_z,'LineWidth',2);
            %LE-LARP
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[CycAvg.ll_cycavg(s),fliplr((CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s)))],colors.l_l_s)
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[CycAvg.ll_cycavg(s),fliplr((CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s)))],colors.l_l_s)
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) + CycAvg.ll_cycstd(s),'Color',colors.l_l)
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s) - CycAvg.ll_cycstd(s),'Color',colors.l_l)
            plot(CycAvg.t(s),CycAvg.ll_cycavg(s),'Color',colors.l_l,'LineWidth',2);
            %RE-LARP
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[CycAvg.rl_cycavg(s),fliplr((CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s)))],colors.r_l_s)
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[CycAvg.rl_cycavg(s),fliplr((CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s)))],colors.r_l_s)
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) + CycAvg.rl_cycstd(s),'Color',colors.r_l)
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s) - CycAvg.rl_cycstd(s),'Color',colors.r_l)
            plot(CycAvg.t(s),CycAvg.rl_cycavg(s),'Color',colors.r_l,'LineWidth',2);
            %LE_RALP
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[CycAvg.lr_cycavg(s),fliplr((CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s)))],colors.l_r_s)
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[CycAvg.lr_cycavg(s),fliplr((CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s)))],colors.l_r_s)
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) + CycAvg.lr_cycstd(s),'Color',colors.l_r)
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s) - CycAvg.lr_cycstd(s),'Color',colors.l_r)
            plot(CycAvg.t(s),CycAvg.lr_cycavg(s),'Color',colors.l_r,'LineWidth',2);
            %RE-RALP
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[CycAvg.rr_cycavg(s),fliplr((CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s)))],colors.r_r_s)
            fill([CycAvg.t(s),fliplr(CycAvg.t(s))],[CycAvg.rr_cycavg(s),fliplr((CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s)))],colors.r_r_s)
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) + CycAvg.rr_cycstd(s),'Color',colors.r_r)
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s) - CycAvg.rr_cycstd(s),'Color',colors.r_r)
            plot(CycAvg.t(s),CycAvg.rr_cycavg(s),'Color',colors.r_r,'LineWidth',2)
            hold off
        case 2
            disp('No Cycle Average defined for this type')
        case 3
            axes(ax)
            if contains(CycAvg.name,{'LH','RH'})
                eye_c = 'z';
            elseif contains(CycAvg.name,{'LA','RP'})
                eye_c = 'l';
            elseif contains(CycAvg.name,{'RA','LP'})
                eye_c = 'r';
            end
            hold on
            cla;
            %Now add the fills and standard deviations and means
            %Head
            fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.head_cycavg,fliplr((CycAvg.head_cycavg + CycAvg.head_cycstd))],0.75*[1 1 1])
            fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.head_cycavg,fliplr((CycAvg.head_cycavg - CycAvg.head_cycstd))],0.75*[1 1 1])
            plot(CycAvg.t,CycAvg.head_cycavg + CycAvg.head_cycstd,'Color','k')
            plot(CycAvg.t,CycAvg.head_cycavg - CycAvg.head_cycstd,'Color','k')
            plot(CycAvg.t,CycAvg.head_cycavg,'Color','k','LineWidth',2);
            %Eye
            fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.eye_cycavg,fliplr((CycAvg.eye_cycavg + CycAvg.eye_cycstd))],colors.(['l_',eye_c,'_s']))
            fill([CycAvg.t,fliplr(CycAvg.t)],[CycAvg.eye_cycavg,fliplr((CycAvg.eye_cycavg - CycAvg.eye_cycstd))],colors.(['l_',eye_c,'_s']))
            plot(CycAvg.t,CycAvg.eye_cycavg + CycAvg.eye_cycstd,'Color',colors.(['l_',eye_c]))
            plot(CycAvg.t,CycAvg.eye_cycavg - CycAvg.eye_cycstd,'Color',colors.(['l_',eye_c]))
            plot(CycAvg.t,CycAvg.eye_cycavg,'Color',colors.(['l_',eye_c]),'LineWidth',2);
            hold off
     end
end