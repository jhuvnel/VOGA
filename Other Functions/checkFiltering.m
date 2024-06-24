% This will plot each filtered accepted cycle in front of its non-filtered
% velocity trace. Press "enter" to continue and any other key+enter to save
% it as a file to redo.

load('VNELcolors.mat')
Cyc_Path = [cd,filesep,'Cycle Averages'];
cycavg_files = extractfield(dir([Cyc_Path,filesep,'CycAvg_*.mat']),'name');
disp([num2str(length(cycavg_files)),' Cycle Average files found'])
fig = figure(1);
set(fig,'Color','w','Units','inches','Position',[0.5 0.5 10 4])
redo = false(length(cycavg_files),1);
for i = 1:length(cycavg_files)
    load([Cyc_Path,filesep,cycavg_files{i}])
    if isfield(CycAvg,'Data_filtvel')
        keep_inds = CycAvg.Data_allcyc.keep_inds(:,CycAvg.keep_tr);
        clf;
        annotation('textbox',[0 0.8 1 0.2],'String',strrep(cycavg_files{i},'_','-'),'HorizontalAlignment','center','EdgeColor','none')
        if contains(cycavg_files{i},{'X','Y'})
            traces = {'x','Y','Z'};
        else
            traces = {'LARP','RALP','Z'};
        end
        for t = 1:length(traces)
            filt_tr = [CycAvg.Data_filtvel.(['LE_Vel_',traces{t}])(keep_inds),CycAvg.Data_filtvel.(['RE_Vel_',traces{t}])(keep_inds)];
            YLim = [25*floor(min(min(filt_tr))/25),25*ceil(max(max(filt_tr))/25)];
            subplot(1,length(traces),t)
            hold on
            plot(CycAvg.t,CycAvg.Data_rawvel.(['LE_Vel_',traces{t}])(keep_inds),'Color',colors.(['l_',lower(traces{t}(1)),'_s']))
            plot(CycAvg.t,CycAvg.Data_rawvel.(['RE_Vel_',traces{t}])(keep_inds),'Color',colors.(['r_',lower(traces{t}(1)),'_s']))
            plot(CycAvg.t,CycAvg.Data_filtvel.(['LE_Vel_',traces{t}])(keep_inds),'Color',colors.(['l_',lower(traces{t}(1))]))
            plot(CycAvg.t,CycAvg.Data_filtvel.(['RE_Vel_',traces{t}])(keep_inds),'Color',colors.(['r_',lower(traces{t}(1))]))
            hold off
            set(gca,'YLim',YLim)
            xlabel('Time (s)')
            ylabel('Angular Velocity (deg/s)')
        end
        redo = input('','s');
        if ~strcmp(redo,'')
            redo(i) = true;
        end
    end
end
disp(cycavg_files(redo))