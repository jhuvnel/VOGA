%% plotRotaryChairVisitSummary
% just gain/phase over freq for all exp type, will be moved to new function soon

function plotRotaryChairVisitSummary(path,Cyc_Path,code_Path,version,Experimenter)
%% Colors and Normative Data
%Implanted ear--current as of sub 9

%Normative Rotary Chair Data (100dps sine)
%Wall et. al. 1984 on patients 50-69
freq = [0.005 0.010 0.020 0.050 0.100 0.200 0.500 1.000];
norm_gain_m = [0.2175	0.3576	0.4803	0.5634	0.56	0.5445	0.5678	0.698];
norm_gain_std = [0.0757	0.1235	0.1474	0.2243	0.2363	0.2509	0.2121	0.2777];
norm_phase_m = [70.9478	43.6413	25.9262	10.7906	4.084	-4.8461	-11.4084 -11.3412];
norm_phase_std = [7.5044	6.526	6.714	8.0519	4.5136	7.6279	6.2448	11.3167];
%% Plot stimulation half-cycle gain over frequency for each of the conditions
    %This can only be done if there is a RotaryChairResults.mat file
    %See if one already exists
    res_file = extractfield(dir([path,filesep,'*RotaryChairResults.mat']),'name');
    if isempty(res_file)
        disp('First, make the table with the cycle sine fits...')
        MakeCycleSineFitTable(path,Cyc_Path);
        res_file = extractfield(dir([path,filesep,'*RotaryChairResults.mat']),'name');
    end
    load(res_file{end},'all_results')
    all_results(~contains(all_results.Condition,'Sine'),:) = [];
    %Get the stim ear
    sub = str2double(all_results.Subject{end}(4:6)); %must be a string with MVI***
    ear = Ears{sub};
    %Figure out which frequenicies you need
    [freq_ax,indf] = sort(cellfun(@str2double,strrep(unique(all_results.Frequency)','Hz','')));
    freqs = unique(all_results.Frequency);
    freqs = freqs(indf);
    fnum = length(freqs);
    %Figure out which experiments to plot
    conds = [all_results.Goggle,cellstr(datestr(all_results.Date)),split(all_results.Condition,' ')];
    conds(:,any(contains(conds,{'Hz','dps','Sine','LHRH'}))) = [];
    conds = join(conds,{' '});
    exp_name = unique(conds,'stable')';
    enum = length(exp_name);    
    %Create arrays for graphing
    excludedatawithgainbelow=0.015;
    %Make the items to plot
    gain = NaN(enum,fnum);
    gain_sd = NaN(enum,fnum);
    phase = NaN(enum,fnum);
    phase_sd = NaN(enum,fnum);
    for i = 1:enum
        eparts = split(exp_name(i),' ');
        inds = find(contains(all_results.Goggle,eparts{1})&all_results.Date==datetime(eparts{2})&contains(all_results.Condition,eparts{3}));
        [~,f_ind] = ismember(all_results.Frequency(inds),freqs);
        if strcmp(ear,'L')
            [sub_gain,eye] = max([all_results.Gain_LZ_HIGH(inds)';all_results.Gain_RZ_HIGH(inds)']);
            sub_gain_sd = all_results.Gain_LZ_HIGH_sd(inds)';
            sub_gain_sd(eye==2) = all_results.Gain_RZ_HIGH_sd(inds(eye==2));
        elseif strcmp(ear,'R')
            [sub_gain,eye] = max([all_results.Gain_LZ_LOW(inds)';all_results.Gain_RZ_LOW(inds)']);
            sub_gain_sd = all_results.Gain_LZ_LOW_sd(inds)';
            sub_gain_sd(eye==2) = all_results.Gain_RZ_LOW_sd(inds(eye==2));
        end
        sub_phase = all_results.Phase_L(inds);
        sub_phase(eye==2) = all_results.Phase_R(inds(eye==2));
        sub_phase_sd = all_results.Phase_L_sd(inds);
        sub_phase_sd(eye==2) = all_results.Phase_R_sd(inds(eye==2));
        sub_phase(sub_gain < excludedatawithgainbelow) = NaN;
        sub_phase_sd(sub_gain < excludedatawithgainbelow) = NaN;
        %remove SD for n=1
        n = all_results.Cycles(inds);
        sub_gain_sd(n==1) = NaN;
        sub_phase_sd(n==1) = NaN;
        %Put into mats
        gain(i,f_ind) = sub_gain;
        gain_sd(i,f_ind) = sub_gain_sd;
        phase(i,f_ind) = sub_phase;
        phase_sd(i,f_ind) = sub_phase_sd;
    end
    % Plot
    h=gobjects(1,enum+3);
    ha = gobjects(1,2);
    logxshift=1.03; %how much to multiply the x value to offset its marker rightward (divide to move left)
    %markerbig=5;
    %markersmall=4;
    linethick=2;
    linethin=1;
    errorbarcapsize=1;
    graydark=0.85;
    graylight=0.95;
    figsizeinches=[7,6];
    %figsizeinchesBoxplot=[2.3,4];
    plot_offset = [logxshift^-1, logxshift^-2, 1, logxshift, logxshift^2];
    figure('Units','inch','Position',[2 2 figsizeinches],'Color',[1,1,1]);%CDS083119a
    annotation('textbox',[0 0 1 1],'String',[Cyc_Path,newline,code_Path,filesep,...
            'plotRotaryChairVisitSummary.m',newline,...
            'VOGA',version,newline,Experimenter],'FontSize',5,...
        'EdgeColor','none','interpreter','none');
    ha(1) = subplot(2,1,1);
    ha(1).Position = [0.1,0.55,0.85,0.4];
    ha(2) = subplot(2,1,2);
    ha(2).Position = [0.1,0.10,0.85,0.4];
    axes(ha(1))
    hold on
    h(enum+3) = fill([freq,fliplr(freq)],[norm_gain_m+2*norm_gain_std,fliplr(norm_gain_m-2*norm_gain_std)],graylight*[1 1 1],'LineStyle','none');
    h(enum+2) = fill([freq,fliplr(freq)],[norm_gain_m+1*norm_gain_std,fliplr(norm_gain_m-1*norm_gain_std)],graydark*[1 1 1],'LineStyle','none');
    h(enum+1) = plot(freq,norm_gain_m,'k--','LineWidth',linethick);
    for i = 1:enum
        h(i) = plot(plot_offset(i)*freq_ax,gain(i,:),'LineWidth',linethick);
        errorbar(plot_offset(i)*freq_ax,gain(i,:),gain_sd(i,:),'Color',h(i).Color,'LineStyle','none','LineWidth',linethin,'CapSize',errorbarcapsize)
    end
    hold off
    legend(h,[exp_name,{'Normal mean','Normal±1SD','Normal±2SD'}],'Location','NorthWest','NumColumns',2)
    title([all_results.Subject{1},' ',all_results.Visit{1},' Rotary Chair Sinusoids'])
    ylabel('Horizontal VOR Gain')
    set(ha(1),'XTick',freq_ax,'XTickLabel',[])
    axes(ha(2))
    hold on
    fill([freq,fliplr(freq)],[norm_phase_m+2*norm_phase_std,fliplr(norm_phase_m-2*norm_phase_std)],graylight*[1 1 1],'LineStyle','none')
    fill([freq,fliplr(freq)],[norm_phase_m+1*norm_phase_std,fliplr(norm_phase_m-1*norm_phase_std)],graydark*[1 1 1],'LineStyle','none')
    plot(freq,norm_phase_m,'k--','LineWidth',linethick)
    for i = 1:enum
        plot(plot_offset(i)*freq_ax,phase(i,:),'Color',h(i).Color,'LineWidth',linethick);
        errorbar(plot_offset(i)*freq_ax,phase(i,:),phase_sd(i,:),'Color',h(i).Color,'LineStyle','none','LineWidth',linethin,'CapSize',errorbarcapsize)
    end
    %plot([0.01,3],[0,0],'k:')
    hold off
    ylabel('Phase Lead (deg)')
    set(ha(2),'XTick',freq_ax,'XTickLabel',freq_ax)
    xlabel('Frequency [Hz]')
    set(ha,'XLim',[0.9*freq_ax(1),2.2],'YTickLabelMode','auto','XScale','log')
    set(ha(1),'YLim',[0 0.85],'YTick',0.1:0.1:0.8)
    set(ha(2),'YLim',[-30,100],'YTick',-20:20:100)
    savefig([path,filesep,all_results.Subject{1},'-',all_results.Visit{1},'-RotaryChair-Sine-LHRH-OverFreq','.fig'])
end