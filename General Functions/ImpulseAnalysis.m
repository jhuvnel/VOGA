%% Create struct with values of interest
list = extractfield(dir([cd,filesep,'Cycle Averages',filesep,'CycAvg_*.mat']),'name');
[indx,tf] = nmlistdlg('PromptString','Select CycAvg files to analyze:',...
                       'SelectionMode','multiple',...
                       'ListSize',[500 600],...
                       'ListString',list); 
if tf                   
    files = list(indx);
    t = 0:0.0001:0.75;
    n = 20; %take first 20
    h_vel = NaN(length(files),n);
    h_wid = NaN(length(files),n);
    h_AUC = NaN(length(files),n); %trapezoidal sum
    e_lat = NaN(length(files),n);
    e_AUC = NaN(length(files),n);
    h_cyc = NaN(length(t),n,length(files));
    e_cyc = NaN(length(t),n,length(files));
    for f = 1:length(files)
        %Assumes a CycAvg alredy exists from using VOGA
        load([cd,filesep,'Cycle Averages',filesep,files{f}],'CycAvg')    
        for i = 1:n %Run per cycle
            if i<=size(CycAvg.head_cyc,2) %if less than 20 cyc, leave as NaN
                %define start and end as 10dps
                head = spline(CycAvg.t,CycAvg.head_cyc(:,i),t);
                eye = spline(CycAvg.t,CycAvg.eye_cyc(:,i),t);
                %head max velocity, head start time, head end time, head impulse width
                [h_vel(f,i),head_max] = max(head);
                h_start = find(head<10&t<t(head_max),1,'last');
                h_end = find(head<10&t>t(head_max),1,'first');
                h_wid(f,i) = 1000*(t(h_end)-t(h_start)); %in ms
                h_AUC(f,i) = trapz(head(h_start:h_end))*median(diff(t));
                h_cyc(:,i,f) = head;
                %eye start time, eye latency, AUC eye
                e_start = find(eye>10&t>=t(h_start),1,'first');
                e_lat(f,i) = 1000*(t(e_start)-t(h_start)); %in ms
                e_AUC(f,i) = trapz(eye(h_start:h_end))*median(diff(t)); 
                e_cyc(:,i,f) = eye;
            end
        end
    end
    all_results.files = files;
    all_results.head_maxvel = h_vel;
    all_results.head_impulsewidth = h_wid;
    all_results.eye_gain = e_AUC./h_AUC;
    all_results.eye_latency = e_lat;
    all_results.head_cyc = h_cyc;
    all_results.eye_cyc = e_cyc;
    save([cd,filesep,datestr(now,'yyyymmdd_HHMMSS'),'_','ImpulseResults.mat'],'all_results')
end
%% Load table in question
res_file = extractfield(dir([cd,filesep,'*Results.mat']),'name')';
if isempty(res_file)
    error('No table with cycle parameters found on this path.')
end
load(res_file{end},'all_results') 
files = all_results.files;
notes = [files,files];
for i = 1:length(files)
    parts = strrep(strrep(split(files{i},'-'),'.mat',''),'_',' ');
    notes{i,1} = parts{6};
    notes{i,2} = [parts{4},' ',parts{end}];
end
%% Plots for each canal
canals = unique(notes(:,1),'stable');
for c = 1:length(canals)
    rel_i = contains(notes(:,1),canals(c));
    rel_notes = notes(rel_i,2);
    rel_head_vel = all_results.head_maxvel(rel_i,:);
    rel_head_wid = all_results.head_impulsewidth(rel_i,:);
    rel_eye_gain = all_results.eye_gain(rel_i,:);
    rel_eye_lat = all_results.eye_latency(rel_i,:); 
    spread = 0.5*rand(size(rel_head_wid,2),1)-0.25;
    x = repmat(1:length(rel_notes),length(spread),1)+spread;    
    %Head metrics
    figure;
    ha(1) = subplot(2,1,1);
    boxplot(rel_head_vel')
    hold on
    plot(x,rel_head_vel','k*')
    hold off
    title([canals{c},' Head Motion Metrics'],'FontSize',15)
    ylabel('Maximum Head Velocity (dps)')
    set(gca,'XTickLabel',[])
    ha(1).Position = [0.13,0.56,0.85,0.39];
    ha(2) = subplot(2,1,2);
    boxplot(rel_head_wid')
    hold on
    plot(x,rel_head_wid','k*')
    hold off
    ylabel('Head Impulse Width (ms)')
    set(gca,'XTickLabel',rel_notes,'XTickLabelRotation',15)
    ha(2).Position = [0.13,0.15,0.85,0.39];
    savefig([canals{c},' Head Motion Metrics.fig'])
    close;
    %Eye metrics
    figure;
    ha(1) = subplot(2,1,1);
    boxplot(rel_eye_gain')
    hold on
    plot(x,rel_eye_gain','k*')
    hold off
    title([canals{c},' VOR Metrics'],'FontSize',15)
    ylabel('Gain (Eye AUC/Head AUC)')
    set(gca,'XTickLabel',[])
    ha(1).Position = [0.13,0.56,0.85,0.39];
    ha(2) = subplot(2,1,2);
    boxplot(rel_eye_lat')
    hold on
    plot(x,rel_eye_lat','k*')
    hold off
    ylabel('Eye Latency (ms)')
    set(gca,'XTickLabel',rel_notes,'XTickLabelRotation',15)
    ha(2).Position = [0.13,0.15,0.85,0.39];
    savefig([canals{c},' VOR Metrics.fig'])
    close;
end
%% Plots with all exps
spread = 0.5*rand(size(all_results.head_maxvel,2),1)-0.25;
x = repmat(1:length(notes),length(spread),1)+spread;
%Head metrics
figure;
boxplot(all_results.head_maxvel')
hold on
plot(x,all_results.head_maxvel','k*')
hold off
title('Maximum Head Velocity During Impulse','FontSize',15)
ylabel('Head Velocity (dps)')
set(gca,'XTickLabel',join(notes,' '),'XTickLabelRotation',30)
savefig('Maximum Head Velocity.fig')
close;
figure;
boxplot(all_results.head_impulsewidth')
hold on
plot(x,all_results.head_impulsewidth','k*')
hold off
title('Head Impulse Width','FontSize',15)
ylabel('Time (ms)')
set(gca,'XTickLabel',join(notes,' '),'XTickLabelRotation',30)
savefig('Head Impulse Width.fig')
close;
%Eye metrics
figure;
boxplot(all_results.eye_gain')
hold on
plot(x,all_results.eye_gain','k*')
hold off
title('VOR Gain','FontSize',15)
ylabel('Eye AUC/Head AUC')
set(gca,'XTickLabel',join(notes,' '),'XTickLabelRotation',30)
savefig('VOR Gain.fig')
close;
figure;
boxplot(all_results.eye_latency')
hold on
plot(x,all_results.eye_latency','k*')
hold off
title('VOR Latency','FontSize',15)
ylabel('Time (ms)')
set(gca,'XTickLabel',join(notes,' '),'XTickLabelRotation',30)
savefig('VOR Latency.fig')
close;
%% Custom Plots

%All manual vs. all aHIT (manual hands on head, aHIT pseudorandom)
sub_notes = {'manual','aHIT'};
aHIT_i = repmat(contains(notes(:,2),'aHIT stacked pseudorandom'),1,size(all_results.head_maxvel,2));
man_i = repmat(contains(notes(:,2),'Manual hands on head'),1,size(all_results.head_maxvel,2));

spread = 0.5*rand(sum(sum(man_i)),1)-0.25;
x = repmat(1:length(sub_notes),length(spread),1)+spread;
%Head metrics
figure;
subplot(1,2,1)
boxplot([all_results.head_maxvel(man_i),all_results.head_maxvel(aHIT_i)])
hold on
plot(x,[all_results.head_maxvel(man_i),all_results.head_maxvel(aHIT_i)],'k*')
hold off
title({'Maximum Head Velocity','During Impulse'},'FontSize',12)
ylabel('Head Velocity (dps)')
set(gca,'XTickLabel',sub_notes,'XTickLabelRotation',0)
subplot(1,2,2)
boxplot([all_results.head_impulsewidth(man_i),all_results.head_impulsewidth(aHIT_i)])
hold on
plot(x,[all_results.head_impulsewidth(man_i),all_results.head_impulsewidth(aHIT_i)],'k*')
hold off
title('Head Impulse Width','FontSize',12)
ylabel('Time (ms)')
set(gca,'XTickLabel',sub_notes,'XTickLabelRotation',0)
savefig('Head Motion Metrics aHIT vs Manual.fig')
close;

[~,var_p1] = vartest2(all_results.head_maxvel(man_i),all_results.head_maxvel(aHIT_i));
disp(['Difference in variance of maximum velocity p-val: ',num2str(var_p1)])
[~,var_p2] = vartest2(all_results.head_impulsewidth(man_i),all_results.head_impulsewidth(aHIT_i));
disp(['Difference in variance of impulse width p-val: ',num2str(var_p2)])

%% Combined aHIV vs. Manual CycAvg
sub_notes = {'manual','aHIT'};
aHIT_files = contains(files,'aHIT')&contains(files,'pseudorandom');
manual_files = contains(files,'Manual')&contains(files,'head');
t = 0:0.0001:0.75;
n = 20;
aHIT_cyc_head = reshape(all_results.head_cyc(:,:,aHIT_files),length(t),[]);
aHIT_cyc_eye = reshape(all_results.eye_cyc(:,:,aHIT_files),length(t),[]);
manual_cyc_head = reshape(all_results.head_cyc(:,:,manual_files),length(t),[]);
manual_cyc_eye = reshape(all_results.eye_cyc(:,:,manual_files),length(t),[]);
figure;
subplot(1,2,1)
plot(1000*t,aHIT_cyc_head,'k',1000*t,aHIT_cyc_eye,'b')
ylabel('Velocity (dps)')
xlabel('Time (ms)')
title('aHIT')
axis([0 500 -100 300])
subplot(1,2,2)
plot(1000*t,manual_cyc_head,'k',1000*t,manual_cyc_eye,'b')
hold on
h(1) = plot(NaN,NaN,'k');
h(2) = plot(NaN,NaN,'b');
hold off
legend(h,{'Head','Eye'})
title('Manual')
ylabel('Velocity (dps)')
xlabel('Time (ms)')
axis([0 500 -100 300])
savefig('Combined Head and Eye Traces aHIT vs Manual.fig')
close;
% spread = 0.5*rand(size(manual_cyc_head,2),1)-0.25;
% x = repmat(1:length(sub_notes),length(spread),1)+spread;
% aHIT_head_SSE = sum((aHIT_cyc_head-mean(aHIT_cyc_head,2)).^2)';
% manual_head_SSE = sum((manual_cyc_head-mean(manual_cyc_head,2)).^2)';
% figure;
% boxplot([manual_head_SSE,aHIT_head_SSE])
% hold on
% plot(x,[manual_head_SSE,aHIT_head_SSE],'k*')
% hold off
% set(gca,'XTickLabel',sub_notes,'XTickLabelRotation',0)
% title('Sum of Squared Errors from Mean Trace Fit')
% 
% [~,var_p3] = vartest2(manual_head_SSE,aHIT_head_SSE);
% disp(['Difference in variance of head trace SSE p-val: ',num2str(var_p3)])