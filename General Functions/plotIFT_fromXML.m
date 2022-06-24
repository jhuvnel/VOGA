%% Assumes files are named with MVI# and visit
clear;
clc;
files = extractfield(dir('*.xml'),'name');
fname = files{end};
sub = fname(strfind(fname,'MVI'):(strfind(fname,'-Vis')-1));
vis = fname(strfind(fname,'Visit'):(strfind(fname,'-IFT')-1));
fdata = cellstr(readlines(fname));
exp_start = find(contains(fdata,'IFT'));
exp_end = find(contains(fdata,'/Result'));
impedance_data = [cell2table(repmat([{sub},{vis}],length(exp_start),1),'VariableNames',{'Subject','Visit'}),...
    array2table(NaT(length(exp_start),1),'VariableNames',{'Date'}),...
    array2table(NaN(length(exp_start),11),'VariableNames',{'Current_cu','Duration_us','E3','E4','E5','E6','E7','E8','E9','E10','E11'})];
for i = 1:length(exp_start)
    rel_text = fdata(exp_start(i):exp_end(i));
    date_str = rel_text{1};    
    date = datetime(date_str(strfind(date_str,'CreateDate')+(12:21)),'Format','yyyy-MM-dd');
    impedance_data{i,3} = date;
    curr_str = rel_text{contains(rel_text,'StimulationCurrent')};
    curr = curr_str((strfind(curr_str,'">')+2):(strfind(curr_str,'</')-1));
    impedance_data{i,4} = str2num(curr);
    dur_str = rel_text{contains(rel_text,'PulseDuration')};
    dur = dur_str((strfind(dur_str,'">')+2):(strfind(dur_str,'</')-1));
    impedance_data{i,5} = str2num(dur)*10^6;
    for j = 3:11
        imp_str = rel_text{find(contains(rel_text,['Number="',num2str(j),'"']))+1};
        imp = imp_str((strfind(imp_str,'">')+2):(strfind(imp_str,'</')-1));
        impedance_data{i,j+3} = str2num(imp);
    end
end
save(strrep(fname,'xml','mat'),'impedance_data')
% Make plots with the values
labs = cellfun(@num2str,num2cell(impedance_data{:,4:5}),'UniformOutput',false);
labs = strcat(labs(:,1),{'(cu), '},labs(:,2),{'(us)'});
fig = figure(1);
set(fig,'Color',[1,1,1]);
plot(3:11,impedance_data{:,6:end}/1000,'-o')
set(gca,'YLim',[0 15],'Xlim',[2.5 11.5])
legend(labs)
xlabel('Electrode Number')
ylabel('Impedance (kOhms)')
title([sub,' ',vis,' Maestro Testing'])
savefig(fig,strrep(fname,'xml','fig'))
saveas(fig,strrep(fname,'xml','png'))
close;