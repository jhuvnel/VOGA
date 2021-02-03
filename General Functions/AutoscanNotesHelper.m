%% Paste current values here (as rounded whole numbers)
currs = [552
576
599
623
651
680
699];
%% Adjust electrode names and other parameters
elec = [repmat({'"RAE11'},7,1)];%;...
    %repmat({'"RPE4'},7,1);....
    %repmat({'"RPE5'},7,1)];
burstpulse = repmat({'200pps'},length(elec),1);
phasedur = repmat({'50us'},length(elec),1);
%% Put it all together
str_curr = strcat(cellstr(num2str(currs)),'uA"');
final = strcat(elec,'-',burstpulse,'-',phasedur,'-',str_curr);
clc;
fprintf('%s\n', final{:})