function filt_tab = MakeCycAvg__filterParameterCheck(filt_tab)
for i = 1:size(filt_tab,2)
    type = filt_tab.Properties.VariableNames{i};
    param_out = filt_tab.(type);
    switch type
        case 'median'
            param_out = abs(floor(param_out)); %positive integer
            param_out(mod(param_out,2)==0) = param_out(mod(param_out,2)==0)+1; %make odd
            param_out(param_out==0) = NaN; %rm if 0
       case {'accel','lowpass'}
            param_out = abs(param_out);
            param_out(param_out==0) = NaN;     
        case 'spline'
            param_out(param_out>1) = 1;
            param_out(param_out<0) = 0;
            param_out(param_out==1) = NaN;
        case 'irlssmooth'
            param_out = abs(floor(param_out));
            param_out(param_out==0) = NaN;
        case 'sgolay1'
            param_out = abs(floor(param_out)); %positive integers
            param_out(param_out==0) = NaN; %remove 0's
        case 'sgolay2'
            param_out = abs(floor(param_out)); %positive integers
            param_out(param_out==0|isnan(filt_tab.sgolay1),:) = NaN; %remove 0's
            param_out(param_out<filt_tab.sgolay1) = filt_tab.sgolay1(param_out<filt_tab.sgolay1)+1;
            param_out(mod(param_out,2)==0) = param_out(mod(param_out,2)==0)+1; %Must be odd and bigger than sgolay1
    end
    filt_tab.(type) = param_out;
end
end