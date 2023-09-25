function filt_tab = MakeCycAvg__filterParameterCheck(filt_tab)
for i = 1:size(filt_tab,2)
    type = filt_tab.Properties.VariableNames{i};
    param_out = filt_tab.(type);
    switch type
        case {'median','mean'}
            param_out = abs(floor(param_out)); %positive integer
            param_out(mod(param_out,2)==0) = param_out(mod(param_out,2)==0)+1; %make odd
            param_out(param_out<2) = NaN; %rm if 0 or 1
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
        case 'sgolay'
            param_out = abs(floor(param_out)); %odd integers greater than 3 (order of 2)
            param_out(mod(param_out,2)==0) = param_out(mod(param_out,2)==0)+1; %make odd
            param_out(param_out<4) = NaN; %remove if 0 or 3
    end
    filt_tab.(type) = param_out;
end
end