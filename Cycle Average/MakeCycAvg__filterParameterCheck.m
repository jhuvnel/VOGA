function param_out = MakeCycAvg__filterParameterCheck(param_in,type)
    param_out = param_in;
    switch type
        case {'Median','median','med','Med'}
            param_out = abs(floor(param_out)); %positive integer
            param_out(mod(param_out,2)==0) = param_out(mod(param_out,2)==0)+1; %make out
            param_out(param_out==0) = NaN; %rm if 0
        case {'Spline','spline'}
            param_out(param_out>1) = 1;
            param_out(param_out<0) = 0;   
            param_out(param_out==1) = NaN;
        case {'Irrlsmooth','smooth','Smooth'}
            param_out = abs(floor(param_out));
            param_out(param_out==0) = NaN;
        case {'sgolay'}
            if size(param_out,2) ~=2
                error('Two parameters needed for Savitzky-Golay')
            end
            param_out = abs(floor(param_out)); %positive integers
            param_out(param_out==0,:) = NaN; %remove 0's 
            param_out(mod(param_out(:,2),2)==0,2) = param_out(mod(param_out(:,2),2)==0,2)+1; %p2 must be odd
            param_out(param_out(:,2)==param_out(:,1),2) = param_out(param_out(:,2)==param_out(:,1),1)+1; %If p1 == p2 p2-> p2+1
            param_out(param_out(:,2)<param_out(:,1),:) = fliplr(param_out(param_out(:,2)<param_out(:,1),:)); %Make p2 > p1
         case {'AccelerationThreshold','acc','accel','vel','Vel'}
            param_out = abs(param_in);
            param_out(param_out==0) = NaN;
    end
end