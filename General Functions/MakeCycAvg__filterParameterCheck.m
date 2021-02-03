function param_out = MakeCycAvg__filterParameterCheck(param_in,type)
    if ~any(isnan(param_in))
        switch type
            case {'Median','median','med','Med'}
                param_in = abs(floor(param_in)); %positive integer
                if param_in==0
                    param_out = [];
                elseif mod(param_in,2)==0 %had to be odd
                    param_out = param_in+1;
                else
                    param_out = param_in;
                end
            case {'Spline','spline'}
                if param_in>1||param_in<0
                    param_out = [];
                else
                    param_out = param_in;
                end
            case {'Irrlsmooth','smooth','Smooth'}
                param_out = abs(floor(param_in));
            case {'AccelerationThreshold','acc','accel','vel','Vel'}
                param_out = abs(param_in);
            case {'sgolay'}
                param_in = abs(floor(param_in)); %positive integers
                if param_in(1)==param_in(2)||any(param_in==0)
                    param_out = []; %illegal case
                    return;
                elseif param_in(1) > param_in(2) %frame len (2nd param) must be bigger than order (1st param), switch if needed
                    param_in = fliplr(param_in);
                end
                if mod(param_in(2),2)==0 %framelen has to be odd
                    param_in(2) = param_in(2)+1;
                end
                param_out = param_in;
        end
    else
        param_out = [];
    end
end