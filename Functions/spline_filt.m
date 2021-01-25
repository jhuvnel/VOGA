function trace_out = spline_filt(t_in,trace_in,t_out,param)
    if nargin < 4 || length(param)~=1
        trace_out = trace_in;
    else
        %Spline filter and resample
        t_in = t_in(~isnan(trace_in));
        trace_in(isnan(trace_in)) = [];
        pp = csaps(t_in,trace_in,param);
        nums = fnplt(pp);
        q = [true,diff(nums(1,:))~=0]; 
        trace_out = spline(nums(1,q),nums(2,q),t_out);
    end
end

