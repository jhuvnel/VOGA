function trace_out = irls_filt(trace_in,param)
    if nargin < 2 || length(param)~=1 || isnan(param)
        trace_out = trace_in;
    else
        trace_out = irlssmooth(trace_in,param);
    end
end

