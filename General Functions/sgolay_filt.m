function trace_out = sgolay_filt(trace_in,param)
    if nargin < 2 || length(param)~=2
        trace_out = trace_in;
    else
        trace_out = sgolayfilt(trace_in,param(1),param(2));
    end
end

