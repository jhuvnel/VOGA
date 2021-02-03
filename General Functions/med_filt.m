function trace_out = med_filt(trace_in,param)
    if nargin < 2 || length(param)~=1
        trace_out = trace_in;
    else
        trace_out = medfilt1(trace_in,param,'omitnan');
    end
end

