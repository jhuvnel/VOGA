function trace_out = accel_QPR(t_in,trace_in,param)
    if nargin < 3 || length(param)~=1 || isnan(param) || all(isnan(trace_in))
        trace_out = trace_in;
    else
        %Find parts of the traces above/below the threshold and linearly
        %interpolate over that
        trace_out = trace_in;
        over_thresh = [abs(diff(reshape(trace_in,[],1))/median(diff(t_in)))>abs(param);false]&[false;abs(diff(reshape(trace_in,[],1))/median(diff(t_in)))>abs(param)];
        trace_out(over_thresh) = NaN;
        nan_locs = find(isnan(trace_out));
        small_gaps = find(diff(nan_locs)>1&diff(nan_locs)<10);
        for i = 1:length(small_gaps)
            trace_out(nan_locs(small_gaps(i)):nan_locs(small_gaps(i)+1)) = NaN;
        end
        t_nonan = t_in(~isnan(trace_out));
        if ~isempty(t_nonan)
            trace_out = interp1(t_nonan,trace_out(~isnan(trace_out)),t_in);
        else
            trace_out = trace_in;
        end
        %plot(t_in,trace_in,'k.',t_in,trace_out,'b*')
    end
end

