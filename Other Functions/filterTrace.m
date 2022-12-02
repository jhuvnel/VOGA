function trace_out = filterTrace(type,trace_in,param,t_in,t_out)
%Add more filter types as needed
switch type
    case 'lowpass'
        if nargin < 4 || length(param)~=1 || isnan(param) || all(isnan(trace_in))
            trace_out = trace_in;
        else
            trace_out = lowpass(trace_in,param,1/median(diff(t_in)));
        end        
    case 'spline'
        if nargin < 5 || length(param)~=1 || isnan(param) || all(isnan(trace_in))
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
    case 'median'
        if nargin < 3 || length(param)~=1 || isnan(param) || all(isnan(trace_in))
            trace_out = trace_in;
        else
            trace_out = medfilt1(trace_in,param,'omitnan');
        end
    case 'sgolay'
        if nargin < 3 || length(param)~=2 || any(isnan(param)) || all(isnan(trace_in))
            trace_out = trace_in;
        else
            trace_out = sgolayfilt(trace_in,param(1),param(2));
        end
    case 'irls'
        if nargin < 3 || length(param)~=1 || isnan(param) || all(isnan(trace_in))
            trace_out = trace_in;
        else
            trace_out = irlssmooth(trace_in,param);
        end
    case 'accel'
        if nargin < 4 || length(param)~=1 || isnan(param) || all(isnan(trace_in))
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
        end
    case 'manual_interp'
        if nargin < 5 || all(isnan(trace_in))
            trace_out = trace_in;
        else
            %Set the t_in values to NaN over the range of t_out and
            %interpolate
            trace_out = trace_in;
            trace_out(t_in) = NaN; %Now remove trace if needed
            trace_out = interp1(t_out(~isnan(trace_out)),trace_out(~isnan(trace_out)),t_out);
            trace_out(isnan(trace_out)) = 0; %trailing zeros
        end
end
end