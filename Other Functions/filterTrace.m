function trace_out = filterTrace(type,trace_in,param,t_out,t_in)
%No filtering needed/possible
if all(isempty(param))||any(isnan(param))||all(isnan(trace_in))
    trace_out = trace_in;
    return;
end
%Error handling: nargin
inval_nargin = nargin<3||nargin<4&&contains(type,{'accel','manual_interp','lowpass'})||...
        nargin<5&&contains(type,{'spline'});
if inval_nargin
    error(['Not enough input arguments to use ',type,' filter.'])
end
%Error handling: parameter number
inval_pnum = contains(type,'sgolay')&&length(param)~=2||...
        contains(type,'manual_interp')&&length(param)<=0||...
        ~contains(type,{'sgolay','manual_interp'})&&length(param)~=1;
if inval_pnum
    error('Incorrect number of filter parameters entered.')
end
%Add more filter types as needed
switch type
    case 'lowpass'
        trace_out = lowpass(trace_in,param,1/median(diff(t_out)));
    case 'spline'
        t_in = t_in(~isnan(trace_in));
        trace_in(isnan(trace_in)) = [];
        pp = csaps(t_in,trace_in,param);
        nums = fnplt(pp);
        q = [true,diff(nums(1,:))~=0];
        trace_out = spline(nums(1,q),nums(2,q),t_out);
    case 'median'
        trace_out = medfilt1(trace_in,param,'omitnan');
    case 'sgolay'
        trace_out = sgolayfilt(trace_in,param(1),param(2));
    case 'irls'
        trace_out = irlssmooth(trace_in,param);
    case 'accel'
        trace_out = trace_in;
        vel_thresh = 30; %Set threshold for saccade
        Fs = 1/median(diff(t_out));
        g_trac = gradient(reshape(trace_in,[],1))*Fs;
        over_thresh = abs(g_trac)>abs(param)&abs(trace_out)>vel_thresh;
        under_thresh = abs(g_trac)<abs(param)&abs(trace_out)<2*vel_thresh;
        o_start = find(diff(over_thresh)>0);
        for i = 1:length(o_start)
            over_thresh(find(under_thresh(1:o_start(i)),1,'last'):find(under_thresh(o_start(i):end),1,'first')-1+o_start(i)) = true;
        end
        trace_out(over_thresh) = NaN;
        t_nonan = t_out(~isnan(trace_out));
        if ~isempty(t_nonan)
            trace_out = interp1(t_nonan,trace_out(~isnan(trace_out)),t_out);
            trace_out(isnan(trace_out)) = 0;
        else
            trace_out = trace_in;
        end
    case 'manual_interp'
        %Set the t_in values to NaN over the range of t_out and interpolate
        trace_out = trace_in;
        trace_out(param) = NaN; %Now remove trace if needed
        trace_out = interp1(t_out(~isnan(trace_out)),trace_out(~isnan(trace_out)),t_out);
        trace_out(isnan(trace_out)) = 0; %trailing zeros
end
end