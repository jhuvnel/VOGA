%% sine_fit 
%This function creates a sine wave with different amplitudes for each
%half-cycle. It also works for sums of sine waves (which each have their
%own amplitude and phase parameters).
%This is a function that is used in fminsearch to find parameters that give
%an optimal fit.
%For n frequencies, p will need to be 3*n x 1.

function vals = sine_fit(t,freq,p)
    all_vals = zeros(length(freq),length(t));
    for i = 1:length(freq)
        pos = p(3*i-2)*sin(2*pi*freq(i)*t+p(3*i)*pi/180+pi);
        pos(pos > 0) = 0; %positive half-cycle
        neg = p(3*i-1)*sin(2*pi*freq(i)*t+p(3*i)*pi/180+pi);
        neg(neg < 0) = 0; %negative half-cycle
        all_vals(i,:) = pos+neg;
    end
    vals = sum(all_vals,1);
end