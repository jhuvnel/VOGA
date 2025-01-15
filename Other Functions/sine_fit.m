%% sine_fit 
%This is a function that is used with fminsearch to find parameters that give
%an optimal fit.
%This function creates a sine wave with different half-cycle amplitudes and phase with offset. 
%It can accommadate multiple frequencies. For n frequencies, p will need to be (3*n+1) x 1.
function vals = sine_fit(t,freq,p)
    p0 = reshape(p(1:end-1),3,length(freq)); %reshape for each frequency
    all_vals = zeros(length(freq),length(t));
    for i = 1:length(freq)
        t0 = sin(2*pi*freq(i)*t+p0(3,i)*pi/180)>=0; %pos half-cycle
        all_vals(i,t0) = p0(1,i)*sin(2*pi*freq(i)*t(t0)+p0(3,i)*pi/180);
        all_vals(i,~t0) = p0(2,i)*sin(2*pi*freq(i)*t(~t0)+p0(3,i)*pi/180+pi);
    end
    vals = sum(all_vals,1)+p(end);
end