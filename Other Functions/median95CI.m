%% median95CI 
% Allows the ability to query items needed for calculation of the 95%CI of the
% median using the method described in *******
% % Inputs
% vals is a 1 x n or n x 1 vector of values
% ouptut allows the user to select what item they want
% 'lowCI'
% 'highCI'
% 'n' - number of non NaN values
% 'r' - lower index of sorted values used for 95%CI 
% 's' - upper index of sorted values used for 95%CI
function out_num = median95CI(vals,output)
vals = reshape(sort(vals(~isnan(vals))),[],1);
out.n = length(vals);
out.r = round(out.n/2-1.96*sqrt(out.n)/2);
out.s = round(1+out.n/2+1.96*sqrt(out.n)/2);
if out.n > 5 % Can be calculated
    out.lowCI = vals(out.r);
    out.highCI = vals(out.s); 
else % n too small for calculations   
    out.lowCI = NaN;
    out.highCI = NaN; 
end
out_num = out.(output);
end