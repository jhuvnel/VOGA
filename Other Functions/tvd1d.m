function y = tvd1d(x,lambda)
if nargin<2
    lambda = 1;
end
[na,nb] = size(x);
x = reshape(x,[],1);
y = x;
repval = @(y,yn,n) [y(1:(n-1));yn;y((n+1):length(y))];
costfun = @(y) sum((x-y).^2,'omitnan')/length(x)+lambda*sum(abs(diff(y)),'omitnan');
opts = optimset('Display','off');
for i = 1:length(x) 
    y(i) = fminsearch(@(yn) costfun(repval(y,yn,i)),y(i),opts);
end 
y = reshape(y,na,nb);
end