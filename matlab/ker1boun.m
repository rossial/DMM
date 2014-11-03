% Usage: function [x,fx] = ker1boun(serie,h,step,minb,maxb)
% A. Rossi, November 2014
function [x,fx] = ker1boun(serie,h,step,minb,maxb)
ll  = length(serie);
llh = ll*h;
x   = (minb:step:maxb)';
for ii=1:length(x)
    g = (1/sqrt(2*pi))*exp(-.5*(((serie-x(ii))/h).^2));
    fx(ii,1) = sum(g)/llh;
end