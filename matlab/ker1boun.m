function [x,fx] = ker1boun(serie,h,step,minb,maxb)
% ------------------------------------------------------------------------
% Usage: function [x,fx] = ker1boun(serie,h,step,minb,maxb)
%
% Copyright (C) 2010-2014 European Commission
%
% This file is part of Program DMM
%
% DMM is free software developed at the Joint Research Centre of the
% European Commission: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% DMM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with DMM.  If not, see <http://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------

ll  = length(serie);
llh = ll*h;
x   = (minb:step:maxb)';
for ii=1:length(x)
    g = (1/sqrt(2*pi))*exp(-.5*(((serie-x(ii))/h).^2));
    fx(ii,1) = sum(g)/llh;
end
