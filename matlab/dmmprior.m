function y = dmmprior(x,a,b,tipo)
% ------------------------------------------------------------------------
% tipo: NT normal, BE beta, IG inverted gamma 2 pdf
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

  y = zeros(size(x));
if nargin ~= 4,
    error('dmmprior requires 4 input arguments');
end
if tipo == 'NT'
    if (b<=0)
        error('Variance b must be positive');
    end
    y = (2*pi*b)^(-.5).*exp(-.5*(x-a).^2/b);
elseif tipo == 'BE'
    if any(any((a<=0)|(b<=0)))
        error('Parameter a and b must be positive');
    end
    I = find((x<0)|(x>1));
    y = x.^(a-1).*(1-x).^(b-1)./beta(a,b);
    y(I) = 0;
elseif tipo == 'IG'
    if any(any((a<=0)|(b<=0)))
        error('Parameter a and b must be positive');
    end
    I = find(x>0);
    y(I) = (-.5*b - 1).*log(x(I))-a./(2*x(I))-gammaln(.5*b)-.5*b.*log(2./a);
    y(I) = exp(y(I));
end
