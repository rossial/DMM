function [prc,iqr] = myprc(Y,NHIST,m,M)
% ------------------------------------------------------------------------
% Percentiles at m and M /% and iqr
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

[N,X] = hist(Y,NHIST);
NC = cumsum(N)/length(Y);
aa = find(NC<m);
if isempty(aa) == 1
    prc(1,1) = min(X);
else
    prc(1,1) = X(max(aa));
end
aa = find(NC>M);
if isempty(aa) == 1
    prc(2,1) = max(X);
else
    prc(2,1) = X(min(aa));
end
iqr = X(min(find(NC>.75)))-X(max(find(NC<.25)));
