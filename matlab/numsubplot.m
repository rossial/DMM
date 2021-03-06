function [nr,nc] = numsubplot(n)
% ------------------------------------------------------------------------
% nr number of rows, nc number of cols
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

if n <= 64
    if n == 1
        nr = 1;
        nc = 1;
    elseif n == 2
        nr = 2;
        nc = 1;
    elseif n == 2
        nr = 2;
        nc = 2;
    elseif n == 3
        nr = 3;
        nc = 1;
    elseif n == 4
        nr = 2;
        nc = 2;
    elseif n == 5 | n == 6
        nr = 3;
        nc = 2;
    elseif n == 7 | n == 8 | n == 9
        nr = 3;
        nc = 3;
    elseif n == 10 | n == 11 | n == 12
        nr = 4;
        nc = 3;
    elseif n == 13 | n == 14 | n == 15 | n == 16
        nr = 4;
        nc = 4;
    elseif n == 17 | n == 18 | n == 19 | n == 20
        nr = 5;
        nc = 4;
    elseif n == 21 | n == 22 | n == 23 | n == 24 | n == 25
        nr = 5;
        nc = 5;
    elseif n >25 & n <31
        nr = 6;
        nc = 5;
    elseif n >30 & n <=36
        nr = 6;
        nc = 6;
    elseif n >36 & n <=42
        nr = 7;
        nc = 6;
    else
        nr = 8;
        nc = 8;
    end
    figure
    set(0,'defaultlinelinewidth',2);
    set(0,'defaultaxesfontweight','bold');
    set(0,'defaulttextfontweight','bold');
    set(0,'defaulttextfontsize',12);
else
    disp('n is too large')
end
