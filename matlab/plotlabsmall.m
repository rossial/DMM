function [lab] =plotlabsmall(annoi,stperi,nobs,freq)
% ------------------------------------------------------------------------
% Create labels for plotting
% Usage [lab] =plotlab(annoi,stperi,nobs,freq);
% INPUT:  annoi  = starting year;
%         stperi = starting period;
%         nobs   = # of obs;
%         freq   = frequency (1 4 12).
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

anno=annoi;
stper=stperi;
lab=cell(nobs,1);
switch freq
case 1
    for i=1:nobs
        appo=num2str(anno);
        lab(i)=cellstr(appo(3:end));
        anno =anno+1;
    end
case 4
    for i=1:nobs
        appo=num2str(anno);
        appo1=strcat(appo(1:4),'-',num2str(stper));
        lab(i)=cellstr(appo1(3:end));
        if stper < 4
            stper=stper+1;
        else
            stper=1;
            anno=anno+1;
        end
    end
case 12
    for i=1:nobs
        appo=num2str(anno);
        if stper < 10
            appo1=strcat(appo(1:4),'-0',num2str(stper));
        else
            appo1=strcat(appo(1:4),'-',num2str(stper));
        end
        lab(i)=cellstr(appo1(3:end));
        if stper < 12
            stper=stper+1;
        else
            stper=1;
            anno=anno+1;
        end
    end
end
