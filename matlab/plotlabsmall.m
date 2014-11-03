% Create labels for plotting
% By A.Rossi, 2001
% Usage [lab] =plotlab(annoi,stperi,nobs,freq);
% INPUT:  annoi  = starting year;
%         stperi = starting period;
%         nobs   = # of obs;
%         freq   = frequency (1 4 12).
function [lab] =plotlabsmall(annoi,stperi,nobs,freq)
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
