% Percentiles at m and M /% and iqr
% A. Rossi, November 2014
function [prc,iqr] = myprc(Y,NHIST,m,M)
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