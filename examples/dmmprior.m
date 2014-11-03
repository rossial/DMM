function y = dmmprior(x,a,b,tipo)
% tipo: NT normal, BE beta, IG inverted gamma 2 pdf 
% A.Rossi (November 2014)
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