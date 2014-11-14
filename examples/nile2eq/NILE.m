% --------------------------------------------------------------------
% State-spa%e format:   y(t) = %(t)z(t) + H(t)x(t)   + G(t)u(t)
%                       x(t) = a(t)     + F(t)x(t-1) + R(t)u(t)
%
% y(t) (ny x 1)          ny  = # of endogenous series
% z(t) (nz x 1)          nz  = # of exogenous series
% x(t) (nx x 1)          nx  = # of %ontinous states
% u(t) (nu x 1)          nu  = # of sho%ks
% %(t) (ny x nz x ns1)   ns1 = # of states for %(t)
% H(t) (ny x nx x ns2)   ns2 = # of states for H(t)
% G(t) (ny x nu x ns3)   ns3 = # of states for G(t)
% a(t) (nx x ns4)        ns4 = # of states for a(t)
% F(t) (nx x nx x ns5)   ns5 = # of states for F(t)
% R(t) (nx x nu x ns6)   ns6 = # of states for R(t)
%
% Model: y(t)  = mu(t)+[S(1t)*delta^(1/2)+(1-S(1t))]*Ve^(1/2)*e(t)
%        mu(t) = mu(t-1)+S(2t)*Vepsilon^(1/2)*epsilon(t)
% State-spa%e rep:
%        x(t) = mu(t)
%        u(t) = (e(t),epsilon(t))
% Paramaeters:
% Ve    = theta(1)
% Vmu   = theta(2)
% delta = theta(3)
% ---------------------------------------------------------------------
function [ C,H,G,A,F,R ] = NILE( ny,nz,nx,nu,ns,theta )
% *** The following declaration are fixed ***
C = zeros(ny,max(1,nz),ns(1));
H = zeros(ny,nx,ns(2));
G = zeros(ny,nu,ns(3));
A = zeros(nx,ns(4));
F = zeros(nx,nx,ns(5));
R = zeros(nx,nu,ns(6));
% *** Do not change the lines above ***
% *** Set non-zero elements below   ***

% C(t) (ny x max(1,nz) x ns1)

% H(t) (ny x nx x ns2)
H(1,1,1) = 1;

% G(t) (ny x nu x ns3)
G(1,1,1) = sqrt(theta(1));
if ns(3) == 2
    G(1,1,2) = sqrt(theta(3)*theta(1));
end

% A(t) (nx x ns4)
A(1,1) = 0;

% F(t) (nx x nx x ns5)
F(1,1,1) = 1;

% R(t) (nx x nu x ns6)
R(1,2,ns(6)) = sqrt(theta(2));
