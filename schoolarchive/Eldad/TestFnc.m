clear all
close all

% a = 10000;
% b = 3;
% 
% N = 50;
% x = linspace(0,1,N+1);
% [X,Y] = meshgrid(x,x);
% 
% 
% u = @(t)(atan(a*t - 1/2*a).*t - 1/2*atan(a*t - 1/2*a) ...
%         -1/(2*a)*log((a*t - 1/2*a).^2 + 1) + b*t);
% up = @(t)(atan(a*(t-0.5)) + b);
% upp = @(t)(a./(a^2 * (t - 0.5).^2 + 1) );
% m = @(t,c,x,y)(1./(up(x) + up(y)) .* (c*exp(t)) );   % Allow for anisotropy with c
% J = @(t)(up(t).*m(t));
% 
% 
% mesh(X,Y,m(X+Y,1,X,Y).*(up(X) + up(Y)))


L = 1;
n = 40;
h = L/n;

[xc,yc]   = ndgrid(h/2:h:L-h/2, h/2:h:L-h/2); % Cell centres
[xdx,ydx] = ndgrid(h:h:L-h, h/2:h:L-h/2); % x walls 
[xdy,ydy] = ndgrid(h/2:h:L-h/2, h:h:L-h); % y walls 




a = 0.02;
a2 = 1e9; 
b = 1;

Gd = atan(a2*(xdx + ydx - 1) + b);
x = xdx; y = ydx;
G = (exp( -( ((x - 0.5).^2)./a + ((y - 0.5).^2)./a  ))) .*Gd;
    
   % m = @(t,x,y,z)(1./(up(t,x,y,z))  ); 
mesh(xdx,ydx,G)