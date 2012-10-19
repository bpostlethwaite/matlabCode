function [q,u] = funrhs(xx,yy,t)
%
% function [q,u] = funrhs(xx,yy,t)
%
% if isw = 1 then fig 9.3 from ascherbook
%
% Generate rhs and u frm known solution u
% u = sin pi* x (1 - exp(-y)) exp (-t) on unit square
% PDE u_t = (au_x)_x + (au_y)_y + q
%
% if isw = 2 then fig 16.13 from ascher-greif

isw = 2;

Jx = length(xx) - 2;
Jxp2 = Jx+2;
Jy = length(yy) - 2;
Jyp2 = Jy+2;

if isw == 1
    
  u = sin(pi*xx) * (1 - exp(-yy)) *  exp(-t);

  a = 1 + xx*ones(1,Jyp2) + ones(Jxp2,1)*yy;
  ux = pi*cos(pi*xx) * (1 - exp(-yy)) * exp(-t);
  uy = sin(pi*xx) * exp(-yy)  * exp(-t);
  
  auxx = ux - pi^2*a .* u;
  auyy = uy - a .* uy;

  q = -(u + auxx + auyy);
  
elseif isw == 2
    
    u = 25 * ones(Jxp2,Jxp2);
    q = zeros(Jxp2,Jxp2);
    
end
