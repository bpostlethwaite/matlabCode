function [gg, dg] = gnfun(x,t)
% function and gradient for Example 9.8

t = t(:);

gg = x(1) * exp(x(2)*t) .* cos(x(3)*t);

g1 = exp(x(2)*t) .* cos(x(3)*t);
g2 = x(1) * t .* exp(x(2)*t) .* cos(x(3)*t);
g3 = -x(1) * t .* exp(x(2)*t) .* sin(x(3)*t);
dg = [g1,g2,g3];
  
  
    