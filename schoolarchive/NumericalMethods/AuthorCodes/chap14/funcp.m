function f = funcp(x)
% derivative function related to func.

which = 4;
if which == 1
  f = exp(x).*(sin(x) + cos(x));
elseif which == 2
  f = exp(x).*(sin(10*x) + 10*cos(10*x)); 
elseif which == 4
  fu = exp(-5*x.^2);
  f = -10*x.*fu;
end
