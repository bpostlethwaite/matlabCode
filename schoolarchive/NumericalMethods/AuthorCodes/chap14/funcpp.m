function f = funcpp(x)
% 2nd derivative function related to func.

which = 4;
if which == 1
  f = exp(x).*(2*cos(x));
elseif which == 2
  f = -exp(x).*(99*sin(10*x)-20*cos(10*x));
elseif which == 4
  fuu = exp(-5*x.^2);
  f = -10*fuu + 100*x.^2 .*fuu;
end
