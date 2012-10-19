function f = func(x)

which = 1;
if which == 1
  f = 2*sqrt(1-x.^2);
elseif which == 2
  f = 200./(2*x.^3-x.^2) .* (5*sin(20./x)).^2;
elseif which == 3
  f = exp(-x.^2);  
end