function f = func(x)
% function to be differentiated.

which = 4;
if which == 1
  f = exp(x).*sin(x);
elseif which == 2 % for simple 2pbvp
  f = exp(x).*sin(10*x);
elseif which == 4 % fft should work on this
    f = exp(-5*x.^2);
end
