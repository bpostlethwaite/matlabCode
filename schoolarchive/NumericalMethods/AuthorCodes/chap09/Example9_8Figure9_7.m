% Example 9.8 -- Figure 9.7 : data fitting using Gauss-Newton

clear all
close all
tol = 1.e-6; itmax = 20;
t = 0:.005:1;
xe = [1,-2,20]';
ge = xe(1)*exp(xe(2)*t) .* cos(xe(3)*t);

% generate data
td = 0:.02:1;
be = xe(1)*exp(xe(2)*td) .* cos(xe(3)*td); % exact values
noise = randn(1,length(be));
b = (1 + .2*noise).*be; b = b(:);  % polluted by random noise

% initial guess
%x = ones(3,1); 
x = [1.2,-1.9,18]';
[gg, dg] = gnfun(x,td);
r = b - gg;
count = 0;

while (norm(r) > tol)* (count < itmax)
  p = dg \ r;
  x = x + p;
  [gg, dg] = gnfun(x,td);
  r = b - gg;
  count = count + 1;
end

x
xx = x(1)*exp(x(2)*td) .* cos(x(3)*td);
plot (t,ge,td,b','go',td,xx,'rd')




