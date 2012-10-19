% Example 16.7 -- Table 16.2 

clear all
format short e
h = [.2,.1,.05,.02,.01,.005,.002]; nh = length(h);

y0 = 1; t0 = 1; tf = 10; ye = 1/tf;

% forward Euler
y(1) = y0; t = t0;
for j = 1:nh
  k = h(j); N = (tf-t0) / k;
  for i = 1:N, y(i+1) = y(i) - k*y(i)^2; end
  err1(j) = abs(ye - y(N+1));
end
rate1 = log2(err1(2:nh)./err1(1:nh-1)) ./ log2(h(2:nh)./h(1:nh-1));
err1
rate1

% explicit midpoint
y(1) = y0; t = t0;
for j = 1:nh
  k = h(j); N = (tf-t0) / k;
  for i = 1:N, ym =  y(i) - .5*k*y(i)^2; y(i+1) = y(i) - k*ym^2; end
  err2(j) = abs(ye - y(N+1));
end
rate2 = log2(err2(2:nh)./err2(1:nh-1)) ./ log2(h(2:nh)./h(1:nh-1));
err2
rate2

% rk4
y(1) = y0; t = t0;
for j = 1:nh
  k = h(j); N = (tf-t0) / k;
  [tt,y] =  rk4(@func,[t0,tf],y0,k);
  err4(j) = abs(ye - y(end));
end
rate4 = log2(err4(2:nh)./err4(1:nh-1)) ./ log2(h(2:nh)./h(1:nh-1));
err4
rate4
  
