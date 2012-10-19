% experiment with myfft

isw = 2;
if isw == 1    
   y = [0,pi/2,pi,pi/2];
elseif isw == 2
  l = 16; % input parameter
  n = 2*l-1; pi2 = 2*pi;
  x = 0 : pi/l : n*pi/l; % abscissae on [0,pi2]
  t = 3/pi2*x - 1;             % abscissae on [-1,2]
  y = t.^2 .* (t+1).^2 .* (t-2).^2 - ...
      exp(-t.^2) .* sin(t+1).^2 .* sin(t-2).^2;  
end

c = myfft(y);

% define coefficients for dft2e
m = length(y);
l = m/2;
a0 = real(c(1))/l; 
a = real((c(2:l)+c(m:-1:l+2)))/m;
b = real(i*(c(2:l)-c(m:-1:l+2)))/m;
a(l) = real(c(l+1))/l;

if isw == 2
  % interpolate on fine mesh, plot and find max interpolation error
  xx = 0:.01*pi:pi2; tt = 3/pi2*xx - 1;
  yexact = tt.^2 .* (tt+1).^2 .* (tt-2).^2 - ...
           exp(-tt.^2) .* sin(tt+1).^2 .* sin(tt-2).^2;
  yapprox = dft2e(xx,a0,a,b);
 
  % display results
  plot(tt,yexact,'b--',tt,yapprox)
  axis([-1 2 -2 6])
  legend('exact','approx')
  xlabel('t')
  err_max = max(abs(yexact-yapprox)./(abs(yexact)+1))
end

