% Examples 15.18--15.19 : find area of unit circle (=\pi)

% method 1:  I = int_{-1}^1 dx [int_c^d dy] 
%        c(x) = -sqrt(1-x^2), d(x) = sqrt(1-x^2).

ll = @(x) -sqrt(1-x.^2);
uu = @(x) sqrt(1-x.^2);
fun1d = @(x) 2*sqrt(1-x.^2);
I1 = quad(fun1d,-1,1)
err1 = abs(I1-pi)

% method 2: Monte Carlo

limit = 1000000; sum = 0;
for i = 1:limit
  x = 2*rand(1,2) - 1;
  if x(1)^2+x(2)^2 < 1
      sum = sum + 1;
  end
end
I2 = 4* sum / limit
err2 = abs(I2-pi)

% method 3: just to check use 1D integration in this special case.

[I3,mesh,fevals] = quads(-1, 1, 1.e-6);
fevals
err3 = abs(I3-pi)