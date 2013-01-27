function yy = dft2e(xx,a0,a,b)
%
%  function yy = dft2e(xx,a0,a,b)
%
%  Given Fourier coefficients, evaluate triginometric polynomial at xx.

l = size(a,2);
yy = a0/2 * ones(size(xx)); 
for k=1:l-1
  yy = yy + a(k)*cos(k*xx) + b(k)*sin(k*xx);
end
yy = yy + a(l)/2 * cos(l*xx);



 