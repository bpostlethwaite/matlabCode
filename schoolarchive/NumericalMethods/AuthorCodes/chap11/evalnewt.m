function p = evalnewt (x, xi, coef)
%
%  function p = evalnewt (x, xi, coef)
%
%  evaluate at x interpolating polynomial in Newton form
%  based on interpolation points xi and coefficients coef

np1 = length(xi);
p = coef(np1)*ones(size(x));
for j=np1-1:-1:1
  p = p.*(x - xi(j)) + coef(j);
end
