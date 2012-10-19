function y = roundc(x,n)
%
% function y = roundc(x,n)
%
% Returns x rounded to n decimal digits
%         (take n < number of significant digits in x )
% Note x can be an array; y would have the same size
   
f = zeros(size(x)); e = f;     % initialize fractions, exponents
nz = find(x ~= 0);             % find nonzero components of x
xx = abs(x(nz));
e(nz) = ceil(log10(xx));       % decimal exponent of nonzero x
f(nz) = xx ./ 10.^e(nz);       % fraction of nonzero x
   
s = 1;
frac = zeros(size(x));
for j = 1:n
  s = s*10;
  d = floor(f*s + 0.5);        % extract digit
  f = f - d/s;
  frac = frac + d/s;           % add digit to rounded fraction
end

y = sign(x) .* (frac .* 10.^e);

