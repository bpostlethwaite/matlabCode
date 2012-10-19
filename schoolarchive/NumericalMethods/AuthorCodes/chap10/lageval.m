function p = lageval (x, xx, yy, rho)
%
%  function newteval (x, xx, table)
%
%  evaluate value of interpolating polynomial at x
%  based on interpolation points xx and divided
%  difference table  table.

np1 = length(xx); 
psi = ones(size(x));

for i=1:np1
    psi = psi.*(x-xx(i));
end
    
p = zeros(size(x));
for j=1:np1
  p = p + yy(j)./ ((x-xx(j))*rho(j));
end
p = p .* psi;