function [coef,table] = divdif (xi, yi)
%
% function [coef,table] = divdif (xi, yi)
%
% Construct a divided difference table based on data points (xi,yi).
% Upon return, the Newton interpolation coefficients are in coef

np1 = length(xi); n = np1-1;
table = zeros(np1,np1); xi = shiftdim(xi); yi = shiftdim(yi);
% construct divided difference table one column at a time 
table(1:np1,1) = yi;
for k = 2:np1
    table(k:np1,k) = (table(k:np1,k-1) - table(k-1:n,k-1)) ./ ...
                     (xi(k:np1) - xi(1:np1-k+1));
end
coef = diag(table);  % the diagonal elements of table
