function [xc, yc, rho] = chebl (f,a,b,n)
%
% given a function f on an interval [a,b], find interpolating
% Chebyshev polynomial of degree n. 
% Use Lagrange form.

% Construct abscissae and data
ni = 0:1:n;
ti = cos(pi/(2*(n+1)).* (2*ni+1));
xc = a + (b-a)/2*(ti+1);
yc = f(xc);

% Coefficients for Lagrange form
rho = ones(n+1,1);
for j=1:n+1
    for i=1:n+1
        if i ~= j
            rho(j) = rho(j) * (xc(j)-xc(i));
        end
    end
end