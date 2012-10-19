function [x] = tridiagSolv(l,c,r,b)
%% Uses gaussian elimination to solve tridiagonal systems
% Requires equal length vectors a1,a2,a3 and b (length n)
n = length(c); % size of matrix

for k = 1:n-1
    g = l(k)/c(k);
    c(k+1) = c(k+1) - g*r(k);
    b(k+1) = b(k+1) - g*b(k);
end

% Now back substitute.
x(n) = b(n) / c(n);
for k = n-1:-1:1
    x(k) = (b(k) - r(k)*x(k+1))/c(k);
end
