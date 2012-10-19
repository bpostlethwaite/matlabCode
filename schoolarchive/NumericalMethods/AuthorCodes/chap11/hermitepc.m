function y = hermitepc (xi, fi, fip, x)
%
% function y = hermitepc (xi, fi, fip, x)
%
% piecewise cubic hermite interpolation
% evaluate at x Hermite cubic interpolant based on xi, fi , fip

y = zeros(size(x));
np1 = length(xi); n = np1-1;
index = 1;

for i = 1:n
    % the ith subinterval
    h = xi(i+1) - xi(i);
    inxxi = find ( (xi(i) <= x).*(x <= xi(i+1)) ); xxi = x(inxxi);
    if length(xxi) > 0
        t = (xxi-xi(i))/h;
        yyi = fi(i) + h*fip(i)*t + (3*(fi(i+1)-fi(i)) - h*(fip(i+1)+2*fip(i)))*t.^2 ...
             + (h*(fip(i+1)+fip(i))-2*(fi(i+1)-fi(i)))*t.^3;
        y(index:index+length(xxi)-1) = yyi(1:length(xxi));
        index = index + length(xxi);
    end
end