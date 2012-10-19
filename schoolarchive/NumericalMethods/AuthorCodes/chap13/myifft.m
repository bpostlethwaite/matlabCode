function y = myifft(c)
%
%  Apply IFFT, assuming the length m of c is a power of 2 and even.
%  

m = length(c);
if m == 2
    y(1) = .5*(c(1)+c(2)); y(2) = .5*(c(1)-c(2));
else
    % recursion step
    l = m/2;
    omegam = exp(i*2*pi/m); omegap = omegam.^(0:1:l-1);
    ceven = c(1:2:m-1); codd = c(2:2:m);
    yeven = myifft(ceven); yodd = myifft(codd);
    y(1:l) = .5*( yeven + omegap .* yodd);
    y(l+1:m) = .5*( yeven - omegap .* yodd);
end