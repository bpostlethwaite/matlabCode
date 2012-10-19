function yh = myfft(y)
%
%  function yh = myfft(y)
%
%  Apply FFT, assuming the length m of y is a power of 2 and even.
m = length(y);
if m == 2
    % directly define coefficients
    yh(1) = .5*(y(1)+y(2)); yh(2) = .5*(y(1)-y(2));
    yh(1) = y(1)+y(2); yh(2) = y(1)-y(2);
else
    % recursion step
    l = m/2;
    yeven = y(1:2:m-1); yodd = y(2:2:m);
    ceven = myfft(yeven); codd = myfft(yodd);
    omegam = exp(-i*2*pi/m); omegap = omegam.^(0:1:l-1);
    yh(1:l) = ceven + omegap .* codd;
    yh(l+1:m) = ceven - omegap .* codd;
end