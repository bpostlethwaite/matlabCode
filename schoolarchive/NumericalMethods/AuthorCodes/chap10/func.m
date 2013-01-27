function f = func(x)

which = 2;
if which == 1 % Runge example
    f = 1./(1+25*x.^2);
    
elseif which == 2  % oscillatory
    f = exp(3*x).*sin(200*x.^2)./(1+20*x.^2);
    
elseif which == 3  % hat function
    f = x;
    I = find(x>pi);
    f(I) = 2*pi - x(I);
    
 elseif which == 4  % examp13.3
    f = x.^2 .* (x+1).^2 .* (x-2).^2 - ...
      exp(-x.^2) .* sin(x+1).^2 .* sin(x-2).^2;
     
 elseif which == 5  % examp13.4
    r = 2;
    f = x.^2 .* (x+1).^r .* (x-2).^r - ...
      exp(-x.^2) .* sin(x+1).^r .* sin(x-2).^r;
  
  elseif which == 6  % examp13.1
     f = cos(3*x) - .5*sin(5*x) + .05*cos(104*x);
end