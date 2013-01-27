function f = ffun(t,y)
%
% f = ffun(t,y)


isw = 2;

if isw == 1
    
  % ODE function for Henon-Heiles  

  f = zeros(4,1);
  f(1) = y(3);
  f(2) = y(4);
  f(3) = -y(1) - 2*y(1)*y(2);
  f(4) = -y(2) - y(1)^2 + y(2)^2;

elseif isw == 2
    
  % ODE function for FPU
  
  m = 3; omega = 100;
  f = zeros(4*m,1);
  f(1:2*m) = y(2*m+1:4*m);
  f(3*m+1:4*m) = -omega^2 * y(m+1:2*m);
  f(2*m+1) = f(2*m+1) - (y(1)-y(m+1))^3;
  f(3*m+1) = f(3*m+1) + (y(1)-y(m+1))^3;
  f(3*m) = f(3*m) - (y(m)+y(2*m))^3;
  f(4*m) = f(4*m) - (y(m)+y(2*m))^3;
  f(2*m+2:3*m) = f(2*m+2:3*m) -...
      (y(2:m)-y(1:m-1)-y(m+2:2*m)-y(m+1:2*m-1)).^3;
  f(2*m+1:3*m-1) = f(2*m+1:3*m-1) + ...
      (y(2:m)-y(1:m-1)-y(m+2:2*m)-y(m+1:2*m-1)).^3;
  f(3*m+1:4*m-1) = f(3*m+1:4*m-1) + ...
      (y(2:m)-y(1:m-1)-y(m+2:2*m)-y(m+1:2*m-1)).^3;
  f(3*m+2:4*m) = f(3*m+2:4*m) + ...
      (y(2:m)-y(1:m-1)-y(m+2:2*m)-y(m+1:2*m-1)).^3;
end
