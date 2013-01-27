function f = phi(x)
% objective function for Examples 9.5, 9.6

number = 3;

if number == 2
  f = .5 * ((1.5 - x(1)*(1-x(2)))^2 + (2.25 - x(1)*(1-x(2)^2))^2 + ...
       (2.625 - x(1)*(1-x(2)^3))^2 );
elseif number == 3
  f = x(1)^4 + x(1)*x(2) + (1+x(2))^2;
end