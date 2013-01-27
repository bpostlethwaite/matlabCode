function [f,J] = func(x)
%
%  Function and Jacobian for nonlinear system

f(1) = x(1)^2 - 2*x(1) - x(2) + 1;
f(2) = x(1)^2 + x(2)^2 - 1;
f = shiftdim(f);

J(1,1) = 2*x(1) - 2;
J(1,2) = -1;
J(2,1) = 2*x(1);
J(2,2) = 2*x(2);