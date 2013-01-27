function [f,Jac] = funcode(u)
%
% function [f,Jac] = funcode(u)
%
% function and Jacobian for boundary value ode example
% size and extension for u
J = length(u);
h = 1/(J+1); h2 = h^2;
u = shiftdim(u);
ue = [0;u;0];

f = (ue(1:J) - 2*ue(2:J+1) + ue(3:J+2))/h2 + exp(ue(2:J+1));

% Do Jacobian in a non-sparse way for clarity
Jac = zeros(J,J);
for j = 1:J
    Jac(j,j) = -2/h2 + exp(u(j));
    if j > 1 , Jac(j,j-1) = 1/h2; end
    if j < J , Jac(j,j+1) = 1/h2; end
end
        