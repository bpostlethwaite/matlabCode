% Question 2

function [F,J] = buildJacoII(u,D,N,h)
%J = D - h^2*spdiags(exp(u),0,(N)^2,(N)^2);
J = D - h^2*spdiags(exp(u),0,(N-2)^2,(N-2)^2);
F = D*u - h^2*exp(u); 
end


