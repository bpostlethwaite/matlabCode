function y = forsub (A,b,p)
%
% function y = forsub (A,b,p)
%
% Given a unit lower triangular, nonsingular n by n matrix A, 
% an n-vector b, and a permutation p,
% return vector y which solves Ay = Pb

n = length(b); 

% permute b according to p
b = b(p);

%forward substitution
y = b;
for k = 2:n
  y(k) = b(k) - A(k,1:k-1) * y(1:k-1);
end
