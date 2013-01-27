function x = backsub (A,b)
%
% function x = backsub (A,b)
%
% Given an upper triangular, nonsingular n by n matrix A and
% an n-vector b, return vector x which solves Ax = b

n = length(b); x = b;
x(n) = b(n) / A(n,n);
for k = n-1:-1:1
  x(k) = ( b(k) - A(k,k+1:n)*x(k+1:n) ) / A(k,k);
end
