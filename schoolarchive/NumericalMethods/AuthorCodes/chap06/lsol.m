function [x,nr] = lsol(b,AA,p)
%
% function [x,nr] = lsol(b,AA,p)
%
% Given the output AA and p of house, which contain Q and R of matrix A,
% and given a right-hand-side vector b,  solve min || b - Ax ||.
% Return also the norm of the residual, nr = ||r|| = ||b - Ax||.

y = b(:); [m,n] = size(AA);
% transform b
for k=1:n
  u = [p(k);AA(k+1:m,k)];
  y(k:m) = y(k:m) - 2*u *(u'*y(k:m));
end
% form upper triangular R and solve
R = triu(AA(1:n,:));
x = R \ y(1:n); nr = norm(y(n+1:m));