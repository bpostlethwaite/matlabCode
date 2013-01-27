function [x,iter,res] = pcgmg (A,b,x0,tol)
%
% function [x,iter,res] = pcgmg (A,b,x0,tol)
% 
% CG with V-cycle preconditioning; assume x0=0

N = sqrt(length(b));
flevel = log2(N+1);
tol2 = tol^2;
x = x0;
r = b - A*x0;
h = poismg(A,r,x0,flevel);
d = r'*h; bb = d;
p = h;

iter = 0;
while d > tol2 * bb
  do = d;
  s = A*p;
  alfa = d / (p'*s);
  x = x + alfa*p;
  r = r - alfa*s;
  h = poismg(A,r,x0,flevel);
  d = r'*h;
  p = h + d/do*p;
  iter = iter + 1;
  res(iter) = norm(r);
end