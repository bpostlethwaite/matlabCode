function [x,k] = newtons(f,x,tol,nmax)
%
% function [x,k] = newtons(f,x,tol,nmax)
%
% This function returns in x a column vector x_k such that
%      || x_k - x_{k-1} || < tol (1 + ||x_k||)
% and in k the number of iterations (Jacobian evaluations) required.
% On entry, x contains an initial guess.
% If k equals nmax then no convergence has been reached.
% 
% The iterates ||f(x_k)|| are recorded. This option
% can be easily turned off

%Initialize
x = x(:); % ensure x is a column vector
fprintf ('k       ||f(x_k)|| \n')
format long g

%Newton
for k=1:nmax
  [fx,Jx] = feval(f,x);
  fprintf ('%d     %e     \n',k-1,norm(fx) )
  p = -Jx \ fx;
  x = x + p;  
  if norm(p) < tol*(1+norm(x))
    fx = feval(f,x);   
    fprintf ('%d     %e     \n',k,norm(fx) ) 
    return
  end
end
k = nmax;