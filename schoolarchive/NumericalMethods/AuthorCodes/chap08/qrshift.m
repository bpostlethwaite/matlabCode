function [lam,iter,A] = qrshift (A,tol)
%
% function [lam,iter,A] = qrshift (A,tol)
%
% Find one eigenvalue lam of A in upper Hessenberg form, 
% return iteration count, too. Also improve A for future 

m = size(A,1); lam = A(m,m); iter=0; I = eye(m);
if m == 1, return, end
while (iter < 100) % max number of iterations
  if (abs(A(m,m-1)) < tol), return, end    % check convergence  
  iter=iter+1;
  [Q,R]=qr(A-lam*I);           % compute the QR decomposition
  A=R*Q+lam*I;                 % find the next iterate 
  lam = A(m,m);                % next shift
end