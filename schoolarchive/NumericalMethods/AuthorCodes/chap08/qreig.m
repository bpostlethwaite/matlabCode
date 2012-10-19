function [lambda,itn] = qreig (A,tol)
%
% function [lambda,itn] = qreig (A,Tol)
%
% Find all real eigenvalues lambda of A 
% Return also iteration counters in itn

% First stage, bring to upper Hessenberg form
A = houseeig(A);

% second stage: deflation loop
n = size(A,1); lambda = []; itn = [];
for j = n:-1:1
  % find jth eigenvalue
  [lambda(j),itn(j),A] = qrshift (A(1:j,1:j),tol);
end