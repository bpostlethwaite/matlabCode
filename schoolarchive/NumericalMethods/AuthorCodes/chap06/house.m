function [A,p] = house(A)
% function [A,p] = house(A)
%
% perform QR decomposition using Householder reflections
% Transformations are of the form P_k = I - 2u_k(u_k^T), so 
% store effecting vector u_k in p(k) + A(k+1:m,k). Assume m > n. 

[m,n]=size(A); p = zeros(1,n);
for k = 1:n
  % define u of length = m-k+1
  z = A(k:m,k);
  e1 = [1; zeros(m-k,1)];
  u = z+sign(z(1))*norm(z)*e1; u = u/norm(u);
  % update nonzero part of A by I-2uu^T
  A(k:m,k:n) = A(k:m,k:n)-2*u*(u'*A(k:m,k:n));
  % store u
  p(k) = u(1);
  A(k+1:m,k) = u(2:m-k+1);
end