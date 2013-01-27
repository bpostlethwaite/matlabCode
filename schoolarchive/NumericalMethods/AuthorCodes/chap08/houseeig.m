function A = houseeig(A)
%
% function A = houseeig(A)
%
% Reduce A to upper Hessenberg form using Householder reflections
n = size(A,1);
for k = 1:n-2
  z=A(k+1:n,k);
  e1=[1; zeros(n-k-1,1)];
  u=z+sign(z(1))*norm(z)*e1; 
  u = u/norm(u);
  % multiply from left and from right by Q = eye(n-k)-2*u*u';
  A(k+1:n,k:n) = A(k+1:n,k:n) - 2*u*(u'*A(k+1:n,k:n));    
  A(1:n,k+1:n) = A(1:n,k+1:n) - 2*(A(1:n,k+1:n)*u)*u';
end