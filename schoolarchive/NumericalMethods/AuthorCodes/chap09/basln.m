function [x,y,s] = basln(A,b,c,sort)
%
% function [x,y,s] = basln(A,b,c,sort)
%
% given a vector of indices, the first l indicate a basis out of A.
% Find corresponding basic solution.

[l,m] = size(A);
B = zeros(l,l); cb = zeros(l,1);

% construct basis
for j=1:l
  B(:,j) = A(:,sort(j));
  cb (j) = c(sort(j));
end 

xb = B \ b;
x = zeros(m,1);
for j=1:l
  x(sort(j)) = xb(j);
end

y = B' \ cb;
s = c - A'*y;