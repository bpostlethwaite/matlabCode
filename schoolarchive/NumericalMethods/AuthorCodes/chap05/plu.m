function [p,A] = plu (A)
%
% function [p,A] = plu (A)
%
% Perform LU decomposition with partial pivoting.
% Upon return the coefficients of L and U replace those 
% of the input n-by-n nonsingular matrix A. The row interchanges
% performed are recorded in the 1D array p.

n = size(A,1);

% initialize permutation vector p
p = 1:n;

% LU decomposition with partial pivoting
for k = 1:n-1

  % find row index of relative maximum in column k
  [val,q] = max(abs(A(k:n,k))); 
  q = q + k-1; 

  % interchange rows k and q and record this in p 
  A([k,q],:)=A([q,k],:); 
  p([k,q])=p([q,k]);

  % compute the corresponding column of L
  J=k+1:n;
  A(J,k) = A(J,k) / A(k,k);
  
  % update submatrix by outer product
  A(J,J) =  A(J,J) - A(J,k) * A(k,J);
  
end

