function [x, trun_index] = TSVDb(U, s, V, b, tol)
%
%           [x, trun_index] = TSVDb(U, s, V, b, tol);
%
%  Compute truncated SVD solution.
%
%  This is experimental, to see if we can make L-curve interactive.
%
%  On Entry:
%     U, V  - Unitary matrices
%        s  - column vector such that
%              K = U*diag(s)*V'
%        b  - right hand side IMAGE (not vector)
%      tol  - tolerance to specify how much regularization to
%             use, 0 <= tol <= 1.  
%
%             Truncation will occur if s values are < trun_tol
%             If trun_tol is not specified then it is set to 0 (i.e,
%             no regularization, which means x = A\b is returned).
%             You can also use the character string 'help' for trun_tol,
%             in which case the generalized cross validation method
%             will be used to pick the truncation tolerance.
%  On Exit:
%           x - tsvd solution of Ax=b
%  trun_index - index of singular value where truncation occurs.
%

%  J. Nagy  03-05-02

%  Modifications:
%  5/25/02, J. Nagy
%           1. Fixed problems with 'help' input.
%           2. Also fixed scaling problems with tol.
%  6/2/02,  J. Nagy
%           Major modifications to allow for using spectral
%           factorizations, in addition to standard singular
%           value decomposition.

[m, n] = size(b);

bhat = U'*b;

if nargin == 4
  s = 1 ./ s;
  s = reshape(s, size(b));
else
  if (ischar(tol))
    if length(tol) > 4
      trun_tol = GCVforSVDb(s, bhat(:));
    else
      trun_tol = GCVforSVD2(s, bhat(:));
    end
  else
    trun_tol = max(s(:)) * tol;
  end
  s = reshape(s, size(b));
  mask = (abs(s) >= trun_tol);
  trun_index = sum(sum(mask));
  mask2 = ~mask;
  s = ( mask .* (1./(s + mask2)) );
end

x = V * (s .* bhat);

if isreal(b)
  x = real(x);
end








