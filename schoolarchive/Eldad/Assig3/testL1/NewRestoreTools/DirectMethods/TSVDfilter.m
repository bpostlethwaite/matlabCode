function x = TSVDfilter(U, s, V, b, tol, delta)
%
%           x = TSVDfilter(U, s, V, b, tol, delta);
%
%  Compute truncated SVD solution.
%
%  On Entry:
%     U, V  - Unitary matrices
%        s  - column vector such that
%              A = U*diag(s)*V'
%        b  - right hand side IMAGE (not vector)
%      tol  - can be either:
%              * number between 0 and 1, in which case truncation will
%                occur if |s| < tol.
%              * 'gcv' in which case the GCV method is used to compute
%                 the regularization parameter.
%              * 'dp' in which case the discrepancy principle is used to
%                 compute the regularization parameter.
%             the default is to use 'gcv'
%     delta - norm of the noise, if it is known (must be specified if
%             'dp' is used).  Otherwise it is not used.
%  On Exit:
%           x - tsvd solution of Ax=b
%

if nargin == 4
  tol = 'gcv';
end

bhat = U'*b;

switch tol

case 'gcv'
  tol = GCVforSVD2(s, bhat(:));

case 'dp'
  error('Discrepance principle has not yet been implemented.')

otherwise
  if ~isa(tol, 'double')
    error('Incorrect input argument for tol.')
  end
  tol = max(s(:))*tol;

end

s = reshape(s, size(b));
mask = (abs(s) >= tol);
mask2 = ~mask;
s = (mask .* (1./(s + mask2)) );

x = V * (s .* bhat);

if isreal(b)
  x = real(x);
end








