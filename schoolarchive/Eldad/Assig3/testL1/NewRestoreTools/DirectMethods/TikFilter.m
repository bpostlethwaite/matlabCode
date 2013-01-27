function x = TikFilter(U, S, V, b, alpha, delta)
%
%         x = TikFilter(U, S, V, b, alpha, delta);
%
%  This function computes a Tikhonov regularized LS solution
%  using the identity matrix as the regularization operator.
%
%  Input: 
%    U, S, V - SVD of A, where S is a column vector containing
%              singular (or spectral) values of A.
%      alpha - can be either:
%               *  number between 0 and 1, which is used for
%                  the regularization parameter.
%               * 'gcv' in which case the GCV method is used
%                  to compute the regularization parameter.
%               *  'dp' in which case the discrepancy principle
%                  is used to compute the regularization parametner.
%              the default is to use 'gcv'
%      delta - norm of the noise, if it is known (must be specified
%              if 'dp' is used).  Otherwise it is not used.
%
%  Output: x - restored image
%

if nargin == 4
  alpha = 'gcv';
end

[m, n] = size(b);

bhat = U'*b;
bhat = bhat(:);

switch alpha

case 'gcv'
  alpha = fminbnd('TikGCVfun', 0, 1, [], bhat, S);

case 'dp'
  error('Discrepancy principle has not yet been implemented.')

otherwise
  if ~isa(alpha, 'double')
    error('Incorrect input argument for alpha.')
  end

end

D = abs(S).^2 + alpha^2;
bhat = conj(S) .* bhat;
xhat = bhat ./ D;
xhat = reshape(xhat, m, n);
x = V * xhat;
