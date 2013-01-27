function [f, alpha] = Tikhonov(varargin)
%
%         f = Tikhonov(K, g, alpha, 'deriv');
%
%         f = Tikhonov(U, s, V, g, alpha, 'deriv');
%
%  You can also use [f, alpha] = ... to return the regularization parameter.
%
%  This function computes a Tikhonov regularized LS solution
%  using the identity matrix as the regularization operator.
%
%  Input: 
%          K - psfMatrix
%          g - blurred image
%    U, s, V - SVD of K, where s is a column vector containing
%              singular (or spectral) values of K.
%
%  Optional Input:
%      alpha - regularization parameter.  If it is not specified,
%              then the GCV method is used to compute a parameter.
%
%  Output: f - restored image
%
%  Optional Output: alpha - regularization parameter
%

%  J. Nagy, 6/4/02
%           Note:  This function replaces TLS.m in version 1.0.

deriv = 0;
if (length(varargin) == 2)
  K = varargin{1};
  g = varargin{2};
  [U, s, V] = svd(K, g);
  alpha = [];
elseif (length(varargin) == 3)
  K = varargin{1};
  g = varargin{2};
  if ischar(varargin{3})
    deriv = 1;
    alpha = [];
  else
    alpha = varargin{3};
  end
  [U, s, V] = svd(K, g);
elseif (length(varargin) == 4)
  if ischar(varargin{4})
    deriv = 1;
    K = varargin{1};
    g = varargin{2};
    alpha = varargin{3};
    [U, s, V] = svd(K, g);
  else
    U = varargin{1};
    s = varargin{2};
    V = varargin{3};
    g = varargin{4};
    alpha = [];
  end
elseif (length(varargin) == 5)
  U = varargin{1};
  s = varargin{2};
  V = varargin{3};
  g = varargin{4};
  if ischar(varargin{5})
    deriv = 1;
    alpha = [];
  else
    alpha = varargin{5};
  end
elseif (length(varargin) == 6)
  U = varargin{1};
  s = varargin{2};
  V = varargin{3};
  g = varargin{4};
  alpha = varargin{5};
  deriv = 1;
else
  error('incorrect number of inputs')
end

[m, n] = size(g);

ghat = U'*g;
ghat = ghat(:);

if deriv
  if isa(U, 'transformMatrix')
    switch U.transform
    case 'fft'
      delta = svd(psfMatrix([0 -1 0; -1 4 -1; 0 -1 0], 'periodic'), g);
    case 'dct'
      delta = svd(psfMatrix([0 -1 0; -1 4 -1; 0 -1 0], 'reflexive'), g);
    otherwise
      error('Incorrect field type')
    end
  else
    error('The derivative regularization operator cannot be used with this blurring function. Use an iterative method.')
  end
else
  delta = ones(size(ghat));
end

if isempty(alpha)
  alpha = fminbnd('GCVfun', 0, 1, [], ghat, s, delta);
%  disp(sprintf('GCV chooses reg. param., alpha = %f', alpha))
end

fhat = normalDLS(s, ghat, alpha, delta);
fhat = reshape(fhat, m, n);
f = V * fhat;