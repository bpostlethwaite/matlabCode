function [x, alpha] = Tikhonovsolver(U, S, V, b, options)
%
%         [x, alpha] = Tikhonovsolver(U, S, V, b, options)
%
%  This function computes a Tikhonov regularized LS solution to the
%  PROJECTED problem, using the identity matrix as the regularization operator.
%
%  Input: 
%    U, S, V - SVD of A, where S is a column vector containing
%              singular (or spectral) values of A.
%          b - noisy image
%    options - structure with fields.  The following fields are needed:
%         RegPar: [value | GCV | {modGCV} | optimal] where
%                     (a) a value for the regularization parameter
%                     (b) a method for choosing a reg. parameter
%                                (i) 'GCV' - standard GCV
%                                (ii) 'modGCV' - modified GCV (default)
%                     (c) 'optimal' - find the optimal reg. parameter
%                                (requires x_true)
%         Omega: if RegPar is 'modGCV', then omega must be 
%                           [value | {adapt}]
%                     (a) a value for omega
%                     (b) 'adapt' - use the adaptive method (default)
%       Vx_true: V_k * x_true (where V_k comes from Lanczos bidiagonlization)
%                 This vector is needed to compute the optimal parameter.
%
%  Output: x - restored image
%      alpha - regularization parameter

bhat = U'*b;
bhat = bhat(:);

alpha = HyBRget(options,'RegPar',[],'fast');
omega = HyBRget(options,'Omega',[],'fast');

switch alpha
  case 'gcv'
    alpha = fminbnd('TikGCV', 0, 1, [], bhat, S);
  case 'modgcv'
    alpha = fminbnd('TikGCV', 0, 1, [], bhat, S, omega);
  case 'optimal'
    Vx_true = HyBRget(options, 'Vx');
    [e, alpha] = OptTikTol(U, S, V, b, Vx_true);
  otherwise
    if ~isa(alpha, 'double')
      error('Incorrect input argument for alpha.')
    end
end

D = abs(S).^2 + alpha^2;
bhat = conj(S) .* bhat(1:length(S));
xhat = bhat ./ D;
x = V * xhat;


%% -----------------------SUBFUNCTIONS----------------------------------
function [e, alpha] = OptTikTol(U, S, V, b, x)
%
%   [e, alpha] = OptTolTol(U, S, V, b, x)
%
%   This function computes the optimal regularization parameter
%       for Tikhonov regularization by minimizing:
%           || filtered solution - true solution ||_2
%   
%   Assume the problem has the form b = Ax + noise, with [U, S, V] = svd(A).
%
%   Input:
%       [U, S, V] - SVD of A, where S is a column vector containing
%                       singular or spectral values of A
%               b - noisy data
%               x - true data
%
%   Output:
%               e - relative error
%           alpha - optimal regularization parameter

bhat = U'*b;
k = length(S);
alpha = fminbnd('TikRelErrors',0,1,[],bhat(1:k),V'*x,S);
e = TikRelErrors(alpha, bhat(1:k), V'*x, S);

