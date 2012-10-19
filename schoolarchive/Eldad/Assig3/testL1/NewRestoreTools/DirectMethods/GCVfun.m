function V = GCVfun(alpha, ghat, s, delta)
%
%    V = GCVfun(alpha, ghat, s);
%
%  This function evaluates the GCV function for Tikhonov
%  regularization.  
%
%  Input:  alpha -  regularization parameter
%           ghat -  vector U'*b, where U = left singular vectors
%              s -  vector containing the singular values
%          delta -  vector containing singular values of the regularization
%                   operator (we assume that the matrix and regularization 
%                   operator have the same singular vectors)
%
%  Output:     V -  the scalar V(alpha).
%

% J. Nagy,  5/15/02

% Modifications:
% 6/6/02, J. Nagy
%         Modifed so that Tikhonov will (in some cases) work for 
%         derivative regularization operators

n = length(ghat);

s2 = s .^ 2;
delta2 = delta .^2;
alpha2 = alpha*alpha*delta2;

t1 = delta2 ./ (s2 + alpha2);
t2 = (ghat .* t1).^2;

V = ( sum(t2) / n ) / ( sum(t1) / n )^2;