function G = TikGCVfun(alpha, bhat, s)
%
%    G = TikGCVfun(alpha, bhat, s);
%
%  This function evaluates the GCV function for Tikhonov
%  regularization.  
%
%  Input:  alpha -  regularization parameter
%           bhat -  vector U'*b, where U = left singular vectors
%              s -  vector containing the singular values
%
%  Output:     G -  the scalar G(alpha).
%

n = length(bhat);

s2 = abs(s) .^ 2;
alpha2 = alpha^2;

t1 = 1 ./ (s2 + alpha2);
t2 = abs(bhat .* t1) .^2;

G = ( sum(t2) / n ) / ( sum(t1) / n )^2;
