function fhat = normalDLS(s, ghat, alpha, delta)
%
%     fhat = normalDLS(s, ghat, alpha, delta);
%
%  This function solves the damped least squares problem:
%
%    minimize  || / ghat \ _ /    S   \ fhat ||
%      fhat    || \  0   /   \ alpha*D/      ||
%                                              2
%  where S and D are diagonal matrices.  We use normal equations
%  to solve this.
%
%  Input:  s  -  column vector containing the diagonal entries
%                of S.
%       ghat  -  column vector
%      alpha  -  scalar
%      delta  -  column vector containing the diagonal entries of D.
%
%  Output: 
%       fhat  -  solution vector
%

%  J. Nagy   01-11-03

D = conj(s).*s + alpha*alpha*( conj(delta).*delta );
ghat = conj(s) .* ghat;

fhat = ghat ./ D;