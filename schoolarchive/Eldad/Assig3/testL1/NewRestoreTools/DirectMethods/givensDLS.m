function fhat = givensDLS(s, ghat, alpha, delta)
%
%     fhat = givensDLS(s, ghat, alpha, delta);
%
%  This function solves the damped least squares problem:
%
%    minimize  || / ghat \ _ /    S   \ fhat ||
%      fhat    || \  0   /   \ alpha*D/      ||
%                                              2
%  where S and D are diagonal matrices.  We use Givens rotations
%  to do this efficiently.
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

%  J. Nagy   06-05-02

n = length(ghat);
delta = alpha*delta;
for i = 1:n
  [cn, sn] = givens(s(i), delta(i));
  s(i) = cn*s(i) - sn*delta(i);
  ghat(i) = cn*ghat(i);
end
fhat = ghat ./ s;