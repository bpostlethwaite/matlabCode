function [c, s] = givens(a, b)
%
%   [c, s] = givens(a, b);
%
%  Compute a Givens rotation to zero an element:
%     / c  -s \ / a \ _ / r \
%     \ s   c / \ b / - \ 0 /
%

%  J. Nagy  03-05-02

if b == 0
  c = 1;
  s = 0;
else
  if abs(b) > abs(a)
    tau = -a/b;
    s = 1 / sqrt(1 + tau*tau);
    c = s * tau;
  else
    tau = -b/a;
    c = 1 / sqrt(1 + tau*tau);
    s = c * tau;
  end
end