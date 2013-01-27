function [p,n] = bisect(func,a,b,fa,fb,atol)
%
% function [p,n] = bisect(func,a,b,fa,fb,atol)
%
% Assuming fa = func(a), fb = func(b), and fa*fb < 0,
% there is a value root in (a,b) such that func(root) = 0.
% This function returns in p a value such that 
%      | p - root | < atol
% and in n the number of iterations required.
%

% check input
if (a >= b) || (fa*fb >= 0) || (atol <= 0)
  disp(' something wrong with the input: quitting');
  p = NaN; n=NaN;
  return
end

n = ceil ( log2 (b-a) - log2 (2*atol));
for k=1:n
  p = (a+b)/2;
  fp = feval(func,p);
  if abs(fp) < eps, n = k; return, end
  if fa * fp < 0
    b = p;
  else
    a = p;
    fa = fp;
  end
end
p = (a+b)/2;
