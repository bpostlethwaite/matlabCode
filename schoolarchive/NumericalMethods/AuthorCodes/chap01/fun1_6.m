function f = fun1_6(x,n)
%
% integrand for evaluating the integral of Example 1.6 stably
f = x.^n ./ (x+10);